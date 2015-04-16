/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.hammerlab.guacamole.{ DelayedMessages, Concordance, DistributedUtil, Common, SparkCommand }
import org.hammerlab.guacamole.Common.Arguments.StructuralVariantArgs
import org.hammerlab.guacamole.ReadSet
import org.hammerlab.guacamole.reads.MappedRead
import org.hammerlab.guacamole.likelihood.Likelihood
import org.hammerlab.guacamole.variants.{ AlleleEvidence, Genotype, AlleleConversions, CalledAllele }

import org.hammerlab.guacamole.filters.GenotypeFilter.GenotypeFilterArguments
import org.hammerlab.guacamole.filters.PileupFilter.PileupFilterArguments
import org.hammerlab.guacamole.filters.{ GenotypeFilter, QualityAlignedReadsFilter }
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.Read
import org.kohsuke.args4j.Option

import scala.math

/**
 * Simple Bayesian variant caller implementation that uses the base and read quality score
 */

class AlignmentPair(
    val firstRead: MappedRead,
    val secondRead: MappedRead,
    lociInterval: Long = 25) {

  def insertSize(): Long = secondRead.end - firstRead.start
  def alignmentScore(): Int = firstRead.alignmentQuality + secondRead.alignmentQuality

  def loci(): Iterable[Long] = {
    val start = firstRead.start
    val end = secondRead.end

    var lociStart = (start / lociInterval) * lociInterval
    if ((start - lociStart) > lociInterval / 2) lociStart += lociInterval

    var lociEnd = ((end / lociInterval) + 1) * lociInterval
    if ((lociEnd - end) > lociInterval / 2) lociEnd -= lociInterval

    // val lociRange: NumericRange = new NumericRange(lociStart, lociEnd, lociInterval.toLong, true)
    lociStart to lociEnd by lociInterval.toLong
  }

  def lociAlignmentPairs: Iterable[(Long, AlignmentPair)] = {
    loci.map(location => (location, this))
  }

}

class AlignmentPairList(
    readName: String,
    firstReadList: Iterable[MappedRead],
    secondReadList: Iterable[MappedRead],
    filter: (MappedRead, MappedRead) => Boolean = ((m1, m2) => true)) {

  def alignmentPairs(): Iterable[AlignmentPair] = {
    for {
      read1 <- firstReadList
      read2 <- secondReadList
      if (filter(read1, read2))
    } yield {
      if (read1.start < read2.start) new AlignmentPair(read1, read2)
      else new AlignmentPair(read2, read1)
    }
  }

  def lociAlignmentPairs(): Iterable[(Long, AlignmentPair)] =
    alignmentPairs
      .flatMap(alignmentPair => alignmentPair.lociAlignmentPairs)

}

object StructuralVariant {

  protected class Arguments extends StructuralVariantArgs {

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "structural-variant"
    override val description = "TBD"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      // val readSet: ReadSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(mapped = true, nonDuplicate = true, passedVendorQualityChecks = true, isPaired = true))
      val readSet: ReadSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(mapped = true, nonDuplicate = true))
      readSet.mappedReads.persist()

      Common.progress(
        "Loaded %,d mapped non-duplicate reads into %,d partitions.".format(readSet.mappedReads.count, readSet.mappedReads.partitions.length))

      val mappedReads: RDD[MappedRead] = readSet.mappedReads
      // mappedReads.take(100).foreach(println)

      val readPairs: RDD[(String, Iterable[MappedRead])] =
        readSet
          .mappedReads
          .groupBy(read => read.readName)

      // readPairs.take(10).foreach(pair => {
      //   println("\nPrinting " + pair._1)
      //   pair._2.foreach(read => {
      //     println(read)
      //     read.matePropertiesOpt match {
      //       case Some(mate) => {
      //         println(mate)
      //       }
      //       case None => {}
      //     }

      //   })
      // })

      val readAlignmentPairLists: RDD[(String, AlignmentPairList)] =
        readPairs
          .map(pair => {
            (pair._1, new AlignmentPairList(pair._1,
              pair._2.filter(read => read.matePropertiesOpt.get.isFirstInPair),
              pair._2.filter(read => !read.matePropertiesOpt.get.isFirstInPair),
              ((m1, m2) => {
                val insertSize = scala.math.abs(m2.end - m1.start)

                if (insertSize > 25000) false
                else if (!m1.isPositiveStrand || m2.isPositiveStrand) false
                else true
              })
            ))

          })

      readAlignmentPairLists.take(10).foreach(pair => {

        println("\nPrinting " + pair._1)
        pair._2.alignmentPairs.foreach(alignmentPair => {
          println(alignmentPair.firstRead)
          println(alignmentPair.secondRead)
          println(alignmentPair.loci)
        })

      })

      println("Filtered values = " + (readPairs.count - readAlignmentPairLists.count))

      val lociAlignmentPairsRDD: RDD[(Long, AlignmentPair)] =
        readAlignmentPairLists
          .map(pair => pair._2)
          .flatMap(alignmentPairList => alignmentPairList.lociAlignmentPairs)

      val groupedLociAlignmentRDD: RDD[(Long, Iterable[AlignmentPair])] =
        lociAlignmentPairsRDD
          .groupBy(pair => pair._1)
          .map(pair => (pair._1, pair._2.map(lap => lap._2)))

      groupedLociAlignmentRDD.take(10).foreach(lap => {

        println("Loci = " + lap._1)
        lap._2.foreach(ap => {
          println("******************************************")
          println(ap.firstRead)
          println(ap.secondRead)
        })

      })
    }

  }
}
