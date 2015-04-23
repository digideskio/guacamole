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
import org.hammerlab.guacamole.AlignmentPairList
import org.hammerlab.guacamole.AlignmentPair
import org.hammerlab.guacamole.SVFeatures
//https://commons.apache.org/proper/commons-math/javadocs/api-3.3/org/apache/commons/math3/distribution/fitting/MultivariateNormalMixtureExpectationMaximization.html
// import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization

/**
 * Simple Bayesian variant caller implementation that uses the base and read quality score
 */

// object SVFeatures {

//   def apply(
//     alignmentPairs: Iterable[AlignmentPair], 
//     insertSizeMean: Long = 300, 
//     insertSizeSD: Long = 30, 
//     maxIterations: Int = 10,
//     convergenceThreshold: Double = 0.0001) {

//     def lognormal: Double(y: Double, mu: Double, sigma: Double) = {
//       -1.0 * Math.log(sigma + Math.sqrt(2.0 * Math.PI)) + -1.0 / 2.0 * Math.pow(((y - mu) / sigma), 2)
//     }

//     def logsumexp: Double(xArray: Array[Double]) = {
//       val m = xArray.max
//       val sumexp = 
//         xArray
//         .map(x => Math.exp(x-m))
//         .sum
//       m + Math.log(sumexp)
//     }

//     def likelihood(
//       yArray: Array[Double], 
//       wArray: Array[Double],
//       muArray: Array[Double],
//       sigma: Double): Double {

//       val sumWeightedLikelihoods: Double = 
//         yArray
//         .map(y => {
//           val weightedPointLikelihoods = 
//             wArray
//             .zip(muArray)
//             //w + lognormal(y, mu, sigma)
//             .map(pair => pair._1 + lognormal(y, pair._2, sigma))

//           logsumexp(weightedPointLikelihoods)
//         })
//         .sum

//       sumWeightedLikelihoods / yArray.length
//     }

//     def gamma(
//       yArray: Array[Double], 
//       wArray: Array[Double],
//       muArray: Array[Double],
//       sigma: Double) Array[Array[Double]] = 
//         yArray
//         .map( y => {
//           val likelihoodRow = 
//             wArray
//             .zip(muArray)
//             .map(pair => pair._1 + lognormal(y, pair._2, sigma))

//           val total = logsumexp(likelihoodRow)

//           likelihoodRow.map(l => l - total)
//         })

//     def updateW(
//       yArray: Array[Double],
//       nArray: Array[Double]): Array[Double] =
//       nArray.map(n => n - Math.log(yArray.length))

//     def updateMuForComponent(
//       gArray: Array[Array[Double]], 
//       yArray: Array[Double],
//       nArray: Array[Double],
//       component: Int): Double = {

//       val lse = 
//         gArray
//         .zip(yArray)
//         .map(pair => pair._1(component) + Math.log(pair._2))

//       val numerator = logsumexp(lse)
//       Math.exp(numerator - n(component))
//     }

//     def calculateN(gArray: Array[Array[Double]]): Array[Double] = {
//       //transpose
//       gArray
//       .foldLeft(Array[Array[Double]]())((acc, row) => {
//         acc.zip(row).map( pair => {
//           pair._1 :+ pair._2
//         })
//       })
//       .map(a => logsumexp(a))
//     }

//     def emStep(
//       yArray: Array[Double], 
//       wArray: Array[Double],
//       muArray: Array[Double],
//       sigma: Double) {
//       val g = gamma(yArray, wArray, muArray, sigma)
//       val nArray = calculateN(g)

//       val wPrime = nArray.map(n => n - Math.log(yArray.length))
//       val mu1Prime = updateMuForComponent
//       val muPrime = muArray.clone
//       muPrime(1) = mu1Prime

//       (wPrime, muPrime)
//     }

//     // Default mean fragment size and standard deviation of library are
//     // INSERT_SIZE=300
//     // INSERT_SIZE_SD=30
//     // These should be command line parameters

//     //Let Y = y1,y2...ym be the observed insert sizes
//     //Library has mean fragement size of mu with standard deviation sigma
//     //we initialize the two components to have means mu and yBar with standard devs of sigma
//     //set alpha to .5

//     //In each E step

//     //fix this later 
//     //filters include nearest neighbors and adaptave quality
//     // val filteredAlignmentPairs: Array[AlignmentPair] = alignmentPairs.toArray

//     def alignmentPairFilter(ap: AlignmentPair): Boolean = {
//       true
//     }

//     val yArray:Array[Double] = 
//       alignmentPairs
//       .filter(alignmentPairFilter)
//       .toArray
//       .map(alignmentPair => alignmentPair.insertSize.toDouble)

//     var i=0;
//     var w:Array[Double] = Array(Math.log(.5), Math.log(.5))
//     var mu:Array[Double] = Array(insertSizeMean, (yArray.sum / yArray.length))
//     val sigma = insertSizeSD
//     var lprime:Double = likelihood(yArray, w, mu, sigma)
//     var l: Double = Double.PositiveInfinity
//     while (Math.abs(l - lprime) < convergenceThreshold || i < maxIterations) {

//       l = lprime

//       val pair = emStep(yArray, w, mu, sigma)

//       w = pair._1
//       mu = pair._2

//       lprime = likelihood(yArray, w, mu, sigma) 

//       i+=1
//     }

//     val nodelOneComponentLikelihood = likelihood(yArray, Array[Double](Math.log(1)), Array[Double](initialMu1), sigma)
//     new SVFeatures(Math.exp(w(0)), mu(1), l - nodelOneComponentLikelihood)
//   }
// }

// class SVFeatures(
//   val w0: Double,
//   val mu2: Double,
//   val lrHeterozygous: Double) {

// }

// class AlignmentPair(
//     val firstRead: MappedRead,
//     val secondRead: MappedRead,
//     lociInterval: Long = 25) {

//   def insertSize(): Long = secondRead.end - firstRead.start
//   def alignmentScore(): Int = firstRead.alignmentQuality + secondRead.alignmentQuality

//   def loci(): Iterable[Long] = {
//     val start = firstRead.start
//     val end = secondRead.end

//     var lociStart = (start / lociInterval) * lociInterval
//     if ((start - lociStart) > lociInterval / 2) lociStart += lociInterval

//     var lociEnd = ((end / lociInterval) + 1) * lociInterval
//     if ((lociEnd - end) > lociInterval / 2) lociEnd -= lociInterval

//     // val lociRange: NumericRange = new NumericRange(lociStart, lociEnd, lociInterval.toLong, true)
//     lociStart to lociEnd by lociInterval.toLong
//   }

//   def lociAlignmentPairs: Iterable[(Long, AlignmentPair)] = {
//     loci.map(locus => (locus, this))
//   }

// }

// class AlignmentPairList(
//     readName: String,
//     firstReadList: Iterable[MappedRead],
//     secondReadList: Iterable[MappedRead],
//     filter: (MappedRead, MappedRead) => Boolean = ((m1, m2) => true)) {

//   def alignmentPairs(): Iterable[AlignmentPair] = {
//     for {
//       read1 <- firstReadList
//       read2 <- secondReadList
//       if (filter(read1, read2))
//     } yield {
//       if (read1.start < read2.start) new AlignmentPair(read1, read2)
//       else new AlignmentPair(read2, read1)
//     }
//   }

//   def lociAlignmentPairs(): Iterable[(Long, AlignmentPair)] =
//     alignmentPairs
//       .flatMap(alignmentPair => alignmentPair.lociAlignmentPairs)

// }

object StructuralVariant {

  //Insert Size mean
  //Insert Size StdDev
  //Read Groups?? - only support single

  protected class Arguments extends StructuralVariantArgs {

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "structural-variant"
    override val description = "TBD"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      // val readSet: ReadSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(mapped = true, nonDuplicate = true, passedVendorQualityChecks = true, isPaired = true))
      // val readSet: ReadSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(mapped = true, nonDuplicate = true))
      val readSet: ReadSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(mapped = true))
      readSet.mappedReads.persist()

      // Common.progress(
      //   "Loaded %,d mapped non-duplicate reads into %,d partitions.".format(readSet.mappedReads.count, readSet.mappedReads.partitions.length))

      val mappedReads: RDD[MappedRead] = readSet.mappedReads
      // mappedReads.take(100).foreach(println)

      val readPairs: RDD[(String, Iterable[MappedRead])] =
        readSet
          .mappedReads
          .groupBy(read => read.readName)

      // readPairs.take(100).foreach(pair => {
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

      // readAlignmentPairLists.take(10).foreach(pair => {

      //   println("\nPrinting " + pair._1)
      //   pair._2.alignmentPairs.foreach(alignmentPair => {
      //     println(alignmentPair.firstRead)
      //     println(alignmentPair.secondRead)
      //     println(alignmentPair.loci)
      //   })

      // })

      // println("Filtered values = " + (readPairs.count - readAlignmentPairLists.count))

      val lociAlignmentPairsRDD: RDD[(Long, AlignmentPair)] =
        readAlignmentPairLists
          .map(pair => pair._2)
          .flatMap(alignmentPairList => alignmentPairList.lociAlignmentPairs)

      val groupedLociAlignmentRDD: RDD[(Long, Iterable[AlignmentPair])] =
        lociAlignmentPairsRDD
          .groupBy(pair => pair._1)
          .map(pair => (pair._1, pair._2.map(lap => lap._2)))

      Common.progress("Computed Loci")

      // val groupedLociAlignmentArray: Array[(Long, Iterable[AlignmentPair])] =
      //   groupedLociAlignmentRDD.take(10)

      // Common.progress("Collected Loci Array")

      // groupedLociAlignmentRDD
      //   .sortBy(pair => pair._2.size, false)
      //   .take(10).foreach(lap => {

      //     println("Loci = " + lap._1)
      //     lap._2.foreach(ap => {
      //       println("******************************************")
      //       println(ap.firstRead)
      //       println(ap.secondRead)
      //     })

      //   })

      // Reduce(GenomicLocation, ReadPairInfos)

      val lociSVFeaturesRDD: RDD[(Long, SVFeatures)] =
        groupedLociAlignmentRDD
          .map(pair => (pair._1, SVFeatures(pair._1, pair._2)))

      // val lociSVFeaturesArray: Array[(Long, SVFeatures)] =
      //   groupedLociAlignmentRDD
      //     .take(1000)
      //     .map(pair => (pair._1, SVFeatures(pair._1, pair._2)))

      Common.progress("Computed Loci Features")

      lociSVFeaturesRDD.foreach(pair => println(pair._1 + ": " + pair._2))

    }

  }
}
