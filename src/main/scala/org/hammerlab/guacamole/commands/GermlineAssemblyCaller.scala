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
import org.hammerlab.guacamole.Common.Arguments.{GermlineCallerArgs, SomaticCallerArgs}
import org.hammerlab.guacamole.DistributedUtil.PerSample
import org.hammerlab.guacamole.alignment.DeBrujinGraph
import org.hammerlab.guacamole.filters.GenotypeFilter.GenotypeFilterArguments
import org.hammerlab.guacamole.{ SparkCommand, DelayedMessages, Common, DistributedUtil }
import org.hammerlab.guacamole.likelihood.Likelihood
import org.hammerlab.guacamole.filters.PileupFilter.PileupFilterArguments
import org.hammerlab.guacamole.filters.SomaticGenotypeFilter.SomaticGenotypeFilterArguments
import org.hammerlab.guacamole.filters._
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.{MappedRead, Read}
import org.hammerlab.guacamole.variants.{CalledAllele, AlleleConversions, AlleleEvidence, CalledSomaticAllele}
import org.hammerlab.guacamole.windowing.SlidingWindow
import org.kohsuke.args4j.{ Option => Opt }

/**
 * Simple assembly based germline variant caller
 *
 */
object GermlineAssemblyCaller {

  protected class Arguments extends GermlineCallerArgs with GenotypeFilterArguments {

    @Opt(name = "--kmer-size", usage = "Length of kmer used for DeBrujin Graph assembly")
    var kmerSize: Int = 85

    @Opt(name = "--snv-window-range", usage = "Number of bases before and after to check for additional matches or deletions")
    var snvWindowRange: Int = 20

    @Opt(name = "--min-average-base-quality", usage = "Minimum average of base qualities in the read")
    var minAverageBaseQuality: Int = 20

    @Opt(name = "--min-alignment-quality", usage = "Minimum alignment qualities of the read")
    var minAlignmentQuality: Int = 30

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "germline-assembly"
    override val description = "call germline variants by assembling the surrounding region of reads"

    def discoverHaplotypes(graph: Option[DeBrujinGraph],
                           windows: PerSample[SlidingWindow[MappedRead]],
                           kmerSize: Int,
                           minAlignmentQuality: Int,
                           minAverageBaseQuality: Int): (Option[DeBrujinGraph], Iterator[CalledAllele]) = {

      val newReads = windows(0).newRegions
      val filteredReads = newReads.filter(_.alignmentQuality > minAlignmentQuality)

      val currentGraph: DeBrujinGraph = ???
      val currentReference: Array[Byte] = ???
      val referenceKmerStart = currentReference.slice(kmerSize)

      currentGraph.depthFirstSearch(referenceKmer)

    }

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val readSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(mapped = true, nonDuplicate = true))
      readSet.mappedReads.persist()
      Common.progress(
        "Loaded %,d mapped non-duplicate reads into %,d partitions.".format(readSet.mappedReads.count, readSet.mappedReads.partitions.length))

      val loci = Common.loci(args, readSet)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(args, loci, readSet.mappedReads)

      val genotypes: RDD[CalledAllele] = DistributedUtil.windowFlatMapWithState[MappedRead, CalledAllele, Option[DeBrujinGraph]](
        readSet.mappedReads,
        lociPartitions,
        skipEmpty = true,
        halfWindowSize = args.snvWindowRange,
        initialState = None,
        (graph, window) =>
          discoverHaplotypes(
            graph,
            window,
            args.kmerSize,
            args.minAlignmentQuality,
            args.minAverageBaseQuality)

      )

      readSet.mappedReads.unpersist()

      val filteredGenotypes = GenotypeFilter(genotypes, args).flatMap(AlleleConversions.calledAlleleToADAMGenotype)
      Common.writeVariantsFromArguments(args, filteredGenotypes)
      DelayedMessages.default.print()
    }


  }
}
