package org.hammerlab.guacamole

import org.hammerlab.guacamole.reads.MappedRead

class AlignmentPairList(
    readName: String,
    firstReadList: Iterable[MappedRead],
    secondReadList: Iterable[MappedRead],
    filter: (MappedRead, MappedRead) => Boolean = ((m1, m2) => true)) {

  def sumMismatchScores(mappedReads: Iterable[MappedRead]): Double = {
    mappedReads.foldLeft(0.0)((acc, mr) => acc + mr.mismatchScore)
  }

  // println("Creating " + readName)
  // firstReadList.foreach(println)
  // secondReadList.foreach(println)

  lazy val sumMismatchScores1: Double = sumMismatchScores(firstReadList)
  lazy val sumMismatchScores2: Double = sumMismatchScores(secondReadList)

  //note that this is not exectued on object creation, only when
  //lociAlignmentPairs are
  def alignmentPairs(): Iterable[AlignmentPair] = {
    for {
      read1 <- firstReadList
      read2 <- secondReadList
      if (filter(read1, read2))
    } yield {
      if (read1.start < read2.start) new AlignmentPair(read1, read2, sumMismatchScores1, sumMismatchScores2)
      else new AlignmentPair(read2, read1, sumMismatchScores2, sumMismatchScores1)
    }
  }

  def lociAlignmentPairs(): Iterable[(GenomicLocation, AlignmentPair)] =
    alignmentPairs
      .flatMap(alignmentPair => alignmentPair.lociAlignmentPairs)

  // def lociAlignmentPairs(): Iterable[(GenomicLocationWithQuality, AlignmentPair)] =
  //   alignmentPairs
  //     .flatMap(alignmentPair => alignmentPair.lociAlignmentPairsWithQuality)

}