package org.hammerlab.guacamole

import org.hammerlab.guacamole.reads.MappedRead

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