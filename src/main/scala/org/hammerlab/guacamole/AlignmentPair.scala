package org.hammerlab.guacamole

import org.hammerlab.guacamole.reads.MappedRead

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
    loci.map(locus => (locus, this))
  }
}