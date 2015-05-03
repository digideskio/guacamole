package org.hammerlab.guacamole

import org.hammerlab.guacamole.reads.MappedRead

class AlignmentPair(
    val firstRead: MappedRead,
    val secondRead: MappedRead,
    sumMismatchScores1: Double,
    sumMismatchScores2: Double,
    lociInterval: Long = 25) {

  require(firstRead.chr == secondRead.chr, "chr values do not match")
  def insertSize(): Long = secondRead.end - firstRead.start
  def alignmentScore(): Int = firstRead.alignmentQuality + secondRead.alignmentQuality

  val probabilityMappingIsCorrect: Double = {
    Math.log(firstRead.mismatchScore / sumMismatchScores1) + Math.log(secondRead.mismatchScore / sumMismatchScores2)
  }

  def loci(): Iterable[GenomicLocation] = {
    val start = firstRead.start
    val end = secondRead.end
    val regex = """chr(\d+)""".r
    val chr: Int = firstRead.chr match {
      case regex(ch) => ch.toInt
    }

    var lociStart = (start / lociInterval) * lociInterval
    if ((start - lociStart) > lociInterval / 2) lociStart += lociInterval

    var lociEnd = ((end / lociInterval) + 1) * lociInterval
    if ((lociEnd - end) > lociInterval / 2) lociEnd -= lociInterval

    // val lociRange: NumericRange = new NumericRange(lociStart, lociEnd, lociInterval.toLong, true)
    (lociStart to lociEnd by lociInterval.toLong)
      .map(locus => new GenomicLocation(chr, locus))

  }

  // def lociAlignmentPairsWithQuality: Iterable[(GenomicLocationWithQuality, AlignmentPair)] = {
  //   loci.map(locus => (GenomicLocationWithQualicy(locus, probabilityMappingIsCorrect), this))
  // }

  def lociAlignmentPairs: Iterable[(GenomicLocation, AlignmentPair)] = {
    loci.map(locus => (locus, this))
  }

  override def toString(): String =
    "AlignmentPair(%d, %d, %f)".format(
      insertSize,
      alignmentScore,
      probabilityMappingIsCorrect)

}