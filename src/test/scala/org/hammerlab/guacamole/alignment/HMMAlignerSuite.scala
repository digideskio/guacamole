package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.TestUtil.Implicits._
import org.scalatest.{Matchers, FunSuite}


class HMMAlignerSuite extends FunSuite with Matchers {

  test("score alignment: exact match") {
    val (alignmentStates, alignmentScores) = HMMAligner.scoreAlignmentPaths("TCGA", "TCGA")
    alignmentScores(4, 4).toInt should be(0)
  }


  test ("score alignment: single mismatch") {
    val (alignmentStates, alignmentScores) = HMMAligner.scoreAlignmentPaths("TCGA", "TCCA", mismatchProbability = 1e-2)
    alignmentScores(4, 4).toInt should be(5)
  }

  test ("align exact match") {
    val alignment = HMMAligner.align("TCGA", "TCGA")
    alignment.toCigar should be("4=")
  }

  test ("align: single mismatch") {
    val alignment = HMMAligner.align("TCGA", "TCCA", mismatchProbability = 1e-2)

    alignment.toCigar should be("2=1X1=")
  }

  test ("align long exact match") {
    val sequence = "TCGATGATCTGAGA"
    val alignment = HMMAligner.align(sequence, sequence)
    alignment.toCigar  should be(sequence.length.toString + "=")
  }

  test ("short align with insertion") {
    val alignment = HMMAligner.align("TCCGA", "TCGA")
    alignment.toCigar should be("2=1I2=")
  }

  test ("long align with insertion") {
    val alignment = HMMAligner.align("TCGACCCTCTGA", "TCGATCTGA")
    alignment.toCigar should be("4=3I5=")
  }

  test ("long align with deletion") {
    val alignment = HMMAligner.align("TCGATCTGA", "TCGACCCTCTGA")
    alignment.toCigar should be("4=3D5=")
  }

  test ("mixed mismatch and insertion") {
    val alignment = HMMAligner.align("TCGACCCTCTTA", "TCGATCTGA")
    alignment.toCigar should be("4=3I3=1X1=")
  }
}
