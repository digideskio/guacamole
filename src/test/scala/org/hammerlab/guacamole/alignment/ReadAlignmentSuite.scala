package org.hammerlab.guacamole.alignment

import org.scalatest.{FunSuite, Matchers}


class ReadAlignmentSuite  extends FunSuite with Matchers {

  test ("test cigar string: all match") {
    val alignment = ReadAlignment(
      Seq(
        AlignmentState.Match,
        AlignmentState.Match,
        AlignmentState.Match,
        AlignmentState.Match,
        AlignmentState.Match,
        AlignmentState.Match), 60)

    alignment.toCigar should be("6M")
  }

  test ("test cigar string: mixed match/insertion") {
    val alignment = ReadAlignment(
      Seq(AlignmentState.Match,
        AlignmentState.Match,
        AlignmentState.Match,
        AlignmentState.Insertion,
        AlignmentState.Insertion,
        AlignmentState.Match), 60)

    alignment.toCigar should be("3M2I1M")
  }


  test ("test cigar string: start with single match") {
    val alignment = ReadAlignment(
      Seq(
        AlignmentState.Match,
        AlignmentState.Insertion,
        AlignmentState.Insertion,
        AlignmentState.Insertion,
        AlignmentState.Insertion,
        AlignmentState.Match), 60)

    alignment.toCigar should be("1M4I1M")
  }
}
