package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases
import org.scalatest.{FunSuite, Matchers}

class KmerSuite extends FunSuite with Matchers {


  test("kmer last base") {
    val tcga = Kmer("TCGA")
    tcga.lastBase should be(Bases.A)

    val acgt = Kmer("ACGT")
    acgt.lastBase should be(Bases.T)
  }

  test("kmer tail") {
    val tcga = Kmer("TCGA")
    tcga.tailKmer should be(Kmer("CGA"))

    val agct = Kmer("AGCT")
    agct.tailKmer should be(Kmer("GCT"))
  }

  test("kmer prefix") {
    val tcga = Kmer("TCGA")
    tcga.prefixKmer should be(Kmer("TCG"))

    val acgt = Kmer("ACGT")
    acgt.prefixKmer should be(Kmer("ACG"))
  }

  test("build sequence from set of kmers") {
    val orignalSequence = "TTTTATATGGGGCGCG"

    val fourMers = orignalSequence.sliding(4).map(Kmer(_)).toSeq
    Bases.basesToString(Kmer.buildSequence(fourMers)) should be(orignalSequence)

    val twoMers = orignalSequence.sliding(2).map(Kmer(_)).toSeq
    Bases.basesToString(Kmer.buildSequence(twoMers)) should be(orignalSequence)
  }

}
