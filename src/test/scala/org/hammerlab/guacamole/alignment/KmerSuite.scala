package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases
import org.hammerlab.guacamole.alignment.Kmer
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

}
