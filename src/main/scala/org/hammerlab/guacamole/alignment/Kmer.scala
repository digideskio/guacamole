package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases

import scala.collection.mutable.ArrayBuffer

case class Kmer(bases: Seq[Byte], kmerSize: Int) {
  def prefix: Seq[Byte] = bases.slice(0, kmerSize)
  def prefixKmer = Kmer(prefix, kmerSize)

  def tail = Kmer(bases.slice(bases.length - kmerSize, bases.length + 1), kmerSize)

  def lastBase = bases.last

  def possibleNext: Seq[Kmer] = {
    Seq(
      Kmer(bases.tail :+ Bases.A, kmerSize),
      Kmer(bases.tail :+ Bases.T, kmerSize),
      Kmer(bases.tail :+ Bases.C, kmerSize),
      Kmer(bases.tail :+ Bases.G, kmerSize)
    )
  }


  def possiblePrevious: Seq[Kmer] = {
    Seq(
      Kmer(Seq(Bases.A) ++ prefix, kmerSize),
      Kmer(Seq(Bases.T) ++ prefix, kmerSize),
      Kmer(Seq(Bases.C) ++ prefix, kmerSize),
      Kmer(Seq(Bases.G) ++ prefix, kmerSize)
    )
  }

  def :+(base: Byte): Kmer = {
      Kmer(bases :+ base, kmerSize)
  }

  override def toString: String = Bases.basesToString(bases)
}

object Kmer {
  def apply(seq: String): Kmer = {
    Kmer(Bases.stringToBases(seq).toArray, kmerSize = seq.length)
  }

  def buildSequence(kmers: Seq[Kmer]): Array[Byte] = {
    val seq = new ArrayBuffer[Byte]
    kmers.headOption.map(seq ++= _.prefix)
    kmers.foreach(seq += _.lastBase)
    seq.toArray
  }
}
