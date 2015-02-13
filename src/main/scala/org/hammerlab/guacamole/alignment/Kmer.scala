package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases

import scala.collection.mutable.ArrayBuffer

case class Kmer(bases: Seq[Byte]) {
  def prefix: Seq[Byte] = bases.slice(0, bases.length - 1)
  def prefixKmer = Kmer(prefix)

  def tailKmer = Kmer(bases.tail)
  def tail = Kmer(bases.tail)

  def lastBase = bases.last

  def possibleNext: Seq[Kmer] = {
    Seq(
      Kmer(bases.tail :+ Bases.A),
      Kmer(bases.tail :+ Bases.T),
      Kmer(bases.tail :+ Bases.C),
      Kmer(bases.tail :+ Bases.G)
    )
  }


  def possiblePrevious: Seq[Kmer] = {
    Seq(
      Kmer(Seq(Bases.A) ++ prefix),
      Kmer(Seq(Bases.T) ++ prefix),
      Kmer(Seq(Bases.C) ++ prefix),
      Kmer(Seq(Bases.G) ++ prefix)
    )
  }

  def :+(base: Byte): Kmer = {
      Kmer(bases :+ base)
  }

  override def toString: String = Bases.basesToString(bases)
}

object Kmer {
  def apply(seq: String): Kmer = {
    Kmer(Bases.stringToBases(seq).toArray)
  }

  def buildSequence(kmers: Seq[Kmer]): Array[Byte] = {
    val seq = new ArrayBuffer[Byte]
    kmers.headOption.map(seq ++= _.prefix)
    kmers.foreach(seq += _.lastBase)
    seq.toArray
  }
}
