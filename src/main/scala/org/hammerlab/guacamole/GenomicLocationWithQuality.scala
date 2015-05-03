package org.hammerlab.guacamole

class GenomicLocationWithQuality(
    val chromosome: Int,
    val position: Long,
    val pMappingCorrect: Double = 1.0) {

  override def equals(other: Any) = other match {
    case that: GenomicLocationWithQuality => {
      (this.chromosome == that.chromosome) &&
        (this.position == that.position) &&
        (this.pMappingCorrect == that.pMappingCorrect)
    }
    case _ => false
  }

  override def hashCode = 41 * (pMappingCorrect.hashCode * (position.hashCode * (41 * (chromosome.hashCode + 41))))

  override def toString(): String =
    "GenomicLocationWithQuality(%d, %d, %f)".format(
      chromosome,
      position,
      pMappingCorrect)

  def toGenomicLocation(): GenomicLocation = new GenomicLocation(chromosome, position)

}

object GenomicLocationWithQuality {
  def apply(
    chromosome: Int,
    position: Long,
    pMappingCorrect: Double): GenomicLocationWithQuality = {
    new GenomicLocationWithQuality(chromosome, position, pMappingCorrect)
  }

  def apply(
    location: GenomicLocation,
    pMappingCorrect: Double): GenomicLocationWithQuality = {
    new GenomicLocationWithQuality(location.chromosome, location.position, pMappingCorrect)
  }
}