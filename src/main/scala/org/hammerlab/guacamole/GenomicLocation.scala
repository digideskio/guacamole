package org.hammerlab.guacamole

class GenomicLocation(
    val chromosome: String,
    val position: Long) {

  override def equals(other: Any) = other match {
    case that: GenomicLocation => {
      (this.chromosome == that.chromosome) &&
        (this.position == that.position)
    }
    case _ => false
  }

  override def hashCode = 41 * (position.hashCode * (41 * (chromosome.hashCode + 41)))

  override def toString(): String =
    "GenomicLocation(%s, %d)".format(
      chromosome,
      position)

}

object GenomicLocation {
  def apply(
    chromosome: String,
    position: Long): GenomicLocation = {
    new GenomicLocation(chromosome, position)
  }
}