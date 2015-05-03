package org.hammerlab.guacamole

import com.esotericsoftware.kryo.io.{ Input, Output }
import com.esotericsoftware.kryo.{ Kryo, Serializer }

class GenomicLocation(
    val chromosome: Int,
    val position: Long) extends Serializable {

  override def equals(other: Any) = other match {
    case that: GenomicLocation => {
      (this.chromosome == that.chromosome) &&
        (this.position == that.position)
    }
    case _ => false
  }

  override def hashCode = 41 * (position.hashCode * (41 * (chromosome.hashCode + 41)))

  override def toString(): String =
    "GenomicLocation(%d, %d)".format(
      chromosome,
      position)

}

class GenomicLocationSerializer
    extends Serializer[GenomicLocation] {
  def write(kryo: Kryo, output: Output, obj: GenomicLocation) = {
    // println("Serializing Genomic Location")
    output.writeInt(obj.chromosome)
    output.writeLong(obj.position)
  }

  def read(kryo: Kryo, input: Input, klass: Class[GenomicLocation]): GenomicLocation = {
    // println("Deserializing Genomic Location")
    val c = input.readInt
    val p = input.readLong
    new GenomicLocation(c, p)
  }
}

// object GenomicLocation {
//   def apply(
//     chromosome: Int,
//     position: Long): GenomicLocation = {
//     new GenomicLocation(chromosome, position)
//   }
// }