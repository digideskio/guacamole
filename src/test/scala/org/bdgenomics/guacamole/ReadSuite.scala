/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.guacamole

import org.scalatest.matchers.ShouldMatchers
import com.google.common.collect._
import com.esotericsoftware.kryo.Kryo
import com.twitter.chill.{ IKryoRegistrar, KryoInstantiator, KryoPool }
import net.sf.samtools.{ Cigar, TextCigarCodec }
import org.bdgenomics.adam.util.MdTag

class ReadSuite extends TestUtil.SparkFunSuite with ShouldMatchers {
  test("serialize / deserialize mapped read") {
    val read = MappedRead(
      Bases.stringToBases("TCGACCCTCGA"),
      Array[Byte]((10 to 20).map(_.toByte): _*),
      true,
      "some sample name",
      "some reference contig",
      50,
      325352323,
      TextCigarCodec.getSingleton.decode(""),
      None, // mdtag
      false)

    val serialized = TestUtil.serialize(read)
    val deserialized = TestUtil.deserialize[MappedRead](serialized)

    // We *should* be able to just use MappedRead's equality implementation, since Scala should implement the equals
    // method for case classes. Somehow, something goes wrong though, and this fails:

    // deserialized should equal(read)

    // So, instead, we'll compare each field ourselves:
    deserialized.sequence should equal(read.sequence)
    deserialized.baseQualities should equal(read.baseQualities)
    deserialized.isDuplicate should equal(read.isDuplicate)
    deserialized.sampleName should equal(read.sampleName)
    deserialized.referenceContig should equal(read.referenceContig)
    deserialized.alignmentQuality should equal(read.alignmentQuality)
    deserialized.start should equal(read.start)
    deserialized.cigar should equal(read.cigar)
    deserialized.mdTag should equal(read.mdTag)
    deserialized.failedVendorQualityChecks should equal(read.failedVendorQualityChecks)
  }

  test("serialize / deserialize unmapped read") {
    val read = UnmappedRead(
      Bases.stringToBases("TCGACCCTCGA"),
      Array[Byte]((10 to 20).map(_.toByte): _*),
      true,
      "some sample name",
      false)

    val serialized = TestUtil.serialize(read)
    val deserialized = TestUtil.deserialize[UnmappedRead](serialized)

    // We *should* be able to just use MappedRead's equality implementation, since Scala should implement the equals
    // method for case classes. Somehow, something goes wrong though, and this fails:

    // deserialized should equal(read)

    // So, instead, we'll compare each field ourselves:
    deserialized.sequence should equal(read.sequence)
    deserialized.baseQualities should equal(read.baseQualities)
    deserialized.isDuplicate should equal(read.isDuplicate)
    deserialized.sampleName should equal(read.sampleName)
    deserialized.failedVendorQualityChecks should equal(read.failedVendorQualityChecks)
  }
}