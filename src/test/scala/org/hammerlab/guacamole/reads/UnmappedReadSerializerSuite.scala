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

package org.hammerlab.guacamole.reads

import org.hammerlab.guacamole.TestUtil
import org.hammerlab.guacamole.TestUtil.Implicits._
import org.scalatest.Matchers

class UnmappedReadSerializerSuite extends TestUtil.SparkFunSuite with Matchers {

  test("serialize / deserialize unmapped read") {
    val read = UnmappedRead(
      22, // token
      "TCGACCCTCGA",
      Array[Byte]((10 to 20).map(_.toByte): _*),
      true,
      "some sample name",
      false,
      isPositiveStrand = true,
      matePropertiesOpt = Some(
        MateProperties(
          isFirstInPair = true,
          inferredInsertSize = Some(300),
          isMateMapped = true,
          Some("chr5"),
          Some(100L),
          false
        )
      )
    )

    val serialized = TestUtil.serialize(read)
    val deserialized = TestUtil.deserialize[UnmappedRead](serialized)

    // We *should* be able to just use MappedRead's equality implementation, since Scala should implement the equals
    // method for case classes. Somehow, something goes wrong though, and this fails:

    // deserialized should equal(read)

    // So, instead, we'll compare each field ourselves:
    deserialized.token should equal(read.token)
    deserialized.sequence should equal(read.sequence)
    deserialized.baseQualities should equal(read.baseQualities)
    deserialized.isDuplicate should equal(read.isDuplicate)
    deserialized.sampleName should equal(read.sampleName)
    deserialized.failedVendorQualityChecks should equal(read.failedVendorQualityChecks)
    deserialized.isPositiveStrand should equal(read.isPositiveStrand)
    deserialized.matePropertiesOpt should equal(read.matePropertiesOpt)
  }

}
