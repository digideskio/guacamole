package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.alignment.AlignmentState.AlignmentState

private[alignment] object AlignmentState extends Enumeration {
  type AlignmentState = Value
  val Match, Mismatch, Insertion, Deletion = Value
}

case class ReadAlignment(alignments: Seq[AlignmentState], alignmentScore: Int) {
  def cigarKey(alignmentOperator: AlignmentState): String = {
    alignmentOperator match {
      case AlignmentState.Match | AlignmentState.Mismatch => "M"
      case AlignmentState.Insertion => "I"
      case AlignmentState.Deletion => "D"
    }
  }

  def toCigar: String = {

    def runLengthEncode(operators: Seq[String]): String = {
      var lastOperator = operators.head
      var i = 1
      val rle = new StringBuffer()
      var currentRun = 1
      while (i < operators.size) {
        if (operators(i) == lastOperator) {
          currentRun += 1
        } else {
          rle.append(currentRun.toString + lastOperator)
          currentRun = 1
        }
        lastOperator = operators(i)
        i += 1
      }
      rle.append(currentRun.toString + lastOperator)
      rle.toString
    }

    runLengthEncode(alignments.map(cigarKey(_)))
  }
}