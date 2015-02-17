package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases
import scala.collection.mutable

class DeBrujinGraph(kmerSize: Int,
                    val kmerCounts: mutable.Map[Seq[Byte], Int],
                    val nodes: mutable.Map[Seq[Byte], Kmer]) {

  type Path = List[Kmer]
  type PathScore = Int
  /**
   * Remove kmers that are not in at least minSupport reads
   * @param minSupport minimum of reads a kmer should appear in
   */
  private def pruneKmers(minSupport: Int) = {
    kmerCounts
      .filter(_._2 < minSupport)
      .foreach( {case (kmer, count) => kmerCounts.remove(kmer) })
  }

//  def mergeNodes(): Unit = {
//    kmerCounts.keys.foreach( kmer => {
//      var mergedKmer = kmer
//      var childNodes = children(kmer)
//      while (childNodes.length == 1) {
//        val child = childNodes(0)
//        val totalCount = kmerCounts(kmer) + kmerCounts(child)
//        kmerCounts.remove(kmer)
//        kmerCounts.remove(child)
//        kmerCounts(Kmer(kmer.))
//
//      }
//    }
//    )
//  }

  def depthFirstSearch(source: Kmer,
                       sink: Kmer,
                       kmerScore: (Kmer => Int) = kmerCounts.getOrElse(_, 0),
                       minPathLength: Int = 1,
                       maxPathLength: Int = Int.MaxValue,
                       maxPaths: Int = 10): List[(Path, PathScore)] = {
    var paths: List[(Path, PathScore)] = List.empty[(Path, PathScore)]
    var frontier: mutable.Stack[Kmer] = mutable.Stack(source)
    var currentPath: Path = List.empty
    var pathScore: Int = 0
    // explore branches until we accumulate the maximum number of appropriate length paths
    while (frontier.nonEmpty && paths.size < maxPaths) {
      val next = frontier.pop()
      println(next)
      pathScore += kmerScore(next)
      currentPath = next :: currentPath
      if (next != sink && currentPath.size < maxPathLength) {
        // Keep searching down tree
        val nextChildren = children(next)
        if (nextChildren.nonEmpty && currentPath.size < maxPathLength) {
          frontier ++= nextChildren
        } else { //found sink or too long
          if (next == sink && currentPath.size + 1 >= minPathLength) {
            // Found legitimate path to sink, save path
            paths = (currentPath, pathScore) :: paths
          } // else degenerate path, too long or too short
          currentPath = List.empty
          pathScore = 0
          println("+++")

        }
      }
    }
    paths
  }

  def roots: Iterable[Kmer] = {
     kmerCounts.keys.filter(parents(_).size == 0)
  }

  def children(node: Kmer): Seq[Kmer] = {
    node.possibleNext.filter(kmerCounts.contains).toSeq
  }

  def parents(node: Kmer): Seq[Kmer] = {
    node.possiblePrevious.filter(kmerCounts.contains).toSeq
  }

}

object DeBrujinGraph {
  type Sequence = Seq[Byte]
  def apply(sequences: Seq[Sequence],
            kmerSize: Int,
            minOccurrence: Int = 1): DeBrujinGraph = {
    val kmerCounts: mutable.Map[Kmer, Int] = mutable.Map.empty[Kmer, Int]

    sequences.filter(Bases.allStandardBases(_))
      .foreach(
        _.sliding(kmerSize)
          .foreach(seq => {
          val kmer = Kmer(seq, kmerSize)
          val count = kmerCounts.getOrElse(kmer, 0)
          kmerCounts.update(kmer, count + 1)
        })
      )

    val graph = new DeBrujinGraph(kmerCounts, kmerSize)
    //TODO(arahuja) Only add in kmers once they minOccurrence rather than post-pruning
    graph.pruneKmers(minOccurrence)
    //graph.mergeNodes()

    graph
  }

  type Sequences = Seq[String]

  def fromString(sequences: Sequences,
            kmerSize: Int,
            minOccurrence: Int = 1,
            isString: Boolean = true): DeBrujinGraph = {
    DeBrujinGraph(sequences.view.map(Bases.stringToBases(_)), kmerSize, minOccurrence)
  }

}
