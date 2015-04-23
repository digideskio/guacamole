package org.hammerlab.guacamole

import scala.math

object SVFeatures {

  def apply(
    locus: Long,
    alignmentPairs: Iterable[AlignmentPair],
    insertSizeMean: Long = 300,
    insertSizeSD: Long = 30,
    maxIterations: Int = 100,
    convergenceThreshold: Double = 0.0001): SVFeatures = {

    def lognormal(y: Double, mu: Double, sigma: Double): Double = {
      -1.0 * math.log(sigma + math.sqrt(2.0 * math.Pi)) + -1.0 / 2.0 * math.pow(((y - mu) / sigma), 2)
    }

    def logsumexp(xArray: Array[Double]): Double = {
      assert(xArray.length > 0)
      val m = xArray.max
      val sumexp =
        xArray
          .map(x => math.exp(x - m))
          .sum
      m + math.log(sumexp)
    }

    def likelihood(
      yArray: Array[Double],
      wArray: Array[Double],
      muArray: Array[Double],
      sigma: Double): Double = {

      assert(yArray.length > 0)
      assert(wArray.length > 0)
      assert(muArray.length > 0)
      val sumWeightedLikelihoods: Double =
        yArray
          .map(y => {
            val weightedPointLikelihoods =
              wArray
                .zip(muArray)
                //w + lognormal(y, mu, sigma)
                .map(pair => pair._1 + lognormal(y, pair._2, sigma))

            logsumexp(weightedPointLikelihoods)
          })
          .sum

      sumWeightedLikelihoods / yArray.length
    }

    def gamma(
      yArray: Array[Double],
      wArray: Array[Double],
      muArray: Array[Double],
      sigma: Double): Array[Array[Double]] =
      yArray
        .map(y => {
          val likelihoodRow =
            wArray
              .zip(muArray)
              .map(pair => pair._1 + lognormal(y, pair._2, sigma))

          val total = logsumexp(likelihoodRow)

          likelihoodRow.map(l => l - total)
        })

    def updateW(
      yArray: Array[Double],
      nArray: Array[Double]): Array[Double] =
      nArray.map(n => n - math.log(yArray.length))

    def updateMuForComponent(
      gArray: Array[Array[Double]],
      yArray: Array[Double],
      nArray: Array[Double],
      component: Int): Double = {

      val lse =
        gArray
          .zip(yArray)
          .map(pair => pair._1(component) + math.log(pair._2))

      val numerator = logsumexp(lse)
      math.exp(numerator - nArray(component))
    }

    def calculateN(gArray: Array[Array[Double]]): Array[Double] = {
      //transpose
      assert(gArray.length > 0)
      val gSeq: Seq[Seq[Double]] = gArray.map(a => a.toSeq).toSeq
      val empty: Seq[Seq[Double]] = Seq()

      gSeq
        .foldLeft(empty)((acc, row) => {

          if (acc.length == 0) row.map(v => Seq(v))
          else {
            acc.zip(row).map(pair => {
              pair._1 :+ pair._2
            })
          }

        })
        .map(a => logsumexp(a.toArray)).toArray
    }

    def emStep(
      yArray: Array[Double],
      wArray: Array[Double],
      muArray: Array[Double],
      sigma: Double): (Array[Double], Array[Double]) = {

      assert(yArray.length > 0)
      assert(wArray.length > 0)
      assert(muArray.length > 0)

      val gArray = gamma(yArray, wArray, muArray, sigma)

      assert(gArray.length > 0)

      val nArray = calculateN(gArray)

      assert(nArray.length > 0)

      val wPrime = nArray.map(n => n - math.log(yArray.length))
      val muPrime = muArray.clone

      assert(wPrime.length > 0)
      assert(muPrime.length > 0)

      // val mu1Prime = updateMuForComponent(gArray, yArray, nArray, 1)
      // muPrime(1) = mu1Prime

      (wPrime, muPrime)
    }

    // Default mean fragment size and standard deviation of library are
    // INSERT_SIZE=300
    // INSERT_SIZE_SD=30
    // These should be command line parameters

    //Let Y = y1,y2...ym be the observed insert sizes
    //Library has mean fragement size of mu with standard deviation sigma
    //we initialize the two components to have means mu and yBar with standard devs of sigma
    //set alpha to .5

    //In each E step

    //fix this later 
    //filters include nearest neighbors and adaptave quality
    // val filteredAlignmentPairs: Array[AlignmentPair] = alignmentPairs.toArray

    def alignmentPairFilter(ap: AlignmentPair): Boolean = {
      true
    }

    val yArray: Array[Double] =
      alignmentPairs
        .filter(alignmentPairFilter)
        .toArray
        .map(alignmentPair => alignmentPair.insertSize.toDouble)

    var i = 0;
    var w: Array[Double] = Array(math.log(.5), math.log(.5))
    var mu: Array[Double] = Array(insertSizeMean, (yArray.sum / yArray.length))
    val sigma = insertSizeSD
    var lprime: Double = likelihood(yArray, w, mu, sigma)
    var l: Double = Double.PositiveInfinity

    assert(w.length > 0)
    assert(mu.length > 0)
    assert(yArray.length > 0)

    while (math.abs(l - lprime) > convergenceThreshold && i < maxIterations) {

      val threshold: Double = math.abs(l - lprime)
      // Common.progress("Locus: %d, Iteration: %d".format(locus, i))
      // Common.progress("Threshold : %f".format(threshold))
      l = lprime

      val pair = emStep(yArray, w, mu, sigma)

      w = pair._1
      mu = pair._2

      assert(w.length > 0)
      assert(mu.length > 0)
      assert(yArray.length > 0)

      lprime = likelihood(yArray, w, mu, sigma)

      i += 1
    }

    val nodelOneComponentLikelihood = likelihood(yArray, Array[Double](math.log(1)), Array[Double](insertSizeMean), sigma)
    // println(l)
    // println(nodelOneComponentLikelihood)

    new SVFeatures(math.exp(w(0)), mu(1), l - nodelOneComponentLikelihood)
  }
}

class SVFeatures(
    val w0: Double,
    val mu2: Double,
    val lrHeterozygous: Double) {

  override def toString(): String =
    "SVFeatures(%f, %f, %f)".format(
      w0,
      mu2,
      lrHeterozygous)

}