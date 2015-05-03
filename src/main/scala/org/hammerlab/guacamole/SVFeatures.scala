package org.hammerlab.guacamole

import scala.math

object SVFeatures {

  def apply(
    locus: GenomicLocation,
    // locus: Long,
    alignmentPairs: Iterable[AlignmentPair],
    insertSizeMean: Long = 300,
    insertSizeSD: Long = 30,
    maxIterations: Int = 100,
    convergenceThreshold: Double = 0.0001,
    maxLogMapQDiff: Double = 6,
    kthNearestNeighbor: Int = 2,
    cutoffInSDs: Int = 5): Option[SVFeatures] = {

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

            logsumexp(weightedPointLikelihoods) // The llh of the point in the GMM
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
      muPrime(1) = updateMuForComponent(gArray, yArray, nArray, 1)

      assert(wPrime.length > 0)
      assert(muPrime.length > 0)

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

    lazy val bestMappingQuality: Double =
      alignmentPairs.map(ap => ap.probabilityMappingIsCorrect).max

    def adaptiveQualityFilter(ap: AlignmentPair): Boolean = {
      (bestMappingQuality - ap.probabilityMappingIsCorrect) <= maxLogMapQDiff
    }

    //naive kNN implementation. Could be optimized, but assuming that 
    // EM is main driver of computational compexity
    // Compute distances to all other alignment pairs,
    // Sort distances into array
    // Get distance of kth nearest neighbor 
    // (keep in mind that distance to itself will be 0 at index 0)
    // take top K+1 values
    // ensure that there does not exists a distance > cutoffInSDs * insertSizeSD
    def outlierDetectionFilter(ap: AlignmentPair): Boolean = {

      val sortedDistances: Array[Long] =
        alignmentPairs
          .map(ap2 => Math.abs(ap.insertSize - ap2.insertSize))
          .toArray
          .sorted
          .take(kthNearestNeighbor + 1)

      !sortedDistances.exists(distance => distance > cutoffInSDs * insertSizeSD)
    }

    val postFilter1: Iterable[AlignmentPair] =
      alignmentPairs
        .filter(adaptiveQualityFilter)

    // Common.progress("Post Filter 1 Size = " + postFilter1.size)

    val postFilter2: Iterable[AlignmentPair] =
      postFilter1
        .filter(outlierDetectionFilter)

    // Common.progress("Post Filter 2 Size = " + postFilter2.size)

    val yArray: Array[Double] =
      postFilter2
        .toArray
        .map(alignmentPair => alignmentPair.insertSize.toDouble)

    //can we get ok data with only one read at a particular locus?
    if (yArray.length >= 1) {
      var i = 0; // The amount of iterations
      var w: Array[Double] = Array(math.log(.5), math.log(.5)) // Array with the logs of the weights of both Gaussians
      var mu: Array[Double] = Array(insertSizeMean, (yArray.sum / yArray.length)) // The means of both Gaussians
      val sigma = insertSizeSD // The stdev of both Gaussians (should be equal)
      var lprime: Double = likelihood(yArray, w, mu, sigma) // LLH of GMM fit
      var l: Double = Double.PositiveInfinity // The previous lprime

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
      // LLH of N(mu, sigma) fit
      val nodelOneComponentLikelihood = likelihood(yArray, Array[Double](math.log(1)), Array[Double](insertSizeMean), sigma)
      // println(l)
      // println(nodelOneComponentLikelihood)

      Some(new SVFeatures(math.exp(w(0)), mu(1), l - nodelOneComponentLikelihood))
    } else None

  }

  def avg(sv1: SVFeatures, sv2: SVFeatures): SVFeatures = {
    new SVFeatures((sv1.w0 + sv2.w0) / 2, (sv1.mu2 + sv2.mu2) / 2, (sv1.lrGMMFit + sv2.lrGMMFit) / 2)
  }
}

class SVFeatures(
    val w0: Double,
    val mu2: Double,
    val lrGMMFit: Double) {

  override def toString(): String =
    "SVFeatures(%f, %f, %f)".format(
      w0,
      mu2,
      lrGMMFit)

}