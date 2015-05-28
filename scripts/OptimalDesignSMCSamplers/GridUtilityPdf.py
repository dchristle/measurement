#!/usr/bin/env python
from __future__ import division
from scipy import *
from pylab import *
import sineExample
from StaticDiscreteDistribution import StaticDiscreteDistribution
import cPickle

class GridUtilityPdf:

	def __init__(self):
		samplesFile = open("gridUtilities.pydata",'r')
		(self.drange, self.expUtilities) = cPickle.load(samplesFile)
		self.expUtilities = array(self.expUtilities)
		samplesFile.close()
		self.normalizingConstants = {}

	def computeNormalizingConstant(self, exponent = 1):
		return 4.0 * average(self.expUtilities ** exponent)

	def evaluate(self, design, exponent = 1):
		assert isinstance(exponent, int)
		if self.normalizingConstants.has_key(exponent):
			normConst = self.normalizingConstants[exponent]
		else:
			normConst = self.computeNormalizingConstant(exponent)
			self.normalizingConstants[exponent] = normConst
		# linear interpolation
		# search for the first index at which expUtilities is larger or equal to the provided design value
		minLargerEqIndex = int(self.drange.searchsorted(design))
		maxSmallerIndex = minLargerEqIndex - 1
		assert maxSmallerIndex >= 0
		maxSmallerD = self.drange[maxSmallerIndex]
		minLargerEqD = self.drange[minLargerEqIndex]
		assert design > maxSmallerD
		assert design <= minLargerEqD
		interpolationFactor = (design - maxSmallerD) / (minLargerEqD - maxSmallerD)
		# we are linearly interpolating between U(d)^exponent as computed at the grid location.
		# this is in contrast to linearly interpolating U(d) and taking the result to a power of exponent
		maxSmallerDExpectedUtility = self.expUtilities[maxSmallerIndex] ** exponent
		minLargerEqDExpectedUtility = self.expUtilities[minLargerEqIndex] ** exponent
		interpolatedExpectedUtility = maxSmallerDExpectedUtility + interpolationFactor * (minLargerEqDExpectedUtility - maxSmallerDExpectedUtility)
		return interpolatedExpectedUtility / normConst


def main():
	testPdf = GridUtilityPdf()
	exponent = 50
	numSamples = 10000
	samples = []
	sampledIntegral = 0.0
	uniformSampler = stats.uniform(-2, 4)
	for i in xrange(numSamples):
		d = uniformSampler.rvs()
		U = testPdf.evaluate(d, exponent)
		samples.append( (d,U) )
		sampledIntegral += 4.0 * U / numSamples
	print "approximate integral evaluation: " , sampledIntegral
	samplesArray = array(samples)
	figure()
	plot(samplesArray[:,0], samplesArray[:,1],'.')
	show()
	

if __name__ == '__main__':
	main()