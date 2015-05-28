from __future__ import division
from pylab import *
from scipy import *
import cPickle
import math
from GridUtilityPdf import GridUtilityPdf
import dateutil
import matplotlib.pyplot as plot
import matplotlib.figure as figure

def plotHist( samples, bins, min, max):
	"Plot a histogram with at least the range from min to max"
	# HACK: we add an artifical sample point at both min and max
	plotsamples = r_[ samples, array([min,max])]
	hist(plotsamples, bins)

def plotHistograms(data, requestedNumHistograms = 10):
	numStoredSamples = data['numStoredSamples']
	numParticles = data['numParticles']
	sampleStore = data['sampleStore']
	numHistograms = min(requestedNumHistograms, numStoredSamples)
	stepSize = numStoredSamples // numHistograms
	for i in xrange(0,numHistograms):
		subplot(numHistograms,1,i + 1)
		plotHist(sampleStore[numStoredSamples-(i*stepSize)-1],100,-2.0, 2.0)

def pointPlots(data):
	numStoredSamples = data['numStoredSamples']
	numParticles = data['numParticles']
	sampleStore = data['sampleStore']
	hold (True)
	ylim(0, numStoredSamples)
	for i in xrange(0,numStoredSamples):
		plot(sampleStore[i], ones(numParticles) * i + 0.0001 * rand(numParticles), 'bx')
	hold (False)

def plotIncrementalWeightVariances(data):
	incrementalWeightStore = data['incrementalWeightStore']
	weightVariances = []
	for incrementalWeights in incrementalWeightStore:
		normalizedWeights = incrementalWeights / sum(incrementalWeights)
		weightVariances.append(var(normalizedWeights))
	plot(weightVariances)

def plotTraces(data, requestedNumTraces = 10):
	numStoredSamples = data['numStoredSamples']
	numParticles = data['numParticles']
	sampleStore = data['sampleStore']
	numTraces = min(numParticles, requestedNumTraces)
	tracesStepSize = numParticles // numTraces
	hold(True)
	for i in xrange(0,numParticles,tracesStepSize):
		particleTrace = [sampleStore[j][i] for j in xrange(numStoredSamples)]
		plot(particleTrace, range(numStoredSamples))
	hold(False)

def plotLastHistogram(data, saveFigureFilename = ""):
	numStoredSamples = data['numStoredSamples']
	sampleStore = data['sampleStore']
	if len(saveFigureFilename) > 0:
		figure(figsize=[10,3])
	plotHist(sampleStore[numStoredSamples-1],100,-2.0, 2.0)
	if len(saveFigureFilename) > 0:
		savefig(saveFigureFilename, dpi=150)

def plotWeightedHistogram(locations, weights, numberOfBins, min, max):
	bins = linspace(min,max, numberOfBins)
	stepSize = bins[1] - bins[0]
	binValues = zeros(numberOfBins, float)
	for i in xrange(len(locations)):
		index = bins.searchsorted(locations[i])
		binValues[index] += weights[i]
	bar(bins, binValues)

def computeAvgSampleLikelihood(particles, exponent, gridPdf):
	likelihoods = [gridPdf.evaluate(d, exponent) for d in particles]
	return mean(likelihoods)

def plotAvgSampleLiklihoods(data):
	n_tStore = data['n_tStore']
	sampleStore = data['sampleStore']
	gridPdf = GridUtilityPdf()
	plotPoints =[]
	for i in xrange(data['numStoredSamples']):
		exponent = n_tStore[i]
		if exponent > 0:
			avgSampleLikelihood = computeAvgSampleLikelihood(sampleStore[i], exponent, gridPdf)
			plotPoints.append((exponent, avgSampleLikelihood))
	plotPoints = array(plotPoints)
	bar(plotPoints[:,0],plotPoints[:,1])

def plotESSLog(data):
	essLog = data['essLog']
	essLog = array(essLog)
	essValues = essLog[:,0]
	iterations = essLog[:,1]
	plot(iterations, essValues)


def loadData(filename):
	samplesFile = open(filename,'r')
	dataArray = cPickle.load(samplesFile)
	samplesFile.close()
	# just for debugging
	numStoredSamples= len(dataArray['sampleStore'])
	numParticles = len(dataArray['sampleStore'][0])
	print "loaded",  filename, "with", numStoredSamples, "samples and", numParticles, "particles."
	print "resampled", dataArray['resampleCount'], "times"
	dataArray['numStoredSamples'] = numStoredSamples
	dataArray['numParticles'] = numParticles
	return dataArray

def plotItAll(filename, requestedNumHistograms, requestedNumTraces, saveFigureFilename = ""):
	data = loadData(filename)
	figure.Figure()
	plotHistograms(data, requestedNumHistograms)
	#figure()
	#plotAvgSampleLiklihoods(data)

	figure.Figure()
	pointPlots(data)
	#figure()
	#plotTraces(data, requestedNumTraces)
	#figure()
	#plotLastHistogram(data, saveFigureFilename)
	figure()
	plotIncrementalWeightVariances(data)
	figure()
	plotESSLog(data)
	show()


if __name__ == '__main__':
	filename = "SMCAnnealingAlgorithm1_samples.pydata"
	saveFigFilename = ""
##	if len(sys.argv) > 1:
##		filename = sys.argv[1]
##	if len(sys.argv) > 2:
##		saveFigFilename = sys.argv[2]
	plotItAll(filename, 10, 5, saveFigFilename)

