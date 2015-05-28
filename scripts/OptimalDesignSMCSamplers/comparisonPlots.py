#!/usr/bin/env python
# encoding: utf-8
"""
comparisonPlots.py

Created by hendrik on 2006-05-26.
Copyright (c) 2006 __MyCompanyName__. All rights reserved.
"""
from __future__ import division
from pylab import *
from scipy import *
import cPickle
import math
from plotSMCsamples import *
from GridUtilityPdf import GridUtilityPdf
import sys
import os
import dateutil
import matplotlib.pyplot as plot
import matplotlib.figure as figure

runsDir = "./runs/"

def main():
	evalExponent = 50
	gridPdf = GridUtilityPdf()
	alltheData = loadRunResults("20parts")
    	comparativeAvgLikelihoodValues = getSpecificAvgSampleLikelihoodValues(alltheData, evalExponent, gridPdf)
	for i in xrange(len(alltheData)):
		print alltheData[i][1], " at n_t = ", evalExponent, ": ", comparativeAvgLikelihoodValues[i]
	comparitivePlots(alltheData, plotESSLog)
	comparitivePlots(alltheData, pointPlots)
	comparitivePlots(alltheData, plotIncrementalWeightVariances)
	compareHistograms(alltheData, evalExponent, gridPdf)
	show()


indices = xrange(sys.maxint)
def Indexed(sequence):
	return zip(sequence, indices)


def loadRunResults(nameComponent):
    dataSets = []
    for runFile in os.listdir(runsDir):
        if not runFile.find(nameComponent) == -1 :
            fullPath = os.path.join(runsDir, runFile)
            print "Loading ", fullPath
            data = loadData(fullPath)
            dataSets.append([fullPath, data['samplerName'] + " - " + data['scheduleName'], data])
    return dataSets

def getSpecificAvgSampleLikelihoodValues(alltheData, desired_n_t, gridPdf):
	avgSampleLikelihoods = []
	for d in alltheData:
		dataSet = d[2]
		index = argmax( array(dataSet['n_tStore']) == desired_n_t)
		particles = dataSet['sampleStore'][index]
		value = computeAvgSampleLikelihood(particles, desired_n_t, gridPdf)
		avgSampleLikelihoods.append(value)
	return avgSampleLikelihoods


def plotSpecificIntegerHistogram(dataSet, desired_n_t):
	index = argmax( array(dataSet['n_tStore']) == desired_n_t)
	print "the index at which n_t = " ,desired_n_t, " is " , index
	particles = dataSet['sampleStore'][index]
	plotHist(particles,100,-2.0, 2.0)

def compareHistograms(alltheData, desired_n_t, gridPdf):
	numPlots = len(alltheData) + 1
	figure()
	grid = linspace(-1.99, 1.99, 2000)
	gridEvals = array([gridPdf.evaluate(v,desired_n_t) for v in grid])
	subplot(numPlots, 1, 1)
	plot(grid, gridEvals)
	title('Grid evaluations')
	for d, i in Indexed(alltheData):
		subplot(numPlots, 1, i+2)
		plotSpecificIntegerHistogram(d[2], desired_n_t)
		title(d[1])


def comparitivePlots(alltheData, plotFct):
	numPlots = len(alltheData)
	figure()
	for d, i in Indexed(alltheData):
		subplot(numPlots, 1, i+1)
		plotFct(d[2])
		title(d[1])


if __name__ == '__main__':
    main()

