from muellerSampler import *
import cPickle

numParallelRuns = 100
chainLength = 200
stepSize = 4

runs = []
endpoints = []
random.seed()

for i in xrange(numParallelRuns):
	print 'Mueller MCMC run ', i
	mcmcRun = MuellerSamplerRun(chainLength, 1, stepSize)
	runs.append(mcmcRun)
	endpoints.append(mcmcRun.finalSamples[0])

samplesFile = open('multiMuellerData.pydata','w')
cPickle.dump((endpoints, runs), samplesFile)
samplesFile.close()
