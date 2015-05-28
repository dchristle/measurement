from scipy import *
from pylab import *
import sineExample
import cPickle

observedData = array( [ (-1.6, 5.0) , ( 1.4, -9.0) ] )  # the assumed 2 first data points

def rejectionSample( data, example, prior, numSamples ):
	upperBoundOnLikelihood = example.noiseDistribution.pdf(0.0) ** data.shape[0]
	samples = []
	rejects = 0
	while len(samples) < numSamples:
		# draw a sample from the prior
		thetaSample = prior.rvs()
		# evaluate Likelihood for this sample
		likelihood = example.evaluateLikelihood(data[:,0],data[:,1], thetaSample)
		u = random.random()
		assert(likelihood < upperBoundOnLikelihood)
		if u * upperBoundOnLikelihood < likelihood:
			# accept the sample
			samples.append(thetaSample)
			print "Got sample ", len(samples)
		else:
			rejects += 1
	acceptanceRatio = numSamples / float(numSamples + rejects) 
	return (samples, acceptanceRatio)

def main():
    example = sineExample.SineExample()
    prior = example.parameterPrior
    # Parameters for computing the posterior (generate samples from it)
    numberSamples = 300 # the number of samples from the posterior to gather

    rejSamples, aratio = rejectionSample(observedData, example, prior, numberSamples)
	
    print "Acceptance Ratio = 1.0 / " , 1.0/aratio 
    
    # Save the samples to disk
    samplesFile = open("rejectionSamples.pydata",'w')
    cPickle.dump(rejSamples, samplesFile)
    samplesFile.close()

if __name__ == '__main__':
    main()

