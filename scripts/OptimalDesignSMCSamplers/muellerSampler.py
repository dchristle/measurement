from __future__ import division
from scipy import *
from pylab import *
import sineExample
from StaticDiscreteDistribution import StaticDiscreteDistribution
from Proposals import *
import cPickle
from plotSMCsamples import plotHist


class LinearSchedule:
	def __init__(self, stepsize):
		self.stepsize = stepsize
		
	def __call__(self,	iteration ):
		return (iteration // self.stepsize) + 1

# global variables, bad, but I dont care.
ex = sineExample.SineExamplePosterior()

# IMPORTANT: This is where the proposal is chosen, currently just a local proposal.
# This proposal also is used by the SMC samplers currently
designProposal = DesignLocalProposal()


class MuellerSamplerRun:

	def __init__(self, numBurnInIterations, numSampleIterations, scheduleStepSize):
		
		designPointInit = random.uniform(-2.0,2.0)

		self.numBurnInIterations = numBurnInIterations
		self.numSampleIterations = numSampleIterations

		self.burninSamples = empty(numBurnInIterations, float)
		self.finalSamples = empty(numSampleIterations, float)

		self.schedule = LinearSchedule(scheduleStepSize)
		
		# Initializations
		curDesignPoint = designPointInit
		curUtility = 1e-20
		oldNumExperiments = 1

		accepts = 0
		rejects = 0
		for i in xrange(numBurnInIterations):
			numExperiments = self.schedule(i)
			if numExperiments != oldNumExperiments:
				# whenever we up the dimensionality, we need to compensate by changing the utility of the previous
				# thanks to division from the future no conversion to float necessary here.This was broken before
				curUtility = curUtility ** (numExperiments / oldNumExperiments)
				oldNumExperiments = numExperiments

			curDesignPoint, curUtility, acceptance = self.MetropolisHastingsStep(curDesignPoint, curUtility, numExperiments)
			if acceptance:
				accepts += 1
			else:
				rejects += 1

			self.burninSamples[i] = curDesignPoint

		self.burnInAcceptanceRate = accepts / (accepts + rejects)
		
		accepts = 0
		rejects = 0
		for i in xrange(numSampleIterations):
			curDesignPoint, curUtility, acceptance = self.MetropolisHastingsStep(curDesignPoint, curUtility, numExperiments)
			if acceptance:
				accepts += 1
			else:
				rejects += 1

			self.finalSamples[i] = curDesignPoint

		self.finalAcceptanceRate = accepts / (accepts + rejects)

	def evaluateUtility( self, designPoint, outcomes ):
		if designPoint < -2 or designPoint > 2:
			return 0.0
		# evaluate all the sine functions corresponding to the parameter samples from the prior
		# at the design point
		functionValues = array([ f(designPoint) for f in ex.postSampleFunctions])

		# now evaluate all the marginal likelihoods p(y | d) for all experiment outcomes y
		marginalLikelihoods = [ average( ex.noiseDistribution.pdf( functionValues - y ) ) for y in outcomes]
		# compute the offset we need to ensure that - log(p(y|d)) > 0
		offset = log( ex.noiseDistribution.pdf(0.0))
		utility = prod(offset - log(marginalLikelihoods))
		return utility.item()

	def MetropolisHastingsStep( self, lastDesignPoint, lastUtility, numExperiments):
		# propose new design
		newDesignPoint = designProposal.propose(lastDesignPoint)
		# now simulate an experiment at that design point
		# draw sample from prior (well, the posterior after the previously observed data)
	
		newOutputSamples = empty(numExperiments,float)
		for k in xrange(numExperiments):
			sampleIndex = random.randint(len(ex.postSampleFunctions))
			newSampledSineFunction = ex.postSampleFunctions[sampleIndex] 
			newOutputSamples[k] = newSampledSineFunction( newDesignPoint ) + ex.noiseDistribution.rvs()

		# now evaluate the acceptance probability for the new sample
		newUtility = self.evaluateUtility( newDesignPoint, newOutputSamples)
	
		A = min( 1.0, newUtility / lastUtility)
		# print "acceptance ratio = ", A
	
		if random.random() < A:
			#accept proposal
			return (newDesignPoint, newUtility, True)
		else:
			# stay at prevsious design
			return (lastDesignPoint, lastUtility, False)


def main():
	random.seed()
	run1 = MuellerSamplerRun(1000, 10, 5)

	#subplot(maxJ, 1, j)
	figure()
	subplot(2,1,1)
	plotHist(run1.burninSamples,50,-2.0, 2.0)
	subplot(2,1,2)
	plotHist(run1.finalSamples,50,-2.0, 2.0)
	
	figure()
	plot (run1.burninSamples, range(run1.numBurnInIterations))
	show()

if __name__ == '__main__':
	main()

