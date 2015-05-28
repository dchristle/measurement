# This Python file uses the following encoding: utf-8
"""
Algorithm 2 in 
H. Kueck, N. de Freitas and Arnaud Doucet
SMC Samplers for Bayesian Optimal Nonlinear Design
Nonlinear Statistical Signal Processing Workshop (NSSPW), 2006

The last outcome is sampled according to  p(y_last|d) (C - log p(y|d))^nu_t
That is, the newly introduced outcome is first sampled from p(y_last|d). Since it is possible to sample
from this directly, no importance sampling is necessary. 
Then the distribution is slowly interpolated from p(y_last|d) to p(y_last|d) (C - log p(y|d)).
The inbetween distributions are fairly close to begin with and the variance of the importance weights therefor 
stays relatively small as a result.
"""

from __future__ import division
from scipy import *
from particleFunctions import *
import Schedules
from GenericSMCSamplerParticleState import *
from GenericSMCSampler import *
import random 

class SMCAnnealingAlgorithm2ParticleState(GenericSMCSamplerParticleState):	
	def __init__(self, design):
		GenericSMCSamplerParticleState.__init__(self,design)
								
	def evaluateLastUtilityComponent(self, marginalLikelihood, utilityFactor, nu_t):
		assert nu_t != 0.0
		return float(marginalLikelihood * ( (self.__class__.ex.maxLogLikelihood - log(marginalLikelihood)) ** nu_t) )
	
	def importanceSampleOutcome(self,outcomeIndex, nu_t ):
		assert nu_t == 0.0
		ex = self.__class__.ex
		newOutcome = random.choice(self.postFctsEvaluations) + ex.noiseDistribution.rvs()
		
		# evaluate marginal likelihood for the new outcome. 
		# this is not actually needed here. Just precomputation for later use.
		newMarginalLikelihood = float(average(ex.noiseDistribution.pdf( self.postFctsEvaluations - newOutcome )))
		# ditto for the full utility factor
		if newMarginalLikelihood == 0.0:
			newUtilityFactor = 0.0
		else:
			newUtilityFactor = float(newMarginalLikelihood * (ex.maxLogLikelihood - log(newMarginalLikelihood)))

		assert outcomeIndex == len(self.outcomes)
		# this is a brand new (not previously sampled) outcome, so we need to extend the vectors
		self.outcomes.append(newOutcome)
		self.marginalLikelihoods.append(newMarginalLikelihood)
		self.utilityFactors.append(newUtilityFactor)
		# we are sampling exactly from the target distribution in this case, so no incremental weight is necessary
		return 1.0
	
	
class SMCAnnealingAlgorithm2(GenericSMCSampler):

	ParticleType = SMCAnnealingAlgorithm2ParticleState

	def __init__(self):
		GenericSMCSampler.__init__(self)
				
	def getMCMCStaticDistributionParameters(self, n_t, nu_t, previous_n_t, previous_nu_t):
		if n_t == previous_n_t:
			# MCMC kernel is defined with invariant distribution pi_{n_t, nu_t}
			return (n_t, nu_t)
		else:
			assert n_t == previous_n_t + 1
			assert previous_nu_t == 1.0
			assert nu_t == 0.0
			return (previous_n_t, 1.0)

	def computePreMCMCStepIncrementalWeight(self, particle, n_t, nu_t, previous_n_t, previous_nu_t):
		
		if n_t == previous_n_t:
			return (SMCAnnealingAlgorithm2ParticleState.ex.maxLogLikelihood - log(particle.marginalLikelihoods[n_t - 1] )) ** (nu_t - previous_nu_t)
		else:
			assert n_t == previous_n_t + 1
			assert previous_nu_t == 1.0
			return 1.0
	
if __name__ == '__main__':
	saveName = "SMCAnnealingAlgorithm2_samples.pydata"
	numParticles = 50
	numIterations = 50
	stepsPerInt = 4
	sampler = SMCAnnealingAlgorithm2()
	schedule = Schedules.Fractional0to1Schedule(stepsPerInt)
	sampler.runSampler(numParticles, numIterations, schedule)
	sampler.saveSamples(saveName)



