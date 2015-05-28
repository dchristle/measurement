# This Python file uses the following encoding: utf-8
"""
This is an improved version of SMCAnnealingAlgorithm1.
This version of the algorithm was not actually described in our workshop paper due to space constraints.

The difference with SMCAnnealingAlgorithm1 is that this improved version uses a MCMC kernel 
with static distribution defined on all outcomes sampled at the previous iteration, 
whereas Algorithm 1 as described in our paper and implemented in SMCAnnealingAlgorithm1 only defines the 
MCMC sampler on all fully faded in outcomes at the previous timestep (that is, it includes 
outcomes that were sampled using importance sampling with previous_nu_t != 1.0 )

As a result, only when an integer boundary is crossed will the improved sampler use importance sampling 
for the new outcome. Algorithm 1 on the other hand will always use importance sampling for 
the last outcome. This results in way higher variance of the importance weights.
"""

from __future__ import division
from scipy import *
from particleFunctions import *
import Schedules
from GenericSMCSamplerParticleState import *
from GenericSMCSampler import *


class SMCAnnealingAlgorithm1_improvedParticleState(GenericSMCSamplerParticleState):	
	def __init__(self, design):
		GenericSMCSamplerParticleState.__init__(self,design)
								
	def evaluateLastUtilityComponent(self, marginalLikelihood, utilityFactor, nu_t):
		return utilityFactor ** nu_t
	
class SMCAnnealingAlgorithm1_improved (GenericSMCSampler):

	ParticleType = SMCAnnealingAlgorithm1_improvedParticleState

	def __init__(self):
		GenericSMCSampler.__init__(self)
				
	def getMCMCStaticDistributionParameters(self, n_t, nu_t, previous_n_t, previous_nu_t):
		if n_t == previous_n_t:
			# MCMC kernel is defined with invariant distribution pi_{n_t, nu_t}
			return (n_t, nu_t)
		else:
			assert n_t == previous_n_t + 1
			assert previous_nu_t == 1.0
			return (previous_n_t, 1.0)

	def computePreMCMCStepIncrementalWeight(self, particle, n_t, nu_t, previous_n_t, previous_nu_t):
		
		if n_t == previous_n_t:
			return particle.utilityFactors[n_t - 1] ** (nu_t - previous_nu_t)
		else:
			assert n_t == previous_n_t + 1
			assert previous_nu_t == 1.0
			#return particle.utilityFactors[previous_n_t - 1] ** (1.0 - previous_nu_t)
			return 1.0
	
if __name__ == '__main__':
	saveName = "SMCAnnealingAlgorithm1_improved_samples.pydata"
	numParticles = 100
	numIterations = 200
	stepsPerInt = 4
	sampler = SMCAnnealingAlgorithm1_improved()
	schedule = Schedules.LinearFractionalSchedule(stepsPerInt)
	sampler.runSampler(numParticles, numIterations, schedule)
	sampler.saveSamples(saveName)



