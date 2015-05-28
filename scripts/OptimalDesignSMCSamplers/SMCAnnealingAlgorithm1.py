# This Python file uses the following encoding: utf-8
"""
Implementation of the SMC sampler described in algorithm 3 in 
A. M. Johansen, A. Doucet, and M. Davy. 
Maximum likelihood parameter estimation for maximum likelihood models using sequential Monte Carlo.
In Proceedings of ICASSP, 2006.

This also corresponds to Algorithm 1 in 
H. Kueck, N. de Freitas and Arnaud Doucet
SMC Samplers for Bayesian Optimal Nonlinear Design
Nonlinear Statistical Signal Processing Workshop (NSSPW), 2006
"""

from __future__ import division
from scipy import *
from particleFunctions import *
import Schedules
from GenericSMCSamplerParticleState import *
from GenericSMCSampler import *


class SMCAnnealingAlgorithm1ParticleState(GenericSMCSamplerParticleState):	
	def __init__(self, design):
		GenericSMCSamplerParticleState.__init__(self,design)
		
	def evaluateJointFctPi(self, marginalLikelihoods, utilityFactors, instancesToSample, nu_t):
		# in this algorithm, MCMC is only defined on the fully active outcomes
		assert nu_t == 1
		return prod(utilityFactors[0:instancesToSample])
		
	def evaluateLastUtilityComponent(self, marginalLikelihood, utilityFactor, nu_t):
		return utilityFactor ** nu_t
	
class SMCAnnealingAlgorithm1 (GenericSMCSampler):

	ParticleType = SMCAnnealingAlgorithm1ParticleState

	def __init__(self):
		GenericSMCSampler.__init__(self)
			
	def getMCMCStaticDistributionParameters(self, n_t, nu_t, previous_n_t, previous_nu_t):
		if previous_nu_t == 1.0:
			mcmc_n_t = previous_n_t
		else: 
			mcmc_n_t = previous_n_t - 1
		mcmc_nu_t = 1
		return (mcmc_n_t, mcmc_nu_t)
	
	
if __name__ == '__main__':
	saveName = "SMCAnnealingAlgorithm1_samples.pydata"
	numParticles = 100
	numIterations = 200
	stepsPerInt = 4
	sampler = SMCAnnealingAlgorithm1()
	schedule = Schedules.LinearFractionalSchedule(stepsPerInt)
	sampler.runSampler(numParticles, numIterations, schedule)
	sampler.saveSamples(saveName)



