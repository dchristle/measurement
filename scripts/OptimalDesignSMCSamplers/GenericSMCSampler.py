from __future__ import division
from scipy import *
from particleFunctions import *
import cPickle
from GenericSMCSamplerParticleState import GenericSMCSamplerParticleState
import scipy.stats as stats

class GenericSMCSampler:

	def __init__(self):
		self.sampleStore = []
		self.weightStore = []
		self.n_tStore = []
		self.essLog = []
		self.particles = []
		self.weights = []
		self.incrementalWeightStore = []
		self.resampleCount = 0
		# setting error handling for scipy
		# seterr(under = 'warn', over = 'warn', divide = 'warn', invalid = 'warn', where = 2)

	def initParticles(self, numParticles, initial_n_t, initial_nu_t):
		assert initial_n_t == 1
		ParticleType = self.__class__.ParticleType
		# draw uniformly in the interval [-2,2]
		positions = list(stats.uniform(-2,4).rvs(size = numParticles))
		for p in positions:
			newParticle = ParticleType(p)
			weight = newParticle.importanceSampleOutcome(0, initial_nu_t)
			self.particles.append(newParticle)
			self.weights.append(weight)

	def getMCMCStaticDistributionParameters(self, n_t, nu_t, previous_n_t, previous_nu_t):
		print "getMCMCStaticDistributionParameters needs to be defined in derived classes"
		assert 0

	def computePreMCMCStepIncrementalWeight(self, particle, n_t, nu_t, previous_n_t, previous_nu_t):
		return 1.0

	def conditionalResample(self, essThreshold):
		ess = computeESS(self.weights)
		print "ESS = ", ess
		if ess < essThreshold:
			print "resampling"
			# do a first resampling based on the weights
			self.particles, self.weights = resample(self.particles, self.weights)
			self.resampleCount += 1
			return (ess, True)
		return (ess, False)

	def runSampler(self,  numParticles, numIterations, schedule, relativeESSThreshold = 0.6):
		self.nameOfScheduleUsed = schedule.__class__.__name__
		essThreshold = relativeESSThreshold * numParticles
		random.seed()
		# initialize particles (defined in subclasses)
		initial_n_t, initial_nu_t = schedule(0)
		self.initParticles(numParticles, initial_n_t, initial_nu_t )
		ess, resampled = self.conditionalResample(essThreshold)
		self.essLog.append( (ess, 0 ))
		if resampled:
			self.essLog.append( (numParticles, 0+0.001 ))

		# now starts the SMC sampling algorithm itself
		for i in xrange(1,numIterations+1):
			# The schedule determines n_t and nu_t
			# n_t is the number of outcomes to be sampled
			# nu_t is the annealing factor for the last component.
			# the last component is slowly faded in by increasing nu_t from 0 to 1
			n_t, nu_t = schedule(i)
			previous_n_t, previous_nu_t = schedule(i-1)

			# get the parameters of the static distribution that the MCMC kernel is defined on
			mcmc_n_t, mcmc_nu_t = self.getMCMCStaticDistributionParameters(n_t, nu_t, previous_n_t, previous_nu_t)

			# make sure that all except for maybe the last component are updated using the MCMC kernel
			assert mcmc_n_t + 1 >= n_t

			# just for debugging
			print " "
			print "iteration ", i, ":  n_t =" ,n_t , "  nu_t =", nu_t
			print "Using MCMC to sample first ",mcmc_n_t," indices  nu_t for MCMC chain = ", mcmc_nu_t
			if n_t > mcmc_n_t:
				print "importance sampling last index ",n_t - 1, " with annealing factor ", nu_t

			incrementalWeights = []
			for p in xrange(len(self.particles)):
				incrementalWeight = 1.0
				particle = self.particles[p]
				# TODO compute the weight compensating for MCMC kernels
				incrementalWeight *= self.computePreMCMCStepIncrementalWeight(particle, n_t, nu_t, previous_n_t, previous_nu_t)
				# do a MH step for each particle (for all pre-existing experiments)
				particle.fullMHStep(mcmc_n_t, mcmc_nu_t)
				if n_t > mcmc_n_t:
					# sample the last component using importance sampling
					incrementalWeight *= particle.importanceSampleOutcome(n_t-1, nu_t)

				self.weights[p] *= incrementalWeight
				incrementalWeights.append(incrementalWeight)
			self.incrementalWeightStore.append(incrementalWeights)

			# output and reset acceptance stats
			print "design   MH acceptance ratio = ", GenericSMCSamplerParticleState.designMHStats.getAcceptanceRatio()
			print "outcomes MH acceptance ratio = ", GenericSMCSamplerParticleState.outcomeMHStats.getAcceptanceRatio()
			GenericSMCSamplerParticleState.designMHStats.reset()
			GenericSMCSamplerParticleState.outcomeMHStats.reset()

			ess, resampled = self.conditionalResample(essThreshold)
			self.essLog.append( (ess, i ))
			if resampled:
				self.essLog.append( (numParticles, i+0.001 ))

			# only store samples when we hit an integer exactly
			# these are samples from the marginal distributions we are interested in,
			# U(d)^J for increasing integer powers J
			# the second condition makes sure that if we stay at the same n_t, nu_t values for several
			# iterations (as is the case with the integer step schedule) we only store samples for the
			# last iteration before moving on to a different target distribution
			if nu_t == 1.0 and schedule(i) != schedule(i+1):
				self.sampleStore.append( array([float(p.design) for p in self.particles]) )
				self.weightStore.append( array( self.weights ) )
				self.n_tStore.append(n_t)


	def saveSamples(self, samplesFileName ):
		samplesFile = open(samplesFileName,'w')
		resultsDict = {}
		resultsDict['sampleStore'] = self.sampleStore
		resultsDict['weightStore'] = self.weightStore
		resultsDict['n_tStore'] = self.n_tStore
		resultsDict['essLog'] = self.essLog
		resultsDict['resampleCount'] = self.resampleCount
		resultsDict['incrementalWeightStore'] = self.incrementalWeightStore
		resultsDict['samplerName'] = self.__class__.__name__
		resultsDict['scheduleName'] = self.nameOfScheduleUsed
		cPickle.dump(resultsDict, samplesFile)
		samplesFile.close()


