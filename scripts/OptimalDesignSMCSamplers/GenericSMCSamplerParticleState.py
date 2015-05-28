from __future__ import division
import sineExample
from Proposals import *
from AcceptanceStats import AcceptanceStats
import sineExample
from scipy import *
import random

class GenericSMCSamplerParticleState:
	# class variables
	designProposal =  DesignLocalProposal()
	outcomeProposal = OutcomeProposal()
	designMHStats = AcceptanceStats()
	outcomeMHStats = AcceptanceStats()
	ex = sineExample.SineExamplePosterior()

	def __init__(self, design):
		self.design = design
		self.outcomes = []
		self.postFctsEvaluations = array([ float(f(self.design)) for f in self.__class__.ex.postSampleFunctions])
		self.marginalLikelihoods = []
		self.utilityFactors = []

	# HACK Simple handcrafted heuristic for fitting a Student-t distribution
	def fitProposal(self, postFctsEvaluations, nu_t):
		"""This function uses a heuristic to come up with a decent standard distribution to use as a 
		   importance sampling proposal for sampling from [p(y|d) log p(y|d)]^nu_t.
		   The postFctsEvaluations are the means of the normal distributions whose sum is p(y|d) """
		fctMean = mean(postFctsEvaluations)
		fctVar= var(postFctsEvaluations)
		# find parameters of proposal distribution
		# variance to aim for
		propvar = 2.0 * fctVar
		propvar += ( 2.0 / nu_t ) 
		dof = 4
		propscale = sqrt( propvar * (dof - 2.0) / dof )
		proposal = stats.t(dof, loc = fctMean, scale = propscale)
		# proposal = stats.norm(loc = fctMean, scale = sqrt(propvar))
		return proposal

	def outcomeMHStep(self,outcomeIndex, nu_t = 1):
		"""performs an MH move for one of the outcomes / auxiliary variables"""
		# propose a new outcome
		proposedOutcome = float(self.outcomeProposal.propose(self.outcomes[outcomeIndex]))
		# evaluate
		proposalMarginalLikelihood = float(average( self.__class__.ex.noiseDistribution.pdf( self.postFctsEvaluations - proposedOutcome )))
		if proposalMarginalLikelihood == 0.0:
			self.__class__.outcomeMHStats.addReject()
			return
		else: 
			proposalUtilityFactor = float(proposalMarginalLikelihood * (self.__class__.ex.maxLogLikelihood - log(proposalMarginalLikelihood)))
		
		# compute acceptance ratio
		if nu_t == 1:
			acceptanceRatio = min( 1.0, (proposalUtilityFactor * self.outcomeProposal.evaluate(self.outcomes[outcomeIndex], proposedOutcome))\
								/ (self.utilityFactors[outcomeIndex] * self.outcomeProposal.evaluate( proposedOutcome, self.outcomes[outcomeIndex])) )
		else:
			proposalUtilityEval = self.evaluateLastUtilityComponent(proposalMarginalLikelihood, proposalUtilityFactor, nu_t)
			currentUtilityEval = self.evaluateLastUtilityComponent(self.marginalLikelihoods[outcomeIndex], self.utilityFactors[outcomeIndex], nu_t)
			acceptanceRatio = min( 1.0, (proposalUtilityEval * self.outcomeProposal.evaluate(self.outcomes[outcomeIndex], proposedOutcome))\
								/ (currentUtilityEval * self.outcomeProposal.evaluate( proposedOutcome, self.outcomes[outcomeIndex])) )
			
		if random.random() < acceptanceRatio:
			#accept proposal
			self.outcomes[outcomeIndex] = proposedOutcome
			self.marginalLikelihoods[outcomeIndex] = proposalMarginalLikelihood
			self.utilityFactors[outcomeIndex] = proposalUtilityFactor
			self.__class__.outcomeMHStats.addAccept()
		else:
			self.__class__.outcomeMHStats.addReject()

	def designMHStep(self,instancesToSample, nu_t):
		# propose new design
		proposedDesignPoint = self.designProposal.propose(self.design)
		if proposedDesignPoint < -2.0 or proposedDesignPoint > 2.0:
			self.__class__.designMHStats.addReject()
			return

		# evaluate all the sine functions corresponding to the parameter samples from the prior
		# at the design point
		proposedPostFctsEvaluations = array([ float(f(proposedDesignPoint)) for f in self.__class__.ex.postSampleFunctions])

		# compute the vector p(Y|d) by using the sampled thetas to approximately integrate over theta
		proposedMarginalLikelihoods = array([ float(average( self.__class__.ex.noiseDistribution.pdf( proposedPostFctsEvaluations - y ))) \
												for y in self.outcomes])
		if any(proposedMarginalLikelihoods[0:instancesToSample] == 0.0):
			self.__class__.designMHStats.addReject()
			return			

		# compute a vector of p(y_j|d) log p(y_j|d)	 for j = 1 ... numExperiments
		proposedUtilityFactors = proposedMarginalLikelihoods * (self.__class__.ex.maxLogLikelihood - log(proposedMarginalLikelihoods))

		# we just take the product over the experiments/ auxiliary variables that the MCMC move is operating on
		proposalUtility = self.evaluateJointFctPi(proposedMarginalLikelihoods, proposedUtilityFactors, instancesToSample, nu_t)
		currentUtility = self.evaluateJointFctPi(self.marginalLikelihoods, self.utilityFactors, instancesToSample, nu_t)

		# compute acceptance probability
		acceptanceRatio = min( 1.0, (proposalUtility * self.designProposal.evaluate(self.design, proposedDesignPoint))
		/ (currentUtility * self.designProposal.evaluate( proposedDesignPoint, self.design) ) )
		if random.random() < acceptanceRatio:
			#accept proposal
			self.design = proposedDesignPoint
			self.postFctsEvaluations = proposedPostFctsEvaluations
			self.marginalLikelihoods = list(proposedMarginalLikelihoods)
			self.utilityFactors = list(proposedUtilityFactors)
			self.__class__.designMHStats.addAccept()
		else:
			# stay at previous design
			self.__class__.designMHStats.addReject()
			
	def importanceSampleOutcome(self,outcomeIndex, nu_t = 1):
		if nu_t == 1:
			return self.importanceSampleFullOutcome(outcomeIndex)
		ex = self.__class__.ex
		# Choose a student-t distribution based on the shape of p(y|d) and the exponent
		proposal = self.fitProposal(self.postFctsEvaluations, nu_t)
		# generate new proposal
		newOutcome = float( proposal.rvs() )

		# evaluate p(y_new|d, theta_i) for all theta samples
		likelihoodEvaluations = ex.noiseDistribution.pdf( self.postFctsEvaluations - newOutcome )

		# evaluate p(y_new|d), that is, evaluate the proposal distribution
		newMarginalLikelihood = float(average( likelihoodEvaluations ))
		# evaluate target
		if newMarginalLikelihood == 0.0:
			newUtilityFactor = 0.0
		else:
			newUtilityFactor = float(newMarginalLikelihood * (ex.maxLogLikelihood - log(newMarginalLikelihood)))

				
		proposalEvaluation = float(proposal.pdf(newOutcome))
		if nu_t == 1:
			targetEval = newUtilityFactor
		else:
			targetEval = self.evaluateLastUtilityComponent(newMarginalLikelihood, newUtilityFactor, nu_t)
		importanceWeight = targetEval / proposalEvaluation

		if len(self.outcomes) == outcomeIndex:
			# this is a brand new (not previously sampled) outcome, so we need to extend the vectors
			self.outcomes.append(newOutcome)
			self.marginalLikelihoods.append(newMarginalLikelihood)
			self.utilityFactors.append(newUtilityFactor)
		else:
			self.outcomes[outcomeIndex] = newOutcome
			self.marginalLikelihoods[outcomeIndex] = newMarginalLikelihood
			self.utilityFactors[outcomeIndex] = newUtilityFactor
		return importanceWeight
	
	def importanceSampleFullOutcome(self,outcomeIndex):
		ex = self.__class__.ex
		# propose from marginal p(y|d)
		newOutcome = random.choice(self.postFctsEvaluations) + ex.noiseDistribution.rvs()
		
		# evaluate marginal likelihood for the new outcome. 
		newMarginalLikelihood = float(average(ex.noiseDistribution.pdf( self.postFctsEvaluations - newOutcome )))
		# ditto for the full utility factor
		if newMarginalLikelihood == 0.0:
			newUtilityFactor = 0.0
			importanceWeight = 0.0
		else:
			newUtilityFactor = float(newMarginalLikelihood * (ex.maxLogLikelihood - log(newMarginalLikelihood)))
			importanceWeight = newUtilityFactor / newMarginalLikelihood
		
		if len(self.outcomes) == outcomeIndex:
			# this is a brand new (not previously sampled) outcome, so we need to extend the vectors
			self.outcomes.append(newOutcome)
			self.marginalLikelihoods.append(newMarginalLikelihood)
			self.utilityFactors.append(newUtilityFactor)
		else:
			self.outcomes[outcomeIndex] = newOutcome
			self.marginalLikelihoods[outcomeIndex] = newMarginalLikelihood
			self.utilityFactors[outcomeIndex] = newUtilityFactor
		return importanceWeight		

	def evaluateJointFctPi(self, marginalLikelihoods, utilityFactors, instancesToSample, nu_t):
		jointEval = prod(utilityFactors[0:instancesToSample-1])
		jointEval *= self.evaluateLastUtilityComponent(marginalLikelihoods[instancesToSample-1], utilityFactors[instancesToSample-1], nu_t)
		return jointEval

	def evaluateLastUtilityComponent(self, marginalLikelihood, utilityFactor, nu_t):
		print "Function evaluateLastUtilityComponent not defined in base class. Needs to be overridden in derived classses."
		assert False
	
	def fullMHStep(self, instancesToSample, nu_t = 1):
		"""This performs a Metropolis Hastings move with invariant distribution h_gamma_t(y,d), 
		where gamma_t=instancesToSample - 1 + nu_t.  instancesToSample is 
		the number of examples / auxiliary variables in the distribution to sample from.
		nu_t is in [0,1] and is the annealing factor for the last outcome"""
		self.designMHStep(instancesToSample, nu_t)
		for k in xrange(instancesToSample -1 ):
			self.outcomeMHStep(k)
		# sample last outcome
		self.outcomeMHStep(instancesToSample-1, nu_t)
	