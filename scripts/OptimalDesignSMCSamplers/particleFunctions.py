from __future__ import division
from scipy import *
from copy import deepcopy

def resample( particles, weights ):
	# draw a single random value in [0,1)
	randVal = random.random()
	weightCDF = cumsum( weights )
	weightSum = weightCDF[-1]
	
	numParticles = len(particles)
	newParticles = []
	stdWeight = 1.0 / numParticles
	newWeights = ones(numParticles) * stdWeight
	
	curParticleIdx	= 0
	for i in xrange(numParticles):
		# u is regularly sampled within the CDF range, with random offset
		u = weightSum * (( i + randVal) / numParticles )
		while ( u > weightCDF[curParticleIdx] ):
			curParticleIdx += 1
		# possible optimization: Only deep copy when duplicate particle
		# but for now, play is save
		# FIXME This doesnt seem to work
		newParticles.append(deepcopy(particles[curParticleIdx]))
	return (newParticles, newWeights)

def computeESS(weights):
	normalizedWeights = weights / sum(weights)
	ess = 1.0 / sum (normalizedWeights ** 2)
	return ess


def countDifferentParticles( particles ):
	count = 0
	lastParticle = -inf
	for p in particles:
		if p != lastParticle:
			lastParticle = p
			count += 1
	return count

