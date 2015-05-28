#!/usr/bin/env python
# encoding: utf-8
"""
compareSamplers.py

Created by hendrik on 2006-05-24.
Copyright (c) 2006. All rights reserved.
"""

from __future__ import division
from scipy import *
from SMCAnnealingAlgorithm1  import SMCAnnealingAlgorithm1
from SMCAnnealingAlgorithm1_improved  import SMCAnnealingAlgorithm1_improved
from SMCAnnealingAlgorithm2  import SMCAnnealingAlgorithm2

import Schedules
import random 

def main():
    numParticles = 100
    numIterations = 509
    stepsPerInt = 10

	baseName = "runs/" + str(numParticles) + "parts_" + str(numIterations) + "iter_" + str(stepsPerInt) + "stepsize_"
	samplerSchedulePairs = [
	(SMCAnnealingAlgorithm1(), Schedules.LinearFractionalSchedule(stepsPerInt)),
#	(SMCAnnealingAlgorithm1_improved(), Schedules.LinearFractionalSchedule(stepsPerInt)),
	(SMCAnnealingAlgorithm2(),Schedules.Fractional0to1Schedule(stepsPerInt)),
	(SMCAnnealingAlgorithm1(), Schedules.IntegerStepsSchedule(stepsPerInt))
	]
	
	for sampler, schedule in samplerSchedulePairs:
		samplerName = sampler.__class__.__name__
		scheduleName = schedule.__class__.__name__
		print "\n========================================================"
		print "Running SMC sampler ", samplerName
		print "========================================================"
		sampler.runSampler(numParticles, numIterations, schedule)
		saveName = baseName + samplerName + "_" + scheduleName + ".pydata"
		sampler.saveSamples(saveName)
		

if __name__ == '__main__':
    main()

