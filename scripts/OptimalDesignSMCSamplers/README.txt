The code included in this distribution implements the algorithms described in 
    H. Kueck, N. de Freitas and Arnaud Doucet
    SMC Samplers for Bayesian Optimal Nonlinear Design
    Nonlinear Statistical Signal Processing Workshop (NSSPW), 2006

Take a look at GenericSMCSampler.py for the implementation of the basic SMC samplers algorithm.
SMCAnnealingAlgorithm1, SMCAnnealingAlgorithm1_improved and SMCAnnealingAlgorithm2 are
specializations (implemented as subclasses) of this algorithm. The headers of these 
files contain short descriptions of the algorithms.

Proposals.py contains the problem specific proposal distributions used by the different
samplers.

Schedules.py provides different types of annealing schedules used by the SMC sampler
algorithms.

Please take a look at compareSamplers.py for how to use the SMC sampler algorithms. 

comparisonPlots.py provides different ways of plotting the results of different runs.

muellerSampler.py implements the MCMC annealing approach proposed in 
    P.Mueller, B. Sanso, and M. de Iorio, 
    “Optimal Bayesian design by inhomogeneous Markov chain simulation,” 
    Journal of the American Statistical Association, vol. 99, pp. 788–798, 2004.

plotPriors.py shows the priors on the parameters for the toy problem described in our paper.
plotPostSamples.py plots the sampled posterior distribution (both in parameter space as well
as the corresponding sine functions) after observing 2 initial data points. The experiment
design problem in this simple example then is to find the optimal x location for the 3rd
measurement.

I apologize for the rather incomplete documentation, the messy code and the lack of comments 
in many parts. Please do not hesitate to contact me at kueck at cs.ubc.ca with any questions 
you might have.
