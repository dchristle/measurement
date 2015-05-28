from scipy import *
import sineExample
from StaticDiscreteDistribution import StaticDiscreteDistribution
import cPickle

samplesFile = open("rejectionSamples.pydata",'r')
rejSamples = cPickle.load(samplesFile)
samplesFile.close()
# turn it into a N x 3 matrix
#rejSamples = array(rejSamples).squeeze()

example = sineExample.SineExample()
prior = example.parameterPrior
numSamples = 10000
expUtilities = []
drange = arange(-2,2.00001,0.01)
offset = log( example.noiseDistribution.pdf(0.0) )
    
for d in drange:
    utility = 0
    # evaluate the sine functions for all parameter values
    functionValues = array([ sineExample.SineFunction(*theta)(d) for theta in rejSamples])
    for i in xrange(numSamples):
        # generate a proposal sample from the kernel machine
        # pick a location
        sampleIndex = random.randint(len(rejSamples))
        chosenSample = rejSamples[sampleIndex]
        sampledFunction = sineExample.SineFunction(*chosenSample) 
        # sample a measurement y for the sampled theta
        noisyY = sampledFunction(d) + example.noiseDistribution.rvs()
        # now evaluate p(noisyY | d) by summing over all theta samples
        marginalLikelihood = mean(example.noiseDistribution.pdf(functionValues - noisyY)) 
#       print "plogp = " , plogp
        utility += (offset - log(marginalLikelihood)) / numSamples
    expUtilities.append(utility)

samplesFile = open("gridUtilities.pydata",'w')
cPickle.dump( (drange, expUtilities), samplesFile)
samplesFile.close()

