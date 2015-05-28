from scipy import *
    
class StaticDiscreteDistribution:
    def __init__(self, weights, normalized = False):
        self.numElements = len(weights)
        # we choose to normalize such that the sum of the weights
        # is N, the number of elements
        if normalized:
            self.weights = len(weights) * array(weights, float)
        else:
            weightSum = sum(weights)
            self.weights = (len(weights) / weightSum) * array(weights,float)
        
        # self.mainIndices = array(range(self.numElements))
        self.aliasIndices = array(range(self.numElements))
        
        # the average weight of an element. As a result of the normalization
        # this is 1.0
        averageWeight = 1.0 
        
        # split indices into those that have above and below average weights
        smallerSet = list(nonzero(less_equal(self.weights, averageWeight)))
        largerSet = list(nonzero(greater(self.weights, averageWeight)))
        
        # now redistribute the weights such that each index carries
        # exactly weight 1.0
        while len(smallerSet) > 0 and len(largerSet) > 0:
            # grab one large and one small index
            smallIdx = smallerSet.pop()
            largeIdx = largerSet.pop()
            # fill up the smaller element
            self.aliasIndices[smallIdx] = largeIdx
            # by taking from the larger one
            self.weights[largeIdx] -= averageWeight - self.weights[smallIdx]
            # now the smaller index is done with. It carries
            # exactly 1/N probability mass
            # Now figure out what to do with the previously larger element
            # which now has lost probability mass
            if self.weights[largeIdx] > averageWeight:
                largerSet.append(largeIdx)
            else:
                smallerSet.append(largeIdx)
        
        # everything has been distributed. The remaining entries thus are neither
        # too small nor too big => they all have exactly weight 1.0
        # => they have no alias. 
        self.weights[smallerSet] = 1.0
        self.weights[largerSet] = 1.0
    def sample(self):
        # choose index
        chosenIndex = random.randint(self.numElements)
        # lookup the decision boundary stored for that element
        thres = self.weights[chosenIndex]
        if random.random() > thres:
            return self.aliasIndices[chosenIndex]
        else:
            return chosenIndex
    
    def rvs(self, size = 1):
        return array([ self.sample() for i in xrange(size)])


if __name__ == "__main__":    
    from pylab import *    
    # test code
    N = 8
    numSamples = 10000
    for i in xrange(10):
        randWeights = exp(random.uniform(0, 10, size=N))
        #randWeights = ones(N,float)
        #randWeights[5] = 10
        sdist = StaticDiscreteDistribution(randWeights)
        #draw samples
        samples = sdist.rvs(size=numSamples)
        counts = [ sum( equal(samples, j)) for j in xrange(N)]
        expectedCounts = randWeights / sum(randWeights) * numSamples
        print array2string(expectedCounts, suppress_small=1)
        print counts
        print ""
    
    