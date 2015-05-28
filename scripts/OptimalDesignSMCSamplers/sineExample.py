from scipy import *
from pylab import *
import scipy.stats as stats
import cPickle
import math
import dateutil
import matplotlib.pyplot as plot
import matplotlib.figure as figure




class SineParameterPrior:
	def __init__(self):
		self.freqPrior = stats.gamma(100, loc=0, scale = 0.01)
		self.ampPrior = stats.gamma(10, loc=0, scale = 1)
		self.phasePrior = stats.uniform(loc=0,scale = 1)
	def rvs(self,**kwargs):
		return ( self.freqPrior.rvs(**kwargs).item(), self.ampPrior.rvs(**kwargs).item(),self.phasePrior.rvs(**kwargs).item() )
	def pdf(self, theta):
		return self.freqPrior.pdf(theta[0]) * self.ampPrior.pdf(theta[1]) * self.phasePrior.pdf(theta[2])

class SineLocalExplorationKernel:
	def __init__(self):
		self.freqKernel = stats.norm(loc = 0, scale = 0.05)
		self.ampKernel = stats.norm(loc = 0, scale = 0.5)
		self.phaseKernel = stats.norm(loc = 0, scale = 0.05)
	def rvs(self,**kwargs):
		return ( self.freqKernel.rvs(**kwargs), self.ampKernel.rvs(**kwargs),self.phaseKernel.rvs(**kwargs))
	def pdf(self, theta):
		return self.freqKernel.pdf(theta[0]) * self.ampKernel.pdf(theta[1]) * self.phaseKernel.pdf(theta[2])
	def evaluate(self, theta1, theta2):
		return self.pdf(theta1-theta2)


class SineExample:
	def __init__(self):
		self.parameterPrior = SineParameterPrior()
		self.noiseDistribution = stats.norm(loc = 0, scale = 0.8)
	def sampleFunction(self):
		(f,a,p) = self.parameterPrior.rvs()
		return SineFunction(f,a,p)
	def evaluateLikelihood( self, x, y, theta):
		"""compute the likelihood of data for parameters theta
		   the data consists of input(s) x and output(s) y"""
		fctValues = SineFunction(*theta)(x)
		if theta[2] < 0.0 or theta[2] > 1.0:
			return 0.0
		return product( self.noiseDistribution.pdf( fctValues - y) )
	def evaluateLogLikelihood( self, x, y, theta):
		"""compute the log likelihood of data for parameters theta
		   the data consists of input(s) x and output(s) y"""
		if theta[2] < 0.0 or theta[2] > 1.0:
			return 0.0
		fctValues = SineFunction(*theta)(x)
		return sum( log(self.noiseDistribution.pdf( fctValues - y) ) )

class SineExamplePosterior:
	def __init__(self):
		self.noiseDistribution = stats.norm(loc = 0, scale = 0.8)
		self.maxLogLikelihood = log(self.noiseDistribution.pdf(0.0))
		self.loadPosteriorSamples()

	def loadPosteriorSamples(self):
		samplesFile = open("rejectionSamples.pydata",'r')
		self.postSamples = cPickle.load(samplesFile)
		samplesFile.close()
		self.postSampleFunctions = [ SineFunction(*theta) for theta in self.postSamples]


class SineFunction:
	def __init__(self, freq, amp, phase):
		self.f = freq
		self.a = amp
		self.p = phase
	def	 __call__(self, x):
		return self.a * sin( ((x * self.f) + self.p) * 2 * pi )
	def __str__(self):
		return str(self.a) + " * sin( " + str(self.f) + " * " +str( 2 * pi) + " * x )"


def plotPrior():
	figure()
	subplots_adjust(hspace = 0.4)
	s = SineExample( )
	thingsToPlot = []
	xr = mgrid[-2:2:1000j]
	hold(True)
	for i in xrange(100):
		theta = s.parameterPrior.rvs()
		f = SineFunction(*theta)
		yr = f(xr)
		subplot(3,1,1)
		plot(xr,yr)
		subplot(3,1,2)
		plot((theta[0],), (theta[1],),'kx')
		xlabel('Frequency')
		ylabel('Amplitude')
		subplot(3,1,3)
		plot((theta[0],), (theta[2],),'kx')
		xlabel('Frequency')
		ylabel('Phase')
	hold(False)
	show()

def plotPriorSingleWindows():
	s = SineExample( )
	thingsToPlot = []
	xr = mgrid[-2:2:1000j]
	figure(1, figsize=[10,3])
	figure(2, figsize=[6,6])
	figure(3, figsize=[6,6])
	hold(True)
	for i in xrange(100):
		theta = s.parameterPrior.rvs()
		f = SineFunction(*theta)
		yr = f(xr)
		figure(1)
		plot(xr,yr)
		figure(2)
		plot((theta[0],), (theta[1],),'kx')
		figure(3)
		plot((theta[0],), (theta[2],),'kx')
	figure(2)
	xlabel('Frequency')
	ylabel('Amplitude')
	figure(3)
	xlabel('Frequency')
	ylabel('Phase')
	hold(False)
	saveFigures()
	show()

def saveFigures():
	figDPI = 300
	figure(1)
	savefig('figures/priorFcts.png', dpi=figDPI)
	figure(2)
	savefig('figures/priorParams1.png', dpi=figDPI)
	figure(3)
	savefig('figures/priorParams2.png', dpi=figDPI)

if __name__ == '__main__':
	plotPriorSingleWindows()
