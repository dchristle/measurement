from pylab import *
from scipy import *
import sineExample

ex = sineExample.SineExamplePosterior()
rejSamples = ex.postSamples
from rejPostSampling import observedData

def plotPosterior():
	figure()
	hold(True)
	xr = mgrid[-2:2:1000j]
	for theta in rejSamples:
	    f = sineExample.SineFunction(*theta)
	    yr = f(xr)
	    subplot(3,1,1)
	    plot(xr,yr)
	    plot(observedData[:,0], observedData[:,1], 'ro')
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


def plotPosteriorSingleWindows():
	close('all')
	thingsToPlot = []
	xr = mgrid[-2:2:1000j]
	figure(1, figsize=[10,3])
	figure(2, figsize=[6,6])
	figure(3, figsize=[6,6])
	hold(True)
	for theta in rejSamples:
		f = sineExample.SineFunction(*theta)
		yr = f(xr)
		figure(1)
		plot(xr,yr)
		figure(2)
		plot((theta[0],), (theta[1],),'kx')
		figure(3)
		plot((theta[0],), (theta[2],),'kx')
	figure(1)
	plot(observedData[:,0], observedData[:,1], 'ro')
	figure(2)
	xlabel('Frequency')
	ylabel('Amplitude')
	figure(3)
	xlabel('Frequency')
	ylabel('Phase')
	hold(False)
	# saveFigures()
	show()

def saveFigures():
	figDPI = 150
	figure(1)
	savefig('figures/posteriorFcts.png', dpi=figDPI)
	figure(2)
	savefig('figures/posteriorParams1.png', dpi=figDPI)
	figure(3)
	savefig('figures/posteriorParams2.png', dpi=figDPI)
	

if __name__ == '__main__':	  
	plotPosteriorSingleWindows()
