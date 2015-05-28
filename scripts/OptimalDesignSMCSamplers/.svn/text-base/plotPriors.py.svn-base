from pylab import *
from scipy import *

from sineExample import SineParameterPrior

def plotDistribution( rng ):
    [mu, var] = rng.stats()
    stdd = sqrt(var)
    minx = mu - 4 * stdd
    maxx = mu + 4 * stdd
    xr = mgrid[minx:maxx:1000j]
    yr = rng.pdf(xr)
    plot(xr,yr)
    

if __name__ == '__main__':
    prior = SineParameterPrior();
    plotsize = [6,3]
    figure(1, figsize=plotsize)
    plotDistribution(prior.freqPrior)
    title('Frequency prior')
    figure(2, figsize=plotsize)
    plotDistribution(prior.ampPrior)
    title('Amplitude prior')
    figure(3, figsize=plotsize)
    plotDistribution(prior.phasePrior)
    title('Phase prior')
    
    #figDPI = 300
    #figure(1)
    #savefig('figures/freqPrior.png', dpi=figDPI)
    #figure(2)
    #savefig('figures/ampPrior.png', dpi=figDPI)
    show()
