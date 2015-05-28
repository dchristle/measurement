#!/usr/bin/env python
from pylab import *
import cPickle
from plotSMCsamples import plotHist

samplesFile = open("multiMuellerData.pydata",'r')
(endpoints, runs) = cPickle.load(samplesFile)
samplesFile.close()

figure(1, figsize=[10,3])
plotHist(endpoints,100,-2.0, 2.0)
# savefig('figures/mueller200_exp50.png', dpi=150)
show()
