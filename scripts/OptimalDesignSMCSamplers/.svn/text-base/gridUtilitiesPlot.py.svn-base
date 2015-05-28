#!/usr/bin/env python

from scipy import *
from pylab import *
import cPickle

samplesFile = open("gridUtilities.pydata",'r')
(drange, expUtilities) = cPickle.load(samplesFile)
samplesFile.close()

# print expUtilities
exponents = [1.0, 5.0, 10.0, 20.0 , 50.0, 100.0]
figure()
for e in range(len(exponents)):
    subplot(len(exponents),1, e+1)
    plot(drange, array(expUtilities) ** exponents[e])


#figure(figsize=[10,3])
#plot(drange, array(expUtilities), '--')
#xlim(-2.0,2.0)
#yticks([])
#savefig('figures/gridUtility_exp_1_HQ.png', dpi=150)

#figure(figsize=[10,3])
#plot(drange, array(expUtilities) ** 50, 'r')
#xlim(-2.0,2.0)
#yticks([])
#savefig('figures/gridUtility_exp_50_HQ.png', dpi=150)

show()
