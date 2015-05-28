import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
import measurement.lib.sic.SiC_ESR_Master as esr_m

# measurement parameters

xsettings = {
        'fbl_time' : 45.0, # seconds
        'ctr_term' : 'PFI0',
        'power' : -10.0, # dbM
        'f_low' : 1.26, #GHz
        'f_high' : 1.42, #Ghz
        'f_step' : 4*1.25e-4, #Ghz
        'dwell_time' : 100.0, # ms
        'temperature_tolerance' : 3.0, # Kelvin
        'MeasCycles' : 550
        }




# Create a measurement object m
m = esr_m.SiC_ESR_Master('default')

# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)


if m.review_params():
    print 'Proceeding with measurement ...'
    m.prepare()
    m.measure()
    m.save_params()
    m.save_stack()
else:
    print 'Measurement aborted!'

# important! hdf5 data must be closed, otherwise will not be readable!
# (can also be done by hand, of course)
m.finish()