# Photoluminescence measurement script
# David J. Christle <christle@uchicago.edu>
# 2015/01/19
#
# This script sets up a double pulse measurement, with or without the pi pulse,
# across an array of delay times between the pulses. The purpose of this experiment
# is to study both the recovery time (implying something about the singlet) and
# the shape of the recovery curve, which depends in a complicated way on the
# various rates in the system.

# measurement parameters

import numpy as np
import qt
import time
import msvcrt
from measurement.lib.sic import singlespin
reload(singlespin)

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 60.0, # seconds
        'step_low' : 0, # step
        'step_high' : 300000, # step
        'step_step' : 2000, # step
        'speed' : 2000, # standa speed
        'ctr_term' : 'PFI0', # terminal
        'dwell_time' : 1000, # ms
        'temperature_tolerance' : 2.0, # Kelvin
        'MeasCycles' : 1,
        'random' : 1
        }


# Create a measurement object m
print 'About to proceed -- waiting 5 s for quit (press q to quit)'
time.sleep(5.0)
name_string = '3C_polarization'
m = singlespin.SiC_PLPolarizationBasic_Master(name_string)
#xsettings['exposure_time'] = 90.0 # seconds
xsettings['MeasCycles'] = 1
# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)
m.sequence(upload=False,program=False,clear=False)


if True:

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

# Alert that measurement has finished
ea_t = qt.instruments['ea']
ls332_t = qt.instruments['ls332']
cur_temp = ls332_t.get_kelvinA()
msg_string = 'Doublepulse measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)

##ps = qt.instruments['xps']
##ps.set_abs_positionZ(12.0)

track_on = True
fbl_t = qt.instruments['fbl']
track_iter = 0
while track_on == True and track_iter < 10:
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    fbl_t.optimize()
    time.sleep(5.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break