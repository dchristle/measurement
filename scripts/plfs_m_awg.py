import numpy as np
import qt
import time
import msvcrt
from measurement.lib.sic import singlespin
reload(singlespin)
# Photoluminescence measurement script

# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 95.0, # seconds
        'M_low' : 15, # in mm
        'M_high' : 91.0, # in mm
        'M_step' : 0.074, # in mm
        'ctr_term' : 'PFI0',
        'dwell_time' : 4000.0, # in ms
        'temperature_tolerance' : 2.0, # Kelvin
        'random' : 1,
        'fast' : True,
        }


# Create a measurement object m
print 'About to proceed -- waiting 5 s for quit (press q to quit)'
time.sleep(5.0)
name_string = 'no_microwave'
m = singlespin.SiC_PLFieldsweep_Master(name_string)

xsettings['MeasCycles'] = 1
xsettings['random'] = False
# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)
m.sequence(upload=False,program=False)


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
#ea_t = qt.instruments['ea']
#ls332_t = qt.instruments['ls332']
#cur_temp = ls332_t.get_kelvinA()
#msg_string = 'PL fieldsweep measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
#ea_t.email_alert(msg_string)

##ps = qt.instruments['xps']
##ps.set_abs_positionZ(12.0)

# Now just keep tracking until 'q' is pressed
track_on = True
fbl_t = qt.instruments['fbl']
track_iter = 0
while track_on == True:
    time.sleep(1.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    fbl_t.optimize()
    time.sleep(1.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break