import numpy as np
import qt
import time
import PL_calibrate
reload(PL_calibrate)
# Photoluminescence measurement script

# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 65.0, # seconds
        'frequency' : 1.233, # GHz
        'detuning' : -0.15, # GHz
        'wavelength' : 1150.0, # nm
        'grating' : 2, # grating number
        'exposure_time' : 1.0, # seconds
        'microwaves' : False, # Boolean
        'power' : 5.0, # dBm
        'constant_attenuation' : 6.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -9.0, # dBm
        'x_displacement' : -1.7, # microns
        'y_displacement' : -1.7, # microns
        'ctr_term' : 'PFI0',
        'temperature_tolerance' : 3.0, # Kelvin
        'random' : 1
        }


# Create a measurement object m
print 'About to proceed -- waiting 5 s for quit (press q to quit)'
time.sleep(5.0)
name_string = 'power %.2f dBm' % (p_array[rr])
m = PL_calibrate.SiC_SpectrumCalibrate_Master(name_string)
xsettings['exposure_time'] = 60.0 # seconds
xsettings['MeasCycles'] = 20
# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)
m.sequence(upload=True,program=True)


if False:
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
#msg_string = 'PL measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
#ea_t.email_alert(msg_string)
#
###ps = qt.instruments['xps']
###ps.set_abs_positionZ(12.0)
#
#track_on = True
#fbl_t = qt.instruments['fbl']
#track_iter = 0
#while track_on == True and track_iter < 50:
#    track_iter = track_iter + 1
#    print 'Tracking for %d iteration.' % track_iter
#    fbl_t.optimize()
#    time.sleep(5.0)
#    if msvcrt.kbhit() or track_on == False:
#                kb_char=msvcrt.getch()
#                if kb_char == "q" or track_on == False: break