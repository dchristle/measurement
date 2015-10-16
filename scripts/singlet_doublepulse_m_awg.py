# Doublepulse PL measurement script
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
        'fbl_time' : 120.0, # seconds
        'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
        'AOM_light_delay' : 655.0, # ns
        'AOM_start_buffer' : 155.0, # ns
        'AOM_end_buffer' : 20.0, # ns
        'AOM_init_length' : 2000.0, # ns
        'AOM_readout_length' : 1600.0, # ns
        'PH_trigger_time' : 0.0, #ns
        'PH_trigger_length' : 50.0, #ns
        'power' : 5.0, # dBm
        'tau_length_start' : 0.0, # ns
        'tau_length_end' : 800.0, # ns
        'tau_length_step' : 15, # ns
        'microwaves' : False, # Boolean
        'frequency' : 1.3358, #GHz
        'desired_power' : -9.0, # dBm
        'pi_length' : 51.0, # ns
        'CFDLevel0' : 125,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 110,
        'CFDZeroCross1' : 10,
        'Binning' : 5,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : -10000,
        'acquisition_time' : 60.0, # s
        'temperature_tolerance' : 2.0, # Kelvin
        'MeasCycles' : 300,
        'random' : 1
        }


tlca = np.array((0.22,))


for cur in tlca:

    tli = qt.instruments['tl']
    tli.set_current(cur)
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)
    if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
    name_string = '3C_microwave_%.2f A' % cur
    m = singlespin.SiC_DoublePulse_Master(name_string)
    xsettings['MeasCycles'] = 170
    # since params is not just a dictionary, it's easy to incrementally load
    # parameters from multiple dictionaries
    # this could be very helpful to load various sets of settings from a global
    # configuration manager!
    m.params.from_dict(xsettings)
    m.sequence(upload=True,program=True,clear=True)


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
    print 'Quit now before next measurement'
    time.sleep(5.0)
    if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break

# Alert that measurement has finished
ea_t = qt.instruments['ea']
ls332_t = qt.instruments['ls332']
cur_temp = ls332_t.get_kelvinA()
msg_string = 'Doublepulse measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)


xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 120.0, # seconds
        'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
        'AOM_light_delay' : 655.0, # ns
        'AOM_start_buffer' : 155.0, # ns
        'AOM_end_buffer' : 40.0, # ns
        'AOM_init_length' : 2000.0, # ns
        'AOM_readout_length' : 1600.0, # ns
        'PH_trigger_time' : 0.0, #ns
        'PH_trigger_length' : 50.0, #ns
        'power' : 5.0, # dBm
        'tau_length_start' : 1800.0, # ns
        'tau_length_end' : 1800.0, # ns
        'tau_length_step' : 1, # ns
        'microwaves' : True, # Boolean
        'frequency' : 1.3358, #GHz
        'desired_power' : -9.0, # dBm
        'pi_length' : 51.0, # ns
        'CFDLevel0' : 125,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 110,
        'CFDZeroCross1' : 10,
        'Binning' : 5,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : -10000,
        'acquisition_time' : 60.0, # s
        'temperature_tolerance' : 2.0, # Kelvin
        'MeasCycles' : 300,
        'random' : 1
        }


tlca = np.array((0.22,))


for cur in tlca:

    tli = qt.instruments['tl']
    tli.set_current(cur)
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)
    if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
    name_string = '3C_microwave_%.2f A' % cur
    m = singlespin.SiC_DoublePulse_Master(name_string)
    xsettings['MeasCycles'] = 170
    # since params is not just a dictionary, it's easy to incrementally load
    # parameters from multiple dictionaries
    # this could be very helpful to load various sets of settings from a global
    # configuration manager!
    m.params.from_dict(xsettings)
    m.sequence(upload=True,program=True,clear=True)


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
    print 'Quit now before next measurement'
    time.sleep(5.0)
    if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break

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
while track_on == True and track_iter < 50:
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    time.sleep(1.0)
    if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
    fbl_t.optimize()
    time.sleep(5.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break