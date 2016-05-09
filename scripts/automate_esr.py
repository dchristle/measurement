# automate defect measurement
import time
import numpy as np
import esr_m_awg
import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
import random
from esr_m_awg import SiC_ESR_Master
import qt
import msvcrt
reload(esr_m_awg)
#reload(SiC_ESR_Master)


fbl = qt.instruments['fbl']
fsm = qt.instruments['fsm']
xps = qt.instruments['xps']
awg = qt.instruments['awg']
reference_defect =  [1.3553020999614673, 1.933735268258758, 0.813615]
updated_reference_defect = []
defect_list = [[12.787952370203088, 7.85261580315423, 0.814255],
               [7.422792977907454, 11.56063093404396, 0.813615],
               [19.656582882155963, 10.397318296305027, 0.81301],
               [17.460559117730522, 9.077662529145742, 0.81486],
               [25.350622160909673, 3.99359998574577, 0.814445],
               [25.6125976571713, 1.9260933703033514, 0.814095],
               [17.294928565627853, 5.133293889722071, 0.81415],
               [17.49421783258949, 9.051636598225887, 0.81498],
               [24.713653195060704, -0.5840262662864592, 0.81309],
               [21.722377239105732, -0.19352430763757097, 0.81222]
               ]


if updated_reference_defect != []:
    drift_vector = [updated_reference_defect[0] - reference_defect[0], updated_reference_defect[1] - reference_defect[1], updated_reference_defect[2] - reference_defect[2]]
else:
    drift_vector = [0, 0, 0]

fsm.move(reference_defect[0] + drift_vector[0],reference_defect[1] + drift_vector[1])
xps.set_abs_positionZ(reference_defect[2] + drift_vector[2])
fbl.optimize()

data = qt.Data(name='automated_position_data')
data.add_coordinate('defect index')
data.add_value('Drift Vector X')
data.add_value('Drift Vector Y')
data.add_value('Drift Vector Z')
data.add_value('FSM X')
data.add_value('FSM Y')
data.add_value('XPS Z')
data.create_file()

new_reference_position = [fsm.get_abs_positionX(), fsm.get_abs_positionY(), xps.get_abs_positionZ()]
drift_vector = [new_reference_position[0] - reference_defect[0],
                    new_reference_position[1] - reference_defect[1],
                    new_reference_position[2] - reference_defect[2]]

for i, defect in enumerate(defect_list):
    print 'Starting new defect %d in 5s...' % i
    try:
        msg_string = [__name__ + ': starting defect index %d, cryostat temperature %.2f K, heater power %.1f percent.' % (i, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
        slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
    except:
        pass
    time.sleep(5.0)
    if msvcrt.kbhit():
        kb_char=msvcrt.getch()
        if kb_char == "q":
            break
    fsm.move(defect[0]+drift_vector[0],defect[1]+drift_vector[1])
    xps.set_abs_positionZ(defect[2]+drift_vector[2])
    fbl.optimize()
    begin_vector = [fsm.get_abs_positionX(), fsm.get_abs_positionY(), xps.get_abs_positionZ()]
    print 'Do a measurement at %.2f %.2f %.2f' % (begin_vector[0], begin_vector[1], begin_vector[2])
    ##########################

    xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 55.0, # seconds
        'ctr_term' : 'PFI0',
        'power' : 5.0, # dBm
        'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -19.0, # dBm
        'f_low' : 1.2, #GHz
        'f_high' : 1.42, #Ghz
        'f_step' : 1.5e-3, #Ghz
        'dwell_time' : 1550.0, # ms
        'temperature_tolerance' : 5.0, # Kelvin
        'MeasCycles' : 6,
        'trigger_period' : 100000.0, #ns
        'dropout' : False,
        'dropout_low' : 1.20, # GHz
        'dropout_high' : 1.35, # GHz
        'Imod' : 0.35
    }

    #
    # qt.instruments['motdl'].reference_search()
    # data.add_data_point(i, drift_vector[0], drift_vector[1], drift_vector[2], begin_vector[0], begin_vector[1], begin_vector[2])
    # xsettings = {
    #         'focus_limit_displacement' : 15, # microns inward
    #         'fbl_time' : 180.0, # seconds
    #         'AOM_start_buffer' : 50.0, # ns
    #         'AOM_length' : 1600.0, # ns
    #         'AOM_light_delay' : 655.0, # ns
    #         'AOM_end_buffer' : 1155.0, # ns
    #         'Sacher_AOM_start_buffer' : 150.0, #ns
    #         'Sacher_AOM_length' : 1000.0, # ns
    #         'Sacher_AOM_light_delay' : 960.0, # ns
    #         'Sacher_AOM_end_buffer' : 1155.0, # ns
    #         'RF_start_buffer' : 300.0, # ns
    #         'readout_length' : 1000.0, # ns
    #         'readout_buffer' : 10.0, # ns
    #         'ctr_term' : 'PFI2',
    #         'piezo_start' : 0, #volts
    #         'piezo_end' : 90, #volts
    #         'piezo_step_size' : 0.2, # volts (dispersion is roughly ~0.4 GHz/V)
    #         'bin_size' : 0.25, # GHz, should be same order of magnitude as (step_size * .1 GHz)
    #         'microwaves' : False, # modulate with microwaves on or off
    #         'microwaves_CW' : True, # are the microwaves CW? i.e. ignore pi pulse length
    #         'pi_length' : 180.0, # ns
    #         'off_resonant_laser' : True, # cycle between resonant and off-resonant
    #         'power' : 5.0, # dBm
    #         'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
    #         'desired_power' : -9.0, # dBm
    #         'freq' : [1.336], #GHz
    #         'dwell_time' : 900.0, # ms
    #         #'filter_set' : ( (270850, 270870), (270950, 270970)),(270810, 270940),
    #         'filter_set' : [(270740,271165)],
    #         'temperature_tolerance' : 6.0, # Kelvin
    #         'MeasCycles' : 1,
    #         'Imod' : 0.5,
    #         }
    #
    name_string = 'defectautomate_%d' % i
    # m = SiC_Toptica_Search_Piezo_Sweep(name_string)
    # xsettings['desired_power'] = -19.0
    #name_string = 'power %.2f dBm' % (p_array[rr])
    # Create a measurement object m with a name we just made indicating the
    # power it's taken at.
    m = esr_m_awg.SiC_ESR_Master(name_string)
    # Change the xsettings dictionary above's entry for 'power' to the desired
    # power.
    xsettings['desired_power'] = -19.0


    # Load all the parameters in the slightly modified xsettings dictionary
    # into the measurement object 'm' that we just made, which will now have
    # the new power
    m.params.from_dict(xsettings)

    # The if/then here is just leftover from previous code -- since True is
    # always True, it will always execute.
    if True:
        print 'Proceeding with measurement ...'
        try:
            msg_string = [__name__ + ': CW ESR measurement %s started. Cryostat temperature %.2f K, heater power %.1f percent.' % (name_string, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
            slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
        except:
            pass
        do_awg_stuff =True
        m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)
        m.prepare()
        m.measure()
        # Save params and save stack I think just save the entire set of parameters
        # in the entire setup somewhere and also save a copy of this file every
        # time a measurement executes alongside the data, so that if it gets
        # modified, we can still go back and look at the original one.
        m.save_params()
        m.save_stack()
    else:
        print 'Measurement aborted!'

    # important! hdf5 data must be closed, otherwise will not be readable!
    # (can also be done by hand, of course)
    # m.finish() will close the HDF5 and end the measurement.
    try:
        msg_string = [__name__ + ': CW ESR measurement %s stopped. Cryostat temperature %.2f K, heater power %.1f percent.' % (name_string, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
        slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
    except:
        pass
    qt.plots['esr_avg'].save_png()
    m.finish()
    #
    #
    # m.params.from_dict(xsettings)
    # do_awg_stuff = True
    # m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)
    #
    # m.prepare()
    # time.sleep(2.0)
    # if msvcrt.kbhit():
    #     kb_char=msvcrt.getch()
    #     if kb_char == "q":
    #         break
    # print 'Proceeding with measurement ...'
    # m.calibrate_laser()
    # m.measure()
    # m.save_params()
    # m.save_stack()
    # qt.plots['piezoscan_single_sweep'].save_png()
    # m.finish()


    ##########################
    fbl.optimize()
    end_vector = [fsm.get_abs_positionX(), fsm.get_abs_positionY(), xps.get_abs_positionZ()]

    reference_defect_estimate = [reference_defect[0] + drift_vector[0] + end_vector[0] - begin_vector[0],
                                 reference_defect[1] + drift_vector[1] + end_vector[1] - begin_vector[1],
                                 reference_defect[2] + drift_vector[2] + end_vector[2] - begin_vector[2]]
    fsm.move(reference_defect_estimate[0], reference_defect_estimate[1])
    xps.set_abs_positionZ(reference_defect_estimate[2])
    fbl.optimize()
    new_reference_position = [fsm.get_abs_positionX(), fsm.get_abs_positionY(), xps.get_abs_positionZ()]
    drift_vector = [new_reference_position[0] - reference_defect[0],
                    new_reference_position[1] - reference_defect[1],
                    new_reference_position[2] - reference_defect[2]]
    data.add_data_point(i, drift_vector[0], drift_vector[1], drift_vector[2], end_vector[0], end_vector[1], end_vector[2])
    print 'Drift vector is now %.2f %.2f %.2f' % (drift_vector[0], drift_vector[1], drift_vector[2])
    try:
        msg_string = [__name__ + ': finished index %d, cryostat temperature %.2f K, heater power %.1f percent.' % (i, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
        slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
        msg_string = [__name__ + ': drift vector is %.2f um %.2f um %.2f um' % (drift_vector[0], drift_vector[1], 1000.0*drift_vector[2]) ]
        slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
    except:
        pass


data.close_file()
