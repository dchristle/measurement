# automate defect measurement
import time
import numpy as np
import topticab_m_awg
from topticab_m_awg import SiC_Toptica_Search_Piezo_Sweep
import qt
import msvcrt
import esr_m_awg
import esr_classifier
reload(esr_classifier)
reload(topticab_m_awg)
reload(esr_m_awg)
#reload(SiC_Toptica_Search_Piezo_Sweep)


fbl = qt.instruments['fbl']
fsm = qt.instruments['fsm']
xps = qt.instruments['xps']
awg = qt.instruments['awg']
reference_defect = [-8.073212356770348, 3.4745444712606366, 0.6299]
updated_reference_defect = [-6.412349875374224, 4.091744224256404, 0.63112]
defect_list = [
                [3.463065610105171, 9.489066635947196, 0.63048],
                [15.982791057710967, 5.714561701386428, 0.63075],
                [8.270766361381119, 10.847944341475827, 0.63148],
                [2.09966849307239, 1.335163411403997, 0.63348],


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

use_defect_names = True
defect_name_list = [ '1K', '1N', '1Q','1G', ]

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
    #
    #qt.instruments['motdl'].reference_search()


    # NOW DO THE TOPTICA MEASUREMENT IF CLASSIFIER OUTPUTS 1
    time.sleep(5.0)
    fbl.optimize()
    if True:
        print 'Starting Toptica scan on defect %d in 5s...' % i
        try:
            msg_string = [__name__ + ': starting defect index %d Toptica scan, cryostat temperature %.2f K, heater power %.1f percent.' % (i, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
            slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
        except:
            pass
        time.sleep(5.0)
        fbl.optimize()
        if msvcrt.kbhit():
            kb_char=msvcrt.getch()
            if kb_char == "q":
                break

        xsettings = {
                'focus_limit_displacement' : 15, # microns inward
                'fbl_time' : 120.0, # seconds
                'AOM_start_buffer' : 50.0, # ns
                'AOM_length' : 1600.0, # ns
                'AOM_light_delay' : 655.0, # ns
                'AOM_end_buffer' : 1155.0, # ns
                'Sacher_AOM_start_buffer' : 150.0, #ns
                'Sacher_AOM_length' : 1000.0, # ns
                'Sacher_AOM_light_delay' : 625.0, # ns
                'Sacher_AOM_end_buffer' : 1155.0, # ns
                'RF_start_buffer' : 300.0, # ns
                'readout_length' : 1000.0, # ns
                'readout_buffer' : 10.0, # ns
                'ctr_term' : 'PFI2',
                'piezo_start' : 0, #volts
                'piezo_end' : 90, #volts
                'piezo_step_size' : 0.25, # volts (dispersion is roughly ~0.4 GHz/V)
                'bin_size' : 0.18, # GHz, should be same order of magnitude as (step_size * .1 GHz)
                'filter_threshold' : 1.5, # GHz
                'microwaves' : False, # modulate with microwaves on or off
                'microwaves_CW' : False, # are the microwaves CW? i.e. ignore pi pulse length
                'pi_length' : 30.0, # ns
                'off_resonant_laser' : True, # cycle between resonant and off-resonant
                'power' : 5.0, # dBm
                'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
                'desired_power' : -19.0, # dBm
                'freq' : [1.3820,], #GHz
                'dwell_time' : 400.0, # ms
                #'filter_set' : ( (270850, 270870), (270950, 270970)),(270810, 270940),
                'filter_set' : [(265150,265490.0),],#, (270951,270974)],
                'current_range' : [0.275, 0.305], # A
                'temperature_tolerance' : 3.0, # Kelvin
                'MeasCycles' : 1,
                'Imod' : 0.3,
                'stabilize_laser' : False,
                }

        if use_defect_names:
            name_string = 'defectautomatePLE_%s' % defect_name_list[i]
        else:
            name_string = 'defectautomatePLE_%d' % i

        m = SiC_Toptica_Search_Piezo_Sweep(name_string)
        xsettings['desired_power'] = -19.0


        m.params.from_dict(xsettings)
        do_awg_stuff = True
        m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)

        m.prepare()
        time.sleep(2.0)
        if msvcrt.kbhit():
            kb_char=msvcrt.getch()
            if kb_char == "q":
                break
        print 'Proceeding with measurement ...'
        #m.calibrate_laser()
        m.measure()
        m.save_params()
        m.save_stack()
        qt.plots['piezoscan_single_sweep'].save_png()
        m.finish()


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
