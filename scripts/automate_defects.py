# automate defect measurement
import time
import numpy as np
import topticab_m_awg
from topticab_m_awg import SiC_Toptica_Search_Piezo_Sweep
import qt
import msvcrt
reload(topticab_m_awg)


fbl = qt.instruments['fbl']
fsm = qt.instruments['fsm']
xps = qt.instruments['xps']
awg = qt.instruments['awg']
reference_defect = [11.501953306181669, -0.5241993237070317, 0.9115]
defect_list = [[7.253176343561596, -0.7637003620202452, 0.913425],
               [3.175042527348811, -12.989455308869934, 0.91207],
               [4.288815677621219, -11.175531099750394, 0.911295],
               [3.0075476079999373, -5.629977213453592, 0.912605],
               [-4.274150483002216, -2.1571366809596637, 0.9121],
               [-5.867956866564423, -3.953027018230962, 0.911305],
               [-5.484955768981131, 4.651952274168486, 0.9119],
               [-5.281175700009079, 10.054927742893181, 0.91136],
               [-0.6067575861894013, 10.738420593014405, 0.912175],
               [-2.5662912369540685, 13.633409610753695, 0.91305]]

fsm.move(reference_defect[0],reference_defect[1])
xps.set_abs_positionZ(reference_defect[2])
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

drift_vector = [0, 0, 0]
reference_defect = [fsm.get_abs_positionX(), fsm.get_abs_positionY(), xps.get_abs_positionZ()]
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
    qt.instruments['motdl'].reference_search()
    data.add_data_point(i, drift_vector[0], drift_vector[1], drift_vector[2], begin_vector[0], begin_vector[1], begin_vector[2])
    xsettings = {
            'focus_limit_displacement' : 15, # microns inward
            'fbl_time' : 180.0, # seconds
            'AOM_start_buffer' : 50.0, # ns
            'AOM_length' : 1600.0, # ns
            'AOM_light_delay' : 655.0, # ns
            'AOM_end_buffer' : 1155.0, # ns
            'Sacher_AOM_start_buffer' : 150.0, #ns
            'Sacher_AOM_length' : 1000.0, # ns
            'Sacher_AOM_light_delay' : 960.0, # ns
            'Sacher_AOM_end_buffer' : 1155.0, # ns
            'RF_start_buffer' : 300.0, # ns
            'readout_length' : 1000.0, # ns
            'readout_buffer' : 10.0, # ns
            'ctr_term' : 'PFI2',
            'piezo_start' : 0, #volts
            'piezo_end' : 90, #volts
            'piezo_step_size' : 0.2, # volts (dispersion is roughly ~0.4 GHz/V)
            'bin_size' : 0.25, # GHz, should be same order of magnitude as (step_size * .1 GHz)
            'microwaves' : False, # modulate with microwaves on or off
            'microwaves_CW' : True, # are the microwaves CW? i.e. ignore pi pulse length
            'pi_length' : 180.0, # ns
            'off_resonant_laser' : True, # cycle between resonant and off-resonant
            'power' : 5.0, # dBm
            'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
            'desired_power' : -9.0, # dBm
            'freq' : [1.336], #GHz
            'dwell_time' : 900.0, # ms
            #'filter_set' : ( (270850, 270870), (270950, 270970)),(270810, 270940),
            'filter_set' : [(270740,271165)],
            'temperature_tolerance' : 6.0, # Kelvin
            'MeasCycles' : 1,
            'Imod' : 0.5,
            }

    name_string = 'defectautomate_%d' % i
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
    m.calibrate_laser()
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
