# automate defect measurement
import time
import numpy as np
import topticab_m_awg
from topticab_m_awg import SiC_Toptica_Search_Piezo_Sweep
import qt
import msvcrt
reload(topticab_m_awg)
reload(SiC_Toptica_Search_Piezo_Sweep)


fbl = qt.instruments['fbl']
fsm = qt.instruments['fsm']
xps = qt.instruments['xps']
awg = qt.instruments['awg']
reference_defect = [11.290182757551214, -2.255799414662634, 0.912165]
updated_reference_defect = []
defect_list = [[14.627732795919837, -1.8826998812055322, 0.914025],
               [8.472612166707194, 10.451051036769988, 0.912405],
               [5.198175350836979, 9.02641732016868, 0.912855],
               [3.4489063824270048, 9.119021594983915, 0.913545],
               [4.540505114378166, 12.317735822056965, 0.912515],
               [1.1053366176145076, 11.657985811046823, 0.911505],
               [-11.06842109939018, 10.982440465412877, 0.912685],
               [-12.223746592000941, 10.289880423288599, 0.9137],
               [-9.040777331305296, 9.739555916079818, 0.913305],
               [-8.876878219650559, 8.330024022909406, 0.91122]
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
