import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
import pyvisa as visa
import gc
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

class SiC_Rabi_Master(m2.Measurement):

    mprefix = 'rabi'

    def sequence(self, upload=True, program=True, clear=False):
        gc.collect()


        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')


        self.params['pts'] = np.uint32(1 + np.ceil(np.abs(self.params['RF_length_end'] - self.params['RF_length_start'])/self.params['RF_length_step']))
        self.params['MW_pulse_durations'] = 1.0e-9*np.linspace(self.params['RF_length_start'], self.params['RF_length_end'], self.params['pts'])
        self._awg = qt.instruments['awg']
        self._awg.stop()
        time.sleep(5.0)
        for i in range(10):
            time.sleep(1.0)
            state = ''
            try:
                state = self._awg.get_state()
            except(visa.visa.VI_ERROR_TMO):
                print 'Waiting for AWG to stop...'
            if state == 'Idle':
                print 'AWG stopped OK.'
                break
        if clear:
            self._awg.clear_waveforms()
            print 'AWG waveforms cleared.'
        elements = []
        # First create a waveform that keeps the AOM and photon counting switch on all the time
        # but leaves the microwave switch off
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='lasercw')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6),
        name='photoncountpulsecw')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        # Add a microwave pulse to allow microwave energy to reach the sample even while tracking.
        # This will give a much more stable measurement for higher powers.
        e.add(pulse.cp(sq_pulseMW, length = self.params['MW_pulse_durations'][int(np.floor(self.params['pts']/2.0))]*1e-9, amplitude = 1.0), name='microwave pulse', start=self.params['RF_delay']*1.0e-9)
        elements.append(e)
        # find the maximum pulse length
        total_rf_pulses = self.params['RF_delay'] + self.params['RF_length_end'] + self.params['RF_buffer']
        AOM_start_time = total_rf_pulses - self.params['AOM_light_delay']
        readout_start_time = AOM_start_time + self.params['AOM_light_delay']
        trigger_period = AOM_start_time + self.params['AOM_length'] + self.params['AOM_light_delay'] + self.params['AOM_end_buffer']
        print 'Total trigger period is %d ns.' % trigger_period
        # Now create the Rabi pulses
        for i in range(self.params['pts']):

            e = element.Element('ElectronRabi_pt-%d' % i, pulsar=qt.pulsar)




            e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_length']*1.0e-9), name='laser init', start=AOM_start_time*1.0e-9)

            e.add(pulse.cp(sq_pulseMW, length = self.params['MW_pulse_durations'][i], amplitude = 1.0), name='microwave pulse', start=self.params['RF_delay']*1.0e-9)

            e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=self.params['readout_length']*1.0e-9),
            name='photoncountpulse', start=readout_start_time*1.0e-9)

            e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=trigger_period*1.0e-9),
            name='MWimodpulse', start=0e-9)

            e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=trigger_period*1.0e-9),
            name='MWqmodpulse', start=0e-9)
            elements.append(e)


        # create a sequence from the pulses
        seq = pulsar.Sequence('ElectronRabi sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
        time.sleep(4.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

        return


    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.1)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.2)
            if q == 4:
                print 'AWG not jumping... clearing VISA.'
                self._awg.clear_visa()

            if q >= 19:
                print 'AWG did not jump to proper waveform!'
        return
    def prepare(self):
        self.start_keystroke_monitor('abort')
        self._stop_measurement = False
        # Set up some instruments
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        #self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        if self._awg.get_state() == 'Idle':
            self._awg.start()
            # set the AWG to CW mode
            print 'Waiting 30 s for AWG to start...'
            time.sleep(30.0)

        for i in range(20):
            time.sleep(5.0)
            state = ''
            print 'Waiting for AWG to start...'
            try:
                state = self._awg.get_state()
            except(visa.visa.VI_ERROR_TMO):
                print 'Still waiting for AWG after timeout...'
            if state == 'Running':
                    print 'AWG started OK...Clearing VISA interface.'
                    self._awg.clear_visa()
                    break
            if state == 'Idle':
                self._awg.start()

        self._awg.sq_forced_jump(1)
        time.sleep(1)
        self.awg_confirm(1)


        self._fbl.optimize()
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))
        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))
        if self.params['RF_length_end'] < self.params['RF_length_start']:
            print 'delay_end is lower than delay_start!'


        # Set the DAQ counter dwell time, units milliseconds
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counter prepared.'
        # Reset the RFSG
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()

        # Now set the power and initial frequency
        self._pxi.set_power(self.params['power'])
        self._pxi.set_frequency(self.params['freq']*1.0e9) # GHz units
        print 'PXI prepared, power and frequency set.'
        # Now set the proper attenuation
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten


        return
    def measure(self):
        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        self._stop_measurement = False
        # Wall time
        t0 = time.time()

        # Populate some arrays

        print '--Rabi meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps)--' % (self.params['RF_length_start'], self.params['RF_length_end'], self.params['RF_length_step'], self.params['pts'] )

        total_count_data = np.zeros(self.params['pts'] , dtype='uint32')
        average_count_data = np.zeros(self.params['pts'] , dtype='float')
        intermediate_total_data = np.zeros( (1,self.params['pts'] ), dtype='uint32')
        signal = np.zeros(1)
        time.sleep(1.0)
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return
        # Set the PXI status to 'on', i.e. generate microwaves
        self._pxi.set_status('on')
        N_cmeas = 0


        # Now set the AWG into CW mode for tracking
        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)
        time.sleep(0.1)
        if self._fbl.optimize() == False:
            if self._fbl.optimize() == False:
                print 'FBL failed twice, breaking.'
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        scan_on = True
        for i in range(self.params['MeasCycles']):

            # Create an index of the waveforms so that we can modify it
            seq_index = range(self.params['pts'])
            if self.params['random'] == 1:
                # Now shuffle the array in place
                shuffle(seq_index)
            # Create an array for the single-sweep data
            temp_count_data = np.zeros(self.params['pts'] , dtype='uint32')


            # Enter the loop for measurement
            t1 = time.time()
            for j in range(int(self.params['pts'])):

                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q']:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    scan_on = False
                    self._stop_measurement = True
                    break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:

                    # set the AWG into CW mode for tracking
                    self._awg.sq_forced_jump(1)
                    self.awg_confirm(1)

                    time.sleep(0.1)
                    # Re-optimize
                    fbl.optimize()

                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()


                # Set the new RF pulse length
                self._awg.sq_forced_jump(seq_index[j]+2) # the +2 is because the indices start at 1, and the first sequence is CW mode
                time.sleep(0.1)
                if j < 5 or (j > 5 and np.random.random() < 0.01):
                    self.awg_confirm(seq_index[j]+2)

                temp_count_data[j] = self._ni63.get('ctr1')
                qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.
                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q'] or scan_on == False:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._stop_measurement = True
                    scan_on = False
                    break
                if msvcrt.kbhit() or scan_on == False or self._stop_measurement == True:
                    kb_char=msvcrt.getch()
                    self._stop_measurement = True
                    if kb_char == "q" or scan_on == False or self._stop_measurement == True:
                        print 'Measurement aborted.'
                        self._stop_measurement = True
                        break
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1

            print 'Cycle %d/%d total time is %.3f, efficiency of %.2f percent. Heater output is at %.1f. ' % (i+1, int(self.params['MeasCycles']), tt, (self.params['pts'] *self.params['dwell_time']/1000.0)/tt*100.0, self._ls332.get_heater_output())

            # Converting to numpy and back is sort of a hack, but it's one line.
            sorted_temp_data = temp_count_data[np.array(seq_index).argsort().tolist()]

            plot2d_0 = qt.Plot2D(1e9*self.params['MW_pulse_durations'],sorted_temp_data, name='rabi_single_sweep', clear=True)
            qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.
            if msvcrt.kbhit() or scan_on == False or self._stop_measurement == True:
                kb_char=msvcrt.getch()
                self._stop_measurement = True
                if kb_char == "q" or scan_on == False or self._stop_measurement == True:
                    print 'Measurement aborted.'
                    self._stop_measurement = True
                    break
            self._keystroke_check('abort')
            if self.keystroke('abort') in ['q','Q']:
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                self._stop_measurement = True
                break
            # Now start checking for other issues. If present, stop.
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
##            if self._snspd.check() == False:
##                print 'SNSPD went normal and could not restore, breaking.'
##                break
            # Checks have all passed, so proceed...

            # Now add the sorted data array to the total array
            # Use the argsort functionality to sort the count data by the frequnecy
            # it was taken at.
            total_count_data = total_count_data + temp_count_data[np.array(seq_index).argsort().tolist()]
            if i == 0:
                intermediate_total_data[0,:] = total_count_data
            elif np.mod(i,10) == 0:
                intermediate_temp_data = np.zeros( (i/10+1,self.params['pts'] ), dtype='uint32')
                intermediate_temp_data[:-1,:] = intermediate_total_data
                intermediate_temp_data[i/10,:] = total_count_data
                intermediate_total_data = np.copy(intermediate_temp_data)
                #print 'size is %s' % (np.size(intermediate_total_data))
                #intermediate_total_data = np.vstack((intermediate_total_data,total_count_data))
            if i == 0:
                signal[0] = self._ni63.get('ctr1')
            elif np.mod(i,10):
                signal = np.hstack((signal,self._ni63.get('ctr1')))
            plot2d_1 = qt.Plot2D(1e9*self.params['MW_pulse_durations'],total_count_data, name='rabi_avg', clear=True)
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)



        # Stop PXI sig gen
        self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Rabi_data', self.h5data, base=self.h5base)
        grp.add('length', data=1e9*self.params['MW_pulse_durations'], unit='ns', note='frequency')
        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')
        grp.add('intermediate', data=intermediate_total_data, unit='', note='intermediate total count data')
        grp.add('signal', data=signal, unit='counts', note='signal rate per N iterations')


        return



# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 150.0, # seconds
        'AOM_length' : 1600.0, # ns
        'AOM_light_delay' : 655.0, # ns
        'AOM_end_buffer' : 1155.0, # ns
        'RF_delay' : 10.0, # ns
        'RF_buffer' : 150.0, # ns
        'readout_length' : 130.0, # ns
        'ctr_term' : 'PFI2',
        'power' : 5.0, # dBm
        'constant_attenuation' : 28.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -9.0, # dBm
        'RF_length_start' : 0.0, # ns
        'RF_length_end' : 1350.0, # ns
        'RF_length_step' : 71.0, # ns
        'freq' : 1.30122, #GHz
        'dwell_time' : 800.0, # ms
        'temperature_tolerance' : 2.0, # Kelvin
        'MeasCycles' : 1200,
        'random' : 1
        }

p_low = -38
p_high = -38
p_nstep = 1

p_array = np.linspace(p_low,p_high,p_nstep)


for rr in range(np.size(p_array)):
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
    name_string = 'power %.2f dBm' % (p_array[rr])
    m = SiC_Rabi_Master(name_string)
    xsettings['readout_length'] = 130.0
    xsettings['desired_power'] = p_array[rr]
    # since params is not just a dictionary, it's easy to incrementally load
    # parameters from multiple dictionaries
    # this could be very helpful to load various sets of settings from a global
    # configuration manager!
    m.params.from_dict(xsettings)
    do_awg_stuff = True
    m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)


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
msg_string = 'Rabi measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)

##ps = qt.instruments['xps']
##ps.set_abs_positionZ(12.0)

track_on = True
fbl_t = qt.instruments['fbl']
track_iter = 0
print 'About to track...'
do_track = True
time.sleep(2.0)
if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q":
                    do_track = False
while track_on == True and track_iter < 50 and do_track == True:
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    fbl_t.optimize()
    time.sleep(5.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break