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

class SiC_Toptica_Piezo_Sweep(m2.Measurement):

    mprefix = 'topticapiezo'

    def sequence(self, upload=True, program=True, clear=False):
        gc.collect()


        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseRES = pulse.SquarePulse(channel='Sacher1160AOM', name='A square pulse on Sacher AOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

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
        # Add a microwave pulse to allow microwave energy to reach the sample even while tracking (if microwaves are enabled)
        # This will give a much more stable measurement for higher powers.
        if self.params['microwaves']:
            e.add(pulse.cp(sq_pulseMW, length=100e-6, amplitude = 1.0), name='microwave pulse', start=0e-9)
        elements.append(e)


##        # find the maximum pulse length
##        total_rf_pulses = self.params['RF_delay'] + self.params['RF_length_end'] + self.params['RF_buffer']
##        AOM_start_time = total_rf_pulses - self.params['AOM_light_delay']
##        readout_start_time = AOM_start_time + self.params['AOM_light_delay']
##        trigger_period = AOM_start_time + self.params['AOM_length'] + self.params['AOM_light_delay'] + self.params['AOM_end_buffer']
##        print 'Total trigger period is %d ns.' % trigger_period


        # Now create the resonant sequence
        e = element.Element('Resonant_mode', pulsar=qt.pulsar)

        if self.params['off_resonant_laser']:
            # Add off-resonant pulse
            e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_length']*1.0e-9), name='offresonant laser init', start=self.params['AOM_start_buffer']*1.0e-9)
            # Determine the resonant laser start time
            resonant_laser_start_time = (self.params['AOM_start_buffer'] + self.params['AOM_length'] + self.params['AOM_end_buffer'] + self.params['AOM_light_delay'] + self.params['Sacher_AOM_start_buffer'])
        else:
            resonant_laser_start_time = self.params['AOM_start_buffer']
        # Now actually add the resonant laser pulse, which depends on the previous if/then statement
        e.add(pulse.cp(sq_pulseRES, amplitude=1, length=self.params['Sacher_AOM_length']*1.0e-9), name='resonant laser', start=resonant_laser_start_time*1.0e-9)
        # Now add its readout pulse
        resonant_readout_start_time = resonant_laser_start_time + self.params['Sacher_AOM_light_delay']
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=self.params['readout_length']*1.0e-9), name='photoncountpulse', start=resonant_readout_start_time*1.0e-9)
        e.add(pulse.cp(sq_pulsePC, amplitude=0.0, length=self.params['readout_buffer']*1.0e-9), name='photonbufferend', start=(resonant_readout_start_time + self.params['readout_length'])*1.0e-9)
        # Now add a microwave pulse, if microwaves are enabled


        microwave_start_time = self.params['AOM_start_buffer'] + self.params['AOM_length'] + self.params['AOM_light_delay'] + self.params['RF_start_buffer']
        trigger_period = resonant_readout_start_time + self.params['Sacher_AOM_length']
        #total_microwave_length = resonant_laser_start_time + self.params['Sacher_AOM_length']+ self.params['Sacher_AOM_end_buffer']
        #e.add(pulse.cp(sq_pulseMW, length = total_microwave_length*1.0e-9, amplitude = 1.0), name='microwave pulse', start=0.0*1.0e-9)
        if self.params['microwaves']:
            e.add(pulse.cp(sq_pulseMW, length=self.params['pi_length']*1e-9, amplitude = 1.0), name='microwave pulse', start=microwave_start_time*1.0e-9)
        else:
            e.add(pulse.cp(sq_pulseMW, length=self.params['pi_length']*1e-9, amplitude = 0.0), name='microwave pulse', start=microwave_start_time*1.0e-9)
        if self.params['microwaves_CW']:
            e.add(pulse.cp(sq_pulseMW, length=trigger_period*1.0e-9, amplitude = 1.0), name='microwave CW pulse', start=0*1.0e-9)
        # Add the I/Q modulator pulses
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=self.params['Imod'], length=trigger_period*1.0e-9,start=0.0e-9),
        name='MWimodpulse', start=0e-9)

        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=trigger_period*1.0e-9,start=0.0e-9),
        name='MWqmodpulse', start=0e-9)

        elements.append(e)


        # Now create a CW version of the waveform where the resonant laser is kept on without the off-resonant laser.
        # The purpose of this is to allow us to measure the resonant laser power only.


        e = element.Element('ResonantCW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseRES, amplitude=1, length=100e-6), name='lasercw')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6),
        name='photoncountpulsecw')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        # Add a microwave pulse to allow microwave energy to reach the sample even while tracking (if microwaves are enabled)
        # This will give a much more stable measurement for higher powers.
        if self.params['microwaves']:
            e.add(pulse.cp(sq_pulseMW, length=self.params['pi_length']*1e-9, amplitude = 1.0), name='microwave pulse', start=0e-9)
        elements.append(e)


        seq = pulsar.Sequence('FineMotor sequence')
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
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']
        self._motdl = qt.instruments['motdl']
        self._toptica = qt.instruments['topt']
        self._fp = qt.instruments['fp']
        self._wvm = qt.instruments['bristol']
        self._flip = qt.instruments['flip']
        self._pm = qt.instruments['pm']
        self._fp = qt.instruments['fp']

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

        if self.params['motor_start'] > self.params['motor_end']:
            logging.warning('Start motor position is greater than end motor position!!!!')
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
        self._pxi.set_frequency(self.params['freq'][0]*1.0e9) # GHz units
        print 'PXI prepared, power and frequency set.'
        # Now set the proper attenuation
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten
        full_attenuation = self.params['power'] - self.params['constant_attenuation'] - np.max((0,np.min((desired_atten,15.5)))) + np.log(self.params['Imod'])/np.log(10.0)*20.0
        print 'Fully attenuated power is %.2f dBm' % full_attenuation

        # Check if wavemeter returns a valid wavelength
        wl_cur = self._wvm.get_wavelength()
        if wl_cur < 1000.0 or wl_cur > 1170.0:
            logging.error(__name__ + ': wavelength %.2f from wavemeter is not in range')
            print 'Wavelength is %.2f nm' % wl_cur

        return
    def filter_update(self, felem, seenelem):
        # first input is a list consistent of a lower bound and upper bound, from the existing filter set
        # second input is a list of list objects, which are lower/upper bounds of ranges of frequencies already
        # observed in experiment
        #
        # output is a list of lists, with each list element consistent of a lower and upper bound of frequencies
        # included in the first input element but NOT included in the second.
        for seen in seenelem:
            # deal with cases
            if felem[0] > self.params['filter_set'][l][0] and felem[0] < self.params['filter_set'][l][1]:
                if felem[1] < self.params['filter_set'][l][1]:
                    # filled in segment is contained within an existing bound, so split it up
                    #
                    # append the second segment
                    new_filter_set.append( [felem[1],self.params['filter_set'][l][1]] )
                    # change the upper bound of the first segment to be the lower bound of
                    # the range we've already scanned
                    new_filter_set.append([ self.params['filter_set'][l][0], felem[0] ])
                    del self.params['filter_set'][l]
                    self.params['filter_set'].append(temp0)
                    self.params['filter_set'].append(temp1)
                else:
                    # the upper bound of the range we've seen exceeds the existing upper bound
                    # so just change the existing upper bound to the lower bound
                    self.params['filter_set'][l] = [self.params['filter_set'][l][0],felem[0]]

            elif felem[0] < self.params['filter_set'][l][0] and felem[1] > self.params['filter_set'][l][0]:
                # the lower bound is not within the set, but the upper bound is either within the set or
                # above it
                if felem[1] < self.params['filter_set'][l][1]:
                    # the upper bound of the set we have scanned is within the existing filter range
                    #
                    # so we set the lower bound of the existing range to this upper bound
                    self.params['filter_set'][l] = [felem[1], self.params['filter_set'][l][1]]
                else:
                    # the upper bound of the set is above the filter range, so we have seen this entire
                    # filter range already
                    self.params['filter_set'][l] = []
        return
    def seek_frequency(self, frequency, step_low, step_high, motor_step = 50):
        fp = qt.instruments['fp']
        motdl = qt.instruments['motdl']
        wvm = qt.instruments['bristol']
        # this is a replacement for the root-finding based algorithm for coarse frequency tuning, e.g. to within
        # a few GHz. The idea here is to replace that algorithm with a simpler 1D line search that steeps from
        # low to high and then moves to the frequency closest to the desired one (stored in frequency).
        #
        print 'Seeking to %.2f GHz with 1D line search...' % frequency
        init_motor_position = motdl.get_position()
        motor_array = np.arange(step_low, step_high, np.round(motor_step))
        motor_freq_list = []
        for mot_pos in motor_array:
            motdl.high_precision_move(mot_pos)
            time.sleep(0.5)
            fp_out = fp.check_stabilization()
            if fp_out == 1:
                # laser is stable at this motor position, so get the frequency
                motor_freq_list.append((mot_pos, wvm.get_frequency()))
                print 'Pos: %d , Freq: %.2f GHz' % (motor_freq_list[-1][0], motor_freq_list[-1][1])
        if not motor_freq_list:
            logging.error(__name__ + '1D line search failed -- no frequencies were stable. Returning to original frequency.')
            motdl.high_precision_move(init_motor_position)
        # find the motor position nearest to the desired frequency
        nearest_pos = motor_freq_list[0][0]
        nearest_dist = motor_freq_list[0][1] - frequency
        for elem in motor_freq_list:
            if np.abs(elem[1] - frequency) < np.abs(nearest_dist):
                nearest_pos = elem[0]
                nearest_dist = elem[1] - frequency
        print 'Found motor step %d to be %.2f GHz away -- seeking.' % (nearest_pos, nearest_dist)
        motdl.high_precision_move(nearest_pos)
        print 'Distance is now %.2f GHz away from %.2f GHz.' % (wvm.get_frequency() - frequency, frequency)
        return

    def measure(self):
        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        self._stop_measurement = False
        # Wall time
        t0 = time.time()

        data = qt.Data(name='wavemotor_sweep')

        data.add_coordinate('frq (GHz)')
        for idx in range(np.size(self.params['freq'])):
            data.add_value('counts')

        plot2d_0 = qt.Plot2D(data, name='piezoscan_single_sweep', clear=True)
        data.create_file()
        # Populate some arrays
        self.params['motor_pts'] = np.uint32(1 + np.ceil(np.abs(self.params['motor_end']-self.params['motor_start'])/self.params['motor_step_size']))
        self.params['motor_array'] = np.uint32(np.linspace(self.params['motor_start'], self.params['motor_start'] + (self.params['motor_pts']-1)*self.params['motor_step_size'], self.params['motor_pts']))

        # Overwrite those arrays
        #self.params['motor_array'] = np.array(( 97610, 97760, 98650),dtype='uint32')
        #self.params['motor_array'] = np.array((97650, 98630-150, 98630, 98630+150),dtype='uint32')
        self.params['motor_pts'] = np.uint32(self.params['motor_array'].size)
        self.params['piezo_pts'] = np.uint32(1 + np.ceil(np.abs(self.params['piezo_end']-self.params['piezo_start'])/self.params['piezo_step_size']))
        #b = np.linspace(self.params['piezo_end'] + (self.params['piezo_pts']-1)*self.params['piezo_step_size'], self.params['piezo_start'], self.params['piezo_pts'])
        self.params['piezo_array'] = np.linspace(self.params['piezo_start'],self.params['piezo_end'], self.params['piezo_pts'])
        #print 'piezo array is %s' %self.params['piezo_array']
		#print '--Toptica motor/piezo scan meas. from %.3f nm to %.3f nm in %.3f nm steps (%f steps)--' % (self.params['wavelength_start'], self.params['wavelength_end'], self.params['wavelength_step_size'], self.params['motor_pts'])



        time.sleep(1.0)
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return
        if self.params['microwaves']:
            # Set the PXI status to 'on', i.e. generate microwaves
            self._pxi.set_status('on')
        N_cmeas = 0

        # Measure the CW power
        temp = np.zeros(10)
        self._awg.sq_forced_jump(3)
        self.awg_confirm(3)
        for k in range(10):
            temp[k] = self._pm.get_power()
            time.sleep(0.5)
        bg_power = np.mean(temp)
        self._flip.flip()
        time.sleep(0.5)
        for k in range(10):
            temp[k] = self._pm.get_power()
            time.sleep(0.5)
        sigbg_power = np.mean(temp)
        time.sleep(0.5)
        self._flip.flip()
        power_meas = sigbg_power-bg_power
        print 'Measured resonant CW power is %.2f nW.' % (power_meas*1.0e9)
        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)

        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

        scan_on = True
        # Start measurement cycle, so go to proper waveform.
        self._awg.sq_forced_jump(2)
        self.awg_confirm(2)

        #make sure motor_step_start is less than motor_step_end
        if self.params['motor_start'] > self.params['motor_end']:
            print 'motor_start is greater than motor_end'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return

		#determine rough frequency range
        if self.params['motor_array'][0] > 10000:
			self._motdl.set_position(self.params['motor_array'][0]-10000)
        else:
			self._motdl.set_position(0)

        self._motdl.set_position(self.params['motor_array'][0])
        self._toptica.set_piezo_voltage(self.params['piezo_array'][0])

        # Monitor laser frequency
        frq_recent = np.zeros(3)
        for zz in range(3):
            time.sleep(1.0)
            frq_recent[zz] = 299792458/self._wvm.get_wavelength()
        # Check if laser is stable
        fp_out = self._fp.check_stabilization()
        if fp_out == 1:
            print 'Laser stable according to FP.'
        elif fp_out == 2:
            print 'Laser multimode according to FP.'
        elif fp_out == 3:
            print 'No signal recorded on the FP.'
        else:
            print 'Unknown output from FP.'


		#This is the high end frequency limit, 100 GHz above the first wavemeter reading at motor_position_start we will limit our array of stored data to
        frq2 = (299792458.0/self._wvm.get_wavelength()) + 100.0 #GHz

        if self.params['motor_array'][0] > 10000:
			self._motdl.set_position(self.params['motor_array'][self.params['motor_pts']-1]-10000)
        else:
			self._motdl.set_position(0)
        self._motdl.set_position(self.params['motor_array'][self.params['motor_pts']-1])

		#This is the reference frequency we will store in the data file, 100 GHz below the wavemeter reading at motor_position_end
        time.sleep(2.0)
        # Check if laser is stable
        fp_out = self._fp.check_stabilization()
        if fp_out == 1:
            print 'Laser stable according to FP.'
        elif fp_out == 2:
            print 'Laser multimode according to FP.'
        elif fp_out == 3:
            print 'No signal recorded on the FP.'
        else:
            print 'Unknown output from FP.'

        frq1 = (299792458.0/self._wvm.get_wavelength()) - 100.0 #Ghz
        print 'Start (reference) frequency %.2f GHz / %.2f nm -- End frequency %.2f GHz / %.2f nm' % (frq1 + 100.0,299792458.0/(frq1 + 100.0), frq2 - 100.0, 299792458.0/(frq2-100.0))
        self.params['bins'] = np.uint32(1 + np.ceil(np.absolute(frq2-frq1)/self.params['bin_size']))
        #column 1 of the data set, i.e. relative frequency
        frq_array = np.linspace(-100.0, np.absolute(frq2-frq1), self.params['bins'])
        #column 2, number of counts
        total_count_data = np.zeros((np.size(frq_array),np.size(self.params['freq'])), dtype='uint32')
        #column 3, number of hits in each bin
        total_hits_data = np.zeros((np.size(frq_array),np.size(self.params['freq'])), dtype='uint32')
        if self.params['motor_array'][0] > 10000:
            self._motdl.set_position(self.params['motor_array'][0]-10000)
        else:
            self._motdl.set_position(0)

        self._motdl.set_position(self.params['motor_array'][0])
        for i in range(self.params['MeasCycles']):
            #Compensate for mechanical backlash by moving first 10000 steps behind desired motor position or to motor position 0, whichever happens first
            if self.params['motor_array'][0] > 10000:
                self._motdl.set_position(self.params['motor_array'][0]-10000)
            else:
                self._motdl.set_position(0)

            # Enter the loop for measurement
            t1 = time.time()
            for j in range(int(self.params['motor_pts'])):

                if self._stop_measurement == True:
                    print 'Measurement aborted.'
                    self._stop_measurement = True
                    break

                # Set the new motor position
                self._toptica.set_piezo_voltage(self.params['piezo_array'][0])
                if (self.params['motor_array'][j]-10000) > 0:
                    self._motdl.set_position(self.params['motor_array'][j]-10000)
                else:
                    self._motdl.set_position(0)
                time.sleep(0.2)
                print 'Motor position set to %d' % (self.params['motor_array'][j])
                self._motdl.set_position(self.params['motor_array'][j])
                time.sleep(1)



                temp_count_data = np.zeros(self.params['piezo_pts'] , dtype='uint32')

                #sweep through the piezo_array voltages, which should go from low to high and then back to low
                for k in range(np.size(self.params['piezo_array'])):

                    #Set the new piezo voltage
                    self._toptica.set_piezo_voltage(self.params['piezo_array'][k])
                    time.sleep(0.05)

                    # Measure frequency and counts
                    cur_frq = 299792458.0/self._wvm.get_wavelength()
                    offset_frq = cur_frq - frq1

                    # Check for superconductivity
                    self._snspd.check()



                    # use filter logic, if enabled
                    filter_inc = False
                    if self.params['filter']:
                        for i, elem in enumerate(self.params['filter_set']):
                            lo, hi = elem
                            if cur_frq < hi and cur_frq > lo:
                                filter_inc = True
                    else:
                        filter_inc = True
                    # Determine if we should measure in the logic statement here
                    # find all nonzero frequencies
                    if cur_frq > frq1 and cur_frq < frq2 and (np.sum(total_hits_data) < 10 or np.min(np.abs( offset_frq - frq_array[np.nonzero(total_hits_data[:,0])])) > self.params['bin_size']) and filter_inc:
                        cts_array_temp = np.zeros(np.size(self.params['freq']))
                        if self.params['stabilize_laser']:
                            fp_out = self._fp.check_stabilization()
                        else:
                            fp_out = 1

                        if (self.params['stabilize_laser'] and fp_out == 1) or not self.params['stabilize_laser']:
                            cur_frq = self._wvm.get_frequency()
                            offset_frq = cur_frq - frq1
                            for idx in range(np.size(self.params['freq'])):
                                self._pxi.set_frequency(self.params['freq'][idx]*1.0e9)
                                time.sleep(0.05)
                                cts = self._ni63.get('ctr1')
                                cts_array_temp[idx] = cts



                                #find where in the 3 column data structure to add counts
                                index = np.searchsorted(frq_array, offset_frq)

                                #update the appropriate two columns keeping track of total counts and total hits
                                total_count_data[index,idx] = total_count_data[index,idx] + cts # temp_count_data[j]
                                total_hits_data[index,idx] = total_hits_data[index,idx] + 1


                                qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.
                            cts_list = list(cts_array_temp)
                            #Live Plot
                            data.add_data_point(cur_frq,*cts_list)

                        # check snspd
                        self._snspd.check()

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
                        self._awg.sq_forced_jump(2)
                        self.awg_confirm(2)
                        time.sleep(0.1)
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


            tt = time.time() - t1

            print 'Cycle %d/%d total time is %.3f, efficiency of %.2f percent. Heater output is at %.1f. ' % (i+1, int(self.params['MeasCycles']), tt, (self.params['motor_pts'] *self.params['dwell_time']/1000.0)/tt*100.0, self._ls332.get_heater_output())

            # Sum along all sweeps so far for the y values, and just use the last frequency displacement measurement
            # for the x-axis. This is an approximation assuming the repeatability is good.
            frq_array_non0 = frq_array[np.nonzero(total_hits_data[:,0])]
            cts_array_non0 = total_count_data[np.nonzero(total_hits_data[:,0]),:]
            hits_array_non0 = total_hits_data[np.nonzero(total_hits_data[:,0]),:]
            avg_cts_array_non0 = np.divide(cts_array_non0.astype(float64),hits_array_non0.astype(float64))

            #plot2d_1 = qt.Plot2D(frq_array_non0,avg_cts_array_non0, name='topticap_avg', clear=True)
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)


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
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break
            # Checks have all passed, so proceed...



        # Set the piezo voltage back to 0
        self._toptica.set_piezo_voltage(0.0)
        # Stop PXI sig gen
        self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)
        print 'Size of freq is %d, counts is %d, avg is %d, init is %d' % (np.size(frq_array), np.size(total_count_data), np.size(total_hits_data), np.size(frq1))
        data.close_file()
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Resonant_data', self.h5data, base=self.h5base)
        grp.add('frequency', data=frq_array_non0, unit='GHz', note='frequency')
        grp.add('counts', data=cts_array_non0, unit='counts', note='total counts per sweep')
        grp.add('total_hits_array', data=hits_array_non0, unit='hits', note='how many times each bin was populated')
        grp.add('initial_measured_frequency', data=frq1, unit='hits', note='base frequency')
        grp.add('power', data=power_meas, unit='W', note='total power array')


        return



# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 180.0, # seconds
        'AOM_start_buffer' : 50.0, # ns
        'AOM_length' : 1600.0, # ns
        'AOM_light_delay' : 655.0, # ns
        'AOM_end_buffer' : 1155.0, # ns
        'Sacher_AOM_start_buffer' : 150.0, #ns
        'Sacher_AOM_length' : 3000.0, # ns
        'Sacher_AOM_light_delay' : 960.0, # ns
        'Sacher_AOM_end_buffer' : 1155.0, # ns
        'RF_start_buffer' : 300.0, # ns
        'readout_length' : 3000.0, # ns
        'readout_buffer' : 10.0, # ns
        'ctr_term' : 'PFI2',
        'motor_start' : 154350, # steps, should be lower than motor_end
        'motor_end' : 154950, # steps
        'motor_step_size' : 50, # steps
        'piezo_start' : 0, #volts
        'piezo_end' : 90, #volts
        'piezo_step_size' : 0.1, # volts (dispersion is roughly ~0.4 GHz/V)
        'bin_size' : 0.15, # GHz, should be same order of magnitude as (step_size * .1 GHz)
        'microwaves' : True, # modulate with microwaves on or off
        'microwaves_CW' : True, # are the microwaves CW? i.e. ignore pi pulse length
        'pi_length' : 180.0, # ns
        'off_resonant_laser' : True, # cycle between resonant and off-resonant
        'power' : 5.0, # dBm
        'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -9.0, # dBm
        'freq' : list([1.3357,]), #GHz
        'dwell_time' : 1000.0, # ms
        'filter' : True,
        #'filter_set' : ( (270850, 270870), (270950, 270970)),
        'filter_set' : [[264850,264960]],
        'temperature_tolerance' : 4.0, # Kelvin
        'MeasCycles' : 1,
        'Imod' : 0.2778,
        'stabilize_laser' : True
        }
def main():
    p_low = -19
    p_high = -19
    p_nstep = 1

    p_array = np.linspace(p_low,p_high,p_nstep)

    # Create a measurement object mqqq

    time.sleep(2.0)
    if msvcrt.kbhit():
        kb_char=msvcrt.getch()
        if kb_char == "q":
            do_track = False

    name_string = 'PL1_nomw'
    m = SiC_Toptica_Piezo_Sweep(name_string)
    #xsettings['desired_power'] = -19.0
    #xsettings['frequency'] = 1.45
    #xsettings['microwaves'] = False


    m.params.from_dict(xsettings)
    do_awg_stuff = True
    m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)

    m.prepare()
    time.sleep(2.0)
    if msvcrt.kbhit():
        kb_char=msvcrt.getch()
        if kb_char == "q":
            return
    print 'Proceeding with measurement ...'
    m.measure()
    m.save_params()
    m.save_stack()
    m.finish()

    #
    # name_string = 'abovetriodefect_1.42GHz_cw'
    # m = SiC_Toptica_Piezo_Sweep(name_string)
    xsettings['desired_power'] = -16
    # xsettings['frequency'] = 1.42
    # xsettings['microwaves'] = True
    #
    #
    # m.params.from_dict(xsettings)
    # do_awg_stuff = True
    # m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)
    # time.sleep(2.0)
    # if msvcrt.kbhit():
    #     kb_char=msvcrt.getch()
    #     if kb_char == "q":
    #         do_track = False
    # print 'Proceeding with measurement ...'
    # m.prepare()
    # m.measure()
    # m.save_params()
    # m.save_stack()
    # m.finish()
    #
    # name_string = 'abovetriodefect_1.3553GHz_cw'
    # m = SiC_Toptica_Piezo_Sweep(name_string)
    # xsettings['desired_power'] = -9.0
    # xsettings['frequency'] = 1.3553
    # xsettings['microwaves'] = True
    #
    #
    # m.params.from_dict(xsettings)
    # do_awg_stuff = True
    # m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)
    # time.sleep(2.0)
    # if msvcrt.kbhit():
    #     kb_char=msvcrt.getch()
    #     if kb_char == "q":
    #         do_track = False
    # print 'Proceeding with measurement ...'
    # m.prepare()
    # m.measure()
    # m.save_params()
    # m.save_stack()
    # m.finish()
    #


    # Alert that measurement has finished
    ea_t = qt.instruments['ea']
    ls332_t = qt.instruments['ls332']
    cur_temp = ls332_t.get_kelvinA()
    msg_string = 'Toptica piezo measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
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



if __name__ == "__main__":
    #Run as main program
    main()