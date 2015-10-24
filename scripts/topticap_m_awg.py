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
            resonant_laser_start_time = (self.params['AOM_start_buffer'] + self.params['AOM_length'] + self.params['AOM_end_buffer'])
        else:
            resonant_laser_start_time = self.params['AOM_start_buffer']
        # Now actually add the resonant laser pulse, which depends on the previous if/then statement
        e.add(pulse.cp(sq_pulseRES, amplitude=1, length=self.params['Sacher_AOM_length']*1.0e-9), name='resonant laser', start=resonant_laser_start_time*1.0e-9)
        # Now add its readout pulse
        resonant_readout_start_time = resonant_laser_start_time + self.params['Sacher_AOM_light_delay']
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=self.params['readout_length']*1.0e-9), name='photoncountpulse', start=resonant_readout_start_time*1.0e-9)
        # Now add a microwave pulse, if microwaves are enabled
        if self.params['microwaves']:
            total_microwave_length = resonant_laser_start_time + self.params['Sacher_AOM_length']+ self.params['Sacher_AOM_end_buffer']
            e.add(pulse.cp(sq_pulseMW, length = total_microwave_length*1.0e-9, amplitude = 1.0), name='microwave pulse', start=0.0*1.0e-9)
            # Add the I/Q modulator pulses
            e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=total_microwave_length*1.0e-9),
            name='MWimodpulse', start=0e-9)

            e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=total_microwave_length*1.0e-9),
            name='MWqmodpulse', start=0e-9)

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

        data = qt.Data(name='wavemotor_sweep')

        data.add_coordinate('frq (GHz)')
        data.add_value('counts')
        plot2d_0 = qt.Plot2D(data, name='piezoscan_single_sweep', clear=True)
        data.create_file()
        # Populate some arrays
        self.params['motor_pts'] = np.uint32(1 + np.ceil(np.abs(self.params['motor_end']-self.params['motor_start'])/self.params['motor_step_size']))
        self.params['motor_array'] = np.uint32(np.linspace(self.params['motor_start'], self.params['motor_start'] + (self.params['motor_pts']-1)*self.params['motor_step_size'], self.params['motor_pts']))
        self.params['piezo_pts'] = np.uint32(1 + np.ceil(np.abs(self.params['piezo_end']-self.params['piezo_start'])/self.params['piezo_step_size']))
        a = np.linspace(self.params['piezo_start'],self.params['piezo_end'], self.params['piezo_pts'])
        #b = np.linspace(self.params['piezo_end'] + (self.params['piezo_pts']-1)*self.params['piezo_step_size'], self.params['piezo_start'], self.params['piezo_pts'])
        self.params['piezo_array'] = a
        print 'piezo array is %s' %self.params['piezo_array']
		#print '--Toptica motor/piezo scan meas. from %.3f nm to %.3f nm in %.3f nm steps (%f steps)--' % (self.params['wavelength_start'], self.params['wavelength_end'], self.params['wavelength_step_size'], self.params['motor_pts'])

        #total_count_data = np.zeros(self.params['motor_pts'] , dtype='uint32')
        #average_count_data = np.zeros(self.params['motor_pts'] , dtype='float')
        #current_count_data = np.zeros(self.params['motor_pts'] ,dtype='uint32')
        #frequency_displacement = np.zeros(self.params['motor_pts'] ,dtype='float')

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


        # Now set the AWG into CW mode for tracking
        #self._awg.sq_forced_jump(1)
        #self.awg_confirm(1)
        #time.sleep(0.1)
        #if self._fbl.optimize() == False:
        #    if self._fbl.optimize() == False:
        #        print 'FBL failed twice, breaking.'
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
        # Check if laser is stable, if not, wait
        for zz in range(30):
            if (np.max(frq_recent) - np.min(frq_recent)) < 1.0:
                print 'Laser stabilized. Dispersion %.2f MHz, range %.2f MHz' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                break
            time.sleep(1.0)
            np.roll(frq_recent,1)
            frq_recent[0] = 299792458/self._wvm.get_wavelength()
            if zz >=29:
                print 'Laser frequency readout on wavemeter did not stabilize after 30 seconds!'

		#This is the high end frequency limit, 100 GHz above the first wavemeter reading at motor_position_start we will limit our array of stored data to
        frq2 = (299792458.0/self._wvm.get_wavelength()) + 100.0 #GHz

        if self.params['motor_array'][0] > 10000:
			self._motdl.set_position(self.params['motor_array'][self.params['motor_pts']-1]-10000)
        else:
			self._motdl.set_position(0)
        self._motdl.set_position(self.params['motor_array'][self.params['motor_pts']-1])

		#This is the reference frequency we will store in the data file, 100 GHz below the wavemeter reading at motor_position_end
        time.sleep(2.0)
        print 'Checking for laser stabilization...'
        # Monitor laser frequency
        frq_recent = np.zeros(5)
        for zz in range(5):
            time.sleep(1.0)
            frq_recent[zz] = 299792458.0/self._wvm.get_wavelength()
        # Check if laser is stable, if not, wait
        for zz in range(30):
            if (np.max(frq_recent) - np.min(frq_recent)) < 1.0:
                print 'Laser stabilized. Dispersion %.2f MHz, range %.2f MHz' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                break
            time.sleep(1.0)
            np.roll(frq_recent,1)
            frq_recent[0] = 299792458.0/self._wvm.get_wavelength()
            if zz >=29:
                print 'Laser frequency readout on wavemeter did not stabilize after 30 seconds! Dispersion %.2f GHz, range %.2f GHz.' % (np.std(frq_recent), np.max(frq_recent)-np.min(frq_recent))

        frq1 = (299792458.0/self._wvm.get_wavelength()) - 100.0 #Ghz
        print 'Start (reference) frequency %.2f GHz / %.2f nm -- End frequency %.2f GHz / %.2f nm' % (frq1 + 100.0,299792458.0/(frq1 + 100.0), frq2 - 100.0, 299792458.0/(frq2-100.0))
        self.params['bins'] = np.uint32(1 + np.ceil(np.absolute(frq2-frq1)/self.params['bin_size']))
        #column 1 of the data set, i.e. relative frequency
        self.params['frq_array'] = np.linspace(-100.0, np.absolute(frq2-frq1), self.params['bins'])
        #column 2, number of counts
        total_count_data = np.zeros(np.size(self.params['frq_array']), dtype='uint32')
        #column 3, number of hits in each bin
        total_hits_data = np.zeros(np.size(self.params['frq_array']), dtype='uint32')
        if self.params['motor_array'][0] > 10000:
			self._motdl.set_position(self.params['motor_array'][0]-10000)
        else:
			self._motdl.set_position(0)

        self._motdl.set_position(self.params['motor_array'][0])
        for i in range(self.params['MeasCycles']):

            # Create array for the single-sweep data
            #temp_count_data = np.zeros(self.params['motor_pts'] , dtype='uint32')
            # Set Sacher wavelength
            #self._epos.set_wavelength(self.params['wavelength_array'][0])
            #time.sleep(1.0)

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
                self._motdl.set_position(self.params['motor_array'][j])
                time.sleep(2)
                # Monitor laser frequency
                frq_recent = np.zeros(5)
                for zz in range(5):
                    time.sleep(1.0)
                    frq_recent[zz] = 299792458/self._wvm.get_wavelength()
                # Check if laser is stable, if not, wait
                for zz in range(30):
                    if (np.max(frq_recent) - np.min(frq_recent)) < 1.0:
                        print 'Laser stabilized. Dispersion %.2f MHz, range %.2f MHz' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                        break
                    time.sleep(1.0)
                    np.roll(frq_recent,1)
                    frq_recent[0] = 299792458/self._wvm.get_wavelength()
                    if zz >=29:
                        print 'Laser frequency readout on wavemeter did not stabilize after 30 seconds!'


                temp_count_data = np.zeros(self.params['piezo_pts'] , dtype='uint32')

                #sweep through the piezo_array voltages, which should go from low to high and then back to low
                for k in range(np.size(self.params['piezo_array'])):
                    #Set the new piezo voltage
                    self._toptica.set_piezo_voltage(self.params['piezo_array'][k])
                    time.sleep(0.05)
                    #Measure frequency and counts
                    cur_frq = 299792458.0/self._wvm.get_wavelength()
                    frq = cur_frq - (frq1 + 100.0) #Ghz
                    self._snspd.check()
                    # this will ensure points aren't retaken -- but it will also cause subsequent scans to skip,
                    # so it's a bit of a hack, just to try this out for now. the hack is that if the current frequency
                    # is nearby any previously sampled frequency to within 1 bin width, we skip it.

                    if cur_frq > frq1 and cur_frq < frq2:
                        cts = self._ni63.get('ctr1')

    				    #Live Plot
                        data.add_data_point(frq,cts)
                        # Keep track of the last frequency for settling purposes
                        last_frq = frq

    				    #find where in the 3 column data structure to add counts
                        index = np.searchsorted(self.params['frq_array'], frq)

    				    #update the appropriate two columns keeping track of total counts and total hits
                        total_count_data[index] = total_count_data[index] + temp_count_data[j]
                        total_hits_data[index] = total_hits_data[index] + 1


                        qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.

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
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1

            print 'Cycle %d/%d total time is %.3f, efficiency of %.2f percent. Heater output is at %.1f. ' % (i+1, int(self.params['MeasCycles']), tt, (self.params['motor_pts'] *self.params['dwell_time']/1000.0)/tt*100.0, self._ls332.get_heater_output())

            # Sum along all sweeps so far for the y values, and just use the last frequency displacement measurement
            # for the x-axis. This is an approximation assuming the repeatability is good.
            frq_array_non0 = self.params['frq_array'][np.nonzero(total_hits_data)]
            cts_array_non0 = total_count_data[np.nonzero(total_hits_data)]
            hits_array_non0 = total_hits_data[np.nonzero(total_hits_data)]
            avg_cts_array_non0 = np.divide(cts_array_non0.astype(float64),hits_array_non0.astype(float64))

            plot2d_1 = qt.Plot2D(frq_array_non0,avg_cts_array_non0, name='topticap_avg', clear=True)
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





        # Stop PXI sig gen
        self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)
        print 'Size of freq is %d, counts is %d, avg is %d, init is %d' % (np.size(self.params['frq_array']), np.size(total_count_data), np.size(total_hits_data), np.size(frq1))
        data.close_file()
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_WaveMotor_data', self.h5data, base=self.h5base)
        grp.add('frequency', data=frq_array_non0, unit='GHz', note='frequency')
        grp.add('counts', data=cts_array_non0, unit='counts', note='total counts per sweep')
        grp.add('total_hits_array', data=hits_array_non0, unit='hits', note='how many times each bin was populated')
        grp.add('initial_measured_frequency', data=frq1, unit='hits', note='base frequency')


        return



# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 50.0, # seconds
        'AOM_start_buffer' : 50.0, # ns
        'AOM_length' : 1600.0, # ns
        'AOM_light_delay' : 655.0, # ns
        'AOM_end_buffer' : 1155.0, # ns
        'Sacher_AOM_length' : 3000.0, # ns
        'Sacher_AOM_light_delay' : 950.0, # ns
        'Sacher_AOM_end_buffer' : 1155.0, # ns
        'readout_length' : 3000.0, # ns
        'ctr_term' : 'PFI2',
        'motor_start' : 92000, # steps, should be lower than motor_end
        'motor_end' : 98000, # steps
        'motor_step_size' : 150, # steps
        'piezo_start' : 0, #volts
        'piezo_end' : 90, #volts
        'piezo_step_size' : 0.16, # volts (dispersion is roughly ~0.4 GHz/V)
        'bin_size' : 0.03, # GHz, should be same order of magnitude as (step_size * .1 GHz)
        'microwaves' : False, # modulate with microwaves on or off
        'off_resonant_laser' : True, # cycle between resonant and off-resonant
        'power' : 5.0, # dBm
        'constant_attenuation' : 28.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -28.0, # dBm
        'freq' : 1.3358, #GHz
        'dwell_time' : 2200.0, # ms
        'temperature_tolerance' : 4.0, # Kelvin
        'MeasCycles' : 1,
        }

p_low = -25
p_high = -25
p_nstep = 1

p_array = np.linspace(p_low,p_high,p_nstep)


for rr in range(np.size(p_array)):
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(1.0)
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
    name_string = 'power %.2f dBm' % (p_array[rr])
    m = SiC_Toptica_Piezo_Sweep(name_string)
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
msg_string = 'Wavemotor measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)

##ps = qt.instruments['xps']
##ps.set_abs_positionZ(12.0)

##track_on = True
##fbl_t = qt.instruments['fbl']
##track_iter = 0
##print 'About to track...'
##do_track = True
##time.sleep(2.0)
##if msvcrt.kbhit():
##                kb_char=msvcrt.getch()
##                if kb_char == "q":
##                    do_track = False
##while track_on == True and track_iter < 50 and do_track == True:
##    track_iter = track_iter + 1
##    print 'Tracking for %d iteration.' % track_iter
##    fbl_t.optimize()
##    time.sleep(5.0)
##    if msvcrt.kbhit() or track_on == False:
##                kb_char=msvcrt.getch()
##                if kb_char == "q" or track_on == False: break