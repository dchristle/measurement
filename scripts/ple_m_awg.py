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

class SiC_PLE_Master(m2.Measurement):

    mprefix = 'ple'

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
        seq = pulsar.Sequence('PLE sequence')
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
        self._epos = qt.instruments['epos']
        self._schr = qt.instruments['schr2']
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
        # Set Sacher temperature to 22.2 degrees.
        self._schr.set_temperature(22.2)

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

        if self.params['wavelength_start'] > self.params['wavelength_end']:
            logging.warning('Start wavelength is greater than end wavelength!!!!')
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


        print '--PLE scan meas. from %.3f nm with %d motor steps in %.3f nm steps (%d steps)--' % (self.params['start_wavelength'], self.params['wavelength_steps_array'][-1], self.params['wavelength_step_size'], self.params['pts'])

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
        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)
        time.sleep(0.1)
        #if self._fbl.optimize() == False:
        #    if self._fbl.optimize() == False:
        #        print 'FBL failed twice, breaking.'
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        scan_on = True
        # Start measurement cycle, so go to proper waveform.
        self._awg.sq_forced_jump(2)
        self.awg_confirm(2)
        self._schr.set_current_coupling_gain(self.params['current_coupling'])
        self._schr.set_current_coupling(self.params['current_coupling_enabled'])
        self._schr.set_current_coupling_direction(0)
        self._schr.set_piezo_offset(0)
        self._schr.set_piezo_status(1)




        data = qt.Data(name='ple_sweep')

        data.add_coordinate('wavelength (nm)')
        data.add_coordinate('piezo voltage (V)')
        data.add_value('frequency (GHz)')
        data.add_value('counts')

        #data.create_file()
        plot2d_0 = qt.Plot2D(data, name='finemotor_single_sweep', clear=True, coorddim=2, valdim=3)

        # The idea here is to do "shotgun spectroscopy" where we just scan over the entire
        # piezo range, step the motor, and repeat the scan. The hope is that we end up
        # getting a complete scan from taking the data and stitching it.

        # The goal is to sweep upward, find a mode hope, then sweep down, find a mode hop
        # and then return to the center of these values. The assumption is that this
        # point is where the diode is matched with the external cavity "the most"
        # and that, then, we can attempt to figure out what the best current feed-forward
        # parameters are to maintain a relatively mode hop free scanning capability.

        # First, let's sweep upward.

        # Rough conversion factor is 245 GHz/nm
        wavelength_start = self.params['wavelength_start']
        wavelength_steps_array = self.params['wavelength_steps_array'] # np.arange(0,2800,90)

        self._epos.set_wavelength(wavelength_start)
        self._schr.set_piezo_offset(0.0)
        time.sleep(0.1)
        print 'Setting wavelength to %.3f nm. Current temperature is %.2f C.' % (wavelength_start, self._schr.get_temperature())
        time.sleep(2)
        piezo_base_delta = np.linspace(0,self.params['piezo_high'],self.params['piezo_steps'])
        # This step controls entire up/back, down/back sweeping of the piezo in one array.
        piezo_delta = np.hstack((piezo_base_delta,piezo_base_delta[::-1],-1.0*piezo_base_delta,-1.0*piezo_base_delta[::-1]))

        #temperature_alternate = [22.2, 22.2-0.5385]


        # Get the current peaks from the FP cavity. These serve as the reference for the
        # rest of the experiment.

        prev_sweep = self._fp.read_sweep_peaks_lorentzian(500,10000,'ai1')

        # We can calculate the delta frequency for small deltas from this sweep.
        cur_sweep_plot = qt.Plot2D(data,'r.-',name='tcouple_plot',coorddim=1, valdim=2)
        cur_sweep_counts_plot = qt.Plot2D(data,'b.',name='ple_live_plot',coorddim=2, valdim=3)
        # Count the peaks
        N_peaks = np.size(self._fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))

        if N_peaks != 4 and N_peaks != 3:
            print '3 or 4 peaks not found -- probably not single mode to start!'
            # This checks if there are 3 or 4 peaks. If there are, it means the
            # diode is probably single mode (but not guaranteed - you have to look
            # at the actual sweep to tell, but this is a good indicator)

        # Let's start sweeping the current upward.
        scan_on = True
        last_wavelength_change_sweep = np.copy(prev_sweep)
        last_wavelength_change_frequency = 0.0
        init_motor_position = self._epos.get_motor_position()
        current_frequency = 0.0
        t1 = time.time()
        for ij in range(np.size(wavelength_steps_array)):
            self._schr.set_piezo_offset(0.0)
            print 'On iteration %d of wavelengths' % ij
            # We need to gently step the wavelength and track the frequency while we do it.
            current_motor_position = self._epos.get_motor_position()
            motor_deltas = int(np.round(((init_motor_position+wavelength_steps_array[ij])-current_motor_position)/10.0))
            temp_deltas = current_temperature + np.linspace(0,self.params['temperature_shift_per_grating_step'],motor_deltas)
            #schr.set_temperature(temperature_alternate[np.mod(ij,2)])
            #current_temperature = temperature_alternate[np.mod(ij,2)]
            for kj in range(motor_deltas):
                self._schr.set_temperature(temp_deltas[kj])


                self._epos.fine_tuning_steps(10)
                time.sleep(2)
                current_temperature = temp_deltas[kj]
                print 'Current temperature is %.3f' % current_temperature

                new_sweep = self._fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
                N_peaks = np.size(new_sweep)
                if N_peaks == 3 or N_peaks == 4:
                    df = self._fp.delta_freq_tp(prev_sweep,new_sweep)
                    prev_sweep = np.copy(new_sweep)
                    current_frequency = current_frequency + df
                    if np.abs(df) > 3.4:
                        print 'df is %.3f' % df
                    qt.msleep(0.005)
                    if msvcrt.kbhit() or scan_on == False:
                        kb_char=msvcrt.getch()
                        if kb_char == "q" or scan_on == False:
                            scan_on = False
                            break
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
                        #fbl.optimize()

                        # Set new track time
                        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
                        self._awg.sq_forced_jump(2)
                        self.awg_confirm(2)
                        time.sleep(0.1)
                #else:
                    #print 'Multimode behavior detected -- ignoring point.'
            # Now, we need to bump the motor up a bit if we're still multimode
            attempts = 0
            #Should only run if N_peaks from previous movement are not 3 or 4.
            while N_peaks != 3 and N_peaks != 4 and attempts < 15:
                self._epos.fine_tuning_steps(10)
                new_sweep = self._fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
                N_peaks = np.size(new_sweep)
                if N_peaks == 3 or N_peaks == 4:
                    df = self._fp.delta_freq_tp(prev_sweep,new_sweep)
                    prev_sweep = np.copy(new_sweep)
                    current_frequency = current_frequency + df
                    if np.abs(df) > 3.4:
                        print 'df is %.3f' % df
                    qt.msleep(0.005)
                    if msvcrt.kbhit() or scan_on == False:
                        kb_char=msvcrt.getch()
                        if kb_char == "q" or scan_on == False:
                            scan_on = False
                            break
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
                        #fbl.optimize()

                        # Set new track time
                        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
                        self._awg.sq_forced_jump(2)
                        self.awg_confirm(2)
                        time.sleep(0.1)
                #else:
                #    print 'Multimode behavior detected -- ignoring point. Attempts beyond wavelength tune are %d' % attempts



            time.sleep(5)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break


            for jj in range(np.size(piezo_delta)):
                self._schr.set_piezo_offset(piezo_delta[jj])
                time.sleep(0.5)
                N_peaks = np.size(self._fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
                new_sweep = self._fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
                if N_peaks == 3 or N_peaks == 4:
                    df = self._fp.delta_freq_tp(prev_sweep,new_sweep)
                    prev_sweep = np.copy(new_sweep)
                    current_frequency = current_frequency + df
                    if np.abs(df) > 3.4:
                        print 'df is %.3f' % df
                    counts = self._ni63.get('ctr1')
                    data.add_data_point(wavelength_steps_array[ij],piezo_delta[jj],current_frequency,counts)
                    qt.msleep(0.005)
                    if msvcrt.kbhit() or scan_on == False:
                        kb_char=msvcrt.getch()
                        if kb_char == "q" or scan_on == False:
                            scan_on = False
                            break
                #else:
                #    print 'Multimode behavior detected -- ignoring point.'

            current_wavelength_change_sweep = np.copy(prev_sweep)
            # Now let's compare the sweeps to get a delta frequency
            df = self._fp.delta_freq_tp(last_wavelength_change_sweep,current_wavelength_change_sweep)
            # Calculate the freq. discrepancy as the 'current frequency', which is integrated
            # from step-to-step, with the 'df', which we just calculated from the FP peak lcoations
            # measured at the last wavelength grating change to this one. In theory,
            # these quantities should be equal, but because are constantly adding small errors
            # together, in practice they will differ by some amount.
            freq_discrepancy = ((current_frequency - last_wavelength_change_frequency) - df)
            print 'Frequency discrepancy from accumulated error is %.4f GHz' % freq_discrepancy
            last_wavelength_change_sweep = np.copy(prev_sweep)
            # Let's now update the current frequency to be the more-accurate version
            # based on a single difference versus the many integrated differences
            current_frequency = last_wavelength_change_frequency + df
            last_wavelength_change_frequency = current_frequency

            # Now start checking for other issues. If present, stop.
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break

            tt = time.time() - t1
            print 'Total sweep time was %.2f seconds.'
            ab = data.get_data()
            frequency_axis = np.sort(ab[:,2])
            counts_axis = ab[np.argsort(ab[:,2]),3]
            plot2d_0 = qt.plot(frequency_axis,counts_axis,name='ple_sorted_sweep', clear=True)
##            if self._snspd.check() == False:
##                print 'SNSPD went normal and could not restore, breaking.'
##                break



        #data.close_file()
        self._schr.set_temperature(22.2)










        ab = data.get_data()
        # Stop PXI sig gen
        self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_PLE_data', self.h5data, base=self.h5base)
        grp.add('wavelength', data=ab[:,0], unit='step', note='wavelength steps')
        grp.add('piezo_voltage', data=ab[:,1], unit='V', note='piezo voltage applied')
        grp.add('frequency', data=ab[:,2], unit='GHz', note='frequency from FP')
        grp.add('counts', data=qb[:,3], unit='counts', note='counts at a particular frequency')



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
        'Sacher_AOM_light_delay' : 655.0, # ns
        'Sacher_AOM_end_buffer' : 1155.0, # ns
        'readout_length' : 3000.0, # ns
        'ctr_term' : 'PFI2',
        'microwaves' : False, # modulate with microwaves on or off
        'off_resonant_laser' : True, # cycle between resonant and off-resonant
        'power' : 5.0, # dBm
        'constant_attenuation' : 28.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -58.0, # dBm
        'freq' : 1.30122, #GHz
        'dwell_time' : 2000.0, # ms
        'temperature_tolerance' : 2.0, # Kelvin
        'wavelength_start' : 1106.100, # nm
        'wavelength_steps_array' : np.arange(0,5800,90), # motor steps array
        'piezo_high' : 7, # volts
        'piezo_steps' : 35, # number of steps
        'temperature_shift_per_grating_step' : -0.1795*1.5, # C/wavelength change
        }

p_low = -58
p_high = -58
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
    m = SiC_PLE_Master(name_string)
    xsettings['desired_power'] = p_array[rr]
    # since params is not just a dictionary, it's easy to incrementally load
    # parameters from multiple dictionaries
    # this could be very helpful to load various sets of settings from a global
    # configuration manager!
    m.params.from_dict(xsettings)
    do_awg_stuff = False
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
msg_string = 'PLE measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
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