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

class SiC_Toptica_TCSPC(m2.Measurement):

    mprefix = 'topticatcspc'

    def sequence(self, upload=True, program=True, clear=False):
        gc.collect()


        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseRES = pulse.SquarePulse(channel='Sacher1160AOM', name='A square pulse on Sacher AOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')
        sq_pulsePH = pulse.SquarePulse(channel='phtrigger', name='A square pulse on phtrigger')

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
        e.add(pulse.cp(sq_pulsePH, amplitude=-0.7, length=self.params['PH_trigger_length']*1.0e-9), name='picoharp trigger', start=self.params['PH_trigger_time']*1.0e-9)
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
            e.add(pulse.cp(sq_pulseMW_Imod, amplitude=self.params['Imod'], length=total_microwave_length*1.0e-9),
            name='MWimodpulse', start=0e-9)

            e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=total_microwave_length*1.0e-9),
            name='MWqmodpulse', start=0e-9)

        elements.append(e)

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


        seq = pulsar.Sequence('TopticaTCSPC sequence')
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
        self._ph = qt.instruments['ph']
        self._ls = qt.instruments['ls']

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        if self._awg.get_state() == 'Idle':
            self._awg.start()
            # set the AWG to CW mode
            for i in range(20):
                time.sleep(10.0)
                state = ''
                print 'Waiting for AWG to start...'
                try:
                    state = self._awg.get_state()
                except(visa.visa.VI_ERROR_TMO):
                    print 'Still waiting for AWG after timeout...'

                if state == 'Running':
                        print 'AWG started OK.'
                        break

                if state == 'Idle':
                    self._awg.start()
        elif self._awg.get_state() == 'Running':
            print 'AWG already started'

        print 'Checking AWG interface status...'
        init_status = self._awg.get_ch3_status()
        awg_is_ok = True
        # check a randomized sequence of 20 commands and see if the response agrees with what we set the channel
        # status to.
        for ij in range(20):
            test_val = random.randint(0,1)
            if test_val == 0:
                self._awg.set_ch3_status('off')
                time.sleep(0.05)
                output = self._awg.get_ch3_status()
                if output == 'on':
                    awg_is_ok == False
            if test_val == 1:
                self._awg.set_ch3_status('on')
                time.sleep(0.05)
                output = self._awg.get_ch3_status()
                if output == 'off':
                    awg_is_ok == False
        self._awg.set_ch3_status(init_status)
        # all the randomized commands should match. If they don't, the VISA interface is clogged and we should
        # clear it
        if not awg_is_ok:
            print 'AWG interface was clogged, clearing it.'
            self._awg.clear_visa()
        else:
            print 'AWG interface is OK.'

        # feedback on laser resonance
        self._awg.sq_forced_jump(3)
        time.sleep(1.0)
        self._ls.optimize()

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


        # Set the DAQ counter dwell time, units milliseconds
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counter prepared.'
        # Configure PH for histogram mode, initialize
        self._ph.start_histogram_mode()

        self._ph.set_Binning(self.params['Binning'])
        self._ph.set_InputCFD0(self.params['CFDLevel0'],self.params['CFDZeroCross0'])
        self._ph.set_InputCFD1(self.params['CFDLevel1'],self.params['CFDZeroCross1'])
        self._ph.set_SyncOffset(self.params['SyncOffset'])
        print 'PicoHarp settings configured.'
        # Check for counts on PH
        if self._ph.get_CountRate1() > 0:
            print 'PH 1 channel receiving counts...'
        else:
            print 'PH 1 not receiving counts!'

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
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten
        full_attenuation = self.params['power'] - self.params['constant_attenuation'] - np.max((0,np.min((desired_atten,15.5)))) + np.log(self.params['Imod'])/np.log(10.0)*10
        print 'Fully attenuated power is %.2f dBm' % full_attenuation

        # Check if wavemeter returns a valid wavelength
        wl_cur = self._wvm.get_wavelength()
        if wl_cur < 1000.0 or wl_cur > 1170.0:
            logging.error(__name__ + ': wavelength %.2f from wavemeter is not in range')
            print 'Wavelength is %.2f nm' % wl_cur


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
            self._stop_measurement = False
            return
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))

        if self._snspd.check():
            print 'SNSPD is superconducting.'
        else:
            print 'SNSPD is not superconducting!'
            return
        print 'Press q now to abort.'
        time.sleep(3.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return
        return
    def measure(self):
        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        self._stop_measurement = False
        # Wall time
        t0 = time.time()


        histogram_array = np.zeros(65536)
        # draw a blank plot, which we will populate later

        time.sleep(0.05)
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return
        if self.params['microwaves']:
            # Set the PXI status to 'on', i.e. generate microwaves
            self._pxi.set_status('on')
        N_cmeas = 0


        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        scan_on = True
        # Start measurement cycle, so go to proper waveform.
        self._awg.sq_forced_jump(2)
        self.awg_confirm(2)

        for jj in range(self.params['MeasCycles']):
            t1 = time.time()
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" :
                    scan_on = False
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
                self._awg.sq_forced_jump(3)
                self.awg_confirm(3)
                # feedback on laser
                self._ls.optimize()
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

            signal_begin = self._ni63.get_ctr0()

            self._ph.ClearHistMem()
            self._ph.StartMeas(int(self.params['acquisition_time']*1000))
            print 'Acquiring signal for %s seconds. ' % (self.params['acquisition_time'])
            time.sleep(self.params['acquisition_time']+0.25)

            n = 0
            while self._ph.get_MeasRunning() and n < 10:
                time.sleep(0.5)
                n = n + 1
            if self._ph.get_MeasRunning():
                print 'Measurement did not finish!'
                break
            # Retrieve measurement
            current_data = self._ph.get_Histogram()

            print 'Got data for sweep %d -- total sum is %d/%d' % (jj, np.sum(current_data),np.sum(histogram_array))
            histogram_array = np.copy(histogram_array) + current_data
            signal_end = self._ni63.get_ctr0()
            # First clear the plot

            # Now just add those arrays to the now-cleared plot

            # Make sure to update it

            # qt.msleep() is also good to do for making sure everything is up to date
            qt.msleep(0.002)
            #plot2d.add(qt.Data(np.array((np.linspace(0,65535*np.power(2,self.params['Binning'])*0.004,65536),histogram_array))))
            data_array = np.zeros((65536,2))
            data_array[:,0] = np.linspace(0,65535*np.power(2,self.params['Binning'])*0.004,65536)
            data_array[:,1] = np.copy(histogram_array)

            plot2d = qt.Plot2D(data_array,name='toptica_tcspc',maxpoints=70000)
            plot2d.update()
            qt.msleep(0.005)
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1
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

            print 'Cycle %d/%d total time is %.3f, efficiency of %.2f percent. Heater output is at %.1f. ' % (jj+1, int(self.params['MeasCycles']), tt, (self.params['acquisition_time'])/tt*100.0, self._ls332.get_heater_output())



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
        #self._toptica.set_piezo_voltage(0.0)
        # Stop PXI sig gen
        self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)

        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_TCSPC', self.h5data, base=self.h5base)
        #grp.add('frequency', data=frq_array_non0, unit='GHz', note='frequency')
        grp.add('counts', data=histogram_array, unit='counts', note='total counts')



        return



# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 150.0, # seconds
        'AOM_start_buffer' : 50.0, # ns
        'AOM_length' : 1600.0, # ns
        'AOM_light_delay' : 655.0, # ns
        'AOM_end_buffer' : 1155.0, # ns
        'Sacher_AOM_length' : 65000.0, # ns
        'Sacher_AOM_light_delay' : 950.0, # ns
        'Sacher_AOM_end_buffer' : 1155.0, # ns
        'readout_length' : 1000.0, # ns
        'ctr_term' : 'PFI2',
        'CFDLevel0' : 320,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 110,
        'CFDZeroCross1' : 10,
        'Binning' : 7,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : -10000,
        'acquisition_time' : 60.0, # s
        'piezo_start' : 0, #volts
        'piezo_end' : 90, #volts
        'piezo_step_size' : 0.08, # volts (dispersion is roughly ~0.4 GHz/V)
        'PH_trigger_time' : 50.0,
        'PH_trigger_length' : 50.0,
        'bin_size' : 0.05, # GHz, should be same order of magnitude as (step_size * .1 GHz)
        'microwaves' : True, # modulate with microwaves on or off
        'off_resonant_laser' : True, # cycle between resonant and off-resonant
        'power' : 5.0, # dBm
        'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -9.0, # dBm
        'freq' : 1.3569, #GHz
        'dwell_time' : 500.0, # ms
        'temperature_tolerance' : 4.0, # Kelvin
        'MeasCycles' : 30,
        'Imod' : 1.0
        }

p_low = -9
p_high = -9
p_nstep = 1

p_array = np.linspace(p_low,p_high,p_nstep)


for rr in range(np.size(p_array)):
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(1.0)
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
    name_string = '3cdefect'
    m = SiC_Toptica_TCSPC(name_string)
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
msg_string = 'Toptica TCSPC measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
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