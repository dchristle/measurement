import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)
import gc

class SiC_Hahn_Master(m2.Measurement):

    mprefix = 'hahn'

    def sequence(self, upload=True, program=True):
        gc.collect()

        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')
        sq_delay = pulse.SquarePulse(channel='MW_pulsemod', name='delay',
                    length = 200e-9, amplitude = 0.)
        # Here is a manual insertion of a tau delay array that I made by hand for this measurement
        #self.params['tau_delay'] = np.array((3500, 5000, 10000, 32000, 72500, 80500, 90000, 160000, 205000, 287000, 345000, 356000, 365500, 380000, 395000),dtype=np.float64)
        #self.params['tau_delay'] = np.array((3500, 5000, 8000, 347500, 351400, 355360, 359290, 363210, 367100, 371070, 375000),dtype=np.float64)
        #self.params['tau_delay'] = np.array((3500, 5000, 8000, 694000, 710000, 726000, 734000),dtype=np.float64)


        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['tau_length_end'] - self.params['tau_length_start'])/self.params['tau_length_step']))
        self.params['tau_delay'] = np.linspace(self.params['tau_length_start'], self.params['tau_length_end'], n_steps)

        self.params['pts'] = np.uint32(self.params['tau_delay'].size)
        self.params['trigger_periods'] = np.zeros(int(self.params['pts']*2))
        self._awg = qt.instruments['awg']
        self._awg.stop()
        time.sleep(1)
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
        elements.append(e)

        # Now create the Hahn pulses
        for i in range(self.params['pts']):
            # we want to calculate the minimum physical duty cycle, and then take the maximum between the minimum physical
            # and the minimum allowable (by a parameter in self.params) duty cycle.
            for j in range(2):
                # precalculate some of the times
                total_rf_length = self.params['RF_delay'] + self.params['pi_length'] + self.params['pi2_length']*2 + self.params['tau_delay'][i] # just for the RF to play out completely
                total_rf_pulses = total_rf_length + self.params['RF_buffer']
                AOM_start_time = total_rf_pulses - self.params['AOM_light_delay']
                readout_start_time = AOM_start_time + self.params['AOM_light_delay']
                #print '%s %s %s %s' % (total_rf_length, total_rf_pulses, AOM_start_time, readout_start_time)
                if j == 0:
                    # this is a positive pulse
                    e = element.Element('ElectronHahn_pos_pt-%d' % i, pulsar=qt.pulsar)
                    IQ_pImod_length = self.params['RF_delay'] + self.params['pi2_length']*2 + self.params['pi_length'] + self.params['tau_delay'][i] + self.params['RF_buffer'] + self.params['AOM_length'] + self.params['AOM_light_delay']
                    IQ_nImod_length = 0.0
                if j == 1:
                    # this is a negative pulse
                    e = element.Element('ElectronHahn_neg_pt-%d' % i, pulsar=qt.pulsar)
                    # negative phase final pulse IQ sequence
                    IQ_pImod_length = self.params['RF_delay'] + self.params['tau_delay'][i]/2.0 + self.params['pi2_length'] + self.params['pi_length'] + 0.5*self.params['tau_delay'][i]/2.0
                    IQ_nImod_length = 0.5*self.params['tau_delay'][i]/2.0 + self.params['pi2_length'] + self.params['RF_buffer'] + self.params['AOM_length'] + self.params['AOM_light_delay']


                e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_length']*1.0e-9), name='laser init', start=AOM_start_time*1.0e-9)
                # compute the explicit start times of the first pi/2 pulse, the pi pulse, and the second pi/2 pulse
                fpi2_start_time = self.params['RF_delay']
                pi_start_time = fpi2_start_time + self.params['pi2_length'] + self.params['tau_delay'][i]/2.0
                spi2_start_time = pi_start_time + self.params['pi_length']  + self.params['tau_delay'][i]/2.0
                #print '%s %s %s' % (fpi2_start_time, pi_start_time, spi2_start_time)
                e.add(pulse.cp(sq_pulseMW, length = self.params['pi2_length']*1.0e-9, amplitude = 1.0), name='first pi2 microwave pulse', start=fpi2_start_time*1.0e-9)

                e.add(pulse.cp(sq_pulseMW, length = self.params['pi_length']*1.0e-9, amplitude = 1.0), name='pi pulse', start=pi_start_time*1.0e-9)

                e.add(pulse.cp(sq_pulseMW, length = self.params['pi2_length']*1.0e-9, amplitude = 1.0), name='second pi2 microwave pulse', start=spi2_start_time*1.0e-9)

                e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=self.params['readout_length']*1.0e-9),
                name='photoncountpulse', start=readout_start_time*1.0e-9)

                e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=IQ_pImod_length*1.0e-9),
                name='MWimodpulse_plus', start=0e-9)
                e.add(pulse.cp(sq_pulseMW_Imod, amplitude=-1.0, length=IQ_nImod_length*1.0e-9),
                name='MWimodpulse_negative', start=IQ_pImod_length*1.0e-9)

                e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=(IQ_pImod_length + IQ_nImod_length)*1.0e-9),
                name='MWqmodpulse', start=0e-9)

                elements.append(e)
                # Now set the total pulse lengths into the array so we can adjust the dwell times so that each point gets the same number of cycles per dwell time

                self.params['trigger_periods'][2*i+j] = (IQ_pImod_length + IQ_nImod_length)





        print 'Trigger period array: %s' % self.params['trigger_periods']
        # create a sequence from the pulses
        seq = pulsar.Sequence('Electron Hahn Echo sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)
        gc.collect()
        if upload:
            qt.pulsar.upload(*elements)
        time.sleep(3.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)
            print 'AWG sequencing complete. Waiting for 3 s...'
            time.sleep(3.0)
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
        for i in range(20):
            try:
                if self._awg.get_state() == 'Idle':
                    self._awg.start()
                break
            except(visa.visa.VI_ERROR_TMO):
                print 'AWG still busy -- trying again...'

        # set the AWG to CW mode
        print 'Waiting 15 s for AWG to start...'
        time.sleep(10.0)

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
        if self.params['tau_length_end'] < self.params['tau_length_start']:
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
        # Wall time
        t0 = time.time()

        # Populate some arrays

        print '--Hahn echo meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps), %.5f GHz --' % (self.params['tau_length_start'], self.params['tau_length_end'], self.params['tau_length_step'], self.params['pts'], self.params['freq'])

        total_count_data = np.zeros(int(self.params['pts']*2.0), dtype='uint32')
        average_count_data = np.zeros(int(self.params['pts']*2.0), dtype='float')
        intermediate_total_data = np.zeros( (1,int(self.params['pts']*2.0)), dtype='uint32')
        signal = np.zeros(1)


        # Set the PXI status to 'on', i.e. generate microwaves
        self._pxi.set_status('on')
        # Start the AWG sequencing
        self._awg.start()
        time.sleep(3)
        print 'Arbitrary waveform generator booted up.'
        N_cmeas = 0
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process
        # Optimize

        # Get the current sequence position index
        prev_awg_sq_position = int(self._awg.get_sq_position())
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
            print 'Starting iteration %d of %d' % ((i+1),self.params['MeasCycles'])

            # Create an index of the waveforms so that we can modify it
            seq_index = range(2*self.params['pts'])
            if self.params['random'] == 1:
                # Now shuffle the array in place
                shuffle(seq_index)
            # Create an array for the single-sweep data
            temp_count_data = np.zeros(int(2*self.params['pts']), dtype='uint32')


            # Enter the loop for measurement
            t1 = time.time()
            for j in range(int(2*self.params['pts'])):

                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:

                    # Now set the AWG into CW mode for tracking
                    self._awg.sq_forced_jump(1)
                    self.awg_confirm(1)
##                    if self._snspd.check() == False:
##                        print 'SNSPD went normal and could not restore, breaking.'
##                        break
                    time.sleep(0.1)
                    # Re-optimize
                    fbl.optimize()


                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()


                # Set the new RF pulse waveform
                self._awg.sq_forced_jump(seq_index[j]+2) # the +2 is because the indices start at 1, and the first sequence is CW mode
                time.sleep(0.01)
                if random.random() < 0.05:
                    # 5% of the time, check to see if the AWG is functioning properly -- this should
                    # still allow us to pick up any sequence jumping errors without causing so much overhead from constant checking
                    self.awg_confirm(seq_index[j]+2)
                # adjust the dwell time to compensate for longer or shorter duty cycles, relative to the shortest
                self._ni63.set_count_time(float(self.params['dwell_time'])/1000.0*float(self.params['trigger_periods'][seq_index[j]])/float(self.params['trigger_periods'][0]))
##                if self._snspd.check() == False:
##                    print 'SNSPD went normal and could not restore, breaking.'
##                    break

                temp_count_data[j] = self._ni63.get('ctr1')
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1
            total_dwell_time = self.params['dwell_time']/1000.0 * np.sum(self.params['trigger_periods'])/self.params['trigger_periods'][0]
            print 'Cycle %d/%d: Total time is %.3f, efficiency of %.2f percent. Heater output is at %.1f.' % (i+1, int(self.params['MeasCycles']), tt, total_dwell_time/tt*100.0, self._ls332.get_heater_output())

            # Converting to numpy and back is sort of a hack, but it's one line.
            sorted_temp_data = temp_count_data[np.array(seq_index).argsort().tolist()]
            #print 'sorted temp data %s' % sorted_temp_data
            plot2d_0 = qt.Plot2D(np.array(range(int(self.params['pts']*2.0))).astype(np.float64),sorted_temp_data.astype(np.float64), name='hahn_single_sweep', clear=True)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False: break
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
            try:
                subtracted_array = np.zeros(self.params['pts'],dtype=np.float64)
                for kk in range(self.params['pts']):
                    subtracted_array[kk] = total_count_data.astype(np.float64)[2*kk+1] - total_count_data.astype(np.float64)[2*kk]
                plot2d_3 = qt.Plot2D(np.array(self.params['tau_delay']).astype(np.float64),subtracted_array.astype(np.float64), name='hahn_subtracted', clear=True)
                if i == 0:
                    intermediate_total_data[0,:] = total_count_data.astype(np.float64)
                elif np.mod(i,10) == 0:
                    intermediate_temp_data = np.zeros( (i/10+1,int(2*self.params['pts'])), dtype=np.float64)
                    intermediate_temp_data[:-1,:] = intermediate_total_data
                    intermediate_temp_data[i/10,:] = total_count_data.astype(np.float64)
                    intermediate_total_data = np.copy(intermediate_temp_data)
                    #print 'size is %s' % (np.size(intermediate_total_data))
                    #intermediate_total_data = np.vstack((intermediate_total_data,total_count_data))
                if i == 0:
                    signal[0] = self._ni63.get('ctr1')
                elif np.mod(i,10):
                    signal = np.hstack((signal,self._ni63.get('ctr1')))
                #print 'range array is %s' % (np.array(range(int(self.params['pts']*2.0))))
                plot2d_1 = qt.Plot2D(np.array(range(int(self.params['pts']*2.0))).astype(np.float64),total_count_data.astype(np.float64), name='hahn_avg', clear=True)
                N_cmeas = N_cmeas + 1
                average_count_data = total_count_data/float(N_cmeas)
            except Exception as e:
                print 'Encountered exception %s in after-measurement data processing, passing...' % e
                pass



        # Stop PXI sig gen
        self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Hahn_data', self.h5data, base=self.h5base)
        grp.add('length', data=self.params['tau_delay'], unit='ns', note='frequency')
        grp.add('trigger_periods', data=self.params['trigger_periods'], unit='ns',note='trigger period lengths')
        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')
        grp.add('intermediate', data=intermediate_total_data, unit='', note='intermediate total count data')
        grp.add('signal', data=signal, unit='counts', note='signal rate per N iterations')


        return



# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 80.0, # seconds
        'AOM_length' : 2000.0, # ns
        'AOM_light_delay' : 655.0, # ns
        'RF_delay' : 450.0, # ns
        'RF_buffer' : 400.0, # ns
        'readout_length' : 210.0, # ns
        'ctr_term' : 'PFI2',
        'power' : 5.0, # dBm
        'constant_attenuation' : 3.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -16.0, # dBm
        'tau_length_start' : 0.0, # ns
        'tau_length_end' : 2000.0, # ns
        'tau_length_step' : 7.0, # ns
        'minimum_duty_cycle' : 3600.0, # ns
        'freq' : 1.23194, #GHz
        'pi2_length' : 18.0, # ns
        'pi_length' : 36.0, # ns
        'dwell_time' : 1000.0, # ms
        'temperature_tolerance' : 2.0, # Kelvin
        'MeasCycles' : 220,
        'random' : 1
        }

p_array = np.array([-0.0])

for rr in range(np.size(p_array)):
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
    name_string = 'power %.2f dBm' % (p_array[rr])
    m = SiC_Hahn_Master(name_string)
    xsettings['readout_length'] = 220.0
    xsettings['desired_power'] = p_array[rr]
    # since params is not just a dictionary, it's easy to incrementally load
    # parameters from multiple dictionaries
    # this could be very helpful to load various sets of settings from a global
    # configuration manager!
    m.params.from_dict(xsettings)
    m.sequence(upload=False,program=False)


    if False:
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
#ea_t = qt.instruments['ea']
#ls332_t = qt.instruments['ls332']
#cur_temp = ls332_t.get_kelvinA()
#msg_string = 'Hahn measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
#ea_t.email_alert(msg_string)

##ps = qt.instruments['xps']
##ps.set_abs_positionZ(12.0)
#
#track_on = True
#fbl_t = qt.instruments['fbl']
#track_iter = 0
#while track_on == True and track_iter < 50:
#    track_iter = track_iter + 1
#    print 'Tracking for %d iteration.' % track_iter
#    fbl_t.optimize()
#    time.sleep(5.0)
#    if msvcrt.kbhit() or track_on == False:
#                kb_char=msvcrt.getch()
#                if kb_char == "q" or track_on == False: break