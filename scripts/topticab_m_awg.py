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
import statsmodels.api as sm
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

class SiC_Toptica_Search_Piezo_Sweep(m2.Measurement):

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
        if upload or program or clear:
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
    def set_motor_and_measure(self, motor_pos):
        self._motdl.high_precision_move(motor_pos)
        self.check_laser_stabilization()
        wl = self._wvm.get_wavelength()

        if wl != 0.0:
            freq = 299792458.0/wl
        else:
            logging.error(__name__ + ': Bristol wavelength returned 0 -- check alignment?')
            freq = 0.0

        return freq
    def calibrate_laser(self, motor_low = 88000, motor_high = 108000):
        # purpose of this method is to do a regression of the wavelength versus motor step.
        # the usefulness of this is to allow the laser frequency seek to be recalibrated on the fly.
        # the necessity of recalibration was evident when even using large bracket widths of 100 GHz
        # would begin to produce motor ranges that did not contain the desired frequency, which causes
        # the root finding search to fail.
        t0 = time.time()
        print 'Running laser calibration from motor %d to %d.' % (motor_low, motor_high)
        # estimate a step size so that we take about 12 steps across the entire range, max.
        n_desired = 22
        # assume we have about 12 steps, come up with a step size, and then round it to the nearest multiple
        # of 20.
        step_guess = np.round((float(motor_high)-float(motor_low))/(n_desired-1) * 1.0/20.0)*20.0
        # in case the step guess is too small, make sure it's at least 20.
        step_size = np.max([20, step_guess])
        motor_array = np.arange(motor_low, motor_high+step_size, step_size)
        # now, check if we have enough steps -- taking only a few points could lead to a badly behaved calibration
        if np.size(motor_array) < 5:
            logging.warning(__name__ + ': calibrate laser method motor sweep had less than 5 steps. Skipped calibration.')
            return False
        # we know now the calibration has at least 5 points, so let's begin taking a sweep
        y = np.zeros(np.size(motor_array))
        for idx, pos in np.ndenumerate(motor_array):
            self._motdl.high_precision_move(pos)
            self.check_laser_stabilization()
            y[idx] = self._wvm.get_frequency()
            print 'motor %d, freq %.0f' % (pos, y[idx])

        # we have two allocated arrays of data, so regress them with a simple quadratic
        # preprocess the data for better fitting/robustness -- commented out
        motor_array_norm = motor_array # (motor_array-np.mean(motor_array))/np.std(motor_array)
        motor_array_mean = np.mean(motor_array)
        motor_array_std = np.std(motor_array)
        ## quadratic
        #X = np.column_stack((motor_array_norm*motor_array_norm, motor_array_norm, np.ones(np.size(motor_array))))
        # linear
        X = np.column_stack((motor_array_norm, np.ones(np.size(motor_array))))
        y_mean = np.mean(y)
        y_std = np.std(y)
        y_n = (y-y_mean)/y_std

        res = sm.OLS(y, X).fit()
        print(res.summary())
        t1 = time.time()
        self._b = res.params[0]
        self._c = res.params[1]
        #self._c = res.params[2]
        print 'Toptica laser calibration took %.0f seconds. Updated calibration constants to b = %.8f c = %.8f' % (t1-t0, self._b, self._c)
        return
    def laser_frequency_seek(self, frequency):
        bracket_init = 60.0
        f_low = frequency - bracket_init
        f_high = frequency + bracket_init

        ## quadratic
        #m_low, its = self.brent_search(lambda x: self._a*x*x + self._b*x + self._c - f_low, 60000, 120000, 2, 60)
        #m_high, its = self.brent_search(lambda x: self._a*x*x + self._b*x + self._c - f_high, 60000, 120000, 2, 60)

        # linear
        m_low, its = self.brent_search(lambda x: self._b*x + self._c - f_low, 60000, 120000, 2, 60)
        m_high, its = self.brent_search(lambda x: self._b*x + self._c - f_high, 60000, 120000, 2, 60)
        #if m_low > m_high:  ## I think this swapping might not serve a function; should be removed if the program
        ## still works.
        #    m_low, m_high = m_high, m_low
        #freq_discrep = lambda mot_pos: set_motor_and_measure(int(np.round(mot_pos))) - frequency
        #print 'Going to try to get to %.2f GHz...' % dfreq
        m_target, its = self.brent_search((lambda mot_pos: self.set_motor_and_measure(int(np.round(mot_pos))) - frequency), m_low, m_high, 20, 14)
        #found_frequency = self.set_motor_and_measure(m_target)
        found_frequency = 299792458.0/self._wvm.get_wavelength()
        print('f_delta = {:.1f}'.format(found_frequency-frequency))
        return found_frequency

    def flatten_ranges(self, range_list):
        # this function takes a list of frequency ranges and flattens them,
        # i.e. if two ranges overlap, they are merged
        saved = list(range_list[0])
        for st, en in sorted(range_list):
            if st <= saved[1]:
                saved[1] = max(saved[1], en)
            else:
                yield tuple(saved)
                saved[0] = st
                saved[1] = en
        yield tuple(saved)

    def difference_ranges(self, a_range, b_range):
        # subtract the elements one by one
        a_saved = a_range

        for b in b_range:
            a_saved = self.difference_single_element(a_saved, b)

        return a_saved

    def difference_single_element(self, a_list, b):
        true_list = []
        for i in range(len(a_list)):
            # pick out one tuple from a_list and save it
            st, en = list(a_list[i]) # also unpack the tuple into start and end

            if b[0] > st and en > b[0] and en <= b[1]:
                #   |st---------------------en|
                #       |b0---------------------b1|
                true_list.append((st,b[0]))

            elif b[0] <= st and b[1] < en and b[1] > st:
                #     |st-----------------------en|
                # |b0---------------------b1|
                true_list.append((b[1],en))

            elif b[0] <= st and b[1] >= en:
                #     |st-----------------------en|
                # |b0-------------------------------b1|
                # the range a is completely contained within the range of b
                # so we return nothing
                pass
            elif b[0] > st and b[1] < en:
                #     |st-----------------------en|
                #          |b0-------------b1|
                # the range of b is completely contained within a, so we
                # return two new sub-ranges of a that weren't in b
                true_list.append((st, b[0]))
                true_list.append((b[1], en))
            else:
                # in this case, b doesn't overlap at all with a, so
                # we just return a unchanged
                #     |st-----------------------en|
                #                                     |b0------------------b1|
                true_list.append((st,en))
        return true_list

    def brent_search(self, f, a, b, tol, maxiters):
        # This method is a fairly general implementation of Brent's search. The
        # usage goal for it is to take advantage of the relatively smooth nature
        # of the relation between grating angle and laser frequency, and allow
        # the program to hunt for a desired laser frequency.
        fa = f(a)
        fb = f(b)
        if fa*fb >= 0:
            print('Root is not bracketed -- exiting root search. a(%d) = %.2f, b(%d) = %.2f.' % (a, fa, b, fb))
            return (a+b)/2.0, -1
        if (np.abs(a) < np.abs(b)):
            # swap a and b
            a, b = b, a
            fa, fb = fb, fa

        c = a
        fc = fa
        fs = -1.0
        mflag = True
        i = 0

        while (fb != 0.0 and fs != 0.0 and np.abs(b - a) > np.abs(tol)) and i < maxiters:
            if fa != fc and fb != fc:
                # do inverse quadratic interpolation
                s = a*fb*fc/((fa - fb)*(fa - fc)) + b*fa*fc/((fb - fa)*(fb - fc)) + c*fa*fb/((fc - fa)*(fc - fb))
            else:
                # do secant method
                s = b - fb*(b - a)/(fb - fa)
            # check conditions to see if we do the bisection method instead
            tmp = (3*a+b)/4.0
            if (not ((tmp < s and s < b) or (b < s and s < tmp)) or \
                (mflag and np.abs(s-b) >= np.abs(b-c)/2.0) or \
                ((not mflag) and np.abs(s-b) >= np.abs(c-d)/2.0)):
                # do bisection instead
                s = (a+b)/2.0
                mflag = True

            else:
                if ((mflag and np.abs(b-c) < np.abs(tol)) or \
                     ((not mflag) and np.abs(c-d) < np.abs(tol))):
                    # do bisection instead
                    s = (a+b)/2.0
                    mflag = True
                else:
                    mflag = False

            fs = f(s)
            d = c
            c = b
            fc = fb
            if fa*fs < 0.0:
                b = s
                fb = fs
            else:
                a = s
                fa = fs

            if np.abs(fa) < np.abs(fb):
                a, b = b, a
                fa, fb = fb, fa
            i = i+1

        return b, i

    def check_laser_stabilization(self):
        # Monitor laser frequency
        frq_recent = np.zeros(3)
        for zz in range(3):
            time.sleep(1.0)
            frq_recent[zz] = 299792458.0/self._wvm.get_wavelength()
        # Check if laser is stable, if not, wait
        for zz in range(30):
            if (np.max(frq_recent) - np.min(frq_recent)) < 1.0:
                #print 'Laser stabilized. Dispersion %.2f MHz, range %.2f MHz' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                break
            time.sleep(1.0)
            np.roll(frq_recent,1)
            frq_recent[0] = 299792458.0/self._wvm.get_wavelength()
            if zz >=29:
                print 'Laser frequency readout on wavemeter did not stabilize after 30 seconds. Dispersion %.2f MHz, range %.2f MHz.' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                return False

        return True
    def filter_laser_frequencies(self, filter_set, frq_array, total_hits_data, frq1):
        print 'Old filter set %s' % filter_set

        frq_idx = np.nonzero(total_hits_data[:,0])
        threshold = 6.0*self.params['bin_size']
        seen_freqs = []
        if np.any(frq_idx):
            # select which ones we've seen, then scale them back to absolute frequency
            frq_subset = frq_array[frq_idx] + frq1
            print('Frq subset shape {}'.format(frq_subset.shape))
            # this array is already sorted from low to high by construction, so we can search for gaps easily
            frq_diffs = np.diff(frq_subset)
            # now lets segment the already swept frequency range into a list of boundaries
            kk = 0
            lb = frq_subset[0]
            ub = frq_subset[0]
            run_length = 0

            while kk < np.size(frq_subset)-1:
                cur_diff = frq_subset[kk+1] - frq_subset[kk]
                if cur_diff < threshold:
                    ub = frq_subset[kk+1]
                    run_length = run_length + 1
                else:
                    # the skip between points is above threshold, but we need to be sure we saw at least 2 points
                    if run_length > 0:
                        seen_freqs.append((lb, ub))
                        run_length = 0
                        lb = frq_subset[kk+1]
                        print('lb is {}'.format(lb))
                        ub = frq_subset[kk+1]
                        if kk + 3 < np.size(frq_subset):
                            lb = frq_subset[kk+1]
                            #print('lb is {}'.format(lb))
                        else:
                            # there are enough points left to make a new set, so just break out of the loop
                            break
                    else:
                        # we didn't see two points, so don't append to the array

                        lb = frq_subset[kk+1]
                        #print('lb is {}'.format(lb))
                        run_length = 0
                if kk == np.size(frq_subset)-2 and run_length > 0:
                    ub = frq_subset[kk+1]
                    seen_freqs.append((lb, ub))
                kk = kk + 1


        print 'Frequencies seen: %s' % seen_freqs
        # now that we've computed what intervals we've observed, find the difference of the sets, and return it.
        new_filter_set = self.difference_ranges(list(filter_set),list(seen_freqs))

        # finally, check for tuples whose total length is less than 3x the bin size and remove them
        print 'New filter set %s' % new_filter_set
        rem_idx = []
        for k in new_filter_set:
            if (k[1] - k[0]) < 3.0 * self.params['bin_size']:
                rem_idx.append(k)

        for k in rem_idx:
            del new_filter_set[new_filter_set.index(k)]

        print 'New filter (small sets removed) is %s' % new_filter_set

        return new_filter_set


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
    def measure_cw_power(self):
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
        return power_meas
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
        # Set laser calibration constants here
        self._a = 2.424e-7
        self._b = -0.104258
        self._c = 281167.0
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
        elif self._awg.get_state() == 'Running':
            print 'AWG already started -- clearing VISA interface anyway.'
            self._awg.clear_visa()

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
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten
        full_attenuation = self.params['power'] - self.params['constant_attenuation'] - np.max((0,np.min((desired_atten,15.5)))) + np.log(self.params['Imod'])/np.log(10.0)*10.0
        print 'Fully attenuated power is %.2f dBm' % full_attenuation

        # Check if wavemeter returns a valid wavelength
        wl_cur = self._wvm.get_wavelength()
        if wl_cur < 1000.0 or wl_cur > 1170.0:
            logging.error(__name__ + ': wavelength %.2f from wavemeter is not in range')
            print 'Wavelength is %.2f nm' % wl_cur

        return

    def measure(self):
        qt.mstart()
        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        self._stop_measurement = False
        # Wall time
        t0 = time.time()

        data = qt.Data(name='wavemotor_sweep')
        data.add_coordinate('laser frequency (GHz)')
        data.add_coordinate('microwave frequency (GHz)')
        data.add_value('counts')
        data.create_file()

        data2d = qt.Data(name='wavemotor_2d')
        data2d.add_coordinate('laser frequency (GHz)')
        data2d.add_value('counts')




        plot2d_0 = qt.Plot2D(data2d, name='piezoscan_single_sweep', coorddim=0, valdim=1)
        plot3d_0 = qt.Plot3D(data, name='piezoscan_full', coorddims=(0,1), valdim=2, style='image')

        # Populate some arrays
        self.params['piezo_pts'] = np.uint32(1 + np.ceil(np.abs(self.params['piezo_end']-self.params['piezo_start'])/self.params['piezo_step_size']))
        self.params['piezo_array'] = np.linspace(self.params['piezo_start'],self.params['piezo_end'], self.params['piezo_pts'])
        self.params['filter_set'].sort()
        self.params['filter_set'] = list(self.flatten_ranges(self.params['filter_set']))
        self.params['filter_set'].sort()
        self.params['working_set'] = self.params['filter_set']


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



        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)

        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

        scan_on = True
        # Start measurement cycle, so go to proper waveform.
        self._awg.sq_forced_jump(2)
        self.awg_confirm(2)


        self._toptica.set_piezo_voltage(self.params['piezo_array'][0])
        # locate the second frequency to seek to
        self.laser_frequency_seek(self.params['working_set'][-1][1]+8.0)
        self.check_laser_stabilization()
        frq2 = self._wvm.get_frequency() + 100.0 #Ghz
        power_meas = self.measure_cw_power()
        # locate the first frequency to seek to
        self.laser_frequency_seek(self.params['working_set'][0][0]-8.0)
        self.check_laser_stabilization()

        frq1 = self._wvm.get_frequency() - 100.0 #Ghz
        power_meas = self.measure_cw_power()

        print 'Start (reference) frequency %.2f GHz / %.2f nm -- End frequency %.2f GHz / %.2f nm' % (frq1 + 100.0, 299792458.0/(frq1 + 100.0), frq2 - 100.0, 299792458.0/(frq2-100.0))
        self.params['bins'] = np.uint32(1 + np.ceil((100.0+np.absolute(frq2-frq1))/self.params['bin_size']))
        #column 1 of the data set, i.e. relative frequency
        frq_array = np.linspace(-100.0, np.absolute(frq2-frq1), self.params['bins'])
        #column 2, number of counts
        total_count_data = np.zeros((np.size(frq_array),np.size(self.params['freq'])), dtype='uint32')
        #column 3, number of hits in each bin
        total_hits_data = np.zeros((np.size(frq_array),np.size(self.params['freq'])), dtype='uint32')

        for i in range(self.params['MeasCycles']):
            # Enter the loop for measurement
            t1 = time.time()
            seek_iterations = 0
            while np.size(self.params['working_set']) != 0:
                self._toptica.set_piezo_voltage(self.params['piezo_array'][0])
                if self._stop_measurement == True:
                    print 'Measurement aborted.'
                    self._stop_measurement = True
                    break


                target_freq = self.params['working_set'][0][0]-6.0
                remedial_target_freq = target_freq
                self.laser_frequency_seek(target_freq)
                cur_frq = self._wvm.get_frequency()
                ik = 0
                while not (cur_frq > self.params['working_set'][0][0]-20.0 and cur_frq < self.params['working_set'][0][0]) and ik < 6:
                    print 'Current frequency is %.2f GHz, versus desired %.2f GHz. Re-scanning.' % (cur_frq, target_freq)
                    # keep decrementing the target frequency by 1 GHz until we get strictly below but not less than 20 GHz away
                    # dispersion is approximately -0.104 GHz/step, so move by about 10 steps upward to decrease by about 1 GHz.
                    current_motor_position = self._motdl.get_position()
                    self._motdl.high_precision_move(current_motor_position+10)
                    time.sleep(2.0)
                    cur_frq = self._wvm.get_frequency()
                    if ik == 5 and not (cur_frq > self.params['working_set'][0][0]-20.0 and cur_frq < self.params['working_set'][0][0]):
                        # if we've failed to home in by stepping the motor in one direction, just try to seek again and hope it works.
                        self.laser_frequency_seek(self.params['working_set'][0][0]-10.0)
                        cur_frq = self._wvm.get_frequency()
                    ik = ik + 1

                #self.measure_cw_power()
                print 'Laser frequency set to %.2f GHz (target %.2f GHz)' % (cur_frq, self.params['working_set'][0][0])
                # Laser frequency seek already checks the stabilization, so this isn't necessary.
                # self.check_laser_stabilization()

                temp_count_data = np.zeros(self.params['piezo_pts'] , dtype='uint32')

                #sweep through the piezo_array voltages, which should go from low to high and then back to low
                for k in range(np.size(self.params['piezo_array'])):

                    #Set the new piezo voltage
                    self._toptica.set_piezo_voltage(self.params['piezo_array'][k])
                    qt.msleep(0.05)

                    # Measure frequency and counts
                    cur_frq = self._wvm.get_frequency()
                    offset_frq = cur_frq - frq1

                    # Check for superconductivity
                    self._snspd.check()



                    # use filter logic
                    filter_inc = False
                    for i, elem in enumerate(self.params['filter_set']):
                        lo, hi = elem
                        if cur_frq < hi and cur_frq > lo:
                            filter_inc = True
                    # Determine if we should measure in the logic statement here
                    # find all nonzero frequencies
                    #(np.sum(total_hits_data) < 10 or np.min(np.abs( offset_frq - frq_array[np.nonzero(total_hits_data[:,0])])) > self.params['bin_size']
                    if (cur_frq > self.params['working_set'][0][0] and cur_frq < self.params['working_set'][0][1]) and not np.any(self.params['bin_size']*0.51 > np.abs(offset_frq - frq_array[np.nonzero(total_hits_data[:,0])])):
                        cts_array_temp = np.zeros(np.size(self.params['freq']))

                        for idx in range(np.size(self.params['freq'])):
                            self._pxi.set_frequency(self.params['freq'][idx]*1.0e9)
                            time.sleep(0.01)
                            cts = self._ni63.get('ctr1')

                            #find where in the 3 column data structure to add counts
                            index = np.searchsorted(frq_array, offset_frq)

                            #update the appropriate two columns keeping track of total counts and total hits
                            total_count_data[index,idx] = total_count_data[index,idx] + cts # temp_count_data[j]
                            total_hits_data[index,idx] = total_hits_data[index,idx] + 1

                            # live plot -- unpack cts_list into args of add_data_point
                            data.add_data_point(cur_frq,self.params['freq'][idx],cts)
                            if idx == 0:
                                data2d.add_data_point(cur_frq, cts)
                                plot2d_0.update()

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

                            qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.

                        data.new_block()
                        plot3d_0.update()
                        qt.msleep(0.002)
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

                    # determine if we are beyond the endpoint of the current desired sweep range
                    if cur_frq > 1.0 + self.params['working_set'][0][1]:
                        print 'Current frequency %.1f too far beyond end of working set (%.2f, %.2f); breaking' % (cur_frq, self.params['working_set'][0][0], self.params['working_set'][0][1])
                        break
                self.params['working_set'] = self.filter_laser_frequencies(self.params['working_set'], frq_array, total_hits_data, frq1)

            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1

            print 'Cycle %d/%d total time is %.3f. Heater output is at %.1f. ' % (i+1, int(self.params['MeasCycles']), tt, self._ls332.get_heater_output())

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
            seek_iterations += 1



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

        qt.mend()
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
        'Sacher_AOM_length' : 1000.0, # ns
        'Sacher_AOM_light_delay' : 960.0, # ns
        'Sacher_AOM_end_buffer' : 1155.0, # ns
        'RF_start_buffer' : 300.0, # ns
        'readout_length' : 1000.0, # ns
        'readout_buffer' : 10.0, # ns
        'ctr_term' : 'PFI2',
        'piezo_start' : 0, #volts
        'piezo_end' : 90, #volts
        'piezo_step_size' : 0.12, # volts (dispersion is roughly ~0.4 GHz/V)
        'bin_size' : 0.12, # GHz, should be same order of magnitude as (step_size * .1 GHz)
        'microwaves' : False, # modulate with microwaves on or off
        'microwaves_CW' : True, # are the microwaves CW? i.e. ignore pi pulse length
        'pi_length' : 180.0, # ns
        'off_resonant_laser' : False, # cycle between resonant and off-resonant
        'power' : 5.0, # dBm
        'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -9.0, # dBm
        'freq' : [1.45], #GHz
        'dwell_time' : 10000.0, # ms
        #'filter_set' : ( (270850, 270870), (270950, 270970)),(270810, 270940),
        'filter_set' : [(270948,271030)],
        'temperature_tolerance' : 12.0, # Kelvin
        'MeasCycles' : 1,
        'Imod' : 0.5,
        }
def main():
    #p_low = -19
    #p_high = -19
    #p_nstep = 1

    #p_array = np.linspace(p_low,p_high,p_nstep)

    # Create a measurement object m

    time.sleep(2.0)
    if msvcrt.kbhit():
        kb_char=msvcrt.getch()
        if kb_char == "q":
            do_track = False

    name_string = 'randomdefect'
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
            return
    print 'Proceeding with measurement ...'
    #m.calibrate_laser()
    m.measure()
    m.save_params()
    m.save_stack()
    m.finish()


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