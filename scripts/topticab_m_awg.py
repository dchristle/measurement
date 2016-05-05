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
import progressbar
import random
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
        self.find_stable_frequency()
        self.check_laser_stabilization()
        wl = self._wvm.get_wavelength()

        if wl != 0.0:
            freq = 299792458.0/wl
        else:
            logging.error(__name__ + ': Bristol wavelength returned 0 -- check alignment?')
            freq = 0.0

        return freq
    def calibrate_laser(self, motor_low = 145000, motor_high = 165000):
        # purpose of this method is to do a regression of the wavelength versus motor step.
        # the usefulness of this is to allow the laser frequency seek to be recalibrated on the fly.
        # the necessity of recalibration was evident when even using large bracket widths of 100 GHz
        # would begin to produce motor ranges that did not contain the desired frequency, which causes
        # the root finding search to fail.
        t0 = time.time()
        print 'Running laser calibration from motor %d to %d.' % (motor_low, motor_high)
        # estimate a step size so that we take about 12 steps across the entire range, max.
        n_desired = 23
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
        motor_freq_list = []
        for idx in range(np.size(motor_array)):
            pos = motor_array[idx]
            self._motdl.high_precision_move(pos)
            fp_out = self._fp.check_stabilization()
            if fp_out == 1:
                motor_freq_list.append((motor_array[idx], self._wvm.get_frequency()))

        motor_array_np = np.array(motor_freq_list)
        motor_array = motor_array_np[:,0]
        y = motor_array_np[:,1]
        # we have two allocated arrays of data, so regress them with a simple quadratic
        # shifting and scaling the motor positions to have a mean of 0 and std = 1 eliminates problems
        # with the condition number of the OLS regression being large, i.e. it's more numerically stable

        motor_array_norm = (motor_array-np.mean(motor_array))/np.std(motor_array)
        self._motor_array_mean = np.mean(motor_array)
        self._motor_array_std = np.std(motor_array)
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
        print 'Toptica laser calibration took %.0f seconds. Updated calibration constants to b = %.8f c = %.8f, motor_array_mean = %.2f, motor_array_std = %.2f' % (t1-t0, self._b, self._c, self._motor_array_mean, self._motor_array_std)
        return

    def coarse_seek_frequency(self, frequency, step_low, step_high, motor_step = 50, tolerance = 0):
        fp = qt.instruments['fp']
        motdl = qt.instruments['motdl']
        wvm = qt.instruments['bristol']
        # this is a replacement for the root-finding based algorithm for coarse frequency tuning, e.g. to within
        # a few GHz. The idea here is to replace that algorithm with a simpler 1D line search that steeps from
        # low to high and then moves to the frequency closest to the desired one (stored in frequency).
        #
        print 'Seeking to %.2f GHz with 1D line search...' % frequency
        init_motor_position = motdl.get_position()
        low_step = np.min(np.array((step_low, step_high)))
        high_step = np.max(np.array((step_low, step_high)))
        print 'low step %d high step %d' % (low_step, high_step)
        motor_array = np.arange(low_step, high_step, np.abs(np.round(motor_step)))
        motor_freq_list = []
        with progressbar.ProgressBar(max_value=np.size(motor_array)) as bar:
            for i in range(np.size(motor_array)):
                mot_pos = motor_array[i]
                bar.update(i)
                motdl.high_precision_move(mot_pos)
                #print 'mot pos %d' % mot_pos
                time.sleep(0.5)
                fp_out = fp.check_stabilization()
                if fp_out == 1:
                    # laser is stable at this motor position, so get the frequency
                    motor_freq_list.append((mot_pos, wvm.get_frequency()))
                    # check if we are within tolerance -- this can help stop the search faster
                    if np.abs(motor_freq_list[-1][1] - frequency) < tolerance:
                        # we are within the allowed tolerance, so return
                        bar.update(np.size(motor_array)-1)
                        final_frequency = wvm.get_frequency()
                        print 'Distance is now %.2f GHz away from %.2f GHz.' % (final_frequency - frequency, frequency)
                        return True

        if not motor_freq_list:
            logging.warning(__name__ + '1D line search failed -- no frequencies were stable. Returning to original frequency.')
            motor_array = np.arange(step_low-1000, step_high+1000, np.round(motor_step)*2)
            motor_freq_list = []
            with progressbar.ProgressBar(max_value=np.size(motor_array)) as bar:
                for i in range(np.size(motor_array)):
                    mot_pos = motor_array[i]
                    bar.update(i)
                    motdl.high_precision_move(mot_pos)
                    time.sleep(0.5)
                    fp_out = fp.check_stabilization()
                    if fp_out == 1:
                        # laser is stable at this motor position, so get the frequency
                        motor_freq_list.append((mot_pos, wvm.get_frequency()))
                        #print 'Pos: %d , Freq: %.2f GHz' % (motor_freq_list[-1][0], motor_freq_list[-1][1])
            if not motor_freq_list:
                logging.warning(__name__ + '1D line search failed again -- no frequencies were stable. Returning to original frequency.')
                motdl.high_precision_move(init_motor_position)
                return False
        # find the motor position nearest to the desired frequency
        nearest_pos = motor_freq_list[0][0]
        nearest_dist = motor_freq_list[0][1] - frequency
        for elem in motor_freq_list:
            if np.abs(elem[1] - frequency) < np.abs(nearest_dist):
                nearest_pos = elem[0]
                nearest_dist = elem[1] - frequency
        #print 'Found motor step %d to be %.2f GHz away -- seeking.' % (nearest_pos, nearest_dist)
        motdl.high_precision_move(nearest_pos)
        final_frequency = wvm.get_frequency()
        print 'Distance is now %.2f GHz away from %.2f GHz.' % (final_frequency - frequency, frequency)
        if np.abs(final_frequency - frequency) > 15.0:
            print(motor_freq_list)
        return True
    def check_if_reference_is_needed(self):
        current_frequency = self._wvm.get_frequency()
        current_motor = self._motdl.get_position()

        motor_array = current_motor + np.arange(100,1100,200)
        masize = np.size(motor_array)
        freq_list = []
        with progressbar.ProgressBar(max_value=masize) as bar:
            for idx in range(masize):
                bar.update(idx)
                self._motdl.high_precision_move(motor_array[idx])
                time.sleep(0.25)
                fp_out = self._fp.check_stabilization()
                if fp_out == 1:
                    freq_list.append((motor_array[idx], self._wvm.get_frequency()))
                    # find the min and max
                    min = freq_list[0][1]
                    max = freq_list[0][1]
                    for elem in freq_list:
                        if elem[1] < min:
                            min = elem[1]
                        if elem[1] > max:
                            max = elem[1]
                    if (max-min) > 15.0:
                        # we've clearly moved, so we can stop taking data
                        self._motdl.high_precision_move(current_motor)
                        bar.update(masize-1)
                        return

        if not freq_list:
            # no frequencies were seen with stable modes - might indicate that the motor is in a multimode
            # position and not moving
            print 'Reference search needed -- did not observe any single mode behavior, so motor might be stuck at a multimode position.'
            self._motdl.reference_search()
            self._motdl.high_precision_move(current_motor)
            return
        # we saw some stable frequencies, but need to determine if they're difference enough to say for sure that
        # the motor still functions properly

        # find the min and max
        min = freq_list[0][1]
        max = freq_list[0][1]
        for elem in freq_list:
            if elem[1] < min:
                min = elem[1]
            if elem[1] > max:
                max = elem[1]
        # we found the max and min elements of the frequencies we saw, so if this is not greater than 15 GHz,
        # we need a reference search
        if max-min < 15.0:
            # need a reference search
            print 'Reference search running - dispersion was %.2f GHz (less than 15 GHz). ' % (max-min)
            print 'Array list:'
            print(freq_list)
            self._motdl.reference_search()

        # presumably, if we have gotten to this point, the motor should be OK and respond to the command to move
        # it back to its original position
        self._motdl.high_precision_move(current_motor)
        return
    def remove_motor_backlash(self, frequency):
        initial_motor_position = self._motdl.get_position()
        current_frequency = self._wvm.get_frequency()
        initial_discrepancy = np.abs(current_frequency-frequency)


        # high precision move is done here, which means that we move 10000 steps left or right and then back again
        # in order to bring the backlash out of the motor


        if (current_frequency-frequency) > 0.0:
            # we are above the desired frequency, so we want to move positive steps. before doing that, we will move
            # 10000 steps negative, and then 10000 steps positive.
            self._motdl.high_precision_move(initial_motor_position,1) # the 1 specifies the final move is increasing
            # re-measure the discrepancy
            its = 0
            fp_out = self._fp.check_stabilization()
            current_frequency = self._wvm.get_frequency()
            while fp_out != 1 and its < 9:
                fp_out = self._fp.check_stabilization()
                if fp_out == 1:
                    current_frequency = self._wvm.get_frequency()
                else:
                    self._motdl.set_position(self._motdl.get_position() + 100)
                if its == 8:
                    logging.warning(__name__ + ': motor did not find a stable frequency after moving. proceeding anyway')
                its = its + 1

            current_frequency = self._wvm.get_frequency()
            if (current_frequency-frequency) < 0.0:
                # the sign of the needed movement has changed. now operate in reverse to try to get the backlash out
                self._motdl.high_precision_move(initial_motor_position,0) # the 0 specifies the final move is decreasing
                # re-measure the discrepancy
                its = 0
                fp_out = self._fp.check_stabilization()
                current_frequency = self._wvm.get_frequency()
                while fp_out != 1 and its < 9:
                    fp_out = self._fp.check_stabilization()
                    if fp_out == 1:
                        current_frequency = self._wvm.get_frequency()
                    else:
                        self._motdl.set_position(self._motdl.get_position() - 100)
                    if its == 8:
                        logging.warning(__name__ + ': motor did not find a stable frequency after moving. proceeding anyway')
                    its = its + 1
                if its == 9:
                    # didn't find a stable frequency, attempt a frequency measurement anyway
                    current_frequency = self._wvm.get_frequency()
                    return False
                # at this point, we likely have an OK measure of frequency, and have tried our best to remove the
                # backlash. if it failed, we just continue with the standard frequency seeking method.
                return True
            else:
                # we still need to move in the positive step direction, and we've removed the backlash, so we're done
                return True
        else:
            # we are below the desired frequency, so we want to move negative steps. before doing that, we will move
            # 10000 steps positive, and then 10000 steps negative.
            self._motdl.high_precision_move(initial_motor_position,0) # the 0 specifies the final move is decreasing
            # re-measure the discrepancy
            its = 0
            fp_out = self._fp.check_stabilization()
            current_frequency = self._wvm.get_frequency()
            while fp_out != 1 and its < 9:
                fp_out = self._fp.check_stabilization()
                if fp_out == 1:
                    current_frequency = self._wvm.get_frequency()
                else:
                    self._motdl.set_position(self._motdl.get_position() - 100)
                if its == 8:
                    logging.warning(__name__ + ': motor did not find a stable frequency after moving. proceeding anyway')
                its = its + 1

            current_frequency = self._wvm.get_frequency()
            if (current_frequency-frequency) > 0.0:
                # the sign of the needed movement has changed. now operate in reverse to try to get the backlash out
                self._motdl.high_precision_move(initial_motor_position,1) # the 1 specifies the final move is increasing
                # re-measure the discrepancy
                its = 0
                fp_out = self._fp.check_stabilization()
                current_frequency = self._wvm.get_frequency()
                while fp_out != 1 and its < 9:
                    fp_out = self._fp.check_stabilization()
                    if fp_out == 1:
                        current_frequency = self._wvm.get_frequency()
                    else:
                        self._motdl.set_position(self._motdl.get_position() + 100)
                    if its == 8:
                        logging.warning(__name__ + ': motor did not find a stable frequency after moving. proceeding anyway')
                    its = its + 1
                if its == 9:
                    # didn't find a stable frequency, attempt a frequency measurement anyway
                    current_frequency = self._wvm.get_frequency()
                    return False
                # at this point, we likely have an OK measure of frequency, and have tried our best to remove the
                # backlash. if it failed, we just continue with the standard frequency seeking method.
                return True
            else:
                # we still need to move in the positive step direction, and we've removed the backlash, so we're done
                return True
    def laser_frequency_seek(self, frequency):
        # this seek is based on a simpler idea of just using the set position only. it suffers from backlash of the
        # motor but should eventually get quite close to the correct setting and fairly quickly

        threshold = 3.0
        self._toptica.set_current(self.params['current_range'][0])
        time.sleep(0.2)
        self.find_single_mode()
        current_frequency = self._wvm.get_frequency()
        if np.abs(current_frequency-frequency) < threshold:
            print 'Frequency is already within threshold -- current %.2f GHz/target %.2f GHz' % (current_frequency,frequency)
            return True



        # try to remove the backlash
        if not self.remove_motor_backlash(frequency):
            # if the first backlash direction doesn't succeed, try again once more
            self.remove_motor_backlash(frequency)

        self.find_single_mode()
        initial_motor_position = self._motdl.get_position()
        current_frequency = self._wvm.get_frequency()
        initial_discrepancy = np.abs(current_frequency-frequency)

        print 'Motor pos: %d, FP %d, GHz: %.2f, cur %.2f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency(), self._toptica.get_current())
        # setup a progress bar based on the relative distance from the start
        with progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ],max_value=1.0) as bar:
            its = 0

            while np.abs(current_frequency-frequency) > threshold and its < 100:
                delta = current_frequency - frequency

                # determine an appropriate motor step size based on the distance away
                if np.abs(delta) > 100.0:
                    step_size = 1000 + random.randint(0,35)
                elif 50 < np.abs(delta) <= 100:
                    step_size = 450 + random.randint(0,40)
                elif 10 < np.abs(delta) <= 50:
                    step_size = 200 + random.randint(0,20)
                else:
                    step_size = 50 - 5 + random.randint(0,10) # adds an element of jitter to the tweaking process

                current_motor_position = self._motdl.get_position()
                if delta >= 0:
                    # current frequency is too high, so we should increase the motor position
                    dm = 1.0*step_size
                else:
                    # current frequenzy is too low, so we should decrease the motor position
                    dm = -1.0*step_size

                if 170000 < (current_motor_position+dm) or current_motor_position+dm < 80000:
                    logging.warning(__name__ + ': Motor tried to move beyond boundaries -- breaking')
                    return False

                self._motdl.set_position(current_motor_position+dm)
                print 'Motor pos: %d, FP %d, GHz: %.2f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency())
                # if we are close enough (within 10 GHz), we want to start sweeping the current to see if we can find
                # a stable point
                # if np.abs(delta) < 10.0:
                #     cur_sweep_array = []
                #     for i, cur in np.ndenumerate(np.arange(self.params['current_range'][0],self.params['current_range'][0],0.001)):
                #         self._toptica.set_current(cur)
                #         time.sleep(0.25)
                #         print 'Motor pos: %d, FP %d, GHz: %.2f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency())
                #         if self._fp.check_stabilization() == 1:
                #             cur_sweep_array.append((current_motor_position+dm, cur, self._wvm.get_frequency()))
                #     # once we've swept the current, let's iterate through and find the closest value
                #     if len(cur_sweep_array) > 0:
                #         # search the array for any combinations that produced a frequency within the threshold
                #         min_dist = 100.0
                #         min_dist_cur = (self.params['current_range'][0]+self.params['current_range'][1])/2.0
                #         for elem in cur_sweep_array:
                #             if np.abs(elem[2]-frequency) < min_dist:
                #                 min_dist = np.abs(elem[2]-frequency)
                #                 min_dist_cur = elem[1]
                #         # if the min_dist satisfies the threshold, then set the current there, verify it, and break
                #         # from the loop
                #         if min_dist < threshold:
                #             self._toptica.set_current(min_dist_cur)
                #             print 'Motor pos: %d, FP %d, GHz: %.2f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency())
                #             time.sleep(1.0)
                #
                #             if fp.check_stabilization() == 1 and np.abs(self._wvm.get_frequency()-frequency) < threshold:
                #                 # yes, it's a good frequency and within the threshold, so break
                #                 current_frequency = self._wvm.get_frequency()
                #                 bar.update(np.clip(np.abs(1.0-np.abs((current_frequency-frequency)/(initial_discrepancy))),0,1))
                #                 break
                #         else:
                #             # min_dist didn't satisfy the threshold, so go back to the center current
                #             self._toptica.set_current((self.params['current_range'][0]+self.params['current_range'][1])/2.0)
                #
                #     else:
                #         # failed to find any stable points -- set the current back to the middle and proceed with the
                #         # standard routine.
                #         self._toptica.set_current((self.params['current_range'][0]+self.params['current_range'][1])/2.0)

                time.sleep(0.25)
                fp_out = self._fp.check_stabilization()
                if fp_out == 1:
                    current_frequency = self._wvm.get_frequency()
                    bar.update(np.clip(np.abs(1.0-np.abs((current_frequency-frequency)/(initial_discrepancy))),0,1))
                    #print 'Motor set to %d, frequency is now %.2f GHz' % (current_motor_position+dm, current_frequency)
                elif fp_out == 2:
                    # do a quick current sweep
                    for cur_idx, cur in np.ndenumerate(np.arange(self.params['current_range'][0],self.params['current_range'][1],0.0005)):
                        self._toptica.set_current(cur)
                        time.sleep(0.25)
                        fp_out = self._fp.check_stabilization()
                        if fp_out == 1:
                            current_frequency = self._wvm.get_frequency()
                            bar.update(np.clip(np.abs(1.0-np.abs((current_frequency-frequency)/(initial_discrepancy))),0,1))
                            print 'Motor pos: %d, FP %d, GHz: %.2f, cur %.2f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency(), self._toptica.get_current())
                            break
                    pass
                    #print 'Motor set to %d, but FP was not single mode.' % (current_motor_position+dm)
                elif fp_out == 3:
                    logging.error(__name__ + ': no signal retrieved from FP, cannot scan wavelength. Breaking.')
                    return False

                its = its + 1

        return True
    def laser_frequency_seek_new(self, frequency):
        # determine the low and high motor positions
        bracket_init = 120.0
        f_low = frequency - bracket_init
        f_high = frequency + bracket_init

        # linear
        m_low, its = self.brent_search(lambda x: self._b*(x - self._motor_array_mean)/self._motor_array_std + self._c - f_low, 120000, 190000, 0.1, 60)
        m_high, its = self.brent_search(lambda x: self._b*(x - self._motor_array_mean)/self._motor_array_std + self._c - f_high, 120000, 190000, 0.1, 60)
        if its == -1:
            print logging.error(__name__ + ': Calibration constants are off by too much to find proper motor boundaries.')
            return False
        # check if a reference of the motor is needed. this can happen because the motor occasionally stops responding
        self.check_if_reference_is_needed()

        # create an array of motor points to sweep - maybe 10 -- with nice round numbers for the motor positions
        # first round to the nearest 100
        mceil_high = np.ceil(m_high/100.0)*100
        mceil_low = np.floor(m_low/100.0)*100


        # create a numpy array with scan points that are fairly round numbers
        scan_step = 100
        print 'Scan step: %d' % scan_step
        # use the coarse frequency seek to scan
        final_frequency_coarse = self.coarse_seek_frequency(frequency, np.min(np.array((mceil_low,mceil_high))), np.max(np.array((mceil_low,mceil_high)))+100, scan_step, tolerance = 0)
        coarse_motor_position = self._motdl.get_position()
        if final_frequency_coarse:
            # coarse step successful, doing fine step
            #print 'Fine scan step: %d' % (np.round(scan_step/2.0/10.0)*10)
            # estimate based on a linear dispersion measure how many steps we need to move
            delta = (self._wvm.get_frequency() - frequency)/(9.96524)*100.0 # +100 steps = -9.96 GHz of tune
            delta_round = np.round(delta/10.0)*10.0
            final_frequency_fine = self.coarse_seek_frequency(frequency, coarse_motor_position + delta - 140, coarse_motor_position + delta + 140,10,0)
            if final_frequency_fine:
                return True
            else:
                self._motdl.high_precision_move(coarse_motor_position)
                print 'Fine scan failed, went back to coarse position. Frequency is %.2f GHz.' % (self._wvm.get_frequency())
        else:
            logging.error(__name__ + 'Coarse scan failed.')
            return False
    def laser_frequency_seek_old(self, frequency):
        bracket_init = 90.0
        f_low = frequency - bracket_init
        f_high = frequency + bracket_init

        ## quadratic
        #m_low, its = self.brent_search(lambda x: self._a*x*x + self._b*x + self._c - f_low, 60000, 120000, 2, 60)
        #m_high, its = self.brent_search(lambda x: self._a*x*x + self._b*x + self._c - f_high, 60000, 120000, 2, 60)

        # linear
        m_low, its = self.brent_search(lambda x: self._b*(x - self._motor_array_mean)/self._motor_array_std + self._c - f_low, 120000, 190000, 2, 60)
        m_high, its = self.brent_search(lambda x: self._b*(x - self._motor_array_mean)/self._motor_array_std + self._c - f_high, 120000, 190000, 2, 60)

        m_target, its = self.brent_search((lambda mot_pos: self.set_motor_and_measure(int(np.round(mot_pos))) - frequency), m_low, m_high, 20, 14)
        if its == -1:
            print 'Since root was not bracketed, recalibrating the laser and trying again...'
            self._motdl.reference_search()
            try:
                msg_string = [__name__ + ': root not bracketed for brents search in laser_frequency seek. re-referenced laser and recalibrating.']
                slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
            except:
                pass

            self.calibrate_laser()
            m_low, its = self.brent_search(lambda x: self._b*(x - self._motor_array_mean)/self._motor_array_std + self._c - f_low, 120000, 190000, 2, 60)
            m_high, its = self.brent_search(lambda x: self._b*(x - self._motor_array_mean)/self._motor_array_std + self._c - f_high, 120000, 190000, 2, 60)
            m_target, its = self.brent_search((lambda mot_pos: self.set_motor_and_measure(int(np.round(mot_pos))) - frequency), m_low, m_high, 20, 14)

        found_frequency = self._wvm.get_frequency()
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
            # find the closest non-root, assuming the function is smooth
            r_idx = np.argmin(np.array((fa,fb)))
            r = np.array((a,b))
            return r[r_idx], -1
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

    def find_single_mode(self):
        fp_out = self._fp.check_stabilization()
        if fp_out == 2:

            for cur_idx, cur in np.ndenumerate(np.arange(self.params['current_range'][0],self.params['current_range'][1],0.0005)):
                self._toptica.set_current(cur)
                time.sleep(0.25)
                fp_out = self._fp.check_stabilization()
                if fp_out == 1:
                    print 'Motor pos: %d, FP %d, GHz: %.2f, cur %.2f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency(), self._toptica.get_current())
                    break
        return

    def find_stable_frequency(self):
        # purpose of this algorithm is to locate a stable frequency nearest to the current motor and voltage settings
        #
        # the algorithm proceeds as:
        # 1) measure the current dispersion using a 3 point measurement
        # 2) if the laser is not stable, increment the voltage by 1 V units; make sure to bound check
        # 3) repeat the dispersion measurement
        # 4) loop through until the laser is stable. if bound of 90 V has been achieved without stability being found,
        # we decrement the motor position by ~10 GHz. In motor terms, it's about -0.104 GHz per step, so we take about
        # 100 steps in the positive direction after resetting the voltage to 0. If this second sweep doesn't work,
        # set the voltage to 0, print a warning, and return.
        disp_max = 0.3 # GHz
        cur_disp = self.measure_laser_dispersion()
        motor_changes = 0
        while cur_disp > disp_max:
            cur_voltage = self._toptica.get_piezo_voltage()
            proposed_voltage = cur_voltage + 1.0
            if proposed_voltage >= 0.0 and proposed_voltage < 90.0:
                # proposed voltage is OK, so move there.
                self._toptica.set_piezo_voltage(proposed_voltage)
                time.sleep(0.2)
                cur_disp = self.measure_laser_dispersion()
            elif proposed_voltage > 90 and motor_changes == 0:
                self._toptica.set_piezo_voltage(0.0)
                # increment the motor by 100 steps
                self._motdl.high_precision_move(self._motdl.get_position() + 100)
                motor_changes = motor_changes +  1
                time.sleep(0.2)
                cur_disp = self.measure_laser_dispersion()
                print 'Executing motor shift to find stable frequency -- found %.2f MHz disp and %.1f GHz frequency. ' % (cur_disp*1000.0, self._wvm.get_frequency())
            elif proposed_voltage > 90 and motor_changes > 0:
                cur_disp = self.measure_laser_dispersion()
                print 'Voltage and motor changes failed to locate a stable frequency. Dispersion: %.2f MHz' % (cur_disp*1000.00)
            else:
                print 'Voltage proposal out of bounds'
                # this case should never happen. voltage proposal will always be greater than or equal to 0.
        print 'Laser dispersion is now %.2f MHz' % (cur_disp*1000.0)
        return

    def measure_laser_dispersion(self):
        # Monitor laser frequency
        frq_recent = np.zeros(3)
        for zz in range(3):
            time.sleep(0.5)
            frq_recent[zz] = self._wvm.get_frequency()
        cur_disp = np.std(frq_recent)*2.920
        return cur_disp

    def check_laser_stabilization(self):
        # Monitor laser frequency
        frq_recent = np.zeros(3)
        for zz in range(3):
            time.sleep(0.5)
            frq_recent[zz] = self._wvm.get_frequency()
        # Check if laser is stable, if not, wait
        for zz in range(30):
            if (np.max(frq_recent) - np.min(frq_recent)) < 1.0:
                #print 'Laser stabilized. Dispersion %.2f MHz, range %.2f MHz' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                break
            time.sleep(1.0)
            np.roll(frq_recent,1)
            frq_recent[0] = self._wvm.get_frequency()
            if zz >=29:
                print 'Laser frequency readout on wavemeter did not stabilize after 30 seconds. Dispersion %.2f MHz, range %.2f MHz.' % (np.std(frq_recent)*1000.0, 1000.0*(np.max(frq_recent)-np.min(frq_recent)))
                return False

        return True

    def filter_laser_frequencies(self, filter_set, frq_array, total_hits_data, frq1):
        print 'Old filter set %s' % filter_set

        frq_idx = np.nonzero(total_hits_data[:,0])
        threshold = self.params['filter_threshold']
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
        self._b = -297.1101
        self._c = 264841.49
        self._motor_array_mean = 155264.44
        self._motor_array_std = 3086.62
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
        full_attenuation = self.params['power'] - self.params['constant_attenuation'] - np.max((0,np.min((desired_atten,15.5)))) + np.log(self.params['Imod'])/np.log(10.0)*20.0
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
        data.add_value('wavelength dispersion (GHz)')
        data.add_value('wavemeter power (W)')
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
        #self.laser_frequency_seek(self.params['working_set'][-1][1]+8.0)
        frq2 = self.params['working_set'][-1][1]+8.0 + 100.0 #Ghz

        # locate the first frequency to seek to
        self.laser_frequency_seek(self.params['working_set'][0][0]-1.0)

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
            #frequency_sentinel = 0.0
            # Enter the loop for measurement
            t1 = time.time()
            seek_iterations = 0
            while np.size(self.params['working_set']) != 0:
                self._toptica.set_piezo_voltage(self.params['piezo_array'][0])
                if self._stop_measurement == True:
                    print 'Measurement aborted.'
                    self._stop_measurement = True
                    break

                # get the desired target frequency
                #
                target_freq = self.params['working_set'][0][0]
                # add a uniform unit random variable to the frequency, for jitter purposes.
                seek_freq = target_freq-1.0*np.random.uniform()
                print 'Seeking to %.2f GHz.' % seek_freq
                self.laser_frequency_seek(seek_freq)
                # keep track of this position, in case we need to return to it after a reference search
                seeked_position = self._motdl.get_position()
                self.find_single_mode()
                cur_frq = self._wvm.get_frequency()

                # now check if we are strictly lower but not more than 20 GHz lower than the target frequency
                ik = 0
                df = cur_frq-target_freq
                print 'df is %.2f GHz.' % (df)
                if not (df > -10 and df < -0.0):
                    # try a seek once more?
                    self.laser_frequency_seek(seek_freq)

                while (not (df > -10.0 and df < -0.25)) and (ik < 15):
                    # determine the discrepancy between the current and target frequencies
                    self.find_single_mode()
                    cur_frq = self._wvm.get_frequency()
                    df = cur_frq-target_freq
                    print 'df is %.2f GHz...' % (df)
                    if ik == 0:
                        initial_df = df
                    # keep decrementing/incrementing the target frequency by 1 GHz until we get strictly below but not less than 20 GHz away
                    # dispersion is approximately -0.104 GHz/step, so move by about 10 steps upward to decrease by about 1 GHz.
                    current_motor_position = self._motdl.get_position()
                    if df > -0.25:
                        # decrement by ~1.5 GHz. don't go smaller than about 10-15 steps because the motor's control
                        # circuitry won't always respond to small moves (c.f. "deadband").
                        print 'df > -3.0... moving motor to %d' % (current_motor_position+50)
                        self._motdl.set_position(current_motor_position+50)

                        time.sleep(1.0)

                    elif df < -10.0:
                        # increment by ~1.5 GHz
                        print 'df < -20.0... moving motor to %d' % (current_motor_position-50)
                        self._motdl.set_position(current_motor_position-50)
                        time.sleep(1.0)

                    if ik == 10 and not (cur_frq > target_freq-10.0 and cur_frq < target_freq-0.0):
                        # if we've failed to home in by stepping the motor in one direction, just try to seek again after
                        # doing a motor reference search and then re-calibration (time intensive - ~4 minutes).
                        print 'Could not reduce frequency delta to within the desired range -- initial df %.2f GHz, final df %.2f GHz. Doing a reference search, calibration, and re-seek.' % (initial_df, df)
                        self._motdl.reference_search()
                        self._motdl.set_position(seeked_position)
                        #self.calibrate_laser()
                        self.laser_frequency_seek(self.params['working_set'][0][0]-3)

                    ik = ik + 1

                cur_frq = self._wvm.get_frequency()
                print 'Laser frequency set to %.2f GHz (target %.2f GHz)' % (cur_frq, target_freq)
                # Laser frequency seek already checks the stabilization, so this isn't necessary.
                # self.check_laser_stabilization()

                temp_count_data = np.zeros(self.params['piezo_pts'] , dtype='uint32')

                #sweep through the piezo_array voltages, which should go from low to high and then back to low
                for k in range(np.size(self.params['piezo_array'])):

                    #Set the new piezo voltage
                    self._toptica.set_piezo_voltage(self.params['piezo_array'][k])
                    qt.msleep(0.5)

                    # The laser can become unstable, so we want to measure the dispersion of its readout; an increase
                    # in the dispersion will tell us when we're near a mode hop in the cavity and should skip measurement.
                    if self.params['stabilize_laser']:
                        fp_out = self._fp.check_stabilization()
                        # check if Fabry Perot check passed. if not,
                        # try to recover the single-mode behavior
                        if not fp_out == 1:
                            self.find_single_mode()
                        # now re-check the FP, since we might be stable now
                        fp_out = self._fp.check_stabilization()
                        if fp_out == 1:
                            # Fabry Perot check passed, but the wavelength might be wrong.
                            wm_disp = np.zeros(3)
                            for i in range(3):
                                wm_disp[i] = self._wvm.get_frequency()
                                if i != 2:
                                    qt.msleep(0.45) # don't wait if this is the last measurement

                            # now check if it's in bounds of where we actually want to measure
                            cur_frq = wm_disp[-1]
                            if (cur_frq > self.params['working_set'][0][0] and cur_frq < self.params['working_set'][0][1]):
                                # Fabry Perot check and frequency check passed, but the final check is whether the dispersion
                                # is low enough. The final dispersion check occurs in a logic statement below.
                                qt.msleep(0.45)

                                # additionally, get the laser power measured in the wavemeter -- this could also be
                                # correlated to instability in the laser cavity mode
                                cur_power = self._wvm.get_power()
                                cur_frq = wm_disp[-1] # self._wvm.get_frequency()
                                cur_disp = np.std(wm_disp)
                                offset_frq = cur_frq - frq1
                            else:
                                # Frequency is not in bounds, so set the values to nonsense ones and the logical
                                # statement below will cause this point to be skipped
                                cur_power = -1.0
                                offset_frq = cur_frq - frq1
                                # the dispersion being reported as "10" will trigger the logic statement not to measure
                                # here, as well.
                                cur_disp = 10.0
                        else:
                            cur_power = -1.0
                            offset_frq = cur_frq - frq1
                            # the dispersion being reported as "10" will trigger the logic statement not to measure here
                            cur_disp = 10.0
                    else:
                        cur_disp = -1.0
                        cur_frq = self._wvm.get_frequency()
                        offset_frq = cur_frq - frq1
                        cur_power = -1.0


                    # Check for superconductivity every 10 points
                    if np.mod(k,10) == 0:
                        self._snspd.check()

                    # use filter logic
                    filter_inc = False
                    for i, elem in enumerate(self.params['filter_set']):
                        lo, hi = elem
                        if cur_frq < hi and cur_frq > lo:
                            filter_inc = True

                    # every ~10% of the complete sweep, print the current frequency
                    if np.mod(k,np.floor(np.size(self.params['piezo_array'])/100.0)) == 0:
                        print 'Motor pos: %d, piezo %.2f V, FP %d, GHz: %.2f, cur %.2f, cur_disp %.2f' % (self._motdl.get_position(), self._toptica.get_piezo_voltage(), self._fp.check_stabilization(), self._wvm.get_frequency(), self._toptica.get_current(), cur_disp)

                    # Determine if we should measure in the logic statement here
                    # find all nonzero frequencies

                    # Dispersion check of 300 MHz, added 2016-03-23. This ensures anything near a mode hop is rejected.
                    if (cur_frq > self.params['working_set'][0][0] and cur_frq < self.params['working_set'][0][1] and cur_disp < 0.3 ) and not np.any(self.params['bin_size']*0.51 > np.abs(offset_frq - frq_array[np.nonzero(total_hits_data[:,0])])):
                        cts_array_temp = np.zeros(np.size(self.params['freq']))

                        for idx in range(np.size(self.params['freq'])):
                            self._pxi.set_frequency(self.params['freq'][idx]*1.0e9)
                            # use explicit settle method instead of a fixed time delay
                            #self._pxi.wait_until_settled(50.0) # argument is maximum time in milliseconds
                            time.sleep(0.01)
                            cts = self._ni63.get('ctr1')

                            #find where in the 3 column data structure to add counts
                            index = np.searchsorted(frq_array, offset_frq)

                            #update the appropriate two columns keeping track of total counts and total hits
                            total_count_data[index,idx] = total_count_data[index,idx] + cts # temp_count_data[j]
                            total_hits_data[index,idx] = total_hits_data[index,idx] + 1

                            # live plot -- unpack cts_list into args of add_data_point
                            data.add_data_point(cur_frq,self.params['freq'][idx],cts,cur_disp,cur_power)
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
                            plot2d_0.update()


                        plot3d_0.update()
                        qt.msleep(0.002)
                        # check snspd
                        self._snspd.check()

                    # Check if a track should occur. If so, track.
                    if time.time() > track_time:

                        # set the AWG into CW mode for tracking
                        self._awg.sq_forced_jump(1)
                        self.awg_confirm(1)

                        qt.msleep(0.1)
                        # Re-optimize
                        self._fbl.optimize()

                        # Set new track time
                        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
                        self._awg.sq_forced_jump(2)
                        self.awg_confirm(2)
                        qt.msleep(0.1)
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
                        print 'Current frequency %.2f too far beyond end of working set (%.2f, %.2f); breaking' % (cur_frq, self.params['working_set'][0][0], self.params['working_set'][0][1])
                        break
                data.new_block()
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
            avg_cts_array_non0 = np.divide(cts_array_non0.astype(np.float64),hits_array_non0.astype(np.float64))

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
        qt.plots['piezoscan_single_sweep'].save_png()
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




def main():
    # measurement parameters

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
            'readout_length' : 200.0, # ns
            'readout_buffer' : 10.0, # ns
            'ctr_term' : 'PFI2',
            'piezo_start' : 0, #volts
            'piezo_end' : 90, #volts
            'piezo_step_size' : 0.20, # volts (dispersion is roughly ~0.4 GHz/V)
            'bin_size' : 0.1, # GHz, should be same order of magnitude as (step_size * .1 GHz)
            'filter_threshold' : 2.5, # GHz
            'microwaves' : True, # modulate with microwaves on or off
            'microwaves_CW' : False, # are the microwaves CW? i.e. ignore pi pulse length
            'pi_length' : 30.0, # ns
            'off_resonant_laser' : True, # cycle between resonant and off-resonant
            'power' : 5.0, # dBm
            'constant_attenuation' : 14.0, # dBm -- set by the fixed attenuators in setup
            'desired_power' : -19.0, # dBm
            'freq' : [1.3776,], #GHz
            'dwell_time' : 500.0, # ms
            #'filter_set' : ( (270850, 270870), (270950, 270970)),(270810, 270940),
            'filter_set' : [(264885,264950)],#, (270951,270974)],
            'current_range' : [0.275, 0.302], # A
            'temperature_tolerance' : 2.0, # Kelvin
            'MeasCycles' : 1,
            'Imod' : 0.533396, #0.2778,
            'stabilize_laser' : True,
            }





    if msvcrt.kbhit():
        kb_char=msvcrt.getch()
        if kb_char == "q":
            do_track = False

    #topt.set_current(ab(atten_array[i]))
    name_string = 'defect_PL1_nomw_%.2f K' % (ls332.get_kelvinA())
    m = SiC_Toptica_Search_Piezo_Sweep(name_string)
    xsettings['desired_power'] = -19.0
    #xsettings['dwell_time'] = base_dwell*np.power(10.0,atten_array[i]/10.0)


    m.params.from_dict(xsettings)

    do_awg_stuff = True


    m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)

    m.prepare()
    print('Press q to quit...')
    time.sleep(2.0)
    if msvcrt.kbhit():
        kb_char=msvcrt.getch()
        if kb_char == "q":
            return
    print 'Proceeding with measurement ...'
    try:
        msg_string = [__name__ + ': Toptica (Brent search) %s started. Cryostat temperature %.2f K, heater power %.1f percent.' % (name_string, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
        slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
    except:
        pass
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

    try:
        msg_string = [__name__ + ': Toptica (Brent search) %s stopped. Cryostat temperature %.2f K, heater power %.1f percent.' % (name_string, qt.instruments['ls332'].get_kelvinA(), qt.instruments['ls332'].get_heater_output() ) ]
        slack.chat.post_message('#singledefectlab', msg_string, as_user=True)
    except:
        pass

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