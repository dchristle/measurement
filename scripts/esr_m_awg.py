import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
import gc
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

# Defines a class called SiC_ESR_Master, and uses class inhereitance to
# derive itself from Measurement class defined in the file measurement.lib
# .measurement2.measurement, which I relabeled as "m2".

class SiC_ESR_Master(m2.Measurement):
    # This prefix is used later for filenames; just indicates that it's an esr
    # measurement
    mprefix = 'esr'
    # Define a 'prepare' function, which prepares the setup to proceed with
    # a measurement by doing things like defining references to the instruments
    # it will use, checking
    def sequence(self, upload=True, program=True, clear=False):
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')
        awg.stop()
        gc.collect()
        time.sleep(1.0)
        if clear:
            awg.clear_waveforms()
        elements = []
        # Create waveform that has laser, microwaves, photon counting, and 1/0 I/Q modulation on
        # all the time for a long period of time (~100 us).
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='laser')
        e.add(pulse.cp(sq_pulseMW, amplitude=1.0, length=100e-6), name='microwaves')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)

        # create a sequence from the pulses -- only one in this case
        seq = pulsar.Sequence('CW ESR Sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(2.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

    def prepare(self):
        track_on = True
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


        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.5)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))

        # Check if the absolute value of the current temperature (kelvinA) and
        # its setpoint is greater than 3.0 K. If it is, print a message saying
        # the temperature is out of range.
        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        if self.params['f_high'] < self.params['f_low']:
            print 'f_high is lower than f_low!'

        # Set the DAQ counter dwell time, units milliseconds
        # This is how long the counter will count for when we give it a command
        # to tell us the counts.
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        # Reset the RFSG - close it, initialize it, then perform a reset.
        # This is just done to clean things up in case it wasn't closed out
        # gracefully before.
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()

        # Now set the power and initial frequency, but don't turn it on yet.
        self._pxi.set_power(self.params['power']-10.0)
        self._pxi.set_frequency((self.params['f_low']+self.params['f_high'])/2.0*1.0e9) # GHz units

        # Now set the proper attenuation
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        if desired_atten > 15.5:
            print 'Cannot attenuate that high -- setting to 15.5 dB attenuation!'
            self._va.set_attenuation(15.5)
            # Should figure out a way to stop the measurement here.
        else:
            self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten
        # Check if the SNSPD is still superconducting
        #if self._snspd.check() == False:
        #    print 'SNSPD went normal and could not restore!'
        # Start the AWG
        self._awg.start()
        time.sleep(1)
        # Set the AWG to the only waveform for this measurement
        self._awg.sq_forced_jump(1)
        time.sleep(0.1)

        return
    # Now define a new method of this particular ESR measurement class called
    # "measure" that does the actual measuring.
    def measure(self):
        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        # Take note of the time the measurement starts.
        t0 = time.time()
        track_on = True
        # Populate some arrays
        # Figure out how many steps we'll take, given a stepsize and freq. range
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['f_high'] - self.params['f_low'])/self.params['f_step']))
        # Make an array of frequencies with that number of steps
        freq = np.linspace(self.params['f_low'], self.params['f_high'], n_steps)
        # If dropout is enabled, dropout the elements and re-count the array
        if self.params['dropout']:
            print 'Dropping out'
            didx = np.logical_and( (freq > self.params['dropout_low']) , (freq < self.params['dropout_high']) )
            print 'Freq size is %s' % freq.size
            freq = np.delete(freq,np.where(didx))
            print 'Now freq size is %s' % freq.size
            n_steps = freq.size

        # Make an empty array to store our total count data in
        total_count_data = np.zeros(n_steps, dtype='uint32')
        # Make an empty array to store tha average count rate data in
        average_count_data = np.zeros(n_steps, dtype='float')

        # Set the PXI status to 'on', i.e. generate microwaves
        print 'Power set to %.2f.' % (self.params['power']-5.0)
        self._pxi.set_power(self.params['power']-5.0)
        self._pxi.set_status('on')
        # Wait 1 second
        time.sleep(1.0)
        # Set the number of completed measurements to 0. This is just a variable
        # that keeps track of how many are done so we can divide the total counts
        # by this number to get an average.
        N_cmeas = 0
        # Optimize the FBL loop. Nominally it should return false if it fails
        # but I haven't implemented that yet.
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return
        # Since microwaves are on, there may be some heating and drifting of the
        # location. So I wait 15 s, then FBL again.
        print 'Waiting 5 s for temperature stabilization...'
        time.sleep(5.0)
        if self._fbl.optimize() == False:
            if self._fbl.optimize() == False:
                print 'FBL failed twice, breaking.'
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return
        # And do it again...
        if self.params['desired_power'] >= -12.0:
            print 'Final power is high -- going to ramp slowly.'
            print 'Waiting 2 s for temperature stabilization, power to %.2f dBm' % (self.params['power']-5.0)
            self._pxi.set_power(self.params['power']-3.0)
            time.sleep(2.0)
            self._keystroke_check('abort')
            if self.keystroke('abort') in ['q','Q']:
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                self._pxi.set_status('off')
                return
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            print 'Waiting 2 s for temperature stabilization, power to %.2f dBm' % (self.params['power']-2.0)
            time.sleep(2.0)
            self._keystroke_check('abort')
            if self.keystroke('abort') in ['q','Q']:
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                self._pxi.set_status('off')
                return
            self._pxi.set_power(self.params['power']-1.5)
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            self._keystroke_check('abort')
            if self.keystroke('abort') in ['q','Q']:
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                self._pxi.set_status('off')
                return
        # Set to final power
        print 'Setting final power to %.2f' % (self.params['power'])

        current_setpoint = self._ls332.get_setpoint1()
        mm = 0
        kk = 0
        while kk < 22 and mm < 1:
            print 'Stabilization loop index is %d' % kk
            time.sleep(2.0)
            self._fbl.optimize()
            current_temp = self._ls332.get_kelvinA()
            if np.abs(current_temp-current_setpoint) < 0.15:
                mm = mm + 1
            else:
                mm = 0
            kk = kk+1

            if kk == 21:
                print 'Temperature did not stabilize!!'
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                self._pxi.set_status('off')
                return

        print 'Stabilized.'
        time.sleep(3.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return

        # OK now just print the measurement we're doing with power, freq. range,
        # and step size for the user to see.
        print '--ESR meas. at %.3f dBm from %.4f to %.4f in %.4f MHz steps (%.2f steps)--' % (self.params['desired_power'],self.params['f_low'], self.params['f_high'], self.params['f_step']*1000.0, n_steps)
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process, too (5 seconds) so the track
        # time is uncorrelated with the measurement.
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

        # This is the main measurement loop that will execute a certain number of
        # times.
        for i in range(self.params['MeasCycles']):
            # This variable is used to break out of several loops, i.e. if one
            # loop breaks it should set scan on to False, and then others will
            # check for scan_on being False and break, too.
            scan_on = True
            # Create a data object just for plotting
##            data = qt.Data(name='esr_measurement')
##
##            data.add_coordinate('frequency (GHz)')
##            data.add_value('counts')

            # Create a copy of the frequency array, so we can modify it
            freq_temp = np.copy(freq)
            # Now shuffle the array in place -- this will make it so we take the
            # data in random order each sweep.
            np.random.shuffle(freq_temp)
            # Create an array for the single-sweep data
            temp_count_data = np.zeros(n_steps, dtype='uint32')

            # Enter the loop for actually sweeping the ESR frequencies
            t1 = time.time() # Take note of the time
            for j in range(n_steps):
                # This code will check for a keyboard key press of "q", and if
                # it detects it, sets scan_on = False and then breaks out of
                # this loop. If scan_on is false, the outer loop checks scan_on,
                # discovers it's false, and then breaks that, too, to end the
                # measurement gracefully. If we didn't use this scan_on variable
                # this inner loop would break, but the outer loop would finish
                # normally and then start over (i.e. not break) and although
                # one sweep would have stopped from our keyboard press of "q",
                # the script would just continue to the next sweep, which we
                # don't want.
                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q']:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    scan_on = False
                    break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:
                    print 'Tracking!'

                    self._fbl.optimize()
                    # Set new track time, fbl_time into the future plus a small
                    # random time.
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

                # Set the frequency to the new desired frequency
                self._pxi.set_frequency(freq_temp[j]*1.0e9) # Remember GHz
                setVal = self._pxi.wait_until_settled() # default wait is 50 ms
                # Check if the frequency settled using an internal routine of
                # the RFSG
                if setVal != 0:
                    print 'Frequency did not settle!'
                    break
                # Read the counts on ctr0 (PFI0) that accumulate in one dwell
                # time (set above).
                temp_count_data[j] = self._ni63.get('ctr0')
                #data.add_data_point(freq_temp[j],temp_count_data[j])


            tt = time.time() - t1
            print 'Cycle %d/%d, time is %.3f, efficiency of %.2f percent. Heater output at %.1f.' % (i+1, self.params['MeasCycles'], tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0, self._ls332.get_heater_output())
            # Sort the count data versus frequency we just took so that it's in
            # regular order, i.e. lowest to highest frequency.
            sorted_temp_data = temp_count_data[freq_temp.argsort()]
##            for h in range(n_steps):
##                data.add_data_point(freq[h], sorted_temp_data[h])
##            plot2d_0 = qt.Plot2D(data, name='esr_single_sweep', clear=True)
            # Make a plot of the single sweep data we just took
            plot2d_0 = qt.Plot2D(freq,sorted_temp_data, name='esr_single_sweep', clear=True)
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            if scan_on == False:
                break
            self._keystroke_check('abort')
            if self.keystroke('abort') in ['q','Q']:
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                break
            # Now start checking for other issues. If present, stop.
            # Check if setpoint and actual temperature are within a tolerance
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
            # Check if the SNSPD is still superconducting
            #if self._snspd.check() == False:
            #    print 'SNSPD went normal and could not restore, breaking.'
            #    break
            # Checks have all passed, so proceed...

            # Now add the sorted data array to the total array
            # Use the argsort functionality to sort the count data by the frequnecy
            # it was taken at.
            total_count_data = total_count_data + temp_count_data[freq_temp.argsort()]
            # Successfully completed another sweep, so increment the variable that
            # keeps track of this.
            N_cmeas = N_cmeas + 1
            # Update the average_count_data to be the new total divided by the
            # number of successful measurements.
            average_count_data = total_count_data/N_cmeas

            # Plot the total counts (even though it's mislabeled as esr_avg here)
            plot2d_1 = qt.Plot2D(freq,total_count_data, name='esr_avg', clear=True)
            qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.



        # Once all sweeps are complete, turn off the microwaves.
        self._pxi.set_status('off')
        # Measurement has ended, so start saving data in an HDF5 object.

        grp = h5.DataGroup('SiC_ESR_data', self.h5data, base=self.h5base)
        grp.add('freq', data=freq, unit='GHz', note='frequency')
        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
        grp.add('avgcounts', data=average_count_data, unit='counts', note='average counts per measurement cycle')
        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')

        # Now the measure method of the class is complete, so return.
        return


# Above we defined a new ESR measurement class. Here we're going to create a
# a dictionary of measurement-specific parameters, create some of these ESR
# measurement class objects, and run them.


# measurement parameters - these are all references above.

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 155.0, # seconds
        'ctr_term' : 'PFI0',
        'power' : 5.0, # dBm
        'constant_attenuation' : 28.0, # dBm -- set by the fixed attenuators in setup
        'desired_power' : -7.0, # dBm
        'f_low' : 1.29, #GHz
        'f_high' : 1.379, #Ghz
        'f_step' : 2*4*1.25e-4, #Ghz
        'dwell_time' : 550.0, # ms
        'temperature_tolerance' : 3.0, # Kelvin
        'MeasCycles' : 800,
        'trigger_period' : 100000.0, #ns
        'dropout' : False,
        'dropout_low' : 1.325, # GHz
        'dropout_high' : 1.342, # GHz
        }




# Generate array of powers -- in this case, just one power.

p_low = -29
p_high = -29
p_nstep = 1

p_array = np.linspace(p_low,p_high,p_nstep)



for rr in range(p_nstep):

    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
    name_string = 'power %.2f dBm' % (p_array[rr])
    # Create a measurement object m with a name we just made indicating the
    # power it's taken at.
    m = SiC_ESR_Master(name_string)
    # Change the xsettings dictionary above's entry for 'power' to the desired
    # power.
    xsettings['desired_power'] = p_array[rr]


    # Load all the parameters in the slightly modified xsettings dictionary
    # into the measurement object 'm' that we just made, which will now have
    # the new power
    m.params.from_dict(xsettings)

    # The if/then here is just leftover from previous code -- since True is
    # always True, it will always execute.
    if True:
        print 'Proceeding with measurement ...'
        do_awg_stuff = True
        m.sequence(upload=do_awg_stuff, program=do_awg_stuff, clear=do_awg_stuff)
        m.prepare()
        m.measure()
        # Save params and save stack I think just save the entire set of parameters
        # in the entire setup somewhere and also save a copy of this file every
        # time a measurement executes alongside the data, so that if it gets
        # modified, we can still go back and look at the original one.
        m.save_params()
        m.save_stack()
    else:
        print 'Measurement aborted!'

    # important! hdf5 data must be closed, otherwise will not be readable!
    # (can also be done by hand, of course)
    # m.finish() will close the HDF5 and end the measurement.
    m.finish()
    # I think that save_params, save_stack, and finish are all methods that are
    # inherited from the measurement class m2.Measurement done at the beginning
    # of the class definition.

# Alert that measurement has finished
ea_t = qt.instruments['ea']
ls332_t = qt.instruments['ls332']
cur_temp = ls332_t.get_kelvinA()
msg_string = 'ESR measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)

# Now just keep tracking until 'q' is pressed
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