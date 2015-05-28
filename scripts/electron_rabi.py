import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar

reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

# Defines a class called SiC_ESR_Master, and uses class inhereitance to
# derive itself from Measurement class defined in the file measurement.lib
# .measurement2.measurement, which I relabeled as "m2".

class ElectronRabi(m2.Measurement):
    # This prefix is used later for filenames; just indicates that it's an esr
    # measurement
    mprefix = 'ElectronRabi'

    mprefix = 'rabi'

    def prepare(self):
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

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        self._ddg.set_trig_source('internal')
        # Compute the overall rate of the measurement sequence
        self._trigger_rate = np.round(1.0/(1.0e-9 * self.params['trigger_period']),5)
        self._ddg.set_trig_rate(self._trigger_rate)

        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-500.0*1.0e-9)
        self._ddg.set_delayC(0.0)
        self._ddg.set_delayD(1.0/self._trigger_rate-500.0*1.0e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-500.0*1.0e-9)
        self._ddg.set_delayG(0.0)
        self._ddg.set_delayH(1.0/self._trigger_rate-500.0*1.0e-9)


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

        return

    def rabi_sequence(self, upload=True):


        # define the necessary pulses
        X = pulselib.MW_IQmod_pulse('Weak pi-pulse',
            I_channel='MW_Imod', Q_channel='MW_Qmod',
            PM_channel='MW_pulsemod',
            frequency = self.params['MW_pulse_frequency'],
            PM_risetime = self.params['MW_pulse_mod_risetime'])

        T = pulse.SquarePulse(channel='MW_Imod', name='delay',
            length = 200e-9, amplitude = 0.)

        # make the elements - one for each ssb frequency
        elements = []
        for i in range(self.params['pts']):

            e = element.Element('ElectronRabi_pt-%d' % i, pulsar=qt.pulsar)
            e.append(T)

            e.append(pulse.cp(X,
                length = self.params['MW_pulse_durations'][i],
                amplitude = self.params['MW_pulse_amplitudes'][i]))

            elements.append(e)


        # create a sequence from the pulses
        seq = pulsar.Sequence('ElectronRabi sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=True)

        # upload the waveforms to the AWG
        if upload:
            qt.pulsar.upload(*elements)

        # program the AWG
        qt.pulsar.program_sequence(seq)

    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['RF_length_end'] - self.params['RF_length_start'])/self.params['RF_length_step']))
        print '--Rabi meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps)--' % (self.params['RF_length_start'], self.params['RF_length_end'], self.params['RF_length_step'], n_steps)
        RF_lengths = np.linspace(self.params['RF_length_start'], self.params['RF_length_end'], n_steps)
        total_count_data = np.zeros(n_steps, dtype='uint32')
        average_count_data = np.zeros(n_steps, dtype='float')
        intermediate_total_data = np.zeros( (1,n_steps), dtype='uint32')
        signal = np.zeros(1)


        # Set the PXI status to 'on', i.e. generate microwaves
        self._pxi.set_status('on')
        N_cmeas = 0
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process
        # Optimize
        prev_aom_delay = self._ddg.get_delayA()
        prev_aom_length = self._ddg.get_delayB()
        prev_readout_delay = self._ddg.get_delayE()
        prev_readout_length = self._ddg.get_delayF()
        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-500.0*1.0e-9)

        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-500.0*1.0e-9)

        self._ddg.set_delayG(0.0)
        self._ddg.set_delayH(1.0/self._trigger_rate-500.0*1.0e-9)

        time.sleep(self._ddgsleep)
        if self._fbl.optimize() == False:
            if self._fbl.optimize() == False:
                print 'FBL failed twice, breaking.'
        self._ddg.set_delayA(prev_aom_delay)
        self._ddg.set_delayB(prev_aom_length)
        self._ddg.set_delayE(prev_readout_delay)
        self._ddg.set_delayF(prev_readout_length)
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        scan_on = True
        for i in range(self.params['MeasCycles']):






            # Create a copy of the frequency array, so we can modify it
            RF_lengths_temp = np.copy(RF_lengths)
            if self.params['random'] == 1:
                # Now shuffle the array in place
                np.random.shuffle(RF_lengths_temp)
            # Create an array for the single-sweep data
            temp_count_data = np.zeros(n_steps, dtype='uint32')


            # Enter the loop for measurement
            t1 = time.time()
            for j in range(n_steps):

                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:
                    # Maybe should check if optimize is successful once that's robust
                    prev_aom_delay = self._ddg.get_delayA()
                    prev_aom_length = self._ddg.get_delayB()
                    prev_readout_delay = self._ddg.get_delayE()
                    prev_readout_length = self._ddg.get_delayF()
                    self._ddg.set_delayA(0.0)
                    self._ddg.set_delayB(1.0/self._trigger_rate-500.0*1e-9)
                    self._ddg.set_delayE(0.0)
                    self._ddg.set_delayF(1.0/self._trigger_rate-500.0*1e-9)

                    time.sleep(self._ddgsleep)
                    fbl.optimize()


                    self._ddg.set_delayA(prev_aom_delay)
                    self._ddg.set_delayB(prev_aom_length)
                    self._ddg.set_delayE(prev_readout_delay)
                    self._ddg.set_delayF(prev_readout_length)
                    time.sleep(self._ddgsleep)
                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()


                # Set the new RF pulse length
                self._ddg.set_delayD(RF_lengths_temp[j]*1.0e-9) # units of ns
                time.sleep(self._ddgsleep)
                self._ni63.set_count_time(self.params['dwell_time']/1000.0)

                temp_count_data[j] = self._ni63.get('ctr1')
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1

            print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0)

            sorted_temp_data = temp_count_data[RF_lengths_temp.argsort()]

            plot2d_0 = qt.Plot2D(RF_lengths,sorted_temp_data, name='rabi_single_sweep', clear=True)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False: break
            # Now start checking for other issues. If present, stop.
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break
            # Checks have all passed, so proceed...

            # Now add the sorted data array to the total array
            # Use the argsort functionality to sort the count data by the frequnecy
            # it was taken at.
            total_count_data = total_count_data + temp_count_data[RF_lengths_temp.argsort()]
            if i == 0:
                intermediate_total_data[0,:] = total_count_data
            elif np.mod(i,10) == 0:
                intermediate_temp_data = np.zeros( (i/10+1,n_steps), dtype='uint32')
                intermediate_temp_data[:-1,:] = intermediate_total_data
                intermediate_temp_data[i/10,:] = total_count_data
                intermediate_total_data = np.copy(intermediate_temp_data)
                #print 'size is %s' % (np.size(intermediate_total_data))
                #intermediate_total_data = np.vstack((intermediate_total_data,total_count_data))
            if i == 0:
                signal[0] = self._ni63.get('ctr1')
            elif np.mod(i,10):
                signal = np.hstack((signal,self._ni63.get('ctr1')))
            plot2d_1 = qt.Plot2D(RF_lengths,total_count_data, name='rabi_avg', clear=True)
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)




        self._pxi.set_status('off')
        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-200.0*1.0e-9)
        self._ddg.set_delayC(0.0)
        self._ddg.set_delayD(1.0/self._trigger_rate-200.0*1.0e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-200.0*1.0e-9)
        self._ddg.set_delayG(0.0)
        self._ddg.set_delayH(1.0/self._trigger_rate-200.0*1.0e-9)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Rabi_data', self.h5data, base=self.h5base)
        grp.add('length', data=RF_lengths, unit='ns', note='frequency')
        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')
        grp.add('intermediate', data=intermediate_total_data, unit='', note='intermediate total count data')
        grp.add('signal', data=signal, unit='counts', note='signal rate per N iterations')


        return









# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'trigger_period' : 4255.0, # ns
        'fbl_time' : 55.0, # seconds
        'AOM_delay' : 1200.0, # ns
        'AOM_length' : 1600.0, # ns
        'AOM_amplitude' : 2.5, # V
        'RF_delay' : 50.0, # ns
        'RF_length' : 10.0,
        'RF_amplitude' : 2.5, # V
        'readout_amplitude' : 2.5, #V
        'readout_delay' : 1855.0,
        'readout_length' : 210.0, # ns
        'ctr_term' : 'PFI2',
        'power' : -4.0, # dbM
        'RF_length_start' : 0.0, # ns
        'RF_length_end' : 1200.0, # ns
        'RF_length_step' : 5, # ns
        'freq' : 1.2349, #GHz
        'dwell_time' : 500.0, # ms
        'temperature_tolerance' : 2.0, # Kelvin
        'MeasCycles' : 800,
        'random' : 1
        }

# Generate array of powers -- in this case, just one power.

p_low = -19
p_high = -19
p_nstep = 1

p_array = np.linspace(p_low,p_high,p_nstep)


##name_string = 'esr power %.3f dBm' % (p_array[0])
##m = SiC_ESR_Master(name_string)
##xsettings['power'] = p_array[0]
### since params is not just a dictionary, it's easy to incrementally load
### parameters from multiple dictionaries
### this could be very helpful to load various sets of settings from a global
### configuration manager!
##m.params.from_dict(xsettings)


##if m.review_params():
##    print 'Proceeding with measurement ...'
##    m.prepare()
##    m.measure()
##    m.save_params()
##    m.save_stack()
##else:
##    print 'Measurement aborted!'
##m.finish()

# This for loop is for measuring at an array of powers, but there's only one,
# so it will only execute once.

for rr in range(p_nstep):

    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)

    name_string = 'power %.2f dBm' % (p_array[rr])
    # Create a measurement object m with a name we just made indicating the
    # power it's taken at.
    m = SiC_ESR_Master(name_string)
    # Change the xsettings dictionary above's entry for 'power' to the desired
    # power.
    xsettings['power'] = p_array[rr]


    # Load all the parameters in the slightly modified xsettings dictionary
    # into the measurement object 'm' that we just made, which will now have
    # the new power
    m.params.from_dict(xsettings)

    # The if/then here is just leftover from previous code -- since True is
    # always True, it will always execute.
    if True:
        print 'Proceeding with measurement ...'
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



### CONSTANTS AND FLAGS
RAMP = True                 # whether to ramp power up to final value if above threshold
RAMP_POWER = -9             # dBm threshold for power ramping
DO_SEQUENCES = True         # if False, we won't resequence the AWG

DO_READOUT = True ##don't change this!


### Tool functions

def measurement_alert():
    # Alert that measurement has finished
    ea_t = qt.instruments['ea']
    ls332_t = qt.instruments['ls332']
    cur_temp = ls332_t.get_kelvinA()
    msg_string = 'ESR measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
    ea_t.email_alert(msg_string)
    return

def continuous_track():
    track_on = True
    fbl_t = qt.instruments['fbl']
    track_iter = 0
    while track_on == True:
        track_iter = track_iter + 1
        print 'Tracking for %d iteration.' % track_iter
        fbl_t.optimize()
        time.sleep(1.0)
        if msvcrt.kbhit() or track_on == False:
                    kb_char=msvcrt.getch()
                    if kb_char == "q" or track_on == False: break
    return

def start_msmt(m):
    m.autoconfig()
    m.update_definitions()
    m.setup()
    m.run()

def finish_msmnt():
    qt.instruments['awg'].stop()
    qt.msleep(1)
    qt.instruments['awg'].set_runmode('CONT')


def default_msmt(name):

    # setup the master measurement
    m = setup_msmt(''+name)



    # load default settings from xsettings dictionary
    m.params.from_dict(xsettings)

    m.params_lt1['max_CR_starts'] = 1000
    m.params_lt1['teleportation_repetitions'] = -1
    m.params['power'] = 40*60 # seconds -- will actually stop 10 sec earlier.

    if DO_SEQUENCES:
        m.rabi_sequence()

    qt.instruments['ZPLServo'].move_out()
    qt.instruments['ZPLServo_lt1'].move_out()

    start_msmt(m)
    finish_msmnt()

if __name__ == '__main__':
    #name = raw_input('Name?')
    default_msmt('PL2_defect3_meas1')
