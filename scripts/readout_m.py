import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt


class SiC_Readout_Master(m2.Measurement):

    mprefix = 'readout'

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

        # Prepare instruments for measurement and verify FBL output

        # Set the trigger source to internal
        self._ddg.set_trig_source('internal')
        # Compute the overall rate of the measurement sequence
        trigger_rate = np.round(1.0/(1.0e-9 * self.params['trigger_period']),5)
        self._ddg.set_trig_rate(trigger_rate)

        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(self.params['trigger_period']*1e-9-10.0*1e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(self.params['trigger_period']*1e-9-10.0*1e-9)


        self._fbl.optimize()

        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))
        if self.params['readout_delay_end'] < self.params['readout_delay_start']:
            print 'delay_end is lower than delay_start!'


        # Set the DAQ counter dwell time, units milliseconds
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        # Set the gated counter source terminal
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counters prepared.'
        # Reset the RFSG
##        self._pxi.close()
##        self._pxi.init_device()
##        self._pxi.reset_device()
##
##        # Now set the power and initial frequency
##        self._pxi.set_power(self.params['power'])
##        self._pxi.set_frequency(self.params['freq']*1.0e9) # GHz units
##        print 'PXI prepared, power and frequency set.'

        # Measure the DDG speed of set/get.
        tdi = time.time()
        self._ddg.set_delayD(1.0e-6)
        self._ddg.get_delayD()
        tdd = time.time() - tdi
        print 'DDG timing is approx. %.3f ms, using this value for sleeps.' % (tdd*1000.0)
        self._ddgsleep = tdd
        # Set all DDG references
        self._ddg.set_referenceA('T0')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceB('A')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceC('T0')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceD('C')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceE('T0')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceF('E')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceG('T0')
        time.sleep(self._ddgsleep)
        self._ddg.set_referenceH('G')
        time.sleep(self._ddgsleep)




        # Set the AOM delay, offset, polarity, and length
        self._ddg.set_delayA(self.params['AOM_delay']*1.0e-9)
        time.sleep(self._ddgsleep)
        self._ddg.set_polarityAB('pos')
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetAB(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeAB(self.params['AOM_amplitude'])
        time.sleep(self._ddgsleep)
        self._ddg.set_delayB(self.params['AOM_length']*1.0e-9)
        time.sleep(self._ddgsleep)
        # Set the RF switching delay

        self._ddg.set_delayC(self.params['AOM_delay']*1.0e-9 + self.params['AOM_length']*1.0e-9 + self.params['RF_delay']*1.0e-9)
        time.sleep(self._ddgsleep)
        self._ddg.set_polarityCD('pos')
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetCD(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeCD(self.params['RF_amplitude'])
        time.sleep(self._ddgsleep)
        self._ddg.set_delayD(1.0e-9) # Just for initialization
        time.sleep(self._ddgsleep)
        print 'DDG references, delays, polarities, and offsets set.'

        # Set the switch photon readout delay

        self._ddg.set_delayE(self.params['readout_delay']*1.0e-9)
        time.sleep(self._ddgsleep)
        self._ddg.set_polarityEF('pos')
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetEF(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeEF(self.params['readout_amplitude'])
        time.sleep(self._ddgsleep)
        self._ddg.set_delayF(self.params['readout_length']*1.0e-9)
        time.sleep(self._ddgsleep)

        # Set offsets and amplitudes
        self._ddg.set_offsetAB(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeAB(2.5)
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetCD(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeCD(2.5)
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetEF(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeEF(2.5)
        time.sleep(self._ddgsleep)






        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['readout_delay_end'] - self.params['readout_delay_start'])/self.params['readout_delay_step']))
        print '--Readout meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps, %.2f window size)--' % (self.params['readout_delay_start'], self.params['readout_delay_end'], self.params['readout_delay_step'], n_steps, self.params['readout_length'])
        delays = np.linspace(self.params['readout_delay_start'], self.params['readout_delay_end'], n_steps)
        total_count_data = np.zeros(n_steps, dtype='uint32')
        average_count_data = np.zeros(n_steps, dtype='float')
        # Create a data object just for plotting
        data = qt.Data(name='readout_measurement')
        data.add_coordinate('delay (ns)')
        data.add_value('counts')

        N_cmeas = 0
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

        for i in range(self.params['MeasCycles']):
            # Optimize
            prev_aom_delay = self._ddg.get_delayA()
            prev_aom_length = self._ddg.get_delayB()
            prev_readout_delay = self._ddg.get_delayE()
            prev_readout_length = self._ddg.get_delayF()
            self._ddg.set_delayA(0.0)
            self._ddg.set_delayB(self.params['trigger_period']*1e-9-50.0*1e-9)
            self._ddg.set_delayE(0.0)
            self._ddg.set_delayF(self.params['trigger_period']*1e-9-50.0*1e-9)
            time.sleep(self._ddgsleep)
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            self._ddg.set_delayA(prev_aom_delay)
            self._ddg.set_delayB(prev_aom_length)
            self._ddg.set_delayE(prev_readout_delay)
            self._ddg.set_delayF(prev_readout_length)
            time.sleep(self._ddgsleep)


            plot2d_0 = qt.Plot2D(data, name='readout_single_sweep', clear=True)
            plot2d_1 = qt.Plot2D(total_count_data, name='readout_avg', clear=True)

            # Create a copy of the frequency array, so we can modify it
            delays_temp = np.copy(delays)

            # Create an array for the single-sweep data
            temp_count_data = np.zeros(n_steps, dtype='uint32')

            # Enter the loop for measurement
            t1 = time.time()
            for j in range(n_steps):

                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" : break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:
                    # Maybe should check if optimize is successful once that's robust
                    prev_aom_delay = self._ddg.get_delayA()
                    prev_aom_length = self._ddg.get_delayB()
                    prev_readout_delay = self._ddg.get_delayE()
                    prev_readout_length = self._ddg.get_delayF()
                    self._ddg.set_delayA(0.0)
                    self._ddg.set_delayB(self.params['trigger_period']*1e-9-50.0*1e-9)
                    self._ddg.set_delayE(0.0)
                    self._ddg.set_delayF(self.params['trigger_period']*1e-9-50.0*1e-9)
                    time.sleep(self._ddgsleep)
                    fbl.optimize()
                    self._ddg.set_delayA(prev_aom_delay)
                    self._ddg.set_delayB(prev_aom_length)
                    self._ddg.set_delayE(prev_readout_delay)
                    self._ddg.set_delayF(prev_readout_length)
                    time.sleep(self._ddgsleep)

                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

                # Set the new readout delay length
                self._ddg.set_delayE(delays_temp[j]*1.0e-9) # units of ns
                time.sleep(self._ddgsleep)

                temp_count_data[j] = self._ni63.get('ctr1')
                data.add_data_point(delays_temp[j],temp_count_data[j])
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1
            if j == 0:
                # During the first iteration, estimate the complete wall time efficiency
                print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
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
            total_count_data = total_count_data + temp_count_data[delays_temp.argsort()]
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)




        self._pxi.set_status('off')
        # Measurement has ended, so start saving data
##        grp = h5.DataGroup('SiC_Readout_data', self.h5data, base=self.h5base)
##        grp.add('freq', data=freq, unit='GHz', note='frequency')
##        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
##        grp.add('avgcounts', data=average_count_data, unit='counts', note='average counts per measurement cycle')
##        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')


        return



# measurement parameters

xsettings = {
        'fbl_time' : 30.0, # seconds
        'trigger_period' : 4000.0, #ns
        'AOM_delay' : 500.0, # ns
        'AOM_length' : 1500.0, # ns
        'AOM_amplitude' : 2.5, # V
        'RF_delay' : 1445.0, # ns, 945 ns for AOM, 500 ns for singlet
        'RF_length' : 10.0,
        'RF_amplitude' : 2.5, # V
        'readout_amplitude' : 2.5, #V
        'readout_delay' : 100.0,
        'readout_length' : 10.0, # ns
        'readout_delay_start' : 100.0, #ns
        'readout_delay_end' : 1400.0, #ns
        'readout_delay_step' : 5.0, #ns
        'ctr_term' : 'PFI2',
        'power' : -57.0, # dbM
        'readout_length_start' : 0.0, # ns
        'readout_length_end' : 1500.0, # ns
        'readout_length_step' : 20.0, # ns
        'freq' : 1.28, #GHz
        'dwell_time' : 1000.0, # ms
        'temperature_tolerance' : 3.0, # Kelvin
        'MeasCycles' : 350,
        'random' : 0
        }


# Create a measurement object m
m = SiC_Readout_Master('default')

# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)


if m.review_params():
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