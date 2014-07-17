import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt


class SiC_PESR_Master(m2.Measurement):

    mprefix = 'pesr'

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


        self._fbl.optimize()

        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))
        if self.params['delay_end'] < self.params['delay_start']:
            print 'delay_end is lower than delay_start!'


        # Set the DAQ counter dwell time, units milliseconds
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        print 'Counter prepared.'
        # Reset the RFSG
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()

        # Now set the power and initial frequency
        self._pxi.set_power(self.params['power'])
        self._pxi.set_frequency(self.params['freq']*1.0e9) # GHz units
        print 'PXI prepared, power and frequency set.'

        # Measure the DDG speed of set/get.
        tdi = time.time()
        self._ddg.set_delayD(6.0e-6)
        self._ddg.get_delayD()
        tdd = time.time() - tdi
        print 'Roundtrip DDG timing is approx. %.3f ms, using this value for sleeps.' % (tdd*1000.0)
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


        # Set the trigger source to internal
        self._ddg.set_trig_source('internal')
        # Compute the overall rate of the measurement sequence
        trigger_rate = np.round(1.0/(1.0e-9 * self.params['trigger_period']),5)
        self._ddg.set_trig_rate(trigger_rate)

        # Set the AOM delay, offset, polarity, and length
        self._ddg.set_delayA(self.params['A_delay']*1.0e-9)
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

        self._ddg.set_delayC(self.params['A_delay']*1.0e-9 + self.params['AOM_length']*1.0e-9 + self.params['RF_delta']*1.0e-9)
        time.sleep(self._ddgsleep)
        self._ddg.set_polarityCD('pos')
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetCD(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeCD(self.params['RF_amplitude'])
        time.sleep(self._ddgsleep)
        self._ddg.set_delayD(2.0e-9) # Just for initialization
        time.sleep(self._ddgsleep)
        print 'DDG references, delays, polarities, and offsets set.'

        # Set the switch photon readout delay

        self._ddg.set_delayG(self.params['A_delay']*1.0e-9 + self.params['readout_delay']*1.0e-9)
        time.sleep(self._ddgsleep)
        self._ddg.set_polarityGH('pos')
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetGH(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeGH(self.params['readout_amplitude'])
        time.sleep(self._ddgsleep)
        self._ddg.set_delayH(self.params['readout_length']*1.0e-9)
        time.sleep(self._ddgsleep)

        # Set the period of the entire measurement







        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['delay_end'] - self.params['delay_start'])/self.params['delay_step']))
        print '--Rabi meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps)--' % (self.params['delay_start'], self.params['delay_end'], self.params['delay_step'], n_steps)
        delays = np.linspace(self.params['delay_start'], self.params['delay_end'], n_steps)
        total_count_data = np.zeros(n_steps, dtype='uint32')
        average_count_data = np.zeros(n_steps, dtype='float')
        # Create a data object just for plotting
        data = qt.Data(name='rabi_measurement')

        data.add_coordinate('delay (ns)')
        data.add_value('counts')
        # Set the PXI status to 'on', i.e. generate microwaves
        self._pxi.set_status('on')
        N_cmeas = 0

        for i in range(self.params['MeasCycles']):
            # Optimize
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'


            plot2d_0 = qt.Plot2D(total_count_data, name='rabi_single_sweep', clear=True)
            plot2d_1 = qt.Plot2D(total_count_data, name='rabi_avg', clear=True)

            # Create a copy of the frequency array, so we can modify it
            delays_temp = np.copy(delays)
            if self.params['random'] == 1:
                # Now shuffle the array in place
                np.random.shuffle(delays_temp)
            # Create an array for the single-sweep data
            temp_count_data = np.zeros(n_steps, dtype='uint32')
            # Set a time that controls when the next feedback occurs
            # Add a bit of randomness to this process
            track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
            # Enter the loop for measurement
            t1 = time.time()
            for j in range(n_steps):

                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" : break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:
                    # Maybe should check if optimize is successful once that's robust
                    fbl.optimize()
                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

                # Set the new RF pulse length
                self._ddg.set_delayD(delays_temp[j]*1.0e-9) # units of ns
                time.sleep(self._ddgsleep)

                temp_count_data[j] = self._ni63.get('ctr0')
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
            total_count_data = total_count_data + temp_count_data(delays_temp.argsort())
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)




        self._pxi.set_status('off')
        # Measurement has ended, so start saving data
##        grp = h5.DataGroup('SiC_ESR_data', self.h5data, base=self.h5base)
##        grp.add('freq', data=freq, unit='GHz', note='frequency')
##        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
##        grp.add('avgcounts', data=average_count_data, unit='counts', note='average counts per measurement cycle')
##        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')


        return



# measurement parameters

xsettings = {
        'fbl_time' : 30.0, # seconds
        'A_delay' : 10.0, # ns
        'AOM_length' : 1200.0, # ns
        'AOM_amplitude' : 2.5, # V
        'RF_delta' : 700.0, # ns
        'RF_length' : 1445.0, # ns, 945 ns for AOM, 500 ns for singlet
        'RF_amplitude' : 2.5, # V
        'readout_amplitude' : 2.5, #V
        'readout_delay' : 400.0, # ns
        'readout_length' : 400.0, # ns
        'trigger_period' : 4000.0, # ns
        'ctr_term' : 'PFI0',
        'power' : -25.0, # dbM
        'delay_start' : 0.0, # ns
        'delay_end' : 400.0, # ns
        'delay_step' : 10.0, # ns
        'freq' : 1.28, #GHz
        'dwell_time' : 100.0, # ms
        'temperature_tolerance' : 3.0, # Kelvin
        'MeasCycles' : 35,
        'random' : 1
        }


# Create a measurement object m
m = SiC_PESR_Master('default')

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