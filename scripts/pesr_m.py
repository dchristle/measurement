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

        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))
        if self.params['freq_high'] < self.params['freq_low']:
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
        self._pxi.set_frequency(1.2*1.0e9) # GHz units
        print 'PXI prepared, power and frequency set.'

        # Measure the DDG speed of set/get.
        tdi = time.time()
        self._ddg.get_delayD()
        tdd = time.time() - tdi
        print 'DDG timing is approx. %.3f ms, using this value for sleeps.' % (tdd*1000.0)
        self._ddgsleep = 2*tdd
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
        #self._ddg.set_trig_source('internal')
        # Compute the overall rate of the measurement sequence
        #trigger_rate = np.round(1.0/(1.0e-9 * self.params['trigger_period']),5)
        #self._ddg.set_trig_rate(trigger_rate)

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

        self._ddg.set_delayC(self.params['RF_delay']*1.0e-9)
        time.sleep(self._ddgsleep)
        self._ddg.set_polarityCD('pos')
        time.sleep(self._ddgsleep)
        self._ddg.set_offsetCD(0.0)
        time.sleep(self._ddgsleep)
        self._ddg.set_amplitudeCD(self.params['RF_amplitude'])
        time.sleep(self._ddgsleep)
        self._ddg.set_delayD(self.params['RF_length']*1.0e-9) # Pi pulse length set
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

        self._ddg.set_polarityGH('pos')
        self._ddg.set_offsetGH(0.0)
        self._ddg.set_delayG(0.0)
        self._ddg.set_delayH(1.0/self._trigger_rate-500.0*1.0e-9)





        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['freq_high'] - self.params['freq_low'])/self.params['freq_step']))
        print '--Pulsed spin resonance meas. from %.4f GHz to %.4f GHz in %.4f MHz steps (%.2f steps)--' % (self.params['freq_low'], self.params['freq_high'], 1000.0*self.params['freq_step'], n_steps)
        freqs = np.linspace(self.params['freq_low'], self.params['freq_high'], n_steps)
        total_count_data = np.zeros(n_steps, dtype='uint32')
        average_count_data = np.zeros(n_steps, dtype='float')
        intermediate_total_data = np.zeros( (n_steps, 1), dtype='uint32')


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
            freqs_temp = np.copy(freqs)
            if self.params['random'] == 1:
                # Now shuffle the array in place
                np.random.shuffle(freqs_temp)
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
                self._pxi.set_frequency(freqs_temp[j]*1.0e9) # units of ns
                time.sleep(self._ddgsleep)
                self._ni63.set_count_time(self.params['dwell_time']/1000.0)

                temp_count_data[j] = self._ni63.get('ctr1')
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1
            if i == 0:
                # During the first iteration, estimate the complete wall time efficiency
                print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0)
            sorted_temp_data = temp_count_data[freqs_temp.argsort()]

            plot2d_0 = qt.Plot2D(freqs,sorted_temp_data, name='pesr_single_sweep', clear=True)
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
            total_count_data = total_count_data + temp_count_data[freqs_temp.argsort()]
##            if i == 0:
##                intermediate_total_data[:,0] = total_count_data
##            elif np.mod(i,10) == 0:
##                intermediate_total_data = np.vstack((intermediate_total_data,total_count_data))
            plot2d_1 = qt.Plot2D(freqs,total_count_data, name='pesr_avg', clear=True)
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)




        self._pxi.set_status('off')
        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-500.0*1.0e-9)
        self._ddg.set_delayC(0.0)
        self._ddg.set_delayD(1.0/self._trigger_rate-500.0*1.0e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-500.0*1.0e-9)
        self._ddg.set_delayG(0.0)
        self._ddg.set_delayH(1.0/self._trigger_rate-500.0*1.0e-9)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_PESR_data', self.h5data, base=self.h5base)
        grp.add('frequency', data=freqs, unit='GHz', note='frequency')
        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')
        #grp.add('intermediatecounts', data=intermediate_total_data, unit='', note='intermediate total count data')


        return



# measurement parameters

xsettings = {
        'trigger_period' : 1165.0, # ns
        'fbl_time' : 55.0, # seconds
        'AOM_delay' : 10.0, # ns
        'AOM_length' : 500.0, # ns
        'AOM_amplitude' : 2.5, # V
        'RF_delay' : 10.0, # ns
        'RF_length' : 57.7, # pi pulse length
        'RF_amplitude' : 2.5, # V
        'readout_amplitude' : 2.5, #V
        'readout_delay' : 665.0,
        'readout_length' : 180.0, # ns
        'ctr_term' : 'PFI2',
        'power' : -5.0, # dBm
        'freq_low' : 1.16, # GHz
        'freq_high' : 1.53, # GHz
        'freq_step' : 8*1.25e-4, # GHz
        'freq' : 1.2329, #GHz
        'dwell_time' : 250.0, # ms
        'temperature_tolerance' : 4.0, # Kelvin
        'MeasCycles' : 450,
        'random' : 1
        }


### Create a measurement object m
##m = SiC_Rabi_Master('default')
##
### since params is not just a dictionary, it's easy to incrementally load
### parameters from multiple dictionaries
### this could be very helpful to load various sets of settings from a global
### configuration manager!
##m.params.from_dict(xsettings)
##
##
##if m.review_params():
##    print 'Proceeding with measurement ...'
##    m.prepare()
##    m.measure()
##    m.save_params()
##    m.save_stack()
##else:
##    print 'Measurement aborted!'
##
### important! hdf5 data must be closed, otherwise will not be readable!
### (can also be done by hand, of course)
##m.finish()

p_array = np.array([-5])

for rr in range(np.size(p_array)):
    # Create a measurement object m
    print 'About to proceed -- waiting 5 s for quit (press q to quit)'
    time.sleep(5.0)
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
    name_string = 'power %.2f dBm' % (p_array[rr])
    m = SiC_PESR_Master(name_string)
    xsettings['readout_length'] = 180.0
    xsettings['power'] = p_array[rr]
    # since params is not just a dictionary, it's easy to incrementally load
    # parameters from multiple dictionaries
    # this could be very helpful to load various sets of settings from a global
    # configuration manager!
    m.params.from_dict(xsettings)


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


track_on = True
fbl_t = qt.instruments['fbl']
track_iter = 0
while track_on == True:
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    fbl_t.optimize()
    time.sleep(5.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break