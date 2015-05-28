import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt


class SiC_BasicPol_Master(m2.Measurement):

    mprefix = 'basicpol'

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
        self._flip = qt.instruments['flip']
        self._st0 = qt.instruments['Standa0']
        self._fm = qt.instruments['fm']

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        self._ddg.set_trig_source('internal')
        # Compute the overall rate of the measurement sequence
        self._trigger_rate = np.round(1.0/(1.0e-9 * self.params['trigger_period']),5)
        self._ddg.set_trig_rate(self._trigger_rate)

        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-100.0*1.0e-9)
        self._ddg.set_delayC(0.0)
        self._ddg.set_delayD(1.0/self._trigger_rate-100.0*1.0e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-100.0*1.0e-9)


        self._fbl.optimize()

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
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        print 'Counter prepared.'


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
        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-100.0*1e-9)

        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-100.0*1e-9)
##        self._ddg.set_delayA(self.params['AOM_delay']*1.0e-9)
##        time.sleep(self._ddgsleep)
##        self._ddg.set_polarityAB('pos')
##        time.sleep(self._ddgsleep)
##        self._ddg.set_offsetAB(0.0)
##        time.sleep(self._ddgsleep)
##        self._ddg.set_amplitudeAB(self.params['AOM_amplitude'])
##        time.sleep(self._ddgsleep)
##        self._ddg.set_delayB(self.params['AOM_length']*1.0e-9)
##        time.sleep(self._ddgsleep)
##        # Set the RF switching delay
##
##        self._ddg.set_delayC(self.params['RF_delay']*1.0e-9)
##        time.sleep(self._ddgsleep)
##        self._ddg.set_polarityCD('pos')
##        time.sleep(self._ddgsleep)
##        self._ddg.set_offsetCD(0.0)
##        time.sleep(self._ddgsleep)
##        self._ddg.set_amplitudeCD(self.params['RF_amplitude'])
##        time.sleep(self._ddgsleep)
##        self._ddg.set_delayD(self.params['RF_length']*1.0e-9) # Just for initialization
##        time.sleep(self._ddgsleep)
##        print 'DDG references, delays, polarities, and offsets set.'
##
##        # Set the switch photon readout delay
##
##        self._ddg.set_delayE(self.params['readout_delay']*1.0e-9)
##        time.sleep(self._ddgsleep)
##        self._ddg.set_polarityEF('pos')
##        time.sleep(self._ddgsleep)
##        self._ddg.set_offsetEF(0.0)
##        time.sleep(self._ddgsleep)
##        self._ddg.set_amplitudeEF(self.params['readout_amplitude'])
##        time.sleep(self._ddgsleep)
##        self._ddg.set_delayF(self.params['readout_length']*1.0e-9)
##        time.sleep(self._ddgsleep)





        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['hwp_step_end'] - self.params['hwp_step_start'])/self.params['hwp_step_step']))
        print '--Basic polarization meas. from %.2f step to %.2f step in %.2f step steps (%.2f steps)--' % (self.params['hwp_step_start'], self.params['hwp_step_end'], self.params['hwp_step_step'], n_steps)
        step_positions = np.linspace(self.params['hwp_step_start'], self.params['hwp_step_end'], n_steps)
        total_count_data = np.zeros(n_steps, dtype='uint32')
        average_count_data = np.zeros(n_steps, dtype='float')
        powers = np.zeros(n_steps, dtype='float')



        self._st0.set_speed(5000)
        self._st0.move(int(np.round(self.params['hwp_step_start'])))
        n = 0
        while n < 50:
            cur_pos = self._st0.get_position()
            if cur_pos == np.round(self.params['hwp_step_start']):
                break
            else:
                n = n + 1
                qt.msleep(0.5)
        N_cmeas = 0
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process
        # Optimize
        prev_aom_delay = self._ddg.get_delayA()
        prev_aom_length = self._ddg.get_delayB()
        prev_readout_delay = self._ddg.get_delayE()
        prev_readout_length = self._ddg.get_delayF()
        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-100.0*1e-9)

        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-100.0*1e-9)
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
        data = qt.Data(name='testmeasurement')

        # Now you provide the information of what data will be saved in the
        # datafile. A distinction is made between 'coordinates', and 'values'.
        # Coordinates are the parameters that you sweep, values are the
        # parameters that you readout (the result of an experiment). This
        # information is used later for plotting purposes.
        # Adding coordinate and value info is optional, but recommended.
        # If you don't supply it, the data class will guess your data format.
        data.add_coordinate('step')
        data.add_value('PL')
        data.add_value('power')
        plot2d = qt.Plot2D(data, name='measure_cvsp', coorddim=0, valdim=1)
        plot2d0 = qt.Plot2D(data, name='measure_pvsp', coorddim=0, valdim=2)
        # The next command will actually create the dirs and files, based
        # on the information provided above. Additionally a settingsfile
        # is created containing the current settings of all the instruments.
        data.create_file()
        for i in range(self.params['MeasCycles']):








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
                    self._ddg.set_delayB(1.0/self._trigger_rate-100.0*1e-9)
                    self._ddg.set_delayE(0.0)
                    self._ddg.set_delayF(1.0/self._trigger_rate-100.0*1e-9)
                    time.sleep(self._ddgsleep)
                    fbl.optimize()

                    self._ddg.set_delayA(prev_aom_delay)
                    self._ddg.set_delayB(prev_aom_length)
                    self._ddg.set_delayE(prev_readout_delay)
                    self._ddg.set_delayF(prev_readout_length)
                    time.sleep(self._ddgsleep)
                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

                s = step_positions[j]
                self._st0.move(int(np.round(s)))
                n = 0
                while n < 50:
                    cur_pos = self._st0.get_position()
                    if cur_pos == np.round(s):
                        break
                    else:
                        n = n + 1
                        qt.msleep(0.5)

                self._ni63.set_count_time(self.params['dwell_time']/1000.0)

                temp_count_data[j] = self._ni63.get('ctr0')

                cur_pow = self._flip.measure_power()
                powers[j] = cur_pow
                data.add_data_point(s,temp_count_data[j],powers[j])
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1
            if j == 0:
                # During the first iteration, estimate the complete wall time efficiency
                print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0)
            sorted_temp_data = temp_count_data

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
            total_count_data = total_count_data + temp_count_data
            plot2d_1 = qt.Plot2D(RF_lengths,total_count_data, name='pol_avg', clear=True)
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)




        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-100.0*1.0e-9)
        self._ddg.set_delayC(0.0)
        self._ddg.set_delayD(1.0/self._trigger_rate-100.0*1.0e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-100.0*1.0e-9)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_BasicPol_data', self.h5data, base=self.h5base)
        grp.add('position', data=step_positions, unit='steps', note='stepper motor angle')
        grp.add('power', data=powers, units='W')
        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')


        return



# measurement parameters

xsettings = {
        'trigger_period' : 5000.0, # ns
        'fbl_time' : 50.0, # seconds
        'AOM_delay' : 1000.0, # ns
        'AOM_length' : 2500.0, # ns
        'AOM_amplitude' : 2.5, # V
        'RF_delay' : 150.0, # ns
        'RF_length' : 10.0,
        'RF_amplitude' : 2.5, # V
        'readout_amplitude' : 2.5, #V

        'ctr_term' : 'PFI0',
        'hwp_step_start' : 0, # ns
        'hwp_step_end' : 200000, # ns
        'hwp_step_step' : 1000, # ns
        'dwell_time' : 500.0, # ms
        'temperature_tolerance' : 3.0, # Kelvin
        'MeasCycles' : 1,
        'random' : 1
        }





m = SiC_BasicPol_Master('default')

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