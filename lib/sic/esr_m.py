import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt


class SiC_ESR_Master(m2.Measurement):

    mprefix = 'esr'

    def prepare(self):
        # Set up some instruments
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']

        # Prepare instruments for measurement and verify FBL output


        self._fbl.optimize()

        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        if self.params['f_high'] < self.params['f_low']:
            print 'f_high is lower than f_low!'

        # Set the DAQ counter dwell time, units milliseconds
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        # Reset the RFSG
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()

        # Now set the power and initial frequency
        self._pxi.set_power(self.params['power'])
        self._pxi.set_frequency(self.params['f_low']*1.0e9) # GHz units



        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['f_high'] - self.params['f_low'])/self.params['f_step']))

        freq = np.linspace(self.params['f_low'], self.params['f_high'], n_steps)
        total_count_data = np.zeros(n_steps, dtype='uint32')
        average_count_data = np.zeros(n_steps, dtype='float')

        # Set the PXI status to 'on', i.e. generate microwaves
        self._pxi.set_status('on')
        time.sleep(1.0)
        N_cmeas = 0
        # Optimize
        if self._fbl.optimize() == False:
            if self._fbl.optimize() == False:
                print 'FBL failed twice, breaking.'
        print '--ESR meas. from %.4f to %.4f in %.4f MHz steps (%.2f steps)--' % (self.params['f_low'], self.params['f_high'], self.params['f_step']*1000.0, n_steps)
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()


        for i in range(self.params['MeasCycles']):

            scan_on = True
            # Create a data object just for plotting
            data = qt.Data(name='esr_measurement', infile=False,inmem=True,tempfile=False)

            data.add_coordinate('frequency (GHz)')
            data.add_value('counts')

            # Create a copy of the frequency array, so we can modify it
            freq_temp = np.copy(freq)
            # Now shuffle the array in place
            np.random.shuffle(freq_temp)
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
                    print 'Tracking!'
                    # Maybe should check if optimize is successful once that's robust
                    self._fbl.optimize()
                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()

                # Set the frequency
                self._pxi.set_frequency(freq_temp[j]*1.0e9) # Remember GHz
                setVal = self._pxi.wait_until_settled() # default wait is 50 ms

                if setVal != 0:
                    print 'Frequency did not settle!'
                    break
                temp_count_data[j] = self._ni63.get('ctr0')
                #data.add_data_point(freq_temp[j],temp_count_data[j])
            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1
            print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0)
            sorted_temp_data = temp_count_data[freq_temp.argsort()]
            for h in range(n_steps):
                data.add_data_point(freq[h], sorted_temp_data[h])
            plot2d_0 = qt.Plot2D(data, name='esr_single_sweep', clear=True)
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
            total_count_data = total_count_data + temp_count_data[freq_temp.argsort()]
            N_cmeas = N_cmeas + 1
            average_count_data = total_count_data/float(N_cmeas)
            avgdata = qt.Data(name='avg_esr_measurement',infile=False,inmem=True,tempfile=False)

            avgdata.add_coordinate('frequency (GHz)')
            avgdata.add_value('counts')
            for h in range(n_steps):
                avgdata.add_data_point(freq[h], total_count_data[h])


            plot2d_1 = qt.Plot2D(avgdata, name='esr_avg', clear=True)




        self._pxi.set_status('off')
        # Measurement has ended, so start saving data
##        grp = h5.DataGroup('SiC_ESR_data', self.h5data, base=self.h5base)
##        grp.add('freq', data=freq, unit='GHz', note='frequency')
##        grp.add('counts', data=total_count_data, unit='counts', note='total counts')
##        grp.add('avgcounts', data=average_count_data, unit='counts', note='average counts per measurement cycle')
##        grp.add('N_cmeas', data=N_cmeas, unit='', note='total completed measurement cycles')


        return