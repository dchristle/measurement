import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt


class SiC_Saturation_Master(m2.Measurement):

    mprefix = 'saturation'

    def prepare(self):
        # Set up some instruments
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._flip = qt.instruments['flip']

        # Prepare instruments for measurement and verify FBL output


        self._fbl.optimize()

        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        if self.params['current_high'] < self.params['current_low']:
            print 'f_high is lower than f_low!'


        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        # Set the DAQ counter dwell time, units milliseconds
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)






        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['current_high'] - self.params['current_low'])/self.params['current_step']))

        tl_current = np.linspace(self.params['current_high'], self.params['current_low'], n_steps)
        signal_count_data = np.zeros(n_steps, dtype='uint32')
        background_count_data = np.zeros(n_steps, dtype='uint32')
        power_data = np.zeros(n_steps, dtype='float')
        # Start at the highest current
        self._tl.set_current(tl_current[0])


        time.sleep(1.0)
        N_cmeas = 0


        print '--Saturation meas. from %.4f to %.4f in %.4f A steps (%.2f steps)--' % (self.params['current_low'], self.params['current_high'], self.params['current_step']*1000.0, n_steps)
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process





        scan_on = True
        # Create a data object just for plotting
##            data = qt.Data(name='esr_measurement')
##
##            data.add_coordinate('frequency (GHz)')
##            data.add_value('counts')

        # Create a copy of the current array, so we can modify it
        current_temp = np.copy(tl_current)
        # Now shuffle the array in place
        #np.random.shuffle(freq_temp)
        # Create an array for the single-sweep data

        # Enter the loop for measurement
        t1 = time.time()
        for j in range(n_steps):

            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" :
                    scan_on = False
                    break
            print 'Tracking before setting current.'
            self._fbl.optimize()
            tl.set_current(tl_current[j])
            print 'Current set to %.3f A, tracking.' % (tl_current[j])
            time.sleep(1.0)
            self._fbl.optimize()
            # Measure power
            power_data[j] = self._flip.measure_power()


            signal_count_data[j] = self._ni63.get('ctr0')
            x0 = self._fsm.get_abs_positionX()
            self._fsm.set_abs_positionX(x0 + self.params['xdisp'])
            background_count_data[j] = self._ni63.get('ctr0')
            self._fsm.set_abs_positionX(x0)
            print 'Measured %d counts over %d counts (%.2f ratio, %.2f defect counts)' % (signal_count_data[j],background_count_data[j],float(signal_count_data[j])/float(background_count_data[j]), signal_count_data[j]-background_count_data[j])

            #data.add_data_point(freq_temp[j],temp_count_data[j])
        # Check for a break, and break out of this loop as well.
        # It's important to check here, before we add the array to the total
        # since doing it the other way risks adding incomplete data to the
        # total array.
        tt = time.time() - t1
        print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['dwell_time']/1000.0)/tt*100.0)

##            for h in range(n_steps):
##                data.add_data_point(freq[h], sorted_temp_data[h])
##            plot2d_0 = qt.Plot2D(data, name='esr_single_sweep', clear=True)


        # Now start checking for other issues. If present, stop.
        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
            print 'Temperature out of bounds, breaking.'

        if self._snspd.check() == False:
            print 'SNSPD went normal and could not restore, breaking.'

        # Checks have all passed, so proceed...

        # Now add the sorted data array to the total array
        # Use the argsort functionality to sort the count data by the frequnecy
        # it was taken at.

##            avgdata = qt.Data(name='avg_esr_measurement')
##
##            avgdata.add_coordinate('frequency (GHz)')
##            avgdata.add_value('counts')
##            for h in range(n_steps):
##                avgdata.add_data_point(freq[h], total_count_data[h])

        subtracted_data = signal_count_data - background_count_data
        plot2d_1 = qt.Plot2D(power_data,signal_count_data,power_data,background_count_data,name='saturation S/B curves', clear=True)
        plot2d_0 = qt.Plot2D(power_data,subtracted_data, name='saturation subtracted data', clear=True)





        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_saturation_data', self.h5data, base=self.h5base)
        grp.add('current', data=tl_current, unit='A')
        grp.add('power', data=power_data, unit='A')
        grp.add('signal_counts', data=signal_count_data, unit='counts', note='signal counts')
        grp.add('background_counts', data=background_count_data, unit='counts', note='background_counts')


        return



# measurement parameters

xsettings = {
        'ctr_term' : 'PFI0',
        'current_low' : 0.08, # A
        'current_high' : 0.671, # A
        'current_step' : 0.05, # A
        'dwell_time' : 2000.0, # ms
        'temperature_tolerance' : 3.0, # Kelvin
        'MeasCycles' : 1,
        'xdisp' : -2.5
        }





# Create a measurement object m

name_string = 'PL saturation'
m = SiC_Saturation_Master(name_string)

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