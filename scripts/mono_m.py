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

# Defines a class called SiC_Monochromator_Master, DC 2015

class SiC_Monochromator_Master(m2.Measurement):
    # This prefix is used later for filenames; just indicates that it's an esr
    # measurement
    mprefix = 'mono'
    # Define a 'prepare' function, which prepares the setup to proceed with
    # a measurement by doing things like defining references to the instruments
    # it will use, checking the SNSPD, and setting the initial monochromator
    # wavelength.

    def prepare(self):
        # you indicate that a measurement is about to start and other
        # processes should stop (like batterycheckers, or temperature
        # monitors)
        qt.mstart()
        # Set a few other variables
        track_on = True
        # Monitor keyboard for abort
        self.start_keystroke_monitor('abort')
        # Set up some instruments
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._ph = qt.instruments['ph']
        self._mc = qt.instruments['mc']
        # Set the wavelength speed and initial wavelength
        self._mc.set_wavelength_speed(self.params['wavelength_speed'])
        print 'Monochromator wavelength speed set to %.2f nm per minute' % self.params['wavelength_speed']
        time.sleep(0.2)
        self._mc.set_wavelength(self.params['wavelength_start'])
        print 'Monochromator initial wavelength set to %.3f nm' % self.params['wavelength_start']
        print 'Retrieving wavelength...'
        current_wl = self._mc.get_wavelength()
        if np.abs(current_wl - self.params['wavelength_start']) < 1.0:
            print 'Current wavelength is %.3f and within 1 nm, continuing...' % current_wl
        else:
            print 'Current wavelength is %.3f and out of bounds!!'


        # Set the DAQ counter dwell time, units milliseconds
        # This is how long the counter will count for when we give it a command
        # to tell us the counts.
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        print 'Counter prepared.'

        # Check if the SNSPD is still superconducting
        if self._snspd.check() == False:
            print 'SNSPD went normal and could not restore!'
        time.sleep(0.1)
        # Create a data structure to store the data in
        self._data = qt.Data(name='monochromator')
        self._data.add_coordinate('Wavelength (nm)')
        self._data.add_value('Counts per second')
        self._data.create_file()
        # Now define a plot for plotting the data as we take it
        self._plot2d = qt.Plot2D(self._data, name='monochromator_sweep', coorddim=0, valdim=1, traceofs=10)
        return
    # Now define a new method of this particular monochromator measurement class called
    # "measure" that does the actual measuring.
    def measure(self):
        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        # Take note of the time the measurement starts.
        t0 = time.time()
        scan_on = True
        # Define a linear array of wavelengths to sweep
        # First compute the number of steps
        n_steps = np.ceil((self.params['wavelength_end'] - self.params['wavelength_start'])/self.params['wavelength_step'])
        # Now define the actual array of wavelengths
        wl_array = self.params['wavelength_start'] + self.params['wavelength_step']*np.arange(n_steps)
        # Now enter the main measurement loop, since the MC should already be at the initial wavelength
        for wl in wl_array:
            if scan_on == False:
                break
            self._keystroke_check('abort')
            if self.keystroke('abort') in ['q','Q']:
                print 'Measurement aborted.'
                self.stop_keystroke_monitor('abort')
                break
            # Now start checking for other issues. If present, stop.
            # Check if the SNSPD is still superconducting
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break
            # Set the wavelength
            self._mc.set_wavelength(wl)
            time.sleep(0.25)
            cps = self._ni63.get('ctr0')*1.0/(self.params['dwell_time']/1000.0)
            # Add the new data point to the data structure
            self._data.add_data_point(wl,cps)

            # Plot the total counts (even though it's mislabeled as esr_avg here)

            qt.msleep(0.003) # keeps GUI responsive and checks if plot needs updating.





        # Now the measure method of the class is complete, so close the file
        self._data.close_file()
        self._plot2d.save_png()
        # lastly tell the secondary processes (if any) that they are allowed to start again.
        qt.mend()
        return


# Above we defined a new ESR measurement class. Here we're going to create a
# a dictionary of measurement-specific parameters, create some of these ESR
# measurement class objects, and run them.


# measurement parameters - these are all references above.

xsettings = {
        'ctr_term' : 'PFI0',
        'wavelength_start' : 1000.0, # nm
        'wavelength_end' : 1070.0, # nm
        'wavelength_step' : 10.0, # nm
        'wavelength_speed' : 300.0, # nm/min grating speed
        'dwell_time' : 500.0, # ms
        'ctr_term' : 'PFI0', # counter terminal for counting
        'MeasCycles' : 20,
        }






name_string = 'monochromator_sweep'
# Create a measurement object m with a name we just made indicating the
# power it's taken at.
m = SiC_Monochromator_Master(name_string)



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
    print 'Measurement completed!'
else:
    print 'Measurement aborted!'


# m.finish() will close the HDF5 and end the measurement.
m.finish()


