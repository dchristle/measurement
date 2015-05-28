import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
import pylab as plt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
import gc
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

# Defines a class called SiC_PL6PH_Master, DC 2015

class SiC_PL6PHDouble_Master(m2.Measurement):
    # This prefix is used later for filenames; just indicates that it's a PicoHarp
    # measurement on PL6
    mprefix = 'pl6phdouble'
    # Define a 'prepare' function, which prepares the setup to proceed with
    # a measurement by doing things like defining references to the instruments
    # it will use, checking the SNSPD, and setting the initial monochromator
    # wavelength.

    def prepare(self):
        # you indicate that a measurement is about to start and other
        # processes should stop (like batterycheckers, or temperature
        # monitors)
        qt.mstart()
        self._scan_on = True
        # Set up some instruments
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._ph = qt.instruments['ph']
        self._mc = qt.instruments['mc']
        if self.params['cycle_waveforms']:
            self._awg = qt.instruments['awg2']

        # Set the wavelength speed and initial wavelength
        self._mc.set_wavelength_speed(self.params['wavelength_speed'])
        time.sleep(0.2)
        self._mc.set_wavelength(self.params['wavelength'][0])
        current_wl = self._mc.get_wavelength()
        if np.abs(current_wl - self.params['wavelength'][0]) < 1.0:
            print 'Current wavelength is %.3f and within 1 nm, continuing...' % current_wl
        else:
            print 'Current wavelength is %.3f and out of bounds!!'

        # Configure PH for histogram mode, initialize it with the correct settings
        self._ph.start_histogram_mode()

        self._ph.set_Binning(self.params['Binning'])
        self._ph.set_InputCFD0(self.params['CFDLevel0'],self.params['CFDZeroCross0'])
        self._ph.set_InputCFD1(self.params['CFDLevel1'],self.params['CFDZeroCross1'])
        self._ph.set_SyncOffset(self.params['SyncOffset'])
        print 'PicoHarp settings configured.'
        # Set the DAQ counter dwell time, units milliseconds
        # This is how long the counter will count for when we give it a command
        # to tell us the counts.
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        self._ni63.set_ctr0_src(self.params['ctr_term'])


        # Check if the SNSPD is still superconducting
        if self._snspd.check() == False:
            print 'SNSPD went normal and could not restore!'
            self._scan_on = False
            return
        time.sleep(0.1)


        if self.params['cycle_waveforms']:
            self._awg = qt.instruments['awg2']
            self._awg.stop()

            awg_ready = False
            i = 0
            while awg_ready == False and i < 30:
                qt.msleep(2.0)
                if self._awg.get_state() == 'Idle':
                    awg_ready = True
                i = i + 1

            self._awg.set_runmode('SEQ') # set to sequence runmode
            qt.msleep(1.5)
            if self._awg.get_runmode() == 'SEQ':
                print 'AWG successfully set to sequence mode -- OK to proceed.'
                self._awg.start()
                awg_ready = False
                i = 0
                while awg_ready == False and i < 30:
                    qt.msleep(2.0)
                    current_state = self._awg.get_state()
                    if current_state == 'Waiting for trigger' or current_state == 'Running':
                        awg_ready = True
                    else:
                        print 'AWG not ready yet, %d/30' % (i+1)
                    i = i + 1
            else:
                print 'AWG not set to sequence mode! Check if AWG is properly connected to qtlab environment before proceeding.'
                self._scan_on = False
                return

        # Check for counts on PH
        if not self._ph.get_CountRate0() > 0:
            print 'PH 0 not receiving counts!'
            self._scan_on = False
            return
        if not self._ph.get_CountRate1() > 0:
            print 'PH 1 not receiving counts!'
            self._scan_on = False
            return



        return
    # Now define a new method of this particular monochromator measurement class called
    # "measure" that does the actual measuring.
    def measure(self):

        # Start keystroke monitor
        self.start_keystroke_monitor('abort')
        # Take note of the time the measurement starts.
        t0 = time.time()
        signal_0_data = np.zeros(self.params['MeasCycles'], dtype='uint32')
        n_wls = self.params['wavelength'].size
        if self.params['cycle_waveforms']:
            n_wfms = self.params['waveform_list'].size
            average_s_data = np.zeros((n_wfms*n_wls, 65536), dtype='uint32')
            print 'Created array for %d waveforms at %d wavelengths.' % (n_wfms, n_wls)
        else:
            print 'Running with existing waveform -- will not cycle. Just wavelengths.'
            n_wfms = 1
            average_s_data = np.zeros(65536*n_wls, dtype='uint32')

        # Enter the main measurement loop
        for i in np.arange(self.params['MeasCycles']):
            for k in np.arange(n_wls):
                self._mc.set_wavelength(self.params['wavelength'][k])
                qt.msleep(0.5)
                gw = self._mc.get_wavelength()
                print 'Monochromator set to %.3f, indicated wavelength is %.3f' % (self.params['wavelength'][k], gw)
                for j in np.arange(n_wfms):
                    if self._scan_on == False:
                        break
                    if self.params['cycle_waveforms']:
                        qt.msleep(0.1)
                        self._awg.sq_forced_jump(int(self.params['waveform_list']))
                        qt.msleep(0.1)
                        on_correct_waveform = False
                        n = 0
                        while on_correct_waveform == False and n < 10:
                            cur_wfm = int(self._awg.get_sq_position())
                            if cur_wfm == int(self.params['waveform_list']):
                                on_correct_waveform == True
                                break
                            else:
                                self._awg.sq_forced_jump(int(self.params['waveform_list']))
                                qt.msleep(0.5)
                            n = n+1
                            qt.msleep(0.2)
                        if n == 10:
                            print 'Waveform did not change to the correct one in the 10 tries! Stopping measurement.'
                            self._scan_on = False

                        print 'Waveform index %d -- waveform %d/%d, wavelength %d/%d.' % (int(self.params['waveform_list']), j+1, n_wfms, k+1, n_wls)
                    # Now start checking for other issues. If present, stop.
                    # Check if the SNSPD is still superconducting
                    if self._snspd.check() == False:
                        print 'SNSPD went normal and could not restore, breaking.'
                        self._scan_on = False
                        break

                    # Start the PicoHarp acquisition, wait, then retrieve the histogram
                    temp_countA = int(self._ph.get_CountRate1())
                    self._ph.ClearHistMem()
                    self._ph.StartMeas(self.params['AcqTime']*1000) # AcqTime in s, arg in ms
                    print 'Acquiring signal for %s s' % (self.params['AcqTime'])
                    # Wait an extra 0.25 seconds
                    time.sleep(self.params['AcqTime']+0.25)

                    n = 0
                    while self._ph.get_MeasRunning() and n < 10:
                        time.sleep(0.5)
                        n = n + 1
                    if self._ph.get_MeasRunning():
                        print 'Measurement did not finish!'
                        self._scan_on = False
                        break
                    # Retrieve measurement
                    current_data = self._ph.get_Histogram()
                    qt.msleep(0.1)
                    temp_countB = int(self._ph.get_CountRate1())
                    print 'Signal count beginning %d, signal count end %d' % (temp_countA, temp_countB)

                    if self.params['cycle_waveforms']:
                        average_s_data[j+k*n_wfms,:] = average_s_data[j+k*n_wfms,:] + current_data
                        sad = np.sum(average_s_data[j+k*n_wfms,:])
                    else:
                        average_s_data = average_s_data + current_data
                        sad = np.sum(average_s_data)

                    # Average the two counts and save it for later as a diagostic
                    signal_0_data[i] = (temp_countA+temp_countB)/2.0
                    scd = np.sum(current_data)
                    print 'Measured %.1f counts during this acqusition, average sum is %.1f counts, on waveform %d of %d' % (scd, sad, j+1, n_wfms)
                    # Plot the total counts (even though it's mislabeled as esr_avg here)
                    self._average_s_data = average_s_data;
                    self._signal_0_data = signal_0_data;
                    if self.params['cycle_waveforms'] and j == 0:
                        plot2dlog0 = qt.Plot2D(np.log(1.0+np.double(average_s_data[0,:])), name='sicpl6ph_logarithmstart', clear=True)
                        #plot2dlog0 = qt.Plot2D(np.log(1.0+np.double(average_s_data[0+k*n_wfms,:])), name='sicpl6ph_logarithmstart', clear=True)
                    if self.params['cycle_waveforms'] and j == n_wfms-1:
                        plot2dlog1 = qt.Plot2D(np.log(1.0+np.double(average_s_data[1,:])), name='sicpl6ph_logarithmend', clear=True)
                        #plot2dlog1 = qt.Plot2D(np.log(1.0+np.double(average_s_data[1+k*n_wfms,:])), name='sicpl6ph_logarithmend', clear=True)
                    if not self.params['cycle_waveforms']:
                        plot2dlog = qt.Plot2D(np.log(1.0+np.double(average_s_data[0,:])), name='sicpl6ph_logarithm', clear=True)
                    # Check again for an abort
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self._scan_on = False
                        self.stop_keystroke_monitor('abort')
                        break
                    # Now start checking for other issues. If present, stop.
                    # Check if the SNSPD is still superconducting
                    if self._snspd.check() == False:
                        print 'SNSPD went normal and could not restore, breaking.'
                        self._scan_on = False
                        break
                    qt.msleep(0.003) # keeps GUI responsive and checks if plot needs updating.
                if self._scan_on == False:
                        break
            if self._scan_on == False:
                        break






        # Now the measure method of the class is complete, so close the file
        grp = h5.DataGroup('SiCpl6ph', self.h5data, base=self.h5base)
        grp.add('average_signal', data=average_s_data, unit='counts', note='total signal count histogram array')
        if self.params['cycle_waveforms']:
            grp.add('n_wfms', data=n_wfms, note='total number of waveforms')
            grp.add('waveform_list', data=self.params['waveform_list'], note='actual waveform indexes used')
        grp.add('binning', data=self.params['Binning'], note='Power of 2 used for binning in the PicoHarp')
        grp.add('signal0', data=signal_0_data, unit='counts/s', note='countrates measured in each sweep of waveforms')
        # lastly tell the secondary processes (if any) that they are allowed to start again.
        qt.mend()
        return


# Above we defined a new ESR measurement class. Here we're going to create a
# a dictionary of measurement-specific parameters, create some of these ESR
# measurement class objects, and run them.


# measurement parameters - these are all references above.

xsettings = {
        'cycle_waveforms' : True, # set True to cycle waveform_list waveforms
        'waveform_list' : np.array((1)),
        'ctr_term' : 'PFI0',
        'wavelength' : np.array((1036.5,1034.0)), # nm
        'wavelength_speed' : 300.0, # nm/min grating speed
        'dwell_time' : 500.0, # ms
        'ctr_term' : 'PFI0', # counter terminal for counting
        'MeasCycles' : 1500,
        'CFDLevel0' : 125,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 125,
        'CFDZeroCross1' : 10,
        'Binning' : 5,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : 0,
        'AcqTime' : 120, # PicoHarp acquisition time in seconds
        }






name_string = 'PL6dbl_picoharp'
# Create a measurement object m with a name we just made indicating the
# power it's taken at.
m = SiC_PL6PHDouble_Master(name_string)



# Load all the parameters in the slightly modified xsettings dictionary
# into the measurement object 'm' that we just made, which will now have
# the new power
m.params.from_dict(xsettings)

# The if/then here is just leftover from previous code -- since True is
# always True, it will always execute.
if m.review_params():
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


