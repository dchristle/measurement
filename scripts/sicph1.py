import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt

class SiCPH_Master(m2.Measurement):

    mprefix = 'SiC_Son'

    def prepare(self):
        # Set up some instruments
        self._ph = qt.instruments['ph']
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']

        # Prepare instruments for measurement and verify FBL output

        # Configure PH for histogram mode, initialize
        self._ph.start_histogram_mode()

        self._ph.set_Binning(self.params['Binning'])
        self._ph.set_InputCFD0(self.params['CFDLevel0'],self.params['CFDZeroCross0'])
        self._ph.set_InputCFD1(self.params['CFDLevel1'],self.params['CFDZeroCross1'])
        self._ph.set_SyncOffset(self.params['SyncOffset'])
        print 'PicoHarp settings configured...'
        self._fbl.optimize()
        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'
        # Check for counts on PH
        if self._ph.get_CountRate0() > 0:
            print 'PH 0 channel receiving counts...'
        else:
            print 'PH 0 not receiving counts!'
        if self._ph.get_CountRate1() > 0:
            print 'PH 1 channel receiving counts...'
        else:
            print 'PH 0 not receiving counts!'
        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        #ret = self._snspd.isenabled()
        #if ret[0] and ret[1] and self._snspd.check():
        #    print 'SNSPD channels enabled, superconducting.'
        #else:
        #    print 'Both SNSPD channels not enabled or not superconducting.'

        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays
        average_b_data = np.zeros(65536, dtype='uint32')
        average_s_data = np.zeros(65536, dtype='uint32')
        for i in range(self.params['MeasCycles']):
            # Optimize
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'

            self._ph.ClearHistMem()
            self._ph.StartMeas(self.params['AcqTime']*1000) # AcqTime in s, arg in ms
            print 'Acquiring signal for %s s' % (self.params['AcqTime'])
            time.sleep(self.params['AcqTime'])

            n = 0
            while self._ph.get_MeasRunning() and n < 10:
                time.sleep(1.5)
                n = n + 1
            if self._ph.get_MeasRunning():
                print 'Measurement did not finish!'
                break
            # Retrieve measurement
            current_data = self._ph.get_Histogram()
            average_s_data = average_s_data + current_data
            # Optimize
            self._ph.ClearHistMem()
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            cur_X = self._fsm.get_abs_positionX()
            print 'Displacing to %s' % (cur_X + self.params['pos_displacement'])
            self._fsm.set_abs_positionX(cur_X + self.params['pos_displacement'])
            self._ph.StartMeas(self.params['AcqTime']*1000) # AcqTime in s, arg in ms
            print 'Acquiring background, iteration %s of %s' % (i+1, self.params['MeasCycles'])
            time.sleep(self.params['AcqTime'])
            n = 0
            while self._ph.get_MeasRunning() and n < 10:
                time.sleep(1.5)
                n = n + 1
            if self._ph.get_MeasRunning():
                print 'Measurement did not finish!'
                break
            # Retrieve measurement
            current_data = self._ph.get_Histogram()
            average_b_data = average_b_data + current_data

            cur_X = self._fsm.get_abs_positionX()
            print 'Returning to spot.'
            self._fsm.set_abs_positionX(cur_X - self.params['pos_displacement'])
            plot2d = qt.Plot2D(average_s_data, name='measureab', clear=True)
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
            if self._ph.get_CountRate0() == 0 or self._ph.get_CountRate1() == 0:
                print 'Count rates are zero on one of the channels, breaking.'
                break


        # Start saving data
        grp = h5.DataGroup('SiCphdata', self.h5data, base=self.h5base)
        grp.add('s', data=average_s_data, unit='counts', note='total s counts')
        grp.add('b', data=average_b_data, unit='counts', note='total b counts')

        return



# measurement parameters

xsettings = {
        'temperature_tolerance' : 3,
        'pos_displacement' : -2.0,
        'fbl_time' : 30.0,
        'CFDLevel0' : 500,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 60,
        'CFDZeroCross1' : 10,
        'Binning' : 0,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : -84000,
        'AcqTime' : 45,
        'MeasCycles' : 345
        }


# Create a measurement object m
m = SiCPH_Master('first measurements')

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