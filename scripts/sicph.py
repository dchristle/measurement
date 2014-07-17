import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt

class SiCPH_Master(m2.Measurement):

    mprefix = 'antibunching'

    def prepare(self):
        # Set up some instruments
        self._ph = qt.instruments['ph']
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._xps = qt.instruments['xps']

        # Prepare instruments for measurement and verify FBL output

        # Configure PH for histogram mode, initialize
        self._ph.start_histogram_mode()

        self._ph.set_Binning(self.params['Binning'])
        self._ph.set_InputCFD0(self.params['CFDLevel0'],self.params['CFDZeroCross0'])
        self._ph.set_InputCFD1(self.params['CFDLevel1'],self.params['CFDZeroCross1'])
        self._ph.set_SyncOffset(self.params['SyncOffset'])
        print 'PicoHarp settings configured...'
        self._fbl.optimize()
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))

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
        T_save = 10
        Nsaves = np.ceil(self.params['MeasCycles']/float(T_save))
        intermediate_b_data = np.zeros( (Nsaves, 65536), dtype='uint32')
        intermediate_s_data = np.zeros( (Nsaves, 65536), dtype='uint32')
        average_b_data = np.zeros(65536, dtype='uint32')
        average_s_data = np.zeros(65536, dtype='uint32')
        background_0_data = np.zeros(self.params['MeasCycles'], dtype='uint32')
        background_1_data = np.zeros(self.params['MeasCycles'], dtype='uint32')
        signal_0_data = np.zeros(self.params['MeasCycles'], dtype='uint32')
        signal_1_data = np.zeros(self.params['MeasCycles'], dtype='uint32')
        signal_0_data_daq = np.zeros(self.params['MeasCycles'], dtype='uint32')
        background_0_data_daq = np.zeros(self.params['MeasCycles'], dtype='uint32')
        temperature_data = np.zeros(self.params['MeasCycles'], dtype='float')
        heater_data = np.zeros(self.params['MeasCycles'], dtype='float')


        for i in range(self.params['MeasCycles']):
            # Optimize
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'
            # Get signal counts in both channels
            time.sleep(0.25)
            self._ni63.set_count_time(1.0)
            signal_0_data_daq[i] = self._ni63.get_ctr0()
            signal_0_data[i] = int(self._ph.get_CountRate0())
            signal_1_data[i] = int(self._ph.get_CountRate1())
            self._ph.ClearHistMem()
            self._ph.StartMeas(self.params['AcqTime']*1000) # AcqTime in s, arg in ms
            print 'Acquiring signal for %s s' % (self.params['AcqTime'])
            time.sleep(self.params['AcqTime']+0.25)

            n = 0
            while self._ph.get_MeasRunning() and n < 10:
                time.sleep(0.5)
                n = n + 1
            if self._ph.get_MeasRunning():
                print 'Measurement did not finish!'
                break
            # Retrieve measurement
            current_data = self._ph.get_Histogram()
            temp_count = self._ni63.get_ctr0()
            print 'Signal count beginning %d, signal count end %d' % (signal_0_data_daq[i], temp_count)
            signal_0_data_daq[i] = (signal_0_data_daq[i] + temp_count)/2.0
            average_s_data = average_s_data + current_data

##            if i == 0:
##                intermediate_s_data[0,:] = np.transpose(average_s_data)
##            elif np.mod(i,20) == 0:
##                intermediate_temp_data = np.zeros( (i/20+1,65536), dtype='uint32')
##                intermediate_temp_data[:-1,:] = intermediate_s_data
##                intermediate_temp_data[i/20,:] = average_s_data
##                intermediate_s_data = np.copy(intermediate_temp_data)
            if np.mod(i,T_save) == 0:
                intermediate_s_data[i/T_save,:] = average_s_data
            # Optimize
            self._ph.ClearHistMem()
            if self._fbl.optimize() == False:
                if self._fbl.optimize() == False:
                    print 'FBL failed twice, breaking.'

            cur_X = self._fsm.get_abs_positionX()
            cur_Y = self._fsm.get_abs_positionY()
            bg_cts = 0.0
            self._ni63.set_count_time(1.0)
            self._fsm.set_abs_positionX(cur_X + self.params['pos_displacement'])
            bg_cts = bg_cts + float(self._ni63.get_ctr0())
            self._fsm.set_abs_positionX(cur_X - self.params['pos_displacement'])
            bg_cts = bg_cts + float(self._ni63.get_ctr0())
            self._fsm.set_abs_positionX(cur_X)
            self._fsm.set_abs_positionY(cur_Y - self.params['pos_displacement'])
            bg_cts = bg_cts + float(self._ni63.get_ctr0())
            self._fsm.set_abs_positionY(cur_Y + self.params['pos_displacement'])
            bg_cts = bg_cts + float(self._ni63.get_ctr0())
            self._fsm.set_abs_positionY(cur_Y)
            print 'Displacing to %s' % (cur_X + self.params['pos_displacement'])
            self._fsm.set_abs_positionX(cur_X + self.params['pos_displacement'])
            # Get background counts in both channels
            time.sleep(0.5)

            background_0_data_daq[i] = bg_cts/4.0
            background_0_data[i] = int(self._ph.get_CountRate0())
            background_1_data[i] = int(self._ph.get_CountRate1())
            avg_bck = float(background_0_data[i]) + float(background_1_data[i])
            avg_sig = float(signal_0_data[i]) + float(signal_1_data[i])
            print 'avg sig %.2f, avg bck %.2f or %.2f 4 point back, S/B %.2f (new S/B is %.2f)' % (avg_sig/2.0, avg_bck/2.0, bg_cts/4.0, (avg_sig-avg_bck)/avg_bck, ((signal_0_data_daq[i]-(bg_cts/4.0))/(bg_cts/4.0)))
            # Now start the background time correlated data acquisition
            #self._ph.StartMeas(self.params['BackTime']*1000) # AcqTime in s, arg in ms
            #print 'Acquiring background for %.2f s, iteration %s of %s' % (self.params['BackTime'], i+1, self.params['MeasCycles'])
            #time.sleep(self.params['BackTime']+0.25)
##            n = 0
##            while self._ph.get_MeasRunning() and n < 10:
##                time.sleep(0.5)
##                n = n + 1
##            if self._ph.get_MeasRunning():
##                print 'Measurement did not finish!'
##                break
            # Retrieve measurement
            #current_data = self._ph.get_Histogram()

            average_b_data = average_b_data + current_data

##            if i == 0:
##                intermediate_b_data[0,:] = average_b_data
##            elif np.mod(i,20) == 0:
##                intermediate_temp_data = np.zeros( (i/20+1,65536), dtype='uint32')
##                intermediate_temp_data[:-1,:] = intermediate_b_data
##                intermediate_temp_data[i/20,:] = average_b_data
##                intermediate_b_data = np.copy(intermediate_temp_data)
            if np.mod(i,T_save) == 0:
                intermediate_b_data[i/T_save,:] = average_b_data

            print 'Returning to spot.'
            self._fsm.set_abs_positionX(cur_X)

            plot2d = qt.Plot2D(average_s_data, name='measureab', clear=True)
            temperature_data[i] = self._ls332.get_kelvinA()
            heater_data[i] = self._ls332.get_heater_output()
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

        N_completed = np.floor(i/float(T_save))
        intermediate_b_data = np.copy(intermediate_b_data[0:N_completed,:])
        intermediate_s_data = np.copy(intermediate_s_data[0:N_completed,:])
        # Start saving data
        grp = h5.DataGroup('SiCphdata', self.h5data, base=self.h5base)
        grp.add('s', data=average_s_data, unit='counts', note='total s counts')
        grp.add('b', data=average_b_data, unit='counts', note='total b counts')
        grp.add('intermediate_s_data', data=intermediate_s_data, unit='counts', note='intermediate signal count data')
        grp.add('intermediate_b_data', data=intermediate_b_data, unit='counts', note='intermediate background count data')
        grp.add('signal0daq', data=signal_0_data_daq, unit='counts/s', note='signal counts in channel 0, measured by daq')
        grp.add('signal0', data=signal_0_data, unit='counts/s', note='signal counts in channel 0')
        grp.add('signal1', data=signal_1_data, unit='counts/s', note='signal counts in channel 1')
        grp.add('background0daq', data=background_0_data_daq, unit='counts/s', note='background counts in channel 0, measured by daq')
        grp.add('background0', data=background_0_data, unit='counts/s', note='background counts in channel 0')
        grp.add('background1', data=background_1_data, unit='counts/s', note='background counts in channel 1')
        grp.add('temperature', data=temperature_data, unit='K', note='temperature')
        grp.add('heater_output', data=heater_data, unit='percent', note='heater output power in percent')

        return

##    def rebin(self, data, binsize):
##        Nbin = np.size(data)/binsize
##        data_reb = np.zeros(Nbin,dtype='uint32')
##        for i in range(Nbin):
##            for j in range(binsize):
##                data_reb[i] = data_reb[i] + data[i*binsize + j]


# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'temperature_tolerance' : 2,
        'pos_displacement' : -1.25,
        'fbl_time' : 55.0,
        'CFDLevel0' : 125,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 400,
        'CFDZeroCross1' : 10,
        'Binning' : 0,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : -84000,
        'AcqTime' : 70,
        'BackTime' : 15,
        'MeasCycles' : 1000
        }


# Create a measurement object m
m = SiCPH_Master('PL4')

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
ea_t = qt.instruments['ea']
ls332_t = qt.instruments['ls332']
cur_temp = ls332_t.get_kelvinA()
msg_string = 'Antibunching PL4 measurement stopped %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)
# Now enter a holding pattern of continuous tracking.
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