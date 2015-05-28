import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt

# David Christle, June 2014 <christle@uchicago.edu>
#
# This class defines a measurement to study the readout/initialization cycle
# of the divacany defects. The idea is to perform a double pulse measurement
# similar to the one done in Manson et al. PRB 2006. There, they were able to
# wait T1 before the first pulse, and then apply the second pulse in short
# succession. Because I only have a DDG right now, I can't wait for different
# times between the pulses -- I have to delay them both equally.
#
# This isn't so bad if we don't care about trying to understand the thermalized
# population of the spin states, however. This type of measurement will be a
# snapshot of a train of equally spaced pulses. We can still get important
# information, hopefully, about the singlet lifetime and the readout rates.

class SiCRO_Master(m2.Measurement):

    mprefix = 'pulse'

    def prepare(self):
        # Set up some instruments
        self._ph = qt.instruments['ph']
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._ddg = qt.instruments['ddg']

        # Prepare instruments for measurement and verify FBL output

        # Prepare DDG
        self._ddg.set_trig_source('internal')
        # Compute the overall rate of the measurement sequence
        self._trigger_rate = np.round(1.0/(1.0e-9 * self.params['trigger_period_start']),5)
        self._ddg.set_trig_rate(self._trigger_rate)

        self._ddg.set_delayA(0.0)
        self._ddg.set_delayB(1.0/self._trigger_rate-300.0*1.0e-9)
        self._ddg.set_delayC(0.0)
        self._ddg.set_delayD(1.0/self._trigger_rate-300.0*1.0e-9)
        self._ddg.set_delayE(0.0)
        self._ddg.set_delayF(1.0/self._trigger_rate-300.0*1.0e-9)
        self._ddg.set_delayG(0.0)
        self._ddg.set_delayH(1.0/self._trigger_rate-300.0*1.0e-9)

        # Configure PH for histogram mode, initialize
        self._ph.start_histogram_mode()

        self._ph.set_Binning(self.params['Binning'])
        self._ph.set_InputCFD0(self.params['CFDLevel0'],self.params['CFDZeroCross0'])
        self._ph.set_InputCFD1(self.params['CFDLevel1'],self.params['CFDZeroCross1'])
        self._ph.set_SyncDiv(self.params['SyncDiv'])
        self._ph.set_SyncOffset(self.params['SyncOffset'])
        print 'PicoHarp settings configured...'
        self._fbl.optimize()

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


        # Measure the DDG speed of set/get.
        tdi = time.time()
        self._ddg.get_delayD()
        tdd = time.time() - tdi
        print 'DDG timing is approx. %.3f ms, using twice this value for sleeps.' % (tdd*1000.0)
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
        self._ddg.set_delayD(self.params['RF_length']*1.0e-9) # Just for initialization
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
        return
    def measure(self):
        # Wall time
        t0 = time.time()

        # Populate some arrays

        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['trigger_period_end'] - self.params['trigger_period_start'])/self.params['trigger_period_step']))

        print '--Readout meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps)--' % (self.params['trigger_period_start'], self.params['trigger_period_end'], self.params['trigger_period_step'], n_steps)

        TP_lengths = np.linspace(self.params['trigger_period_start'], self.params['trigger_period_end'], n_steps)

        average_s_data = np.zeros((n_steps, 65536), dtype='uint32')
        temperature_data = np.zeros(self.params['MeasCycles'], dtype='float')
        heater_data = np.zeros(self.params['MeasCycles'], dtype='float')
        scan_on = True


        for i in range(self.params['MeasCycles']):

            # Enter the loop for measurement
            t1 = time.time()
            for j in range(n_steps):

                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break


                # Track before measurement
                fbl.optimize()
                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break

                # Set the new trigger period
                current_tp_rate = np.round(1.0/(1.0e-9 * TP_lengths[j]),5)
                self._ddg.set_trig_rate(current_tp_rate)
                time.sleep(self._ddgsleep)
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

                average_s_data[j,:] = average_s_data[j,:] + current_data
                plot2d_0 = qt.Plot2D(range(65536),average_s_data[0,:], name='readout_start_average', clear=True)



            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1

            print 'Total time is %.3f, efficiency of %.2f percent.' % (tt, (n_steps*self.params['AcqTime'])/tt*100.0)




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


        # Start saving data
        grp = h5.DataGroup('SiCphdata', self.h5data, base=self.h5base)
        grp.add('s', data=average_s_data, unit='counts', note='total s counts')
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
        'temperature_tolerance' : 1,
        'CFDLevel0' : 200,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 125,
        'CFDZeroCross1' : 10,
        'Binning' : 6,
        'Offset' : 0,
        'SyncDiv' : 2,
        'SyncOffset' : 0,
        'AcqTime' : 50,
        'MeasCycles' : 1000,
        'trigger_period_start' : 1650, # ns,
        'trigger_period_step' : 50.0, # ns,
        'trigger_period_end' : 2100.0, # ns
        'AOM_delay' : 0, # ns
        'AOM_length' : 1600.0, # ns
        'AOM_amplitude' : 2.5, # V
        'RF_delay' : 00.0, # ns
        'RF_length' : 10.0,
        'RF_amplitude' : 2.5, # V
        'readout_amplitude' : 2.5, #V
        'readout_delay' : 0.0,
        'readout_length' : 210.0, # ns


        }


# Create a measurement object m
m = SiCRO_Master('PL2')

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