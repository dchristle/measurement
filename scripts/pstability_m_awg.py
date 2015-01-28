import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

class SiC_Photostability_Master(m2.Measurement):

    mprefix = 'photostability'

    def sequence(self, upload = True):
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

        elements = []
        # Create waveform that has laser, microwaves, photon counting, and 1/0 I/Q modulation on
        # all the time for a long period of time (~100 us).
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='laser')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)

        # create a sequence from the pulses -- only one in this case
        seq = pulsar.Sequence('CW ESR Sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(2.0)
        # program the AWG
        qt.pulsar.program_sequence(seq)

    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.1)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.2)

        if q >= 20:
            print 'AWG did not jump to proper waveform!'
        return
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
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        self._awg.start()
        time.sleep(3.0)
        self._awg.sq_forced_jump(1)
        time.sleep(3.0)
        self.awg_confirm(1)


        self._fbl.optimize()
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))
        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'


        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))




        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counter prepared.'




        return
    def measure(self):
        print 'Measuring counts versus time...'
        time.sleep(0.2)
        self._fbl.optimize()
        # Get the current counter time
        prev_count_time = self._ni63.get_count_time()
        # Now set the new count time to 1/ the rate we want to count at
        self._ni63.set_count_time(1.0/self.params['readout_rate'])
        # Create voltage array by doing the conversion
        dummy_aochannel = 'ao3' # using this ao channel, because it's unused
        # Figure out how many points we need to sample
        n_points = int(np.round(self.params['readout_rate']*self.params['readout_time']))
        # We're going to just write a bunch of 0's at the prescribed rate
        x_V_array = np.zeros(n_points)
        # Execute write and count function with raw voltage values

        carray = self._ni63.write_and_count(x_V_array,
            dummy_aochannel,'ctr0')
        self._ni63.set_count_time(prev_count_time)

        diff_counts = np.diff(carray)
        qt.Plot2D(range(n_points-1),diff_counts)
        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Photostability_data', self.h5data, base=self.h5base)
        grp.add('readout_rate', data=self.params['readout_rate'], unit='Hz', note='frequency')
        grp.add('readout_time', data=self.params['readout_time'], unit='s', note='total duration')
        grp.add('counts', data=diff_counts, unit='counts', note='count array')
        return



# measurement parameters

xsettings = {
        'focus_limit_displacement' : 20, # microns inward
        'fbl_time' : 55.0, # seconds
        'readout_rate' : 100.0, # Hz
        'readout_time' : 220.0, # seconds
        'ctr_term' : 'PFI0',
        'temperature_tolerance' : 2.0, # Kelvin
        'random' : 1
        }


# Create a measurement object m
print 'About to proceed -- waiting 5 s for quit (press q to quit)'
time.sleep(5.0)
name_string = 'power %.2f dBm' % (p_array[rr])
m = SiC_Photostability_Master(name_string)
xsettings['readout_time'] = 200.0 # seconds
xsettings['readout_rate'] = 30.0 # Hz
# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)
m.sequence(upload=False)


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

# Alert that measurement has finished
ea_t = qt.instruments['ea']
ls332_t = qt.instruments['ls332']
cur_temp = ls332_t.get_kelvinA()
msg_string = 'Rabi measurement stopped at %s, temperature is %.2f K' % (time.strftime('%c'), cur_temp)
ea_t.email_alert(msg_string)

##ps = qt.instruments['xps']
##ps.set_abs_positionZ(12.0)

track_on = True
fbl_t = qt.instruments['fbl']
track_iter = 0
while track_on == True and track_iter < 50:
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    fbl_t.optimize()
    time.sleep(5.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break