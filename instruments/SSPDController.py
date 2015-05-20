# SSPDController.py, virtual instrument for SSPD measurements
# Reinier Heeres <reinier@heeres.eu>, 2009 - 2010
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

from instrument import Instrument
import types
import qt
import time

class SSPDController(Instrument):

    def __init__(self, name, ni_ins, resistance=500):
        Instrument.__init__(self, name, tags=['virtual'])

        self._ni_ins = ni_ins
        if self._ni_ins is not None:
            self._ni_ins.set('ctr0_src', 'PFI0')
            self._ni_ins.set('ctr1_src', 'PFI1')

        self.add_parameter('bias', type=types.FloatType, channels=[0, 1],
                flags=Instrument.FLAG_SET | Instrument.FLAG_SOFTGET,
                minval=-10, maxval=10, units='V',
                )
        self.add_parameter('vmeas', type=types.FloatType, channels=[0, 1],
                flags=Instrument.FLAG_GET, units='V',
                )
        self.add_parameter('thresh', type=types.IntType, channels=[0, 1],
                flags=Instrument.FLAG_SET | Instrument.FLAG_SOFTGET,
                minval=0, maxval=31,
                )

        # Read counts of one or both channels
        self.add_parameter('counts0', type=types.IntType,
                flags=Instrument.FLAG_GET)
        self.add_parameter('counts1', type=types.IntType,
                flags=Instrument.FLAG_GET)
        self.add_parameter('counts', type=types.ListType,
                flags=Instrument.FLAG_GET)

        self.add_parameter('resistance', type=types.FloatType,
                flags=Instrument.FLAG_SET | Instrument.FLAG_SOFTGET,
                units='kOhm')

        self.add_parameter('zerosignal', type=types.FloatType,
                flags=Instrument.FLAG_SET | Instrument.FLAG_SOFTGET,
                minval=0.000, maxval=0.5, units='V',
                doc='Absolute voltage level to consider zero')

        self.add_parameter('inttime', type=types.FloatType,
                flags=Instrument.FLAG_SET | Instrument.FLAG_SOFTGET,
                minval=0, maxval=1e3, units='sec',
                doc='Integration time in seconds')

        self.set_inttime(0.2)
        self.set('bias0', 0)
        self.set('bias1', 0)
        self.set('thresh0', 0x12)
        self.set('thresh1', 0x12)
        self.set_zerosignal(0.01)
        self.set_resistance(resistance)

        self.add_function('check0')
        self.add_function('check1')
        self.add_function('check')
#        self.add_function('iv')
#        self.add_function('standard_iv')
        self.add_function('iv_counts0')
        self.add_function('iv_counts1')
        self.add_function('reset_bias')
        self.add_function('ch0on')
        self.add_function('ch0off')
        self.add_function('ch1on')
        self.add_function('ch1off')
        self.add_function('zero')
        self.add_function('isenabled')

    def do_set_bias(self, val, channel=0, check=False):
        if self._ni_ins is None:
            return
        self._bias = val
        self._ni_ins.set('ao%d' % channel, val)
        if check:
            self.check(channel)

    def do_get_vmeas(self, channel=0):
        if self._ni_ins is None:
            return
        return self._ni_ins.get('ai%d' % channel)

    def do_set_thresh(self, val, channel=0):
        if self._ni_ins is None:
            return
        if channel == 0:
            lines = 'port0/line0:7'
        else:
            lines = 'port0/line8:15'
        self._ni_ins.digital_out(lines, val)
    def isenabled(self, channel):
        ret = np.zeros(2)
        ret[0] = (self.do_get_bias0() > 0.05)
        ret[1] = (self.do_get_bias1() > 0.05)
        return ret

    def do_get_counts0(self):
        if self._ni_ins is None:
            return
        return self._ni_ins.get('ctr0')

    def do_get_counts1(self):
        if self._ni_ins is None:
            return
        return self._ni_ins.get('ctr1')

    def do_get_counts(self):
        if self._ni_ins is None:
            return
        return self._ni_ins.read_counters(['ctr0', 'ctr1'])

    def do_set_resistance(self, val):
        self._resistance = val

    def do_set_inttime(self, val):
        if self._ni_ins is None:
            return
        return self._ni_ins.set_count_time(val)

    def do_set_zerosignal(self, val):
        self._zerosignal = val

    def check_chan(self, channel):
        biaspar = 'bias%d' % channel
        vbias = self.get(biaspar)
        vmeaspar = 'vmeas%d' % channel

        n = 0
        val = self.get(vmeaspar)
        while abs(val) > self._zerosignal and n < 10:
            print 'detector switched normal, restoring...'
            was_super = False
            self.set(biaspar, 0)
            time.sleep(0.5)
            self.set(biaspar, vbias)
            time.sleep(0.5)
            val = self.get(vmeaspar)
            n += 1

        if abs(val) > self._zerosignal:
            print 'Unable to restore detector!'
            return False
        else:
            return True

    def check0(self):
        return self.check_chan(0)

    def check1(self):
        return self.check_chan(1)

    def check(self):
        ret = self.check0() & self.check1()
        return ret

    def iv(self, channel, start=0, stop=3, step=0.2, delay=0.2):
        '''
        Take an IV on channel chan:
        - start / stop in V
        - steps
        '''

        biaspar = 'bias%d' % channel
        vmeaspar = 'vmeas%d' % channel

        r = self.get_resistance() / 1000.0 # Not sure why this is dividing. Bug?
        n = (int(abs(stop - start) / step)) + 1
        data = qt.Data(name='iv')
        data.add_coordinate('Vbias', units='V')
        data.add_value('Vmeas', units='V')
        data.create_file()
        for i in range(n):
            v_out = start + i * step
            current = v_out/r
            self.set(biaspar, v_out, check=False)
            qt.msleep(delay)
            v_meas = self.get(vmeaspar)
            print 'v_out, current_out, v_meas: %f, %f, %f' %(v_out, current, v_meas)
            data.add_data_point(v_out, v_meas)

        self.set(biaspar, 0)
        qt.plot(data, name='iv', clear=True)
        return data


    def iv_counts(self, channel, start=0, stop=7.5, step=0.25, delay=0.2):
        '''
        Take an IV on channel chan and measure count rate:
        - start / stop in V
        - steps
        '''

        print ''

        biaspar = 'bias%d' % channel
        vmeaspar = 'vmeas%d' % channel
        ni6216 = qt.instruments['NIDAQ6216'] # Get internal DAQ instrument
        ni6216.set_count_time(1.0)
        r = self.get_resistance() / 1000.0
        n = (int(abs(stop - start) / step)) + 1
        data = qt.Data(name='iv')
        data.add_coordinate('Vbias', units='V')
        data.add_value('Vmeas', units='V')
        data.add_value('Counts', units='')
        data.create_file()
        for i in range(n):
            v_out = start + i * step
            current = v_out/r
            self.set(biaspar, v_out, check=False)
            time.sleep(delay)
            v_meas = self.get(vmeaspar)
            if channel == 0:
                counts = self.get_counts0()
            else:
                counts = self.get_counts1()
            print 'v_out, current_out, v_meas, counts: %f, %f, %f, %d' % (v_out, current, v_meas, counts)
            data.add_data_point(v_out, v_meas, counts)

            if v_meas > 0.5:
                break

        self.set(biaspar, 0)

        qt.plot(data, name='ivcounts', clear=True)
        qt.plot(data, name='ivcounts', valdim=2, right=True)

        return data

    def iv_counts_sr400(self, channel, start=0, stop=10, step=0.25, delay=0.2):
        '''
        Take an IV on channel chan and measure count rate:
        - start / stop in V
        - steps
        '''

        print ''
        sr400 = qt.instruments['sr400']

        biaspar = 'bias%d' % channel
        vmeaspar = 'vmeas%d' % channel

        r = self.get_resistance() / 1000.0
        n = (int(abs(stop - start) / step)) + 1
        data = qt.Data(name='iv')
        data.add_coordinate('Vbias', units='V')
        data.add_value('Vmeas', units='V')
        data.add_value('Counts', units='')
        data.create_file()
        for i in range(n):
            v_out = start + i * step
            current = v_out/r
            self.set(biaspar, v_out, check=False)
            time.sleep(delay)
            v_meas = self.get(vmeaspar)
            if channel == 0:
                counts = sr400.get_countA()
            else:
                counts = sr400.get_countB()
            print 'v_out, current_out, v_meas, counts: %f, %f, %f, %d' % (v_out, current, v_meas, counts)
            data.add_data_point(v_out, v_meas, counts)

            if v_meas > 0.5:
                break

        self.set(biaspar, 0)

        qt.plot(data, name='ivcounts', clear=True)
        qt.plot(data, name='ivcounts', valdim=2, right=True)

        return data

    def iv_counts_disc_sr400(self, channel, start=0, stop=10, step=0.25, delay=0.2, discarray = [0.2]):
        '''
        Take an IV on channel chan and measure count rate:
        - start / stop in V
        - steps
        '''

        print ''
        sr400 = qt.instruments['sr400']

        biaspar = 'bias%d' % channel
        vmeaspar = 'vmeas%d' % channel



        r = self.get_resistance() / 1000.0
        n = (int(abs(stop - start) / step)) + 1
        data = qt.Data(name='iv')
        data.add_coordinate('Vbias', units='V')
        data.add_value('Vmeas', units='V')
        data.add_value('Counts', units='')
        data.create_file()

        plot3d = qt.Plot3D(data, name='3D', style='image', coorddims=(1,0), valdim=2)
        for disc in discarray:
            if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
            sr400.set_disc_levelA(disc)
            sr400.set_disc_levelB(disc)
            for i in range(n):
                if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
                v_out = start + i * step
                current = v_out/r
                self.set(biaspar, v_out, check=False)
                time.sleep(delay)
                v_meas = self.get(vmeaspar)
                if channel == 0:
                    counts = sr400.get_countA()
                else:
                    counts = sr400.get_countB()
                print 'v_out, current_out, v_meas, counts: %f, %f, %f, %d' % (v_out, current, v_meas, counts)
                data.add_data_point(v_out, v_meas, counts)
                data.add_data_point(v_out, disc, counts)

                if v_meas > 0.5:
                    break

        self.set(biaspar, 0)

        qt.plot(data, name='ivcounts', clear=True)
        qt.plot(data, name='ivcounts', valdim=2, right=True)

        return data
    def iv_counts0(self, ret=False):
        val = self.iv_counts(0)
        if ret:
            return val

    def iv_counts1(self, ret=False):
        val = self.iv_counts(1)
        if ret:
            return val

    def standard_iv(self, channel):
        d = self.iv(channel, 0, 3, 0.2, 0.2)
        return d

    def ch0on(self):
        # Constant for this channel,determined by D. Christle 2014/01/31
        # Quantum efficiency: 23.3% @ 1060 nm
        # Dark counts: 5 per second
        default_bias_ch0 = 3.45
        self.set_bias0(default_bias_ch0)
        print 'SNSPD channel 0 enabled, bias set to %s V' % default_bias_ch0
        print 'Checking for superconducting state...'
        time.sleep(2)
        result = self.check()
        print 'Superconducting state: %s' % result
        return result

    def ch1on(self):
        # Constant for this channel,determined by D. Christle 2014/01/31
        # Quantum efficiency: 24.7% @ 1060 nm
        # Dark counts: 20 per second
        default_bias_ch1 = 4.75
        self.set_bias1(default_bias_ch1)
        print 'SNSPD channel 1 enabled, bias set to %s V' % default_bias_ch1
        print 'Checking for superconducting state...'
        time.sleep(2)
        result = self.check()
        print 'Superconducting state: %s' % result
        return result

    def ch0off(self):
        aa = self.set_bias0(0.0)
        print 'SNSPD channel 0 disabled'
        return aa

    def ch1off(self):
        aa = self.set_bias1(0.0)
        print 'SNSPD channel 1 disabled'
        return aa

    def zero(self):
        aa = self.set_bias0(0.0)
        bb = self.set_bias1(0.0)
        print 'Both SNSPD channels disabled'
        return (aa and bb)

    def reset_bias(self):
        b0 = self.get_bias0()
        b1 = self.get_bias1()
        self.set_bias0(0)
        self.set_bias1(0)
        time.sleep(1)
        self.set_bias0(b0)
        self.set_bias1(b1)

