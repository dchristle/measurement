#-------------------------------------------------------------------------------
# Sacher.py driver for Sacher Laser
# Author:      David Christle
# Created:     21/04/2014

from instrument import Instrument
import visa
import types
import logging
import numpy
import time
import qt


class Sacher(Instrument):

#----------------------------------------------
# Initialization
#----------------------------------------------

    def __init__(self, name, address, reset = False):

        Instrument.__init__(self, name, tags = ['physical'])

        self._address = address
        self._visa = visa.instrument(self._address)


        # Add functions

        self.add_function('get_all')

        self.add_parameter('current_limit',
            flags = Instrument.FLAG_GETSET,
            type = types.FloatType,
            units = 'A',
            minval=0.0,maxval=0.6)

        self.add_parameter('current',
            flags = Instrument.FLAG_GETSET,
            type = types.FloatType,
            units = 'A',
            minval=0.0,maxval=0.6)

        self.add_parameter('PDcurr',
            flags = Instrument.FLAG_GET,
            type = types.FloatType,
            units = 'A')

        self.add_parameter('CALtable',
            flags = Instrument.FLAG_GET,
            type = types.IntType)

        self.add_parameter('power',
            flags = Instrument.FLAG_GET,
            type = types.FloatType,
            units = 'W')

        self.add_parameter('TEC_status',
            flags = Instrument.FLAG_GETSET,
            type = types.StringType,
            format_map = {
                'OFF': 0,
                'ON': 1,
              })
        self.add_parameter('laser_status',
            flags = Instrument.FLAG_GETSET,
            type = types.StringType,
            format_map = {
                'OFF': 0,
                'ON': 1,
              })

        self.add_parameter('temperature',
            flags = Instrument.FLAG_GETSET,
            type = types.FloatType,
            units = 'C',
            minval=17.0,maxval=25.0)

        self.add_parameter('piezo_status',
            flags = Instrument.FLAG_GETSET,
            type = types.StringType,
            format_map = {
                'OFF': 0,
                'ON': 1,
            })

        self.add_parameter('piezo_offset',
            flags = Instrument.FLAG_GETSET,
            type = types.FloatType,
            units = 'V')

        self.add_parameter('current_coupling',
            flags = Instrument.FLAG_GETSET,
            type = types.IntType,
            format_map = {
                0: 'OFF',
                1: 'ON',
            })
        self.add_parameter('current_coupling_gain',
            flags = Instrument.FLAG_GETSET,
            type = types.FloatType,
            units = 'V')
##        self.add_parameter('TEC_current',
##            flats = Instrument.FLAG_GET,
##            type = types.FloatType,
##            units = 'mA')
        self.add_parameter('PID_P',
            flats = Instrument.FLAG_GETSET,
            type = types.FloatType)
        self.add_parameter('PID_I',
            flats = Instrument.FLAG_GETSET,
            type = types.FloatType)
        self.add_parameter('PID_D',
            flats = Instrument.FLAG_GETSET,
            type = types.FloatType)
        self.add_parameter('current_coupling_direction',
            flags = Instrument.FLAG_GETSET,
            type = types.IntType,
            format_map = {
                0: 'NEG',
                1: 'POS',
            })





#---------------------------------------------
# Class functions
#---------------------------------------------

    def get_all(self):

        self.get_current_limit()
        self.get_current()
        self.get_power()
        self.get_TEC_status()
        self.get_temperature()
        self.get_piezo_status()
        self.get_piezo_offset()
        return


    def do_get_current_limit(self):

        return float(self._visa.ask('L:ILIM?'))

    def do_set_current_limit(self, cl):
        self._visa.write('L:ILIM %.3f' % cl)
        return True


    def do_set_current(self, cur):
        time.sleep(0.1)
        st = self.get_TEC_status()
        time.sleep(0.1)
        if st == 'ON':

            self._visa.write('L:CURR %s' % cur)
            time.sleep(1)
            I = float(self._visa.ask('L:CURR?'))
            print 'LD current is %f' % I

        else:

            print 'Could not set the current because TEC is off!'
            print 'Turn on TEC with .set_TEC_status(1)'

    def on(self):
        self.set_laser_status(1)
        time.sleep(3.2)
        reply = self.get_laser_status()
        if reply == 'ON':
            return
        else:
            logging.error(__name__ + ': laser did not turn on!')
        return
    def off(self):
        self.set_laser_status(0)
        time.sleep(0.2)
        reply = self.get_laser_status()
        if reply == 'OFF':
            return
        else:
            logging.error(__name__ + ': laser did not turn off!')
        return

    def do_get_current(self):

        return float(self._visa.ask('L:CURR?'))


    def do_get_PDcurr(self):

        return float(self._visa.ask('PD:CURR?'))


    def do_get_CALtable(self):

        return int(self._visa.ask('PD:CALT?'))


    def do_get_power(self):

        return float(self._visa.ask('L:POW?'))


    def do_get_TEC_status(self):

        return '%s' % self._visa.ask('TEC:ENA?')


    def do_set_TEC_status(self, TEC_status):

        self._visa.write('TEC:ENA %s' % TEC_status)

    def do_get_laser_status(self):

        return '%s' % self._visa.ask('L:STAT?')


    def do_set_laser_status(self, laser_status):

        self._visa.write('L:STAT %s' % laser_status)

    def do_set_piezo_status(self, piezo_status):

        self._visa.write('P:ENA %s' % piezo_status)

    def do_get_piezo_status(self):

        return self._visa.ask('P:ENA?')

    def do_get_piezo_offset(self):
        response = self._visa.ask('P:OFFS?')
        return float(response)


    def do_set_piezo_offset(self,po):

        self._visa.write(('P:OFFS %.3f' % po))
        return True


    def do_get_current_coupling(self):
        return self._visa.ask('CC:ENA?')


    def do_set_current_coupling(self,cc):
        self._visa.write('CC:ENA %d' % cc)
        return True


    def do_get_current_coupling_gain(self):
        return self._visa.ask('CC:GAIN?')


    def do_set_current_coupling_gain(self, cg):
        self._visa.write('CC:GAIN %.4f' % cg)
        return True


    def do_set_temperature(self, temp):
        self._visa.write('TEC:TEMP %.3f' % temp)
        return


    def do_get_temperature(self):
        return self._visa.ask('TEC:TEMP?')


    def do_get_TEC_current(self):
        return (1000.0*float(self._visa.ask('TEC:CURR?')))


    def do_get_PID_D(self):
        return float(self._visa.ask('TEC:DP?'))


    def do_set_PID_D(self, dp):
        self._visa.write('TEC:DP %s' % dp)
        return True


    def do_get_PID_P(self):
        return float(self._visa.ask('TEC:PP?'))


    def do_set_PID_P(self, pp):
        self._visa.write('TEC:PP %s' % pp)
        return True


    def do_get_PID_I(self):
        return float(self._visa.ask('TEC:IP?'))


    def do_set_PID_I(self, ip):
        self._visa.write('TEC:IP %s' % ip)
        return True


    def do_get_current_coupling_direction(self):
        return self._visa.ask('CC:DIR?')


    def do_set_current_coupling_direction(self, ccd):
        self._visa.write('CC:DIR %d' % ccd)
        return





