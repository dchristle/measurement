# PI/Acton monochromator driver via RS-232
#
# David Christle <christle@uchicago.edu>, July 2014


from instrument import Instrument
import visa
import types
import logging
import math
import time
import pyvisa
import numpy as np

class Monochromator(Instrument):


    def __init__(self, name, address, reset=False):
        Instrument.__init__(self, name, tags=['physical'])

        self._address = address


        self.add_parameter('wavelength',
            flags=Instrument.FLAG_GETSET,
            type=types.FloatType,
            units='nm')

        self.add_parameter('wavelength_speed',
            flags=Instrument.FLAG_GETSET,
            type=types.FloatType,
            minval=1)

        self.add_parameter('grating',
            flags=Instrument.FLAG_GETSET,
            type=types.IntType,
            minval=1)


        self.add_function('buffer_clear')




        self._open_serial_connection()
        if reset:
            self.reset()
        else:
            self.get_all()


    # Open serial connection
    def _open_serial_connection(self):
        logging.debug(__name__ + ' : Opening serial connection')

        self._visa = pyvisa.visa.SerialInstrument(self._address,
                baud_rate=9600, data_bits=8, stop_bits=1,
                parity=pyvisa.visa.no_parity, term_chars=pyvisa.visa.CR+pyvisa.visa.LF,
                send_end=False,timeout=5)
        # The purpose of the short timeout is so that the buffer_clear()
        # operation that takes place with every command to ensure the proper
        # output doesn't take too long. Each buffer_clear() usually takes one
        # entire timeout period, since most of the time, the buffer is in fact
        # clear.

    # Close serial connection
    def _close_serial_connection(self):
        '''
        Closes the serial connection
        '''
        logging.debug(__name__ + ' : Closing serial connection')
        self._visa.close()

    def buffer_clear(self): # Got this from Zaber code
        navail = pyvisa.vpp43.get_attribute(self._visa.vi, pyvisa.vpp43.VI_ATTR_ASRL_AVAIL_NUM)
        if navail > 0:
            reply = pyvisa.vpp43.read(self._visa.vi, navail)

    def reset(self):
        self._visa.write('*rst')
        time.sleep(3) # Sleep to avoid trying to talk to the device too quickly

    def get_all(self):
        self.get_wavelength()
        self.get_wavelength_speed()

    def do_get_wavelength_speed(self):
        logging.debug(__name__ + 'reading wavelength speed')
        self.buffer_clear()
        tstr = self._visa.ask('?NM/MIN')
        return float(tstr.split(' ')[1])
    def do_set_wavelength_speed(self, wl_speed):
        logging.debug(__name__ + 'setting wavelength speed')
        self.buffer_clear()
        self._visa.write('%.3f NM/MIN' % wl_speed)
        time.sleep(0.5)
        act_speed = self.get_wavelength_speed()
        if np.abs(act_speed-wl_speed) < 1:
            return True
        else:
            print 'Wavelength speed not set properly!'
            return False

    def do_get_grating(self):
        logging.debug(__name__ + 'reading grating')
        self.buffer_clear()
        tstr = self._visa.ask('?GRATING')
        return float(tstr.split(' ')[1])
    def do_set_grating(self, grating):
        logging.debug(__name__ + 'setting grating')
        self.buffer_clear()
        self._visa.write('%d GRATING' % grating)
        time.sleep(15.0)
        self.buffer_clear()
        final_grating = self.get_grating()
        if final_grating == grating:
            return
        else:
            print 'Grating did not set within timeout!'
            return False

    def do_get_wavelength(self):
        logging.debug(__name__ + 'reading wavelength')
        self.buffer_clear()
        tstr = self._visa.ask('?NM')
        return float(tstr.split(' ')[1])

    def do_set_wavelength(self, wavelength):
        logging.debug(__name__ + 'setting wavelength')
        self.buffer_clear()
        wl_speed = self.get_wavelength_speed()
        cur_wl = self.get_wavelength()
        tgt_wl = wavelength
        est_time = np.abs(cur_wl-tgt_wl)/wl_speed*60.0
        self._visa.write('%.3f NM' % wavelength)
        time.sleep(est_time + 2.0)
        self.buffer_clear()
        final_wl = self.get_wavelength()
        if np.abs(final_wl-tgt_wl) < 0.1:
            return
        else:
            return False
