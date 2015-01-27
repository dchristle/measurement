# MiniCircuits variable attenuator instrument
# David Christle <christle@uchicago.edu>, September 2014
# This instrument takes in a configuration of digital out lines on the NI DAQ
# and writes the corresponding high/low values to the attenuator.
from instrument import Instrument
import types
import logging
import time
import qt
import math
import numpy as np




class Variable_Attenuator(Instrument):

    def __init__(self, name):
        Instrument.__init__(self, name)


        self._ni63 = qt.instruments['NIDAQ6363']
        # these values are hardcoded based on the wiring to the DAQ -- could be put into cfgman
        self.digital_lines = {
                    'latch' : 'port0/line1',
                    '0.5db' : 'port0/line3',
                    '1.0db' : 'port0/line2',
                    '2.0db' : 'port0/line6',
                    '4.0db' : 'port0/line4',
                    '8.0db' : 'port0/line5',
                }
                
        self.add_parameter('attenuation',
            type=types.FloatType,
            flags=Instrument.FLAG_SET|Instrument.FLAG_SOFTGET,
            units='dBm',
            format='%.1f',
            minval=0.0, maxval=15.5)
            

        
    def do_set_attenuation(self, atten):
        # multiply by 2 (since the lowest unit is 0.5 dB attenuation) and then round it and
        # convert to a string
        string_bin = self.dec2bin(np.round(atten*2))
        # now zero fill it to a constant length of 5
        string_bin = string_bin.zfill(5)
        #print 'string is %s' % string_bin
        # make a list of the lines, in order
        line_list = [self.digital_lines['8.0db'], self.digital_lines['4.0db'], self.digital_lines['2.0db'], self.digital_lines['1.0db'], self.digital_lines['0.5db']]
        # now turn the latch off
        self._ni63.digital_out(self.digital_lines['latch'],0)
        for i in range(5):
            self._ni63.digital_out(line_list[i],int(string_bin[i]))
        # now re-enable the latch
        self._ni63.digital_out(self.digital_lines['latch'],1)
        time.sleep(0.05)
        # now disable the latch
        self._ni63.digital_out(self.digital_lines['latch'],0)
            
        
        
        
        return True
        
    def dec2bin(self, f):
        if f >= 1:
            g = int(math.log(f, 2))
        else:
            g = -1
        h = g + 1
        ig = math.pow(2, g)
        st = ""    
        while f > 0 or ig >= 1: 
            if f < 1:
                if len(st[h:]) >= 10: # 10 fractional digits max
                       break
            if f >= ig:
                st += "1"
                f -= ig
            else:
                st += "0"
            ig /= 2
        if not st[h:]:
            st = st[:h]
        else:
            st = st[:h] + "." + st[h:]
        if (not st[:h]) and (not st[h:]):
            st = '%s' % 0
        return st