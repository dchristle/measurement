# XPS Sweep back and forth script
# David Christle <christle@uchicago.edu>, 2013
#
# This script just sweeps the FSM back and forth. The FSM should already be
# intialized as a Qt instrument.
#

import time
import logging
import qt
import numpy
import msvcrt

def main(channel = 'Z'):
    ##channel = 'X' # This can be 'X' or 'Y' or 'Z' channels
    min_position = -7.0 # mm
    max_position = 7.4 # m
    rate = 0.5 # number of points written per second to the stage
    density = 2 # number of points across the full scale range
    wait =  6 # Wait time before sweeping back in seconds

    x_mm_array_f = numpy.linspace(min_position,max_position,density)
    x_mm_array_b = numpy.linspace(max_position,min_position,density)
    xps = qt.instruments['xps']

    while 1:
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        print 'Sweeping %s forward...' % channel
        for i in x_mm_array_f:
            xps.set('abs_position%s' % channel, i)
            time.sleep(1.0/rate)
        time.sleep(wait)
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        print 'Sweeping %s backward...' % channel
        for i in x_mm_array_b:
            xps.set('abs_position%s' % channel, i)
            time.sleep(1.0/rate)
        time.sleep(wait)


    f_string = 'set_abs_position' + channel
    # Create function that corresponds to the abs_position set function of
    # the correct channel, then call it to reset the voltage to 0.
    reset_channel = getattr(xps, f_string)
    reset_channel(0.0)
    return

if __name__ == '__main__':
    main(channel='Y')
