# QE of SNSPD


import numpy
from time import time,sleep
import os
import qt
import logging
import msvcrt

def simple(s_vec):
    '''
    this example is based on 'measure_module.py'
    you will need to set up a vector of frequencies and then call the
    measurment script.

    To run the function type in the terminal:

    sv = numpy.arange(0,300000,1000)
    esr_meas.simple(sv)
    '''


    snspd = qt.instruments['snspd']
    sr400 = qt.instruments['sr400']




    qt.mstart()

    data = qt.Data(name='QE calculation, channel 0')
    data.add_coordinate('Bias V')
    data.add_value('Counts')
    data.create_file()
    filename=data.get_filepath()[:-4]

    snspd.set_bias1(0)
    plot2d_1 = qt.Plot2D(data, name='measure2D', valdim=1)

    for s in s_vec:
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        snspd.set_bias1(s)
        logging.debug('step set: %s' % (s))
        print 'bias at %s' % snspd.get_bias0()
        n = 0
        qt.msleep(1.2)
        counts = sr400.get_counterA()


        qt.msleep(0.5)

        data.add_data_point(s,counts)


    data.new_block()


    plot2d_1.save_png(filename+'.png')
    data.close_file()
    qt.mend()
    snspd.set_bias1(0)