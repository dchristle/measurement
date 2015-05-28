# Sweep HWP and measure power from FieldMaster


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


    st0 = qt.instruments['Standa0']


    ni63 = qt.instruments['NIDAQ6363']


    qt.mstart()

    data = qt.Data(name='hwp_sweep')
    data.add_coordinate('Step Position')
    data.add_value('Counts')
    data.create_file()
    filename=data.get_filepath()[:-4]

    plot2d_1 = qt.Plot2D(data, name='measure2D', valdim=1)


    for s in s_vec:
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        st0.move(int(s))
        logging.debug('step set: %s' % (s))
        n = 0
        qt.msleep(0.2)
        while n < 50:
            cur_pos = st0.get_position()
            if cur_pos == s:
                break
            else:
                n = n + 1
                qt.msleep(0.5)

        qt.msleep(0.5)
        counts = ni63.get('ctr0')
        data.add_data_point(s,counts)



    plot2d_1.save_png(filename+'.png')
    data.close_file()
    qt.mend()
