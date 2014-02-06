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
    fm = qt.instruments['fm']
    tl = qt.instruments['tl']

    tl.on()



    qt.mstart()

    data = qt.Data(name='basic_HWP_calibration')
    data.add_coordinate('Step Position')
    data.add_value('Power')
    data.create_file()
    filename=data.get_filepath()[:-4]


    plot2d_1 = qt.Plot2D(data, name='measure2D', valdim=1)
    powers = []
    for s in s_vec:
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        st0.move(s)
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
        power = fm.get_power()
        powers.append(power)
        data.add_data_point(s,power)


    data.new_block()


    powersnp = numpy.array(powers)

    print 'Maximum power measured was:' + '{:.7e}'.format(powersnp.max())
    print 'Setting step to %s' % s_vec[powersnp.argmax()]
    st0.move(s_vec[powersnp.argmax()])
    n = 0
    while n < 50:
        cur_pos = st0.get_position()
        if cur_pos == s_vec[powersnp.argmax()]:
            break
        else:
            n = n + 1
            qt.msleep(2)

    new_power = fm.get_power()
    print 'New power is ' + '{:.7e}'.format(new_power)
    plot2d_1.set_plottitle('Set to step %.05f, max power scanned %.4e, new power %.4e' % (s_vec[powersnp.argmax()], powersnp.max(), new_power) )

    plot2d_1.save_png(filename+'.png')
    data.close_file()
    qt.mend()
    tl.off()