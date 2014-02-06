# Monitor cooldown of SNSPD detector


import numpy
from time import time,sleep
import os
import qt
import logging
import msvcrt

def simple(t_maximum, meas_period):

    # This example will monitor the SNSPD cooldown for a period of up to
    # t_maximum, which is measured in HOURS.



    ls211 = qt.instruments['ls211']
    snspd = qt.instruments['snspd']

    # meas_period is the measurement period denoted in seconds.





    qt.mstart()

    data = qt.Data(name='snspd_cooldown_monitor')
    data.add_coordinate('Measurement Index')
    data.add_value('Temperature')
    data.add_value('Channel 0 Voltage')
    data.add_value('Channel 1 Voltage')
    data.create_file()
    filename=data.get_filepath()[:-4]


    plot2d_1 = qt.Plot2D(data, name='measureSNSPDT', valdim=1)
    plot2d_2 = qt.Plot2D(data, name='measureSNSPDVs', valdim=2)
    plot2d_2.add_data(data, valdim=3)
    meas_idx = 0

    # Set bias on detectors
    snspd.set_bias0(0.2)
    snspd.set_bias1(0.2)
    while True:
        if (meas_idx > t_maximum*3600/15): break
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        currentT = ls211.get_temperature()
        currentV0 = snspd.get_vmeas0()
        currentV1 = snspd.get_vmeas1()

        qt.msleep(meas_period)
        data.add_data_point(meas_idx,currentT,currentV0,currentV1)
        meas_idx = meas_idx + 1



    plot2d_1.set_plottitle('SNSPD cooldown, measurement period %s seconds, measured for %s seconds total' % (meas_period, meas_period*meas_idx))

    plot2d_1.save_png(filename+'.png')
    data.close_file()
    qt.mend()
