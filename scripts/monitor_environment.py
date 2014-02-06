# Monitor temperature and humidity over optical table


import numpy
from time import time,sleep
import os
import qt
import logging
import msvcrt
import datetime

def simple(t_maximum, meas_period):

    # This example will monitor the T/H for a period of up to
    # t_maximum, which is measured in HOURS.



    thmon = qt.instruments['thmon']

    # meas_period is the measurement period denoted in seconds.





    qt.mstart()

    data = qt.Data(name='environment_monitor')
    data.add_coordinate('Time (minutes)')
    data.add_value('Temperature')
    data.add_value('Humidity')
    data.create_file()
    filename=data.get_filepath()[:-4]


    plot2d_1 = qt.Plot2D(data, name='env_temperature', valdim=1)
    plot2d_2 = qt.Plot2D(data, name='env_humidity', valdim=2)

    starttime = datetime.datetime.now()


    timedelta = 0
    while timedelta < t_maximum*3600:
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): break
        currentT = thmon.get_probe1_temperature()*9/5+32
        currentH = thmon.get_humidity()
        currenttime = datetime.datetime.now()
        c = currenttime - starttime
        timedelta = (c.days * 86400 + c.seconds)
        qt.msleep(meas_period)
        data.add_data_point(timedelta/60.0,currentT,currentH)




    plot2d_1.set_plottitle('Environment temperature monitoring, SiC SD Lab, ' + starttime.strftime("%A, %d. %B %Y %I:%M%p"))
    plot2d_2.set_plottitle('Environment humidity monitoring, SiC SD Lab, ' + starttime.strftime("%A, %d. %B %Y %I:%M%p"))

    plot2d_1.save_png(filename+'.png')
    data.close_file()
    qt.mend()
