# Power vs relative motor step for Sacher laser

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time
import msvcrt
import datetime


# Be careful with this program -- it's essential you stay as close as possible
# to the range where the laser will lase properly. If you go too far beyond, you
# can cause the metal rod between the motor and the grating to fall out.
#
# A good step size to play around with is about 10,000 steps. I think about 200,000
# to 600,000 steps is the full range of the grating, depending on the model, so
# the process you should employ is to try to make small relative motions with
# a sufficient current applied and try to observe a power increase.

motor_step_amount = 15000
approx_wavelength = 1140
LD_curr = 0.139
schr = qt.instruments['schru']
pm = qt.instruments['pm']

schrm = qt.instruments['schrmu']
schr.set_laser_status('ON')
schr.set_current(LD_curr)
# you indicate that a measurement is about to start and other
# processes should stop (like batterycheckers, or temperature
# monitors)
qt.mstart()

# Next a new data object is made.
# The file will be placed in the folder:
# <datadir>/<datestamp>/<timestamp>_testmeasurement/
# and will be called:
# <timestamp>_testmeasurement.dat
# to find out what 'datadir' is set to, type: qt.config.get('datadir')
data = qt.Data(name='PvsMotorsacher')
data.add_coordinate('Relative motor position')

data.add_value('Power (mW)')

data.create_file()

plot2d_1 = qt.plot(data, name='p vs motor',linewidth=6,clear=True)

pm.set_wavelength(approx_wavelength*1e-9)
print ('Seting the power meter wavelength to %f' % approx_wavelength)

# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.


# Next plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.


#plot3d_1.set_plottitle('Sacher resonant laser calibration ' + datetime.date.today().isoformat())
# Measurement loop
cur_position = 0
meas_go = True
ii = 0
while meas_go == True:
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" :
                    meas_go = False
                    break
    if meas_go == False:
        break



    if cur_position == 0:
        schrm.set_target_position(motor_step_amount,False,True)
        cur_position = motor_step_amount
    if cur_position == motor_step_amount:
        schrm.set_target_position(-1*motor_step_amount,False,True)
        cur_position = 0

    time.sleep(5)
    powrr = pm.get_power()*1000.0

    data.add_data_point(ii,powrr)
    plot2d_1.update()
    ii = ii+1


# after the measurement ends, you need to close the data file.
plot2d_1.save_png()
data.close_file()
# lastly tell the secondary processes (if any) that they are allowed to start again.
qt.mend()
#schr.set_current(0.00)
#schr.set_laser_status('OFF')
