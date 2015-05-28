# Current vs. Power for Sacher laser

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time
import msvcrt
import datetime


#####################################################
# here is where the actual measurement program starts
#####################################################

# you define vectors of what you want to sweep. In this case
# the LD current (LD_curr). We also create the vector of
# wavelenghts (W_vec) that we`re gonna change manually

W_vec = np.linspace(1139.6, 1139.6,1)
LD_curr = np.linspace(0.12,0.18,30)
schr = qt.instruments['schru']
pm = qt.instruments['pm']
schr.set_current(LD_curr[0])
schrm = qt.instruments['schrmu']
schr.set_laser_status('ON')
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
data = qt.Data(name='PvsIsacher')
data.add_coordinate('Wavelength (nm)')
data.add_coordinate('Current (mA)')
data.add_value('Power (mW)')

data.create_file()

plot3d_1 = qt.Plot3D(data, name='3D_1', style='image')



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
meas_go = True
for wav in W_vec:
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
    if meas_go == False:
        break
    #Set the laser wavelength
    pm.set_wavelength(wav*1e-9)
    print ('Seting the laser wavelength to %f' % wav)
    schr.set_current(LD_curr[0])
    #schrm.set_wavelength(wav)
    print 'Wavelength set to %.2f' % wav

    #raw_input("When ready press Enter to continue with the measurement...")


    #Set LD current and collect power values
    for curr in LD_curr:
        if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" :
                    meas_go = False
                    break
        schr.set_current(curr)

        time.sleep(25)
        powrr = pm.get_power()
        print 'Power is %.2f mW' % (powrr*1000.0)
        data.add_data_point(wav,curr,powrr)
        plot3d_1.update()


# after the measurement ends, you need to close the data file.
plot3d_1.save_png()
data.close_file()
# lastly tell the secondary processes (if any) that they are allowed to start again.
qt.mend()
#schr.set_current(0.00)
#schr.set_laser_status('OFF')
