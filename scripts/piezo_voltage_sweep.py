# Piezo offset sweep and measure peak locations

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time




#####################################################
# here is where the actual measurement program starts
#####################################################

# you define two vectors of what you want to sweep. In this case
# a magnetic field (b_vec) and a frequency (f_vec)
piezo_vec = np.linspace(-12.5,12.5,100)
schr = qt.instruments['schr']
fp = qt.instruments['fp']
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
data = qt.Data(name='piezo_sweep')

# Now you provide the information of what data will be saved in the
# datafile. A distinction is made between 'coordinates', and 'values'.
# Coordinates are the parameters that you sweep, values are the
# parameters that you readout (the result of an experiment). This
# information is used later for plotting purposes.
# Adding coordinate and value info is optional, but recommended.
# If you don't supply it, the data class will guess your data format.
data.add_coordinate('piezo voltage')
data.add_value('lowest peak voltage')

# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.
data.create_file()

# Next two plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.
plot2d = qt.Plot2D(data, name='fp1_lowest', coorddim=0, valdim=1)

# preparation is done, now start the measurement.
# It is actually a simple loop.
for pv in piezo_vec:
    schr.set_piezo_offset(pv)
    time.sleep(1.5)
    peaks = fp.read_sweep_peaks(500,10000,'ai1')
    time.sleep(0.1)
    #fp.read_sweep_plot(500, 10000, 'ai1')
    print 'peaks: %s' % peaks
    data.add_data_point(pv,peaks[0])
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False: break


# after the measurement ends, you need to close the data file.
data.close_file()
# lastly tell the secondary processes (if any) that they are allowed to start again.
qt.mend()
