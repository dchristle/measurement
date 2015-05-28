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
piezo_vec = np.linspace(-12,12.5,150)
curr_gain = np.linspace(0*10**(-3), 5*10**(-3),12)
schr = qt.instruments['schru']
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
data.add_coordinate('piezo voltage (V)')
data.add_coordinate('current gain (mA/V)')
data.add_value('frequency (GHz)')

# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.
data.create_file()

# Next two plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.
plot2d = qt.Plot2D(data, name='fp_freq', coorddim=0, valdim=2, maxtraces=20)
schr.set_piezo_offset(piezo_vec[0])
# Now attempt to find a current where the mode is stable
curr_diff = np.linspace(0,25*10**(-3),25)
curr0 = schr.get_current()
for cdf in curr_diff:
    schr.set_current(curr0+cdf)
    time.sleep(3)
    curr_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
    if np.size(curr_peaks) < 3:
        print 'Not enough peaks to start at current %.4f' % (curr0+cdf)
    else:
        print 'Peaks found at current %.4f, starting.' % (curr0+cdf)
        break




# preparation is done, now start the measurement.
# It is actually a simple loop.
scan_on = True
for cg in curr_gain:
    schr.set_current_coupling_gain(cg)

    time.sleep(1.00)
    curr_freq = 0
    for pv in piezo_vec:
        schr.set_piezo_offset(pv)
        time.sleep(1.0)
        curr_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
        if np.size(curr_peaks) < 3:
            print 'Not enough peaks to start!'
        else:
            print 'Peaks found at piezo %.2f, starting.' % pv
            break
    for pv in piezo_vec:
        schr.set_piezo_offset(pv)
        time.sleep(0.5)
        new_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
        time.sleep(0.1)
        #fp.read_sweep_plot(500, 10000, 'ai1')

        if np.size(new_peaks) == 3 or np.size(new_peaks) == 4:
            df = fp.delta_freq(curr_peaks,new_peaks)
            curr_freq = curr_freq + df
            curr_peaks = new_peaks
            data.add_data_point(pv,cg*1000.0,curr_freq)
            print 'peaks: %s, df %.3f' % (new_peaks, df)
        else:
            print 'peaks: %s' % new_peaks




        if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q":
                        scan_on = False
                        break
    data.new_block()
    if msvcrt.kbhit() or scan_on == False:
                    kb_char=msvcrt.getch()
                    if kb_char == "q" or scan_on == False: break


# after the measurement ends, you need to close the data file.
data.close_file()
# lastly tell the secondary processes (if any) that they are allowed to start again.
qt.mend()
