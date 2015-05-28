# Piezo offset sweep and measure peak locations

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time
import msvcrt




#####################################################
# here is where the actual measurement program starts
#####################################################

# you define two vectors of what you want to sweep. In this case
# a magnetic field (b_vec) and a frequency (f_vec)
piezo_vec = np.linspace(-0.5,0.5,50)
temp_vec = np.linspace(21.90,22.1,10)
schr = qt.instruments['schr']
fp = qt.instruments['fp']
schr.set_temperature(22)

print 'Attempting initial temperature stabilization.'
while n < 12:
        print 'iteration %d for init. temp.' % n
        time.sleep(10.0)
        t_actual = schr.get_temperature()
        print 'Measured temperature %.3f, desired %.3f' % (t_actual, 22)
        if np.abs(22 - t_actual) < 0.02:
            time.sleep(10.0)
            t_actual = schr.get_temperature()
            if np.abs(22 - t_actual) < 0.02:
                print 'Temperature stable. Proceeding.'
                break
        n = n + 1

qt.mstart()


data = qt.Data(name='temp_sweep')


data.add_coordinate('piezo voltage (V)')
data.add_coordinate('temperature (K)')
data.add_value('frequency (GHz)')


data.create_file()


plot2d = qt.Plot2D(data, name='fp_freqvt', coorddim=0, valdim=2, maxtraces=20, traceofs=2)
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
for tt in temp_vec:
    # Reset the piezo to the first voltage step
    schr.set_piezo_offset(piezo_vec[0])
    # Set the temperature to the desired temperature
    schr.set_temperature(tt)
    n = 0
    while n < 12:
        print 'Attempting to get to %.3f C' % tt
        time.sleep(10.0)
        t_actual = schr.get_temperature()
        print 'Actual temperature is %.3f, attempting to get to %.3f' % (t_actual, tt)
        if np.abs(tt - t_actual) < 0.03:
            time.sleep(10.0)
            t_actual = schr.get_temperature()
            if np.abs(tt-t_actual) < 0.03:
                print 'Temperature stabilized at %.3f, proceeding.' % t_actual
                break
        n = n + 1

    time.sleep(1.0)
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
            data.add_data_point(pv,tt,curr_freq)
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
