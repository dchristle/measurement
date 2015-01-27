# Calibrate ECDL current coupling

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time
import msvcrt






schr.set_piezo_offset(0.0)
time.sleep(1.0)
curr0 = schr.get_current()

piezo_diff = np.linspace(-0.35,0.35,20)
curr_diff = np.linspace(-10*10**(-3),10*10**(-3),20)
schr = qt.instruments['schr']
fp = qt.instruments['fp']





# preparation is done, now start the measurement.
# It is actually a simple loop.
scan_on = True



time.sleep(1.00)
curr_freq = 0

time.sleep(1.0)
curr_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
if np.size(curr_peaks) < 3:
    print 'Not enough peaks to start!'
    scan_on = False
else:
    print 'Peaks found, starting.'
    scan_on = True
    qt.mstart()

    data = qt.Data(name='piezo_sweep')
    data.add_coordinate('piezo voltage (V)')
    data.add_value('frequency (GHz)')
    data.create_file()
    plot2d = qt.Plot2D(data, name='fp_freqvsp', coorddim=0, valdim=1, maxtraces=5)
    for pv in piezo_diff:
        schr.set_piezo_offset(pv)
        time.sleep(0.5)
        new_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
        time.sleep(0.1)
        #fp.read_sweep_plot(500, 10000, 'ai1')

        if np.size(new_peaks) == 3 or np.size(new_peaks) == 4:
            df = fp.delta_freq(curr_peaks,new_peaks)
            curr_freq = curr_freq + df
            curr_peaks = new_peaks
            data.add_data_point(pv,curr_freq)
            print 'peaks: %s, df %.3f' % (new_peaks, df)
        else:
            print 'peaks: %s' % new_peaks




        if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q":
                        scan_on = False
                        break
    data.close_file()



print 'Starting current sweep.'
schr.set_piezo_offset(0.0)

time.sleep(3.0)
curr_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
if np.size(curr_peaks) < 3:
    print 'Not enough peaks to start!'
    scan_on = False
else:
    print 'Peaks found at piezo starting.'
    scan_on = True
    qt.mstart()

    data = qt.Data(name='current_sweep')
    data.add_coordinate('current (mA)')
    data.add_value('frequency (GHz)')
    data.create_file()
    plot2d = qt.Plot2D(data, name='fp_freqvscurr', coorddim=0, valdim=1, maxtraces=5)
    for d_curr in curr_diff:
        schr.set_current(curr0 + d_curr)
        time.sleep(2.5)
        new_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
        time.sleep(0.1)
        #fp.read_sweep_plot(500, 10000, 'ai1')

        if np.size(new_peaks) == 3 or np.size(new_peaks) == 4:
            df = fp.delta_freq(curr_peaks,new_peaks)
            curr_freq = curr_freq + df
            curr_peaks = new_peaks
            data.add_data_point(curr0+d_curr,curr_freq)
            print 'peaks: %s, df %.3f' % (new_peaks, df)
        else:
            print 'peaks: %s' % new_peaks




        if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q":
                        scan_on = False
                        break
    data.close_file()
schr.set_current(curr0)
qt.mend()
