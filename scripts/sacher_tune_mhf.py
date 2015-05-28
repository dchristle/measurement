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
piezo_vec = np.linspace(-12,12.5,150)
curr_gain = np.linspace(0*10**(-3), 5*10**(-3),12)
schr = qt.instruments['schr2']
fp = qt.instruments['fp']
schr.set_current_coupling_gain(0.0)
schr.set_current_coupling(0)

qt.mstart()

data = qt.Data(name='piezo_sweep')

data.add_coordinate('piezo voltage (V)')
data.add_coordinate('current gain (mA/V)')
data.add_value('frequency (GHz)')

data.create_file()


plot2d = qt.Plot2D(data, name='fp_freq', coorddim=0, valdim=2, maxtraces=20)
schr.set_piezo_offset(piezo_vec[0])

# Now attempt to find a current where the mode is stable
curr_diff = np.linspace(-20*10**(-3),20*10**(-3),55)
curr0 = schr.get_current()
print 'Starting current is %.4f A.' % curr0
for cdf in curr_diff:
    print 'Setting current to %.4f A.' % (curr0+cdf)
    schr.set_current(curr0+cdf)
    time.sleep(3)
    curr_peaks = np.sort(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
    start_sweep = fp.read_sweep(500,10000,'ai1')
    prev_sweep = np.copy(start_sweep)

    if np.size(curr_peaks) < 3:
        print 'Not enough peaks to start at current %.4f' % (curr0+cdf)
    else:
        print 'Peaks found at current %.4f, continuing to sweep until a mode hop is detected.' % (curr0+cdf)
        mhfcl = curr0+cdf
        for cdfn in curr_diff:
            print 'Setting current to %.4f A.' % (mhfcl+cdfn)
            schr.set_current(mhfcl+cdfn)
            time.sleep(3)
            new_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
            new_sweep = fp.read_sweep(500,10000,'ai1')
            if np.size(new_peaks) == 3 or np.size(new_peaks) == 4:
                #df = fp.delta_freq(curr_peaks,new_peaks)
                df = fp.delta_cross_correlation(np.linspace(0.0,float(np.size(500))/10000,np.size(500)),prev_sweep,new_sweep)
                curr_freq = curr_freq + df
                curr_peaks = new_peaks
                if np.abs(df) > 0.7:
                    print 'Frequency changed by %.2f, mode hop.' % df
                    mhfch = mhfcl+cdf
                    break
            if np.size(new_peaks) < 3:
                print 'Not enough peaks to start at current %.4f' % (mhfcl+cdfn)
                mhfch = mhfcl+cdfn
                break
            else:
                print 'Peaks found at current %.4f, continuing to sweep until a mode hop is detected.' % (mhfcl+cdfn)
                mhfch = mhfcl + cdfn
        break

print 'Mode hop free current range is: %.3f to %.3f A' % (mhfcl, mhfch)

avg_mhf = (mhfcl + mhfch)/2.0
print 'Setting current to average %.3f A.' % avg_mhf
schr.set_current(avg_mhf)
schr.set_current_coupling(1)


# preparation is done, now start the measurement.
# It is actually a simple loop.
scan_on = True
for cg in curr_gain:
    schr.set_current_coupling_gain(cg)
    time.sleep(0.1)
    schr.set_piezo_offset(piezo_vec[0])
    time.sleep(5.00)
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
