# Piezo offset sweep and measure peak locations

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time
import msvcrt


def simple_line(x, A, B):
    return A*x+B


schr = qt.instruments['schr2']
fp = qt.instruments['fp']
schr.set_current_coupling_gain(0.0)
schr.set_current_coupling(0)
schr.set_piezo_offset(0)


qt.mstart()

data = qt.Data(name='piezo_sweep')

data.add_coordinate('current (mA)')
data.add_value('frequency (GHz)')

data.create_file()


# First, we want to get the current value of the injection current

start_current = schr.get_current()

# The goal is to sweep upward, find a mode hope, then sweep down, find a mode hop
# and then return to the center of these values. The assumption is that this
# point is where the diode is matched with the external cavity "the most"
# and that, then, we can attempt to figure out what the best current feed-forward
# parameters are to maintain a relatively mode hop free scanning capability.

# First, let's sweep upward.

cur_delta = np.linspace(0,35e-3,35)

# Get the current peaks from the FP cavity. These serve as the reference for the
# rest of the experiment.

init_sweep = fp.read_sweep(500,10000,'ai1')

# We can calculate the delta frequency for small deltas from this sweep.
cur_sweep_plot = qt.Plot2D(data,name='csweep_plot')
# Count the peaks
N_peaks = np.size(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))

if N_peaks != 4 and N_peaks != 3:
    print '3 or 4 peaks not found -- probably not single mode to start!'
    # This checks if there are 3 or 4 peaks. If there are, it means the
    # diode is probably single mode (but not guaranteed - you have to look
    # at the actual sweep to tell, but this is a good indicator)

# Let's start sweeping the current upward.
scan_on = True
# Initial frequency is 0
cur_frequency = 0.0
last_frequency = 0.0
prev_sweep = np.copy(init_sweep)
for dc in cur_delta:
    schr.set_current(start_current+dc)
    time.sleep(3)
    if msvcrt.kbhit() or scan_on == False:
        kb_char=msvcrt.getch()
        if kb_char == "q" or scan_on == False:
            scan_on = False
            break
    new_sweep = fp.read_sweep(500,10000,'ai1')
    #fp.read_sweep_plot(500,10000,'ai1')
    df = fp.delta_cross_correlation(np.linspace(0.0,float(500)/10000.,500),prev_sweep,new_sweep)
    prev_sweep = np.copy(new_sweep)
    if np.abs(df) < 0.7:
        # No mode hop.

        cur_frequency = cur_frequency + df
        data.add_data_point(start_current+dc,cur_frequency)
        qt.msleep(0.002)
        print 'No mode hop at %.2f A, current frequency is %.2f' % (start_current+dc, cur_frequency)
        current_limit_high = start_current+dc
        freq_limit_high = cur_frequency
    else:
        # Mode hop detected
        current_limit_high = start_current+dc
        print 'Mode hop detected at %.2f A -- resetting to 0.' % current_limit_high
        freq_limit_high = cur_frequency
        break


schr.set_current(start_current)

# Either we never observed a mode hop or we observed one and have now returned the current to zero offset.

time.sleep(3)

cur_delta = np.linspace(0,-35e-3,35)

cur_frequency = 0.0
last_frequency = 0.0
# Update this sweep; only interested in calculating deltas. We can compare at the end.
prev_sweep = fp.read_sweep(500,10000,'ai1')
if scan_on == True:
    for dc in cur_delta:
        schr.set_current(start_current+dc)
        time.sleep(3)
        if msvcrt.kbhit() or scan_on == False:
            kb_char=msvcrt.getch()
            if kb_char == "q" or scan_on == False:
                scan_on = False
                break
        new_sweep = fp.read_sweep(500,10000,'ai1')
        df = fp.delta_cross_correlation(np.linspace(0.0,float(500.)/10000.,500),prev_sweep,new_sweep)
        prev_sweep = np.copy(new_sweep)
        if np.abs(df) < 0.7:
            # No mode hop.

            cur_frequency = cur_frequency + df
            data.add_data_point(start_current+dc,cur_frequency)
            print 'No mode hop at %.2f A, current frequency is %.2f' % (start_current+dc, cur_frequency)
            current_limit_low = start_current+dc
            freq_limit_low = cur_frequency
        else:
            # Mode hop detected
            current_limit_low = start_current+dc
            print 'Mode hop detected at %.2f A -- resetting to 0.' % current_limit_low
            freq_limit_low = cur_frequency
            break

    schr.set_current(start_current)
    time.sleep(3)

    # Now let's ramp the current to the center value, and monitor the df so we
    # can alert the user how far that changes.
    desired_current = (current_limit_high+current_limit_low)/2.0
    slow_ramp_array = np.linspace(start_current,desired_current,35)

    current_frequency = 0.0
    new_sweep = fp.read_sweep(500,10000,'ai1')
    df = fp.delta_cross_correlation(np.linspace(0.0,float(500)/10000,500),init_sweep,new_sweep)
    print 'Have moved %.3f GHz from the init sweep so far.' % df
    current_frequency = current_frequency + df
    prev_sweep = np.copy(new_sweep)
if scan_on == True:
    for cur_val in slow_ramp_array:
        schr.set_current(cur_val)
        time.sleep(1.0)
        new_sweep = fp.read_sweep(500,10000,'ai1')
        df = fp.delta_cross_correlation(np.linspace(0.0,float(500)/10000.,500),prev_sweep,new_sweep)
        prev_sweep = np.copy(new_sweep)
        current_frequency = current_frequency + df

    print 'Final frequency adjustment is %.3f GHz away from the original frequency.' % current_frequency


data_array = data.get_data()
data.close_file()

# Fit a linear curve to the data


x_guess = np.array((34,0))
fit_output = scp.optimize.curve_fit(simple_line, data_array[:,0]-start_current, data_array[:,1], x_guess, maxfev = 20000)

print 'Fit output is %s' % fit_output
cur_slope = fit_output[0][0]/1000.0
print 'Current slope is %.3f GHz/mA' % cur_slope






