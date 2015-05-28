# Piezo offset sweep and measure peak locations

from numpy import pi, random, arange, size
from time import time,sleep
import numpy as np
import time
import msvcrt
import scipy as scp


def simple_line(x, A, B):
    return A*x+B


schr = qt.instruments['schr2']
fp = qt.instruments['fp']
epos = qt.instruments['epos']
schr.set_current_coupling_gain(0.000)
schr.set_current_coupling(0)
schr.set_current_coupling_direction(0)
schr.set_piezo_offset(0)
schr.set_piezo_status(1)


qt.mstart()

data = qt.Data(name='mhf_sweep')

data.add_coordinate('wavelength (nm)')
data.add_coordinate('piezo voltage (V)')
data.add_value('frequency (GHz)')

data.create_file()


# The idea here is to do "shotgun spectroscopy" where we just scan over the entire
# piezo range, step the motor, and repeat the scan. The hope is that we end up
# getting a complete scan from taking the data and stitching it.

# The goal is to sweep upward, find a mode hope, then sweep down, find a mode hop
# and then return to the center of these values. The assumption is that this
# point is where the diode is matched with the external cavity "the most"
# and that, then, we can attempt to figure out what the best current feed-forward
# parameters are to maintain a relatively mode hop free scanning capability.

# First, let's sweep upward.

# Rough conversion factor is 245 GHz/nm
wavelength_start = 1106.480
wavelength_steps_array = np.arange(0,2800,100)

epos.set_wavelength(wavelength_start)
schr.set_piezo_offset(0.0)
print 'Set wavelength to %.3f nm' % wavelength_start
time.sleep(2)
piezo_delta = np.linspace(0,7,35)
schr.set_temperature(22.2)
current_temperature= 22.2
temperature_alternate = [22.2, 22.2-0.5385]


# Get the current peaks from the FP cavity. These serve as the reference for the
# rest of the experiment.

prev_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')

# We can calculate the delta frequency for small deltas from this sweep.
cur_sweep_plot = qt.Plot2D(data,name='tcouple_plot',coorddim=1, valdim=2)
# Count the peaks
N_peaks = np.size(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))

if N_peaks != 4 and N_peaks != 3:
    print '3 or 4 peaks not found -- probably not single mode to start!'
    # This checks if there are 3 or 4 peaks. If there are, it means the
    # diode is probably single mode (but not guaranteed - you have to look
    # at the actual sweep to tell, but this is a good indicator)

# Let's start sweeping the current upward.
scan_on = True
last_wavelength_change_sweep = np.copy(prev_sweep)
last_wavelength_change_frequency = 0.0
init_motor_position = epos.get_motor_position()
current_frequency = 0.0
for ij in range(np.size(wavelength_steps_array)):
    schr.set_piezo_offset(0.0)
    print 'On iteration %d of wavelengths' % ij
    # We need to gently step the wavelength and track the frequency while we do it.
    current_motor_position = epos.get_motor_position()
    motor_deltas = int(np.round(((init_motor_position+wavelength_steps_array[ij])-current_motor_position)/10.0))
    temp_deltas = current_temperature + np.linspace(0,-0.1795*1.5,motor_deltas)
    #schr.set_temperature(temperature_alternate[np.mod(ij,2)])
    #current_temperature = temperature_alternate[np.mod(ij,2)]
    for kj in range(motor_deltas):
        schr.set_temperature(temp_deltas[kj])


        epos.fine_tuning_steps(10)
        time.sleep(2)
        current_temperature = temp_deltas[kj]
        print 'Current temperature is %.3f' % current_temperature

        new_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
        N_peaks = np.size(new_sweep)
        if N_peaks == 3 or N_peaks == 4:
            df = fp.delta_freq_tp(prev_sweep,new_sweep)
            prev_sweep = np.copy(new_sweep)
            current_frequency = current_frequency + df
            if np.abs(df) > 3.4:
                print 'df is %.3f' % df
            qt.msleep(0.005)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break
        #else:
            #print 'Multimode behavior detected -- ignoring point.'
    # Now, we need to bump the motor up a bit if we're still multimode
    attempts = 0
    #Should only run if N_peaks from previous movement are not 3 or 4.
    while N_peaks != 3 and N_peaks != 4 and attempts < 15:
        epos.fine_tuning_steps(10)
        new_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
        N_peaks = np.size(new_sweep)
        if N_peaks == 3 or N_peaks == 4:
            df = fp.delta_freq_tp(prev_sweep,new_sweep)
            prev_sweep = np.copy(new_sweep)
            current_frequency = current_frequency + df
            if np.abs(df) > 3.4:
                print 'df is %.3f' % df
            qt.msleep(0.005)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break
        #else:
        #    print 'Multimode behavior detected -- ignoring point. Attempts beyond wavelength tune are %d' % attempts



    time.sleep(5)
    if msvcrt.kbhit() or scan_on == False:
        kb_char=msvcrt.getch()
        if kb_char == "q" or scan_on == False:
            scan_on = False
            break


    for jj in range(np.size(piezo_delta)):
        schr.set_piezo_offset(piezo_delta[jj])
        time.sleep(0.5)
        N_peaks = np.size(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
        new_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
        if N_peaks == 3 or N_peaks == 4:
            df = fp.delta_freq_tp(prev_sweep,new_sweep)
            prev_sweep = np.copy(new_sweep)
            current_frequency = current_frequency + df
            if np.abs(df) > 3.4:
                print 'df is %.3f' % df
            data.add_data_point(wavelength_steps_array[ij],piezo_delta[jj],current_frequency)
            qt.msleep(0.005)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break
        #else:
        #    print 'Multimode behavior detected -- ignoring point.'
    # Now sweep through the array in reverse.
    for jj in range(np.size(piezo_delta)):
        schr.set_piezo_offset(piezo_delta[-(1+jj)])
        time.sleep(0.5)
        N_peaks = np.size(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
        new_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
        if N_peaks == 3 or N_peaks == 4:
            df = fp.delta_freq_tp(prev_sweep,new_sweep)
            prev_sweep = np.copy(new_sweep)
            current_frequency = current_frequency + df
            if np.abs(df)>3.4:
                print 'df is %.3f' % df
            data.add_data_point(wavelength_steps_array[ij],piezo_delta[-(1+jj)],current_frequency)
            qt.msleep(0.005)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break
        #else:
        #    print 'Multimode behavior detected -- ignoring point.'
    schr.set_piezo_offset(0.0)
    time.sleep(5)
    if msvcrt.kbhit() or scan_on == False:
        kb_char=msvcrt.getch()
        if kb_char == "q" or scan_on == False:
            scan_on = False
            break
    print 'Now scanning in the negative direction'

    for jj in range(np.size(piezo_delta)):
        schr.set_piezo_offset(-1.0*piezo_delta[jj])
        time.sleep(0.5)
        N_peaks = np.size(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
        new_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
        if N_peaks == 3 or N_peaks == 4:
            df = fp.delta_freq_tp(prev_sweep,new_sweep)
            current_frequency = current_frequency + df
            prev_sweep = np.copy(new_sweep)
            if np.abs(df) > 3.4:
                print 'df is %.3f' % df
            data.add_data_point(wavelength_steps_array[ij],-1.0*piezo_delta[jj],current_frequency)
            qt.msleep(0.005)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break
        #else:
        #    print 'Multimode behavior detected -- ignoring point.'
    for jj in range(np.size(piezo_delta)):
        schr.set_piezo_offset(-1.0*piezo_delta[-(1+jj)])
        time.sleep(0.5)
        N_peaks = np.size(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
        new_sweep = fp.read_sweep_peaks_lorentzian(500,10000,'ai1')
        if N_peaks == 3 or N_peaks == 4:
            df = fp.delta_freq_tp(prev_sweep,new_sweep)
            current_frequency = current_frequency + df
            prev_sweep = np.copy(new_sweep)
            if np.abs(df) > 3.4:
                print 'df is %.3f' % df
            data.add_data_point(wavelength_steps_array[ij],-1.0*piezo_delta[-(1+jj)],current_frequency)
            qt.msleep(0.005)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False:
                    scan_on = False
                    break
        #else:
        #    print 'Multimode behavior detected -- ignoring point.'
    current_wavelength_change_sweep = np.copy(prev_sweep)
    # Now let's compare the sweeps to get a delta frequency
    df = fp.delta_freq_tp(last_wavelength_change_sweep,current_wavelength_change_sweep)
    # Calculate the freq. discrepancy as the 'current frequency', which is integrated
    # from step-to-step, with the 'df', which we just calculated from the FP peak lcoations
    # measured at the last wavelength grating change to this one. In theory,
    # these quantities should be equal, but because are constantly adding small errors
    # together, in practice they will differ by some amount.
    freq_discrepancy = ((current_frequency - last_wavelength_change_frequency) - df)
    print 'Frequency discrepancy from accumulated error is %.4f GHz' % freq_discrepancy
    last_wavelength_change_sweep = np.copy(prev_sweep)
    # Let's now update the current frequency to be the more-accurate version
    # based on a single difference versus the many integrated differences
    current_frequency = last_wavelength_change_frequency + df
    last_wavelength_change_frequency = current_frequency



data.close_file()
schr.set_temperature(22.2)
current_temperature= 22.2
