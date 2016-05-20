# Toptica ECDL laser utilities
# David J Christle <christle@uchicago.edu>, 2016/05/13
#
# The purpose of this module is to aggregate routines used for tuning the performance of my Toptica external cavity
# diode laser. The specific performance factors of interest are the stability of the laser into lasing in a single
# frequency mode, the ability to tune the grating using the internal stepper motor to an arbitrary frequency,
# and the adjustment of the current feed-forward parameter to obtain a large mode hop free tuning range.
#
# These routines are an iterative work-in-progress -- I am not precisely sure what routines will end up giving
# the best performance, so this is simply a collection of things I've tried.
#
#
import numpy as np
import qt
import time
# define the instruments we're using; not bad to define these globally since they are usually defined when qtlab
# is first opened, anyway
topt = qt.instruments['topt']
motdl = qt.instruments['motdl']
fp = qt.instruments['fp']
bristol = qt.instruments['bristol']

def tune_mhf():
    data = qt.Data(name='fast_mhf_tune')
    data.add_coordinate('voltage (V)')
    data.add_value('frequency delta (GHz)')
    data.create_file()
    fastmhfplot = qt.Plot2D(data, name='fast_mhf', coorddim=0, valdim=1)

    topt.set_piezo_voltage(0)
    t0 = time.time()
    # define a few currents to check for single mode behavior
    current_array = np.linspace(0.270,0.306,37)
    stable_currents = []

    # now with the set of stable currents, we want to try to adjust the feedforward parameter and sweep the piezo
    # until we detect a modehop. we need to do it fast, so we'll use only the fabry perot and no wavemeter.
    piezo_array = np.linspace(0,90,120)
    for cur in current_array:
        print 'Setting current to %.4f A.' % (cur)
        topt.set_current(cur)
        topt.set_piezo_voltage(0.0)
        time.sleep(1.0)
        fp_out = fp.check_stabilization()
        if fp_out != 1:
            pass

        # FP is stable, so get the peak locations for frequency tracking
        curr_peaks = np.sort(fp.get_peaks())#np.sort(fp.read_sweep_peaks_lorentzian(500,10000,'ai1'))
        start_sweep = fp.read_sweep(500,10000,'ai1')
        prev_sweep = np.copy(start_sweep)
        curr_freq = 0.0
        if np.size(curr_peaks) < 3:
            print 'Not enough peaks to start at current %.4f' % (cur)
            pass

        print 'Peaks found at current %.4f A, continuing to sweep until a mode hop is detected.' % (cur)

        for pv in piezo_array:
            topt.set_piezo_voltage(pv)

            #new_peaks = np.sort(fp.read_sweep_peaks(500,10000,'ai1'))
            new_peaks = np.sort(fp.get_peaks())
            #new_sweep = fp.read_sweep(500,10000,'ai1')
            if np.size(new_peaks) == 3 or np.size(new_peaks) == 4:
                df = fp.delta_freq_tp(curr_peaks,new_peaks)
                #df = fp.delta_cross_correlation(np.linspace(0.0,float(np.size(500))/10000,np.size(500)),prev_sweep,new_sweep)
                curr_freq = curr_freq + df
                data.add_data_point(pv,curr_freq)
                curr_peaks = new_peaks
                if np.abs(df) > 1.5:
                    print 'Frequency changed by %.2f, mode hop. Total range was %.2f ' % (df, curr_freq-df)
                    stable_currents.append([cur,curr_freq-df]) # the last frequency before a mode hop
                    break
                else:
                    data.add_data_point(pv,curr_freq)
            if np.size(new_peaks) < 3:
                print 'Not enough peaks to start at current %.4f' % (cur)
                #print(new_peaks)
                break

    t1 = time.time()
    print 'MHF tuning took %.1f seconds' % (t1-t0)
    return stable_currents

