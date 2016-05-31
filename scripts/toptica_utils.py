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
import math
import progressbar
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


def set_toptica_to_analog_scan():
    topt = qt.instruments['topt']
    topt.set_scan_enabled(False)
    # set the piezo and FF to zero
    topt.get_piezo_voltage()
    topt.set_piezo_voltage(0.1)
    # set the scan output to external channel A (channel A is channel 20, the piezo direct is channel 50)
    topt.set_scan_output(20)
    topt.set_external_input(True)
    # set the analog scaling factor to 17.5 (c.f. with Andrew Grimes, the engineer at Toptica)
    topt.set_external_input_factor(17.5)
    # set the external input signal to "Fine Input 1", which is channel 0
    topt.set_external_input_signal(0)
    # now set the external input to true

    return


def set_toptica_to_regular_scan():
    topt = qt.instruments['topt']
    topt.set_scan_enabled(False)
    # set the piezo to zero
    topt.get_piezo_voltage()
    topt.set_piezo_voltage(0.0)
    # turn the external input off
    topt.set_external_input(False)
    # set analog scaling factor to 0 - probably doesn't do anything
    topt.set_external_input_factor(0)
    # set the scan output to regular piezo channel (channel A is channel 20, the piezo direct is channel 50)
    topt.set_scan_output(50)
    # again ensure the piezo voltage is 0
    topt.get_piezo_voltage()
    topt.set_piezo_voltage(0.0)
    return

def prepare_external_scan():
    topt = qt.instruments['topt']
    # set parameters like the offset/amplitude/frequency
    topt.set_scan_offset(45.0/17.5) # set to the equivalent of 45 V on the piezo
    topt.set_scan_amplitude(60.0/17.5) # set the amplitude to 60, so that we get +- 30 V (60 Vpp) of scan
    topt.set_scan_frequency(5.0)
    return

def read_fp_output():
    topt = qt.instruments['topt']
    scan_frequency = topt.get_scan_frequency()
    scan_rate = 10000.0
    scan_samples = math.floor(1.0/scan_frequency*scan_rate)
    raw_sweep = self.read_sweep(scan_samples,scan_rate,'ai1')
    t_axis = np.linspace(1,490,490)

    # use a threshold to filter out only points where the photodiode records signal
    v_pd_thresh = 0.005
    filt_idx = raw_sweep > v_pd_thresh

    # create a pandas dataframe
    sweep_frame = pd.DataFrame({'time':t_axis[filt_idx], 'voltage':raw_sweep[filt_idx]})
    # now return the classification
    return self.count_fabry_perot(sweep_frame)


def read_sweep(samples, trigchan, rate, threshold_voltage, channel):
    ni63 = qt.instruments['NIDAQ6363']
    devchan = channel
    rsamples = ni63.read_array_analog(samples, trigchan, rate, threshold_voltage, -10, 10, 5, devchan)
    return rsamples


def get_peaks(data):

    if np.max(data[:,1]) - np.min(data[:,1]) < 0.02:
        # no signal case is triggered, return 3
        return -1, -1

    threshold_pct = 0.01

    peakidxs = data[:,1] > (threshold_pct*(np.max(data[:,1])-np.min(data[:,1])) + np.min(data[:,1]))

    if np.size(peakidxs) == 0:
        return -1, -1
    else:
        peak_locations = []
        current_peak_vals = []
        times = data[:,0]
        times = times[peakidxs]
        peak_counter = 1
        for i in range(np.size(times)-1):
            current_peak_vals.append(times[i])
            # count the number of gaps and keep track of points within this peak
            if times[i+1] - times[i] > 0.0001:  # 0.1 milliseconds
                peak_counter = peak_counter + 1
                peak_locations.append(float(np.mean(np.array(current_peak_vals))))
                current_peak_vals = []

    peak_locations.append(float(np.mean(np.array(current_peak_vals))))

    return peak_counter, peak_locations


def get_peaks_and_count(threshold_voltage = 0.1):
    raw_data = read_sweep(40000,'APFI0',200000,threshold_voltage,'ai4')
    cdata = np.column_stack((np.linspace(0,0.2,40000), raw_data))
    n_peaks, peak_locs = get_peaks(cdata)
    return n_peaks, peak_locs

def grid_search_print():
    qt.mstart()
    data = qt.Data(name='grid_search_MHF')
    data.add_coordinate('current (mA)')
    data.add_coordinate('feedforward (mA/V)')
    data.add_value('Mode Hops')
    data.create_file()
    plot3d = qt.Plot3D(data, name='measure_MH3D', style='image')
    topt = qt.instruments['topt']
    set_toptica_to_analog_scan()
    #topt.set_scan_offset(2.0)
    #topt.set_scan_amplitude(4.0)
    topt.set_scan_offset(45.0/17.5) # set to the equivalent of 45 V on the piezo
    topt.set_scan_amplitude(60.0/17.5) # set the amplitude to 60, so that we get +- 30 V (60 Vpp) of scan
    topt.set_scan_frequency(5.0)
    topt.set_scan_enabled(True)
    current_values = np.linspace(0.230,0.300,30)
    feedforward_values = np.linspace(-0.48,-0.60,5)
    laser_performance = np.zeros((np.size(current_values)*np.size(feedforward_values),1))
    for i in range(np.size(current_values)):
        for j in range(np.size(feedforward_values)):
            while k < 10:
                cur = current_values[i]
                ff = feedforward_values[j]
                topt.set_feedforward_factor(ff)
                topt.set_current(cur)
                time.sleep(0.2)
                threshold_voltage = 45.0/17.5 - 0.5*60.0/17.5+0.1
                n_peaks, peak_locs = get_peaks_and_count(threshold_voltage)
                peakdiffs = np.diff(np.array(peak_locs))
                norm_peakdiffs = peakdiffs/np.mean(peakdiffs)
                n_modehops = np.sum(np.abs(norm_peakdiffs) > 0.2)
                #print 'i %s j %s' % (i ,j)
                laser_performance[np.size(feedforward_values)*i+j] = n_modehops

                #print 'Found %d modehops for current %.3f A and FF %.3f' % (n_modehops, cur, ff)
                qt.msleep(0.001)
                if n_modehops > 1:
                    k = 10
                    break
            data.add_data_point(cur,ff,n_modehops)

        data.new_block()
    qt.mend()
    return laser_performance

def grid_search_optimize():
    topt = qt.instruments['topt']
    set_toptica_to_analog_scan()
    topt.set_scan_offset(2.0)
    topt.set_scan_amplitude(4.0)
    topt.set_scan_frequency(5.0)
    topt.set_scan_offset(25.0/17.5) # set to the equivalent of 45 V on the piezo
    topt.set_scan_amplitude(50.0/17.5) # set the amplitude to 60, so that we get +- 30 V (60 Vpp) of scan
    topt.set_scan_frequency(5.0)
    topt.set_scan_enabled(True)
    current_values = np.linspace(0.25,0.283,15)
    feedforward_values = np.linspace(-0.40,-0.51,5)
    laser_performance = np.zeros((np.size(current_values)*np.size(feedforward_values),3))
    masize = np.size(current_values)*np.size(feedforward_values)
    with progressbar.ProgressBar(max_value=masize) as bar:
        for i in range(np.size(current_values)):
            for j in range(np.size(feedforward_values)):
                k = 0
                while k < 10:
                    cur = current_values[i]
                    ff = feedforward_values[j]
                    topt.set_feedforward_factor(ff)
                    topt.set_current(cur)
                    threshold_voltage = 25.0/17.5 - 0.5*50.0/17.5+0.1
                    n_peaks, peak_locs = get_peaks_and_count(threshold_voltage)
                    if np.size(peak_locs) < 2:
                        n_modehops = -1
                        laser_performance[np.size(feedforward_values)*i+j,0] = cur
                        laser_performance[np.size(feedforward_values)*i+j,1] = ff
                        laser_performance[np.size(feedforward_values)*i+j,2] = -1
                    else:
                        peakdiffs = np.diff(np.array(peak_locs))
                        #print('peak locs %s' % peak_locs)
                        norm_peakdiffs = peakdiffs/np.mean(peakdiffs)
                        n_modehops = np.sum(np.abs(norm_peakdiffs) > 0.2)
                        bar.update(np.size(feedforward_values)*i+j)
                        laser_performance[np.size(feedforward_values)*i+j,0] = cur
                        laser_performance[np.size(feedforward_values)*i+j,1] = ff
                        laser_performance[np.size(feedforward_values)*i+j,2] = n_modehops
                        #data.add_data_point(cur,ff,n_modehops)
                        #print 'Found %d modehops for current %.3f A and FF %.3f' % (n_modehops, cur, ff)
                        qt.msleep(0.001)
                    if n_modehops > 0:
                        k = 10
                        break
                    else:
                        print 'Found < 1 modehops -- trying again...'
                        k = k + 1

            #data.new_block()
        #qt.mend()

    # take the points in the array with MHs greater than 20
    lp_thresh_idx = laser_performance[:,2] > 20
    laser_performance_threshed = laser_performance[lp_thresh_idx,:]
    idx = np.argmin(laser_performance_threshed[:,2])
    topt.set_scan_enabled(False)
    set_toptica_to_regular_scan()
    topt.set_scan_enabled(False)
    #best_array = laser_performance_threshed[idx,:]
    topt.set_feedforward_factor(laser_performance_threshed[idx,1])
    topt.set_current(laser_performance_threshed[idx,0])
    return

def brent_search(f, a, b, tol, maxiters):
    # This method is a fairly general implementation of Brent's search. The
    # usage goal for it is to take advantage of the relatively smooth nature
    # of the relation between grating angle and laser frequency, and allow
    # the program to hunt for a desired laser frequency.
    fa = f(a)
    fb = f(b)
    if fa*fb >= 0:
        print('Root is not bracketed -- exiting root search. a(%d) = %.2f, b(%d) = %.2f.' % (a, fa, b, fb))
        # find the closest non-root, assuming the function is smooth
        r_idx = np.argmin(np.array((fa,fb)))
        r = np.array((a,b))
        return r[r_idx], -1
    if (np.abs(a) < np.abs(b)):
        # swap a and b
        a, b = b, a
        fa, fb = fb, fa

    c = a
    fc = fa
    fs = -1.0
    mflag = True
    i = 0

    while (fb != 0.0 and fs != 0.0 and np.abs(b - a) > np.abs(tol)) and i < maxiters:
        if fa != fc and fb != fc:
            # do inverse quadratic interpolation
            s = a*fb*fc/((fa - fb)*(fa - fc)) + b*fa*fc/((fb - fa)*(fb - fc)) + c*fa*fb/((fc - fa)*(fc - fb))
        else:
            # do secant method
            s = b - fb*(b - a)/(fb - fa)
        # check conditions to see if we do the bisection method instead
        tmp = (3*a+b)/4.0
        if (not ((tmp < s and s < b) or (b < s and s < tmp)) or \
            (mflag and np.abs(s-b) >= np.abs(b-c)/2.0) or \
            ((not mflag) and np.abs(s-b) >= np.abs(c-d)/2.0)):
            # do bisection instead
            s = (a+b)/2.0
            mflag = True

        else:
            if ((mflag and np.abs(b-c) < np.abs(tol)) or \
                 ((not mflag) and np.abs(c-d) < np.abs(tol))):
                # do bisection instead
                s = (a+b)/2.0
                mflag = True
            else:
                mflag = False

        fs = f(s)
        d = c
        c = b
        fc = fb
        if fa*fs < 0.0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs

        if np.abs(fa) < np.abs(fb):
            a, b = b, a
            fa, fb = fb, fa
        i = i+1

    return b, i


# def find_single_mode():
#     fp = qt.instruments['fp']
#     fp_out = fp.check_stabilization()
#     if fp_out == 2:
#         # try to wrap around the current list, so we aren't always starting at the bottom
#         current_sweep_list = np.arange(self.params['current_range'][0],self.params['current_range'][1],0.0005)
#         N_currents = np.size(current_sweep_list)
#         nearest_idx = self.find_nearest(current_sweep_list,self._toptica.get_current())
#         # now that we found the nearest, sweep through the list, starting at this index
#         for idx_offset in range(N_currents):
#             self._toptica.set_current(current_sweep_list[(nearest_idx + idx_offset) % N_currents])
#             time.sleep(0.25)
#             fp_out = self._fp.check_stabilization()
#             if fp_out == 1:
#                 print 'FSM: Motor pos: %d, FP %d, GHz: %.2f, cur %.3f' % (self._motdl.get_position(), self._fp.check_stabilization(), self._wvm.get_frequency(), self._toptica.get_current())
#                 break
#     return
#
# def set_motor_and_measure(self, motor_pos):
#     motdl.high_precision_move(motor_pos)
#     time.sleep(0.25)
#     self.find_single_mode()
#     time.sleep(0.2)
#     #self.check_laser_stabilization()
#     wl = bristol.get_wavelength()
#
#
#     if wl != 0.0:
#         freq = 299792458.0/wl
#     else:
#         logging.error(__name__ + ': Bristol wavelength returned 0 -- check alignment?')
#         freq = 0.0
#     print 'Brent search frequency: %.2f GHz...' % freq
#     return freq
# def laser_frequency_seek_rf(frequency):
#     # determine the low and high motor positions
#     bracket_init = 120.0
#     f_low = frequency - bracket_init
#     f_high = frequency + bracket_init
#     a = 2.424e-7
#     b = -297.1101
#     c = 264841.49
#     motor_array_mean = 155264.44
#     motor_array_std = 3086.62
#     # linear
#     m_low, its = brent_search(lambda x: b*(x - motor_array_mean)/motor_array_std + c - f_low, 120000, 190000, 0.1, 60)
#     m_high, its = brent_search(lambda x: b*(x - motor_array_mean)/motor_array_std + c - f_high, 120000, 190000, 0.1, 60)
#     if its == -1:
#         print logging.error(__name__ + ': Calibration constants are off by too much to find proper motor boundaries.')
#         return False
#     # check if a reference of the motor is needed. this can happen because the motor occasionally stops responding
#     #self.check_if_reference_is_needed()
#
#     # create an array of motor points to sweep - maybe 10 -- with nice round numbers for the motor positions
#     # first round to the nearest 100
#     mceil_high = np.ceil(m_high/100.0)*100
#     mceil_low = np.floor(m_low/100.0)*100
#
#
#     m_target, its = self.brent_search((lambda mot_pos: self.set_motor_and_measure(int(np.round(mot_pos))) - frequency), m_low, m_high, 3, 16)