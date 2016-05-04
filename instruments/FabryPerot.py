# Fabry-Perot interferometer
#
# David Christle <christle@uchicago.edu>, April 27, 2014
# The purpose of this code is to calibrate and read out the Fabry Perot interferometer
# using the DAQ. The original settings are a 50 ms period, using a 20 V sweep.
# The DAQ should trigger off a digital tick that occurs on PFI 8, and read at a
# fixed rate for 50 ms. This code is pretty messy considering I'm trying to push
# this paper out pretty quickly.

from instrument import Instrument
import types
import logging
import time
import qt
import math
import numpy as np
import scipy.cluster
import msvcrt
import lib.math.peakfind as pf
import scipy as scp
import pandas as pd



class FabryPerot(Instrument):

    def __init__(self, name, channels=2):
        Instrument.__init__(self, name, tags=['positioner'])
        # Import the DAQ as an instrument object to write to.
        self._ni63 = qt.instruments['NIDAQ6363']

        # Store related constants for the interferometer here.

        self.FP_params = {
                'period' : 50.0, # ms
                'FSR' : 10.0, # GHz
                'sample_rate' : 10000.0, # Hz
                'aiport' : 'ai1',
                'trigger' : 'PFI8',
                'thresholdV' : 0.25,
                'thresholdsep' : 0.006

                }
        # Instrument parameters
##        self.add_parameter('abs_position',
##            type=types.FloatType,
##            channels=('X', 'Y'),
##            flags=Instrument.FLAG_SET|Instrument.FLAG_SOFTGET,
##            units='um',
##            format='%.04f')


        # Instrument functions
        self.add_function('read_sweep')
        self.add_function('read_sweep_plot')
        self.add_function('read_sweep_centroids')

        self.add_function('threshold')







    def read_sweep_plot(self, samples, rate, channel):
        rsamples = self.read_sweep(samples, rate, channel)
        time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
        qt.plot(time_axis,rsamples,name='fpplot1',traceofs=0.2, maxtraces=20)
        return
    def focus_plot(self):
        samples = 500
        rate = 10000
        channel = 'ai1'
        while True:
            rsamples = self.read_sweep(samples, rate, channel)
            time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
            qt.plot(time_axis,rsamples,name='fp_focusplot',traceofs=0.2, clear=True,linewidth=6)
            time.sleep(0.5)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break

        return
    def max_intensity_plot(self):
        samples = 500
        rate = 10000
        channel = 'ai1'
        data = qt.Data(name='max_int_data')
        data.add_coordinate('time')
        data.add_value('intensity')
        t0 = time.time()
        miplot = qt.plot(data,name='maxintplot',traceofs=0.2, coorddim=0, valdim=1, maxpoints=100, linewidth=6, clear=True)
        while True:
            rsamples = self.read_sweep(samples, rate, channel)
            data.add_data_point((time.time()-t0),np.max(rsamples))
            miplot.update()

            time.sleep(0.5)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break

        return
    def avg_max_intensity_plot(self):
        samples = 500
        rate = 10000
        channel = 'ai1'
        data = qt.Data(name='avg_max_int_data')
        data.add_coordinate('time')
        data.add_value('intensity')
        t0 = time.time()
        miplot = qt.plot(data,name='avgmaxintplot',traceofs=0.2, coorddim=0, valdim=1, maxpoints=100,linewidth=7,clear=True)
        while True:
            rsamples = self.read_sweep(samples, rate, channel)
            time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
            p = pf.PeakFinder(time_axis,rsamples, maxpeaks = 4,threshold = 4,fitwidth=35,fitfunc='FIT_LORENTZIAN',amp_threshold=0.2)
            peaks = p.find(sign=1, bgorder=1)
            peak_array = np.array(peaks)
            data.add_data_point((time.time()-t0),np.mean(peak_array[:,1]))
            miplot.update()
            time.sleep(0.5)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break

        return
    def read_sweep_centroids(self, samples, rate, channel):

        rsamples = self.read_sweep(samples, rate, channel)
        time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
        tsamples = self.threshold(time_axis, rsamples)

        return tsamples
    def read_sweep_peaks(self, samples, rate, channel):

        rsamples = self.read_sweep(samples, rate, channel)
        time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
        tsamples = self.threshold(time_axis, rsamples)
        peaks = self.find_peaks(tsamples)

        return np.sort(peaks)
    def lorentzian(self, x, A, gamma, x0, C):
        return A*(gamma*gamma/(gamma*gamma+np.power(x-x0,2)))+C

    def read_sweep_peaks_lorentzian(self, samples, rate, channel):

        rsamples = self.read_sweep(samples, rate, channel)
        time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
        tsamples = self.threshold(time_axis, rsamples)
        peaks = self.find_peaks(tsamples)
        accurate_peaks = []
        for i in range(np.size(peaks)):
            # Find the points around the peak
            distance_array = np.abs(time_axis - peaks[i])
            min_idx = np.argmin(distance_array)
            idx_width = 0.9/625.0
            distance_array = np.abs(time_axis - (peaks[i] + idx_width))
            upper_idx = np.argmin(distance_array)
            distance_array = np.abs(time_axis - (peaks[i] - idx_width))
            lower_idx = np.argmin(distance_array)
            fit_obj = scp.optimize.curve_fit(self.lorentzian, time_axis[lower_idx:upper_idx], rsamples[lower_idx:upper_idx], np.array((np.max(rsamples[lower_idx:upper_idx])-np.min(rsamples[lower_idx:upper_idx]), 0.00025, peaks[i], 0)),maxfev=20000)
            if np.abs(fit_obj[0][2] - peaks[i])*625 > 1.0:
                accurate_peaks.append(peaks[i])
                print 'Lorentzian fit outside range -- returning original threshold peak. %.4f' % (fit_obj[0][2])
            else:
                accurate_peaks.append(fit_obj[0][2])


        return np.sort(np.array(accurate_peaks))
    def delta_freq_tp(self, prev, curr):
        Npoints = 500
        df_array = np.linspace(-4.5/625.0,4.5/625,Npoints)
        mse = np.zeros(Npoints)
        for i in range(Npoints):
            mse[i] = self.peak_obj_fun(prev,curr,df_array[i])
        min_idx = np.argmin(mse)
        pof_anon = lambda x: self.peak_obj_fun(prev,curr,x)
        fit_width = np.abs(df_array[1] - df_array[0])
        #min_obj = scp.optimize.minimize(pof_anon,df_array[min_idx],method='L-BFGS-B', bounds=((df_array[min_idx]-2*fit_width,df_array[min_idx]+2*fit_width)))
        min_obj = scp.optimize.minimize_scalar(pof_anon,df_array[min_idx],method='bounded', bounds=(df_array[min_idx]-2*fit_width,df_array[min_idx]+2*fit_width))
        return float(min_obj.x*625.0)
    def peak_obj_fun(self,prev,curr,df):
        # Find the two smallest peak displacements
        Nprev = np.size(prev)
        Ncurr = np.size(curr)
        displacements = np.zeros(Nprev*Ncurr)
        for i in range(Nprev):
            for j in range(Ncurr):
                new_freq = prev[i]+df
                # This "if" statement is adjusting for the edges of the
                # active-scanned FP readout -- it adjusts for the wrapover
                # of the peaks.
                if new_freq > 0.0495:
                    new_freq = new_freq - 0.0495 + 0.001
                if new_freq < 0.001:
                    new_freq = new_freq + 0.0495 - 0.001
                displacements[i*Ncurr+j] = np.power(new_freq-curr[j],2.0)
        # Now sort the array
        sorted_displacements = np.sort(displacements)
        # And return the two or three shortest displacement distances, squared and summed.
        return sorted_displacements[0] + sorted_displacements[1] + sorted_displacements[2]



        return Q
    def read_sweep_peakfinder(self, samples, rate, channel):

        rsamples = self.read_sweep(samples, rate, channel)
        time_axis = np.linspace(0.0,float(np.size(rsamples))/rate,np.size(rsamples))
        p = pf.PeakFinder(time_axis,rsamples, maxpeaks = 4,threshold = 4,fitwidth=35,fitfunc='FIT_LORENTZIAN',amp_threshold=0.2)
        peaks = p.find(sign=1, bgorder=1)

        return np.sort(np.array(peaks)[:,0])

    def read_sweep(self, samples, rate, channel):
        devchan = channel
        rsamples = self._ni63.readarray(samples, self.FP_params['trigger'], rate, -10.0, 10.0, 10.0, devchan)
        return rsamples

    def threshold(self, time, samples):
        tsamples = time[samples > self.FP_params['thresholdV']]
        #print tsamples
        #print samples[samples > self.FP_params['thresholdV']]
        return tsamples

    def kmeans_centroid(self, samples):

        for k in range(4):
            centroids, labels = scipy.cluster.vq.kmeans2(samples, 5-k, 500)

            sorted_centroids = np.sort(centroids)
            #print sorted_centroids
            diffs = np.diff(sorted_centroids)
            min_separation = np.min(diffs)
            if min_separation > self.FP_params['thresholdsep']:
                break
        print 'Found %s peaks, centroids %s' % (5-k, sorted_centroids)
        return sorted_centroids


    def find_peaks(self, samples):
        sorted_samples = np.sort(samples)
        sorted_samples_diff = np.diff(sorted_samples)
        Npeaks = 1 + np.size(sorted_samples_diff[sorted_samples_diff > self.FP_params['thresholdsep']])

        ac = np.arange(np.size(sorted_samples_diff))

        split_arrays = np.split(sorted_samples,ac[sorted_samples_diff > self.FP_params['thresholdsep']]+1)

        peak_means = np.zeros(Npeaks)
        for i in range(Npeaks):
            peak_means[i] = np.mean(split_arrays[i])

        return peak_means

    def delta_freq(self, prev, curr):
        # prev is an array of previous time values
        # current is an array of the current time values
        # Goal of this function is to compute the change in frequency based on
        # the known free spectral range of the interferometer.
        # Idea is to take only the peaks that "survive" between the difference
        # and then find the minimum frequency difference that produces agreement.

        if np.size(curr) == np.size(prev):
            raw_diff = curr-prev
            if raw_diff.any() == False:
                return 0.0
        deltaF = np.linspace(-4,4,900)
        mses = np.zeros(900)
        for ij in range(900):
            mses[ij] = self.peak_mse(prev,curr,deltaF[ij])

        # Find all mse's that are local minima
        lm = (np.diff(np.sign(np.diff(mses))) > 0).nonzero()[0] + 1
        # Get the absolute magnitudes of lm's
        dfm = np.abs(deltaF[lm])
        lmas = lm[dfm.argsort()]
        simple_est = (curr[0]-prev[0])/0.00156 # set by 10 ghz and time separation
        print 'simple estimate %.3f' % simple_est

        return deltaF[lmas[0]]
    def peak_mse(self, prev, curr, deltaF):
        N_prev = np.size(prev)
        N_curr = np.size(curr)
        #print 'N_prev %d, Ncurr %d' % (N_prev, N_curr)
        M = np.min(np.array([N_prev,N_curr]))
        peakmses = np.zeros(N_prev*N_curr)

        for i in range(N_prev):
            for j in range(N_curr):

                if j == (N_curr-1):
                    slope = (curr[j] - curr[j-1])/10.0
                elif j == 0:
                    slope = (curr[j+1]-curr[j])/10.0
                else:
                    slope = ((curr[j+1]-curr[j])/10.0 + (curr[j] - curr[j-1])/10.0)/2.0
                peakmses[i*N_curr + j] = np.power(curr[j] - np.mod((prev[i] + deltaF*slope),0.5),2)
        sortedpmses = np.sort(peakmses)
        # Compute the sum of the M smallest squared errors
        mse = 0
        for k in range(M):
            mse = mse + sortedpmses[k]


        # Return this sum of squared errors
        return mse

    def check_stabilization(self):
        raw_sweep = self.read_sweep(490,10000,'ai1')
        t_axis = np.linspace(1,490,490)

        # use a threshold to filter out only points where the photodiode records signal
        v_pd_thresh = 0.005
        filt_idx = raw_sweep > v_pd_thresh

        # create a pandas dataframe
        sweep_frame = pd.DataFrame({'time':t_axis[filt_idx], 'voltage':raw_sweep[filt_idx]})
        # now return the classification
        return self.classify_fabry_perot(sweep_frame)

    def classify_fabry_perot(self, dat):
        # classify the no signal condition

        if np.max(dat['voltage'])-np.min(dat['voltage']) < 0.02:
            # no signal case is triggered, return 3
            return 3

        threshold_pct = 0.06

        peakidxs = dat['voltage'].values > (threshold_pct*(np.max(dat['voltage'])-np.min(dat['voltage'])) + np.min(dat['voltage']))
        if np.size(peakidxs) == 0:
            peak_counter = 0
            return 3
        else:
            peak_counter = 1
            times = dat['time'].values
            times = times[peakidxs]
            #print(times)
            for i in range(np.size(times)-1):
                # count the number of gaps

                if times[i+1] - times[i] > 10:
                    peak_counter = peak_counter + 1

        if peak_counter > 6:
            #print('peaks are {:d}'.format(peak_counter))
            return 2
        else:
            #print('peaks are {:d}'.format(peak_counter))
            return 1












