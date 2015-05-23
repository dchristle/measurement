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
                'thresholdV' : 0.1,
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
            idx_width = 0.8/625.0
            distance_array = np.abs(time_axis - (peaks[i] + idx_width))
            upper_idx = np.argmin(distance_array)
            distance_array = np.abs(time_axis - (peaks[i] - idx_width))
            lower_idx = np.argmin(distance_array)
            fit_obj = scp.optimize.curve_fit(self.lorentzian, time_axis[lower_idx:upper_idx], rsamples[lower_idx:upper_idx], np.array((np.max(rsamples[lower_idx:upper_idx])-np.min(rsamples[lower_idx:upper_idx]), 0.08/625.0, peaks[i], np.min(rsamples[lower_idx:upper_idx]))))
            if np.abs(fit_obj[0][2] - peaks[i])*625 > 1.0:
                accurate_peaks.append(peaks[i])
                print 'Lorentzian fit outside range -- returning original threshold peak.'
            else:
                accurate_peaks.append(fit_obj[0][2])


        return np.sort(np.array(accurate_peaks))
    def delta_freq_tp(self, prev, curr):


        return df
    def peak_obj_fun(self,prev,curr):


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













