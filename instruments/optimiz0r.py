# Line search optimizer using FSM/XPS
# David Christle <christle@uchicago.edu>, 2014
#
# Taken from original optimiz0r.py file written by Wolfgang P..
#
# This routine executes a line search along X and Y using the FSM, and along Z
# using the XPS. The X and Y are synchronized reads using the DAQ while the XPS
# is a simple set-position/read-counts routine. The three curves are fit with
# Gaussians. The guesses and final fits are plotted, for a diagnostic to check
# the accuracy of the guesses.

from instrument import Instrument
import types
import qt
import msvcrt
import numpy as np
import scipy
from analysis.lib.fitting import fit,common
reload(fit)
reload(common)
class optimiz0r(Instrument):



    def __init__(self, name, dimension_set='default'):
        Instrument.__init__(self, name)
        self.dimension_sets = {
            'default' : {
                'x' : {
                    'scan_length' : 1.5,
                    'nr_of_points' : 40,
                    'qt_ins' : 'fsm',
                    'channel' : 'X',
                    'sigma' : 0.35,
                    'drift_prior' : 0.1,

                    },
                'y' : {
                    'scan_length' : 1.5,
                    'nr_of_points' : 40,
                    'qt_ins' : 'fsm',
                    'channel' : 'Y',
                    'sigma' : 0.35,
                    'drift_prior' : 0.1,
                    },
                'z' : {
                    'scan_length' : 2.0,
                    'nr_of_points' : 50,
                    'qt_ins' : 'xps',
                    'channel' : 'Z',
                    'sigma' : 1.1/1000.0,
                    'drift_prior' : 0.1/1000.0,
                    },
                'xxps' : {
                    'scan_length' : 2.4,
                    'nr_of_points' : 44,
                    'qt_ins' : 'xps',
                    'ins_attr_set' : 'set_abs_positionX',
                    'ins_attr_get' : 'get_abs_positionX'
                    },
                'yxps' : {
                    'scan_length' : 2.4,
                    'nr_of_points' : 44,
                    'qt_ins' : 'xps',
                    'ins_attr_set' : 'set_abs_positionY',
                    'ins_attr_get' : 'get_abs_positionY'
                    },
                'xyz' : ['x','y','z'],
                'xy': ['y','x'],
                'fullxyz' : ['z','y','x','yxps','xxps']
                },
            }
        # Get the instruments we're going to use -- could remove this from
        # hardcode in the future
        self._fsm = qt.instruments['fsm']
        self._xps = qt.instruments['xps']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._opt_pos_prev = {
                                'x' : 0.0,
                                'y' : 0.0,
                                'z' : 3.0
                                }
        self._opt_sigma_prev = {
                                'x' : 0.35,
                                'y' : 0.35,
                                'z' : 1.1/1000.0
                                }
        self._opt_pos = {'x' : self._fsm.get_abs_positionX(),
                         'y' : self._fsm.get_abs_positionY(),
                         'z' : self._xps.get_abs_positionZ()}
        self.add_function('optimize')
        self.dimensions = self.dimension_sets[dimension_set]


    def optimize(self, plot=0, dims='xyz', cycles=1):
        if 'fbl_plot' in qt.plots:
            qt.plots['fbl_plot'].clear()
        qt.plot(name='fbl_plot',clear=True,needtempfile=False)

        for c in range(cycles):
            ret = True

            for d in self.dimensions[dims]:
                scan_length = self.dimensions[d]['scan_length']
                nr_of_points = self.dimensions[d]['nr_of_points']

                if d == 'x':
                    # Get the current position from the FSM - only works if the
                    # FSM has been written to at least once!
                    cur_pos = self._fsm.get_abs_positionX()
                    if cur_pos == None:
                        print 'Current position is none! FSM not initialized?'
                        break

                    self._opt_pos_prev['x'] = cur_pos
                    # Create the desired array of points to measure

                    point_array = np.linspace(cur_pos-scan_length/2.0,cur_pos+scan_length/2.0,nr_of_points)

                    # Add the beginning point so that it's in the array twice
                    # at the beginning. This is because we're looking for the
                    # difference between counter reads to get the counts per
                    # second at a position.
                    temp_point_array = np.insert(point_array,0,cur_pos-scan_length/2.0,axis=0)
                    # Use AO Smooth Goto code to smoothly go to the beginning point
                    self._fsm.AO_smooth(cur_pos, cur_pos-scan_length/2.0, 'X')
                    # Now write the points and get the counts
                    fsm_rate = 40.0 # Hz

                    counts = self._fsm.sweep_and_count(temp_point_array,fsm_rate, 'ctr0','PFI0','X')
                    # Find the difference between readouts to get the counts measured
                    # at a certain point, then divide by the period to get the estimate
                    # of the counts per second
                    cps = np.diff(counts)*fsm_rate
                    # Call the fitting routine to determine optimal position
                    # and set it as self._opt_pos['x'], in this case.
                    ret = self.process_fit(point_array, np.diff(counts), 'x', fsm_rate)
                    print 'Previous x: %.3f, new optimum: %.3f (delta %.1f nm)' % (self._opt_pos_prev['x'], self._opt_pos['x'], self._opt_pos['x']*1E3-self._opt_pos_prev['x']*1E3)
                    #
                    if np.array(point_array).min() < self._opt_pos['x'] < np.array(point_array).max():
                        self._fsm.set_abs_positionX(self._opt_pos['x'])
                    else:
                        self._fsm.set_abs_positionX(self._opt_pos_prev['x'])
                        self._opt_pos['x'] = self._opt_pos_prev['x']
                        print'Optimum outside scan range: Position is set to previous maximum'
                        ret = False

                    #print 'Position changed %d nm' % (self._opt_pos['x']*1E3-self._opt_pos_prev['x']*1E3)

                if d == 'y':
                    # Get the current position from the FSM - only works if the
                    # FSM has been written to at least once!
                    cur_pos = self._fsm.get_abs_positionY()
                    self._opt_pos_prev['y'] = cur_pos
                    # Create the desired array of points to measure
                    point_array = np.linspace(cur_pos-scan_length/2.0,cur_pos+scan_length/2.0,nr_of_points)

                    # Add the beginning point so that it's in the array twice
                    # at the beginning. This is because we're looking for the
                    # difference between counter reads to get the counts per
                    # second at a position.
                    temp_point_array = np.insert(point_array,0,cur_pos-scan_length/2.0,axis=0)
                    ##print 'y temp point array %s' % temp_point_array
                    # Use AO Smooth Goto code to smoothly go to the beginning point
                    self._fsm.AO_smooth(cur_pos, cur_pos-scan_length/2.0, 'Y')
                    # Now write the points and get the counts
                    fsm_rate = 40.0 # Hz
                    counts = self._fsm.sweep_and_count(temp_point_array,fsm_rate, 'ctr0','PFI0','Y')
                    # Find the difference between readouts to get the counts measured
                    # at a certain point, then divide by the period to get the estimate
                    # of the counts per second
                    cps = np.diff(counts)*fsm_rate
                    # Call the fitting routine to determine optimal position
                    # and set it as self._opt_pos['y'], in this case.
                    ret = self.process_fit(point_array, np.diff(counts), 'y', fsm_rate)
                    print 'Previous y: %.3f, new optimum: %.3f (delta %.1f nm)' % (self._opt_pos_prev['y'], self._opt_pos['y'], self._opt_pos['y']*1.0E3-self._opt_pos_prev['y']*1.0E3)

                    if np.array(point_array).min() < self._opt_pos['y'] < np.array(point_array).max():
                        self._fsm.set_abs_positionY(self._opt_pos['y'])
                    else:
                        self._fsm.set_abs_positionY(self._opt_pos_prev['y'])
                        self._opt_pos['y'] = self._opt_pos_prev['y']
                        print 'Optimum outside scan range: Position is set to previous maximum'
                        ret = False

                    #print 'Position changed %d nm' % (self._opt_pos['y']*1.0E3-self._opt_pos_prev['y']*1.0E3)

                if d == 'z':
                    # Get the current position from the FSM - only works if the
                    # FSM has been written to at least once!
                    cur_pos = self._xps.get_abs_positionZ() # In mm, not um!
                    self._opt_pos_prev['z'] = cur_pos # In mm
                    # Create the desired array of points to measure
                    # Note 0.001 converts from mm to um!
                    point_array = np.linspace(cur_pos-scan_length*0.001/2.0,cur_pos+scan_length*0.001/2.0,nr_of_points)

                    # No need to double up on initial position in the point array
                    # since the counts will be read out step-by-step and not in
                    # a synchronized DAQ read.

                    # Set the position on the XPS to the initial point; sort of redundant
                    self._xps.set_abs_positionZ((cur_pos-0.001*scan_length/2.0))
                    # Now write the points and get the counts
                    counts = np.zeros(nr_of_points)
                    # HARDCODED SETUP FOR COUNT READING
                    xps_rate = 40.0 # Hz
                    xps_settle_time = 10.0*0.001 # 10 ms
                    prev_ctr0_src = self._ni63.get_ctr0_src()
                    self._ni63.set_ctr0_src('PFI0')
                    prev_count_time = self._ni63.get_count_time()
                    self._ni63.set_count_time(1.0/xps_rate)
                    for i in range(nr_of_points):
                        # Point array is in mm, so this should work
                        self._xps.set_abs_positionZ((point_array[i]))
                        qt.msleep(xps_settle_time)
                        counts[i] = self._ni63.get_ctr0()
                    self._ni63.set_count_time(prev_count_time)
                    self._ni63.set_ctr0_src(prev_ctr0_src)
                    # Find the difference between readouts to get the counts measured
                    # at a certain point, then divide by the period to get the estimate
                    # of the counts per second
                    cps = counts*xps_rate

                    # Call the fitting routine to determine optimal position
                    # and set it as self._opt_pos['z'], in this case.
                    ret = self.process_fit(point_array, counts, 'z', xps_rate)
                    print 'Previous z: %.6f, new optimum is %.6f (delta %.1f nm)' % (self._opt_pos_prev['z'], self._opt_pos['z'], self._opt_pos['z']*1.0E6-self._opt_pos_prev['z']*1.0E6)


                    #
                    if np.array(point_array).min() < self._opt_pos['z'] < np.array(point_array).max():
                        self._xps.set_abs_positionZ(self._opt_pos['z'])
                    else:
                        self._xps.set_abs_positionZ(self._opt_pos_prev['z'])
                        print 'Optimum outside scan range: Position is set to previous maximum'
                        self._opt_pos['z'] = self._opt_pos_prev['z']
                        ret = False
                    # Use 10^6 instead of 10^3 to convert from mm to nm
                    qt.msleep(0.1)


                qt.msleep(0.05)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break


        return ret

    def process_fit(self, p, cr, dimension, rate):

        # p is position array
        # cr is countrate array
        # d is a string corresponding to the dimension

        # Get the Gaussian width sigma guess from the dimensions dictionary
        sigma_guess = self.dimensions[dimension]['sigma']
        # Always do a Gaussian fit, for now
        self._gaussian_fit = True

        if self._gaussian_fit:
            # Come up with some robust estimates, for low signal/noise conditions
            a_guess = np.array(cr).min()*0.9
            A_guess = np.array(cr).max()-np.array(cr).min()
            p_size = np.size(p)
            x0_guess = p[np.round(p_size/2.0)]  #np.sum(np.array(cr) * np.array(p))/np.sum(np.array(cr)**2)
            #sigma_guess = 1.0#np.sqrt(np.sum(((np.array(cr)-x0_guess)**2)*(np.array(p))) / np.sum((np.array(p))))
            #print 'Guesses: %r %r %r %r' % (a_guess,  x0_guess,A_guess,sigma_guess)
            # The following statement is that if the old sigma is within a factor of 2 of the initial fixed value
            # then use the old sigma as the guess for the fit.
            if np.abs(self._opt_sigma_prev[dimension]/self.dimensions[dimension]['sigma']) < 2.0 and np.abs(self._opt_sigma_prev[dimension]/self.dimensions[dimension]['sigma']) > 0.3:
                sigma_guess = self._opt_sigma_prev[dimension]
            else:
                sigma_guess = self.dimensions[dimension]['sigma']

            ababfunc = (a_guess + A_guess*np.exp(-(p-x0_guess)**2/(2.0*sigma_guess**2)))*rate

            # New method for fitting
            # Purpose of this new method is to improve the robustness of both the numerical optimization and the 'mistracks'
            # that can occur because of nearby luminescence.
            # I think that the numerical robustness can be improved by using a
            # multiple start technique -- start the optimization at a series of
            # guesses across the frequency range, and note the 'best fit' value
            # of each. Pick the one with the best fit of each of these fits, and
            # return it, subject to the condition that the center of the peak is
            # within the min/max of the scan range. This way, it will be no
            # worse than the existing method and significantly more robust to
            # finding the global minimum.
            # The second modification is to reduce the probability of a fit to
            # a nearby peak that isn't the one we want, which I call a 'mistrack'.
            # This can be done by multiplying the likelihood by a prior, making
            # it a Bayesian method, that incorporates our knowledge that the
            # peak we want will not drift too far from scan to scan. Therefore,
            # given two peaks that fit equally well, the one nearest to the old
            # peak location should be selected preferentially. The exception to
            # this is if the other peak fits exceedingly well; these statements
            # are quantified by the exact prior probability distribution used.
            bf_ret = self.bayesian_fit(p, cr, np.array((A_guess, x0_guess, sigma_guess, a_guess)),self.dimensions[dimension]['drift_prior'])

            # check if fit result is in range, process it accordingly
            if not (bf_ret[1] > np.min(p)) and (bf_ret[1] < np.max(p)):
                # fit out of range, set location to max
                pamax = np.argmax(cr)
                self._opt_pos[dimension] = p[pamax]
                ret = False

                print '(%s) fit failed! Set to maximum.' % dimension


            else:
                self._opt_pos[dimension] = bf_ret[1]
                self._opt_sigma_prev[dimension] = bf_ret[2]
                ret = True

                final_fit_func = (bf_ret[3] + bf_ret[0]*np.exp(-(p-bf_ret[1])**2/(2.0*bf_ret[2]**2)))*rate
                qt.plot(p,cr*rate,p,ababfunc,p,final_fit_func,name='fbl_plot',needtempfile=False)
                #print '(%s) optimize succeeded!' % self.get_name()





        return ret
    def bayesian_fit(self, x, y, guess, dp):
        A_guess = guess[0]
        x0_guess = guess[1]
        sigma_guess = guess[2]
        a_guess = guess[3]

        # use the previous guess for the gaussian's center as our prior
        # probability density parameter for the center. The 300 nm width
        # and degrees of freedom are hardcoded here, for now.
        prior = np.array((guess[1], dp, 4))

        # make an anonymous function to substitute into the optimization
        # routine for minimization
        func = lambda p: -1.0*self.bayesian_gauss(p,x,y,prior)

        # now allocate an array to store the 'best fit' values in
        F = np.zeros((10,1))
        xf = np.zeros((10,4))
        ##gval = func(np.array((A_guess, x0_guess, sigma_guess, a_guess)))
        ##print 'guess params are %s' % (np.array((A_guess, x0_guess, sigma_guess, a_guess)))
        ##print 'guess val is %s' % gval
        # set up the array of starting x guesses; other parameters are the same
        # between runs
        x_trials = np.linspace(np.min(x),np.max(x),10)
        for i in range(10):
            res = scipy.optimize.minimize(func, np.array((A_guess, x_trials[i], sigma_guess, a_guess)), jac=None, bounds=((np.min(y)*0.1,np.max(y)), (np.min(x), np.max(x)), (0.05*sigma_guess,5*sigma_guess), (0,0.95*np.max(y))), method='L-BFGS-B', options={'disp': False, 'eps' : 2.0e-8})
            # put the fit results into the array
            if res.success:
                xf[i,:] = res.x
                F[i] = res.fun
            else:
                # fit failed, set F to inf
                xf[i,:] = np.array((A_guess, x_trials[i], sigma_guess, a_guess))
                F[i] = np.inf
        # now determine what fit was the best fit
        min_index = np.argmin(F)
        #print 'x F %s %s' % (xf, F)
        # select out the best fit parameters
        return xf[min_index,:]


    def bayesian_gauss(self, p, x, y, prior):
        # p is a parameter array, x is the spatial position array, y is the
        # array of counts received at that position, prior is an array of
        # fixed prior parameters
        #
        # The parameter array correspondence is:
        # A = p[0], the amplitude of the gaussian
        # mu = p[1], the location parameter of the gaussian
        # sigma = p[2], the width of the gaussian
        # C = p[3], the background offset
        #
        # For the prior parameter density, I use a student t distribution
        # mu_t = prior[0], the location parameter of the student t
        # sigma_t = prior[1], the scale parameter of the student t
        # v_t = prior[2], the degrees of freedom of the student t

        # compute the ideal curve
        y_ideal = p[3] + p[0]*np.exp(-(x-p[1])**2/(2.0*p[2]**2))

        # note the slight tweak to the error weights -- sqrt(y+1.0) to handle
        # the condition when y = 0; won't affect the statistical result much
        chis = np.sum(-(y_ideal - y)**2.0/(2.0*(y+1.0)) - np.log(np.sqrt(y+1.0)))
        # student t argument definition
        q = (p[1]-prior[0])/prior[1]
        # unnormalized log of student t
        priorlogpdf = np.log((prior[2]/(prior[2]+q**2.0))**((1+prior[2])/2.0))

        return (priorlogpdf + chis)
