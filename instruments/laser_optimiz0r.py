# Line search optimizer using Toptica piezo
# David Christle <christle@uchicago.edu>, May 2016
#

from instrument import Instrument
import types
import qt
import msvcrt
import numpy as np
import scipy
from analysis.lib.fitting import fit,common
reload(fit)
reload(common)
class laser_optimiz0r(Instrument):



    def __init__(self, name, dimension_set='default'):
        Instrument.__init__(self, name)
        self.dimension_sets = {
            'default' : {
                'V' : {
                    'scan_length' : 4,
                    'nr_of_points' : 40,
                    'qt_ins' : 'topt',
                    'sigma' : 1.0,
                    'drift_prior' : 1.0,

                    },
                'Vdim' : ['V'],
                },
        }
        # Get the instruments we're going to use -- could remove this from
        # hardcode in the future
        #self._fsm = qt.instruments['fsm']
        #self._xps = qt.instruments['xps']
        self._topt = qt.instruments['topt']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._opt_pos_prev = {
                                'V' : 0.0,

                                }
        self._opt_sigma_prev = {
                                'V' : 1.0,
                                }
        self._opt_pos = {'V' : self._topt.get_piezo_voltage(),
                        }
        self.add_function('optimize')
        self.dimensions = self.dimension_sets[dimension_set]


    def sweep_laser_and_count(self, point_array):
        self._ni63.set_count_time(100.0/1000.0)
        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src('PFI2')
        count_array = np.zeros(point_array.shape)
        for i, voltage in np.ndenumerate(point_array):
            self._topt.set_piezo_voltage(voltage)
            cts = self._ni63.get('ctr1')
            count_array[i] = cts
        return count_array

    def optimize(self, plot=0, dims='Vdim', cycles=1):
        if 'laserscan_plot' in qt.plots:
            qt.plots['laserscan_plot'].clear()
        qt.plot(name='laserscan_plot',clear=True,needtempfile=False)

        for c in range(cycles):
            ret = True
            #print '%s' % self.dimensions['V']
            for d in self.dimensions[dims]:
                #print 'dims is %s d is %s' % (dims, d)
                scan_length = self.dimensions[d]['scan_length']
                nr_of_points = self.dimensions[d]['nr_of_points']

                if d == 'V':
                    # Get the current voltage from the Toptica
                    cur_V = self._topt.get_piezo_voltage()
                    if cur_V == None:
                        logging.warning(__name__ + 'Current position is none! Toptica voltage not initialized?')
                        break

                    self._opt_pos_prev['V'] = cur_V
                    # Create the desired array of points to measure

                    point_array = np.linspace(cur_V-scan_length/2.0,cur_V+scan_length/2.0,nr_of_points)

                    count_array = self.sweep_laser_and_count(point_array)
                    # Find the difference between readouts to get the counts measured
                    # at a certain point, then divide by the period to get the estimate
                    # of the counts per second

                    # Call the fitting routine to determine optimal position
                    # and set it as self._opt_pos['x'], in this case.
                    ret = self.process_fit(point_array, count_array, 'V', 10.0)
                    print 'Previous V/sigma: %.3f/%.3f, new optimum V: %.3f (delta %.1f V)' % (self._opt_pos_prev['V'], self._opt_sigma_prev['V'], self._opt_pos['V'], self._opt_pos['V']-self._opt_pos_prev['V'])
                    #
                    if np.array(point_array).min() < self._opt_pos['V'] < np.array(point_array).max():
                        self._topt.set_piezo_voltage(self._opt_pos['V'])
                    else:
                        self._topt.set_piezo_voltage(self._opt_pos_prev['V'])
                        self._opt_pos['V'] = self._opt_pos_prev['V']
                        print 'Optimum outside scan range: Position is set to previous maximum'
                        ret = False



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
            if np.abs(self._opt_sigma_prev[dimension]/self.dimensions[dimension]['sigma']) < 3.0 and np.abs(self._opt_sigma_prev[dimension]/self.dimensions[dimension]['sigma']) > 0.3:
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
                qt.plot(p,cr*rate,p,ababfunc,p,final_fit_func,name='laserscan_plot',needtempfile=False)
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
