# generalized pattern search optimizer
# David Christle <christle@uchicago.edu>, August 2014
#
# This virtual instrument collects a dictionary of control functions and an
# output function, and then uses a derivative-free method called a pattern
# search to minimize this function. It also has a sample size control method
# that is based on a hypothesis test, which allows it to converge toward an
# optimum when the change in function value is smaller than the current
import qt
from instrument import Instrument
import numpy as np
import gobject
import types
import time
import logging
import measurement.lib.qsturng as qst
import msvcrt

class dmads_opt(Instrument):
    def __init__(self, name, **kw):
        Instrument.__init__(self, name)





        self.add_function('add_control_variable')
        self.add_function('remove_control_variable')
        self.add_function('get_control_variables')
        self.add_function('set_output_variable')
        self.add_function('iterate')
        self.add_function('optimize')
        self.add_function('get_current_x')
        self.add_function('set_x0')
        self.add_function('set_max_sample_size')
        self.add_function('get_max_sample_size')
        self.add_function('set_debug_on')
        self.add_function('set_debug_off')
        self.add_function('set_iterations')
        self.add_function('get_iterations')

        # Constants for the pattern search
        self._sleeptime = 0.1

        self._debug = False
        self._cur_f = 0.0
        self._nfeval = 0
        self._cur_sse = 0.0
        self._mesh_size = 1.0
        self._sample_size = 2
        self._max_sample_size = 4
        self._history = []
        self._controlv_f = []
        self._controlv_names = []
        self._x0 = []
        self._cur_x = []
        self._lb = []
        self._ub = []

        self._optimization_is_on = False




    #--------------get_set
    def reset(self):
        # resets the instrument completely, getting rid of all variables
        self._cur_f = 0.0
        self._cur_sse = 0.0
        self._mesh_size = 1.0
        self._sample_size = 2
        self._history = []
        self._controlv_f = []
        self._controlv_names = []
        self._x0 = []
        self._cur_x = []
        self._lb = []
        self._ub = []
        self._hypd = qt.Data(name='hypd')
        self._hypd.add_coordinate('iteration')
        self._iterations = 70


        self._optimization_is_on = False
        self._nfeval = 0

        return True
    def initialize(self):
        # keeps the control variables and their properties, but sets up for
        # another optimization

        self._found_better_point = False
        self._cur_f = 0.0
        self._cur_sse = 0.0
        self._n = len(self._controlv_names)
        self._beta = np.ones(self._n)*2.0
        self._r = np.ones(self._n)
        self._D = np.hstack((-np.eye(self._n),np.eye(self._n)))
        self._sample_size = 2
        self._history = []
        self._cur_x = []
        self._hypd = qt.Data(name='hypd')
        self._hypd.add_coordinate('iteration')
        self._hypd.add_value('current f')
        self._hypd.add_value('N fevals')
        self._hypd.add_value('min/max ratio')
        self._fed = qt.Data(name='fed')
        self._fed.add_coordinate('function evals')
        self._fed.add_value('current abs delta f')
        self._optimization_is_on = False
        self._nfeval = 0
        self._prev_r_max = 1.
        self.poll_size_initialization()
        self._cur_x = self._x0
        self._cur_f = self._evaluate_point(self._cur_x,1)[0]
        return

    def poll_size_initialization(self):
        self._delta = np.ones(self._n)
        for i in range(self._n):
            ub = np.array(self._ub)
            lb = np.array(self._lb)
            self._delta[i] = (ub[i] - lb[i])/10.0

        self._deltainit = np.copy(self._delta)
        return
    def update_poll_size_vector(self, success, direction):
        if success:
            print 'Expanding poll size vector'
            d_inf = np.linalg.norm(direction)
            for i in range(self._n):
                if np.abs(direction[i]) > 1.0/self._n * d_inf:
                    self._r[i] = self._r[i] + 1.0


        else:
            print 'Shrinking poll size vector.'
            self._r = self._r - 1.

                    # Otherwise, r_i stays the same
        # Now check if a mesh index is too small relative to the others
##        max_r = np.max(self._r)
##        for i in range(self._n):
##            if self._r[i] < -2.0 and self._r[i] < 2*max_r:
##                self._r[i] = self._r[i] + 1.
        return

    def compute_poll_size_vector(self):
        self._delta = self._deltainit*np.power(self._beta,self._r)
        return

    def iterate(self):
        self.compute_mesh_size_vector()
        if self._found_better_point:
            ret = self.search_step()
        else:
            self.poll_step()
        return

    def compute_mesh_size_vector(self):
        # Compute little delta (d) using Definition 2.2
        d = np.zeros(self._n)

        for i in range(self._n):
            d[i] = np.power(np.min([self._delta[i],self._deltainit[i]]),2.)/(np.sqrt(self._n)*self._deltainit[i])
        self._d = np.copy(d)
        return
    def search_step(self):
        # This step is optional
        # If the last search or poll step was successful, we'll try to move in the same
        # direction again
        self._found_better_point = False

        x_trial = self._cur_x + 3*np.dot(np.diag(self._d),self._last_direction.T).T

        is_feas = self._feasible_points(x_trial[None])

        # Re-evaluate the current point -- important for a noisy function
        self._cur_f = self._evaluate_point(self._cur_x,1)[0]
        print 'SEARCH STEP current x is %s, cur f is %s' % (self._cur_x, self._cur_f)
        # Now evaluate the feasible points only
        self._found_better_point = False
        if is_feas:
            f_evaluated = self._evaluate_point(x_trial,1)
            print 'SEARCH found f %s' % f_evaluated[0]
            if f_evaluated[0] < self._cur_f:
                # Found a successful point
                print 'Successful SEARCH, %f less than %f' % (f_evaluated[0], self._cur_f)
                self._cur_x = x_trial
                self._cur_f = f_evaluated[0]
                direction = self._last_direction
                #self._found_better_point = True
                return True
        return False
    def generate_basis(self):
        # We need to generate an orthonormal basis, so first we can generate a
        # random direction in the unit sphere (Muller, 1952)

        v = np.random.normal(0,1,self._n)
        # Normalize it
        v_norm = v/np.linalg.norm(v)
        # Compute Householder matrix. The [None] makes the array 2D so we can
        # transpose it.
        self._H = np.eye(self._n) - 2.0*np.dot(v_norm[None].T,v_norm[None])
        return

    def poll_step(self):
        # Generate a new basis
        self.generate_basis()
        # Update the poll and mesh size vectors
        self.compute_poll_size_vector()
        self.compute_mesh_size_vector()
        # Now round the normalized directions of the Householder matrix to the
        # grid, for this iteration
        a, b = self._H.shape
        Hdd = np.eye(self._n)
        for i in range(a):
            for j in range(b):
                Hdd[i,j] = np.round(self._delta[j]*self._H[i,j]/self._d[j])*self._d[j]
        # Now construct the maximal positive basis with -Hdd and Hdd concatenated
        H2n = np.array(np.hstack((-Hdd,Hdd)).T)
        # Construct an array of points by adding the maximal basis to the current
        # point.
        x_array = np.zeros((2*self._n,self._n))
        for i in range(2*self._n):
            x_array[i,:] = np.array(self._cur_x) + np.dot(np.diag(self._d),H2n[i,:].T).T
        # Check which points are feasible and create a vector that indicates
        # feasibility, so we don't try to evaluate infeasible (out of bounds)
        # values.

        is_feas = self._feasible_points(x_array)
        # Initialize direction; useful for when poll fails
        direction = H2n[1,:]
        # Re-evaluate the current point -- important for a noisy function
        self._cur_f = self._evaluate_point(self._cur_x,1)[0]
        print 'current x is %s, cur f is %s' % (self._cur_x, self._cur_f)
        # Now evaluate the feasible points only
        self._found_better_point = False
        for i in range(2*self._n):
            if is_feas[i]:
                print 'trying x %s' % x_array[i,:]
                f_evaluated = self._evaluate_point(x_array[i,:],1)
                print 'found f %s' % f_evaluated[0]
                if f_evaluated[0] < self._cur_f:
                    # Found a successful point
                    print 'Successful poll, %f less than %f' % (f_evaluated[0], self._cur_f)
                    self._cur_x = x_array[i,:]
                    self._cur_f = f_evaluated[0]
                    direction = H2n[i,:]
                    self._found_better_point = True
                    break
        print 'r is %s' % self._r
        # Now update the poll size vector, depending on if successful poll step
        self.update_poll_size_vector(self._found_better_point, direction)
        self._last_direction = direction
        return



    def set_max_sample_size(self, mss):
        self._max_sample_size = int(mss)
        return
    def get_max_sample_size(self):
        return self._max_sample_size

    def set_iterations(self, numits):
        self._iterations = int(numits)
        return
    def get_iterations(self):
        return self._iterations

    def get_current_x(self):
        return self._cur_x

    def set_x0(self, x0):
        # use for starting the optimizer from a different point than the one
        # defined when adding the original control variables
        self._x0 = np.float(x0)
        return
    def set_debug_on(self):
        self._debug = True
        print 'Debug set to true'
        return
    def set_debug_off(self):
        self._debug = False
        print 'Debug set to false'
        return


    def add_control_variable(self,name,f_handle, x0, lb, ub):
        self._controlv_names.append(name)
        self._controlv_f.append(f_handle)
        self._x0.append(np.float(x0))
        self._cur_x.append(np.float(x0))
        self._lb.append(np.float(lb))
        self._ub.append(np.float(ub))
        return
    def remove_control_variable(self,name):
        try:
            vidx = self._controlv_names.index(name)
        except:
            logging.warning(__name__ + ('control variable named %s not found' % name))
            return
        del self._controlv_names[vidx]
        del self._controlv_f[vidx]
        del self._x0[vidx]
        np.delete(self._cur_x, vidx)
        del self._lb[vidx]
        del self._ub[vidx]
        return


    def get_control_variables(self):
        return self._controlv_names

    def set_output_variable(self, fhandle):
        self._output_f = fhandle
        return

    def optimize(self):
        self.initialize()
        self.generate_basis()

        for i in range(int(self._iterations)):
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
            print 'Iteration %d, cur_f %.2f' % (i, self._cur_f)
            self.iterate()
            print 'Current x location is: %s, %s' % (np.array(self._cur_x), self._cur_f)
            self._evaluate_point(self._cur_x)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q": break
            time.sleep(0.01)
        print 'Done!'


        print 'Total function evals: %d' % self._nfeval
#        plotff = qt.plot(self._fed, coorddim=0, valdim=1, maxtraces=1, name='fevalp', clear=True)
        return np.array(self._cur_x)

    def _feasible_points(self, x_array):
        ub = np.array(self._ub)
        lb = np.array(self._lb)
        logical_index = np.array([True] * x_array.shape[0])
        for i in np.arange(x_array.shape[0]):
            x = x_array[i,:]
            isnotfeas = np.logical_or(lb >= x, ub <= x)
            if np.any(isnotfeas):
                # point is not feasible
                logical_index[i] = False
            else:
                # point is feasible
                logical_index[i] = True
        return logical_index


    def _evaluate_point(self, x, n=1):
        output_array = np.zeros(n)
        # set each output variable to the desired value
        for i in range(self._n):
                self._controlv_f[i](x[i])
        # sleep before measuring
        time.sleep(self._sleeptime)
        # measure n times -- useful for hypothesis testing
        for q in range(n):
            output_array[q] = self._output_f()
        # return the array of output values measured
        self._nfeval = self._nfeval + n
        return output_array


    def _logscale(self):
        # Log scaling from Numerical Optimization
        # Somewhat incomplete, not all cases considered.
        lb = np.array(self._lb)
        ub = np.array(self._ub)
        eps_l = 1.0e-6
        lf = np.logical_and(np.abs(lb) >= eps_l, np.abs(ub) <= eps_l)
        lu = np.logical_and(np.abs(lb) >= eps_l, np.abs(ub) >= eps_l)
        uf = np.logical_and(np.abs(lb) <= eps_l, np.abs(ub) <= eps_l)

        log_s = np.zeros(lb.size)
        log_s[lf] = np.log2(np.abs(lb[lf]))
        log_s[lu] = (np.log2(np.abs(lb[lu])) + np.log2(np.abs(ub[lu])))/2.0
        log_s[uf] = np.log2(np.abs(ub[uf]))
        self._scale = np.power(2.0,np.floor(log_s))*0.2 # added 0.2 to reduce the size a bit of the initial scale
        self._scale = np.abs(ub-lb)*0.1
        print 'SCALE is %s' % self._scale
        return


