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

class gps_optimiz0r(Instrument):
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


        # Constants for the pattern search
        self._sleeptime = 0.25
        self._alpha = 0.05 # significance test alpha
        self._max_sample_size = 800
        self._mesh_expansion_factor = 2.0
        self._mesh_shrinkage_factor = 0.5
        self._mads = False

        self._debug = False
        self._cur_f = 0.0
        self._nfeval = 0
        self._cur_sse = 0.0
        self._mesh_size = 1.0
        self._sample_size = 2
        self._max_sample_size = 8
        self._history = []
        self._controlv_f = []
        self._controlv_names = []
        self._x0 = []
        self._cur_x = []
        self._lb = []
        self._ub = []
        self._hypd = qt.Data(name='hypd')
        self._hypd.add_coordinate('iteration')
        self._hypd.add_value('T4')
        self._hypd.add_value('T4 crit')
        self._hypd.add_value('sample size')
        self._hypd.add_value('angle 1')
        self._hypd.add_value('angle 2')
        self._hypd.add_value('cmax position')
        self._hypd.add_value('current f')
        self._hypd.add_value('N fevals')

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
        self._hypd.add_value('T4')
        self._hypd.add_value('T4 crit')
        self._hypd.add_value('sample size')
        self._hypd.add_value('angle 1')
        self._hypd.add_value('angle 2')
        self._hypd.add_value('cmax position')
        self._hypd.add_value('current f')
        self._hypd.add_value('N fevals')
        self._optimization_is_on = False
        self._nfeval = 0

        return True
    def initialize(self):
        # keeps the control variables and their properties, but sets up for
        # another optimization
        self._cur_f = 0.0
        self._cur_sse = 0.0
        self._mesh_size = 1.0
        self._sample_size = 2
        self._history = []
        self._cur_x = []
        self._hypd = qt.Data(name='hypd')
        self._hypd.add_coordinate('iteration')
        self._hypd.add_value('T4')
        self._hypd.add_value('T4 crit')
        self._hypd.add_value('sample size')
        self._hypd.add_value('angle 1')
        self._hypd.add_value('angle 2')
        self._hypd.add_value('cmax position')
        self._hypd.add_value('current f')
        self._hypd.add_value('N fevals')
        self._fed = qt.Data(name='fed')
        self._fed.add_coordinate('function evals')
        self._fed.add_value('current abs delta f')
        self._optimization_is_on = False
        self._nfeval = 0

    def set_max_sample_size(self, mss):
        self._max_sample_size = int(mss)
        return
    def get_max_sample_size(self):
        return self._max_sample_size

    def get_current_x(self):
        return self._cur_x

    def set_x0(self, x0):
        # use for starting the optimizer from a different point than the one
        # defined when adding the original control variables
        self._x0 = x0
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
        self._x0.append(x0)
        self._cur_x.append(x0)
        self._lb.append(lb)
        self._ub.append(ub)
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
    def _feasible_points(self, x_array):
        ub = np.array(self._ub)
        lb = np.array(self._lb)
        logical_index = np.array([True] * x_array.shape[0])
        for i in np.arange(x_array.shape[0]):
            x = x_array[i,:]
            isnotfeas = np.logical_or(np.abs(lb) >= x, np.abs(ub) <= x)
            if np.any(isnotfeas):
                # point is not feasible
                logical_index[i] = False
            else:
                # point is feasible
                logical_index[i] = True
        return x_array[logical_index,:]

    def iterate(self):
        # set the overall scale for each dimension
        if not self._optimization_is_on:
            self._logscale()
            self._optimization_is_on = True
            self._cur_x = self._x0
            self._cur_f = self._evaluate_point(self._cur_x,self._sample_size)[0]
            if self._debug:
                print 'Initial f is %.2f' % self._cur_f
            self._iteration_number = 0
        else:
            self._iteration_number = self._iteration_number + 1

        # re-poll the current location -- necessary for unbiasing the
        # hypothesis testing
        cur_output = self._evaluate_point(self._cur_x,self._sample_size)
        self._cur_f = cur_output[0]
        self._cur_sse = cur_output[1]
        # generate a set of points to poll
        x_basis = self._gps2nbasis(self._cur_x)
        # now check if each point is feasible, i.e. within bounds
        # this is important since gps2n can generate trial points outside
        # the bounds
        if self._debug:
            print 'x_basis is %s' % x_basis
        new_x = self._feasible_points(x_basis)

        # now iterate through the array, evaluating each point
        if self._debug:
            print 'current x is: %s' % self._cur_x
            print 'new x shape is: %s' % new_x
        if new_x.shape[0] == 0:
            self._mesh_size = self._mesh_size*self._mesh_shrinkage_factor
            if self._debug:
                print 'no feasible points, shrinking mesh'
                print 'new x is %s' % new_x
            self._history.append(self._cur_f)
            return
        new_f = []
        new_sse = []
        for i in range(new_x.shape[0]):
            point_output = self._evaluate_point(new_x[i,:],self._sample_size)
            new_f.append(point_output[0])
            new_sse.append(point_output[1])
        if self._debug:
            print 'new f array is: %s' % np.array(new_f)
        new_f = np.array(new_f)
        min_f = np.amin(new_f)
        if self._debug:
            print 'min f is %s' % min_f
        if np.array(self._x0).size > 1:
            # test for statistical significance
            sse = np.float(self._cur_sse) + np.sum(np.array(new_sse))
            n = 1.0 + np.float(new_x.shape[0])
            M = (n+1.0)*np.float(self._sample_size)
            mse = sse/(M-(n+1))
            Y1 = np.amin(np.append(new_f,self._cur_f))
            YN1 = np.amax(np.append(new_f,self._cur_f))
            T4 = (YN1-Y1)/np.sqrt(mse/np.float(self._sample_size))
            q_crit = qst.qsturng(1-self._alpha, n+1, M-(n+1))
            self._hypd.add_data_point(self._iteration_number, T4, q_crit, self._sample_size, self._cur_x[0]*180.0/np.pi, self._cur_x[1]*180.0/np.pi, self._cur_x[2], self._cur_f, self._nfeval)
            self._fed.add_data_point(self._nfeval, np.log(np.abs(self._cur_f + 600)))
            if self._debug:
                print 'mse, M, n %.2f %.2f %.2f' % (mse, M, n)
                print 'T4 statistic is: %.3f, crit. is %.2f' % (T4, q_crit)
            if q_crit > T4:
                self._sample_size = np.min(np.array((2*self._sample_size, self._max_sample_size)))
                self._history.append(self._cur_f)
                return False



        if min_f < self._cur_f:
            if self._debug:
                print 'successful iteration, %.2f to %.2f' % (self._cur_f, min_f)
                print 'new x selected %s' % new_x[np.argmin(new_f),:]
            self._cur_x = new_x[np.argmin(new_f),:]
            self._cur_f = new_f[np.argmin(new_f)]
            self._mesh_size = self._mesh_size*self._mesh_expansion_factor
            self._history.append(self._cur_f)
            return True
        else:
            if self._debug:
                print 'no successful polls, shrinking mesh'
            self._mesh_size = self._mesh_size*self._mesh_shrinkage_factor
            self._history.append(self._cur_f)
            return False

        return

    def optimize(self):
        self.initialize()
        for i in range(40):
            print 'iteration %d, cur_f %.2f' % (i, self._cur_f)
            self.iterate()
            print 'Current x location is: %s, %s' % (np.array(self._cur_x), self._cur_f)
            time.sleep(0.1)
        print 'Done!'

        self._optimization_is_on = False
        hist = np.array(self._history)
        qt.plot(np.arange(hist.size),hist,name='optimizegps',clear=True)
        plot2d = qt.plot(self._hypd, coorddim=0, valdim=1,name='optT4',clear=True)
        plot2d.add_data(self._hypd, coorddim=0, valdim=2, maxtraces=1)
        plot2d.add_data(self._hypd, coorddim=0, valdim=3, maxtraces=1)

        plota1 = qt.plot(self._hypd, coorddim=0, valdim=4,name='angle1',clear=True)
        plota2 = qt.plot(self._hypd, coorddim=0, valdim=5,name='angle2',clear=True)
        plotcmax = qt.plot(self._hypd, coorddim=0, valdim=6,name='cmax',clear=True)
        plotfevals = qt.plot(self._hypd, coorddim=0, valdim=8, maxtraces=1, clear=True)
        print 'Total function evals: %d' % self._nfeval
        plotff = qt.plot(self._fed, coorddim=0, valdim=1, maxtraces=1, name='fevalp', clear=True)
        return



    def _evaluate_point(self, x, n=1):
        output_array = np.zeros(n)
        ub = np.array(self._ub)
        # set each output variable to the desired value
        for i in range(ub.size):
                self._controlv_f[i](x[i])
        # sleep before measuring
        time.sleep(self._sleeptime)
        # measure n times -- useful for hypothesis testing
        for q in range(n):
            output_array[q] = self._output_f()
        # return the array of output values measured
        self._nfeval = self._nfeval + n
        return (np.mean(output_array), np.sum((output_array - np.mean(output_array))**2.0))


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
        self._scale = np.power(2.0,np.floor(log_s))
        return

    def _gps2nbasis(self, x):
        # Detect if bounds are active, then create the maximal 2n basis
        # In principle, other bases could also be generated here

        # Define a tolerance -- right now, it's relative to the overall scale
        tol = 0.001*self._scale
        # Detect which bounds are active
        lb = np.array(self._lb)
        ub = np.array(self._ub)
        if self._debug:
            print 'tol %s' % tol
            print 'lb %s' % lb
            print 'ub %s' % ub
            print 'x %s' % x
        lb_a = np.logical_or(np.abs(x - lb) <= tol , np.abs(x-ub) <= tol)
        # Create an identity matrix, and then select which directions are
        # feasible (basis vectors for the search) and which directions are
        # not feasible (they form the tangent cone).
        I_m = np.eye(lb.size)
        self.basis = I_m[:,np.logical_not(lb_a)]
        self.tc = I_m[:,lb_a]
        if self._mads:
            # see Audet -- create a lower triangular of random signs, scaled by
            # the *rounded* poll size parameter. Then create the diagonal
            # scaled by just the poll parameter (not necessarily an integer).
            # This follows section 4.1 in Audet (2004).
            self._poll_size = np.sqrt(self._mesh_size)
            lowerT = np.tril(np.round(self._poll_size+1)*np.random.rand(lb.size,lb.size)-0.5,-1)
            diag = np.diagflat(self._poll_size*np.sign(np.random.rand(lb.size,1)-0.5))
            self.basis = lowerT + diag
            # now randomize the ordering of the basis vectors
            # there has to be a faster way to do this with indexing than the for
            # loop I am currently using
            scramble_vecA = np.random.permutation(lb.size)
            scramble_vecB = np.random.permutation(lb.size)
            output_basis = np.eye(lb.size)
            for i in range(lb.size):
                for j in range(lb.size):
                    output_basis[i,j] = self.basis[scramble_vecA[i],scramble_vecB[j]]
            # output_basis should now be "d" in Audet (2004) -- not yet a complete
            # n+1 set, but we'll do that in a moment
            self._basis = output_basis
        if self._mads:
            n1vector = -1.0*np.sum(self.basis,1)
            self._basis = np.vstack((n1vector,self.basis.T)).T
            Nbasis = self.basis.shape[1]
            Ntangent = np.sum(lb_a)
            Ntotal = Nbasis + Ntangent
            dirMat = np.hstack((self.basis, self.tc))
            indexV = np.concatenate( (np.arange(Nbasis), np.arange(Nbasis,Ntotal), np.arange(Nbasis,Ntotal)) )
            signV = np.concatenate( (np.ones(Nbasis), np.ones(Ntangent), -1*np.ones(Ntangent) ) );
        else:
            # Figure out the number of basis vectors/vectors within the tangent cone
            Nbasis = np.sum(np.logical_not(lb_a))
            Ntangent = np.sum(lb_a)
            Ntotal = Nbasis + Ntangent
            # Assemble the matrix of directions (both basis and tangent cone)
            dirMat = np.hstack((self.basis, self.tc))
            # Index them, take care of signs for 2N basis
            indexV = np.concatenate( (np.arange(Nbasis), np.arange(Nbasis), np.arange(Nbasis,Ntotal), np.arange(Nbasis,Ntotal)) )
            signV = np.concatenate( (np.ones(Nbasis), -1*np.ones(Nbasis), np.ones(Ntangent), -1*np.ones(Ntangent) ) );
            # Create a vector for ordering the polling process -- not relevent
            # yet since we will do complete polling

        orderV = np.arange(indexV.size)
        # np.random.shuffle(orderV)
        # Iterate to create a list of polling locations
        if self._debug:
            print '%s, %s' % (indexV.size, lb.size)
        new_x = np.zeros((indexV.size,lb.size))
        if self._debug:
            print 'mesh scale %s' % (self._mesh_size*self._scale)
        for i in range(indexV.size):
            direction = signV[i]*dirMat[:,indexV[orderV[i]]]
            new_x[i,:] = x + self._mesh_size*self._scale*direction
        return new_x





##    def __init__(self, name, mos_ins=qt.instruments['master_of_space'],
##            adwin_ins=qt.instruments['adwin']):
##        Instrument.__init__(self, name)
##
##        self.add_function('optimize')
##        self.mos = mos_ins
##        self.adwin= adwin_ins
##
##        ins_pars  = {'a'  :   {'type':types.FloatType,'val':1.5,'flags':Instrument.FLAG_GETSET},
##                     'b'  :   {'type':types.FloatType,'val':1.3,'flags':Instrument.FLAG_GETSET},
##                     'c'  :   {'type':types.FloatType,'val':0.5,'flags':Instrument.FLAG_GETSET},
##                     'd'  :   {'type':types.FloatType,'val':0.9,'flags':Instrument.FLAG_GETSET},
##                     }
##        instrument_helper.create_get_set(self,ins_pars)
##
##
##
##    def _measure_at_position(self, pos, int_time, cnt, speed):
##        self.mos.move_to_xyz_pos(('x','y','z'),pos,speed,blocking=True)
##        qt.msleep(0.001)
##        return self.adwin.measure_counts(int_time)[cnt-1]
##
##    def _tetra_volume(self, V):
##        return 1/6.*np.abs(np.linalg.det(V[0:3]-V[3]))
##
##    def optimize(self,xyz_range=[.5,0.5,1.0], xyz_tolerance_factor=0.002, max_cycles=15,
##                  cnt=1, int_time=50, speed=2000,
##                  do_final_countrate_check=True):
##
##        #in principle, the method below works for any dimension D, and number of vertices N!
##        #Only need to change the initial simplex shape to size N and set_functions to dim D
##        D=3
##        N=5
##
##
##        a=self._a
##        b=self._b
##        c=self._c
##        d=self._d
##        search_range=np.array(xyz_range)
##        tolerance=search_range*xyz_tolerance_factor
##        pos = np.array([self.mos.get_x(), self.mos.get_y(), self.mos.get_z()])
##        old_cnt = self.adwin.measure_counts(int_time)[cnt-1]
##        new_pos=pos.copy()
##
##        tetra=np.array([[ 0, 0, 0],
##                        [ 1, 1, 1],
##                        [ 1,-1,-1],
##                        [-1, 1,-1],
##                        [-1,-1, 1],
##                      ],dtype=np.float)
##
##        V=pos+tetra*search_range
##
##        J=np.zeros(N)
##
##        for i in range(1,len(J)):
##          J[i]=self._measure_at_position(V[i],int_time,cnt,speed)
##        print 'Simplex J countrates:',J/(int_time/1000.)
##
##        for j in range(max_cycles):
##          print '\n ============== \n j: ', j
##          J[0]=self._measure_at_position(V[0],int_time,cnt,speed)
##          J_o=np.sort(J)
##          V_o=V[np.argsort(J)]
##
##          v_c=1./(N-1.)*np.sum(V_o[1:],axis=0)
##          v_r=(1-a)*V_o[0]+a*v_c
##
##          j_r=self._measure_at_position(v_r,int_time,cnt,speed)
##          print 'New vertex j_r countrates:',j_r/(int_time/1000.)
##
##          if J_o[1]<=j_r<=J_o[-1]:
##            V_o[0]=v_r
##            print 'case 1: replace'
##          elif j_r>J_o[-1]:
##            v_e=b*v_r+(1-b)*v_c
##            j_e=self._measure_at_position(v_e,int_time,cnt,speed)
##            if j_e>J_o[-1]:
##              V_o[0]=v_e #replace [0] or V[-1] here??
##              print 'case 2.1: expand far'
##            else:
##              V_o[0]=v_r #replace [0] or V[-1] here??
##              print 'case 2.2: expand close'
##          elif j_r<J_o[1]:
##            v_s=c*v_r+(1-c)*v_c
##            j_s=self._measure_at_position(v_s,int_time,cnt,speed)
##            if j_s>=J[0]:
##              V_o[0]=v_s
##              print 'case 3.1: stay close'
##            else:
##              print 'case 3.2: shrink'
##              V_o=d*V_o+(1-d)*V_o[-1]
##          else:
##            print 'j_r comparison error. Check values.'
##
##          V=V_o
##          J=J_o
##
##          var=np.var(V,axis=0)
##          print 'Average xyz_var', np.sum(var)
##          if (var<tolerance).all():
##            new_pos = V[-1] if j_r<J_o[-1] else  V[0]
##            print '========================='
##            print 'Convex success'
##            print '========================='
##            break
##
##          if (msvcrt.kbhit() and (msvcrt.getch() == 'q')):
##            do_final_countrate_check = False
##            print 'User interrupt'
##            break
##
##          #vol=self._tetra_volume(V) #N=4 simplex volume.
##          #print 'vol', vol
##          #if vol[j] < 0.0002:
##          #  break
##
##          #----------------------------------------------------
##
##
##        if do_final_countrate_check:
##          print "Proposed position x change %d nm" % \
##                        (1000*V[-1][0]-1000*pos[0])
##          print "Proposed position y change %d nm" % \
##                        (1000*V[-1][1]-1000*pos[1])
##          print "Proposed position z change %d nm" % \
##                        (1000*V[-1][2]-1000*pos[2])
##          new_cnt = self._measure_at_position(V[-1], int_time, cnt, speed/2.)
##          print 'Old countrates', old_cnt/(int_time/1000.)
##          print 'New countrates', new_cnt/(int_time/1000.)
##          if new_cnt>old_cnt:
##            print 'New position accepted'
##          else:
##            print 'Old position kept'
##            self.mos.move_to_xyz_pos(('x','y','z'),pos,speed,blocking=True)
##
##        else:
##          print "Position x changed %d nm" % \
##                          (1000*new_pos[0]-1000*pos[0])
##          print "Position y changed %d nm" % \
##                          (1000*new_pos[1]-1000*pos[1])
##          print "Position z changed %d nm" % \
##                          (1000*new_pos[2]-1000*pos[2])
##          print 'Old countrates', old_cnt/(int_time/1000.)
##          print  "Countrates at new position: %d" % \
##                      (float(self._measure_at_position(new_pos, int_time,cnt, speed))/(int_time/1000.))
##
