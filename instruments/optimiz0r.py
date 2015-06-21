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
from analysis.lib.fitting import fit,common
reload(fit)
reload(common)
class optimiz0r(Instrument):



    def __init__(self, name, dimension_set='default'):
        Instrument.__init__(self, name)
        self.dimension_sets = {
            'default' : {
                'x' : {
                    'scan_length' : 2.0,
                    'nr_of_points' : 60,
                    'qt_ins' : 'fsm',
                    'channel' : 'X',
                    'sigma' : 0.35

                    },
                'y' : {
                    'scan_length' : 2.0,
                    'nr_of_points' : 60,
                    'qt_ins' : 'fsm',
                    'channel' : 'Y',
                    'sigma' : 0.35
                    },
                'z' : {
                    'scan_length' : 4.0,
                    'nr_of_points' : 50,
                    'qt_ins' : 'xps',
                    'channel' : 'Z',
                    'sigma' : 1.1/1000.0
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
                                'z' : 2.0/1000.0
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
                    ret = self.process_fit(point_array, cps, 'x')
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
                    ret = self.process_fit(point_array, cps, 'y')
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
                    ret = self.process_fit(point_array, cps, 'z')
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


                qt.msleep(1)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break


        return ret

    def process_fit(self, p, cr, dimension):

        # p is position array
        # cr is countrate array
        # d is a string corresponding to the dimension

        # Get the Gaussian width sigma guess from the dimensions dictionary
        sigma_guess = self.dimensions[dimension]['sigma']
        # Always do a Gaussian fit, for now
        self._gaussian_fit = True

        if self._gaussian_fit:
            # Come up with some robust estimates, for low signal/noise conditions
            a_guess = np.array(cr).min()
            A_guess = np.array(cr).max()-np.array(cr).min()
            p_size = np.size(p)
            x0_guess = p[np.round(p_size/2.0)] #np.sum(np.array(cr) * np.array(p))/np.sum(np.array(cr)**2)
            #sigma_guess = 1.0#np.sqrt(np.sum(((np.array(cr)-x0_guess)**2)*(np.array(p))) / np.sum((np.array(p))))
            #print 'Guesses: %r %r %r %r' % (a_guess,  x0_guess,A_guess,sigma_guess)
            # The following statement is that if the old sigma is within a factor of 2 of the initial fixed value
            # then use the old sigma as the guess for the fit.
            if np.abs(self._opt_sigma_prev[dimension]/self.dimensions[dimension]['sigma']) < 2.0 and np.abs(self._opt_sigma_prev[dimension]/self.dimensions[dimension]['sigma']) > 0.3:
                sigma_guess = self._opt_sigma_prev[dimension]
            else:
                sigma_guess = self.dimensions[dimension]['sigma']

            gaussian_fit = fit.fit1d(np.array(p,dtype=float), np.array(cr,dtype=float),common.fit_gauss, a_guess,
                    x0_guess,A_guess, sigma_guess, do_print=False,ret=True)

            ababfunc = a_guess + A_guess*np.exp(-(p-x0_guess)**2/(2.0*sigma_guess**2))

            if type(gaussian_fit) != dict:
                pamax = np.argmax(cr)
                self._opt_pos[dimension] = p[pamax]
                ret = False
                print '(%s) fit failed! Set to maximum.' % dimension


            else:

                if gaussian_fit['success'] != False:
                    self._fit_result = [gaussian_fit['params'][1],
                            gaussian_fit['params'][2],
                            gaussian_fit['params'][3],
                            gaussian_fit['params'][0] ]
                    self._fit_error = [gaussian_fit['error'][1],
                            gaussian_fit['error'][2],
                            gaussian_fit['error'][3],
                            gaussian_fit['error'][0] ]
                    self._opt_pos[dimension] = self._fit_result[0]
                    self._opt_sigma_prev[dimension] = self._fit_result[2]
                    ret = True
                    final_fit_func = gaussian_fit['params'][0] + gaussian_fit['params'][2]*np.exp(-(p-gaussian_fit['params'][1])**2/(2.0*gaussian_fit['params'][3]**2))
                    qt.plot(p,cr,p,ababfunc,p,final_fit_func,name='fbl_plot',needtempfile=False)
                    #print '(%s) optimize succeeded!' % self.get_name()


                else:
                    self.set_data('fit', zeros(len(p)))
                    ret = False
                    print '(%s) optimize failed! Set to maximum.' % dimension
                    self._opt_pos[dimension] = p[np.argmax(cr)]
                    qt.plot(p,cr,p,ababfunc,name='fbl_plot',clear=True,needtempfile=False)




        return ret
