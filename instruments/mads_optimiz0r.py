# Pattern search optimizer using FSM/XPS
# David Christle <christle@uchicago.edu>, September 2014
#
# This is a pattern search based optimizer that specifically uses
# the dynamic mesh-adaptive direct search algorithm to compensate for
# drift in the cryostat. The idea is to poll randomized and scaled directional
# moves in three dimensions and to move if a better optimum is detected. The
# algorithm adjusts to the relative 'scale' of the axes as the algorithm 
# progresses so that each iteration is making meaningful changes to
# the control variables (not too large, not too small) to rapidly
# locate an optimum.

from instrument import Instrument
import types
import qt
import msvcrt
import numpy as np
from analysis.lib.fitting import fit,common
reload(fit)
reload(common)
class mads_optimiz0r(Instrument):



    def __init__(self, name, dimension_set='default'):
        Instrument.__init__(self, name)
        self.dimension_sets = {
            'default' : {
                'x' : {
                    'scan_length' : 1.0,
                    'nr_of_points' : 60,
                    'qt_ins' : 'fsm',
                    'channel' : 'X',
                    'sigma' : 0.5

                    },
                'y' : {
                    'scan_length' : 1.0,
                    'nr_of_points' : 60,
                    'qt_ins' : 'fsm',
                    'channel' : 'Y',
                    'sigma' : 0.5
                    },
                'z' : {
                    'scan_length' : 1.4,
                    'nr_of_points' : 60,
                    'qt_ins' : 'xps',
                    'channel' : 'Z',
                    'sigma' : 2.0/1000.0
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
        self._opt_pos = {'x' : self._fsm.get_abs_positionX(),
                         'y' : self._fsm.get_abs_positionY(),
                         'z' : self._xps.get_abs_positionZ()}
        self.add_function('optimize')
        self.dimensions = self.dimension_sets[dimension_set]


    def optimize(self, plot=0, dims='xyz', cycles=1):
        if 'fbl_plot' in qt.plots:
            qt.plots['fbl_plot'].clear()
        #qt.plot(name='fbl_plot',clear=True,needtempfile=False)

        for c in range(cycles):
            ret = True


            # Get the current position from the FSM - only works if the
            # FSM has been written to at least once!
            cur_pos = self._fsm.get_abs_positionX()
            
            if cur_pos == None:
                print 'Current position is none! FSM not initialized?'
                break

            self._opt_pos_prev['x'] = cur_pos
            cur_pos = self._fsm.get_abs_positionY()
            if cur_pos == None:
                print 'Current position is none! FSM not initialized?'
                break              
            self._opt_pos_prev['y'] = cur_pos                
            cur_pos = self._xps.get_abs_positionZ() # In mm, not um!
            self._opt_pos_prev['z'] = cur_pos # In mm
            
            prev_ctr0_src = self._ni63.get_ctr0_src()
            self._ni63.set_ctr0_src('PFI0')
            prev_count_time = self._ni63.get_count_time()
            count_rate = 10.0 # Hz
            self._ni63.set_count_time(1.0/count_rate)
            #if not 'gps_fbl' in qt.instruments:
            # GPS/MADS optimizer does not exist, create it
            self.gps = qt.instruments.create('gps_fbl','gps_optimiz0r')
            
            self.gps.set_output_variable(lambda : -1.0*count_rate*self._ni63.get_ctr0())
            
            x_low = self._opt_pos_prev['x']-self.dimensions['x']['scan_length']/2.0
            x_high = self._opt_pos_prev['x']+self.dimensions['x']['scan_length']/2.0
            
            self.gps.add_control_variable('fsm_x',self._fsm.set_abs_positionX,self._opt_pos_prev['x'],x_low,x_high)
            
            y_low = self._opt_pos_prev['y']-self.dimensions['y']['scan_length']/2.0
            y_high = self._opt_pos_prev['y']+self.dimensions['y']['scan_length']/2.0
            
            self.gps.add_control_variable('fsm_y',self._fsm.set_abs_positionY,self._opt_pos_prev['y'],y_low,y_high)
            
            z_low = self._opt_pos_prev['z']-0.001*self.dimensions['z']['scan_length']/2.0
            z_high = self._opt_pos_prev['z']+0.001*self.dimensions['z']['scan_length']/2.0
            
            self.gps.add_control_variable('xps_z',self._xps.set_abs_positionZ,self._opt_pos_prev['z'],z_low,z_high)
            self.gps.set_max_sample_size(2)
            self.gps.set_iterations(10)
            new_opt = self.gps.optimize()
            self._opt_pos['x'] = new_opt[0]
            self._opt_pos['y'] = new_opt[1]
            self._opt_pos['z'] = new_opt[2]

            print 'Previous: (%.3f, %.3f, %.3f), new optimum: (%.3f, %.3f, %.3f) (delta (%.1f, %.1f, %.1f) nm)' % (self._opt_pos_prev['x'], self._opt_pos_prev['y'], self._opt_pos_prev['z'], self._opt_pos['x'], self._opt_pos['y'], self._opt_pos['z'], self._opt_pos['x']*1.0E3-self._opt_pos_prev['x']*1.0E3, self._opt_pos['y']*1.0E3-self._opt_pos_prev['y']*1.0E3, self._opt_pos['z']*1.0E6-self._opt_pos_prev['z']*1.0E6)
            self._ni63.set_count_time(prev_count_time)
            self._ni63.set_ctr0_src(prev_ctr0_src)



            qt.msleep(0.1)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break


        return ret
