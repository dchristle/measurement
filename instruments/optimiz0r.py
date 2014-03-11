# Line search optimizer using FSM/XPS
#
#
# Taken from original optimiz0r.py file written by Wolfgang P.; modified
# to use FSM and replace the FBL_MainDC2.vi
#
# Right now, the optimizer just executes a line search along each axis and
# doesn't fit any Gaussians like the existing VI; I actually think in some
# or most cases this might be the simpler, better way.

from instrument import Instrument
import types
import qt
import msvcrt

class optimiz0r(Instrument):

    dimension_sets = {
            'default' : {
                'x' : {
                    'scan_length' : 2.,
                    'nr_of_points' : 31,#99,
#                    'pixel_time' : 50,
                    },
                'y' : {
                    'scan_length' : 2.,
                    'nr_of_points' : 31,#99,
#                    'pixel_time' : 50,
                    },
                'z' : {
                    'scan_length' : 4.,
                    'nr_of_points' : 31,#99,
#                    'pixel_time' : 50,
                    },
                'xxps' : {
                    'scan_length' : 2.,
                    'nr_of_points' : 31,
                    },
                'yxps' : {
                    'scan_length' : 2.,
                    'nr_of_points' : 31,
                    },
                'zyx' : ['z','y','x'],
                'xyonly':['y','x'],
                'fullxyz' : ['z','y','x','yxps','xxps']
                },
            }

    def __init__(self, name, opt1d_ins=qt.instruments['opt1d_counts'],
            fsm_ins=qt.instruments['fsm'],
            xps_ins=qt.instruments['xps'],
            dimension_set='default'):
        Instrument.__init__(self, name)

        self.add_function('optimize')
        self.opt1d_ins = opt1d_ins
        self.dimensions = self.dimension_sets[dimension_set]

        self.mos = mos_ins

    def optimize(self, cycles=1, cnt=1, int_time=50, dims=[], order='zyx'):
        ret=True
        for c in range(cycles):

            if len(dims) == 0:
                dims = self.dimensions[order]

            for d in dims:
                ret=ret and self.opt1d_ins.run(dimension=d, counter = cnt,
                        pixel_time=int_time, **self.dimensions[d])
                qt.msleep(1)
            if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break


        return ret
