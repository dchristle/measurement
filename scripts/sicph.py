import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2


class SiCPH_Master(m2.Measurement):

    mprefix = 'SiC_Son'

    def prepare(self):
        # Set up some instruments
        self._ph = qt.instruments['ph']
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']

        # Prepare instruments for measurement and verify FBL output

        # Configure PH for histogram mode, initialize
        self._ph.start_histogram_mode()

        self._ph.set_Range(self.params['Range'])
        self._ph.set_CFDLevel0(self.params['CFDLevel0'])
        self._ph.set_CFDLevel1(self.params['CFDLevel1'])
        self._ph.set_CFDZeroCross0(self.params['CFDZeroCross0'])
        self._ph.set_CFDZeroCross0(self.params['CFDZeroCross1'])
        print 'PicoHarp settings configured...'
        self._fbl.optimize()
        # Use some logic here to decide what's going on
        # i.e. position, chi sq., signal amp, background amp
        print 'FBL optimized...'
        return
    def measure(self):



##        x = np.linspace(0,self.params['xmax'],self.params['xpts'])
##        y = np.linspace(0,self.params['ymax'],self.params['ypts'])
##        z = np.zeros((self.params['xpts'],self.params['ypts']))
##
##        for i,xval in enumerate(x):
##            print 'linesweep %d / %d ...' % (i+1, self.params['xpts'])
##            for j,yval in enumerate(y):
##                qt.msleep(0.01)
##                z[i,j] = xval*yval
##
##        # save the data into the pre-created group.
##        # note the passed meta-data (optional).
##        # you can have a look at the data with HDFView
##        # (you can get it from hdfgroup.com)
##        grp = h5.DataGroup('xy-scan', self.h5data, base=self.h5base)
##        grp.add('x', data=x, unit='um', note='somewhat inaccurate')
##        grp.add('y', data=y, unit='um')
##        grp.add('z', data=z, unit='counts per second', dimensions='1=x, 2=y')

        return



# measurement parameters

xsettings = {
        'int_time' : 30.0,
        'pos_displacement' : 2.5,
        'fbl_time' : 30.0,
        'CFDLevel0' : 200,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 200,
        'CFDZeroCross1' : 10,
        'Range' : 0,
        'Offset' : 0,
        'SyncDiv' : 1
        }


# Create a measurement object m
m = SiCPH_Master('first measurements')

# since params is not just a dictionary, it's easy to incrementally load
# parameters from multiple dictionaries
# this could be very helpful to load various sets of settings from a global
# configuration manager!
m.params.from_dict(xsettings)


if m.review_params():
    print 'Proceeding with measurement ...'
    m.prepare()
    m.measure()
    m.save_params()
    m.save_stack()
else:
    print 'Measurement aborted!'

# important! hdf5 data must be closed, otherwise will not be readable!
# (can also be done by hand, of course)
m.finish()