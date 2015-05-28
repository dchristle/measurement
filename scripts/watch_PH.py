import time
import msvcrt

ph=qt.instruments['pm']

qt.mstart()
data = qt.Data(name='testmeasurement')

# Now you provide the information of what data will be saved in the
# datafile. A distinction is made between 'coordinates', and 'values'.
# Coordinates are the parameters that you sweep, values are the
# parameters that you readout (the result of an experiment). This
# information is used later for plotting purposes.
# Adding coordinate and value info is optional, but recommended.
# If you don't supply it, the data class will guess your data format.
data.add_coordinate('time')
data.add_value('counts')

# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.

# Next two plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.

import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
import pylab as plt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
import gc
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

# Defines a class called SiC_PL6PH_Master, DC 2015















# Above we defined a new ESR measurement class. Here we're going to create a
# a dictionary of measurement-specific parameters, create some of these ESR
# measurement class objects, and run them.


# measurement parameters - these are all references above.

xsettings = {
        'cycle_waveforms' : True, # set True to cycle waveform_list waveforms
        'waveform_list' : np.array((1,2,3,4,5)),
        'ctr_term' : 'PFI0',
        'wavelength' : 1036.5, # nm
        'wavelength_speed' : 300.0, # nm/min grating speed
        'dwell_time' : 500.0, # ms
        'ctr_term' : 'PFI0', # counter terminal for counting
        'MeasCycles' : 1500,
        'CFDLevel0' : 125,
        'CFDZeroCross0' : 10,
        'CFDLevel1' : 125,
        'CFDZeroCross1' : 10,
        'Binning' : 5,
        'Offset' : 0,
        'SyncDiv' : 1,
        'SyncOffset' : 0,
        'AcqTime' : 1, # PicoHarp acquisition time in seconds
        'bin_low' : 29850,
        'bin_high' : 30000,
        }


snspd = qt.instruments['snspd']
ph = qt.instruments['ph']



# Configure PH for histogram mode, initialize it with the correct settings
ph.start_histogram_mode()

ph.set_Binning(xsettings['Binning'])
ph.set_InputCFD0(xsettings['CFDLevel0'],xsettings['CFDZeroCross0'])
ph.set_InputCFD1(xsettings['CFDLevel1'],xsettings['CFDZeroCross1'])
ph.set_SyncOffset(xsettings['SyncOffset'])
print 'PicoHarp settings configured.'



plot2d = qt.Plot2D(data, 'b-', linewidth=12, name='measure2D1', coorddim=0, valdim=1, maxpoints=70)
cont = True
t0 = time.time()
while cont:
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
    # Check if the SNSPD is still superconducting
    if snspd.check() == False:
        print 'SNSPD went normal and could not restore!'
        scan_on = False
        break

    ph.ClearHistMem()
    ph.StartMeas(xsettings['AcqTime']*1000) # AcqTime in s, arg in ms
    print 'Acquiring signal for %s s' % (xsettings['AcqTime'])
    # Wait an extra 0.25 seconds
    time.sleep(xsettings['AcqTime']+0.05)


    current_data = ph.get_Histogram()

    bin_sum = np.sum(current_data[int(xsettings['bin_low']):int(xsettings['bin_high'])])



    print 'PH: %.3f' % bin_sum
    data.add_data_point(time.time()-t0, bin_sum)
    plot2d.update()
    time.sleep(0.01)

qt.mend()