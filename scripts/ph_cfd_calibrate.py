import time
import msvcrt
import numpy as np


ph=qt.instruments['ph']

qt.mstart()
data0 = qt.Data(name='phcal')

# Now you provide the information of what data will be saved in the
# datafile. A distinction is made between 'coordinates', and 'values'.
# Coordinates are the parameters that you sweep, values are the
# parameters that you readout (the result of an experiment). This
# information is used later for plotting purposes.
# Adding coordinate and value info is optional, but recommended.
# If you don't supply it, the data class will guess your data format.
data0.add_coordinate('CFD')
data0.add_value('counts0')
data0.add_value('counts1')

# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.

# Next two plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.
plot2d = qt.Plot2D(data0, name='phcal', coorddim=0, valdims=1)
plot2d.add_data(data0, coorddim=0, valdim=2)
cont = True
t0 = time.time()
cfd_array = np.linspace(50,800,40)
for i in range(40):
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
    ph.set_InputCFD0(int(cfd_array[i]),10)
    ph.set_InputCFD1(int(cfd_array[i]),10)
    navgs = 20
    c0 = 0
    c1 = 0
    for j in range(navgs):
        time.sleep(0.1)
        c0 = c0 + ph.get_CountRate0()
        c1 = c1 + ph.get_CountRate1()
    c0 = float(c0)/float(navgs)
    c1 = float(c1)/float(navgs)
    time.sleep(0.5)
    data0.add_data_point(cfd_array[i],c0,c1)
    print 'CFD is %.2f, Ch0 counts: %.2f, Ch1 counts: %.2f' % (cfd_array[i], c0,c1)

    plot2d.update()


qt.mend()