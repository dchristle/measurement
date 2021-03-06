import time
import msvcrt

ph=qt.instruments['ph']

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
plot2d = qt.Plot2D(data, name='measure2D', coorddim=0, valdim=1, maxpoints=100)
cont = True
t0 = time.time()
while cont:
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
    print 'Ch0 counts: %s, Ch1 counts: %s' % (ph.get_CountRate0(), ph.get_CountRate1())
    data.add_data_point(time.time()-t0, ph.get_CountRate1())
    plot2d.update()
    time.sleep(0.25)

qt.mend()