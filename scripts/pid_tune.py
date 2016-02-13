import time
import msvcrt
import qt
ls332 = qt.instruments['ls332']

qt.mstart()
data = qt.Data(name='pid')


# Now you provide the information of what data will be saved in the
# datafile. A distinction is made between 'coordinates', and 'values'.
# Coordinates are the parameters that you sweep, values are the
# parameters that you readout (the result of an experiment). This
# information is used later for plotting purposes.
# Adding coordinate and value info is optional, but recommended.
# If you don't supply it, the data class will guess your data format.
data.add_coordinate('time')
data.add_value('temperature')
data.create_file()
# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.

# Next two plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.
plot2d = qt.Plot2D(data, name='measure2D', coorddim=0, valdim=1)
cont = True

ls332.set_cmode1(3)
ls332.set_mout1(20)
time.sleep(20)
t0 = time.time()
ls332.set_mout1(80)
while cont:
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
    data.add_data_point(time.time()-t0, ls332.get_kelvinA())
    plot2d.update()
    time.sleep(1.0)
    if (time.time()-t0) > 180.0:
        t1 = time.time()-t0
        cont = False
ls332.set_mout1(20)
t1 = time.time()
print 'Entering second step.'
cont = True
while cont:
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break

    data.add_data_point(time.time()-t0, ls332.get_kelvinA())
    plot2d.update()
    time.sleep(1.0)
    if (time.time()-t1) > 180.0:
        cont = False
data.close_file()
qt.mend()