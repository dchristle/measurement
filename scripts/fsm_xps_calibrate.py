import time
import msvcrt
import numpy as np


ph=qt.instruments['ph']
fbl = qt.instruments['fbl']
xps = qt.instruments['xps']
fsm = qt.instruments['fsm']
qt.mstart()

data0 = qt.Data(name='fsmcal')

# Now you provide the information of what data will be saved in the
# datafile. A distinction is made between 'coordinates', and 'values'.
# Coordinates are the parameters that you sweep, values are the
# parameters that you readout (the result of an experiment). This
# information is used later for plotting purposes.
# Adding coordinate and value info is optional, but recommended.
# If you don't supply it, the data class will guess your data format.
data0.add_coordinate('XPS displacement')
data0.add_value('X FSM displacement')
data0.add_value('Y FSM displacement')
data0.create_file()

# The next command will actually create the dirs and files, based
# on the information provided above. Additionally a settingsfile
# is created containing the current settings of all the instruments.

# Next two plot-objects are created. First argument is the data object
# that needs to be plotted. To prevent new windows from popping up each
# measurement a 'name' can be provided so that window can be reused.
# If the 'name' doesn't already exists, a new window with that name
# will be created. For 3d plots, a plotting style is set.
plot2d = qt.Plot2D(data0, name='fsmcal', coorddim=0, valdims=1)
plot2d.add_data(data0, coorddim=0, valdim=2)
cont = True
t0 = time.time()
disp_array = np.linspace(0,7,40)
x_orig = xps.get_abs_positionX()
y_orig = xps.get_abs_positionY()
fbl.optimize()
xf_orig = fsm.get_abs_positionX()
yf_orig = fsm.get_abs_positionY()
for i in range(15):
    if msvcrt.kbhit():
                kb_char=msvcrt.getch()
                if kb_char == "q" : break
    xps.set_abs_positionY(y_orig + 0.001*disp_array[i])
    time.sleep(0.1)
    fbl.optimize()
    fbl.optimize()
    xf_cur = fsm.get_abs_positionX()
    xf_diff = xf_cur-xf_orig
    yf_cur = fsm.get_abs_positionY()
    yf_diff = yf_cur - yf_orig
    data0.add_data_point(disp_array[i],xf_diff,yf_diff)


    plot2d.update()
xps.set_abs_positionY(y_orig)
fsm.set_abs_positionX(xf_orig)
fsm.set_abs_positionY(yf_orig)
data0.close_file()
qt.mend()