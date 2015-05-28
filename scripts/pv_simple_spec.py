import numpy as np
import time

pv = qt.instruments['pv_here']

mc = qt.instruments['mc_here']


EXPOSURE_TIME = 1
GRATINGNO = 1
desired_wavelength = 1101.0;
desired_wlspeed = 500.0;

pv.set_exposure_time(EXPOSURE_TIME)


##cur_gno = mc.get_grating()
##if cur_gno != GRATINGNO:
##    mc.set_grating(GRATINGNO)
##mc.set_wavelength_speed(desired_wlspeed,timeout=30.0)
##mc.set_wavelength(desired_wavelength-50.0,timeout=30.0)
##time.sleep(1.0)
##mc.set_wavelength(desired_wavelength,timeout=30.0)

x_array = np.array(range(1024)) + 1

#w_conv = lambda x: 0.5241608 - -2.1974490e-6*(x-1000) - -3.4234693877e-8*(x-1000)*(x-1000)
#w_conv = lambda x: 0.12861 - 1.90909e-6*x-1.5051e-8*x*x



# Spectrometer constants

pwidth = 25.0e-3
center_pixel = 512.5
m = 1.0;



if GRATINGNO == 1:## Recalibrated on 04/25/2015
    d = 1.0/600.0
    gamma=29.3392/180.0*np.pi
    fl =300.6738
    delta = -1.93952/180.0*np.pi
    dpixel = 7.9292
    dl1 = 0.00066976
    dl2 = 0.0000049564
if GRATINGNO == 2:
    d = 1.0/150.0
    gamma=30.2992/180.0*np.pi
    fl = 298.6884
    delta = -2.8663/180.0*np.pi
    dpixel = 6.5353
    dl1 =  0.0075519577896
    dl2 = -0.0000309008512




# Begin wavelength conversion calculation
x_array = x_array - center_pixel + dpixel + dl1*(desired_wavelength-1100.0) + dl2*(desired_wavelength-1100.0)*(desired_wavelength-1100.0)

xi = np.arctan(x_array*pwidth*np.cos(delta)/(fl+x_array*pwidth*np.sin(delta)));

# Generate the psi and gamma values from the linear fit relations
psi = np.arcsin(m*desired_wavelength*1.0e-6/(2.0*d*np.cos(gamma/2.0)));

lambda_s = (d/m)*(np.sin(psi-gamma/2.0)+np.sin(psi+gamma/2.0+xi))*1.0e6;

spectrum_array = pv.acquire_image(timeout=(EXPOSURE_TIME+5.0));
wl_array = lambda_s
#qt.plot(wl_array,spectrum_array,clear=True)
qt.plot(wl_array,spectrum_array[0],clear=True)

