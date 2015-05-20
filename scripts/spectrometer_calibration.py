import numpy as np
import time
ws = qt.instruments['ws_here']

mc = qt.instruments['mc_here']
import scipy as scp


#w_conv = lambda x: 0.5241608 - -2.1974490e-6*(x-1000) - -3.4234693877e-8*(x-1000)*(x-1000)
#w_conv = lambda x: 0.12861 - 1.90909e-6*x-1.5051e-8*x*x

# Alignment wavelength list

wca = np.array((1013.975, 1128.740, 1357.021, 1367.351, 1395.055, 1529.582))

def lorentzian(x, gamma, x0, A, C):
    return A*(gamma*gamma/(gamma*gamma + (x-x0)*(x-x0))) + C
def gaussian(x, sigma, x0, A, C):
    return A*np.exp(-0.5*np.power((x-x0),2)/np.power(sigma,2)) + C
def get_wl_array(desired_wavelength, x_array, GRATINGNO=2):
    # Spectrometer constants

    pwidth = 25.0e-3
    center_pixel = 512.5
    m = 1.0;



    if GRATINGNO == 1: ## These are wrong after 04/25/2015 -- I had to recalibrate the spectrometer because it was torqued, and only calibrated grating #2
        d = 1.0/600.0
        gamma=29.3392/180.0*np.pi
        fl =300.6738
        delta = -1.93952/180.0*np.pi
        dpixel = 7.9292
        dl1 = 0.00066976
        dl2 = 0.0000049564
    if GRATINGNO == 2:
        d = 1.0/150.0
        gamma=27.8040/180.0*np.pi
        fl = 300.4224
        delta = -4.9729/180.0*np.pi
        dpixel = 6.0324
        dl1 =  -0.0015846424354
        dl2 = 0.0000132849551




    # Begin wavelength conversion calculation
    x_array = x_array - center_pixel + dpixel + dl1*(desired_wavelength-1100.0) + dl2*(desired_wavelength-1100.0)*(desired_wavelength-1100.0)

    xi = np.arctan(x_array*pwidth*np.cos(delta)/(fl+x_array*pwidth*np.sin(delta)));

    # Generate the psi and gamma values from the linear fit relations
    psi = np.arcsin(m*desired_wavelength*1.0e-6/(2.0*d*np.cos(gamma/2.0)));

    lambda_s = (d/m)*(np.sin(psi-gamma/2.0)+np.sin(psi+gamma/2.0+xi))*1.0e6;
    return lambda_s




EXPOSURE_TIME = 15
GRATINGNO = 1

desired_wlspeed = 500.0;

ws.set_exposure_time(EXPOSURE_TIME)


cur_gno = mc.get_grating()
if cur_gno != GRATINGNO:
    mc.set_grating(GRATINGNO)
mc.set_wavelength_speed(desired_wlspeed,timeout=10)


x_array = np.array(range(1024)) + 1

mc_steps = np.array((1005.0, 1007.0, 1020.0, 1030.0, 1105.0, 1110.0, 1120.0, 1130.0, 1345.0, 1360.0, 1370.0, 1390.0, 1400.0, 1515, 1520, 1535.0, 1545.0))
cal_array = []
bw = 3.0 # nm
for cur_wl in mc_steps:

    print 'Setting 50 nm before wl.'
    mc.set_wavelength(cur_wl-50.0,timeout=30)

    time.sleep(1.0)
    print 'Setting wl.'
    mc.set_wavelength(cur_wl,timeout=30)
    print 'wl set, getting spectrum'
    spectrum_array = np.array(ws.get_spectrum([],timeout=(EXPOSURE_TIME+6.0))[0,:,0]);
    wl_array = get_wl_array(cur_wl, x_array, GRATINGNO)
    print wl_array
    print 'got spectrum'
    for cal_wl in wca:
        distance_array = np.abs(wl_array-cal_wl)
        print distance_array
        if any(distance_array < bw):
            print 'pixels within area for %s' % cal_wl
            # check the predicted pixel is far enough away from ends
            idx = np.argmin(distance_array)
            print 'idx is %s' % idx
            if idx > 20 and idx < 1024-20:
                print 'idx in bounds'
                # it's far enough away, so use it for calibration
                center_wl = wl_array[idx]
                print 'center wl is %s' % center_wl
                idx_select = np.logical_and((wl_array < (center_wl + bw)),(wl_array > (center_wl - bw)))
                print 'selected idx %s' % idx_select
                x_vals = x_array[idx_select]
                y_vals = spectrum_array[idx_select]
                x_guess = np.array((2, np.mean(x_vals), np.max(y_vals)-np.min(y_vals), np.min(y_vals)))
                print 'doing optimize...'
                scpout =  scp.optimize.curve_fit(gaussian, x_vals, y_vals, x_guess, maxfev = 10000)
                print scpout
                #qt.plot(x_vals,y_vals)
                xvp = np.linspace(x_vals[0],x_vals[-1],100)
                #qt.plot(xvp,gaussian(xvp, scpout[0][0],scpout[0][1],scpout[0][2],scpout[0][3]))
                cal_array.append((cur_wl, scpout[0][1], cal_wl))
                print 'Appended to cal array: %s %s %s' % (cur_wl, scpout[0][1], cal_wl)
#qt.plot(wl_array,spectrum_array,clear=True)


##
##data_array = np.array([x_array+center_pixel, spectrum_array])
##name_for_data = ('spectrum at %.2f' % desired_wavelength)
##data = qt.Data(name=name_for_data)
##data.add_coordinate('wavelength (nm)')
##data.add_value('counts')
##data.create_file()
##for ij in range(wl_array.size):
##    data.add_data_point(wl_array[ij],spectrum_array[ij])
##    #data.add_data_point(x_array[ij],spectrum_array[ij])
##
##
##
##
##data.close_file()
