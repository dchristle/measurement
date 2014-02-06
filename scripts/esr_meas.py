# Take one: ESR measurement using NI-RFSG and Lockin drivers
# File name: ESRmeasurement_test.py
# F. J. Heremans <jhereman@gmail.com>
#
# If you did not already look at 'basic_measure_script.py'
# do that first.

# To load this module type 'import measure_module'
# to run an example function read on.

# Modules are often more convenient then scripts. In a
# module you can define man functions that can be resued
# in other functions or scripts. The functions are accessed
# by first importing the module, and then call a function
# within the module with <module>.<function>()
# Remember that a module has its own namespace, so many
# object that exist in the main namespace (like instruments)
# need to be imported explicitly.

import numpy as np
from time import time,sleep
import os
import qt
import logging
import msvcrt
from scipy.special import erf


def check_for_abort():
    if msvcrt.kbhit() and msvcrt.getch() == "q" :
        return False
    return True

######################################################
# example 1 - basic
######################################################
##def autophase(xa, ya):
##    # Routine will suggest an autophase for the lockin
##
##
##
##    suggestedPhase=0
##
##
##
##    sX = xa
##    sY = ya
##    sigma = errors;
##    A=max(abs(sX)+abs(sY))*2;
##    myllfunc = lambda theta: np.sum(np.log(exp(-(sY*np.cos(theta)-sX*np.sin(theta))**2/(2*sigma**2))*(erf((A+sX*np.cos(theta)+sY*np.sin(theta))/(np.sqrt(2)*sigma))-erf((-A+sX*np.cos(theta)+sY*np.sin(theta))/(sqrt(2)*sigma)))/(4*sqrt(2)*3.14159**1.5*sigma)))
##    myexp = lambda theta: np.exp(np.sum(np.log(np.exp(-(sY*np.cos(theta)-sX*np.sin(theta))**2/(2*sigma**2))*(erf((A+sX*np.cos(theta)+sY*np.sin(theta))/(sqrt(2)*sigma))-erf((-A+sX*np.cos(theta)+sY*np.sin(theta))/(sqrt(2)*sigma)))/(4*sqrt(2)*3.14159**1.5*sigma))))
##    myexptheta = lambda theta: theta*np.exp(np.sum(np.log(np.exp(-(sY*np.cos(theta)-sX*np.sin(theta))**2/(2*sigma**2))*(erf((A+sX*np.cos(theta)+sY*np.sin(theta))/(sqrt(2)*sigma))-erf((-A+sX*np.cos(theta)+sY*np.sin(theta))/(sqrt(2)*sigma)))/(4*sqrt(2)*3.14159**1.5*sigma))))
##
##    npoints = 300
##    thetavals = np.linspace(-3.141592,3.141592,npoints)
##    for ij,theta in enumerate(thetavals):
##        logvals[ij] = myllfunc(thetavals[ij])
##
##    C = np.amax(logvals)
##    I = np.argmax(logvals)
##
##    nDraws = 5000;
##    unifDraws = 2*(np.random.uniform(0,1,size=100)-0.5)*3.1415926535 + thetavals[I];
##    myllMax = C
##    ARfratio = np.zeros(nDraws)
##    for ij in range(nDraws):
##        ARfratio[ij] = np.exp(myllfunc(unifDraws[ij])-myllMax)
##
##
##    acceptDraws = np.random.uniform(0,1,size=nDraws)
##    ARidx = acceptDraws < ARfratio
##    ARaccepted = unifDraws(ARidx);
##
##    intFound = 0;
##    meanEst = mean(ARaccepted);
##    intLength = 0;
##    while intFound == 0
##    intLength = intLength + 0.001;
##    intARidx = (ARaccepted < meanEst + intLength & ARaccepted > meanEst - intLength);
##    if sum(intARidx)/length(ARaccepted) >= 0.95
##        intFound = 1;
##    end
##
##    end
##
##
##    if I > 1 && I < length(thetavals)
##    fdderiv = (logvals(I+1)-2*logvals(I)+logvals(I-1))/stepsize^2;
##    else
##    fdderiv = 1;
##    end
##
##    %suggestedPhase = thetavals(I)*180/pi;
##    %suggestedPhase = mean(RWval)*180/pi;
##    suggestedPhase = mean(ARaccepted)*180/pi;
##    [~, sPpi] = min([abs(suggestedPhase),abs(suggestedPhase+180),abs(suggestedPhase-180)]);
##
##    if sPpi == 2
##    suggestedPhase = suggestedPhase + 180;
##    end
##    if sPpi == 3
##    suggestedPhase = suggestedPhase - 180;
##    end
##
##    %montecarloError = std(RWval*180/pi);
##    %montecarloError = std(ARaccepted)*180/pi;
##    montecarloError = intLength*180/pi;
##
##    suggestedError = 1/sqrt(-1*fdderiv)*180/pi;
##    %suggestedError = countAccepted/9999;
def power_sweep(f_vec,p_vec,TC=11,reps=15):
    for i,p in enumerate(p_vec):
        simple(f_vec,p,TC,reps)
def simple(f_vec,power,TC = '100e-3s',reps = 15):
    '''
    this example is based on 'measure_module.py'
    you will need to set up a vector of frequencies and then call the
    measurment script.

    To run the function type in the terminal:

    fv = numpy.arange(1e9,2e9,50e6)
    esr_meas.simple(fv,power,TC)
    '''


    li = qt.instruments['lockin']
    pxi = qt.instruments['pxi']
    pxi.reset_device()
    qt.msleep(0.2)
    pxi.on()
    li.set_TC(TC)
    pxi.set_power(power)
    TCval = li.get_TC()

    qt.mstart()

    data = qt.Data(name='basic_esr')
    data.add_coordinate('Frequency, NI_RFSG [GHz]')
    data.add_value('Lockin X [mV]')
    data.add_value('Lockin Y [mv]')
    data.create_file()
    filename=data.get_filepath()[:-4]

    is_measuring = True
    qt.msleep(0.2)
    total_signal = np.zeros([np.size(f_vec),2])
    total_reps = 0
    for cur_rep in range(reps):
        for i,f in enumerate(f_vec):
            if (msvcrt.kbhit() and (msvcrt.getch() == 'q')):
                is_measuring = False
                break
            pxi.set_frequency(f)
            print 'frequency set: %s GHz' % (f*1e-9)
            qt.msleep(3*TCval+0.01)
            dx,dy = li.get_XY()
            total_signal[i,0] = total_signal[i,0] + dx
            total_signal[i,1] = total_signal[i,1] + dy
        total_reps = total_reps + 1


        #p_c = qt.Plot2D(f_vec, np.power(total_signal[:,0],2)+numpy.power(total_signal[:,1],2), 'bO-', name='ESR', clear=True)
        p_c = qt.Plot2D(f_vec, total_signal[:,0], 'bO-', name='ESR_X', clear=True)
        #p_c.add_data(f_vec, total_signal[:,1], name='ESR_Y', clear=True)


        if not is_measuring: break
    plot2d_1 = qt.Plot2D(data, name='measure2D', valdim=1)
    plot2d_1.add_data(data, valdim=2)
    data.add_data_point(f_vec,total_signal[:,0]/total_reps,total_signal[:,1]/total_reps)
    pxi.off()
    pxi.reset_device()
    plot2d_1.save_png(filename+'.png')
    data.close_file()
    qt.mend()



