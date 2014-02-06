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

import numpy
from time import time,sleep
import os
import qt


def instr_setup():
    # The purpose of this is to setup the instruments, so that we can just
    # reference them later, instead of constantly creating and destroying
    # the objects.

    # define instruments
    Lockin_address = 'GPIB::17'
    NI_RFSG_address = 'IQ5611'
    qt.instruments.create('pxi','NI_RFSG',resource_name=NI_RFSG_address)
    qt.instruments.create('lockin','Lockin_726x',address = Lockin_address)

######################################################
# example 1 - basic
######################################################

def example1(f_vec):
    '''
    this example is based on 'measure_module.py'
    you will need to set up a vector of frequencies and then call the
    measurment script.

    To run the function type in the terminal:

    fv = numpy.arange(1e9,2e9,50e6)
    ESRmeasurement_test.example1(fv)
    '''


    li = qt.instruments['lockin']
    pxi = qt.instruments['pxi']


    qt.mstart()

    data = qt.Data(name='basic_esr')
    data.add_coordinate('Frequency, NI_RFSG [GHz]')
    data.add_value('Lockin X [mV]')
    data.add_value('Lockin Y [mv]')
    data.create_file()

    plot2d_1 = qt.Plot2D(data, name='measure2D',valdim=1)
    plot2d_2 = qt.Plot2D(data, name='measure2D',valdim=2)
##    plot2d_2 = qt.Plot2D(data, name='measure2D',valdim=2)
##    plot3d = qt.Plot3D(data, name='measure3D', style='image')



    for f in f_vec:
        pxi.set_frequency(f)
        print 'frequency set: %s GHz' % (f*1e-9)

        dx,dy = li.get_XY()
        data.add_data_point(f,dx,dy)

        qt.msleep(3*li.get_TC()+0.01)
    data.new_block()

    data.close_file()
    qt.mend()


#######################
# example 2 - data
#######################

def example2(f_vec, b_vec):
    '''
    This example introduces three new features:
    1) setting format and/or precision of data in the datafile.
       using 'precision' will keep the default scientific notation,
       'format' can be anything you like
        => add_coordinate(precision=<nr>)
        => add_coordinate(format='<format_string>')
    2) specify specific filepath for the data file (in stead of
       automatic filepath)
        => create_file(filepath=<filepath>)
    3) turn off automatic saving of instrument-settings-file.
        => create_file(settings_file=False)

    To run the function type in the terminal:

    fv=numpy.arange(0,10,0.01)
    bv=numpy.arange(-5,5,0.1)
    measure_module.example2(fv, bv)
    '''

    qt.mstart()

    # this shows how to change format of saved data (per column)
    data = qt.Data(name='testmeasurement')
    data.add_coordinate('frequency, mw src 1 [Hz]', precision=3)
    data.add_coordinate('Bfield, ivvi dac 3 [mV]', format='%.12f')
    data.add_value('Psw SQUID', format='%.3e')
    data.create_file()

    # this shows how to save to a specific path and name, and how
    # to avoid a settings file to be created. The directory is first
    # retreived from the previous data object
    dir = data.get_dir()
    maxfilepath = os.path.join(dir, 'maxvals.dat')

    data_max = qt.Data(name='maxvals')
    data_max.add_coordinate('Bfield, ivvi dac 3 [mV]')
    data_max.add_value('resonance frequency [Hz]')
    data_max.create_file(
            filepath=maxfilepath,
            settings_file=False)

    plot2d = qt.Plot2D(data, name='measure2D')
    plot3d = qt.Plot3D(data, name='measure3D', style='image')

    plot2dmax = qt.Plot2D(data_max, name='maxvals')

    for b in b_vec:
        fake_ivvi_set_dac_3(b)

        last_trace = []
        for f in f_vec:
            fake_mw_src_set_freq(f)

            result = fake_readout_psw()
            data.add_data_point(f, b, result)

            last_trace.append(result)

            qt.msleep(0.01)
        data.new_block()

        loc_of_max = numpy.argmax(last_trace)
        freq_at_max = f_vec[loc_of_max]

        data_max.add_data_point(b, freq_at_max)

    data.close_file()
    data_max.close_file()
    qt.mend()


#######################
# example 3 - plotting
#######################

def example3(x_vec=numpy.linspace(0,10,10), y_vec=numpy.linspace(0,10,50)):
    '''
    To run the function type in the terminal:

    measure_module.example3()
    '''

    qt.mstart()

    data = qt.Data(name='testmeasurement')
    data.add_coordinate('x')
    data.add_coordinate('y')
    data.add_value('z1')
    data.add_value('z2')
    data.add_value('z3')
    data.create_file()

    plot2d_1 = qt.Plot2D(data, name='2D_1', coorddim=1, valdim=2)

    plot2d_2 = qt.Plot2D(data, name='2D_2', coorddim=1, valdim=2, maxtraces=1)

    plot2d_3 = qt.Plot2D(data, name='2D_3', coorddim=1, valdim=2, maxtraces=1)
    plot2d_3.add_data(data, coorddim=1, valdim=3, maxtraces=1)
    plot2d_3.add_data(data, coorddim=1, valdim=4, maxtraces=1)

    plot2d_4 = qt.Plot2D(data, name='2D_4', coorddim=1, valdim=2, mintime=0.3)

    plot2d_5 = qt.Plot2D(data, name='2D_5', coorddim=1, valdim=2, autoupdate=False)

    plot3d_1 = qt.Plot3D(data, name='3D_1', style='image')

    plot3d_2 = qt.Plot3D(data, name='3D_2', style='image', coorddims=(1,0), valdim=4)

    for x in x_vec:
        for y in y_vec:

            z1 = numpy.sin(x+y)
            z2 = numpy.cos(x+y)
            z3 = numpy.sin(x+2*y)

            data.add_data_point(x, y, z1, z2, z3)

            if z1>0:
                plot2d_5.update()

            qt.msleep(0.1)
        data.new_block()

    plot2d_1.save_png()
    plot2d_1.save_gp()

    plot3d_2.save_png()
    plot3d_2.save_gp()

    data.close_file()
    qt.mend()

