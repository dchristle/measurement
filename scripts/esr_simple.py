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

def simple(f_vec):
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

    plot2d_1 = qt.Plot2D(data, name='measure2D', valdim=1)
    plot2d_1.add_data(data, valdim=2)
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

