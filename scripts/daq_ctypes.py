#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      DiamondAdmin
#
# Created:     19/11/2013
# Copyright:   (c) DiamondAdmin 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# import system libraries

import ctypes
import numpy
import threading
# load any DLLs
nidaq = ctypes.windll.LoadLibrary("nicaiu.dll") # load the DLL
def CHK(err ):
    """a simple error checking routine"""
    if err < 0:
        buf_size = 100
        buf = ctypes.create_string_buffer('\000' * buf_size)
        nidaq.DAQmxGetErrorString(err,ctypes.byref(buf),buf_size)
        raise RuntimeError('nidaq call failed with error %d: %s'%(err,repr(buf.value)))
    if err > 0:
        buf_size = 100
        buf = ctypes.create_string_buffer('\000' * buf_size)
        nidaq.DAQmxGetErrorString(err,ctypes.byref(buf),buf_size)
        raise RuntimeError('nidaq generated warning %d: %s'%(err,repr(buf.value)))

def run( taskHandle ):
    counter = 0
    CHK(nidaq.DAQmxStartTask( taskHandle ))

def stop(taskHandle):
    nidaq.DAQmxStopTask( taskHandle )
    nidaq.DAQmxClearTask( taskHandle )

def main():

    ##############################
    # Setup some typedefs and constants
    # to correspond with values in
    # C:\Program Files\National Instruments\NI-DAQ\DAQmx ANSI C Dev\include\NIDAQmx.h
    # the typedefs
    int32 = ctypes.c_long
    uInt32 = ctypes.c_ulong
    uInt64 = ctypes.c_ulonglong
    float64 = ctypes.c_double
    TaskHandle = uInt32
    # the constants
    DAQmx_Val_Cfg_Default = int32(-1)
    DAQmx_Val_Volts = 10348
    DAQmx_Val_Rising = 10280
    DAQmx_Val_FiniteSamps = 10178
    DAQmx_Val_ContSamps = 10123
    DAQmx_Val_GroupByChannel = 0
    DAQmx_Val_SampleClock = 12487
    ##############################


    running = True
    sampleRate = 2
    # thanks to Lenard Lindstrom for this one
    DOUBLEPTR = ctypes.POINTER(ctypes.c_double)
    taskHandle = TaskHandle( 0 )
    data = numpy.arange(-1,1,0.05, dtype=numpy.float64)
    # setup the DAQ hardware
    pointsRead = uInt32()
    CHK(nidaq.DAQmxCreateTask("",ctypes.byref(taskHandle)))
    CHK(nidaq.DAQmxCreateAOVoltageChan(taskHandle,
                                   "DAQ6363/ao1",
                                   "",
                                   float64(-10.0),
                                   float64(10.0),
                                   DAQmx_Val_Volts,
                                   None))
    CHK(nidaq.DAQmxCfgSampClkTiming(taskHandle,
                                "",
                                float64(sampleRate),
                                DAQmx_Val_Rising,
                                DAQmx_Val_FiniteSamps,
                                uInt64(data.size)));
    CHK(nidaq.DAQmxWriteAnalogF64(taskHandle,
                              int32(data.size),
                              False,
                              float64(-1),
                              DAQmx_Val_GroupByChannel,
                              ctypes.cast(data.ctypes.data, DOUBLEPTR),
                              ctypes.byref(pointsRead),
                              None))
    CHK(nidaq.DAQmxExportSignal(taskHandle, DAQmx_Val_SampleClock,
    "/DAQ6363/PFI0"))
    CHK(nidaq.DAQmxStartTask (taskHandle))

if __name__ == '__main__':
    main()
