#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      DiamondAdmin
#
# Created:     17/11/2013
# Copyright:   (c) DiamondAdmin 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from PyDAQmx.DAQmxTypes import *
from PyDAQmx import *
import ctypes
import numpy
def main():
    read = int32()
    data = numpy.arange(-1,5.0,1, dtype=numpy.float64)
    analogOut = Task()
    sampleRate = 3
    analogOut.CreateAOVoltageChan("/DAQ6363/ao1",
    "",
    -10.0,
    10.0,
    DAQmx_Val_Volts,
    None)

    analogOut.CfgSampClkTiming(None,
    float64(sampleRate),
    DAQmx_Val_Rising,
    DAQmx_Val_FiniteSamps,
    data.size)

    analogOut.WriteAnalogF64(data.size,
    0,
    float64(-1),
    DAQmx_Val_GroupByChannel,
    data,
    None,
    None)

    ##analogOut.ExportSignal(DAQmx_Val_SampleClock,
    ##"/DAQ6363/PFI0")

    analogOut.StartTask()



##    taskHandle = TaskHandle(0)
    # Create task

##    PyDAQmx.DAQmxCreateTask("",byref(taskHandle))
##    PyDAQmx.DAQmxFunctions.DAQmxCreateAO
##    print '%s' % DAQmx_Val_Cfg_Default

##    print '%s' % data
##    try:
##        PyDAQmx.DAQmxCreateAOVoltageChan(taskHandle,"DAQ6363/ao0",
##        PyDAQmx.DAQmx_Val_Cfg_Default,-10.0,10.0,DAQmx_Val_Volts,None)
##    except PyDAQmx.DAQError as err:
##        print "DAQmx Error: %s"%err
##        DAQmxCreate
##
##    finally:
##        if taskHandle:
##           # DAQmx Stop Code
##          PyDAQmx.DAQmxStopTask(taskHandle)
##          PyDAQmx.DAQmxClearTask(taskHandle)

if __name__ == '__main__':
    main()
