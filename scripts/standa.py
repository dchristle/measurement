#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      DiamondAdmin
#
# Created:     18/12/2013
# Copyright:   (c) DiamondAdmin 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import ctypes
from ctypes import wintypes

class USMCDevices(ctypes.Structure):
    _fields_ = [
      ("NOD", wintypes.DWORD),
      ("Serial", ctypes.POINTER(ctypes.c_char_p)),
      ("Version", ctypes.POINTER(ctypes.c_char_p)),
    ]


def main():
    usmc = ctypes.cdll.USMCDLL
    init = usmc.USMC_Init
    init.restype = wintypes.DWORD
    init.argtypes = [ctypes.POINTER(USMCDevices)]
    dev = USMCDevices()
    init(dev)

    devices = [dev.Serial[i] + b':' + dev.Version[i]
              for i in range(dev.NOD)]
    print('\n'.join(d.decode('ascii') for d in devices))
    pass

if __name__ == '__main__':
    main()
