import time
import qt
from qt import *
from instrument import Instrument




class Thorlabs_MFF001(Instrument):

    def __init__(self, name, line):

        Instrument.__init__(self, name)

        self.line = line

        self.ni63 = qt.instruments['NIDAQ6363']

        self.add_function('flip')




    def flip(self):



        l = self.line


        self.ni63.digital_out(l,True)
        time.sleep(1)
        self.ni63.digital_out(l,False)
    def measure_power(self):
        self.fm = qt.instruments['fm']
        self.fsm = qt.instruments['fsm']
        x0 = self.fsm.get_abs_positionX()
        y0 = self.fsm.get_abs_positionY()
        self.fsm.set_abs_positionX(0.0)
        self.fsm.set_abs_positionY(0.0)
        self.fm.set_numavgs(20)
        time.sleep(3)
        p_init = self.fm.get_power()
        self.flip()
        time.sleep(3)
        p_final = self.fm.get_power()
        self.flip()
        self.fsm.set_abs_positionX(x0)
        self.fsm.set_abs_positionY(y0)

        return (p_final - p_init)







