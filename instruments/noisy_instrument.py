from instrument import Instrument
import types
import numpy as np

class noisy_instrument(Instrument):
    '''this is a dummy noisy instrument'''

    def __init__(self, name, reset=False):
        Instrument.__init__(self,name)

        self.add_parameter('angle1',
                type=types.FloatType,
                flags=Instrument.FLAG_GETSET | \
                Instrument.FLAG_GET_AFTER_SET,
                minval=-np.pi, maxval=np.pi,
                units='rad')
        self.add_parameter('angle2',
                type=types.FloatType,
                flags=Instrument.FLAG_GETSET | \
                Instrument.FLAG_GET_AFTER_SET,
                minval=0., maxval=np.pi,
                units='dBm')
        self.add_parameter('phase',
                type=types.FloatType,
                flags=Instrument.FLAG_GETSET | \
                Instrument.FLAG_GET_AFTER_SET,
                minval=-180., maxval=180.)
        self.add_parameter('cmax_position',
                type=types.FloatType,
                flags=Instrument.FLAG_GETSET | \
                Instrument.FLAG_GET_AFTER_SET,
                minval=1.0, maxval = 7.2)
        self.add_parameter('sigma',
                type=types.FloatType,
                flags=Instrument.FLAG_GETSET | \
                Instrument.FLAG_GET_AFTER_SET,
                minval=0, maxval = 10000)

        self.add_function('get_counts')
        self.add_function('set_optimum')
        self.add_function('set_initial')
        self.set_sigma(10)
        self.set_angle1(32.0/180.*np.pi)
        self.set_angle2(14.0/180*np.pi)
        self.set_phase(45.)
        self.set_cmax_position(5.2)




#### communication with machine

    def do_get_angle1(self):
        return self._angle1

    def do_set_angle1(self, angle1):
        self._angle1 = angle1

    def do_get_angle2(self):
        return self._angle2

    def do_set_angle2(self, angle2):
        self._angle2 = angle2

    def do_get_phase(self):
        return self._phase

    def do_set_phase(self, phase):
        self._phase = phase

    def do_get_cmax_position(self):
        return self._cmax_position

    def do_set_cmax_position(self,status):
        self._cmax_position = status

    def do_get_sigma(self):
        return self._sigma

    def do_set_sigma(self,status):
        self._sigma = status

### here we get a noisy
    def get_counts(self):
        # define the function
        f_ideal = (np.exp(-0.5*(self._angle1-42.8/180.0*np.pi)**2.0/((5.0/180.*np.pi)**2.0)) + \
            np.exp(-0.5*(self._angle2-29.8/180.0*np.pi)**2.0/((12.0/180.0*np.pi)**2.0)) + \
            np.exp(-0.5*(self._phase-2.0)**2.0/(5.0**2.0)))*1.0 - \
            15.0*(self._cmax_position-3.4289)**2.0
        f_inverted = f_ideal*(-1.0)
        noise = np.random.randn(1)[0]*self._sigma
        output = f_ideal + noise
        return output
    def set_optimum(self):
        self.set_angle1(42.8/180.*np.pi)
        self.set_angle2(29.8/180*np.pi)
        self.set_phase(2.0)
        self.set_cmax_position(3.4289)
        return
    def set_initial(self):
        self.set_angle1(32.0/180.*np.pi)
        self.set_angle2(14.0/180*np.pi)
        self.set_phase(45.)
        self.set_cmax_position(5.2)
        return



