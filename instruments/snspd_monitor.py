# SNSPD monitor virtual instrument
#
# The purpose of this instrument is to provide a background monitor for the
# SNSPD temperature, etc. during cooldowns. This instrument works by starting
# a monitor that is updated using the update method of this instrument. This
# update method is called once every update interval using the gobject
# timeout_add method to call the updating routine.

from instrument import Instrument
import types
import qt
import gobject
import pprint
import time
import hdf5_data as h5

class snspd_monitor(Instrument):

    def __init__(self, name):
        Instrument.__init__(self, name)
        self.name = name
        # Set up a few instruments as private subobjects
        self._snspd = qt.instruments['snspd']
        self._ls211 = qt.instruments['ls211']
        # Define time = 0 to be when the instrument is created
        self._t0 = 0
        # Set a parameter called is_running to be used to control the logic
        # of the update routine
        self.add_parameter('is_running',
            type=types.BooleanType,
            flags=Instrument.FLAG_GETSET)
        # Set the update interval here; can be changed
        self.add_parameter('update_interval',
            type=types.FloatType,
            unit='s',
            flags=Instrument.FLAG_GETSET)
        # Expose these functions to the user to start/stop/show the monitors
        self.add_function('start')
        self.add_function('stop')
        self.add_function('show_monitors')
        # Set default values of parameters
        self.set_is_running(False)
        self.set_update_interval(20)

    def do_get_is_running(self):
        return self._is_running

    def do_set_is_running(self, val):
        self._is_running = val

    def do_get_update_interval(self):
        return self._update_interval

    def do_set_update_interval(self, val):
        self._update_interval = val

    def start(self):
        # If running, do nothing. If not already running, start the monitor.
        if not self._is_running:
            # Set t0 to the start time.
            self._t0 = time.time()

            # Create a normal private qtlab data object called monitor_data.
            self._monitor_data = qt.Data(
            name='snspd_monitor_data')
            self._monitor_data.add_coordinate('time', format='%.2f')
            self._monitor_data.add_value('temperature', format='%.2f')

            # Create an HDF5 object for storing; we'll write to it upon stop
            self._monitor_dat = h5.HDF5Data(name='snspd_data')
            self._monitor_grp = h5.DataGroup('snspd_data_group', self._monitor_dat)

            # register some data dimensions
            self._monitor_grp.add_coordinate('time', unit='s', format='%.2f')
            self._monitor_grp.add_value('temperature', unit='K', format='%.2f')


            # Set the is_running variable to true, since we're now running
            self._is_running = True
            try:
                # Get the temperature, set variable _last_temperature to the current
                # temperature. This is a bit of a hack of the original adwin_monit0r
                # that just keeps track of this value so that when the user calls
                # show_monitors, we can print it.
                self._last_temperature = self._ls211.get_temperature()

            except:
                # If some sort of error occurs, print to the screen and cease
                # the gobject.timeout_add functionality by returning False
                print 'Could not get temperature from ls211, will stop now.'
                return False
            # Get temperature was successful
            # Add the just-measured temperature to the data object
            t = time.time() - self._t0
            self._monitor_data.add_data_point(t,
                    self._last_temperature)
            # Use gobject.timeout_add to run self._update every _update_interval
            gobject.timeout_add(int(self._update_interval*1e3), self._update)

    def stop(self):
        self._is_running = False
        data = self._monitor_data.get_data()
        self._monitor_grp['time'] = data[:,0]
        self._monitor_grp['temperature'] = data[:,1]
        self._monitor_dat.close()
    def show_plot(self):
        self._plot_on = True
        # Set up the plot object as a private subobject of the instrument
        self._monitor_plot = qt.Plot2D(
        self._monitor_data, name='snspd monitor', coorddim=0,
                    valdim=1, maxpoints=100, clear=True)
    def stop_plot(self):
        self._plot_on = False
        self._monitor_plot.clear()

    def show_monitors(self):
        print 'Temperature:'
        print '----'
        print ('%s' % self._last_temperature) + ' Kelvin'
        print





    def _update(self):
        # If not running, return False.
        # Returning False to gobject.timeout_add will cause it to stop
        # calling the update routine at all.
        if not self._is_running:
            return False
        # Calculate the current time
        t = time.time() - self._t0

        try:
            # Get the temperature, set variable _last_temperature to the current
            # temperature. This is a bit of a hack of the original adwin_monit0r
            # that just keeps track of this value so that when the user calls
            # show_monitors, we can print it.
            self._last_temperature = self._ls211.get_temperature()
            # Add the just-measured temperature to the data object
            self._monitor_data.add_data_point(t,
                    self._last_temperature)
        except:
            # If some sort of error occurs, print to the screen and cease
            # the gobject.timeout_add functionality by returning False
            print 'Could not get temperature from ls211, will stop now.'
            return False

        # If everything goes okay, return True, so that gobject.timeout_add
        # continutes to work.
        return True