import numpy as np
import logging
import qt
import hdf5_data as h5
import measurement.lib.measurement2.measurement as m2
import time
import msvcrt
from measurement.lib.pulsar import pulse, pulselib, element, pulsar
from random import shuffle
reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

class SiC_Spectrum_Master(m2.Measurement):

    mprefix = 'spectrum'

    def sequence(self, upload = True, program=True):
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

        elements = []
        # Create waveform that has laser, microwaves, photon counting, and 1/0 I/Q modulation on
        # all the time for a long period of time (~100 us).
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='laser')
        e.add(pulse.cp(sq_pulseMW, amplitude=1.0, length=100e-6), name='microwaves')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)

        # Create a sequence that lacks the AOM pulse -- we want to turn the AOM off for background subtraction
        e = element.Element('dark_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=0, length=100e-6), name='laseroff')
        e.add(pulse.cp(sq_pulseMW, amplitude=1.0, length=100e-6), name='microwaves')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)



        seq = pulsar.Sequence('Photoluminescence Sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(3.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.1)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.2)

        if q >= 20:
            print 'AWG did not jump to proper waveform!'
        return
    def prepare(self):
        self.start_keystroke_monitor('abort')
        # Set up some instruments
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']
        self._ws = qt.instruments['ws_here'] # winspec -- this is a remote instrument
        self._mc = qt.instruments['mc'] # monochromator -- this is a remote instrument

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        self._awg.start()
        time.sleep(5.0)
        self._awg.sq_forced_jump(1)
        time.sleep(3.0)
        self.awg_confirm(1)


        self._fbl.optimize()
        print 'FBL optimized...'
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))



        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))




        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counter prepared.'
        # Configure the exposure time in WinSpec
        self._ws.set_exposure_time(self.params['exposure_time'])
        # Set the PXI system
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()
        # Now set the new parameters
        self._pxi.set_power(self.params['power'])
        self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
        # Now set the proper attenuation
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten
        if self.params['microwaves']:
            self._pxi.set_status('on')
        else:
            self._pxi.set_status('off')
        time.sleep(2.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return
        self._fbl.optimize()
        time.sleep(5.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return

        if not self.params['microwaves']:
            self._fbl.optimize()
            time.sleep(10.0)
            self._fbl.optimize()





        return
    def wavelength_array(self,desired_wavelength,GRATINGNO):
        x_array = np.array(range(1024)) + 1
        # Spectrometer constants

        pwidth = 25.0e-3
        center_pixel = 512.5
        m = 1.0;



        if GRATINGNO == 1:
            d = 1.0/600.0
            gamma=30.5852/180.0*np.pi
            fl = 297.5702
            delta = -4.5612/180.0*np.pi
            dpixel = -0.6415
            dl1 = -0.0009530671345
            dl2 = 0.0000041768437
        if GRATINGNO == 2:
            d = 1.0/150.0
            gamma=32.0063/180.0*np.pi
            fl = 296.2841
            delta = -5.2596/180.0*np.pi
            dpixel = -1.6232
            dl1 = -0.0008462231572
            dl2 =  0.0000047503028




        # Begin wavelength conversion calculation
        x_array = x_array - center_pixel + dpixel + dl1*(desired_wavelength-1000.0) + dl2*(desired_wavelength-1000.0)*(desired_wavelength-1000.0)

        xi = np.arctan(x_array*pwidth*np.cos(delta)/(fl+x_array*pwidth*np.sin(delta)));

        # Generate the psi and gamma values from the linear fit relations
        psi = np.arcsin(m*desired_wavelength*1.0e-6/(2.0*d*np.cos(gamma/2.0)));

        lambda_s = (d/m)*(np.sin(psi-gamma/2.0)+np.sin(psi+gamma/2.0+xi))*1.0e6;
        return lambda_s

    def measure(self):
        print 'Photoluminescence spectra measurement sequence initiated.'

        self._stop_measurement = False
        # Take note of the time the measurement starts.
        t0 = time.time()
        lo_spec_total = np.zeros(1024,dtype=np.uint32)
        lr_spec_total = np.zeros(1024,dtype=np.uint32)
        no_spec_total = np.zeros(1024,dtype=np.uint32)
        bg_spec_total = np.zeros(1024,dtype=np.uint32)

        wl_array = self.wavelength_array(1150,2)
        for k in range(int(self.params['MeasCycles'])):
            print 'Starting iteration %d of %d' % ((k+1),self.params['MeasCycles'])
            # Take the spectra in random order
            seq_idx = range(4)
            if self.params['random'] == 1:
                shuffle(seq_idx)

            for j in range(4):
                if seq_idx[j] == 0 and self.params['microwaves']:
                    # laser & off-resonant microwaves
                    print 'PL iteration %d/%d, step %d/%d: Laser with off-resonant microwaves, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    lo_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                    if not type(lo_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        lo_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))




                if seq_idx[j] == 1:
                    # laser & resonant microwaves
                    print 'PL iteration %d/%d, step %d/%d: Laser with resonant microwaves, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    self._pxi.set_frequency((self.params['frequency'])*1.0e9)
                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    lr_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+5.0))
                    if not type(lr_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        lr_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))


                if seq_idx[j] == 2:
                    # no laser & off-resonant microwaves
                    print 'PL iteration %d/%d, step %d/%d: No laser with off-resonant microwaves, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    self._awg.sq_forced_jump(2)
                    time.sleep(0.2)
                    self.awg_confirm(2)
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    no_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                    if not type(no_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        no_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))


                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                if seq_idx[j] == 3:
                    print 'PL iteration %d/%d, step %d/%d: Laser and off-resonant microwaves, displaced to a background spot, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    # laser & off-resonant microwaves, and displaced from the defect -- this is the background signal
                    self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    cur_X = self._fsm.get_abs_positionX()
                    cur_Y = self._fsm.get_abs_positionY()
                    self._fsm.set_abs_positionX(cur_X + self.params['x_displacement'])
                    self._fsm.set_abs_positionY(cur_Y + self.params['y_displacement'])
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    bg_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+5.0))
                    self._fsm.set_abs_positionX(cur_X)
                    self._fsm.set_abs_positionY(cur_Y)
                    if not type(bg_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        self._fsm.set_abs_positionX(cur_X + self.params['x_displacement'])
                        self._fsm.set_abs_positionY(cur_Y + self.params['y_displacement'])
                        bg_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))



            	if not self._stop_measurement:
            	    self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q'] or self._stop_measurement:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._pxi.set_status('off')
                    self._stop_measurement = True
                    print 'Breaking from internal measurement loop...'
                    break
            detector_temperature = self._ws.get_temperature()
            target_temperature = self._ws.get_target_temperature()
            if np.abs(float(detector_temperature) - float(target_temperature)) > self.params['temperature_tolerance']:
                print 'Spectrometer camera warming up -- breaking from main measurement loop...'
                break
            # Now start checking for other issues. If present, stop.
            # Check if setpoint and actual temperature are within a tolerance
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
            # Check if the SNSPD is still superconducting
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break

            if self._stop_measurement:
                print 'Breaking from main measurement loop...'
                break
            if not ((type(lr_spec) == np.ndarray) and (type(no_spec) == np.ndarray) and (type(bg_spec) == np.ndarray)):
                break
            if self.params['microwaves']:
                if not (type(lo_spec) == np.ndarray):
                    break
            if self.params['microwaves']:
                lo_spec_total = lo_spec_total + lo_spec[0].T[0].astype(np.uint32)
            lr_spec_total = lr_spec_total + lr_spec[0].T[0].astype(np.uint32)
            no_spec_total = no_spec_total + no_spec[0].T[0].astype(np.uint32)
            bg_spec_total = bg_spec_total + bg_spec[0].T[0].astype(np.uint32)
            if self.params['microwaves']:
                lo_plot_array = lo_spec_total.astype(np.float64)-no_spec_total.astype(np.float64)
            lr_plot_array = lr_spec_total.astype(np.float64)-no_spec_total.astype(np.float64)
            bg_plot_array = bg_spec_total.astype(np.float64)-no_spec_total.astype(np.float64)
            if self.params['microwaves']:
                qt.plot(wl_array,lo_plot_array,name='lo_spec',clear=True)
            qt.plot(wl_array,lr_plot_array,name='lr_spec',clear=True)
            qt.plot(wl_array,no_spec_total,name='no_spec',clear=True)
            qt.plot(wl_array,bg_plot_array,name='bg_spec',clear=True)
            qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.

        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Photoluminescence_data', self.h5data, base=self.h5base)
        if self.params['microwaves']:
            grp.add('lo_spec_total', data=lo_spec_total, unit='counts', note='lo total count array')
        grp.add('lr_spec_total', data=lr_spec_total, unit='counts', note='lr total count array')
        grp.add('no_spec_total', data=no_spec_total, unit='counts', note='no total count array')
        grp.add('bg_spec_total', data=bg_spec_total, unit='counts', note='bg total count array')
        return

class SiC_PLvsT_Master(m2.Measurement):

    mprefix = 'PLvsT'

    def sequence(self, upload = True, program=True):
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

        elements = []
        # Create waveform that has laser, microwaves, photon counting, and 1/0 I/Q modulation on
        # all the time for a long period of time (~100 us).
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='laser')
        e.add(pulse.cp(sq_pulseMW, amplitude=1.0, length=100e-6), name='microwaves')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)

        # Create a sequence that lacks the AOM pulse -- we want to turn the AOM off for background subtraction
        e = element.Element('dark_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=0, length=100e-6), name='laseroff')
        e.add(pulse.cp(sq_pulseMW, amplitude=1.0, length=100e-6), name='microwaves')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)



        seq = pulsar.Sequence('Photoluminescence Sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(3.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.1)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.2)

        if q >= 20:
            print 'AWG did not jump to proper waveform!'
        return
    def prepare(self):
        self.start_keystroke_monitor('abort')
        # Set up some instruments
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']
        self._ws = qt.instruments['ws_here'] # winspec -- this is a remote instrument
        self._mc = qt.instruments['mc'] # monochromator -- this is a remote instrument

        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        self._awg.start()
        time.sleep(5.0)
        self._awg.sq_forced_jump(1)
        time.sleep(3.0)
        self.awg_confirm(1)


        self._fbl.optimize()
        print 'FBL optimized...'
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))



        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))




        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counter prepared.'
        # Configure the exposure time in WinSpec
        self._ws.set_exposure_time(self.params['exposure_time'])
        # Set the PXI system
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()
        # Now set the new parameters
        self._pxi.set_power(self.params['power'])
        self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
        # Now set the proper attenuation
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation.' % desired_atten
        if self.params['microwaves']:
            self._pxi.set_status('on')
        else:
            self._pxi.set_status('off')
        time.sleep(2.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return
        self._fbl.optimize()
        time.sleep(5.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return

        if not self.params['microwaves']:
            self._fbl.optimize()
            time.sleep(10.0)
            self._fbl.optimize()





        return
    def wavelength_array(self,desired_wavelength,GRATINGNO):
        x_array = np.array(range(1024)) + 1
        # Spectrometer constants

        pwidth = 25.0e-3
        center_pixel = 512.5
        m = 1.0;



        if GRATINGNO == 1:
            d = 1.0/600.0
            gamma=30.5852/180.0*np.pi
            fl = 297.5702
            delta = -4.5612/180.0*np.pi
            dpixel = -0.6415
            dl1 = -0.0009530671345
            dl2 = 0.0000041768437
        if GRATINGNO == 2:
            d = 1.0/150.0
            gamma=32.0063/180.0*np.pi
            fl = 296.2841
            delta = -5.2596/180.0*np.pi
            dpixel = -1.6232
            dl1 = -0.0008462231572
            dl2 =  0.0000047503028




        # Begin wavelength conversion calculation
        x_array = x_array - center_pixel + dpixel + dl1*(desired_wavelength-1000.0) + dl2*(desired_wavelength-1000.0)*(desired_wavelength-1000.0)

        xi = np.arctan(x_array*pwidth*np.cos(delta)/(fl+x_array*pwidth*np.sin(delta)));

        # Generate the psi and gamma values from the linear fit relations
        psi = np.arcsin(m*desired_wavelength*1.0e-6/(2.0*d*np.cos(gamma/2.0)));

        lambda_s = (d/m)*(np.sin(psi-gamma/2.0)+np.sin(psi+gamma/2.0+xi))*1.0e6;
        return lambda_s

    def measure(self):
        print 'Photoluminescence spectra measurement sequence initiated.'

        self._stop_measurement = False
        # Take note of the time the measurement starts.
        t0 = time.time()
        lo_spec_total = np.zeros(1024,dtype=np.uint32)
        lr_spec_total = np.zeros(1024,dtype=np.uint32)
        no_spec_total = np.zeros(1024,dtype=np.uint32)
        bg_spec_total = np.zeros(1024,dtype=np.uint32)

        wl_array = self.wavelength_array(1150,2)
        for k in range(int(self.params['MeasCycles'])):
            print 'Starting iteration %d of %d' % ((k+1),self.params['MeasCycles'])
            # Take the spectra in random order
            seq_idx = range(4)
            if self.params['random'] == 1:
                shuffle(seq_idx)

            for j in range(4):
                if seq_idx[j] == 0 and self.params['microwaves']:
                    # laser & off-resonant microwaves
                    print 'PL iteration %d/%d, step %d/%d: Laser with off-resonant microwaves, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    lo_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                    if not type(lo_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        lo_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))




                if seq_idx[j] == 1:
                    # laser & resonant microwaves
                    print 'PL iteration %d/%d, step %d/%d: Laser with resonant microwaves, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    self._pxi.set_frequency((self.params['frequency'])*1.0e9)
                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    lr_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+5.0))
                    if not type(lr_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        lr_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))


                if seq_idx[j] == 2:
                    # no laser & off-resonant microwaves
                    print 'PL iteration %d/%d, step %d/%d: No laser with off-resonant microwaves, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    self._awg.sq_forced_jump(2)
                    time.sleep(0.2)
                    self.awg_confirm(2)
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    no_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                    if not type(no_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        no_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))


                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                if seq_idx[j] == 3:
                    print 'PL iteration %d/%d, step %d/%d: Laser and off-resonant microwaves, displaced to a background spot, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 4, self.params['exposure_time'])
                    # laser & off-resonant microwaves, and displaced from the defect -- this is the background signal
                    self._pxi.set_frequency((self.params['frequency']+self.params['detuning'])*1.0e9)
                    self._awg.sq_forced_jump(1)
                    time.sleep(0.2)
                    self.awg_confirm(1)
                    self._fbl.optimize()
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    cur_X = self._fsm.get_abs_positionX()
                    cur_Y = self._fsm.get_abs_positionY()
                    self._fsm.set_abs_positionX(cur_X + self.params['x_displacement'])
                    self._fsm.set_abs_positionY(cur_Y + self.params['y_displacement'])
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    bg_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+5.0))
                    self._fsm.set_abs_positionX(cur_X)
                    self._fsm.set_abs_positionY(cur_Y)
                    if not type(bg_spec) == np.ndarray:
                        time.sleep(10.0)
                        self._fbl.optimize()
                        self._fsm.set_abs_positionX(cur_X + self.params['x_displacement'])
                        self._fsm.set_abs_positionY(cur_Y + self.params['y_displacement'])
                        bg_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))



            	if not self._stop_measurement:
            	    self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q'] or self._stop_measurement:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._pxi.set_status('off')
                    self._stop_measurement = True
                    print 'Breaking from internal measurement loop...'
                    break
            detector_temperature = self._ws.get_temperature()
            target_temperature = self._ws.get_target_temperature()
            if np.abs(float(detector_temperature) - float(target_temperature)) > self.params['temperature_tolerance']:
                print 'Spectrometer camera warming up -- breaking from main measurement loop...'
                break
            # Now start checking for other issues. If present, stop.
            # Check if setpoint and actual temperature are within a tolerance
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
            # Check if the SNSPD is still superconducting
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break

            if self._stop_measurement:
                print 'Breaking from main measurement loop...'
                break
            if not ((type(lr_spec) == np.ndarray) and (type(no_spec) == np.ndarray) and (type(bg_spec) == np.ndarray)):
                break
            if self.params['microwaves']:
                if not (type(lo_spec) == np.ndarray):
                    break
            if self.params['microwaves']:
                lo_spec_total = lo_spec_total + lo_spec[0].T[0].astype(np.uint32)
            lr_spec_total = lr_spec_total + lr_spec[0].T[0].astype(np.uint32)
            no_spec_total = no_spec_total + no_spec[0].T[0].astype(np.uint32)
            bg_spec_total = bg_spec_total + bg_spec[0].T[0].astype(np.uint32)
            if self.params['microwaves']:
                lo_plot_array = lo_spec_total.astype(np.float64)-no_spec_total.astype(np.float64)
            lr_plot_array = lr_spec_total.astype(np.float64)-no_spec_total.astype(np.float64)
            bg_plot_array = bg_spec_total.astype(np.float64)-no_spec_total.astype(np.float64)
            if self.params['microwaves']:
                qt.plot(wl_array,lo_plot_array,name='lo_spec',clear=True)
            qt.plot(wl_array,lr_plot_array,name='lr_spec',clear=True)
            qt.plot(wl_array,no_spec_total,name='no_spec',clear=True)
            qt.plot(wl_array,bg_plot_array,name='bg_spec',clear=True)
            qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.

        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Photoluminescence_data', self.h5data, base=self.h5base)
        if self.params['microwaves']:
            grp.add('lo_spec_total', data=lo_spec_total, unit='counts', note='lo total count array')
        grp.add('lr_spec_total', data=lr_spec_total, unit='counts', note='lr total count array')
        grp.add('no_spec_total', data=no_spec_total, unit='counts', note='no total count array')
        grp.add('bg_spec_total', data=bg_spec_total, unit='counts', note='bg total count array')
        return


class SiC_DoublePulse_Master(m2.Measurement):

    mprefix = 'doublepulse'

    def sequence(self, upload = True, program=True, clear=False):
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulsePH = pulse.SquarePulse(channel='phtrigger', name='A square pulse on phtrigger')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

        self.params['pts'] = np.uint32(1 + np.ceil(np.abs(self.params['tau_length_end'] - self.params['tau_length_start'])/self.params['tau_length_step']))
        self.params['tau_delay'] = np.linspace(self.params['tau_length_start'], self.params['tau_length_end'], self.params['pts'])
        self._awg = qt.instruments['awg']
        print 'Stopping AWG.'
        self._awg.stop()
        time.sleep(3)
        if clear:
            self._awg.clear_waveforms()
            print 'AWG waveforms cleared.'
        elements = []

        # First create a waveform that keeps the AOM and photon counting switch on all the time
        # but leaves the microwave switch off
        e = element.Element('CW_mode_wfm', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1.0, length=100e-6,start=0.0e-9), name='lasercw')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6,start=0.0e-9), name='photoncountpulsecw')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        # The following pulse is for testing purposes only and should be commented out.
##        e.add(pulse.cp(sq_pulsePH, amplitude=-0.7, length=50e-9), name='picoharp trigger', start=self.params['PH_trigger_time']*1.0e-9)
        elements.append(e)


        # Now create the Double Pulse pulses
        for i in range(self.params['pts']):

            for j in range(2):
                if j == 0:
                    e = element.Element('DoublePulse_pos_pt-%d' % i, pulsar=qt.pulsar)

                    e.add(pulse.cp(sq_pulsePH, amplitude=-0.7, length=self.params['PH_trigger_length']*1.0e-9), name='picoharp trigger', start=self.params['PH_trigger_time']*1.0e-9)
                    e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_init_length']*1.0e-9), name='laser init', start=self.params['AOM_start_buffer']*1.0e-9)
                    readout_start_time = self.params['AOM_start_buffer'] + self.params['AOM_init_length'] + self.params['tau_delay'][i]
                    e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_readout_length']*1.0e-9), name='laser readout', start=readout_start_time*1.0e-9)

                    trigger_period = self.params['AOM_start_buffer'] + self.params['AOM_init_length'] + self.params['tau_delay'][self.params['pts']-1] + self.params['AOM_readout_length'] + self.params['AOM_light_delay'] + self.params['AOM_end_buffer']
                    e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=trigger_period*1.0e-9),
                    name='MWimodpulse', start=0e-9)

                    e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=trigger_period*1.0e-9),
                    name='MWqmodpulse', start=0e-9)
                    elements.append(e)

                if j == 1 and self.params['microwaves']:
                    e = element.Element('DoublePulse_pos_mwave_pt_-%d' % i, pulsar=qt.pulsar)

                    e.add(pulse.cp(sq_pulsePH, amplitude=-0.7, length=self.params['PH_trigger_length']*1.0e-9), name='picoharp trigger', start=self.params['PH_trigger_time']*1.0e-9)
                    e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_init_length']*1.0e-9), name='laser init', start=self.params['AOM_start_buffer']*1.0e-9)
                    readout_start_time = self.params['AOM_start_buffer'] + self.params['AOM_init_length'] + self.params['tau_delay'][i]
                    e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=self.params['AOM_readout_length']*1.0e-9), name='laser readout', start=readout_start_time*1.0e-9)

                    # Calculate the center between the two AOM pulses, i.e. where to insert the pi-pulse. Split the pi pulse by 2 so its center is exactly at the center.
                    center_time = (self.params['AOM_readout_length'] - self.params['AOM_init_length'])/2.0 + self.params['AOM_start_buffer'] + self.params['AOM_init_length'] + self.params['AOM_light_delay'] - self.params['pi_length']/2.0
                    e.add(pulse.cp(sq_pulseMW, length = self.params['pi_length']*1.0e-9, amplitude = 1.0), name='microwave pi pulse', start=center_time*1.0e-9)

                    trigger_period = self.params['AOM_start_buffer'] + self.params['AOM_init_length'] + self.params['tau_delay'][self.params['pts']-1] + self.params['AOM_readout_length'] + self.params['AOM_light_delay'] + self.params['AOM_end_buffer']
                    e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=trigger_period*1.0e-9),
                    name='MWimodpulse', start=0e-9)

                    e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=trigger_period*1.0e-9),
                    name='MWqmodpulse', start=0e-9)
                    elements.append(e)


        if self.params['microwaves']:
            # Double the points if we're using microwaves.
            self.params['pts'] = 2*self.params['pts']
            ntd = []
            for tau in self.params['tau_delay']:
                for j in range(2):
                    ntd.append(tau)
            self.params['tau_delay'] = np.copy(np.array(ntd))

        seq = pulsar.Sequence('DoublePulse sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(3.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.5)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.3)

        if q >= 20:
            print 'AWG did not jump to proper waveform!'
        return
    def prepare(self):
        self.start_keystroke_monitor('abort')
        self._stop_measurement = False
        # Set up some instruments
        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        #self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']
        self._ph = qt.instruments['ph']

        # Prepare instruments for measurement and verify FBL output

        # Configure PH for histogram mode, initialize
        self._ph.start_histogram_mode()

        self._ph.set_Binning(self.params['Binning'])
        self._ph.set_InputCFD0(self.params['CFDLevel0'],self.params['CFDZeroCross0'])
        self._ph.set_InputCFD1(self.params['CFDLevel1'],self.params['CFDZeroCross1'])
        self._ph.set_SyncOffset(self.params['SyncOffset'])
        print 'PicoHarp settings configured.'
        # Check for counts on PH
        if self._ph.get_CountRate1() > 0:
            print 'PH 1 channel receiving counts...'
        else:
            print 'PH 1 not receiving counts!'

        # Set the PXI system
        self._pxi.close()
        self._pxi.init_device()
        self._pxi.reset_device()
        # Now set the new parameters
        self._pxi.set_power(self.params['power'])
        self._pxi.set_frequency((self.params['frequency'])*1.0e9)
        print 'PXI reset & configured: %.3f dBm at %.3f GHz' % (self.params['power'], self.params['frequency'])
        # Now set the proper attenuation
        desired_atten = self.params['power'] - self.params['constant_attenuation'] - self.params['desired_power']
        self._va.set_attenuation(desired_atten)
        print 'Variable attenuator set to %.1f dB attenuation -- should be %.3f' % (desired_atten, self.params['desired_power'])
        if self.params['microwaves']:
            self._pxi.set_status('on')
            print 'Microwave generation enabled.'
        else:
            self._pxi.set_status('off')
        time.sleep(2.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            self._stop_measurement = False
            return

        # set the AWG to CW mode
        print 'Booting up AWG.'
        self._awg.start()
        time.sleep(5.0)
        print 'Sequence force jump to CW waveform...'
        self._awg.sq_forced_jump(1)
        time.sleep(3.0)
        self.awg_confirm(1)
        print 'CW waveform confirmed.'


        self._fbl.optimize()
        print 'FBL optimized...'
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))



        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
            self._stop_measurement = False
            return

        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))

##        if self._snspd.check():
##            print 'SNSPD is superconducting.'
##        else:
##            print 'SNSPD is not superconducting!'
##            self._stop_measurement = False
##            return


        print 'Press q now to abort.'
        time.sleep(3.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return
        self._fbl.optimize()
        time.sleep(5.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._stop_measurement = True
            return

        return

    def measure(self):
        if self._stop_measurement:
            print 'Breaking from measurement sequence...'
            return

        print 'PicoHarp double pulse measurement sequence initiated.'
        # Wall time
        t0 = time.time()

        # Populate some arrays
        histogram_array =  np.zeros( (self.params['pts'], 65536), dtype='uint32')
        temporary_histogram_array =  np.zeros( (self.params['pts'], 65536), dtype='uint32')
        print '--Double pulse meas. from %.4f ns to %.4f ns in %.4f ns steps (%.2f steps), %.5f GHz --' % (self.params['tau_length_start'], self.params['tau_length_end'], self.params['tau_length_step'], self.params['pts'], self.params['frequency'])


        # Start the AWG sequencing
        self._awg.start()
        time.sleep(3)
        N_cmeas = 0
        # Set a time that controls when the next feedback occurs
        # Add a bit of randomness to this process
        # Optimize

        # Get the current sequence position index
        prev_awg_sq_position = int(self._awg.get_sq_position())
        # Now set the AWG into CW mode for tracking
        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)
        time.sleep(0.1)
        if self._fbl.optimize() == False:
            if self._fbl.optimize() == False:
                print 'FBL failed twice, breaking.'
        # Set the AWG back to the previous sequence position index
        self._awg.sq_forced_jump(prev_awg_sq_position)
        self.awg_confirm(prev_awg_sq_position)
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        scan_on = True
        for i in range(self.params['MeasCycles']):

            # Create an index of the waveforms so that we can modify it
            seq_index = range(self.params['pts'])
            if self.params['random'] == 1:
                # Now shuffle the array in place
                shuffle(seq_index)
            # Create an array for the single-sweep data
            temp_count_data = np.zeros(self.params['pts'], dtype='uint32')


            # Enter the loop for measurement
            t1 = time.time()
            for j in range(int(self.params['pts'])):

                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
                # Check if a track should occur. If so, track.
                if time.time() > track_time:
                    # Maybe should check if optimize is successful once that's robust
                    # Get the current sequence position index
                    prev_awg_sq_position = int(self._awg.get_sq_position())
                    # Now set the AWG into CW mode for tracking
                    self._awg.sq_forced_jump(1)
                    self.awg_confirm(1)

                    time.sleep(0.1)
                    # Re-optimize
                    self._fbl.optimize()

                    # Set the AWG back to the previous sequence position index
                    self._awg.sq_forced_jump(prev_awg_sq_position)
                    self.awg_confirm(prev_awg_sq_position)

                    # Set new track time
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()


                # Set the new double pulse separation
                self._awg.sq_forced_jump(seq_index[j]+2) # the +2 is because the indices start at 1, and the first sequence is CW mode
                self.awg_confirm(seq_index[j]+2)


                # Check signal count
                time.sleep(0.25)
                self._ni63.set_count_time(1.0)
                if msvcrt.kbhit():
                    kb_char=msvcrt.getch()
                    if kb_char == "q" :
                        scan_on = False
                        break
                signal_begin = self._ni63.get_ctr0()

                self._ph.ClearHistMem()
                self._ph.StartMeas(int(self.params['acquisition_time']*1000))
                print 'Acquiring signal for %s seconds, tau %.3f -- index %d of %d' % (self.params['acquisition_time'], self.params['tau_delay'][seq_index[j]], j+1, self.params['pts'])
                time.sleep(self.params['acquisition_time']+0.25)

                n = 0
                while self._ph.get_MeasRunning() and n < 10:
                    time.sleep(0.5)
                    n = n + 1
                if self._ph.get_MeasRunning():
                    print 'Measurement did not finish!'
                    break
                # Retrieve measurement
                current_data = self._ph.get_Histogram()
                print 'Got data for index %d' % (seq_index[j])
                temporary_histogram_array[seq_index[j],:] = current_data
                signal_end = self._ni63.get_ctr0()
                # Accumulate the current histogram
                print 'Signal count beginning %d, signal count end %d. Heater power at %.1f.' % (signal_begin, signal_end, self._ls332.get_heater_output())
                # If this is the last index j, then we've completed a measurement cycle and can accumulate the array into the final output array.
                # Doing this prevents some arrays from averaging more than others if a measurement is stopped before completing every measurement cycle.
                if int(j) == int(self.params['pts']-1):
                    histogram_array = histogram_array + temporary_histogram_array
                if int(j) == 1:
                    plott2d1 = qt.Plot2D(current_data, name='dp1_plot', clear=True)
                qt.msleep(0.02)



            # Check for a break, and break out of this loop as well.
            # It's important to check here, before we add the array to the total
            # since doing it the other way risks adding incomplete data to the
            # total array.
            tt = time.time() - t1

            print 'Cycle %d/%d: Total time is %.3f, efficiency of %.2f percent. Heater output at %.1f.' % (i+1, int(self.params['MeasCycles']), tt, (self.params['pts']*self.params['acquisition_time'])/tt*100.0, self._ls332.get_heater_output())



            plott2d = qt.Plot2D(histogram_array[int(self.params['pts']-1),:], name='dp_plot', clear=True)
            if msvcrt.kbhit() or scan_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or scan_on == False: break
            # Now start checking for other issues. If present, stop.
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
##            if self._snspd.check() == False:
##                print 'SNSPD went normal and could not restore, breaking.'
##                break
            # Checks have all passed, so proceed...
            qt.msleep(0.02)


        if self.params['microwaves'] == True:
            # Stop PXI sig gen
            self._pxi.set_status('off')
        # Set AWG to CW mode
        self._awg.sq_forced_jump(1)
        self.awg_confirm(1)
        # Measurement has ended, so start saving data

        grp = h5.DataGroup('SiC_DoublePulse', self.h5data, base=self.h5base)
        grp.add('tau_array', data=self.params['tau_delay'], unit='ns', note='double pulse separation lengths')
        grp.add('histogram_array', data=histogram_array, unit='counts', note='intermediate signal count data')

        return

class SiC_PLPolarizationBasic_Master(m2.Measurement):

    mprefix = 'PLPolarization'

    def sequence(self, upload = True, program=True, clear=True):
        self._awg = qt.instruments['awg']
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

        elements = []
        # Create waveform that has laser, microwaves, photon counting, and 1/0 I/Q modulation on
        # all the time for a long period of time (~100 us).
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='laser')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)





        seq = pulsar.Sequence('CW Polarization Sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)
        if clear:
            self._awg.clear_waveforms()
            print 'AWG waveforms cleared.'
        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(3.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.1)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.2)

        if q >= 20:
            print 'AWG did not jump to proper waveform!'
        return
    def prepare(self):
        self.start_keystroke_monitor('abort')
        # Set up some instruments

        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']
        self._st0 = qt.instruments['Standa0']
        self._flip = qt.instruments['flip']
        self._pm = qt.instruments['pm']


        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        self._awg.start()
        time.sleep(5.0)
        self._awg.sq_forced_jump(1)
        time.sleep(3.0)
        self.awg_confirm(1)


        self._fbl.optimize()
        print 'FBL optimized...'
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))



        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))


        self._data = qt.Data(name='PLpolarizationsweep')
        self._data.add_coordinate('Standa Step')
        self._data.add_value('PL Counts')
        self._data.add_value('Power (mW)')
        # Set the DAQ counter PFI channel/dwell
        self._ni63.set_count_time(self.params['dwell_time']/1000.0)
        self._ni63.set_ctr0_src(self.params['ctr_term'])
        print 'Counter prepared.'

        # set speed on standa
        self._st0.set_speed(int(self.params['speed']))

        time.sleep(1.0)


        return







    def measure(self):
        print 'Photoluminescence polarization measurement sequence initiated.'
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['step_high'] - self.params['step_low'])/self.params['step_step']))
        # Make an array of Standa polarization step values
        Mvals = np.linspace(self.params['step_low'], self.params['step_high'], n_steps)
        self._st0.move(int(Mvals[0]))

        self._stop_measurement = False
        # Take note of the time the measurement starts.
        t0 = time.time()
        m_pl_array_temp = np.zeros(n_steps)
        m_pl_array_total = np.zeros(n_steps)
        power_array = np.zeros(n_steps)
        deltaM = 0
        plot2d22 = qt.Plot2D(self._data, name='PLpolarization_intensity', coorddim=0, valdim=1)
        plot2d33 = qt.Plot2D(self._data, name='PLpolarization_power', coorddim=0, valdim=2)

        for k in range(int(self.params['MeasCycles'])):
            print 'Starting iteration %d of %d' % ((k+1),self.params['MeasCycles'])


            for ij in range(n_steps):


##                if time.time() > track_time or deltaM > 3000:
##                    print 'Tracking!'
##
##                    self._fbl.optimize()
##                    # Set new track time, fbl_time into the future plus a small
##                    # random time.
##                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
##                    deltaM = 0.0
                current_position = self._st0.get_position()
                self._st0.move(int(Mvals[ij]))
                deltaM = deltaM + (Mvals[ij] - current_position)


                print 'Standa to %d' % float(Mvals[ij])
                nn = 0
                while nn < 30:
                    time.sleep(1)
                    cur_pos = self._st0.get_position()
                    if np.abs(cur_pos - Mvals[ij]) < 2:
                        break
                    nn = nn + 1

                qt.msleep(0.1)
                self._fbl.optimize()
                self._ni63.set_count_time(self.params['dwell_time']/1000.0)
                self._flip.flip()
                curX = self._fsm.get_abs_positionX()
                curY = self._fsm.get_abs_positionY()
                self._fsm.zero()
                time.sleep(1.0)
                current_power = self._pm.get_power()*1.0e3
                self._flip.flip()
                self._fsm.set_abs_positionX(curX)
                self._fsm.set_abs_positionY(curY)
                time.sleep(1.0)
                m_pl_array_temp[ij] = self._ni63.get('ctr0')*1.0/(self.params['dwell_time']/1000.0)
                self._data.add_data_point(Mvals[ij],m_pl_array_temp[ij],current_power)
                qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.
                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q']:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._stop_measurement = True
                    break

            m_pl_array_total = m_pl_array_total + m_pl_array_temp

            if not self._stop_measurement:
                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q'] or self._stop_measurement:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._stop_measurement = True
                    print 'Breaking from internal measurement loop...'
                    break


            # Now start checking for other issues. If present, stop.
            # Check if setpoint and actual temperature are within a tolerance
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
            # Check if the SNSPD is still superconducting
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break

            if self._stop_measurement:
                print 'Breaking from main measurement loop...'
                break



        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_PLpolarization_data', self.h5data, base=self.h5base)
        grp.add('stepvals', data=Mvals, unit='step', note='step position')
        grp.add('pl_array', data=m_pl_array_total, unit='counts', note='pl counts')
        grp.add('power_array', data=power_array, unit='mW', note='power')

        return



class SiC_PLFieldsweep_Master(m2.Measurement):

    mprefix = 'PLfieldsweep'

    def sequence(self, upload = True, program=True):
        # define the pulses we'll use
        sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
        sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
        sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')
        sq_pulseMW_Imod = pulse.SquarePulse(channel='MW_Imod', name='A square pulse on MW I modulation')
        sq_pulseMW_Qmod = pulse.SquarePulse(channel='MW_Qmod', name='A square pulse on MW I modulation')

        elements = []
        # Create waveform that has laser, microwaves, photon counting, and 1/0 I/Q modulation on
        # all the time for a long period of time (~100 us).
        e = element.Element('CW_mode', pulsar=qt.pulsar)
        e.add(pulse.cp(sq_pulseAOM, amplitude=1, length=100e-6), name='laser')
        e.add(pulse.cp(sq_pulseMW, amplitude=1.0, length=100e-6), name='microwaves')
        e.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=100e-6), name='photoncountpulse')
        e.add(pulse.cp(sq_pulseMW_Imod, amplitude=1.0, length=100e-6),
        name='MWimodpulsecw', start=0e-9)
        e.add(pulse.cp(sq_pulseMW_Qmod, amplitude=0.0, length=100e-6),
        name='MWqmodpulsecw', start=0e-9)
        elements.append(e)





        seq = pulsar.Sequence('Fieldsweep Sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=False, repetitions=-1)

        if upload:
            qt.pulsar.upload(*elements)
            time.sleep(3.0)
        # program the AWG
        if program:
            qt.pulsar.program_sequence(seq)

    def awg_confirm(self, seq_el):
        q = 0
        time.sleep(0.1)
        while q < 20:

            cur_pos = int(self._awg.get_sq_position())
            if cur_pos == seq_el:
                break
            else:
                q = q + 1
            time.sleep(0.2)

        if q >= 20:
            print 'AWG did not jump to proper waveform!'
        return
    def prepare(self):
        self.start_keystroke_monitor('abort')
        # Set up some instruments

        self._fbl = qt.instruments['fbl']
        self._tl = qt.instruments['tl']
        self._ni63 = qt.instruments['NIDAQ6363']
        self._snspd = qt.instruments['snspd']
        self._fsm = qt.instruments['fsm']
        self._ls332 = qt.instruments['ls332']
        self._pxi = qt.instruments['pxi']
        self._ddg = qt.instruments['ddg']
        self._xps = qt.instruments['xps']
        self._awg = qt.instruments['awg']
        self._va = qt.instruments['va']


        # Prepare instruments for measurement and verify FBL output
        # Set the trigger source to internal

        # set the AWG to CW mode
        self._awg.start()
        time.sleep(5.0)
        self._awg.sq_forced_jump(1)
        time.sleep(3.0)
        self.awg_confirm(1)


        self._fbl.optimize()
        print 'FBL optimized...'
        # Set focus axis limit
        cur_Z = self._xps.get_abs_positionZ()
        self._xps.set_parameter_bounds('abs_positionZ',cur_Z-(self.params['focus_limit_displacement']*0.001),12.1)
        print 'Current Z is %.4f, focus limit set to %.4f' % (cur_Z, cur_Z-(self.params['focus_limit_displacement']*0.001))



        if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > 3.0:
            print 'Temperature away from setpoint!'
        else:
            print 'Temperature in reference (%.2f from setpoint), proceeding.' % (np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()))




        # Set the DAQ counter PFI channel (default is 'PFI0')
        self._ni63.set_ctr1_src(self.params['ctr_term'])
        print 'Counter prepared.'
        # Configure the exposure time in WinSpec



        time.sleep(1.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return



        self._data = qt.Data(name='PLfieldsweep')

        # Now you provide the information of what data will be saved in the
        # datafile. A distinction is made between 'coordinates', and 'values'.
        # Coordinates are the parameters that you sweep, values are the
        # parameters that you readout (the result of an experiment). This
        # information is used later for plotting purposes.
        # Adding coordinate and value info is optional, but recommended.
        # If you don't supply it, the data class will guess your data format.
        self._data.add_coordinate('Magnet Position (mm)')
        self._data.add_value('PL Counts')








    def measure(self):
        print 'Photoluminescence fieldsweep measurement sequence initiated.'
        track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
        n_steps = np.uint32(1 + np.ceil(np.abs(self.params['M_high'] - self.params['M_low'])/self.params['M_step']))
        # Make an array of magnet step values
        Mvals = np.linspace(self.params['M_low'], self.params['M_high'], n_steps)
        self._xps.set_abs_positionM(Mvals[0])

        self._stop_measurement = False
        # Take note of the time the measurement starts.
        t0 = time.time()
        m_pl_array_temp = np.zeros(n_steps)
        m_pl_array_total = np.zeros(n_steps)

        plot2d22 = qt.Plot2D(self._data, name='PLfs_intensity', coorddim=0, valdim=1)



        for k in range(int(self.params['MeasCycles'])):
            print 'Starting iteration %d of %d' % ((k+1),self.params['MeasCycles'])
            # Take the spectra in random order
            Mvals_temp = np.copy(Mvals)
            if self.params['random'] == 1:
                np.random.shuffle(Mvals_temp)

            for ij in range(n_steps):

                self._awg.sq_forced_jump(1)
                time.sleep(0.5)
                self.awg_confirm(1)

                if time.time() > track_time:
                    print 'Tracking!'

                    self._fbl.optimize()
                    # Set new track time, fbl_time into the future plus a small
                    # random time.
                    track_time = time.time() + self.params['fbl_time'] + 5.0*np.random.uniform()
                self._xps.set_abs_positionM(Mvals_temp[ij])
                if ij == 1:
                    qt.msleep(1.0)
                    self._fbl.optimize()
                print 'Magnet to %.3f' % float(Mvals_temp[ij])
                nn = 0
                while nn < 30:
                    time.sleep(0.5)
                    cur_pos = self._xps.get_abs_positionM()
                    if np.abs(cur_pos - Mvals_temp[ij]) < 0.01:
                        break
                    nn = nn + 1
                self._ni63.set_count_time(self.params['dwell_time']/1000.0)
                m_pl_array_temp[ij] = self._ni63.get('ctr0')*1.0/(self.params['dwell_time']/1000.0)
                self._data.add_data_point(Mvals_temp[ij],m_pl_array_temp[ij])
                qt.msleep(0.002) # keeps GUI responsive and checks if plot needs updating.
                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q']:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._stop_measurement = True
                    break
            sorted_temp_data = m_pl_array_temp[Mvals_temp.argsort()]
            m_pl_array_total = m_pl_array_total + sorted_temp_data

            if not self._stop_measurement:
                self._keystroke_check('abort')
                if self.keystroke('abort') in ['q','Q'] or self._stop_measurement:
                    print 'Measurement aborted.'
                    self.stop_keystroke_monitor('abort')
                    self._stop_measurement = True
                    print 'Breaking from internal measurement loop...'
                    break


            # Now start checking for other issues. If present, stop.
            # Check if setpoint and actual temperature are within a tolerance
            if np.abs(self._ls332.get_kelvinA() - self._ls332.get_setpoint1()) > self.params['temperature_tolerance']:
                print 'Temperature out of bounds, breaking.'
                break
            # Check if the SNSPD is still superconducting
            if self._snspd.check() == False:
                print 'SNSPD went normal and could not restore, breaking.'
                break

            if self._stop_measurement:
                print 'Breaking from main measurement loop...'
                break



        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_PLfieldsweep_data', self.h5data, base=self.h5base)
        grp.add('m_xvals', data=Mvals, unit='mm', note='magnet position array in mm')
        grp.add('m_pl_array_total', data=m_pl_array_total, unit='counts', note='pl counts array total')

        return