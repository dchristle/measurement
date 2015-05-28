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

class SiC_SpectrumCalibrate_Master(m2.Measurement):

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
            time.sleep(2.0)
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
        self._mc = qt.instruments['mc_here'] # monochromator -- this is a remote instrument
        self._flip = qt.instruments['flip']
        
        # Prepare instruments for measurement and verify FBL output

        print 'Counter prepared.'
        # Configure the exposure time in WinSpec
        self._ws.set_exposure_time(self.params['exposure_time'])

        
        self._keystroke_check('abort')
        time.sleep(5.0)
        self._keystroke_check('abort')
        if self.keystroke('abort') in ['q','Q']:
            print 'Measurement aborted.'
            self.stop_keystroke_monitor('abort')
            self._pxi.set_status('off')
            return

        




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
        off_spec_total = np.zeros(1024,dtype=np.uint32)
        on_spec_total = np.zeros(1024,dtype=np.uint32)

        
        wl_array = self.wavelength_array(1150,2)
        for k in range(int(self.params['MeasCycles'])):
            print 'Starting iteration %d of %d' % ((k+1),self.params['MeasCycles'])
            # Take the spectra in random order
            seq_idx = range(2)
            if self.params['random'] == 1:
                shuffle(seq_idx)
                
            for j in range(2):
                if seq_idx[j] == 0:
                    # ON measurement
                    print 'PL iteration %d/%d, step %d/%d: ON spectrum, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 2, self.params['exposure_time'])

                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._pxi.set_status('off')
                        self._stop_measurement = True
                        break
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    on_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                    if not type(on_spec) == np.ndarray:
                        time.sleep(10.0)
                        on_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                        
                    
                    
                    
                if seq_idx[j] == 1:
                    # OFF measurement
                    print 'PL iteration %d/%d, step %d/%d: OFF spectrum, %d seconds' % ((k+1),self.params['MeasCycles'], j+1, 2, self.params['exposure_time'])
                    
                    self._keystroke_check('abort')
                    if self.keystroke('abort') in ['q','Q']:
                        print 'Measurement aborted.'
                        self.stop_keystroke_monitor('abort')
                        self._stop_measurement = True
                        break
                    print 'Flipping mirror...'
                    self._flip.flip()
                    time.sleep(1.5)
                    
                    print 'Exposing for %d seconds...' % self.params['exposure_time']
                    off_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+5.0))
                    if not type(off_spec) == np.ndarray:
                        time.sleep(10.0)
                        off_spec = self._ws.get_spectrum([],timeout=(self.params['exposure_time']+6.0))
                    print 'Flipping mirror back...'
                    self._flip.flip()
                    time.sleep(1.5)
                    
                    
 
                    
                    
                    
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
            if float(detector_temperature) > -80.0:
                print 'Spectrometer camera warming up -- breaking from main measurement loop...'
                break
            # Now start checking for other issues. If present, stop.

            
            if self._stop_measurement:
                print 'Breaking from main measurement loop...'
                break
            if not ((type(on_spec) == np.ndarray) and (type(off_spec) == np.ndarray)):
                break
            off_spec_total = off_spec_total + off_spec[0].T[0].astype(np.uint32)
            on_spec_total = on_spec_total + on_spec[0].T[0].astype(np.uint32)
            print 'size wl %s size on %s' % (wl_array.size, on_spec_total.size)
            qt.plot(wl_array,on_spec_total,name='on_spec',clear=True)
            qt.plot(wl_array,off_spec_total,name='off_spec',clear=True)


        # Measurement has ended, so start saving data
        grp = h5.DataGroup('SiC_Photoluminescence_data', self.h5data, base=self.h5base)
        grp.add('on_spec_total', data=on_spec_total, unit='counts', note='on total count array')
        grp.add('off_spec_total', data=off_spec_total, unit='counts', note='off total count array')

        return
