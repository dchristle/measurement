import msvcrt
import numpy as np
import qt
import hdf5_data as h5
import logging

import measurement.lib.config.adwins as adwins_cfg
import measurement.lib.measurement2.measurement as m2
from measurement.lib.measurement2.adwin_ssro import ssro
from measurement.lib.pulsar import pulse, pulselib, element, pulsar

class PulsarMeasurement(ssro.IntegratedSSRO):
    mprefix = 'PulsarMeasurement'
    awg = None
    mwsrc = None

    def __init__(self, name):
        ssro.IntegratedSSRO.__init__(self, name)

        self.params['measurement_type'] = self.mprefix

    def setup(self, wait_for_awg=True, mw=True):
        ssro.IntegratedSSRO.setup(self)

        if mw:
            self.mwsrc.set_iq('on')
            self.mwsrc.set_pulm('on')
            self.mwsrc.set_frequency(self.params['mw_frq'])
            self.mwsrc.set_power(self.params['mw_power'])
            self.mwsrc.set_status('on')

        self.awg.set_runmode('SEQ')
        self.awg.start()

        if wait_for_awg:
            i=0
            awg_ready = False
            while not awg_ready and i<40:
                if (msvcrt.kbhit() and (msvcrt.getch() == 'x')):
                    raise Exception('User abort while waiting for AWG')

                try:
                    if self.awg.get_state() == 'Waiting for trigger':
                        awg_ready = True
                except:
                    print 'waiting for awg: usually means awg is still busy and doesnt respond'
                    print 'waiting', i, '/ 40'
                    i=i+1

                qt.msleep(0.5)

            if not awg_ready:
                raise Exception('AWG not ready')

    def generate_sequence(self):
        pass

    def stop_sequence(self):
        self.awg.stop()

    def finish(self,**kw):
        ssro.IntegratedSSRO.finish(self,**kw)

        self.awg.stop()
        self.awg.set_runmode('CONT')

        self.mwsrc.set_status('off')
        self.mwsrc.set_iq('off')
        self.mwsrc.set_pulm('off')

# class PulsarMeasurement

class DarkESR(PulsarMeasurement):
    mprefix = 'PulsarDarkESR'

    def autoconfig(self):
        PulsarMeasurement.autoconfig(self)

        self.params['sequence_wait_time'] = \
            int(np.ceil(self.params['pulse_length']*1e6)+15)

        self.params['sweep_name'] = 'MW frq (GHz)'
        self.params['sweep_pts'] = (np.linspace(self.params['ssbmod_frq_start'],
            self.params['ssbmod_frq_stop'], self.params['pts']) + \
                self.params['mw_frq'])*1e-9

    def generate_sequence(self, upload=True):

        # define the necessary pulses
        X = pulselib.MW_IQmod_pulse('Weak pi-pulse',
            I_channel='MW_Imod', Q_channel='MW_Qmod',
            PM_channel='MW_pulsemod',
            amplitude = self.params['ssbmod_amplitude'],
            length = self.params['pulse_length'],
            PM_risetime = self.params['MW_pulse_mod_risetime'])

        T = pulse.SquarePulse(channel='MW_Imod', name='delay')
        T.amplitude = 0.
        T.length = 2e-6

        # make the elements - one for each ssb frequency
        elements = []
        for i, f in enumerate(np.linspace(self.params['ssbmod_frq_start'],
            self.params['ssbmod_frq_stop'], self.params['pts'])):

            e = element.Element('DarkESR_frq-%d' % i, pulsar=qt.pulsar)
            e.add(T, name='wait')
            e.add(X(frequency=f), refpulse='wait')
            elements.append(e)

        # create a sequence from the pulses
        seq = pulsar.Sequence('DarkESR sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=True)

        # upload the waveforms to the AWG
        if upload:
            qt.pulsar.upload(*elements)

        # program the AWG
        qt.pulsar.program_sequence(seq)

        # some debugging:
        # elements[-1].print_overview()


class ElectronRabi(PulsarMeasurement):
    mprefix = 'ElectronRabi'


    def autoconfig(self):
        PulsarMeasurement.autoconfig(self)

        self.params['sequence_wait_time'] = \
            int(np.ceil(np.max(self.params['MW_pulse_durations'])*1e6)+10)


    def generate_sequence(self, upload=True):

        # define the necessary pulses
        X = pulselib.MW_IQmod_pulse('Weak pi-pulse',
            I_channel='MW_Imod', Q_channel='MW_Qmod',
            PM_channel='MW_pulsemod',
            frequency = self.params['MW_pulse_frequency'],
            PM_risetime = self.params['MW_pulse_mod_risetime'])

        T = pulse.SquarePulse(channel='MW_Imod', name='delay',
            length = 200e-9, amplitude = 0.)

        # make the elements - one for each ssb frequency
        elements = []
        for i in range(self.params['pts']):

            e = element.Element('ElectronRabi_pt-%d' % i, pulsar=qt.pulsar)
            e.append(T)

            e.append(pulse.cp(X,
                length = self.params['MW_pulse_durations'][i],
                amplitude = self.params['MW_pulse_amplitudes'][i]))

            elements.append(e)


        # create a sequence from the pulses
        seq = pulsar.Sequence('ElectronRabi sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=True)

        # upload the waveforms to the AWG
        if upload:
            qt.pulsar.upload(*elements)

        # program the AWG
        qt.pulsar.program_sequence(seq)

        # some debugging:
        # elements[-1].print_overview()

class ElectronRamsey(PulsarMeasurement):
    mprefix = 'ElectronRamsey'

    def autoconfig(self):
        PulsarMeasurement.autoconfig(self)

        self.params['sequence_wait_time'] = \
            int(np.ceil(np.max(self.params['evolution_times'])*1e6)+10)

    def generate_sequence(self, upload=True):

        # define the necessary pulses
        X = pulselib.IQ_CORPSE_pi2_pulse('CORPSE pi2-pulse',
            I_channel='MW_Imod',
            Q_channel='MW_Qmod',
            PM_channel='MW_pulsemod',
            PM_risetime = self.params['MW_pulse_mod_risetime'],
            frequency = self.params['CORPSE_pi2_mod_frq'],
            amplitude = self.params['CORPSE_pi2_amp'],
            length_24p3 = self.params['CORPSE_pi2_24p3_duration'],
            length_m318p6 = self.params['CORPSE_pi2_m318p6_duration'],
            length_384p3 = self.params['CORPSE_pi2_384p3_duration'])

        T = pulse.SquarePulse(channel='MW_Imod', name='delay',
            length = 200e-9, amplitude = 0.)

        # make the elements - one for each ssb frequency
        elements = []
        for i in range(self.params['pts']):

            e = element.Element('ElectronRamsey_pt-%d' % i, pulsar=qt.pulsar,
                global_time = True)
            e.append(T)

            e.append(pulse.cp(X,
                amplitude = self.params['CORPSE_pi2_amps'][i],
                phase = self.params['CORPSE_pi2_phases1'][i]))

            e.append(pulse.cp(T,
                length = self.params['evolution_times'][i]))

            e.append(pulse.cp(X,
                amplitude = self.params['CORPSE_pi2_amps'][i],
                phase = self.params['CORPSE_pi2_phases2'][i]))

            elements.append(e)

        # create a sequence from the pulses
        seq = pulsar.Sequence('ElectronRamsey sequence')
        for e in elements:
            seq.append(name=e.name, wfname=e.name, trigger_wait=True)

        # upload the waveforms to the AWG
        if upload:
            qt.pulsar.upload(*elements)

        # program the AWG
        qt.pulsar.program_sequence(seq)

        # some debugging:
        # elements[-1].print_overview()

class MBI(PulsarMeasurement):
    mprefix = 'PulsarMBI'
    adwin_process = 'MBI'

    def autoconfig(self):

        self.params['sweep_length'] = self.params['pts']
        self.params['repetitions'] = \
                self.params['nr_of_ROsequences'] * \
                self.params['pts'] * \
                self.params['reps_per_ROsequence']

        # Calling autoconfig from sequence.SequenceSSRO and thus
        # from ssro.IntegratedSSRO after defining self.params['repetitions'],
        # since the autoconfig of IntegratedSSRO uses this parameter.
        PulsarMeasurement.autoconfig(self)

        self.adwin_process_params['Ex_MBI_voltage'] = \
            self.E_aom.power_to_voltage(
                    self.params['Ex_MBI_amplitude'])

    def run(self, autoconfig=True, setup=True):
        if autoconfig:
            self.autoconfig()

        if setup:
            self.setup()

        for i in range(10):
            self.physical_adwin.Stop_Process(i+1)
            qt.msleep(0.1)

        qt.msleep(1)
        self.adwin.load_MBI()
        qt.msleep(1)

        length = self.params['nr_of_ROsequences']
        self.physical_adwin.Set_Data_Long(
                np.array(self.params['A_SP_durations'], dtype=int), 33, 1, length)
        self.physical_adwin.Set_Data_Long(
                np.array(self.params['E_RO_durations'], dtype=int), 34, 1, length)

        v = [ self.A_aom.power_to_voltage(p) for p in self.params['A_SP_amplitudes'] ]
        self.physical_adwin.Set_Data_Float(np.array(v), 35, 1, length)

        v = [ self.E_aom.power_to_voltage(p) for p in self.params['E_RO_amplitudes'] ]
        self.physical_adwin.Set_Data_Float(np.array(v), 36, 1, length)

        self.physical_adwin.Set_Data_Long(
                np.array(self.params['send_AWG_start'], dtype=int), 37, 1, length)

        self.physical_adwin.Set_Data_Long(
                np.array(self.params['sequence_wait_time'], dtype=int), 38, 1, length)

        self.start_adwin_process(stop_processes=['counter'], load=False)
        qt.msleep(1)
        self.start_keystroke_monitor('abort')

        while self.adwin_process_running():

            if self.keystroke('abort') != '':
                print 'aborted.'
                self.stop_keystroke_monitor('abort')
                break

            reps_completed = self.adwin_var('completed_reps')
            print('completed %s / %s readout repetitions' % \
                    (reps_completed, self.params['repetitions']))
            qt.msleep(1)

        try:
            self.stop_keystroke_monitor('abort')
        except KeyError:
            pass # means it's already stopped

        self.stop_adwin_process()
        reps_completed = self.adwin_var('completed_reps')
        print('completed %s / %s readout repetitions' % \
                (reps_completed, self.params['repetitions']))


    def save(self, name='adwindata'):
        reps = self.adwin_var('completed_reps')
        sweeps = self.params['pts'] * self.params['reps_per_ROsequence']

        self.save_adwin_data(name,
                [   ('CR_before', sweeps),
                    ('CR_after', sweeps),
                    ('MBI_attempts', sweeps),
                    ('statistics', 10),
                    ('ssro_results', reps),
                    ('MBI_cycles', sweeps),
                    ('MBI_time', sweeps),  ])
        return

    def _MBI_element(self):
        # define the necessary pulses
        T = pulse.SquarePulse(channel='MW_pulsemod',
            length = 10e-9, amplitude = 0)

        X = pulselib.MW_IQmod_pulse('MBI MW pulse',
            I_channel = 'MW_Imod', Q_channel = 'MW_Qmod',
            PM_channel = 'MW_pulsemod',
            frequency = self.params['AWG_MBI_MW_pulse_ssbmod_frq'],
            amplitude = self.params['AWG_MBI_MW_pulse_amp'],
            length = self.params['AWG_MBI_MW_pulse_duration'],
            PM_risetime = self.params['MW_pulse_mod_risetime'])

        adwin_sync = pulse.SquarePulse(channel='adwin_sync',
            length = (self.params['AWG_to_adwin_ttl_trigger_duration'] \
                + self.params['AWG_wait_for_adwin_MBI_duration']),
            amplitude = 2)

        # the actual element
        mbi_element = element.Element('MBI CNOT', pulsar=qt.pulsar)
        mbi_element.append(T)
        mbi_element.append(X)
        mbi_element.append(adwin_sync)

        return mbi_element







