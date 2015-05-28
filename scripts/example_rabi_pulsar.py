# Attempt to program a basic Rabi sequence

import qt
import numpy as np

from measurement.lib.pulsar import pulse, pulselib, element, pulsar


reload(pulse)
reload(element)
reload(pulsar)
reload(pulselib)

test_element = element.Element('RabiElement', pulsar=qt.pulsar)
# we copied the channel definition from out global pulsar
# print 'Channel definitions: '
# pprint.pprint(test_element._channels)
# print

# define some bogus pulses.
sq_pulseAOM = pulse.SquarePulse(channel='AOM975', name='A square pulse on ThorLabsAOM')
sq_pulseMW = pulse.SquarePulse(channel='MW_pulsemod', name='A square pulse on MW modulation')
sq_pulsePC = pulse.SquarePulse(channel='photoncount', name='A square pulse on photon counting switch')

sq_delay = pulse.SquarePulse(channel='MW_Imod', name='delay',
            length = 200e-9, amplitude = 0.)
# create a few of those
test_element.add(pulse.cp(sq_pulseAOM, amplitude=1, length=1.8e-6),
    name='laser init')
test_element.add(pulse.cp(sq_delay, amplitude=0.0, length=0.8e-6),
	name='delay1', refpulse='laser init', refpoint='end')
test_element.add(pulse.cp(sq_pulseMW, amplitude=1, length=0.6e-6),
	name='microwave pulse', refpulse='delay1', refpoint='end')
test_element.add(pulse.cp(sq_delay, amplitude=0.0, length=0.3e-6),
	name='delay2', refpulse='microwave pulse', refpoint='end')
test_element.add(pulse.cp(sq_pulsePC, amplitude=1.0, length=0.21e-6, time_offset=400e-9),
    name='photoncountpulse')

#print 'Element overview:'
test_element.print_overview()
#print


# upload waveforms
qt.pulsar.upload(test_element)

# create the sequnce
# note that we re-use the same waveforms (just with different identifier
# names)
seq = pulsar.Sequence('RabiSequence')
seq.append(name='first element', wfname='RabiElement', trigger_wait=False,
    goto_target='first element', repetitions=-1)


# program the Sequence
qt.pulsar.program_sequence(seq)