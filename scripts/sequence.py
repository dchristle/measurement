import qt

from measurement.lib.pulsar import pulsar
reload(pulsar)

pulsar.Pulsar.AWG = qt.instruments['awg']

# FIXME in principle we only want to create that once, at startup
try:
    del qt.pulsar
except:
    pass
qt.pulsar = pulsar.Pulsar()

### channels
# RF



qt.pulsar.define_channel(id='ch1', name='RF', type='analog', high=1.5,
    low=-1.5, offset=0., delay=0*1e-9, active=True)

# MW
qt.pulsar.define_channel(id='ch1_marker1', name='MW_pulsemod', type='marker',
    high=2.5, low=0, offset=0., delay=0.*1.0e-9, active=True)
qt.pulsar.define_channel(id='ch3', name='MW_Imod', type='analog', high=0.95,
    low=-0.95, offset=0., delay=0.*1.0e-9, active=True)
qt.pulsar.define_channel(id='ch4', name='MW_Qmod', type='analog', high=0.95,
    low=-0.95, offset=0., delay=0.*1.0e-9, active=True)

# sync DDG
qt.pulsar.define_channel(id='ch3_marker2', name='ddg_sync', type='marker',
    high=2.5, low=0, offset=0., delay=0., active=True)


# light
qt.pulsar.define_channel(id='ch2_marker1', name='AOM975', type='marker',
    high=2.5, low=0, offset=0., delay=0, active=True)
qt.pulsar.define_channel(id='ch2_marker2', name='Sacher1160AOM', type='marker',
    high=2.5, low=0, offset=0., delay=0, active=True)
# photon counting
qt.pulsar.define_channel(id='ch3_marker1', name='photoncount', type='marker',
    high=2.5, low=0, offset=0., delay=0.0e-9, active=True)
qt.pulsar.define_channel(id='ch4_marker2', name='phtrigger', type='marker',
    high=0.0, low=-0.7, offset=0.0, delay=0.0e-9, active=True)
##qt.pulsar.set_channel_opt('Velocity1AOM','high', qt.instruments['MatisseAOM'].get_sec_V_max())
##qt.pulsar.set_channel_opt('Velocity1AOM','low', qt.instruments['MatisseAOM'].get_sec_V_off())
##qt.pulsar.define_channel(id='ch1_marker2', name='YellowAOM', type='marker',
##    high=0.4, low=0, offset=0., delay=750e-9, active=True)
##qt.pulsar.set_channel_opt('YellowAOM','high', qt.instruments['YellowAOM'].get_sec_V_max())
##qt.pulsar.set_channel_opt('YellowAOM','low', qt.instruments['YellowAOM'].get_sec_V_off())

#qt.pulsar.define_channel(id='ch2', name='Velocity1AOM', type='analog',
#    high=0.4, low=0, offset=0., delay=700e-9, active=True)
#qt.pulsar.define_channel(id='ch2_marker2', name='YellowAOM', type='marker',
#    high=0.4, low=0, offset=0., delay=750e-9, active=True)

# Trigger AWG LT2
#qt.pulsar.define_channel(id='ch3_marker1', name='AWG_LT2_trigger', type='marker',
#    high=2.0, low=0, offset=0., delay=0., active=True)