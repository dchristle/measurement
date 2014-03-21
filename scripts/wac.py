import numpy as np

varr = np.linspace(-3,6,12)

ni63 = qt.instruments['NIDAQ6363']
ni63.set_ctr0_src('PFI0')
ni63.set_count_time(0.2)
carray = ni63.write_and_count(varr,'ao0','ctr0')

print '%s' % carray

print '%s' % (carray[1]-carray[0])