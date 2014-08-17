
ni = qt.instruments.create('ni','noisy_instrument')
ni.set_sigma(1)
gps = qt.instruments.create('gps','gps_optimiz0r')

gps.set_output_variable(lambda : -1.0*ni.get_counts())

gps.add_control_variable('angle1',ni.set_angle1,20.0/180.0*np.pi,10.0/180*np.pi,70.0/180*np.pi)

gps.add_control_variable('angle2',ni.set_angle2,12.0/180.0*np.pi,1.0/180.*np.pi,70.0/180*np.pi)

gps.add_control_variable('cmax',ni.set_cmax_position,6.9,1.0,7.0)

gps.optimize()