##v = np.random.normal(0,1,3)
### Normalize it
##v_norm = v/np.linalg.norm(v)
##print 'v norm is %s' % v_norm
### Compute Householder matrix. The [None] makes the array 2D so we can
### transpose it.
##H = np.eye(3) - 2.0*np.dot(v_norm[None].T,v_norm[None])
##
##print 'H is %s' % H



dmad = qt.instruments.create('dmad','dmads_opt')
dmad.add_control_variable('phase',ni1.set_phase,2,-35.,35.)
dmad.add_control_variable('angle1',ni1.set_angle1,0.6,0.,np.pi)
dmad.add_control_variable('angle2',ni1.set_angle2,0.5,0.,np.pi)
dmad.add_control_variable('cmax',ni1.set_cmax_position,3.1,1.0,5)
anon_f = lambda: -1*ni1.get_counts()
dmad.set_output_variable(anon_f)
dmad.optimize()