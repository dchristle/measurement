
import numpy as np

##float s = atof( argv[1] ); // seed value e.g. 1.0
##int p = atoi( argv[2] ); // sampled # of indices e.g 64
##int z = atoi( argv[3] ); // total # of indices e.g. 256
##int i; // Fourier grid index, e.g. 1 through 256
##int k; // generated gap size
##int n; // temporary # indices
##int *v; // temporary storage vector
##int j; // wrking variable
##float ld = (float) z / (float) p; // establish 1/fraction
##float adj = 2.0*(ld-1); // initial guess of adjustment
##srand48( s );
##v = ( int* ) malloc( z*sizeof( z ) );
##do {
##i = 0; n = 0;
##while ( i < z ) {
##v[n] = i;
##i += 1;
##k =poisson(adj*sin((float)(i+0.5)/(float)(z+1)*1.5707963268) );
##i += k;
##n += 1;
##}
##if ( n > p ) adj *= 1.02; // too many pts created
##if ( n < p ) adj /= 1.02; // too few pts created
##} while ( n != p ); // if not at first, try, try again
##for ( j = 0 ; j < p ; k++ ) printf( "%d\n", v[j] );
##}

seed = 1.0
p = 24
z = 640


v = np.zeros(z)
adj = 2.0*(z/p-1.0)
qq = 0
i = 0.0
n = 0.0
while n != p and qq < 200:
    zz = 0.0
    i = 0.0
    n = 0.0
    v = np.zeros(z)
    while i < z and zz < 300:
        v[n] = i
        i = i+1
        #print 'arg is %.3f, %.2f, %.2f' % (adj*np.sin((float(i)+0.5)/(float(z)+1.0)*np.pi/2.0),i,z)
        k = np.random.poisson(adj*np.sin((float(i)+0.5)/(float(z)+1.0)*np.pi/2.0))
        i = i+k
        n = n+1
        zz = zz+1
    if n > p:
        adj = adj*1.02
    if n < p:
        adj = adj/1.02
    qq = qq+1

print 'qq reached %.2f, zz reached %.2f' %(qq,zz)

print v[0:(p-1)]*10.0