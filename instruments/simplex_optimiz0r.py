# generalized pattern search optimizer
# David Christle <christle@uchicago.edu>, August 2014
#
# This virtual instrument collects a dictionary of control functions
import qt
from instrument import Instrument
import numpy as np
import gobject
import instrument_helper
import types

class gps_optimizer(Instrument):
    def __init__(self, name, **kw):
        Instrument.__init__(self, name)


        ins_pars  = {'read_interval'    :   {'type':types.FloatType,'flags':Instrument.FLAG_GETSET},
                    'is_running'       :   {'type':types.BooleanType,'flags':Instrument.FLAG_GETSET},
                    }
        instrument_helper.create_get_set(self,ins_pars)

        self.add_function('start')
        self.add_function('stop')
        self.add_function('add_optimizer')
        self.add_function('get_optimizers')
        self.add_function('manual_optimize')

        self._optimizers=[]

        self.set_is_running(False)
        self.set_read_interval(60) #s
        self._is_waiting=False
        self._timer=-1

    #--------------get_set

    def get_optimizers(self):
        return self._optimizers

    def add_optimizer(self,name,**kw):
        ins_list=qt.instruments.get_instrument_names()
        try:
            ins_list.index(name)
        except ValueError:
            qt.instruments.create(name,'simple_optimizer',**kw)
        self._optimizers.append(name)

    def remove_optimizer(self,name):
        self._optimizers.remove(name)

    def _optimize(self,*arg):
        pass #implement in child class

    def _check_treshold(self):
        pass

    def manual_optimize(self,optimizer):
        self._optimize(optimizer)

    def start(self):
        if self.get_is_running():
            print 'ALREADY RUNNING'
            return False
        if self._is_waiting:
            print 'still waiting from previous run, will wait 20 sec extra.'
            qt.msleep(20)
            self._stop_waiting()
        self.set_is_running(True)
        self._timer=gobject.timeout_add(int(self.get_read_interval()*1e3),\
                self._check_treshold)
        return True

    def stop(self):
        self.set_is_running(False)
        return gobject.source_remove(self._timer)

    def start_waiting(self,time):
        self._is_waiting=True
        gobject.timeout_add(int(time*1e3),\
                self._stop_waiting)

    def _stop_waiting(self):
        self._is_waiting=False

        return False
    def _logscale(self):
        # Log scaling from Numerical Optimization
        # Somewhat incomplete, not all cases considered.
        lb = self.lb
        ub = self.ub
        eps_l = 1.0e-12
        lf = (np.abs(LB) >= eps_l & np.abs(UB) <= eps_l)
        lu = (np.abs(LB) >= eps_l & np.abs(UB) >= eps_l)
        uf = (np.abs(LB) <= eps_l & np.abs(UB) <= eps_l)

        log_s = np.zeros(lb.size)
        log_s[lf] = np.log2(np.abs(lb[lf]))
        log_s[lu] = (np.log2(np.abs(lb[lu])) + np.log2(np.abs(ub[lu])))/2.0
        log_s[uf] = np.log2(np.abs(ub[uf]))
        self._scale = np.power(2.0,np.round(log_s))
        return

    def _gps2nbasis(self, x):
        # Detect if bounds are active, then create the maximal 2n basis
        # In principle, other bases could also be generated here
        lb_a = (np.abs(x - self.lb) <= tol or np.abs(x-self.ub) <= tol)
        I_m = np.eye(self.lb.size)
        self.basis = I_m[:,not(lb_a)]
        self.tc = I_m[:,lb_a]
        # Figure out the number of basis vectors/vectors within the tangent cone
        Nbasis = np.sum(not(lb_a))
        Ntangent = np.sum(lb_a)
        Ntotal = Nbasis + Ntangent
        # Assemble the matrix of directions
        dirMat = np.hstack((self.basis, self.tc))
        # Index them
        indexV = np.concatenate( (np.arange(Nbasis), np.arange(Nbasis), np.arange(Nbasis,Ntotal), np.arange(Nbasis,Ntotal)) )
        signV = np.concatenate( (np.ones(Nbasis), -1*np.ones(Nbasis), np.ones(Ntangent), -1*np.ones(Ntangent) ) );
        # Create a vector for ordering the polling process -- not relevent
        # yet since we will do complete polling
        orderV = np.arange(indexV.size)
        # Iterate to create a list of polling locations
        new_x = np.zeros(indexV.size,self.lb.size)
        for i in range(indexV.size):
            direction = signV[i]*dirMat[:,indexV[orderV[k]]]
            new_x[k,:] = x + self._meshsize*self._scale*direction
        return new_x





    def __init__(self, name, mos_ins=qt.instruments['master_of_space'],
            adwin_ins=qt.instruments['adwin']):
        Instrument.__init__(self, name)

        self.add_function('optimize')
        self.mos = mos_ins
        self.adwin= adwin_ins

        ins_pars  = {'a'  :   {'type':types.FloatType,'val':1.5,'flags':Instrument.FLAG_GETSET},
                     'b'  :   {'type':types.FloatType,'val':1.3,'flags':Instrument.FLAG_GETSET},
                     'c'  :   {'type':types.FloatType,'val':0.5,'flags':Instrument.FLAG_GETSET},
                     'd'  :   {'type':types.FloatType,'val':0.9,'flags':Instrument.FLAG_GETSET},
                     }
        instrument_helper.create_get_set(self,ins_pars)



    def _measure_at_position(self, pos, int_time, cnt, speed):
        self.mos.move_to_xyz_pos(('x','y','z'),pos,speed,blocking=True)
        qt.msleep(0.001)
        return self.adwin.measure_counts(int_time)[cnt-1]

    def _tetra_volume(self, V):
        return 1/6.*np.abs(np.linalg.det(V[0:3]-V[3]))

    def optimize(self,xyz_range=[.5,0.5,1.0], xyz_tolerance_factor=0.002, max_cycles=15,
                  cnt=1, int_time=50, speed=2000,
                  do_final_countrate_check=True):

        #in principle, the method below works for any dimension D, and number of vertices N!
        #Only need to change the initial simplex shape to size N and set_functions to dim D
        D=3
        N=5


        a=self._a
        b=self._b
        c=self._c
        d=self._d
        search_range=np.array(xyz_range)
        tolerance=search_range*xyz_tolerance_factor
        pos = np.array([self.mos.get_x(), self.mos.get_y(), self.mos.get_z()])
        old_cnt = self.adwin.measure_counts(int_time)[cnt-1]
        new_pos=pos.copy()

        tetra=np.array([[ 0, 0, 0],
                        [ 1, 1, 1],
                        [ 1,-1,-1],
                        [-1, 1,-1],
                        [-1,-1, 1],
                      ],dtype=np.float)

        V=pos+tetra*search_range

        J=np.zeros(N)

        for i in range(1,len(J)):
          J[i]=self._measure_at_position(V[i],int_time,cnt,speed)
        print 'Simplex J countrates:',J/(int_time/1000.)

        for j in range(max_cycles):
          print '\n ============== \n j: ', j
          J[0]=self._measure_at_position(V[0],int_time,cnt,speed)
          J_o=np.sort(J)
          V_o=V[np.argsort(J)]

          v_c=1./(N-1.)*np.sum(V_o[1:],axis=0)
          v_r=(1-a)*V_o[0]+a*v_c

          j_r=self._measure_at_position(v_r,int_time,cnt,speed)
          print 'New vertex j_r countrates:',j_r/(int_time/1000.)

          if J_o[1]<=j_r<=J_o[-1]:
            V_o[0]=v_r
            print 'case 1: replace'
          elif j_r>J_o[-1]:
            v_e=b*v_r+(1-b)*v_c
            j_e=self._measure_at_position(v_e,int_time,cnt,speed)
            if j_e>J_o[-1]:
              V_o[0]=v_e #replace [0] or V[-1] here??
              print 'case 2.1: expand far'
            else:
              V_o[0]=v_r #replace [0] or V[-1] here??
              print 'case 2.2: expand close'
          elif j_r<J_o[1]:
            v_s=c*v_r+(1-c)*v_c
            j_s=self._measure_at_position(v_s,int_time,cnt,speed)
            if j_s>=J[0]:
              V_o[0]=v_s
              print 'case 3.1: stay close'
            else:
              print 'case 3.2: shrink'
              V_o=d*V_o+(1-d)*V_o[-1]
          else:
            print 'j_r comparison error. Check values.'

          V=V_o
          J=J_o

          var=np.var(V,axis=0)
          print 'Average xyz_var', np.sum(var)
          if (var<tolerance).all():
            new_pos = V[-1] if j_r<J_o[-1] else  V[0]
            print '========================='
            print 'Convex success'
            print '========================='
            break

          if (msvcrt.kbhit() and (msvcrt.getch() == 'q')):
            do_final_countrate_check = False
            print 'User interrupt'
            break

          #vol=self._tetra_volume(V) #N=4 simplex volume.
          #print 'vol', vol
          #if vol[j] < 0.0002:
          #  break

          #----------------------------------------------------


        if do_final_countrate_check:
          print "Proposed position x change %d nm" % \
                        (1000*V[-1][0]-1000*pos[0])
          print "Proposed position y change %d nm" % \
                        (1000*V[-1][1]-1000*pos[1])
          print "Proposed position z change %d nm" % \
                        (1000*V[-1][2]-1000*pos[2])
          new_cnt = self._measure_at_position(V[-1], int_time, cnt, speed/2.)
          print 'Old countrates', old_cnt/(int_time/1000.)
          print 'New countrates', new_cnt/(int_time/1000.)
          if new_cnt>old_cnt:
            print 'New position accepted'
          else:
            print 'Old position kept'
            self.mos.move_to_xyz_pos(('x','y','z'),pos,speed,blocking=True)

        else:
          print "Position x changed %d nm" % \
                          (1000*new_pos[0]-1000*pos[0])
          print "Position y changed %d nm" % \
                          (1000*new_pos[1]-1000*pos[1])
          print "Position z changed %d nm" % \
                          (1000*new_pos[2]-1000*pos[2])
          print 'Old countrates', old_cnt/(int_time/1000.)
          print  "Countrates at new position: %d" % \
                      (float(self._measure_at_position(new_pos, int_time,cnt, speed))/(int_time/1000.))

