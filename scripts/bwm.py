import qt
import msvcrt
import numpy
import math
import analysis.lib.fitting.fit as fit
import analysis.lib.fitting.common as common


def array_erf(x):
    result = numpy.zeros(x.shape[0])
    for ij in range(0,x.shape[0]):
        result[ij] = erf(x[ij])

    return result
def erf(x):
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of x
    sign = 1
    if x < 0:
        sign = -1
    x = abs(x)

    # A & S 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)

    return sign*y
def fit_knife_simple(g_a, g_x0, g_w, g_b, *arg):
    """
    fits a knife edge function,
        y(x) = a/2 * (1-erf(sqrt(2)*(x-x0)/w)) + b

    Initial guesses, in this order:
        g_a : amplitude
        g_x0 : center position
        g_w : beam waist
        g_b : background


    """
    fitfunc_str = "a/2 * (1-erf(sqrt(2)*(x-x0)/w)) + b"

    a = fit.Parameter(g_a, 'a')
    x0 = fit.Parameter(g_x0, 'x0')
    w = fit.Parameter(g_w, 'w')
    b = fit.Parameter(g_b, 'b')
    # tau = fit.Parameter(g_tau, 'tau')
    p0 = [a, x0, w, b]

    def fitfunc(x) :

        return a()*0.5*(1 - array_erf(math.sqrt(2)*(x-x0())/w())) + b()
    #*(1 - math.erf(math.sqrt(2)*(x-x0())/w())) + b()
    return p0, fitfunc, fitfunc_str



# end damped rabi

def meas_BW(min_pos, max_pos, steps, name = 'beamwaist'):


#generate list of steps
    x_list = numpy.linspace(min_pos, max_pos, steps)

    ins_xps = qt.instruments['xps']
    ins_fm = qt.instruments['fm']


    # create data object
    qt.mstart()

    qt.msleep(0.2)


    d = qt.Data(name='BeamWaist')
    d.add_coordinate('displacement (mm)')
    d.add_value('power')
    d.create_file()
    filename=d.get_filepath()[:-4]

    plot2d = qt.Plot2D(d, name=name)
    stop_scan = False
    result = numpy.zeros(steps)
    for i,cur_x in enumerate(x_list):

        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): stop_scan=True
        ins_xps.set_abs_positionZ(float(-1*cur_x))

        qt.msleep(0.05)
        result[i] = ins_fm.get_power()
        d.add_data_point(cur_x, result[i])

        if stop_scan: break


    ins_xps.set_abs_positionZ(-1*min_pos)
    i_max = result.tolist().index(max(result))
    i_min = result.tolist().index(min(result))

    print 'result imax is: %s' % result[i_max]
    knife_fit = fit.fit1d(x_list, result,fit_knife_simple, result[i_max],
                    0, 2, result[i_min], do_print=False,ret=True)
    if type(knife_fit) != dict:
                print 'fit failed!'
    else:
        print 'best fit params are: A: %.03f x0: %.03f w: %.03f b: %.03f' % (
        knife_fit['params'][0], knife_fit['params'][1], knife_fit['params'][2],
        knife_fit['params'][3])
    plot2d.set_plottitle('Beam Waist: %.03f +- %.03f' % (knife_fit['params'][2], knife_fit['error'][2]))
    d.close_file()
    plot2d.save_png(filename+'.png')

    qt.mend()
