# Exponential model slice sampling
# David Christle <christle@uchicago.edu>

import numpy
import numpy.random
import time

import matplotlib.pyplot as plt



def main():
##    d = qt.Data()
##    d.add_coordinate('i')
##    d.add_value('A')
##    d.add_value('tau')
##    d.add_value('B')
    domain_min = 0
    domain_max = 25
    search_points = 100
    search_array = numpy.linspace(domain_min,domain_max,search_points)
    data_x = numpy.linspace(0,10,12)

    data_y = 1.2+numpy.exp(-data_x/17.7)+numpy.random.randn(data_x.shape[0])/10

    data_sigma = numpy.ones(data_x.shape[0])*1/10

    llh = exp_llh(numpy.array([1.1,3.2,0.1],float),data_x,data_y,data_sigma)
    Nsamples = 10
    x_out = numpy.zeros([Nsamples,3])
    for i in range(Nsamples):
        xx, last_llh = slicesample(numpy.array([1,4,0.1],float), lambda x: exp_llh(x,data_x,data_y,data_sigma), None, 10.1, False)
        x_out[i,0:3] = xx
##        d.add_data_point(i,xx[0],xx[1],xx[2])



##    p = qt.Plot2D()
##    p.add_data(d, coorddim=0, valdim=1)
##    p.add_data(d, coorddim=0, valdim=2)
##    p.add_data(d, coorddim=0, valdim=3)

    print '%s' % x_out
    print 'Begin search.'
    inf_gain = numpy.zeros(Nsamples,float)
    exp_inf_gain = numpy.zeros(search_points,float)
    for ij in range(search_points):
        for im in range(Nsamples):
            ##print '%s xout' % x_out[im,:]
            sigma = 0.1
            pred = exp_prediction_draw(search_array[ij],x_out[im,:], sigma)
            est = exp_sampling_dist(search_array[ij],pred,x_out, sigma)
            log_est = numpy.log(est)
            inf_gain[im] = numpy.mean(log_est)

        exp_inf_gain[ij] = -numpy.mean(inf_gain)
        print 'est. log average %.3f' % exp_inf_gain[ij]

    max_idx = numpy.argmax(exp_inf_gain)
    print 'Found highest information gain at %.2f of %.2f' % (search_array[max_idx],exp_inf_gain[max_idx])
    x_data = data_x
    i_chosen = numpy.random.randint(Nsamples)
    y1 = exp_prediction(x_data,x_out[i_chosen])
    i_chosen = numpy.random.randint(Nsamples)
    y2 = exp_prediction(x_data,x_out[i_chosen])
    i_chosen = numpy.random.randint(Nsamples)
    y3 = exp_prediction(x_data,x_out[i_chosen])
    plt.ion()
    plt.show()
    for jj in range(50):
        i_chosen = numpy.random.randint(Nsamples)
        y4 = exp_prediction(x_data,x_out[i_chosen])

        plt.plot(x_data, y4, 'r--')
        plt.draw()
        time.sleep(0.05)
    return

def exp_prediction(x_0, params):
    output = params[0]*numpy.exp(-x_0/params[1]) + params[2]
    return output

def exp_prediction_draw(x_0, params, sigma):
    mean = params[0]*numpy.exp(-x_0/params[1]) + params[2]
    output = numpy.random.randn(1)*sigma+mean
    return output

def exp_sampling_dist(x_0, datum, params, noise):
    yd = params[:,0]*numpy.exp(-x_0/params[:,1]) + params[:,2] - numpy.ones(params[:,1].shape[0])*datum
    output = 1/numpy.sqrt(2*numpy.pi*numpy.power(noise,2)) * numpy.exp(-numpy.power(yd,2)/(2*numpy.power(noise,2)))
    return output

def exp_llh(params, data_x, data_y, data_sigma):
    numpy.set_printoptions(precision=3)
    A = params[0]
    tau = params[1]
    B = params[2]

    llh = numpy.sum(numpy.log(data_sigma)) - numpy.sum(numpy.power((data_y-(A*numpy.exp(-data_x/tau)+B))/(2*numpy.power(data_sigma,2)),2))
    return llh
# This slice sample code is originally from:
# http://machinelearningmisc.blogspot.com/2010/01/slice-sampling-in-python.html
def slicesample(xx, llh_func, last_llh=None, sigma=1, step_out=True):
    dims = xx.shape[0]
    perm = range(dims)
    numpy.random.shuffle(perm)

    if (type(sigma).__name__ == 'int') or (type(sigma).__name__ == 'float'):
        sigma = numpy.tile(sigma, dims)

    elif (type(sigma).__name__ == 'tuple') or (type(sigma).__name__ == 'list'):
        sigma = numpy.array(sigma)

    if last_llh is None:
        last_llh = llh_func(xx)
    for d in perm:
        llh0   = last_llh.copy()
        rr     = numpy.random.rand(1)
        x_l    = xx.copy()
        x_l[d] = x_l[d] - rr*sigma[d]
        x_r    = xx.copy()
        x_r[d] = x_r[d] + (1-rr)*sigma[d]

        if step_out:
            llh_l = llh_func(x_l)
            while llh_l > llh0:
                x_l[d] = x_l[d] - sigma[d]
                llh_l  = llh_func(x_l)
                print 'Stepping dim %s left out to %s with ll %s...' % (d, x_l[d], llh_l)
            llh_r = llh_func(x_r)
            while llh_r > llh0:
                x_r[d] = x_r[d] + sigma[d]
                llh_r  = llh_func(x_r)
                print 'Stepping dim %s right out to %s with ll %s, to get below %s and %s...' % (d, x_r[d], llh_r, llh0, last_llh)
        x_cur = xx.copy()
        while True:

            xd = numpy.random.rand(1)*(x_r[d] - x_l[d]) + x_l[d]


            x_cur[d] = xd
            last_llh = llh_func(x_cur)
            if last_llh > llh0:

                xx[d] = xd.copy()

                break
            elif xd > xx[d]:
                x_r[d] = xd.copy()
            elif xd < xx[d]:
                x_l[d] = xd.copy()
            else:
                raise RuntimeException("Slice sampler shrank too far.")
    return xx, last_llh