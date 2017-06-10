from pylab import plot, show, figure, xlabel, ylabel, xlim, ylim, legend
from numpy.fft import fft, fftfreq, fftshift
from math import exp
from numpy import sqrt, diag
from scipy.optimize import curve_fit

def readData(filename):
    t,p1,p2 = [],[],[]
    
    file = open(filename, 'r')
    for line in file:
        currLine = map(float, line.split())
        t.append(currLine[0])
        p1.append(currLine[1])
        p2.append(currLine[2])
    return t, p1, p2

def plotData(t,p1,p2):
    figure()
    plot(t,p1,label='Pendulum 1')
    plot(t,p2,label='Pendulum 2')
    legend(loc=0)
    show()

def fourierTransform(t,p1,p2):
    freq = fftshift(fftfreq(len(t),t[1]-t[0]))
    f1 = fftshift(fft(p1))
    f2 = fftshift(fft(p2))
    
    return freq, f1, f2

# assume that the side peaks are at x0-x1 and x0+x1
# just three gaussians A*exp(-a*(x-x0)**2)
def fitTriplePeak(xData, x0,x1, a0,a1, A0,A1_l,A1_r, C):
    return map(lambda x: A0*exp(-a0*(x-x0)**2)+
                         A1_r*exp(-a1*(x-x0+x1)**2)+
                         A1_l*exp(-a1*(x-x0-x1)**2)+C, xData)

def fitFFT(freq, f, p0):
    params, cov = curve_fit(fitTriplePeak, freq,abs(f), p0=p0)
    errors = sqrt(diag(cov))
    
    frequencies = params[:2]
    errors_freq = errors[:2]
    
    return frequencies, errors_freq

# initial parameters determined by looking at plots
p0_1=[0,0.67, 6900,2400, 133063,51000,51000, 3]
p0_2=[0,0.67, 6000,3080, 64400,54000,54000, 3]