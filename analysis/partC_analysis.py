from pylab import plot, show, figure, xlabel,ylabel, xlim,ylim, errorbar, legend, subplot, savefig
from numpy.fft import fft,fftfreq,fftshift
from math import pi
from scipy.optimize import curve_fit
from numpy import sqrt, diag

k = {1:3.20778210044,
     2:30.4874228538,
     3:165.561223371,
     4:55.4326207947}
l = {1:0.1330,
     2:0.2230,
     3:0.3080,
     4:0.4180}
L = 0.5705
m1,m2 = 2.931,2.946
M = m1*m2/(m1+m2)
lError = 0.003
lError= 0.003


def readData(filename):
    t,p1,p2 = [],[],[]
    
    file = open(filename, 'r')
    for line in file:
        currLine = map(float, line.split())
        t.append(currLine[0])
        p1.append(currLine[1])
        p2.append(currLine[2])
    return t, p1, p2

def naturalFrequency(f):
    return 2*pi*f
    
def fourierTransform(t,p1,p2):
    freq = map(naturalFrequency, fftfreq(len(t),t[1]-t[0]))
    f1 = abs(fft(p1))
    f2 = abs(fft(p2))
    
    # set zero frequency components to 0
    f1[0] = 0
    f2[0] = 0
    
    freq = fftshift(freq)
    f1 = fftshift(f1)
    f2 = fftshift(f2)
    
    # Truncate data to between -50 and +50
    startNdx = findClosest(freq,-50)
    endNdx = findClosest(freq,50)+1
    freq = freq[startNdx:endNdx]
    f1 = f1[startNdx:endNdx]
    f2 = f2[startNdx:endNdx]
    
    return freq, f1, f2

# dx is a list of the errors in the independent variables x
def errorAddition(dx):
    return sum(map(lambda x: x**2, dx))**0.5

# y is the value calculated, x,dx are lists of variables and their errors
def errorMultiplication(y,x,dx):
    return y*(sum([(dx[i]/x[i])**2 for i in range(len(x))])**0.5)

# y=x^n
def errorPower(y,x,dx,n):
    return abs(n)*y*dx/n

# y = ax^n
# dy = n*a*x^(n-1)dx, ax^(n-1)=y/x
# dy = n*y*dx/x

# y = a/x
# dy = a/x^2dx
# dy = a/x *dx/x = y*dx/x

def expectedFrequency(k,l,m,L):
    g = 9.81
    return (g/L)**0.5, ((g/L)+2*k*l**2/(m*L**2))**0.5

def expectedFrequency2(k,l,m,L):
    g = 9.81
    w_0 = (g/L)**0.5
    return w_0, w_0+2*k*l**2/(w_0*m*L**2)

def expectedFrequencyError(k,l,m,L,kError,lError,mError,LError):
    f0,f1 = expectedFrequency(k,l,m,L)
    
    g = 9.81
    
    f1_sqr_1 = 2*k*l**2/(m*L**2)
    f1_sqr_error = errorAddition([errorPower(g/L,L,LError,-1), errorMultiplication(f1_sqr_1,[k,l**2,m,L**2],[kError,errorPower(l**2,l,lError,2),mError,errorPower(L**2,L,LError,2)])])
    
    error_f0 = errorPower(f0,L,LError,0.5)
    error_f1 = errorPower(f1,f1**2,f1_sqr_error,0.5)
    
    return error_f0, error_f1

def findClosest(x,target):
    error = float('inf')
    ans = 0
    for i in range(len(x)):
        if abs(x[i]-target)<error:
            error = abs(x[i]-target)
            ans = i
    return ans

def Lorentz(x, A,x0,w):
    return A*w**2/(4*(x-x0)**2+w**2)

def fourPeakLorentz(xData, A0,A1,x1,x0,w0,w1,C):
    return map(lambda x: Lorentz(x, A0,x0,w0)+Lorentz(x, A0,-x0,w0)+
                         Lorentz(x, A1,x1,w1)+Lorentz(x, A1,-x1,w1)+C, xData)

def findPeaks(x,y, peak_guess):
    A = max(y) # assume the peaks are of comparable heights
    w = 2*(x[1]-x[0]) # since we're fitting a fourier transform
    
    params,cov = curve_fit(fourPeakLorentz, x,y, p0=[A,A,peak_guess[0],peak_guess[1],w,w,min(y)])
    
    return params[2],params[3]

def plotTimeDomain(t,p1,p2,filename):
    figure(num=filename+' Time Domain')
    subplot(2,1,1)
    plot(t,p1,label='Pendulum 1')
    xlabel('Time (s)')
    ylabel('Angle (rad)')
    legend(loc=0)
    
    subplot(2,1,2)
    plot(t,p2,label='Pendulum 2')
    xlabel('Time (s)')
    ylabel('Angle (rad)')
    legend(loc=0)
    show()

def plotFreqDomain(freq,f1,f2,filename):
    figure(num=filename+' Frequency Domain')
    plot(freq,f1,label='Pendulum 1')
    plot(freq,f2,label='Pendulum 2')
    xlabel('Frequency (rad/s)')
    ylabel('Amplitude (arb. units)')
    legend(loc=0)
    xlim(-12,12)
    show()

def outputFinalData(filename, peaks_f1,peaks_f2,beat,beatError,expectedBeat,expectedBeatError):
    pass

# spring is 1,2,3,4, same with position
def analyzeData(filename, spring,position):
    t,p1,p2 = readData(filename)
    freq,f1,f2 = fourierTransform(t,p1,p2)
    
    expectedOmega_0, expectedOmega_1 = expectedFrequency(k[spring],l[position],2*M,L)
    expectedBeat = expectedOmega_1-expectedOmega_0
    #expectedBeatError = errorAddition(expectedFrequencyError(k[spring],l[position],M,L))
    
    peaks_f1 = findPeaks(freq,f1, [expectedOmega_0,expectedOmega_1])
    peaks_f2 = findPeaks(freq,f2, [expectedOmega_0,expectedOmega_1])
    
    f_error = 0.5*(freq[1]-freq[0])
    
    beatFreq_1 = abs(peaks_f1[1]-peaks_f1[0])
    beatFreq_2 = abs(peaks_f2[1]-peaks_f2[0])
    
    beatFreq = 0.5*(beatFreq_1+beatFreq_2)
    beatFreqError = f_error
    
    plotTimeDomain(t,p1,p2,filename)
    savefig(filename.replace('.txt','_time.png'))
    plotFreqDomain(freq,f1,f2,filename)
    savefig(filename.replace('.txt','_freq.png'))

    return peaks_f1,peaks_f2, f_error, beatFreq,beatFreqError

#analyzeData('partc_s1_l4.txt',1,4)
#analyzeData('partc_s2_l1.txt',2,1)
#analyzeData('partc_s3_l1.txt',3,1)
#analyzeData('partc_s3_l3.txt',3,3)


print analyzeData('partc_s2_l3.txt',2,3), 's2l3'
print analyzeData('partc_s3_l2.txt',3,2), 's3l2'
print analyzeData('partc_s3_l3_nxg.txt',3,3), 's3l3_nxg'
print analyzeData('partc_s4_l2.txt',4,2), 's4l2'


f_0Error = 0.009
f_0 = 4.147
f_23,f_23Error = 4.82,0.02
f_32,f_32Error = 5.86,0.03
f_33,f_33Error = 7.07,0.04
f_42,f_42Error = 4.79,0.01

beat_23,beat_23Error = f_23-f_0, errorAddition([f_23Error,f_0Error])
beat_32,beat_32Error = f_32-f_0, errorAddition([f_32Error,f_0Error])
beat_33,beat_33Error = f_33-f_0, errorAddition([f_33Error,f_0Error])
beat_42,beat_42Error = f_42-f_0, errorAddition([f_42Error,f_0Error])

rows = [[2,3,0.64475255856726132, beat_23],
        [3,2,1.4327185498282597, beat_32],
        [3,3,2.4774796324241315, beat_33],
        [4,2,0.58619399836624275, beat_42]]
errors = [[None, None, 0.058768292832781593, beat_23Error],
          [None, None, 0.058152918038720657, beat_32Error],
          [None, None, 0.053144532753088024, beat_33Error],
          [None, None, 0.034440953008981, beat_42Error]]
units = [['' for j in i] for i in rows]

f_1 = [4.7892788764101955, 5.6180963214705448, 6.6742552052830977, 4.7645321364439628]
f_2 = [4.8433631674458191, 5.6231711263560928, 6.6619674824091426, 4.7641814758092673]
f_error = [0.083110916761633291, 0.082240645381929767, 0.075157718985398247, 0.048706862846355392]
f_1_theory = [f_23, f_32, f_33, f_42]
f_1_error = [f_23Error, f_32Error, f_33Error, f_42Error]
k_error = [0.6, 1, 1, 0.2]
l_error = [0.3*0.01, 0.3*0.01, 0.3*0.01, 0.3*0.01]
df_1 = [f_23-f_1[0], -f_1[1]+f_32, -f_1[2]+f_33, -f_1[3]+f_42]
kl2 = [k[2]*l[3]**2, k[3]*l[2]**2, k[3]*l[3]**2, k[4]*l[2]**2]

f0 = [4.1823990104605748, 4.160737916260917, 4.1840476904604111, 4.191782657709707, 4.1868734172599451, 4.1943900055840322, 4.1773682704003949, 4.1789573451203497]

def linearFit(xData, a,b):
    return map(lambda x: a*x+b, xData)

def plotDeviations(kl2, K,L, f_1,f_2,f_t, f_error,f_t_error, k_error,l_error):
    kl2_error = [errorMultiplication(kl2[i], [K[i],L[i]**2],[k_error[i],errorPower(L[i]**2,L[i],l_error[i],2)]) for i in range(len(kl2))]
    avg_f = [0.5*(f_1[i]+f_2[i]) for i in range(len(f_1))]
    avg_error = map(lambda x: x/2**0.5, f_error)
    
    deviations = [100*abs(f_t[i]-avg_f[i])/f_t[i] for i in range(len(avg_f))]
    deviations_error = [errorMultiplication(deviations[i], [abs(f_t[i]-avg_f[i]),f_t[i]], [errorAddition([f_t_error[i], avg_error[i]]),f_t_error[i]]) for i in range(len(avg_error))]
    
    params,cov = curve_fit(linearFit, kl2,deviations, sigma=deviations_error)
    
    figure()
    errorbar(kl2, deviations, xerr=kl2_error, yerr=deviations_error, fmt='.')
    plot(xlim(),linearFit(xlim(),params[0],params[1]))
    xlabel('Coupling Strength, $kl^{2}$, (Nm)')
    ylabel('Percentage Deviation from Theory')
    show()
    
    return params, sqrt(diag(cov))