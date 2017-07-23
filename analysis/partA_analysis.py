from pylab import plot, show, figure, xlabel, ylabel, legend, xlim, ylim, errorbar, savefig, subplot
from scipy.optimize import curve_fit
from numpy import sqrt, diag

def readDataSpring(filename):
    file = open(filename, 'r')
    line = file.readline().split()[5:]
    file.close()
    
    m,x,dm = [],[],[]
    for i in range(0,len(line),3):
        m.append(0.001*float(line[i]))
        x.append(0.01*float(line[i+1]))
        dm.append(0.001*float(line[i+2]))
    return m,x,dm

def plotSpringData(m,x,dm,dx, params):
    figure()
    errorbar(m,x, xerr=dm,yerr=dx, fmt='.')
    xlim(xmin=0.02)
    limits = xlim()
    plot(limits, linear(limits, params[0],params[1]))
    xlim(limits)
    ylabel('Displacement (m)')
    xlabel('Applied Mass (kg)')
    show()

def linear(xData, k,b):
    g = 9.81
    return map(lambda x: (g/k)*x+b, xData)

def analyzeSpring(filename):
    m,x,dm = readDataSpring(filename)
    dx = 0.0002
    
    params,cov = curve_fit(linear, m,x, sigma=dx)
    errors = sqrt(diag(cov))
    
    k = params[0]
    k_error = errors[0]
    
    plotSpringData(m,x,dm,dx, params)
    savefig(filename.replace('.txt','.png'))
    
    #return k,k_error
    return params, m,x,dm,dx, k,k_error

# filenames is a list of filenames for springs 1,2,3,4
def analyzeAllSprings(filenames):
    k = []
    k_error = []
    for filename in filenames:
        currK, currK_error = analyzeSpring(filename)
        k.append(currK)
        k_error.append(currK_error)
    
    file = open('parta_springs_results.txt', 'w')
    file.write('Spring\tk (N/m)\n')
    for s in range(len(k)):
        file.write('{0}\t{1}+/-{2}\n'.format(s+1,k[s],k_error[s]))
    file.close()
    
    return makeRowsFinal([i+1 for i in range(len(k))], k,k_error)

filenames = ['parta_s1.txt', 'parta_s2.txt',
             'parta_s3.txt', 'parta_s4.txt']
#analyzeAllSprings(filenames)


def makeRowsFinal(springs, k,k_error):
    rows = [[springs[i],k[i]] for i in range(len(springs))]
    errors = [[None,k_error[i]] for i in range(len(springs))]
    units = [['' for j in i] for i in rows]
    
    return rows,errors,units

def makeRows(filename):
    rows = [[]]
    error = [[]]
    
    m,x,dm = readDataSpring(filename)
    dx = 0.02
    
    for i in range(len(m)):
        rows.append([1000*m[i],100*x[i]])
        error.append([1000*dm[i],dx])
    
    units = [['' for j in i] for i in rows]
    
    return rows, error, units
        

def readDataPendula(filename):
    file = open(filename, 'r')
    file.readline()
    p1 = map(float, file.readline().split())
    n1, m1,L_t1,L_b1,l1_1,l2_1,l3_1,l4_1,l_p = p1
    L1_1 = l_p - 0.5*(L_t1+L_b1)
    l1_1,l2_1,l3_1,l4_1 = l_p-l1_1, l_p-l2_1, l_p-l3_1, l_p-l4_1
    
    p2 = map(float, file.readline().split())
    n2, m2,L_t2,L_b2,l1_2,l2_2,l3_2,l4_2,l_p = p1
    L1_2 = l_p - 0.5*(L_t2+L_b2)
    l1_2,l2_2,l3_2,l4_2 = l_p-l1_2, l_p-l2_2, l_p-l3_2, l_p-l4_2
    file.close()

params1,m1,x1,dm1,dx1,k1,dk1 = analyzeSpring(filenames[0])
params2,m2,x2,dm2,dx2,k2,dk2 = analyzeSpring(filenames[1])
params3,m3,x3,dm3,dx3,k3,dk3 = analyzeSpring(filenames[2])
params4,m4,x4,dm4,dx4,k4,dk4 = analyzeSpring(filenames[3])

figure()

subplot(211)
errorbar(m1,x1, xerr=dm1,yerr=dx1, fmt='.', label='Spring 1: $k=3.21\\pm0.02Nm^{-1}$')
errorbar(m2,x2, xerr=dm2,yerr=dx2, fmt='.', label='Spring 2: $k=30.5\\pm0.6Nm^{-1}$')
xlim(xmin=0.02)
limits = xlim()
plot(limits, linear(limits, params1[0],params1[1]))
plot(limits, linear(limits, params2[0],params2[1]))
xlim(limits)
ylabel('Displacement (m)')
xlabel('Applied Mass (kg)')
legend(loc=0)

subplot(212)
errorbar(m3,x3, xerr=dm3,yerr=dx3, fmt='.', label='Spring 3: $k=166\\pm1Nm^{-1}$')
errorbar(m4,x4, xerr=dm4,yerr=dx4, fmt='.', label='Spring 4: $k=55.4\\pm0.2Nm^{-1}$')
xlim(xmin=0.02)
limits = xlim()
plot(limits, linear(limits, params3[0],params3[1]))
plot(limits, linear(limits, params4[0],params4[1]))
xlim(limits)
ylabel('Displacement (m)')
xlabel('Applied Mass (kg)')
legend(loc=0)

show()