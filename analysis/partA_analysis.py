from pylab import plot, show, figure, xlabel, ylabel, legend, xlim, ylim, errorbar, savefig
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
    plot(xlim(), linear(xlim(), params[0],params[1]))
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
    
    return k,k_error

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

filenames = ['parta_s1.txt', 'parta_s2.txt',
             'parta_s3.txt', 'parta_s4.txt']
analyzeAllSprings(filenames)


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