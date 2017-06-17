from pylab import figure, plot, show, xlabel, ylabel, errorbar, xlim, ylim, savefig,semilogy,title,legend,clf
from scipy.optimize import curve_fit
from numpy import sqrt, diag, pi, sin, linspace, fft,absolute,amax,argmax,where,abs,log

def readData(filename):
    file = open(filename, 'r')

    t,pend1,pend2 = [],[],[]
    
    for line in file:
        currLine = map(float, line.split())
        t.append(currLine[0])
        pend1.append(currLine[1])
        pend2.append(currLine[2])
    file.close()
    
    return t,pend1,pend2

def readpositions():
    file = open('parta_pendulum.txt', 'r')
    file.readline()
    m,L_t,L_b,l_1,l_2,l_3,l_4,P = [],[],[],[],[],[],[],[]
    vals = []
    avgs = []
    for line in file:
        currLine = map(float, line.split())
        for i in range(1,9):
            if i > 1:
                vals.append(currLine[i]/100)
            else:
                vals.append(currLine[i])
    file.close()
    for j in range(0,8):
        avgs.append((vals[j]+vals[8+j])/2)
    avgs[2] = avgs[7]-abs(avgs[1]-avgs[2])/2
    for k in range(3,7):
        avgs[k] = avgs[7]-avgs[k]
    del avgs[1]
    return avgs

def FFTdata(filename):
    samp_t = 1/5000.
    t,p1,p2 = readData(filename)

    ferr = 1/(2*t[len(t)-1])
    Fk1 = fft.fft(p1,norm="ortho") #Fourier coeff
    nf = fft.fftfreq(len(t),samp_t) #natural freq
    Fk1 = fft.fftshift(Fk1) # 0 at center
    nf = fft.fftshift(nf) # 0 at center
    Fk2 = fft.fft(p2,norm="ortho") #Fourier coeff
    Fk2 = fft.fftshift(Fk2) # 0 at center
    return 2*pi*nf, Fk1, Fk2, 2*pi*ferr

def FindFreq(Fcoef,nfrq,s1,swid):
    maxF = 0
    maxfindex = 0
    frq = 0
    for i in range(Fcoef.shape[0]):
        if nfrq[i] < s1-swid or nfrq[i] > s1+swid:
            continue
        else:
            if maxF < Fcoef[i]:
                maxF = Fcoef[i]
                frq = nfrq[i]
                maxfindex = i
    return frq

def PlotData(sval,lval,mode):
    if mode%2 == 0:
        smode = 'e'
        st = 'Even'
    elif mode%2 != 0:
        smode = 'o'
        st = 'Odd'
    fname = 'partb_s' + str(sval) + '_l' + str(lval) + '_' + smode + '.txt'
    t,p1,p2 = readData(fname)
    clf()
    plot(t,p1,label='Pendulum 1')
    plot(t,p2,label='Pendulum 2')
    xlim([t[0],t[len(t)-1]])
    xlabel('Time $t$ [s]')
    ylabel('Signal Amplitude [V]')
    title('Signal of %s Mode with Spring # %s at Location # %s'%(st,sval,lval))
    legend()
    #show()
    savefig('PartBPlots/signal_'+fname.replace('.txt','.png'))
    return
    
def AnalyzeData(sval,lval,mode):
    consts = readpositions()
    kval = [3.20778210044,30.4874228538,165.561223371,55.4326207947]
    kerr = [0.0178610036282,0.647160877092,1.17279065989,0.248952490127]
    k = kval[sval-1]
    ke = kerr[sval-1]
    lv = consts[1+lval]
    m = consts[0]
    L = consts[1]
    if mode%2 == 0:
        smode = 'e'
        st = 'Even'
        expectf = sqrt(9.81/L)
    elif mode%2 != 0:
        smode = 'o'
        st = 'Odd'
        expectf = sqrt(9.81/L) + k*lv*lv/(sqrt(9.81/L)*m*L*L)
    fname = 'partb_s' + str(sval) + '_l' + str(lval) + '_' + smode + '.txt'
    angf,Fc1,Fc2,frqerr = FFTdata(fname)
    abFF1 = absolute(Fc1)/amax(absolute(Fc1))
    abFF2 = absolute(Fc2)/amax(absolute(Fc2))
    frq1 = FindFreq(abFF1,angf,expectf,2*pi*0.2)
    frq2 = FindFreq(abFF2,angf,expectf,2*pi*0.2)
    print 'spring #',sval,"location #",lval,"mode:",smode
    print "Angular Frequency for P1:",frq1,"+/-",frqerr
    print "Angular Frequency for P2:",frq2,"+/-",frqerr,'\n'
    clf()
    plot(angf,abFF1,label='Pendulum 1')
    plot(angf,abFF2,label='Pendulum 2')
    xlim([-2*pi,2*pi])
    xlabel('Angular Frequency $\omega$ [rad/s]')
    ylim([0.0001,1])
    ylabel('Weights [arb. units]')
    title('%s Mode with Spring # %s at Location # %s'%(st,sval,lval))
    legend()
    #show()
    savefig('PartBPlots/'+fname.replace('.txt','.png'))
    return frq1,frq2,frqerr,k,ke

def writefiles():
    outfile = file('partb_all_results.txt','w')
    outfile.write('spring\tlocation\tmode\tFreq_P1\tFreq_P2\tFreq_err\tspring_const\tspring_consterr\n')
    for i in range(1,5):
        for j in range(1,5):
            for k in range(0,2):
                fp1,fp2,fe,spr,spre = AnalyzeData(i,j,k)
                PlotData(i,j,k)
                outfile.write('\t'.join(map(str,[i,j,k,fp1,fp2,fe,spr,spre])) +'\n')
    outfile.close()
    return
writefiles()
#PlotData(1,1,1)
#AnalyzeData(1,1,1)            
