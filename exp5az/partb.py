from pylab import figure, plot, show, xlabel, ylabel, errorbar, xlim, ylim, savefig,semilogy,title,legend,clf
from scipy.optimize import curve_fit
from numpy import sqrt, diag, pi, sin, linspace, fft,absolute,amax,amin,argmax,where,abs,log

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

def adderr(x1):
    ans = 0
    for i in range(len(x1)):
        ans += x1[i]**2 
    return sqrt(ans)

def multerr(x1,x1e):
    ans = 0
    mult = 1
    for i in range(len(x1)):
        ans += (x1e[i]/x1[i])**2
        mult *= x1[i]
    return mult*ans

def readpositions():
    file = open('parta_pendulum.txt', 'r')
    file.readline()
    merr = 0.5/1000
    lerr = 0.3/100
    m,L_t,L_b,l_1,l_2,l_3,l_4,P = [],[],[],[],[],[],[],[]
    vals = []
    avgs = []
    averr = []
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
        if j == 0:
            averr.append(merr/sqrt(2))
        else:
            averr.append(lerr/sqrt(2))
    avgs[2] = avgs[7]-abs(avgs[1]-avgs[2])/2
    averr[2] = adderr([averr[7],multerr([abs(avgs[1]-avgs[2]),0.5],[adderr(averr[1:3]),0.0])])
    for k in range(3,7):
        avgs[k] = avgs[7]-avgs[k]
        averr[k] = adderr([averr[7],averr[k]])
    del avgs[1]
    del averr[1]
    return avgs,averr

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

def linFN(xarr,m,b):
    return map(lambda x: m*x+b,xarr)

def PlotData(sval,lval,mode,saveplot=0):
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
    if saveplot == 0:
        return
    elif saveplot == 1:
        savefig('PartBPlots/signal_'+fname.replace('.txt','.png'))
    else:
        show()
    return
    
def AnalyzeData(sval,lval,mode,saveplot=0):
    consts,conerr = readpositions()
    kval = [3.20778210044,30.4874228538,165.561223371,55.4326207947]
    kerr = [0.0178610036282,0.647160877092,1.17279065989,0.248952490127]
    k = kval[sval-1]
    ke = kerr[sval-1]
    lv = consts[1+lval]
    lverr = conerr[1+lval]
    m = consts[0]
    L = consts[1]
    if mode%2 == 0:
        smode = 'e'
        st = 'Even'
        expectf = sqrt(9.81/L)
    elif mode%2 != 0:
        smode = 'o'
        st = 'Odd'
        expectf = sqrt(9.81/L + 2*k*lv*lv/(m*L*L))
        expectf1 = sqrt(9.81/L) + k*lv*lv/(sqrt(9.81/L)*m*L*L)
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
    if abs(4*pi) < abs(expectf):
        xlim([-2*pi,2*pi])
    else:
        xlim([-4*pi,4*pi])
    xlabel('Angular Frequency $\omega$ [rad/s]')
    ylim([0.0001,1])
    ylabel('Weights [arb. units]')
    legend()
    if saveplot == 0:
        return frq1,frq2,frqerr,k,ke,lv,lverr
    elif saveplot == 1:
        savefig('PartBPlots/'+fname.replace('.txt','.png'))
    else:
        show()
    return frq1,frq2,frqerr,k,ke,lv,lverr

def writefiles(createplots=0):
    outfile = file('partb_all_results.txt','w')
    outfile.write('spring\tloc\tmode\tFreq_P1\tFreq_P2\tFreq_err\tspring_const\tspring_consterr\tloc_val\tloc_err\n')    
    for i in range(1,5):
        for j in range(1,5):
            for k in range(0,2):
                fp1,fp2,fe,spr,spre,le,lerr = AnalyzeData(i,j,k,createplots)
                PlotData(i,j,k,createplots)
                outfile.write('\t'.join(map(str,[i,j,k,fp1,fp2,fe,spr,spre,le,lerr])) +'\n')
    outfile.close()
    return

def ReadResults():
    bigarr = []
    file = open('partb_all_results.txt','r')
    file.readline()
    for line in file:
        currLine = map(float, line.split())
        bigarr.append(currLine)
    file.close()
    return bigarr

def rearrange(arr):
    sm1 = []
    for n in range(4):
        dm1 = []
        for ni in range(4):
            dm1.append(arr[4*ni+n])
        sm1 += dm1
    return sm1

def FitLine(xar,yar,year):
    pr1,cv1 = curve_fit(linFN,xar,yar,sigma=year)
    return pr1,cv1

def PlotBFigs(location=0,spring=0,mode=1,sfig=0):
    bigarray = ReadResults()
    xval = []
    yval = []
    xerr = []
    yerr = []
    ytl = 'Angular Frequency Squared $\omega^{2}$ [$rad^{2}/s^{2}$]'
    index = []
    for i in range(len(bigarray)):
        if location !=0 and spring !=0:
            return
        if spring == 0 and bigarray[i][2] == mode:
            xtl = 'Squared Distance from Point $l^{2}$ [$m^{2}$]'
            ttl = ''
            for j in range(1,5):
                if j == bigarray[i][1]:
                    index.append(j)
                    xval.append(bigarray[i][8]**2)
                    xerr.append(multerr([bigarray[i][8],bigarray[i][8]],[bigarray[i][9],bigarray[i][9]]))                    
                    yval.append(bigarray[i][3]**2)
                    yerr.append(multerr([bigarray[i][3],bigarray[i][3]],[bigarray[i][5],bigarray[i][5]]))
        if location == 0 and bigarray[i][2] == mode:
            xtl = 'Spring Constant $k$ [N/m]'
            for k in range(1,5):
                if k == bigarray[i][0]:
                    index.append(k)
                    xval.append(bigarray[i][6])
                    xerr.append(bigarray[i][7])
                    yval.append(bigarray[i][3]**2)
                    yerr.append(multerr([bigarray[i][3],bigarray[i][3]],[bigarray[i][5],bigarray[i][5]]))
    clf()  
    if location == 0:
        lbl = 'Location '
        xval = rearrange(xval)
        xerr = rearrange(xerr)
        yval = rearrange(yval)
        yerr = rearrange(yerr)
        index = rearrange(index)
    if spring == 0:
        lbl = 'Spring '
    slope = []
    slopeerr = []
    yint = []
    yinterr = []
    xsss = linspace(amin(xval),amax(xval),10)
    for i in range(4):
        errorbar(xval[4*i:4*i+4],yval[4*i:4*i+4],xerr=xerr[4*i:4*i+4],yerr=yerr[4*i:4*i+4],fmt=".",label= lbl+'%s'%(str(i+1)))
        par,cov = FitLine(xval[4*i:4*i+4],yval[4*i:4*i+4],yerr[4*i:4*i+4])
        print cov
        slope.append(par[0])
        slopeerr.append(sqrt(diag(cov)[0]))
        yint.append(par[1])
        yinterr.append(sqrt(diag(cov)[1]))
        plot(xsss,linFN(xsss,slope[i],yint[i]),"k")
    print sqrt(yint)
    xlabel(xtl)
    ylabel(ytl)
    legend()
    if sfig == 0:
        show()
    if sfig == 1:
        savefig('PartBPlots/linplots_partb_s%s_l%s_m%s.png'%(str(spring),str(location),str(mode)))
    if sfig == 2:
        savefig('PartBPlots/linplots_partb_s%s_l%s_m%s.png'%(str(spring),str(location),str(mode)))
        outfile = file('linresults_partb_s%s_l%s_m%s.txt'%(str(spring),str(location),str(mode)),'w')
        outfile.write('slope\tsloper_err\tyint\tyint\n')
        for j in range(4):
            outfile.write('\t'.join(map(str,[slope[j],slopeerr[j],yint[j],yinterr[j]])) +'\n')
        outfile.close()
    return

#PlotBFigs(1,0)
#readpositions()
#writefiles()
PlotData(1,1,1,2)
#AnalyzeData(1,1,1)            
