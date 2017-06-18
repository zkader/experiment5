from math import log10,floor

def sigfig(val,err):
    return round(val,-int(floor(log10(abs(err)))))

def latexSI(xval,xerr,unit):
    if log10(abs(xerr)) < 0:
        return '\\SI{'+str(xval)+'('+str(xerr)[-1]+')}{'+unit+'}'
    else:
        return '\\SI{'+str(int(xval))+'('+str(int(xerr))+')}{'+unit+'}'

def printLaxtexTable(fname,val_ind,err_ind,unit):
    file = open(fname,'r')
    if len(val_ind) != len(unit) and len(val_ind) != len(err_ind):
        return
    for line1 in file:
        cLine = line1.split()
        numcol = len(cLine)
        break
    arr = []
    for i in range(numcol):
        arr.append([])
    for line in file:
        currLine = map(float, line.split())
        for j in range(numcol):
            arr[j].append(currLine[j])
    file.close()
    for l in range(len(arr[0])):
        dstr = ''
        for k in range(len(arr)):
            if k in val_ind:
                for m in range(len(val_ind)):
                    if k == val_ind[m]:
                        dstr += latexSI(sigfig(arr[val_ind[m]][l],arr[err_ind[m]][l]),sigfig(arr[err_ind[m]][l],arr[err_ind[m]][l]),unit[m]) + ' & '
            elif k not in val_ind and k not in err_ind:
                dstr += latexSI(arr[k][l],0,"") + ' & '
        print dstr +' \\\\\n'
    return arr

#mapping of file
##bvalind = [3,4,6,8]
##bvalerr = [5,5,7,9]
##unitmap = ["rad/s","rad/s","N/m","m"]
##printLaxtexTable('partb_all_results.txt',bvalind,bvalerr,unitmap)

b1 = [0,2]
be = [1,3]
ut = ["?","rad^{2}/s^{2}"]
arr1 = printLaxtexTable('linresults_partb_s0_l1_m1.txt',b1,be,ut)
#print sigfig(0.036455,0.005662)
#print sigfig(0.005662,0.005662)
