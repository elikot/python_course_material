#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../code/')
from tools_2 import peakdetect
from tools_2 import savitzky_golay ## Smoothing algorithm ##
from scipy.interpolate import interp1d ##Interpolation ##
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #

############################ read in datafile:
name="C3A17.dat"
f = open('C3A17.dat', 'r')
#name="C3A18.dat"
#f = open('C3A18.dat', 'r')
datalines=f.readlines()
xd=[]
yd=[]
for el in datalines:
    d=el.split('\t')
    xd.append(float(d[0]))
    yd.append(float(d[1]))
xd=np.array(xd)
yd=np.array(yd)
#### smoothing #############ind_bg_low
yd_smooth=savitzky_golay(yd,55,3)



#### Background subtraction ######
# defining the 'background' part of the spectrum #

### Method 1:
y_bg=np.mean(yd[0:2])
background=[y_bg]*len(xd)

### Method 2:
#ind_bg_low = (xd >= np.min(xd)) &(xd <= np.min(xd)+1) 
#ind_bg_high = (xd >= np.max(xd)-1) & (xd <= np.max(xd)) 
#x_bg = np.concatenate((xd[ind_bg_low],xd[ind_bg_high])) # joins the ends
#y_bg = np.concatenate((yd[ind_bg_low],yd[ind_bg_high]))
### interpolating the background #
#f = interp1d(x_bg,y_bg, kind='linear')
#background = f(xd)

# removing fitted background #
y_bg_corr = yd_smooth - background # corrected intensity data


########################### finding maxima
x_max,x_min=peakdetect(y_bg_corr,xd, lookahead=60, delta=15)

###### Defining function ######################################################################################################################################
def tslfn(x,p):
      A=p[0]
      b=p[1]
      xc=p[2]
      w=p[3]
      k=11604.5
      return (A*b**(b/(b-1))*np.exp(k*w*(x-xc)/(x*xc)) *(1+(b-1)*(x**2)*(1-2*x/(k*w))*np.exp(k*w*(x-xc)/(x*xc))/(xc**2)+2*(b-1)*xc/(k*w))**(-b/(b-1)))

def within_bounds(p):
    if p[1]>0.8 and p[1]<2.2:
       return True
    else:
       return False
      

def residuals(p,y,x):
    if within_bounds(p):
       err=y-tslfn(x,p)
       return err
    else:
       print 'b is not within bounds'
       
       

################################# Fitting
#### initial values p=[A,b,xc,w]
p=[]
for i in range(1): #len(x_max)):kop
    p.append((x_max[i][1])) # A
    p.append(2.0) #b
    p.append((x_max[i][0])) # xc
    p.append(0.3) # w

#print len(p)
#plt.plot(xd,tslfn(xd,p),'g-')
#plt.plot(xd,residuals(p,y_bg_corr,xd), 'r-', label='residuals')

## optimization:
pbest,cov,infodict,mesg,ier=leastsq(residuals,p,args=(y_bg_corr,xd), full_output=1) 
#print pbest
best_parameters=pbest

## fit to data:
fit=tslfn(xd,best_parameters)
#print best_parameters


ssErr=(infodict['fvec']**2).sum()
ssTot=((y_bg_corr-y_bg_corr.mean())**2).sum()
rsquared=1-(ssErr/ssTot)
#print rsquared
#print y_bg_corr.mean()

#### Plotting ###########
plt.rcParams.update({'font.size':10})
plt.plot(xd, yd, 'k-', label='data, {0}'.format(name))
#plt.plot(xd, yd_smooth, 'k', label='smoothed, {0}'.format(name))
plt.plot(xd, background, 'g--', label='background')
plt.plot(xd, y_bg_corr, 'bo', label='smoothed - bg, {0}'.format(name))
plt.plot(xd,fit,'r.', label='fit')
plt.plot(xd,residuals(best_parameters,y_bg_corr,xd), 'r-', label='residuals')
plt.text(xd[0]*1.01, np.max(yd), r'$R^2={0}$ '.format(np.around(rsquared,3)))
plt.text(xd[0]*1.01, np.max(yd)-50, r'$b={0};  Ea={1} $'.format(np.around(best_parameters[1],3),np.around(best_parameters[3],decimals=4)))

plt.ylabel('Intensity')
plt.xlabel('Temperature, K')
plt.axis([xd[0], xd[-1], -25, np.max(yd)*1.1])
plt.legend(loc=0)
plt.show()

###############################################
