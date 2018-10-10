from scipy import *
from scipy.optimize import leastsq
import matplotlib.pylab as plt
from matplotlib import colors
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np

from tkinter import *
from tkinter import ttk
import math
import time
import os

import scipy.optimize as scimin
import matplotlib.pyplot as mpl

from decimal import *

#sum the squares of the values of a list
def sq_sum(l):
    rt = 0
    for i in l:
        rt += i ** 2
    return rt

#given a list of individual peaks, return a single list containing a composite peak
def make_composite(l):
    longest_len = 0
    for i in l:
        if len(i) > longest_len:
            longest_len = len(i)
    c = [0] * longest_len
    for i in range(len(l)):
        for j in range(longest_len):
            try:
                c[j] += l[i][j]
            except:
                pass
    return c

# least square fit with a single peak
#note: y is an extraneous parameter. not used
def one_peak(p,x,y):
    A, mean, sigma, tau = p
    exp_ = (x-mean)/tau + sigma**2/(2*tau**2)
    erfc_ = 1/math.sqrt(2)*((x-mean)/sigma + sigma/tau)
    ret = []
    for i, j in zip(exp_, erfc_):
        #if the number is too big to calculate in a float, i.e. overflow, then we probably didn't want to use this fit anyway
        #replace the number with infinity instead
        try:
            ret.append(A * math.exp(i) * math.erfc(j))
        except:
            ret.append(float("inf"))
    return ret

def one_peak_residuals(p,x,y):
    return y-one_peak(p,x,y)

def two_peak_residuals(p,x,y):
    p1, p2 = p[:int(len(p)/2)], p[int(len(p)/2):]
    return y-one_peak(p1,x,y)-one_peak(p2,x,y)

def three_peak_residuals(p,x,y):
    p1, p2, p3 = p[:int(len(p)/3)], p[int(len(p)/3):2*int(len(p)/3)], p[2*int(len(p)/3):]
    return y-one_peak(p1,x,y)-one_peak(p2,x,y)-one_peak(p3,x,y)


#data = np.loadtxt("alpha_peak_to_fit.txt")
#data = np.loadtxt("Two_alpha_peaks.txt")
data = np.loadtxt("Three_alpha_peaks.txt")
datax = data[:,0]
datay = data[:,1]     

is_one_peak = True
is_two_peaks = False
is_three_peaks = False

#one peak fitting
p1=[400, 80, 7, 5] # initial parameters guess
p,cov,infodict,mesg,ier=scimin.leastsq(one_peak_residuals, p1, (datax, datay),full_output=True) #traditionnal least squares fit
p_one_peak = p[:]
sum_one_peak = sq_sum(one_peak_residuals(p, datax, datay))

#two peak fitting
p2 = [400,82,7,5,400,95,7,5] #guess for two peaks
p,cov,infodict,mesg,ier=scimin.leastsq(two_peak_residuals, p2, (datax, datay),full_output=True)
p_two_peak = p[:]
sum_two_peak = sq_sum(two_peak_residuals(p_two_peak, datax, datay))

if sum_two_peak < sum_one_peak:
    is_two_peaks = True

#three peak fitting
p3= [50,50,7,8,50,75,7,8,50,80,7,8] #guess for three peaks
p,cov,infodict,mesg,ier=scimin.leastsq(three_peak_residuals, p3, (datax, datay),full_output=True)
p_three_peak = p[:]
sum_three_peak = sq_sum(three_peak_residuals(p_three_peak, datax, datay))
   
print("sums: 1", sum_one_peak, "2", sum_two_peak, "3", sum_three_peak)
if sum_three_peak < sum_two_peak and sum_three_peak < sum_one_peak:
    is_one_peak = False
    is_two_peaks = False
    is_three_peaks = True

# plotting
ax=mpl.figure().add_subplot(1,1,1)
ax.plot(datax,datay,ls="",marker="x",color="blue",mew=2.0,label="Data")
plot_x = np.linspace(0,100,1000)

#plot one peak
if is_one_peak:
    print("one peak", p)
    ax.plot(plot_x,one_peak(p_one_peak, plot_x, []),color="red",label="1peak")
   
#plot two peaks
elif is_two_peaks:
    p1, p2 = p_two_peak[:int(len(p_two_peak)/2)], p_two_peak[int(len(p_two_peak)/2):]
    print("two peaks", p1, p2)
    l1=one_peak(p1, plot_x, [])
    l2=one_peak(p2, plot_x, [])
    c = make_composite([l1, l2])

    ax.plot(plot_x,l1,color="red",label="1 of 2")
    ax.plot(plot_x,l2,color="magenta",label="2 of 2")
    ax.plot(plot_x,c,color="black",label="composite")

elif is_three_peaks:
    p1, p2, p3 = p[:int(len(p)/3)], p[int(len(p)/3):2*int(len(p)/3)], p[2*int(len(p)/3):]
    print("three peaks", p1, p2, p3)
    l1=one_peak(p1, plot_x, [])
    l2=one_peak(p2, plot_x, [])
    l3=one_peak(p3, plot_x, [])
    c = make_composite([l1, l2, l3])

    ax.plot(plot_x,l1,color="red",label="1 of 3")
    ax.plot(plot_x,l2,color="magenta",label="2 of 3")
    ax.plot(plot_x,l3,color="cyan",label="3 of 3")
    ax.plot(plot_x,c,color="black",label="composite")
   
ax.legend(loc=2)
mpl.show()
