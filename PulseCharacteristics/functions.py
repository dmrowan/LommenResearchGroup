import numpy as np
import matplotlib.pyplot as plt

"""
Various functions used for data analysis
"""

def gauss(x, a, m, s):  #single peak Gaussian, no shift, do not use if need shift (fitting to pulse profiles)
    return(a*np.exp(-(((x - m)/s)**2)/2)) #a=amplitude; m=mean; s=standard deviation/width

def gauss2(x, a, m, s, b, c, e, d): #double peak Gaussian with a shift, use when fitting to both main pulse and interpulse
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def power(x, a, b): #power law
    return(a*(x**(-b))) #a=scaling factor; b=power

def power2(x, a, b, c, d): #sum of two power laws
    return(a*(x**(-b)) + c*(x**(-d))) 

def lognormal(x, a, m, s): #log normal; a=scaling factor (NOT amplitude); m=mean and s=standard deviation but not really
    #just go to the log normal wikipedia page, this is the PDF log normal distribution for description of parameters
    return(a*(1/x)*(1/(s*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-m)**2/(2*(s**2)))))

def convolve(y, n): #convolve function with itself n times using numpy convolve; use where you need to do a "real" convolution
    template = y
    for n in range(n):
        convo = np.convolve(template, y, mode='full')
        template = convo
    return(template)

def convolve2(y1, y2, n): #convolve two different functions n times; template is different from the array it is convolved with 
    template = y2
    for n in range(n):
        convo = np.convolve(template, y1, mode='full')
        template = convo
    return(template)

def convolve1(y1, y2): #the "convolution" used for finding the location of the peak value of y1; returns the location of the max value of convolution
    convo = []
    template = y2
    for i in range(len(y1)):
        convo.append(np.sum(y1*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    return(convo.index(m))

def hist_to_curve(vals, binwidths): #converts a histogram to a curve
    yvals, xlims = np.histogram(vals, bins = binwidths) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot
    return(xvals, yvals)

def fwhm(x, y): #calculates the full width half max of a curve
    xind = np.where(y >= (max(y)/2))
    l = x[xind[0][0]]
    r = x[xind[0][-1]]
    return(r-l)
