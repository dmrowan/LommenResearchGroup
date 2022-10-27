import numpy as np
import matplotlib.pyplot as plt


def gauss(x, a, m, s):
    return(a*np.exp(-(((x - m)/s)**2)/2))

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def power(x, a, b):
    return(a*(x**(-b)))

def power2(x, a, b, c, d):
    return(a*(x**(-b)) + c*(x**(-d)))

def lognormal(x, a, m, s):
    return(a*(1/x)*(1/(s*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-m)**2/(2*(s**2)))))

def convolve(y, n):
    template = y
    for n in range(n):
        convo = np.convolve(template, y, mode='full')
        template = convo
    return(template)

def convolve2(y1, y2, n):
    template = y2
    for n in range(n):
        convo = np.convolve(template, y1, mode='full')
        template = convo
    return(template)

def convolve1(y1, y2):
    convo = []
    template = y2
    for i in range(len(y1)):
        convo.append(np.sum(y1*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    return(convo.index(m))

def hist_to_curve(vals, binwidths):
    yvals, xlims = np.histogram(vals, bins = binwidths) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot
    return(xvals, yvals)

def fwhm(x, y):
    xind = np.where(y >= (max(y)/2))
    l = x[xind[0][0]]
    r = x[xind[0][-1]]
    return(r-l)
