#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.ticker as tck
mpl.use('PDF')
import matplotlib.pyplot as plt
import scipy.misc as scmi
import scipy.special as scsp

mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', weight='normal')
# mpl.rc('font', size=15)
mpl.rcParams.update({'text.usetex':True, 'text.latex.preamble':[r'\usepackage{amsmath}']})

h = 1.0
hquer = h / (2.*np.pi)
m = 1
w = 1

def H(p,n):
    return scsp.hermite(n)(p)

def mom_distribution(p,n):
    return 1./(2.**n * scmi.factorial(n) * np.sqrt(np.pi)) * np.exp(-p**2) * H(p,n)**2 

def plot_distribution(n, (lo,hi), filename):
    fig = plt.figure(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    ax = plt.gca()
    ax.get_yaxis().set_tick_params(which='both', direction='in', right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in', top=True)

    n_points = 1000
    px = np.linspace(lo, hi, n_points)

    plt.plot(px, [mom_distribution(p,n) for p in px], 'k-')
    plt.xlabel(r'$p / \sqrt{m\omega\hbar}$')
    plt.ylabel(r'$w_{' + str(n) + r'}(p)$')
    plt.title(r'Impulsverteilung im harmonischen Oszillator $(n=' + str(n) +  r')$')

    ax.set_xlim([lo,hi])
    ax.set_ylim(bottom=0)

    plt.savefig(filename)
    print "saved plot in ", filename
    plt.close(fig)

plot_distribution(0  , (-4, 4), 'Impulsverteilungen_harmonischer_Oszillator_n-0.pdf')
plot_distribution(1  , (-4, 4), 'Impulsverteilungen_harmonischer_Oszillator_n-1.pdf')
plot_distribution(5  , (-6, 6), 'Impulsverteilungen_harmonischer_Oszillator_n-5.pdf')
plot_distribution(10 , (-10, 10), 'Impulsverteilungen_harmonischer_Oszillator_n-10.pdf')
