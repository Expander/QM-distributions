#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.ticker as tck
mpl.use('PDF')
import matplotlib.pyplot as plt

mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', weight='normal')
# mpl.rc('font', size=15)
mpl.rcParams.update({'text.usetex':True, 'text.latex.preamble':[r'\usepackage{amsmath}']})

h = 1.0
hquer = h / (2.*np.pi)
L = 1.0

def mom(n, sgn=1):
    return sgn*1.0*n*h/(2.0*L)

p1 = mom(1)

def mom_prefactor(p,n):
    return 4.0*L*n**2 / (hquer*np.pi**3) / ((p/p1)**2 - n**2)**2

def mom_even(p,n):
    return mom_prefactor(p,n) * np.sin(np.pi/2 * p/p1)**2

def mom_odd(p,n):
    return mom_prefactor(p,n) * np.cos(np.pi/2 * p/p1)**2

def mom_distribution(p,n):
    if np.mod(n,2) == 0:
        return mom_even(p,n)
    else:
        return mom_odd(p,n)

def plot_distribution(n, (lo,hi), filename):
    fig = plt.figure(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    ax = plt.gca()
    ax.get_yaxis().set_tick_params(which='both', direction='in', right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in', top=True)

    unit = hquer/L
    n_points = 1000
    px = np.linspace(lo*unit, hi*unit, n_points)

    plt.plot(px/unit, [mom_distribution(p,n)*unit for p in px], 'k-')
    plt.xlabel(r'$p\cdot L/\hbar$')
    plt.ylabel(r'$w(p)\cdot L/\hbar$')
    plt.title(r'Impulsverteilung im Potentialtopf $(n=' + str(n) +  r')$')

    ax.set_xlim([lo,hi])
    ax.set_ylim(bottom=0)

    plt_width  = ax.get_xlim()[1] - ax.get_xlim()[0]
    plt_height = ax.get_ylim()[1] - ax.get_ylim()[0]

    plt.annotate('', xy = (mom(n,-1)/unit, -0.15*plt_height),
                 xytext = (mom(n,-1)/unit, -0.01*plt_height),
                 arrowprops=dict(arrowstyle='<|-', facecolor='black'), annotation_clip=False) 

    plt.annotate('', xy = (mom(n,+1)/unit, -0.15*plt_height),
                 xytext = (mom(n,+1)/unit, -0.01*plt_height),
                 arrowprops=dict(arrowstyle='<|-', facecolor='black'), annotation_clip=False) 

    plt.savefig(filename)
    print "saved plot in ", filename
    plt.close(fig)

plot_distribution(1  , (-20 ,20) , 'Impulsverteilungen_Potentialtopf_n-1.pdf')
plot_distribution(5  , (-40 ,40) , 'Impulsverteilungen_Potentialtopf_n-5.pdf')
plot_distribution(100, (-400,400), 'Impulsverteilungen_Potentialtopf_n-100.pdf')
