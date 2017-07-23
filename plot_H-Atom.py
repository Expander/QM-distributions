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

def psi_p(p,n,l,m):
    if n == 1:
        return (2*np.sqrt(2.)/np.pi) * 1./(p**2 + 1)**2

def distr(p,n,l,m):
    return psi_p(p,n,l,m)**2

def plot_distribution(n,l,m, (lo,hi), filename):
    fig = plt.figure(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    ax = plt.gca()
    ax.get_yaxis().set_tick_params(which='both', direction='in', right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in', top=True)

    unit = hquer/L
    n_points = 1000
    px = np.linspace(lo*unit, hi*unit, n_points)

    plt.plot(px/unit, [distr(p,n,l,m) for p in px], 'k-')
    plt.xlabel(r'$p\cdot a/\hbar$')
    plt.ylabel(r'$w_{' + str(n) + str(l) + str(m) + r'}(p)$')
    plt.title(r'Impulsverteilung im H-Atom $(n=' + str(n) +  r',l=' + str(l) +  r',m=' + str(m) +  r')$')

    ax.set_xlim([lo,hi])
    ax.set_ylim(bottom=0)

    plt_width  = ax.get_xlim()[1] - ax.get_xlim()[0]
    plt_height = ax.get_ylim()[1] - ax.get_ylim()[0]

    # plt.annotate('', xy = (mom(n,-1)/unit, -0.15*plt_height),
    #              xytext = (mom(n,-1)/unit, -0.01*plt_height),
    #              arrowprops=dict(arrowstyle='<|-', facecolor='black'), annotation_clip=False) 

    # plt.annotate('', xy = (mom(n,+1)/unit, -0.15*plt_height),
    #              xytext = (mom(n,+1)/unit, -0.01*plt_height),
    #              arrowprops=dict(arrowstyle='<|-', facecolor='black'), annotation_clip=False) 

    plt.savefig(filename)
    print "saved plot in ", filename
    plt.close(fig)

plot_distribution(1,0,0, (-20 ,20) , 'Impulsverteilungen_H-Atom_n-1.pdf')
