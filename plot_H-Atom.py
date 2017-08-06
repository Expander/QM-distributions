#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.ticker as tck
mpl.use('PDF')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', weight='normal')
# mpl.rc('font', size=15)
mpl.rcParams.update({'text.usetex':True, 'text.latex.preamble':[r'\usepackage{amsmath}']})

h = 1.0
hquer = h / (2.*np.pi)
L = 1.0

def Y(theta,phi,l,m):
    if l == 0:
        return 1./(2.*np.sqrt(np.pi))
    if l == 1:
        return 0.5 * np.sqrt(3./np.pi)*np.cos(theta)

def F(p,n,l):
    # Bethe/Salpeter Eq.(8.10)
    if n == 1 and l == 0:
        return 4.*np.sqrt(2./np.pi) / (p**2 + 1)**2
    if n == 2 and l == 0:
        return 32./np.sqrt(np.pi) * (4.*p**2 - 1) / (4.*p**2 + 1)**3
    if n == 2 and l == 1:
        return 128./np.sqrt(3.*np.pi) * p / (4.*p**2 + 1)**3

def psi_p(p,theta,n,l,m):
    # Bethe/Salpeter Eq.(8.4)
    return F(p,n,l) * Y(theta,0,l,m)

def distr(p,theta,n,l,m):
    if l == 0:
        return 4*np.pi * p**2 * psi_p(p,theta,n,l,m)**2
    if l == 1:
        return 2*np.pi * p**2 * np.sin(theta) * psi_p(p,theta,n,l,m)**2

def plot_distribution_1D(n,l,m, (lo,hi), filename):
    fig = plt.figure(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    ax = plt.gca()
    ax.get_yaxis().set_tick_params(which='both', direction='in', right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in', top=True)

    unit = hquer/L
    n_points = 1000
    px = np.linspace(lo, hi, n_points)

    plt.plot(px, [distr(p,0,n,l,m) for p in px], 'k-')
    plt.xlabel(r'$p\cdot a/\hbar$')
    plt.ylabel(r'$w_{' + str(n) + str(l) + str(m) + r'}(p)$')
    plt.title(r'Impulsverteilung im H-Atom $(n=' + str(n) +  r',l=' + str(l) +  r',m=' + str(m) +  r')$', fontsize=11)

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

def plot_distribution_2D(n,l,m, (lo,hi), filename):
    fig = plt.figure(figsize=(5,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    ax = plt.gca()
    ax.get_yaxis().set_tick_params(which='both', direction='in', right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in', top=True)

    unit = hquer/L
    n_points = 100
    x = np.linspace(lo, hi, n_points)
    y = np.linspace(0, np.pi, n_points)
    x,y = np.meshgrid(x,y)

    z = distr(x,y,n,l,m)

    cax = ax.imshow(z, origin='lower',
                    extent=[x.min(), x.max(), y.min(), y.max()],
                    aspect='auto', cmap=plt.get_cmap('Greys'))

    ax.set_yticks([0., .25*np.pi, .5*np.pi, .75*np.pi, np.pi])
    ax.set_yticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])

    ax.set_xlabel(r'$p\cdot a/\hbar$')
    ax.set_ylabel(r'$\theta$')
    plt.title(r'Impulsverteilung im H-Atom $(n=' + str(n) +  r',l=' + str(l) +  r',m=' + str(m) +  r')$')

    clb = plt.colorbar(cax)
    clb.ax.set_ylabel(r'$w_{' + str(n) + str(l) + str(m) + r'}(p,\theta)$', rotation=-90, labelpad=15)
    clb.ax.tick_params(axis='y', direction='in')

    ax.set_xlim([lo,hi])
    ax.set_ylim([0,np.pi])
    # ax.set_ylim(bottom=0)

    plt.savefig(filename) # bbox_inches='tight'
    print "saved plot in ", filename
    plt.close(fig)

plot_distribution_1D(1,0,0, (0,4) , 'Impulsverteilungen_H-Atom_n-1.pdf')
plot_distribution_1D(2,0,0, (0,2) , 'Impulsverteilungen_H-Atom_n-2_l-0.pdf')
plot_distribution_2D(2,1,0, (0,2) , 'Impulsverteilungen_H-Atom_n-2_l-1.pdf')
