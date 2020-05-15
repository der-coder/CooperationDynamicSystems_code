#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 12:43:08 2020

@author: Isaac Ayala, CINVESTAV RYMA
"""


from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate
from math import exp

def solve_rossler(alpha = 1.0, beta = 1.0, gamma = 1.0, max_time = 10, N = 5):
    """Plot a solution to the Lorenz differential equations."""
    resolution = 72
    
    fig1 = plt.figure(dpi = resolution)
    ax1, ax2, ax3 = fig1.subplots(3, 1)

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    # ax.axis('off')
   
    def rossler_deriv(xi, t0, alpha=alpha, beta = beta, gamma = gamma):
        """Compute the time-derivative of a Lorenz system."""
        
        xi1, xi2, xi3 = xi
        
        xi1_dot = -xi1 - exp(xi3)
        xi2_dot = xi1 + alpha * xi2
        xi3_dot = gamma * exp(-xi3) + xi1 - beta
        
        return [xi1_dot, xi2_dot, xi3_dot]

    # Choose random starting points, uniformly distributed from -15 to 15
    np.random.seed(1)
    xi0 = -15 + 30 * np.random.random((N, 3))
    xi0[:,2] = np.log(abs(xi0[:,2]))

    # Solve for the trajectories
    t = np.linspace(0, max_time, int(250*max_time))
    xi_t = np.asarray([integrate.odeint(rossler_deriv, xi0i, t)
                      for xi0i in xi0])
    
    # choose a different color for each trajectory
    colors = plt.cm.viridis(np.linspace(0, 1, N))

    for i in range(N):
        xi1, xi2, xi3= xi_t[i,:,:].T
        
        ax1.plot(t, xi1, '-', c=colors[i])
        ax1.set_title('xi1')
        
        ax2.plot(t, xi2, '-', c=colors[i])
        ax2.set_title('xi2')
        
        ax3.plot(t, xi3, '-', c=colors[i])
        ax3.set_title('xi3')
        
        lines = ax.plot(xi1, xi2, xi3, '-', c=colors[i])
        plt.setp(lines, linewidth=2)
    angle = 104
    ax.view_init(30, angle)
    plt.show()

    return t, xi_t

if __name__ == '__main__':
    solve_rossler()