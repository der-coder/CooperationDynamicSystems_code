# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:57:51 2020

@author: Isaac Ayala, CINVESTAV RYMA SALTILLO
"""


from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate

def solve_chua(alpha=15.6, Lambda = 25., m1 = -5./7, m2 = -3./7,  N = 5, max_time = 15.0):
    """Plot a solution to the Lorenz differential equations."""
    resolution = 1200
    
    fig1 = plt.figure(dpi = resolution)
    ax1, ax2, ax3 = fig1.subplots(3, 1)

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    ax.axis('off')
    
    def chua_deriv(x, t0, alpha=alpha, Lambda=Lambda, m1=m1, m2=m2):
        """Compute the time-derivative of a Lorenz system."""
        # x1, x2, x1_hat, x2_hat, e1, e2, v = x
        x1, x2, x3 = x
        phi = m1 * x1 + m2 * (abs(x1 + 1) - abs(x1 -1))
             
        x1_dot = alpha * ( -x1 + x2 - phi)
        x2_dot = x1 - x2 + x3
        x3_dot = - Lambda * x2
        
        return [x1_dot, x2_dot, x3_dot]

    # Choose random starting points, uniformly distributed from -15 to 15
    np.random.seed(1)
    x0 = -15 + 30 * np.random.random((N, 3))

    # Solve for the trajectories
    t = np.linspace(0, max_time, int(250*max_time))
    x_t = np.asarray([integrate.odeint(chua_deriv, x0i, t)
                      for x0i in x0])
    
    # choose a different color for each trajectory
    colors = plt.cm.viridis(np.linspace(0, 1, N))

    for i in range(N):
        x1, x2, x3= x_t[i,:,:].T
        
        ax1.plot(t, x1, '-', c=colors[i])
        ax1.set_title('x1')
        
        ax2.plot(t, x2, '-', c=colors[i])
        ax2.set_title('x2')
        
        ax3.plot(t, x3, '-', c=colors[i])
        ax3.set_title('x3')
        
        lines = ax.plot(x1, x2, x3, '-', c=colors[i])
        plt.setp(lines, linewidth=2)
    angle = 104
    ax.view_init(15, angle)
    plt.show()

    return t, x_t

if __name__ == '__main__':
    solve_chua()