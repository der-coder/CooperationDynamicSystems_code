# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:57:51 2020

@author: Isaac Ayala, CINVESTAV RYMA SALTILLO
"""


from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate

def solve_pecora(alpha=1.0, beta=-1.0, k1 = 1.0,  N = 5, max_time = 10.0):
    """Plot a solution to the Lorenz differential equations."""
    k2 = 1 + alpha;
    resolution = 1200
    
    fig1 = plt.figure(dpi = resolution)
    ax1, ax2 = fig1.subplots(1, 2)
    
    fig2 = plt.figure(dpi = resolution)
    ax3, ax4 = fig2.subplots(1, 2)
    
    fig3 = plt.figure(dpi = resolution)
    ax5 = fig3.add_axes([0, 0, 1, 1], projection=None)
    
    def pecora_deriv(x, t0, alpha=alpha, beta=beta, k1=k1, k2=k2):
        """Compute the time-derivative of a Lorenz system."""
        x1, x2, x1_hat, x2_hat, e1, e2, v = x
             
        x1_dot = x2
        x2_dot = alpha * x1 + beta * x2
        
        x1h_dot = x2_hat + k1 * e1
        x2h_dot = alpha * x1_hat + beta * x2_hat + k2 * e1
        
        e1_dot = e2 - k1 * e1
        e2_dot = (alpha - k2) * e1 + beta * e2
        
        v_dot = 2 * (alpha - k2 + 1) * e1 * e2 + 2 * beta * (e2 ** 2) - 2 * k1 * (e1 ** 2)
        return [x1_dot, x2_dot, x1h_dot, x2h_dot, e1_dot, e2_dot, v_dot]

    # Choose random starting points, uniformly distributed from -15 to 15
    np.random.seed(1)
    x0 = -15 + 30 * np.random.random((N, 2))
    x0_h = -18 + 30 * np.random.random((N, 2))
    
    e0 = x0 - x0_h
    v0 = e0[:,0] ** 2 + e0[:,1] ** 2
    
    state0 = np.hstack((x0, x0_h, e0))
    state0 = np.c_[state0, v0]

    # Solve for the trajectories
    t = np.linspace(0, max_time, int(250*max_time))
    x_t = np.asarray([integrate.odeint(pecora_deriv, x0i, t)
                      for x0i in state0])
    
    # choose a different color for each trajectory
    colors = plt.cm.viridis(np.linspace(0, 1, N))

    for i in range(N):
        x1, x2, x1h, x2h, e1, e2, v = x_t[i,:,:].T
        
        ax1.plot(x1, x2, '-', c=colors[i])
        ax1.set_title('Original system')
        
        ax2.plot(x1h, x2h, '-', c=colors[i])
        ax2.set_title('Observer')
        
        ax3.plot(t, e1, '-', c=colors[i])
        ax3.set_title('Observer error 1')
        
        ax4.plot(t, e2, '-', c=colors[i])
        ax4.set_title('Observer error 2')
        
        ax5.plot(t, v, '-', c=colors[i])
        ax5.set_title('Lyapunov candidate function')
        
    plt.show()

    return t, x_t

if __name__ == '__main__':
    solve_pecora()