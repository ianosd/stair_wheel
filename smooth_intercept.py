import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from collections import namedtuple 


WheelGeometry = namedtuple('WheelGeometry', ['r', 'rp', 'R', 'Dtheta', 'curve_length', 'segm_count'])


def make_geometry(r, Dtheta, segm_count):
    assert 1 == r.ndim
    if r.dtype is np.dtype(int):
        r = np.array(r, dtype=float)


    theta_step = Dtheta/(r.size - 1)
    rp = (r[1:] - r[0:-1])/theta_step

    R = np.sqrt(r[0:-1]**2 + rp**2)

    curve_length= np.sum(np.sqrt(r[0:-1]**2 + r[1:]**2 - 2*np.cos(theta_step)*r[0:-1]*r[1:]))

    return WheelGeometry(r, rp, R, Dtheta, curve_length, segm_count)


def compute_constraints(g, stair_width, stair_height):
    c0 = g.curve_length - g.rp[-1]*g.r[-1]/g.R[-1] + g.rp[0]*g.r[0]/g.R[0] - stair_width
    c1 = g.r[-1]**2/g.R[-1] - g.r[0]**2/g.R[0] - stair_height

    v0 = np.array([-g.rp[-1], g.r[-1]])
    v1 = np.array([-g.rp[0], g.r[0]])

    dot = v0.dot(v1)
    det = v0[0]*v1[1] - v1[0]*v0[1]

    Dalpha = math.atan2(det, dot)

    c2 = g.Dtheta + Dalpha - 2*np.pi/g.segm_count

    return np.array([c0, c1, c2])


def plot_wheel(g, *args):
    if not args:
        args = ('b',)

    thetas_one_segm = np.linspace(0, g.Dtheta, g.r.size)

    for i in range(g.segm_count):
        thetas = thetas_one_segm + i*2*np.pi/g.segm_count
        x = g.r*np.cos(thetas)
        y = g.r*np.sin(thetas)

        plt.plot(x, y, *args)
        plt.plot([0, x[0]], [0, y[0]], *args)
        plt.plot([0, x[-1]], [0, y[-1]], *args)

stair_width = 1
stair_height = 0.7
segm_count = 3


def constr(arr):
    return compute_constraints(make_geometry(arr[:-1], arr[-1], segm_count), stair_width, stair_height)

def objective(arr):
    g = make_geometry(arr[:-1], arr[-1], segm_count)
    r = arr[:-1]
    segments = np.sqrt(r[0:-1]**2 + r[1:]**2 - 2*np.cos(theta_step)*r[0:-1]*r[1:])

    segments2b2 = np.sqrt(r[0::2]
    return  np.sum(g.rp**2)

def curve_length(arr):
    g = make_geometry(arr[:-1], arr[-1], segm_count)
    return stair_width - g.curve_length

r0 = np.linspace(1, 1, 10)
Dtheta0 = 1
x0 = np.append(r0, Dtheta0)

cons = {'type': 'eq', 'fun': constr}
cons_length = {'type': 'ineq', 'fun': curve_length}

bounds_constr = optimize.Bounds(0, np.append(np.repeat(100, r0.size), 2*np.pi/segm_count))
x_star = optimize.minimize(objective, x0, bounds=bounds_constr, constraints=(cons, cons_length)).x

# print(constr(x_star))

g_0 = make_geometry(r0, Dtheta0, segm_count)
g_star = make_geometry(x_star[:-1], x_star[-1], segm_count)

# print(g_star)

# plot_wheel(g_0, 'b')
plot_wheel(g_star, 'r')

def plot_stairs():
    x = -1
    y = -4

    for i in range(3):
        plt.plot([x, x+stair_width], [y, y], 'b')
        x += stair_width
        plt.plot([x, x], [y, y+stair_height], 'b')
        y += stair_height

plot_stairs()
plt.axis('equal')
plt.show()


