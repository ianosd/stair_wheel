import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from collections import namedtuple 


WheelGeometry = namedtuple('WheelGeometry', ['r', 'rp', 'R', 'Dtheta', 'curve_length', 'segm_count', 'theta_step', 'chunk_lengths'])


def make_geometry(r, Dtheta, segm_count):
    assert 1 == r.ndim
    if r.dtype is np.dtype(int):
        r = np.array(r, dtype=float)


    theta_step = Dtheta/(r.size - 1)
    rp = (r[1:] - r[0:-1])/theta_step

    R = np.sqrt(r[0:-1]**2 + rp**2)

    chunk_lengths = np.sqrt(r[0:-1]**2 + r[1:]**2 - 2*np.cos(theta_step)*r[0:-1]*r[1:])
    curve_length = np.sum(chunk_lengths)

    return WheelGeometry(r, rp, R, Dtheta, curve_length, segm_count, theta_step, chunk_lengths)


def compute_constraints(g, stair_width, stair_height):
    alpha_0 = compute_backward_angle(g.r[0], g.r[1], g.theta_step)
    alpha_1 = compute_angle(g.r[-2], g.r[-1], g.theta_step)

    Dalpha = alpha_1 - alpha_0
    c0 = g.curve_length - np.cos(alpha_1)*g.r[-1] + np.cos(alpha_0)*g.r[0] - stair_width
    c1 = g.r[-1]*np.sin(alpha_1) - g.r[0]*np.sin(alpha_0) - stair_height

    c2 = g.Dtheta + Dalpha - 2*np.pi/g.segm_count

    return np.array([c0, c1, c2])


def plot_wheel(g, *args, **kwargs):
    if not args:
        args = ('b',)

    if 'transform' in kwargs:
        transform = kwargs['transform']
    else:
        transform = np.eye(3)

    xx, yy = [], []
    thetas_one_segm = np.linspace(0, g.Dtheta, g.r.size)
    for i in range(g.segm_count):
        thetas = thetas_one_segm + i*2*np.pi/g.segm_count
        x = g.r*np.cos(thetas)
        y = g.r*np.sin(thetas)
        
        joined = np.concatenate(([x], [y], np.ones((1, x.size))))
        transformed = transform.dot(joined)
    
        x = transformed[0, :]
        y = transformed[1, :]

        xx.append(x)
        yy.append(y)

        if i == 0:
            stair_args = ([x[0] - stair_width, x[0], x[0], x[0] + stair_width], [y[0]-stair_height, y[0] - stair_height, y[0], y[0]], 'b')

    x = np.concatenate(xx)
    x = np.append(x, x[0])
    y = np.concatenate(yy)
    y = np.append(y, y[0])
    plt.plot(*stair_args)
    plt.plot(x, y, *args)
    plt.plot(0, 0, 'r*')

stair_width = 1
stair_height = 0.5
segm_count = 3


def compute_angle(r0, r1, dtheta):
    r2 = np.sqrt(r0**2 + r1**2 - 2*np.cos(dtheta)*r0*r1)
    sin_angle = np.sin(dtheta) * r0/r2
    cos_angle = (r1**2 + r2**2 - r0**2)/(2*r1*r2)
    return np.arctan2(sin_angle, cos_angle)

def compute_backward_angle(r0, r1, dtheta):
    return dtheta + compute_angle(r0, r1, dtheta)

def constr(arr):
    return compute_constraints(make_geometry(arr[:-1], arr[-1], segm_count), stair_width, stair_height)

def centre_points(g):
    bwd_angles = compute_backward_angle(g.r[:-1], g.r[1:], g.theta_step)
    l = np.insert(np.cumsum(g.chunk_lengths), 0, 0)
    return np.array([g.r[:-1]*np.cos(bwd_angles) + l[:-1], g.r[:-1]*np.sin(bwd_angles)])

def objective(arr):
    g = make_geometry(arr[:-1], arr[-1], segm_count)
    r = arr[:-1]
    theta_step = g.Dtheta/(r.size - 1)
    segments = np.sqrt(r[0:-1]**2 + r[1:]**2 - 2*np.cos(theta_step)*r[0:-1]*r[1:])
    segments2b2 = np.sqrt(r[0:-2]**2 + r[2:]**2 - 2*np.cos(2*theta_step)*r[0:-2]*r[2:])
    angles = np.arccos((segments[0:-1]**2 + segments[1:]**2 - segments2b2**2)/(2*segments[0:-1]*segments[1:]))
    # return -sum(g.r**2)
    return max(np.pi - angles)

def curve_length(arr):
    g = make_geometry(arr[:-1], arr[-1], segm_count)
    return stair_width - g.curve_length

r0 = np.linspace(4, 9, 7)
Dtheta0 = 1
x0 = np.append(r0, Dtheta0)

cons = {'type': 'eq', 'fun': constr}
cons_length = {'type': 'ineq', 'fun': curve_length}

bounds_constr = optimize.Bounds(0.1, np.append(np.repeat(100, r0.size), 2*np.pi/segm_count))
x_star = optimize.minimize(objective, x0, bounds=bounds_constr, constraints=(cons, cons_length)).x

# print(constr(x_star))

g_0 = make_geometry(r0, Dtheta0, segm_count)
g_star = make_geometry(x_star[:-1], x_star[-1], segm_count)

def plot_wheel_and_stair(g, *args):
    if not args:
        args = ('b',)

    angle0 = - compute_backward_angle(g.r[0], g.r[1], g.theta_step)
    print("Angle0", angle0)

    t = np.array([[np.cos(angle0), -np.sin(angle0), 0], [np.sin(angle0), np.cos(angle0), 0], [0, 0, 1]])

    plot_wheel(g, transform=t, *args)
    path = centre_points(g)
    plt.plot(path[0, :], path[1, :])

def plot_stairs():
    x = -1
    y = -4

    for i in range(3):
        plt.plot([x, x+stair_width], [y, y], 'b')
        x += stair_width
        plt.plot([x, x], [y, y+stair_height], 'b')
        y += stair_height

plot_wheel_and_stair(g_star, 'y')
plot_stairs()
plt.axis('equal')
plt.show()

