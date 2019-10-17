import numpy as np
import scipy.optimize as optimize


class StarWheel:
    def __init__(self, r, theta, n):
        self.r = r
        self.theta = theta
        self.n = n
        
        self._cache = {}
        
    # We define some additional useful properties
    
    @property
    def theta_step(self):
        return self.theta/(self.r.size - 1)
    
    @property
    def rp(self):
        return (self.r[1:] - self.r[0:-1])/self.theta_step
    
    @property
    def curve_length(self):
        return np.sum(self.chunk_lengths)
    
    @property
    def chunk_lengths(self):
        return np.sqrt(self.r[0:-1]**2 + self.r[1:]**2 - 2*np.cos(self.theta_step)*self.r[0:-1]*self.r[1:])

    @property
    def radii_and_theta_array(self):
        return np.append(self.r, self.theta)
    
    def get_transform(self, s):
        floor_s = np.floor(s)
        
        dfirst = self._centre_trajectory_data[:, 0]
        dlast = self._centre_trajectory_data[:, -1]
        delta_translation = floor_s*(dlast[0:2] - dfirst[0:2])
        delta_rotation = floor_s*(dlast[2] - dfirst[2])
        
        num_transforms = self._centre_trajectory_data.shape[1]
        transform_index = np.floor((s - floor_s) * num_transforms)
        transform_index = int(min(num_transforms-1, transform_index))
        
        return transform2d(delta_rotation + self._centre_trajectory_data[2, transform_index], delta_translation + self._centre_trajectory_data[0:2, transform_index])
        
    @property
    def _centre_trajectory_data(self):
        try:
            return self._cache['centre_trajectory_data']
        except KeyError:
            pass
        
        g = self
        bwd_angles = compute_backward_angle(g.r[:-1], g.r[1:], g.theta_step)
        fwd_angles = compute_angle(g.r[:-1], g.r[1:], g.theta_step)

        steps = 10
        alphas = np.linspace(0, 1, steps)
        result = np.empty((3, 0), dtype=g.r.dtype)
        l = np.cumsum(g.chunk_lengths)
        rotation_angles = -(
            bwd_angles
            + g.theta_step*np.arange(g.r.size-1)
        )
        
        for i in range(g.r.size-2):
            fwd_angle = fwd_angles[i]
            bwd_angle = bwd_angles[i+1]
            
            pos_angles = np.pi - lerp(fwd_angle, bwd_angle, alphas)
            rotation_angles_interp = lerp(rotation_angles[i], rotation_angles[i+1], alphas)
            data = np.array([g.r[i+1]*np.cos(pos_angles) + l[i], g.r[i+1]*np.sin(pos_angles), rotation_angles_interp])
            result = np.concatenate((result, data), 1)
            
        self._cache['centre_trajectory_data'] = result
        return result
    
    @property
    def centre_trajectory(self):
        return self._centre_trajectory_data[0:2, :]
    
    @property
    def centre_trajectory_length(self):
        trajectory_segments = self.centre_trajectory[0:2, 1:] - self.centre_trajectory[:, :-1]
        return np.sum(np.sqrt(np.sum(trajectory_segments**2, 1)))

    
def transform2d(angle, translation):
    return np.array([
        [np.cos(angle), -np.sin(angle), translation[0]],
        [np.sin(angle), np.cos(angle), translation[1]],
        [0, 0, 1]
    ])
    

def compute_angle(r0, r1, dtheta):
    r2 = np.sqrt(r0**2 + r1**2 - 2*np.cos(dtheta)*r0*r1)
    sin_angle = np.sin(dtheta) * r0/r2
    cos_angle = (r1**2 + r2**2 - r0**2)/(2*r1*r2)
    return np.arctan2(sin_angle, cos_angle)


def compute_backward_angle(r0, r1, dtheta):
    return dtheta + compute_angle(r0, r1, dtheta)


def lerp(a, b, s):
    return a + (b-a)*s


def array_to_wheel_geometry(array, n):
    return StarWheel(array[:-1], array[-1], n)


def decorate_to_take_array(func, n, *args, **kwargs):
    def decorated_func(array):
        return func(array_to_wheel_geometry(array, n), *args, **kwargs)
    return decorated_func


def zero_constraints(g, stair_width, stair_height):
    alpha_0 = compute_backward_angle(g.r[0], g.r[1], g.theta_step)
    alpha_1 = compute_angle(g.r[-2], g.r[-1], g.theta_step)

    Dalpha = alpha_1 - alpha_0
    c0 = g.curve_length - np.cos(alpha_1)*g.r[-1] + np.cos(alpha_0)*g.r[0] - stair_width
    c1 = g.r[-1]*np.sin(alpha_1) - g.r[0]*np.sin(alpha_0) - stair_height

    c2 = g.theta + Dalpha - 2*np.pi/g.n

    return np.array([c0, c1, c2])


def inequality_constraints(g, stair_width):
    curve_length_constraint = stair_width*0.8 - g.curve_length
    
    c = np.cos(g.theta_step)
    s = np.sin(g.theta_step)
    c2 = np.cos(2*g.theta_step)
    s2 = np.sin(2*g.theta_step)

    r = g.r

    curve_is_convex_constraints = r[:-2]*(s*r[1:-1] - s2*r[2:]) + r[1:-1]*r[2:]*(c*s2 - c2*s)
    return np.append(curve_is_convex_constraints, curve_length_constraint)


def get_trajectory_length(g):
    return g.centre_trajectory_length


def get_starting_point(n, number_of_radii, stair_width, stair_height):
    starting_theta = 2*np.pi / n
    r0 = stair_width / starting_theta
    slope = stair_height/stair_width
    return StarWheel(theta=starting_theta, n=n, r=r0*np.exp(slope*np.linspace(0, starting_theta, number_of_radii)))


def compute_star_wheel(n, number_of_radii, stair_width, stair_height):
    x0 = get_starting_point(n, number_of_radii, stair_width, stair_height).radii_and_theta_array

    eq_constr = {'type': 'eq', 'fun': decorate_to_take_array(zero_constraints, n, stair_width, stair_height)}
    ineq_constr = {'type': 'ineq', 'fun': decorate_to_take_array(inequality_constraints, n, stair_width)}
    bounds_constr = optimize.Bounds(0.1, np.append(np.repeat(100, number_of_radii), 2*np.pi/n))

    optimisation_objective = decorate_to_take_array(get_trajectory_length, n);

    opt_res = optimize.minimize(optimisation_objective, x0, bounds=bounds_constr, constraints=(eq_constr, ineq_constr))
    return array_to_wheel_geometry(opt_res.x, n)
