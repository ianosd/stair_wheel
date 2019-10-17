import matplotlib.pyplot as plt
from matplotlib import animation, rc
import numpy as np

def gen_stair_data(stair_width, stair_height, start_x=0, start_y=0, num_climbs=1):
    xs = []
    ys = []
    for i in range(num_climbs):
        xs += [start_x, start_x+stair_width]
        ys += [start_y, start_y]
        start_x += stair_width
        xs += [start_x, start_x]
        ys += [start_y, start_y+stair_height]
        start_y += stair_height
        
    return xs, ys

def plot_stairs(*args, **kwargs):
    plt.plot(*gen_stair_data(*args, **kwargs))

def create_wheel_plot_data(g, *args, **kwargs):
    if 'transform' in kwargs:
        transform = kwargs['transform']
    else:
        transform = np.eye(3)
        
    xx, yy = [], []
    thetas_one_segm = np.linspace(0, g.theta, g.r.size)
    for i in range(g.n):
        thetas = thetas_one_segm + i*2*np.pi/g.n
        x = g.r*np.cos(thetas)
        y = g.r*np.sin(thetas)
        
        joined = np.concatenate(([x], [y], np.ones((1, x.size))))
        transformed = transform.dot(joined)
    
        x = transformed[0, :]
        y = transformed[1, :]

        xx.append(x)
        yy.append(y)

    x = np.concatenate(xx)
    x = np.append(x, x[0])
    y = np.concatenate(yy)
    y = np.append(y, y[0])
    
    return {'contour': (x, y), 'centre': ([transform[0, 2]], [transform[1, 2]]) }


def plot_star_wheel(g, *args, **kwargs):
    if not args:
        args = ('b',)
    
    plot_data = create_wheel_plot_data(g, *args, **kwargs)
    
    p0 = plt.plot(*plot_data['contour'], *args)
    p1 = plt.plot(*plot_data['centre'], 'r*')
    
    return p0 + p1


def create_animation(g_result, stair_width, stair_height, n):
    fig, ax = plt.subplots()
    ax.set_xlim(-stair_width, stair_width*n)
    ax.set_ylim(-stair_height, 2*stair_height*n+30)
    ax.set_aspect('equal')

    lines = (ax.plot([], [], 'b:', label='result shape')[0], ax.plot([], [], 'r-', label="centre trajectory")[0], ax.plot(*gen_stair_data(stair_width, stair_height, -0.5, 0, n), 'g', label="stairs"))

    ax.legend()
    frame_count = n*20
    def animate(frame):
        data = create_wheel_plot_data(g_result, transform=g_result.get_transform(n*frame/frame_count))
        lines[0].set_data(*data['contour'])
        
        cum_data = (create_wheel_plot_data(g_result, transform=g_result.get_transform(n*f/frame_count)) for f in range(frame+1))
        centre_data = list(zip(*((d['centre'][0][0], d['centre'][1][0]) for d in cum_data)))
        lines[1].set_data(centre_data)
        return lines

    return animation.FuncAnimation(fig, animate, frames=frame_count, interval=50)

def export_manufacturing_material():
    import svgwrite
    import math

    from svglib.svglib import svg2rlg
    from reportlab.graphics import renderPDF

    def cm(num):
        if type(num) is tuple:
            return tuple(cm(x) for x in num)
        return str(num) + "cm"

    base_name="result_{}__{}__{}".format(g_result.n, str(stair_width).replace('.', '_'), str(stair_height).replace('.', '_'))
    svg_name = base_name + ".svg"
    pdf_name = base_name + ".pdf"
    d = svgwrite.Drawing(filename=base_name + ".svg", width="21cm", height="29.7cm")

    # draw support for splitting an angle in n pieces, where n is the n in our geometry definition
    support = d.add(d.g(stroke_width="1mm", stroke="black"))
    n = g_result.n
    r = 6
    c_x = r + 1
    c_y = r + 1

    line = d.line

    def add_x(d, point, size=0.5):
        d.add(line((cm(point[0] - size), cm(point[1] - size)), (cm(point[0] + size), cm(point[1] + size))))
        d.add(line((cm(point[0] - size), cm(point[1] + size)), (cm(point[0] + size), cm(point[1] - size))))
        
    add_x(support, (c_x, c_y))
    for i in range(n):
        start_point = (c_x + math.cos(i*2*math.pi/n)*r, c_y + math.sin(i*2*math.pi/n)*r)
        end_point = (c_x + math.cos(i*2*math.pi/n + g_result.Dtheta)*r, c_y + math.sin(i*2*math.pi/n + g_result.Dtheta)*r)
        add_x(support, start_point)
        add_x(support, end_point)
        support.add(line(cm(start_point), cm(end_point)))
        
        
    support.add(d.rect(cm((0, 0)), cm((21, 29.7)), fill_opacity=0))
    def draw_curve(drawing, geometry, tx, ty, phi):
        cos_phi = math.cos(phi)
        sin_phi = math.sin(phi)
        def transform(p):
            return (tx + cos_phi * p[0] - sin_phi * p[1], ty + sin_phi * p[0] + cos_phi * p[1])
        def line(x0, y0, x1, y1):
            p0 = transform((x0, y0))
            p1 = transform((x1, y1))
            return d.line(cm(p0), cm(p1))
        s = 0.90
        drawing.add(line(s*geometry.r[0], 0, (2-s)*geometry.r[0], 0))
        drawing.add(line(s*geometry.r[-1]*math.cos(geometry.Dtheta), s*geometry.r[-1]*math.sin(geometry.Dtheta), (2-s)*geometry.r[-1]*math.cos(geometry.Dtheta), (2-s)*geometry.r[-1]*math.sin(geometry.Dtheta)))
        
        thetas = np.linspace(0, geometry.Dtheta, len(geometry.r))
        for i in range(len(thetas) - 1):
            drawing.add(line(geometry.r[i]*math.cos(thetas[i]), geometry.r[i]*math.sin(thetas[i]), geometry.r[i+1]*math.cos(thetas[i+1]), geometry.r[i+1]*math.sin(thetas[i+1])))
        
    draw_curve(support, g_result, 1, 19 + g_result.r[0], 3*math.pi/2)
        
    d.save()
    renderPDF.drawToFile(svg2rlg(svg_name), pdf_name)