import numpy as np


def circle(triangle):
    x1, y1, x2, y2, x3, y3 = tuple(triangle)

    a = x1 - x2
    b = y1 - y2
    c = x1 - x3
    d = y1 - y3
    a1 = (x1 ** 2 - x2 ** 2 + y1 ** 2 - y2 ** 2) / 2.0
    a2 = (x1 ** 2 - x3 ** 2 + y1 ** 2 - y3 ** 2) / 2.0
    theta = b * c - a * d

    if (abs(theta)) < 1e-7:
        raise RuntimeError('3 points in line')

    x0 = (b * a2 - d * a1) / theta
    y0 = (c * a1 - a * a2) / theta
    r = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

    return x0, y0, r


def load_data():
    grid = np.loadtxt('4RYRbggridt.dat', dtype=np.float64)
    triangles_set = np.loadtxt('4RYRbgnod.dat', dtype=np.int64) - 1

    triangle = np.empty(6)
    circle_info = np.empty((len(triangles_set), 3))
    for i in range(len(triangles_set)):
        for j in range(0, 3):
            triangle[j * 2] = grid[triangles_set[i, j], 0]
            triangle[j * 2 + 1] = grid[triangles_set[i, j], 1]
        circle_info[i] = circle(triangle)

    np.savetxt("circle_info.csv", circle_info)


load_data()
