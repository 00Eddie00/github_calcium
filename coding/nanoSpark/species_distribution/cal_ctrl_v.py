import numpy as np

grid_coordinates = np.loadtxt("../config/grid/grid_coordinates.csv", delimiter=",")  # 点坐标
neighbor = np.loadtxt("../config/grid/neighbor.csv", int, delimiter=",")  # 邻点
two_neighbor = np.loadtxt("../config/grid/two_neighbor.csv", int, delimiter=",")  # 交界处的特殊邻点


def radius(i, r_neighbor):
    if r_neighbor == -2:
        for j in range(len(two_neighbor)):
            if two_neighbor[j][0] == i:
                r_neighbor_1 = two_neighbor[j][1]
                break
        r = grid_coordinates[r_neighbor_1][1]
    else:
        r = grid_coordinates[r_neighbor][1]
    return r


def cal_ctrl_v():
    ctrl_v = np.zeros(3360)

    for i in range(3360):
        z2 = grid_coordinates[i][0]
        r2 = grid_coordinates[i][1]
        p32, p12, p23, p21 = neighbor[i]  # 上下外内

        if p32 == -1:
            z3 = z2
            z1 = grid_coordinates[p12][0]
        else:
            z3 = grid_coordinates[p32][0]
            if p12 == -1:
                z1 = z2
            else:
                z1 = grid_coordinates[p12][0]

        if p23 == -1:
            r3 = r2
            r1 = grid_coordinates[p21][1]
        else:
            r3 = radius(i, p23)
            if p21 == -1:
                r1 = r2
            else:
                r1 = radius(i, p21)

        r1_5 = (r2 + r1) / 2
        r2_5 = (r2 + r3) / 2
        ctrl_v[i] = np.pi * (r2_5 ** 2 - r1_5 ** 2) * (z3 - z1)

    return ctrl_v


if __name__ == '__main__':
    np.save('../config/ctrl_v', cal_ctrl_v())
    print("success")
