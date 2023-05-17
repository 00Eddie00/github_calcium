import numpy as np

grid_coordinates = np.loadtxt("../config/new_grid/grid_coordinates.csv", delimiter=",")  # 点坐标
neighbor = np.loadtxt("../config/new_grid/neighbor.csv", int, delimiter=",")  # 邻点
two_neighbor = np.loadtxt("../config/new_grid/two_neighbor.csv", int, delimiter=",")  # 交界处的特殊邻点


#  确定r轴邻居的r坐标
def radius(i, r_neighbor):
    r_neighbor_1 = -1
    if r_neighbor == -2:
        for j in range(len(two_neighbor)):
            if two_neighbor[j][0] == i:
                r_neighbor_1 = two_neighbor[j][1]
                break
        r = grid_coordinates[r_neighbor_1][1]
    else:
        r = grid_coordinates[r_neighbor][1]
    return r


def calculating_coefficient():
    coefficient = np.empty((len(grid_coordinates), 4))

    for i in range(len(grid_coordinates)):
        z2 = grid_coordinates[i][0]
        r2 = grid_coordinates[i][1]
        p32, p12, p23, p21 = neighbor[i]  # 上、下、外、内
        # 计算顶面系数
        if p32 == -1:  # 不存在上顶面，则一定存在下顶面
            coe1 = 0
            z3 = z2
            z1 = grid_coordinates[p12][0]
            h = z3 - z1
            coe2 = 1 / (2 * h * (z1 - z2))
        else:  # 存在上顶面的情况
            z3 = grid_coordinates[p32][0]
            if p12 == -1:  # 不存在下顶面，则一定存在上顶面
                coe2 = 0
                z1 = z2
                h = z3 - z1
                coe1 = 1 / (2 * h * (z3 - z2))
            else:
                z1 = grid_coordinates[p12][0]
                h = z3 - z1
                coe1 = 1 / (2 * h * (z3 - z2))
                coe2 = 1 / (2 * h * (z1 - z2))
        # 计算侧面系数
        if p23 == -1:  # 不存在外侧面，则一定存在内侧面
            coe3 = 0
            r3 = r2
            r1 = grid_coordinates[p21][1]
            coe4 = r1 / ((r3 ** 2 - r1 ** 2) * (r1 - r2))
        else:  # 存在外侧面的情况
            r3 = radius(i, p23)
            if p21 == -1:
                coe4 = 0
                r1 = r2
                coe3 = r3 / ((r3 ** 2 - r1 ** 2) * (r3 - r2))
            else:
                r1 = radius(i, p21)
                coe3 = r3 / ((r3 ** 2 - r1 ** 2) * (r3 - r2))
                coe4 = r1 / ((r3 ** 2 - r1 ** 2) * (r1 - r2))

        coefficient[i] = coe1, coe2, coe3, coe4

    # np.savetxt("./config/coefficient/coefficient_v2.csv", coefficient, delimiter=',')
    np.save("../config/new_coefficient/coefficient", coefficient)


if __name__ == '__main__':
    calculating_coefficient()
