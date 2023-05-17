import numpy as np

point_coordinate = np.loadtxt("../config/grid/grid_coordinates.csv", delimiter=",")  # 点坐标
point_neighbor = np.loadtxt("../config/grid/neighbor.csv", delimiter=",")  # 邻点
two_neighbor = np.loadtxt("../config/grid/two_neighbor.csv", delimiter=",")  # 交界处的特殊邻点
point_totals = len(point_coordinate)  # 点总数


count_Z = 0
count_P = 0
for i in range(0, point_totals):
    neighbor_z = int(point_neighbor[i][0])
    neighbor_r = int(point_neighbor[i][2])
    if neighbor_z != -1:
        dist = point_coordinate[neighbor_z][1] - point_coordinate[i][1]
        if dist != 0:
            print(i, ":", dist)
            count_Z += 1
    if neighbor_r != -1 and neighbor_r != -2:
        dist = point_coordinate[neighbor_r][0] - point_coordinate[i][0]
        if dist != 0:
            print(i, ":", dist)
            count_P += 1
print("完毕")
print(count_Z)
print(count_P)
