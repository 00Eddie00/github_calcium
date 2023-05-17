import numpy as np

grid_coordinates = np.loadtxt("../config/new_grid/grid_coordinates.csv", delimiter=",")
point_nums = len(grid_coordinates)
# 建立一个二维数组存储各点的邻居；第一维表示点编号，第二维表示四个方向的此点的邻居
# 索引0：z轴正方向，索引1：z轴负方向，索引2：ρ轴正方向，索引3：ρ负方向；
# -1 代表此方向上无邻居；-2代表此方向有两个邻居，需要求中点。
neighbor = np.full((point_nums, 4), -1)
two_neighbor = np.empty((23, 3))

# 纳米区域  z:-7.5_7.5 r:0_296
for i in range(3360):
    # z轴邻居
    if -7.5 < grid_coordinates[i][0] < 7.5:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif -7.5 == grid_coordinates[i][0]:
        neighbor[i][0] = i + 1
    else:
        neighbor[i][1] = i - 1
    # r轴邻居
    if 0 < grid_coordinates[i][1] < 296:
        neighbor[i][2] = i + 24
        neighbor[i][3] = i - 24
    elif 0 == grid_coordinates[i][1]:
        neighbor[i][2] = i + 24
    else:  # r=296处即衔接处附近
        for j in range(3360, 3515):
            if grid_coordinates[j][0] == grid_coordinates[i][0] and grid_coordinates[j][1] == 300:
                neighbor[i][2] = j
        neighbor[i][3] = i - 24

# 交界处 z：-2475_2475 r:300
count = 0
for i in range(3360, 3515):
    # z轴邻居
    if -2475 < grid_coordinates[i][0] < 2475:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif -2475 == grid_coordinates[i][0]:
        neighbor[i][0] = i + 1
    else:
        neighbor[i][1] = i - 1
    # r轴邻居
    if grid_coordinates[i][0] < -7.5:
        neighbor[i][2] = i + 155
        neighbor[i][3] = i + 13951
    elif grid_coordinates[i][0] > 7.5:
        neighbor[i][2] = i + 133
        neighbor[i][3] = i + 15486
    else:  # -7.5 <= x <= 7.5
        if grid_coordinates[i][0] <= -5.5:
            neighbor[i][2] = -2
            two_neighbor[count][0] = i
            two_neighbor[count][1] = 3579
            two_neighbor[count][2] = 3580
            count += 1
            neighbor[i][3] = i - 89
        elif grid_coordinates[i][0] == -5:
            neighbor[i][2] = 3580
            neighbor[i][3] = -2
            two_neighbor[count][0] = i
            two_neighbor[count][1] = 3338
            two_neighbor[count][2] = 3339
            count += 1
        elif grid_coordinates[i][0] <= -0.5:
            neighbor[i][2] = -2
            two_neighbor[count][0] = i
            two_neighbor[count][1] = 3580
            two_neighbor[count][2] = 3581
            count += 1
        elif grid_coordinates[i][0] == 0:
            neighbor[i][2] = 3581
        elif grid_coordinates[i][0] <= 4.5:
            neighbor[i][2] = -2
            two_neighbor[count][0] = i
            two_neighbor[count][1] = 3581
            two_neighbor[count][2] = 3582
            count += 1
        elif grid_coordinates[i][0] == 5:
            neighbor[i][2] = 3582
        else:
            neighbor[i][2] = -2
            two_neighbor[count][0] = i
            two_neighbor[count][1] = 3582
            two_neighbor[count][2] = 3583
            count += 1

    if -4.5 <= grid_coordinates[i][0] <= 7.5:
        neighbor[i][3] = i - 90

# 开放区域 z:-2475_2475 r:305_4946
for i in range(3515, 15751):
    # z轴邻居
    if -2475 < grid_coordinates[i][0] < 2475:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif -2475 == grid_coordinates[i][0]:
        neighbor[i][0] = i + 1
    else:
        neighbor[i][1] = i - 1
    # r轴邻居
    if 305 < grid_coordinates[i][1] < 4946:
        neighbor[i][2] = i + 133
        neighbor[i][3] = i - 133
    elif 305 == grid_coordinates[i][1]:
        neighbor[i][2] = i + 133
        for j in range(3360, 3515):
            if grid_coordinates[j][0] == grid_coordinates[i][0] and grid_coordinates[i][1] == 300:
                neighbor[i][3] = j
    else:
        neighbor[i][3] = i - 133

# 新开放区域2 z：-2475 ~ -11 r: 0 - 295
for i in range(15751, 17376):
    # z轴邻居
    if -2475 < grid_coordinates[i][0] < -11:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif -2475 == grid_coordinates[i][0]:
        neighbor[i][0] = i + 1
    elif -11 == grid_coordinates[i][0]:
        neighbor[i][1] = i - 1
    # r轴邻居
    if 0 < grid_coordinates[i][1] < 295:
        neighbor[i][2] = i + 65
        neighbor[i][3] = i - 65
    elif 0 == grid_coordinates[i][1]:
        neighbor[i][2] = i + 65
    else:  # 295 == grid_coordinates[i][1]:
        neighbor[i][2] = i - 13951
        neighbor[i][3] = i - 65

# 新开放区域1 z：11 ~ 2475 r: 0 - 295
for i in range(17376, point_nums):
    # z轴邻居
    if 11 < grid_coordinates[i][0] < 2475:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif 11 == grid_coordinates[i][0]:
        neighbor[i][0] = i + 1
    elif 2475 == grid_coordinates[i][0]:
        neighbor[i][1] = i - 1
    # r轴邻居
    if 0 < grid_coordinates[i][1] < 295:
        neighbor[i][2] = i + 65
        neighbor[i][3] = i - 65
    elif 0 == grid_coordinates[i][1]:
        neighbor[i][2] = i + 65
    else:  # 295 == grid_coordinates[i][1]:
        neighbor[i][2] = i - 15486
        neighbor[i][3] = i - 65

np.savetxt("../config/new_grid/neighbor.csv", neighbor, fmt='%d', delimiter=",")
np.savetxt("../config/new_grid/two_neighbor.csv", two_neighbor, fmt='%d', delimiter=",")
