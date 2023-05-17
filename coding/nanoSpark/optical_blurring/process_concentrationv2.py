import numpy as np
from math import sqrt
from config.parameters_for_ob import nano_point_totals, xy_length, z_length, C_VAL

# 对原始数据预处理，生成可用的极坐标系下的二维浓度矩阵
# 生成一个存储亚空间r轴坐标的一维数组
def pretreatment(original_data):
    # 提取原始数据中亚空间中的点，并转换为141*24的二维矩阵
    temp_1 = original_data[0:nano_point_totals].reshape((140, 24))
    temp_2 = np.concatenate((original_data[3425:3428], original_data[3429:3450]))  # 交界处（300nm处）的点
    temp = np.concatenate((temp_1, temp_2[np.newaxis, :]))  # 第一维表示r轴，第二维表示z轴
    # z轴按1nm的间隔保留16个点
    polar_concentration = np.empty((len(temp), z_length + 1))
    for i in range(len(temp)):
        for j in range(8):
            polar_concentration[i][j] = temp[i][j]
        for j in range(8, 16):
            polar_concentration[i][j] = temp[i][2 * j - 7]

    # 存储亚空间r轴坐标
    nano_r1 = np.linspace(0, 15, 31)  # 间距为 0.5nm 的部分
    nano_r2 = np.linspace(15, 50, 36)  # 间距为 1nm 的部分
    nano_r3 = np.linspace(50, 100, 26)  # 间距为 2nm 的部分
    nano_r4 = np.linspace(100, 300, 51)  # 间距为 4nm 的部分
    nano_r = np.concatenate((nano_r1[:-1], nano_r2[:-1], nano_r3[:-1], nano_r4))

    return polar_concentration, nano_r


# 计算目标点在ρ轴上的前后邻居的ρ坐标，即半径radius
def cal_neighbor_radius(radius):
    if 0 <= radius < 15:
        border = 0
        spacing = 0.5
    elif 15 <= radius < 50:
        border = 15
        spacing = 1
    elif 50 <= radius < 100:
        border = 50
        spacing = 2
    else:
        border = 100
        spacing = 4

    pre_radius = ((radius - border) // spacing) * spacing + border
    if (radius - border) % spacing != 0:
        next_radius = pre_radius + spacing
    else:
        next_radius = -1
    return pre_radius, next_radius


# 得到目标点的浓度
def get_concentration(radius, polar_concentration, nano_r):
    pre_radius, next_radius = cal_neighbor_radius(radius)
    pre_concentration = polar_concentration[nano_r.searchsorted(pre_radius)]  # 前邻居的浓度
    # 若两个坐标系下的点正好对应，则前邻居的浓度即所求浓度，否则进行插值计算
    if next_radius != -1:
        next_concentration = polar_concentration[nano_r.searchsorted(next_radius)]  # 后邻居的浓度
        # 插值计算
        dist = next_radius - pre_radius
        dist_pre = radius - pre_radius
        dist_next = next_radius - radius
        concentration = (dist_pre * next_concentration + dist_next * pre_concentration) / dist
    else:
        concentration = pre_concentration
    return concentration


# 拼接?
def concatenate_matrix(matrix):
    b_matrix = np.flip(matrix[1:], 0)  # flip：axis=0：上下翻转，意味着把行看成整体，行的顺序发生颠倒，每一行的元素不发生改变
    c_matrix = np.flip(b_matrix[:, 1:], 1)  # axis=1：左右翻转，意味着把列看成整体，列的顺序发生颠倒，每一列的元素不发生改变
    d_matrix = np.flip(matrix[:, 1:], 1)
    result = np.block([[[c_matrix], [b_matrix]], [[d_matrix], [matrix]]])
    return result


# 生成可用的亚空间浓度矩阵
def process_concentration(original_data):
    # 1. 对原始浓度数据进行预处理得到一个二维矩阵（polar_concentration，141*16）：存储极坐标系下亚空间中各点浓度
    #    生成一个一维数组（nano_r，140）：存储亚空间r轴坐标值
    polar_concentration, nano_r = pretreatment(original_data)

    # 2. 坐标转换
    rectangle_concentration = np.full((xy_length + 1, xy_length + 1, z_length + 1), C_VAL, dtype=float)
    for i in range(xy_length + 1):
        for j in range(xy_length + 1):
            # 2.1 计算直角坐标系下各点对应极坐标系下的r轴坐标值，
            radius = sqrt(i ** 2 + j ** 2)
            # 2.2 根据radius确定浓度
            if radius <= 300:  # 如果位于亚空间内部，就需要去 polar_concentration 中得到相应浓度
                result = get_concentration(radius, polar_concentration, nano_r)
                rectangle_concentration[i][j] = result

    # 3. 将rectangle_concentration拼接成完整的浓度矩阵
    processed_con_matrix = concatenate_matrix(rectangle_concentration)
    return processed_con_matrix

if __name__ == "__main__":
    pretreatment(np.loadtxt("..\Result\TS70000_kRyR=311999711.0000832_09-13-14-56-25\CaF\CaF00000100.csv"))