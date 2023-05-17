import numpy as np

from config.parameters import *

grid_coordinates = np.loadtxt("../config/new_grid/grid_coordinates.csv", delimiter=",")  # 点坐标
neighbor = np.loadtxt("../config/new_grid/neighbor.csv", int, delimiter=",")  # 邻点
two_neighbor = np.loadtxt("../config/new_grid/two_neighbor.csv", int, delimiter=",")  # 交界处的特殊邻点
coefficient = np.load("../config/new_coefficient/coefficient.npy")


def cal_j_dye(c_ca: np.ndarray, c_caf: np.ndarray) -> np.ndarray:
    """ 计算J_dye

    :param c_ca: 肌质中的钙离子浓度 f
    :param c_caf: 肌质中的荧光钙离子浓度 g
    :return: J_dye
    """
    j_dye = np.zeros(len(c_caf))
    j_dye[:3360] = -K_PLUS_CAF * c_ca[:3360] * (TOTAL_CAF - c_caf[:3360]) + K_MINUS_CAF * c_caf[:3360]
    j_dye[3425:3450] = -K_PLUS_CAF * c_ca[3425:3450] * (TOTAL_CAF - c_caf[3425:3450]) + K_MINUS_CAF * c_caf[3425:3450]
    return j_dye


def cal_j_buffer(c_ca: np.ndarray, c_cam: np.ndarray, c_trc: np.ndarray, c_srm: np.ndarray, c_slm: np.ndarray):
    """计算J_i

    :param c_ca: 肌质中的钙离子浓度
    :param c_cam:
    :param c_trc:
    :param c_srm:
    :param c_slm:
    :return: j_cam, j_trc, j_srm, j_slm
    """
    j_cam = -K_PLUS_CAM * c_ca * (TOTAL_CAM - c_cam) + K_MINUS_CAM * c_cam
    j_trc = -K_PLUS_TRC * c_ca * (TOTAL_TRC - c_trc) + K_MINUS_TRC * c_trc
    j_srm = -K_PLUS_SRM * c_ca * (TOTAL_SRM - c_srm) + K_MINUS_SRM * c_srm
    j_slm = -K_PLUS_SLM * c_ca * (TOTAL_SLM - c_slm) + K_MINUS_SLM * c_slm
    return j_cam, j_trc, j_srm, j_slm


def cal_j_ryr(k_ryr):
    j_ryr = np.zeros(POINT_TOTALS, float)

    if k_ryr == 0:
        return j_ryr

    for i in range(POINT_TOTALS):
        if 23 <= i <= 263 and (i - 23) % 24 == 0:
            j_ryr[i] = k_ryr * C_CA_JSR

    return j_ryr


def det_f_r_neighbor(i, r_neighbor, iteration):
    z2 = grid_coordinates[i][0]
    f_r_neighbor = iteration[i]

    if r_neighbor != -1 and r_neighbor != -2:
        f_r_neighbor = iteration[r_neighbor]
    elif r_neighbor == -2:  # 没对齐
        for j in range(len(two_neighbor)):
            if two_neighbor[j][0] == i:
                r_neighbor_1 = two_neighbor[j][1]
                r_neighbor_2 = two_neighbor[j][2]
                dist_1 = abs(grid_coordinates[r_neighbor_2][0] - z2)
                dist_2 = abs(grid_coordinates[r_neighbor_1][0] - z2)
                dist = dist_1 + dist_2
                f_r_neighbor = (iteration[r_neighbor_1] * dist_1 + iteration[r_neighbor_2] * dist_2) / dist
                break

    return f_r_neighbor


def cal_f_concentration(f, g, h_1, h_2, h_3, h_4, k_ryr):
    denominator = np.zeros(POINT_TOTALS)
    iteration = np.copy(f)
    temp = np.empty(POINT_TOTALS)  # 存储当次迭代的值

    j_dye = cal_j_dye(iteration, g)
    j_cam, j_trc, j_srm, j_slm = cal_j_buffer(c_ca=iteration, c_cam=h_1, c_trc=h_2, c_srm=h_3, c_slm=h_4)
    j_buffers = j_cam + j_trc + j_srm + j_slm
    j_ryr = cal_j_ryr(k_ryr)
    j_item = j_dye + j_buffers + j_ryr

    for i in range(POINT_TOTALS):
        if 23 <= i <= 263 and (i - 23) % 24 == 0:
            denominator[i] = 1 - D_CA * np.sum(coefficient[i] * np.array([-1, 1, -1, 1])) * DT + k_ryr * DT
        else:
            denominator[i] = 1 - D_CA * np.sum(coefficient[i] * np.array([-1, 1, -1, 1])) * DT
    # 进入迭代
    for k in range(ITERATION_TIMES):
        for i in range(POINT_TOTALS):
            p32, p12, p23, p21 = neighbor[i]
            # 保证此次迭代邻居的值用的是上一次的
            if p32 != -1:  # z轴正方向有
                f32 = iteration[p32]
            else:  # z轴正方向没有
                f32 = iteration[i]
            if p12 != -1:
                f12 = iteration[p12]
            else:
                f12 = iteration[i]
            f23 = det_f_r_neighbor(i, p23, iteration)
            f21 = det_f_r_neighbor(i, p21, iteration)

            temp[i] = (D_CA * np.sum(coefficient[i] * np.array([f32, -f12, f23, -f21])) * DT + j_item[i] * DT + f[
                i]) / denominator[i]
        iteration = np.copy(temp)

    return iteration


# 计算各点的CaF浓度、CaBi浓度
def cal_g_h_concentration(f, g, h_1, h_2, h_3, h_4):
    j_dye = cal_j_dye(c_ca=f, c_caf=g)
    j_cam, j_trc, j_srm, j_slm = cal_j_buffer(c_ca=f, c_cam=h_1, c_trc=h_2, c_srm=h_3, c_slm=h_4)

    g_1_new = DT * (-j_dye) + g

    h_1_new = DT * (-j_cam) + h_1
    h_2_new = DT * (-j_trc) + h_2
    h_3_new = DT * (-j_srm) + h_3
    h_4_new = DT * (-j_slm) + h_4

    return g_1_new, h_1_new, h_2_new, h_3_new, h_4_new
