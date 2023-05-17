import numpy as np
import os
from scipy.ndimage import convolve1d
from process_concentration_matrix import process_concentration


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    matrix = original_matrix.copy()
    # 第一维
    for i in np.arange(matrix.shape[1]):
        for j in np.arange(matrix.shape[2]):
            one_line = matrix[:, i, j]
            matrix[:, i, j] = convolve1d(one_line, kernel[0], mode='constant', cval=xy_c_val)
    # 第二维
    for i in np.arange(matrix.shape[0]):
        for j in np.arange(matrix.shape[2]):
            one_line = matrix[i, :, j]
            matrix[i, :, j] = convolve1d(one_line, kernel[1], mode='constant', cval=xy_c_val)
    # 第三维
    for i in np.arange(matrix.shape[0]):
        for j in np.arange(matrix.shape[1]):
            one_line = matrix[i, j, :]
            matrix[i, j, :] = convolve1d(one_line, kernel[2], mode='constant', cval=z_c_val)
    return matrix


# # 多物质时间分布(不包括CaG)
# def multi_temporal_distribution(result_set, kernel, xy_c_val, z_c_val):
#     ca_dir_name = f"{result_set}/Ca"
#     caf_dir_name = f"{result_set}/CaF"
#     ca_dir_list = os.listdir(ca_dir_name)
#     caf_dir_list = os.listdir(caf_dir_name)
#
#     ca = np.empty(len(ca_dir_list))
#     caf = np.empty(len(caf_dir_list))
#     ca_bc = np.empty(len(caf_dir_list))
#     caf_psf = np.empty(len(caf_dir_list))
#
#     count = 0
#     for i in ca_dir_list:
#         ca[count] = np.loadtxt(f"{ca_dir_name}/{i}")[23]
#         count += 1
#
#     count = 0
#     for i in caf_dir_list:
#         processed_con_matrix = process_concentration(np.loadtxt(f"{caf_dir_name}/{i}"), xy_c_val)
#         caf[count] = processed_con_matrix[300, 300, 15]
#         # back-calculation of [Ca2+] with [CaF]
#         ca_bc[count] = K_DIV_CAF * processed_con_matrix[300, 300, 15] / (TOTAL_CAF - processed_con_matrix[300, 300, 15])
#         # [CaF] convolved with PSF
#         # 0：边缘；100、200：中部；300：中心
#         caf_psf[count] = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)[300, 300, 15]
#         count += 1
#         print(count)
#     ca_bc_psf = K_DIV_CAF * caf_psf / (TOTAL_CAF - caf_psf)  # back-calculation of [Ca2+] with [CaF] convolved with PSF
#
#     return np.stack((ca, caf, caf_psf, ca_bc, ca_bc_psf), axis=-1)


# 荧光染料卷积后空间分布
def spatial_distribution(result_set, species, kernel, xy_c_val, z_c_val):
    fluo_file_name = f"{result_set}\\NANO\\{species}\\{species}00010000.csv"
    processed_con_matrix = process_concentration(np.loadtxt(fluo_file_name), xy_c_val)  # 处理得到荧光钙离子的浓度矩阵
    fluo_psf = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)[:, 300, 15]  # [CaF] convolved with PSF
    return fluo_psf


# 荧光染料卷积后时间分布
def temporal_distribution(result_set, species, kernel, xy_c_val, z_c_val):
    fluo_dir_name = f"{result_set}\\NANO\\{species}"
    fluo_dir_list = os.listdir(fluo_dir_name)
    fluo_dir_list_len = len(fluo_dir_list)
    position_con = np.empty(fluo_dir_list_len)

    for fluo_file_index in range(fluo_dir_list_len):
        print(fluo_file_index)
        # 处理浓度
        processed_con_matrix = process_concentration(np.loadtxt(f"{fluo_dir_name}\\{fluo_dir_list[fluo_file_index]}"),
                                                     xy_c_val)
        # 卷积
        con_result = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
        # 选取位置保存
        position_con[fluo_file_index] = con_result[300, 300, 15]  # 应该是这样吧
    return position_con


# 物质的时间分布
def temporal_distribution_no_conv(result_set, species, xy_c_val, position_list):
    fluo_dir_name = f"{result_set}\\NANO\\{species}"
    fluo_dir_list = os.listdir(fluo_dir_name)
    fluo_dir_list_len = len(fluo_dir_list)

    position_list_len = len(position_list)
    position_con = np.empty((position_list_len, fluo_dir_list_len))

    for fluo_file_index in range(fluo_dir_list_len):
        print(fluo_file_index)
        processed_con_matrix = process_concentration(np.loadtxt(f"{fluo_dir_name}\\{fluo_dir_list[fluo_file_index]}"),
                                                     xy_c_val)
        # 选取位置保存
        for position_index in range(position_list_len):
            i = position_list[position_index][0] + 300
            j = position_list[position_index][1] + 300
            position_con[position_index, fluo_file_index] = processed_con_matrix[i, j, 15]
    return position_con


# 物质的空间分布
def spatial_distribution_no_conv(result_set, species, xy_c_val):
    fluo_file_name = f"{result_set}\\NANO\\{species}\\{species}00010000.csv"
    processed_con = process_concentration(np.loadtxt(fluo_file_name), xy_c_val)[:, 300, 15]  # 处理得到荧光钙离子的浓度矩阵
    return processed_con


def optical_blurring(result_set, species):
    kernel = np.load("", allow_pickle=True)
    C_VAL = 0.0  # 用于纳米钙火花
    xy_c_val = C_VAL

    np.savetxt(f'{result_set}\\{species}_t_dis_conv.csv',
               temporal_distribution(result_set, species, kernel, xy_c_val, 0))


if __name__ == '__main__':
    my_result_set = ""
    my_species = "CaG"
