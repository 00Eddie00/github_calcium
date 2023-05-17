import numpy as np
import os
from scipy.ndimage import convolve1d
from process_concentrationv2 import process_concentration
from config.parameters_for_ob import C_VAL
from config.parameters import K_DIV_CAF, TOTAL_CAF

# 三维卷积
def convolve3d(matrix, kernel, xy_c_val, z_c_val):
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


# 时间分布
def temporal_distribution(dir_name):
    kernel = np.load("../config/optical_blurring/kernel.npy", allow_pickle=True)
    dir_list = os.listdir(dir_name)

    ca_bc = np.empty(len(dir_list))
    caf_psf = np.empty(len(dir_list))

    count = 0
    for i in dir_list:
        processed_con_matrix = process_concentration(np.loadtxt(f"{dir_name}//{i}"))
        # back-calculation of [Ca2+] with [CaF]
        # ca_bc[count] = K_DIV_CAF * processed_con_matrix[300, 300, 15] / (TOTAL_CAF - processed_con_matrix[300, 300, 15])
        # [CaF] convolved with PSF
        caf_psf[count] = convolve3d(processed_con_matrix, kernel, C_VAL, C_VAL)[300, :, :]
        count += 1
        print(count)
    # ca_bc_psf = K_DIV_CAF * caf_psf / (TOTAL_CAF - caf_psf)  # back-calculation of [Ca2+] with [CaF] convolved with PSF

    # return np.stack((ca_bc, caf_psf, ca_bc_psf), axis=-1)

    return caf_psf

def new_temporal_distribution(dir_name):
    kernel = np.load("../config/optical_blurring/kernel.npy", allow_pickle=True)
    dir_list = os.listdir(dir_name)

    count = 0
    for i in dir_list:
        processed_con_matrix = process_concentration(np.loadtxt(f"{dir_name}//{i}"))
        result = convolve3d(processed_con_matrix, kernel, C_VAL, C_VAL)[300, :, :]
        count += 1
        print(count)
        np.savetxt(f"..\\Result\\TS70000_kRyR=311999711.0000832_09-23-11-38-22-new_grid\\convolved_result\\{i}", result)

# 空间分布
def spatial_distribution(file_name):
    kernel = np.load("../config/optical_blurring/kernel.npy", allow_pickle=True)  # 导入卷积核
    processed_con_matrix = process_concentration(np.loadtxt(file_name))  # 处理得到荧光钙离子的浓度矩阵
    # back-calculation of [Ca2+] with [CaF]
    # ca_bc = K_DIV_CAF * processed_con_matrix[300, :, 15] / (TOTAL_CAF - processed_con_matrix[300, :, 15])
    caf_psf = convolve3d(processed_con_matrix, kernel, C_VAL, C_VAL)[300, :, :]  # [CaF] convolved with PSF
    # ca_bc_psf = K_DIV_CAF * caf_psf / (TOTAL_CAF - caf_psf)  # back-calculation of [Ca2+] with [CaF] convolved with PSF

    # return np.stack((ca_bc, caf_psf, ca_bc_psf), axis=-1)
    return caf_psf

if __name__ == '__main__':
    path = "..\Result\TS70000_kRyR=311999711.0000832_09-23-11-38-22-new_grid"
    # np.savetxt(f"{path}/t_dis.csv", temporal_distribution(f"{path}/CaF"))
    # np.save("./x_dis_for_contour", spatial_distribution("../Result/TS70000_kRyR=311999711.0000832_09-13-14-56-25-old_grid/CaF/CaF00010000.csv"))
    new_temporal_distribution(f"{path}/CaF")