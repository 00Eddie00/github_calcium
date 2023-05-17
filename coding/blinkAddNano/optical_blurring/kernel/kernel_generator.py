import numpy as np
from math import exp, sqrt
from optical_blurring.kernel.kernel_parameters import xy_count, z_count, sigma_xy, sigma_z


# 高斯函数
def gaussian_function(w, sigma):
    result = exp(-(w ** 2) / (2 * sigma)) / sqrt(2 * np.pi * sigma)
    return result


# 计算xy_count z_count
def cal_count(sigma):
    base_line = gaussian_function(0, sigma)
    current_value = base_line
    count = 1
    while current_value > base_line * 0.01:
        current_value = gaussian_function(count, sigma)
        count += 1
    print(count)


# 生成卷积核
def generate_kernel():
    # x、y方向上
    xy_parameter = np.array([gaussian_function(i, sigma_xy) for i in range(xy_count + 1)])
    # z方向上
    z_parameter = np.array([gaussian_function(i, sigma_z) for i in range(z_count + 1)])
    # 拼接
    xy_parameter = np.concatenate((xy_parameter[xy_count + 1:0:-1], xy_parameter))
    z_parameter = np.concatenate((z_parameter[z_count + 1:0:-1], z_parameter))

    kernel = np.asarray([xy_parameter, xy_parameter, z_parameter], dtype='object')
    return kernel


if __name__ == '__main__':
    np.save("./kernel_v1", generate_kernel())
