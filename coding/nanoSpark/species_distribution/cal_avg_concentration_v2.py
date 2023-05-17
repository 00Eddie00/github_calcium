import numpy as np
import os


def cal_avg_concentration(dir_name, species):
    filenames = os.listdir(f"{dir_name}{species}/")  # 所有步数的浓度文件
    ctrl_v = np.load("../config/ctrl_v.npy")  # 亚空间每个单元的控制体积
    total_v = np.pi * 298 ** 2 * 15.0  # 亚空间总体积
    avg_c = np.zeros(len(filenames))  # 存放每一步的亚空间平均浓度
    count = 0
    for filename in filenames:
        c_matrix = np.loadtxt(f"{dir_name}{species}/{filename}")[0:3360]
        total_c = np.sum(c_matrix * ctrl_v) / 2
        avg_c[count] = total_c / total_v
        count += 1
    np.savetxt(f"{dir_name}avg_c_{species}.csv", avg_c)


if __name__ == '__main__':
    my_dir_name = "../Result/TS100000_kRyR=311999711.0000832_01-16-20-55-32/"
    cal_avg_concentration(my_dir_name, "Ca")
    cal_avg_concentration(my_dir_name, "CaF")
    print("success!!!")
