import numpy as np


# 读入浓度文件，只取前3360个点，存成数组
def generate_c_matrix(dir_name):
    temp = np.loadtxt(dir_name)
    c_matrix = temp[0:3360]
    return c_matrix


def cal_avg_concentration(dir_name, species):
    avg_c = np.zeros(701)
    ctrl_v = np.load("../config/ctrl_v.npy")
    total_v = np.pi * 298 ** 2 * 15.0
    for i in range(701):
        load_file_number = str(i * 100).zfill(8)
        c_matrix = generate_c_matrix(f"{dir_name}{species}/{species}{load_file_number}.csv")
        total_c = np.sum(c_matrix * ctrl_v) / 2
        avg_c[i] = total_c / total_v
    np.savetxt(f"{dir_name}avg_c_{species}.csv", avg_c)


if __name__ == '__main__':
    my_dir_name = "../Result/TS50000_kRyR=311999711.0000832_08-13-14-41-34/"
    cal_avg_concentration(my_dir_name, "Ca")
    cal_avg_concentration(my_dir_name, "CaF")
    print("success!!!")
