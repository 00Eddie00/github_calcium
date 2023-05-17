import numpy as np
import os


def cal_avg_concentration(result_set, area, species):
    filenames = os.listdir(f"{result_set}/{area}/{species}/")  # 所有步数的浓度文件
    ctrl_v = np.load("./nano_grid/blink_mirrior/ctrl_triangle_area.npy")  # 亚空间每个单元的控制体积
    # total_v = 282662.89171587746 # blink_style
    total_v = 158832.085418743  # 450 nm
    avg_c = np.zeros(len(filenames))  # 存放每一步的亚空间平均浓度
    count = 0
    for filename in filenames:
        c_matrix = np.loadtxt(f"{result_set}/{area}/{species}/{filename}", delimiter=",")
        total_c = np.sum(c_matrix * ctrl_v) / 3
        avg_c[count] = total_c / total_v
        count += 1
    np.savetxt(f"{result_set}/{area}_{species}_avg_c.csv", avg_c)


if __name__ == '__main__':
    my_result_set = "../FixedUnstructuredResultSets/ResultSet_22-0901-2228_情况六-子情况2"
    cal_avg_concentration(my_result_set, "NANO", "Ca")
    # cal_avg_concentration(my_result_set, "CaF")
    # cal_avg_concentration(my_result_set, "CaG")
    print("success!")
