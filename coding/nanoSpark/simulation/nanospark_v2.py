import numpy as np
import time

from utils.tool_mkdir import mkdir
from config.parameters import *
from calculating_iteration_v3 import cal_f_concentration, cal_g_h_concentration


# 初始化各物质浓度
def initialize_parameter():
    ini_ca = INITIAL_C_CA
    ini_caf = ini_ca * TOTAL_CAF / (ini_ca + K_DIV_CAF)
    ini_cam = ini_ca * TOTAL_CAM / (ini_ca + K_DIV_CAM)
    ini_trc = ini_ca * TOTAL_TRC / (ini_ca + K_DIV_TRC)
    ini_srm = ini_ca * TOTAL_SRM / (ini_ca + K_DIV_SRM)
    ini_slm = ini_ca * TOTAL_SLM / (ini_ca + K_DIV_SLM)

    f = np.full(POINT_TOTALS, ini_ca)
    g = np.zeros(POINT_TOTALS)
    g[:3360] = ini_caf
    g[3425:3450] = ini_caf
    h_1 = np.full(POINT_TOTALS, ini_cam)
    h_2 = np.full(POINT_TOTALS, ini_trc)
    h_3 = np.full(POINT_TOTALS, ini_srm)
    h_4 = np.full(POINT_TOTALS, ini_slm)

    return f, g, h_1, h_2, h_3, h_4


# 循环控制
def loop(total_steps, do_save, dir_name, is_continue):
    # 创建用于存储结果的文件夹
    mkdir(dir_name + "Ca")
    mkdir(dir_name + "CaF")
    mkdir(dir_name + "CaB")
    # ryr通道开放时间对应的步数
    release_step = int(RELEASE_TIME / DT)
    # 根据是否断点续算赋不同的值
    if is_continue:
        # 得到已经计算的步数：finished_step  和  接下来应该计算的步数：current_step
        with open(dir_name + "STATE.csv", 'r') as load_state_file:
            finished_step = int(load_state_file.readline())
        current_step = finished_step + 1
        # 控制ryr通道关闭或开启
        if current_step >= release_step:
            k_ryr = 0
        else:
            k_ryr = K_RYR
        # 载入钙离子、荧光钙离子、各缓冲物的浓度
        load_file_number = str(finished_step).zfill(8)
        f = np.loadtxt(f"{dir_name}Ca/Ca{load_file_number}.csv")
        g = np.loadtxt(f"{dir_name}CaF/CaF{load_file_number}.csv")
        cab = np.loadtxt(f"{dir_name}CaB/CaB{load_file_number}.csv", delimiter=',')
        h_1 = cab[:, 0]
        h_2 = cab[:, 1]
        h_3 = cab[:, 2]
        h_4 = cab[:, 3]
    else:
        current_step = 1
        k_ryr = K_RYR  # 开启ryr通道
        f, g, h_1, h_2, h_3, h_4 = initialize_parameter()
        # 记录初始时刻钙离子浓度值
        ca_file_name = f"{dir_name}Ca/Ca00000000.csv"
        np.savetxt(ca_file_name, f)
        print(ca_file_name, "SAVED")
        # 记录初始时刻荧光钙离子浓度值
        caf_file_name = f"{dir_name}CaF/CaF00000000.csv"
        np.savetxt(caf_file_name, g)
        print(caf_file_name, "SAVED")
    # 开始计算
    for i in range(current_step, total_steps + 1):
        print(f"------------ CURRENT_STEP {i} TOTAL_STEPS {total_steps} ------------")
        # 计算钙离子浓度、荧光钙离子浓度、各缓冲物浓度
        temp = cal_f_concentration(f, g, h_1, h_2, h_3, h_4, k_ryr)
        g, h_1, h_2, h_3, h_4 = cal_g_h_concentration(f, g, h_1, h_2, h_3, h_4)
        f = temp
        # 记录钙离子浓度、荧光钙离子浓度、各缓冲物浓度，每100步存一次
        if i % do_save == 0:
            file_number = str(i).zfill(8)
            # 记录钙离子浓度
            ca_file_name = f"{dir_name}Ca/Ca{file_number}.csv"
            np.savetxt(ca_file_name, f)
            print(ca_file_name, "SAVED")
            # 记录荧光钙离子浓度
            caf_file_name = f"{dir_name}CaF/CaF{file_number}.csv"
            np.savetxt(caf_file_name, g)
            print(caf_file_name, "SAVED")
            # 记录各缓冲物浓度
            cab_file_name = f"{dir_name}CaB/CaB{file_number}.csv"
            np.savetxt(cab_file_name,
                       np.concatenate((h_1[:, np.newaxis], h_2[:, np.newaxis], h_3[:, np.newaxis], h_4[:, np.newaxis]),
                                      axis=1), delimiter=',')
            print(cab_file_name, "SAVED")
            # 记录已经计算的步数
            with open(dir_name + "STATE.csv", "w") as state_file:
                state_file.write(str(i))
            print("CURRENT STATE SAVED")
        # 控制ryr通道关闭
        if i == release_step:
            k_ryr = 0


# 程序入口
def main():
    # total_steps = 100000
    # do_save = 100
    total_steps = 100
    do_save = 1
    is_continue = False

    # if is_continue:
    #     note = input("请输入文件目录名：")
    # else:
    #     parameters_info = f"kRyR={K_RYR}"
    #     current_time = time.strftime("%m-%d-%H-%M-%S", time.localtime())
    #     note = f"TS{total_steps}_{parameters_info}_{current_time}"
    # dir_name = f"../Result/{note}/"
    dir_name = ""
    loop(total_steps, do_save, dir_name, is_continue)


if __name__ == '__main__':
    main()
