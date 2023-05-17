import numpy as np
import time

from utils.tool_mkdir import mkdir
from config.parameters import NR, CCACYTREST, DT, KRYR2, BCAF, KDCAF, BCAM, KDCAM, BTRC, KDTRC, BSRM, KDSRM, BSLM, \
    KDSLM, RELEASE_TIME
from simulation.calculating import generate_coefficient_matrix, generate_constant_matrix, cal_f_g_concentration, \
    cal_h_concentration, \
    cal_avg_concentration


# **************************************************导入网格信息**************************************************
def load_grid_info():
    spark_grid = np.loadtxt("../config/spark_mesh.csv", delimiter=",")
    return spark_grid


# **************************************************初始化各物质浓度***************************************************
def initial_parameter():
    # 各物质的初始值
    ini_ca = CCACYTREST
    ini_caf = ini_ca * BCAF / (ini_ca + KDCAF)
    ini_cam = ini_ca * BCAM / (ini_ca + KDCAM)
    ini_trc = ini_ca * BTRC / (ini_ca + KDTRC)
    ini_srm = ini_ca * BSRM / (ini_ca + KDSRM)
    ini_slm = ini_ca * BSLM / (ini_ca + KDSLM)

    f = np.full(NR, ini_ca)
    g = np.full(NR, ini_caf)
    h_1 = np.full(NR, ini_cam)
    h_2 = np.full(NR, ini_trc)
    h_3 = np.full(NR, ini_srm)
    h_4 = np.full(NR, ini_slm)

    return f, g, h_1, h_2, h_3, h_4


# ***************************************************循环控制***************************************************
def not_empty(s):
    return s and s.strip()


# 循环
def loop(total_steps, do_save, dir_name, is_flag):
    f, g, h_1, h_2, h_3, h_4 = initial_parameter()

    # 创建文件夹
    mkdir(dir_name + "Ca")
    mkdir(dir_name + "CaF")
    mkdir(dir_name + "CaB")

    # 导入网格
    spark_grid = load_grid_info()
    k_ryr = KRYR2
    release_time = RELEASE_TIME  # 0.02s,这个就是ryr通道开放的时间，后80毫秒是恢复的
    release_step = int(release_time / DT)  # 释放时间的迭代次数
    current_step = 1  # 目前执行的步数

    if is_flag:
        with open(dir_name + "STATE.dat", 'r') as load_state_file:
            finished_step = int(load_state_file.readline())
        current_step = finished_step + 1
        load_file_number = str(finished_step).zfill(8)

        # 载入Ca.dat
        load_ca_file_name = dir_name + "Ca\\Ca" + load_file_number + '.dat'
        with open(load_ca_file_name, 'r') as load_ca_file:
            count = 0
            for line in load_ca_file.readlines():
                if count < NR:
                    f[count] = float(line)
                count = count + 1
        # 载入CaF.dat
        load_caf_file_name = dir_name + "CaF\\CaF" + load_file_number + '.dat'
        with open(load_caf_file_name, 'r') as load_caf_file:
            count = 0
            for line in load_caf_file.readlines():
                if count < NR:
                    g[count] = float(line)
                count = count + 1
        # 载入CaB.dat
        load_cab_file_name = dir_name + "CaB\\CaB" + load_file_number + '.dat'
        with open(load_cab_file_name, 'r') as load_cab_file:
            count = 0
            for line in load_cab_file.readlines():
                current_line = list(filter(not_empty, line.strip("\n").split(" ")))
                if count < NR:
                    h_1[count] = float(current_line[0])
                    h_2[count] = float(current_line[1])
                    h_3[count] = float(current_line[2])
                    h_4[count] = float(current_line[3])
                count = count + 1

        if current_step >= release_step:
            k_ryr = 0
    else:
        avg_ca = cal_avg_concentration(spark_grid, f)
        avg_caf = cal_avg_concentration(spark_grid, g)
        # Ca.dat文件
        ca_file_name = dir_name + "Ca\\Ca00000000.dat"
        with open(ca_file_name, "w") as ca_file:
            for j in range(0, NR):
                ca_file.write(str(f[j]) + "\n")
            ca_file.write(str(avg_ca) + "\n")
        print(ca_file_name, "SAVED")

        # CaF.dat文件
        caf_file_name = dir_name + "CaF\\CaF00000000.dat"
        with open(caf_file_name, "w") as caf_file:
            for j in range(0, NR):
                caf_file.write(str(g[j]) + "\n")
            caf_file.write(str(avg_caf) + "\n")
            caf_file.write(str(j) + "\n")
        print(caf_file_name, "SAVED")


    # 生成系数矩阵
    coe_spark, coe_caf = generate_coefficient_matrix(spark_grid, k_ryr)

    for i in range(current_step, total_steps + 1):
        print("------------ CURRENT_STEP", i, "TOTAL_STEPS", total_steps, "------------")
        # 生成常数矩阵
        # 要用到各个缓冲物n时刻的值，所以先计算
        const_spark, const_caf = generate_constant_matrix(f, g, h_1, h_2, h_3, h_4, spark_grid, k_ryr)

        # 计算
        # 计算缓冲物浓度时需要使用n时刻的Ca离子浓度，因此先计算缓冲物浓度
        h_1, h_2, h_3, h_4 = cal_h_concentration(f, h_1, h_2, h_3, h_4)
        print("*******************")
        print("h_1：")
        print(h_1)
        print("*******************")
        # 计算钙离子浓度
        f = cal_f_g_concentration(coe_spark, const_spark)
        print("f：")
        print(f)
        print("*******************")
        # 计算荧光钙离子浓度
        g = cal_f_g_concentration(coe_caf, const_caf)

        # 输出数据
        if i % do_save == 0:
            file_number = str(i).zfill(8)

            # 计算平均值
            avg_ca = cal_avg_concentration(spark_grid, f)
            avg_caf = cal_avg_concentration(spark_grid, g)

            # Ca.dat文件
            ca_file_name = dir_name + "Ca\\Ca" + file_number + '.dat'
            with open(ca_file_name, "w") as ca_file:
                for j in range(0, NR):
                    ca_file.write(str(f[j]) + "\n")
                ca_file.write(str(avg_ca) + "\n")
                ca_file.write(str(i) + "\n")
            print(ca_file_name, "SAVED")

            # CaF.dat文件
            caf_file_name = dir_name + "CaF\\CaF" + file_number + '.dat'
            with open(caf_file_name, "w") as caf_file:
                for j in range(0, NR):
                    caf_file.write(str(g[j]) + "\n")
                caf_file.write(str(avg_caf) + "\n")
                caf_file.write(str(i) + "\n")
            print(caf_file_name, "SAVED")

            # CaB.dat文件
            cab_file_name = dir_name + "CaB\\CaB" + file_number + '.dat'
            with open(cab_file_name, "w") as cab_file:
                for j in range(0, NR):
                    cab_file.write(str(h_1[j]) + " " + str(h_2[j]) + " " + str(h_3[j]) + " " + str(h_4[j]) + "\n")
                # cab_file.write(str(avg_cab1) + " " + str(avg_cab2) + " " + str(avg_cab3) + " " + str(avg_cab4) + "\n")
                cab_file.write(str(i) + "\n")
            print(cab_file_name, "SAVED")

            # 记录已经计算的步数
            with open(dir_name + "STATE.dat", "w") as state_file:
                state_file.write(str(i))
            print("CURRENT STATE SAVED")
        if i == release_step:
            k_ryr = 0
            coe_spark, coe_caf = generate_coefficient_matrix(spark_grid, k_ryr)  # 生成系数矩阵

# ***************************************************程序入口***************************************************
def judge(answer_str):
    if answer_str == 'Y' or answer_str == 'y':
        return True
    else:
        return False

def main():
    total_steps = int(input("输入总步数："))
    do_save = int(input("输入保存数据的步数间隔："))
    parameter_info = input("请输入参数情况")
    is_continue = judge(input("是否继上次结束点运行？（Y/N）："))
    if is_continue:
        note = input("请输入文件备注名：")
    else:
        note = "_totalSteps" + str(total_steps) + "_" + time.strftime("%m-%d-%H-%M", time.localtime()) + parameter_info
    dir_name = "..\\Result\\result" + note + "\\"

    initial_parameter()
    loop(total_steps, do_save, dir_name, is_continue)

if __name__ == '__main__':
    main()
