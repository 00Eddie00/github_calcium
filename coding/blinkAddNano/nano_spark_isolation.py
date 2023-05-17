import numpy as np
import time
import shutil

from tool_mkdir import mkdir
from nano_calculator_solve2 import nano_core
from parameters import *
# from blink.blink import blink_core
from blink.constant import *


# 循环控制
def loop(total_steps, do_save, dir_name, is_continue):
    nano_c_ca_prefix = f"{dir_name}/NANO/Ca/Ca"
    nano_c_cag_prefix = f"{dir_name}/NANO/CaG/CaG"
    nano_c_caf_prefix = f"{dir_name}/NANO/CaF/CaF"
    nano_c_cab_prefix = f"{dir_name}/NANO/CaB/CaB"
    # blink_c_ca_prefix = f"{dir_name}/BLINK/Ca/Ca"
    # blink_c_caf_prefix = f"{dir_name}/BLINK/CaF/CaF"
    # ryr通道开放时间对应的步数
    release_step = int(RELEASE_TIME / DT)
    # 根据是否断点续算赋不同的值
    if is_continue:
        # 得到已经计算的步数：finished_step；接下来应该计算的步数：current_step
        with open(f"{dir_name}/STATE.csv", 'r') as load_state_file:
            finished_step = int(load_state_file.readline())
        current_step = finished_step + 1
        # 控制ryr通道关闭或开启
        if current_step >= release_step:
            k_ryr = 0
            k_ryr_blink = 0
        else:
            k_ryr = K_RYR
            k_ryr_blink = DCARYR
        # 载入钙离子、荧光钙离子、各缓冲物的浓度
        load_file_number = str(finished_step).zfill(8)
        # 纳米区域
        nano_c_ca = np.loadtxt(f"{nano_c_ca_prefix}{load_file_number}.csv")
        nano_c_cag = np.loadtxt(f"{nano_c_cag_prefix}{load_file_number}.csv")
        nano_c_caf = np.loadtxt(f"{nano_c_caf_prefix}{load_file_number}.csv")
        nano_c_cab = np.loadtxt(f"{nano_c_cab_prefix}{load_file_number}.csv", delimiter=',')
        nano_c_cam = nano_c_cab[:, 0]
        nano_c_trc = nano_c_cab[:, 1]
        nano_c_srm = nano_c_cab[:, 2]
        nano_c_slm = nano_c_cab[:, 3]
        # 开放区域
        # blink_c_ca = np.loadtxt(f"{blink_c_ca_prefix}{load_file_number}.csv")
        # blink_c_caf = np.loadtxt(f"{blink_c_caf_prefix}{load_file_number}.csv")
    else:
        current_step = 1
        # 开启ryr通道
        k_ryr = K_RYR
        k_ryr_blink = DCARYR
        # 初始化各物质浓度
        ini_ca = INITIAL_C_CA
        ini_cag = ini_ca * TOTAL_CAG / (ini_ca + K_DIV_CAG)
        ini_caf = ini_ca * TOTAL_CAF / (ini_ca + K_DIV_CAF)
        ini_cam = ini_ca * TOTAL_CAM / (ini_ca + K_DIV_CAM)
        ini_trc = ini_ca * TOTAL_TRC / (ini_ca + K_DIV_TRC)
        ini_srm = ini_ca * TOTAL_SRM / (ini_ca + K_DIV_SRM)
        ini_slm = ini_ca * TOTAL_SLM / (ini_ca + K_DIV_SLM)

        nano_c_ca = np.full(NANO_POINT_COUNT, ini_ca)
        nano_c_cag = np.full(NANO_POINT_COUNT, ini_cag)
        nano_c_caf = np.full(NANO_POINT_COUNT, ini_caf)
        nano_c_cam = np.full(NANO_POINT_COUNT, ini_cam)
        nano_c_trc = np.full(NANO_POINT_COUNT, ini_trc)
        nano_c_srm = np.full(NANO_POINT_COUNT, ini_srm)
        nano_c_slm = np.full(NANO_POINT_COUNT, ini_slm)

        blink_c_ca = np.full(BLINK_POINT_COUNT, CCAJSR)
        blink_c_caf = np.full(BLINK_POINT_COUNT, CCAF)
        # 纳米区域
        # 记录初始时刻钙离子浓度值
        ca_file_name = f"{nano_c_ca_prefix}00000000.csv"
        np.savetxt(ca_file_name, nano_c_ca)
        print(ca_file_name, "SAVED")
        # 记录初始时刻荧光（GCaMP6f）钙离子浓度值
        cag_file_name = f"{nano_c_cag_prefix}00000000.csv"
        np.savetxt(cag_file_name, nano_c_cag)
        print(cag_file_name, "SAVED")
        # 记录初始时刻荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{nano_c_caf_prefix}00000000.csv"
        np.savetxt(caf_file_name, nano_c_caf)
        print(caf_file_name, "SAVED")
        # JSR：钙空穴
        # 记录初始时刻钙离子浓度值
        # ca_file_name = f"{blink_c_ca_prefix}00000000.csv"
        # np.savetxt(ca_file_name, blink_c_ca)
        # print(ca_file_name, "SAVED")
        # # 记录初始时刻荧光钙离子浓度值
        # caf_file_name = f"{blink_c_caf_prefix}00000000.csv"
        # np.savetxt(caf_file_name, blink_c_caf)
        # print(caf_file_name, "SAVED")

    # 开始计算
    ryr_cca = 0.0001
    for i in range(current_step, total_steps + 1):
        print(f"------------ CURRENT_STEP {i} TOTAL_STEPS {total_steps} ------------")

        # blink_c_ca, blink_c_caf, c_ca_jsr = blink_core(blink_c_ca, blink_c_caf, ryr_cca, k_ryr_blink)
        c_ca_jsr = C_CA_JSR
        nano_c_ca, nano_c_caf, nano_c_cag, nano_c_cam, nano_c_trc, nano_c_srm, nano_c_slm, ryr_cca = \
            nano_core(nano_c_ca, nano_c_cag, nano_c_caf, nano_c_cam, nano_c_trc, nano_c_srm, nano_c_slm, c_ca_jsr,
                      k_ryr)

        # 记录钙离子浓度、荧光钙离子浓度、各缓冲物浓度
        if i % do_save == 0:
            file_number = str(i).zfill(8)
            # 记录钙离子浓度
            ca_file_name = f"{nano_c_ca_prefix}{file_number}.csv"
            np.savetxt(ca_file_name, nano_c_ca)
            print(ca_file_name, "SAVED")
            # 记录初始时刻荧光（GCaMP6f）钙离子浓度值
            cag_file_name = f"{nano_c_cag_prefix}{file_number}.csv"
            np.savetxt(cag_file_name, nano_c_cag)
            print(cag_file_name, "SAVED")
            # 记录初始时刻荧光（Fluo-3）钙离子浓度值
            caf_file_name = f"{nano_c_caf_prefix}{file_number}.csv"
            np.savetxt(caf_file_name, nano_c_caf)
            print(caf_file_name, "SAVED")
            # 记录各缓冲物浓度
            cab_file_name = f"{nano_c_cab_prefix}{file_number}.csv"
            np.savetxt(cab_file_name,
                       np.concatenate((nano_c_cam[:, np.newaxis], nano_c_trc[:, np.newaxis], nano_c_srm[:, np.newaxis],
                                       nano_c_slm[:, np.newaxis]),
                                      axis=1), delimiter=',')
            print(cab_file_name, "SAVED")

            # 记录钙离子浓度
            # ca_file_name = f"{blink_c_ca_prefix}{file_number}.csv"
            # np.savetxt(ca_file_name, blink_c_ca)
            # print(ca_file_name, "SAVED")
            # # 记录初始时刻荧光（Fluo-3）钙离子浓度值
            # caf_file_name = f"{blink_c_caf_prefix}{file_number}.csv"
            # np.savetxt(caf_file_name, blink_c_caf)
            # print(caf_file_name, "SAVED")

            # 记录已经计算的步数
            with open(f"{dir_name}/STATE.csv", "w") as state_file:
                state_file.write(str(i))
            print("CURRENT STATE SAVED")
        # 控制ryr通道关闭
        if i == release_step:
            k_ryr = 0
            k_ryr_blink = 0


# 程序入口
def main(adjust_item, note):
    total_steps = 200
    do_save = 1
    is_continue = False
    result_sets = "../FixedUnstructuredResultSets/"
    if is_continue:
        dir_name = input("请输入文件目录名：")
        dir_name = f"{result_sets}{dir_name}"
    else:
        current_time = time.strftime("%y-%m%d-%H%M", time.localtime())
        dir_name = f"{result_sets}ResultSet_{current_time}_{adjust_item}"
        # 创建用于存储结果的文件夹
        mkdir(f"{dir_name}/NANO/Ca")
        mkdir(f"{dir_name}/NANO/CaG")
        mkdir(f"{dir_name}/NANO/CaF")
        mkdir(f"{dir_name}/NANO/CaB")
        mkdir(f"{dir_name}/BLINK/Ca")
        mkdir(f"{dir_name}/BLINK/CaF")
        # 保存所有参数
        destination = f"{dir_name}/parameters.md"
        shutil.copy("./parameters.py", destination)
        with open(destination, "a", encoding="utf-8") as des_obj:
            des_obj.write("\n# 备注\n")
            des_obj.write(note)

    loop(total_steps, do_save, dir_name, is_continue)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--adjust', type=str, default="")
    # parser.add_argument('--note', type=str, default="")
    # args = parser.parse_args()
    adjust = input()
    notes = input()
    main(adjust, notes)
