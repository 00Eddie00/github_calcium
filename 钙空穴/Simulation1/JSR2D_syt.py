import numpy as np
import os

global nod, NP, gridt, npoch, NE, noe, fn_avg, gn_avg, total_area
global single_area, control_area, near_triangle, num_in_triangle, nix_multiply_l, niy_multiply_l, b_arr, c_arr, nmax
global B_INNER, B_WALL, B_INFLOW, B_OUTFLOW, B_SYMMETRY, B_SOURCE
global CCAMYO, CCAFSR, BCSQ, KDCSQ, DCAJSR, DF, DCARYR, DCAFSR, DT, RELEASE_TIMES
global CCAJSR, Gn, NEWCCA, NEWF
global K1, K2, F
global DSAVE, NSTEP, CURSTEP


def init_constant():
    global CCAMYO, CCAFSR, BCSQ, KDCSQ, DCAJSR, DF, DCARYR, DCAFSR, DT, RELEASE_TIMES
    global B_INNER, B_WALL, B_INFLOW, B_OUTFLOW, B_SYMMETRY, B_SOURCE
    global CCAJSR, Gn, NEWCCA, NEWF
    global K1, K2, F
    # Jdiffusion扩散项
    DCAJSR = 3.5 * 10 ** 8  # 终池jSR 自由Ca离子扩散系数
    DF = 2 * 10 ** 7  # 荧光指示剂染料 扩散系数

    # Jbuff
    BCSQ = 14.0  # 方程中的B，所有肌钙集蛋白（肌浆网SR上调控钙储存的主要结合蛋白，调节SR钙释放）的量
    KDCSQ = 0.63  # 方程中的K，肌钙集蛋白的Ca离子解离常数

    # Jrefill,Jryr两项
    # 扩散系数
    DCARYR = 0.25 * 6.5 * 10 ** 7  # RyR通道 Ca离子扩散系数 KRyR 6.5 * 10 ** 7
    DCAFSR = 0.7854 * 10 ** 6  # 连接的fSR通道 Ca离子扩散系数 KfSR
    # 浓度
    CCAMYO = 0.0001  # [Ca2+]cyto  胞浆的Ca浓度
    CCAFSR = 1.0  # [Ca2+]fSR  fSR通道中Ca离子浓度

    # Jdye染料
    F = 0.1  # [F]T   fluo-3总浓度
    K1 = 48800  # KF+  Ca离子和fluo-3的结合速率
    K2 = 19520  # KF-  Ca离子和fluo-3的解离速率
    # K1 = 0
    # K2 = 0

    # 关闭RyR通道用到
    DT = 2 * 10 ** -6  # dt
    RELEASE_TIMES = 2 * 10 ** -2  # 0.02s,这个就是ryr通道开放的时间，后80毫秒是恢复的

    # 点的属性
    B_INNER = 0  # 内部
    B_INFLOW = 2  # 入流
    B_WALL = 1  # 壁面
    B_OUTFLOW = 4  # 出流
    B_SYMMETRY = 8
    B_SOURCE = 16

    CCAJSR = np.ones(NP, float)  # fn,肌质网每一点钙的浓度
    Gn = np.ones(NP, float) / 14.0  # gn初始值 1/14  固定
    NEWCCA = np.empty(NP, float)  # 存每一次迭代的中间结果fn
    NEWF = np.empty(NP, float)  # 存每一次迭代的中间结果gn


def load_files():
    global nod, NP, gridt, npoch, NE, noe
    grid_dir = "./GRID"
    i = 0
    with open(grid_dir + "/gridt.dat") as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            i = i + 1
    NP = i
    file_obj.close()

    i = 0
    gridt = np.zeros([2, NP], float)
    with open(grid_dir + "/gridt.dat") as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            line_list = line.split()
            gridt[0, i] = float(line_list[0])  # x
            gridt[1, i] = float(line_list[1])  # y
            i = i + 1
    file_obj.close()

    i = 0
    npoch = np.zeros(NP, int)
    with open(grid_dir + "/npoch.dat") as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            line_list = line.split()
            npoch[i] = int(line_list[0])
            i = i + 1
    file_obj.close()

    i = 0
    with open(grid_dir + "/nod.dat") as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            i = i + 1
    NE = i
    file_obj.close()

    i = 0
    nod = np.zeros([3, NE], int)
    with open(grid_dir + "/nod.dat") as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            line_list = line.split()
            for j in range(0, 3):
                nod[j, i] = int(line_list[j]) - 1
            i = i + 1
    file_obj.close()

    i = 0
    noe = np.zeros([3, NE], int)
    with open(grid_dir + "/noe.dat") as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            line_list = line.split()
            for j in range(0, 3):
                noe[j, i] = int(line_list[j]) - 1
            i = i + 1
    file_obj.close()


def cal_elements():
    global single_area, control_area, near_triangle, num_in_triangle, nix_multiply_l, niy_multiply_l, b_arr, c_arr, nmax, total_area
    single_area = np.zeros(NE, float)
    control_area = np.zeros(NP, float)
    cal_max = np.zeros(NP, int)
    area_matrix = np.zeros([3, 3], float)
    nix_multiply_l = np.zeros([3, NE], float)
    niy_multiply_l = np.zeros([3, NE], float)
    b_arr = np.zeros([3, NE], float)
    c_arr = np.zeros([3, NE], float)
    total_area = 0.0
    nmax = 0
    for i in range(0, NE):
        for j in range(0, 3):
            p = nod[j, i]
            area_matrix[j, 0] = 1.0
            area_matrix[j, 1] = gridt[0, p]
            area_matrix[j, 2] = gridt[1, p]
        area = np.linalg.det(area_matrix) * 0.5
        single_area[i] = area
        total_area = total_area + area
        for v in range(0, 3):
            p = nod[v, i]
            control_area[p] = control_area[p] + area
            cal_max[p] = cal_max[p] + 1
    nmax = int(max(cal_max))
    near_triangle = np.full([nmax, NP], -1)
    num_in_triangle = np.full([nmax, NP], -1)
    num_temp = np.zeros(NP, int)
    for j in range(0, NE):
        for k in range(0, 3):
            i = nod[k, j]
            t = num_temp[i]
            near_triangle[t, i] = j
            num_in_triangle[t, i] = k
            num_temp[i] = num_temp[i] + 1

            nod2 = nod[(k + 1) % 3, j]
            nod3 = nod[(k + 2) % 3, j]
            nix_multiply_l[k, j] = gridt[1, nod2] - gridt[1, nod3]
            niy_multiply_l[k, j] = gridt[0, nod3] - gridt[0, nod2]
            b_arr[k, j] = (gridt[1, nod3] - gridt[1, nod2]) / (2 * single_area[j])
            c_arr[k, j] = (gridt[0, nod2] - gridt[0, nod3]) / (2 * single_area[j])


def judge_point(n):
    in_boundary, out_boundary, inner = [], [], []
    for j in range(0, 3):
        k = nod[j, n]
        if npoch[k] == B_INFLOW:
            in_boundary.append(k)
        elif npoch[k] == B_OUTFLOW:
            out_boundary.append(k)
        elif npoch[k] == B_INNER:
            inner.append(k)
    return in_boundary, out_boundary, inner


def fn(iteration):
    global NEWCCA
    tempcca = np.empty(NP, float)  # 保存每一次迭代后的结果
    for it in range(0, iteration):
        for p in range(0, NP):
            m1, m5, m6 = 0.0, 0.0, 0.0  # 分别对应(1)(5)(6)式
            m2 = DT * (K2 * Gn[p] - K1 * CCAJSR[p] * (F - Gn[p])) * control_area[p]  # 计算（2）式
            n7, n8, n9, n10 = 0.0, 0.0, 0.0, 0.0  # 分别对应(7)(8)(9)(10)式
            in_boundary_list = []  # 保存有入流边界的三角形编号
            out_boundary_list = []  # 保存有出流边界的三角形编号
            for o in range(0, nmax):  # 计算相邻三角形
                n = near_triangle[o, p]  # n为三角形编号
                if n == -1:
                    break
                i = num_in_triangle[o, p]
                nod_2, nod_3 = nod[(i + 1) % 3, n], nod[(i + 2) % 3, n]  # 其余两点的编号及在三角形中的序号
                f1i, f2i, f3i = CCAJSR[p], CCAJSR[nod_2], CCAJSR[nod_3]  # f1i, f2i, f3i取n时刻的值
                if it > 0:  # 判断是否为第一次迭代
                    f1i, f2i, f3i = (CCAJSR[p] + NEWCCA[p]) / 2.0, (CCAJSR[nod_2] + NEWCCA[nod_2]) / 2.0, (
                            CCAJSR[nod_3] + NEWCCA[nod_3]) / 2.0
                avg = (1 / 3) * (f1i + f2i + f3i)  # fni为第 i 个三角形单元的三个顶点 n 时刻的平均值
                m5 = m5 + (1.0 + ((BCSQ * KDCSQ) / (KDCSQ + avg) ** 2)) * single_area[n]  # 计算（5）式
                in_boundary, out_boundary, inner = judge_point(n)
                if len(in_boundary) == 2:
                    in_boundary_list.append(o)  # 记录该点的triangleNumber的下标
                elif len(out_boundary) == 2:
                    out_boundary_list.append(o)
                if npoch[nod_2] > B_INNER and npoch[nod_3] > B_INNER:
                    continue
                m1 = m1 + (nix_multiply_l[i, n] * (f2i * b_arr[(i + 1) % 3, n] + f3i * b_arr[(i + 2) % 3, n]) +
                           niy_multiply_l[i, n] * (f2i * c_arr[(i + 1) % 3, n] + f3i * c_arr[(i + 2) % 3, n]))  # 计算（1）式
                m6 = m6 + (nix_multiply_l[i, n] * b_arr[i, n] + niy_multiply_l[i, n] * c_arr[i, n])  # 计算（6）式
            tempcca[p] = (DT * DCAJSR * m1 + m2 + CCAJSR[p] * m5 + 0.5 * CCAJSR[p] * DT * DCAJSR * m6) / (
                    m5 - 0.5 * DT * DCAJSR * m6)  # 存入临时数组保存
            if len(in_boundary_list) > 0:  # 入流边界数目大于0
                for q in range(0, len(in_boundary_list)):
                    o1 = in_boundary_list[q]
                    in_n = near_triangle[o1, p]  # 取出三角形编号
                    i = num_in_triangle[o1, p]
                    in_boundary, out_boundary, inner = judge_point(in_n)
                    a = in_boundary[0]
                    b = in_boundary[1]
                    nod_2, nod_3 = nod[(i + 1) % 3, in_n], nod[(i + 2) % 3, in_n]
                    f2j, f3j = CCAJSR[nod_2], CCAJSR[nod_3]
                    if it > 0:  # 判断是否为第一次迭代
                        f2j, f3j = NEWCCA[nod_2], NEWCCA[nod_3]
                    Lj = np.sqrt((gridt[0, a] - gridt[0, b]) ** 2 + (gridt[1, a] - gridt[1, b]) ** 2)  # 入流边界长度
                    n7 = n7 + DCAFSR * (CCAFSR - (f2j + f3j) / 3.0) * Lj  # 计算（7）式
                    n8 = n8 + DCAFSR * Lj  # 计算（8）式
                tempcca[p] = (DT * DCAJSR * m1 + m2 + CCAJSR[p] * m5 + 0.5 * CCAJSR[p] * DT * DCAJSR * m6 + n7 * DT) / (
                        m5 - 0.5 * DT * DCAJSR * m6 + (n8 * DT) / 3.0)  # 若有入流边界则修改临时数组
            elif len(out_boundary_list) > 0:  # 出流边界数目大于0
                for r in range(0, len(out_boundary_list)):
                    o2 = out_boundary_list[r]
                    out_n = near_triangle[o2, p]  # 三角形编号
                    i = num_in_triangle[o2, p]
                    in_boundary, out_boundary, inner = judge_point(out_n)
                    a = out_boundary[0]
                    b = out_boundary[1]
                    nod_2, nod_3 = nod[(i + 1) % 3, out_n], nod[(i + 2) % 3, out_n]
                    f2k, f3k = CCAJSR[nod_2], CCAJSR[nod_3]
                    if it > 0:
                        f2k, f3k = NEWCCA[nod_2], NEWCCA[nod_3]
                    Lk = np.sqrt((gridt[0, a] - gridt[0, b]) ** 2 + (gridt[1, a] - gridt[1, b]) ** 2)  # 出流边界长度
                    n9 = n9 + DCARYR * (CCAMYO - (f2k + f3k) / 3.0) * Lk  # 计算（9）式
                    n10 = n10 + DCARYR * Lk  # 计算（10）式
                tempcca[p] = (DT * DCAJSR * m1 + m2 + CCAJSR[p] * m5 + 0.5 * CCAJSR[p] * DT * DCAJSR * m6 + n9 * DT) / (
                        m5 - 0.5 * DT * DCAJSR * m6 + (n10 * DT) / 3.0)  # 若有出流边界则修改临时数组
        for t in range(0, NP):
            NEWCCA[t] = tempcca[t]  # 保存这次迭代的结果
            tempcca[t] = 0.0


def gn(iteration):
    global NEWF
    tempf = np.empty(NP, float)  # 保存每一次迭代后的结果
    for it in range(0, iteration):
        for p in range(0, NP):
            m1, m2, m6 = 0.0, 0.0, 0.0,  # 分别对应(1)(2)(6)式
            for o in range(0, nmax):
                n = near_triangle[o, p]  # 三角形编号
                if n == -1:
                    break
                i = num_in_triangle[o, p]
                nod_2, nod_3 = nod[(i + 1) % 3, n], nod[(i + 2) % 3, n]
                if npoch[nod_2] > B_INNER and npoch[nod_3] > B_INNER:
                    continue
                g2i, g3i = Gn[nod_2], Gn[nod_3]
                if it > 0:  # 判断是否为第一次迭代
                    g2i, g3i = (Gn[nod_2] + NEWF[nod_2]) / 2.0, (Gn[nod_3] + NEWF[nod_3]) / 2.0
                m1 = m1 + DT * DF * (
                        nix_multiply_l[i, n] * (g2i * b_arr[(i + 1) % 3, n] + g3i * b_arr[(i + 2) % 3, n]) +
                        niy_multiply_l[i, n] * (
                                g2i * c_arr[(i + 1) % 3, n] + g3i * c_arr[(i + 2) % 3, n]))
                m6 = m6 + DT * DF * 0.5 * (nix_multiply_l[i, n] * b_arr[i, n] + niy_multiply_l[i, n] * c_arr[i, n])
            m2 = DT * (K1 * CCAJSR[p] * (F - Gn[p]) - K2 * Gn[p]) * control_area[p]
            tempf[p] = (m1 + m2 + Gn[p] * control_area[p] + m6 * Gn[p]) / (
                    control_area[p] - m6)
        for t in range(0, NP):
            NEWF[t] = tempf[t]  # 保存上这次迭代的结果
            tempf[t] = 0.0


def doloop_diffusion():
    global DSAVE, NSTEP, CURSTEP, RELEASE_TIMES, fn_avg, gn_avg, CCAJSR, Gn, DCARYR, NEWF, NEWCCA
    CURSTEP = 1  # 从多少步开始，默认从第1步开始
    DSAVE = 1  # 每做多少步保存一次
    NSTEP = 50000  # 程序到多少步结束
    ITERATION = 10  # 迭代次数
    PATH = "DATA\\DATA_blink"  # 保存路径
    pathfn = PATH + '\\Fn'  # DATA\\DATA_blink\\Fn
    pathgn = PATH + '\\Gn'  # DATA\\DATA_blink\\Gn
    if not os.path.exists(pathfn):
        os.makedirs(pathfn)
    if not os.path.exists(pathgn):
        os.makedirs(pathgn)
    if CURSTEP == 1:
        cal_avg()
        fn_file_name = "Fn00000000.dat"
        write_fn = os.path.join(pathfn, fn_file_name)
        file = open(write_fn, 'w')
        for j in range(0, NP):
            file.write(str(CCAJSR[j]) + '\n')
        file.write(str(fn_avg) + '\n')
        print("数据写入 ", fn_file_name)
        file.flush()
        file.close()
        gn_file_name = "Gn00000000.dat"
        write_gn = os.path.join(pathgn, gn_file_name)
        file = open(write_gn, 'w')
        for j in range(0, NP):
            file.write(str(Gn[j]) + '\n')
        file.write(str(gn_avg) + '\n')
        print("数据写入 ", gn_file_name)
        file.flush()
        file.close()
    else:
        length = len(str(CURSTEP - 1))
        rest_name = "0" * (8 - length) + str(CURSTEP - 1) + ".dat"
        fn_file_name = "Fn" + rest_name
        gn_file_name = "Gn" + rest_name

        i = 0
        data_dir = "./DATA/DATA_blink"
        with open(data_dir + "/Fn/" + fn_file_name) as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                line_list = line.split()
                if i < NP:
                    CCAJSR[i] = float(line_list[0])
                else:
                    fn_avg = float(line_list[0])
                i = i + 1
        file_obj.close()
        i = 0
        with open(data_dir + "/Gn/" + gn_file_name) as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                line_list = line.split()
                if i < NP:
                    Gn[i] = float(line_list[0])
                else:
                    gn_avg = float(line_list[0])
                i = i + 1
        file_obj.close()

    RELEASE_STEP = int((RELEASE_TIMES + (10 ** -6)) / DT)  # 释放时间的迭代次数 RELEASE_TIMES=2 * 10 ** -2s DT=2 * 10 ** -6s
    if CURSTEP >= RELEASE_STEP:
        DCARYR = 0
        # 释放的点就停止释放

    for j in range(CURSTEP, NSTEP + 1):  # NSTEP = 50000程序到多少步结束
        # 每一步需要迭代 5-10 次
        fn(ITERATION)  # ITERATION = 10 迭代次数
        gn(ITERATION)
        for i in range(0, NP):
            CCAJSR[i] = NEWCCA[i]
            Gn[i] = NEWF[i]
        if (j % DSAVE == 0) or (j == NSTEP):
            cal_avg()
            length = len(str(j))
            rest_name = "0" * (8 - length) + str(j) + ".dat"
            fn_file_name = "Fn" + rest_name
            gn_file_name = "Gn" + rest_name
            write_fn = os.path.join(pathfn, fn_file_name)
            file = open(write_fn, 'w')
            for k in range(0, NP):
                file.write(str(CCAJSR[k]) + '\n')
            file.write(str(fn_avg) + '\n')
            print("数据写入 ", fn_file_name)
            file.flush()
            file.close()
            write_gn = os.path.join(pathgn, gn_file_name)
            file = open(write_gn, 'w')
            for k in range(0, NP):
                file.write(str(Gn[k]) + '\n')
            file.write(str(gn_avg) + '\n')
            print("数据写入 ", gn_file_name)
            file.flush()
            file.close()
        if j == RELEASE_STEP:
            DCARYR = 0
        for i in range(0, NP):
            NEWCCA[i] = 0.
            NEWF[i] = 0.


def cal_avg():
    global fn_avg, gn_avg
    fn_total, gn_total = 0.0, 0.0
    for i in range(0, NP):
        fn_total = fn_total + CCAJSR[i] * control_area[i]
        gn_total = gn_total + Gn[i] * control_area[i]
    fn_total = fn_total / 3.0
    gn_total = gn_total / 3.0
    fn_avg = fn_total / total_area
    gn_avg = gn_total / total_area


def main():
    load_files()
    # print("load_files执行\n")
    init_constant()
    # print("init_constant执行\n")
    cal_elements()
    # print("cal_elements执行\n")
    doloop_diffusion()


main()
