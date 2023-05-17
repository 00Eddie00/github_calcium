import numpy as np
from config.parameters import PI, NR, DT, DCACYT, DCAF, CCASTORE, \
    KMCAF, BCAF, KPCAF, \
    KMCAM, BCAM, KPCAM, \
    KMTRC, BTRC, KPTRC, \
    KMSLM, BSLM, KPSLM, \
    KMSRM, BSRM, KPSRM


# 计算J_dye，传入ca浓度和荧光材料与钙离子结合的浓度
# KPCAF为钙离子与染料结合的结合速率常数，KMCAF为钙离子与染料离解的离解速率常数
def cal_dye(c_ca, c_caf):
    j_dye = -KPCAF * c_ca * (BCAF - c_caf) + KMCAF * c_caf
    return j_dye


# 计算J_i，传入的为各缓冲物与钙离子结合产物的浓度
# KPCAM，KPTRC，KPSRM，KPSLM为钙离子与各钙离子缓冲物结合的结合速率常数
# KMCAM，KMTRC，KMSRM，KMSLM为钙离子与各钙离子缓冲物离解的离解速率常数
# BCAM，BTRC，BSRM，BSLM为各钙离子缓冲物的总浓度
def cal_j_i(c_ca, c_cam, c_trc, c_srm, c_slm):
    j_cam = -KPCAM * c_ca * (BCAM - c_cam) + KMCAM * c_cam
    j_trc = -KPTRC * c_ca * (BTRC - c_trc) + KMTRC * c_trc
    j_srm = -KPSRM * c_ca * (BSRM - c_srm) + KMSRM * c_srm
    j_slm = -KPSLM * c_ca * (BSLM - c_slm) + KMSLM * c_slm
    return j_cam, j_trc, j_srm, j_slm


# 计算J_buffers
def cal_j_buffers(c_ca, c_cam, c_trc, c_srm, c_slm):
    j_cam, j_trc, j_srm, j_slm = cal_j_i(c_ca, c_cam, c_trc, c_srm, c_slm)
    j_buffers = j_cam + j_trc + j_srm + j_slm
    return j_buffers


# 产生系数矩阵，即文档中的A矩阵
# k_ryr 应该为ryr通道处钙离子释放的扩散系数，我觉得！！！！！！！！！！
def generate_coefficient_matrix(spark_grid, k_ryr):
    coe_spark = np.full((NR, NR), 0, dtype = float)
    coe_caf = np.full((NR, NR), 0, dtype = float)
    for j in range(0, NR):
        if j == 0:
            r2 = spark_grid[0]
            r3 = spark_grid[1]
            r2_5 = (r2 + r3) / 2.0
            v = (4 * PI * r2_5 ** 3) / 3.0
            f_A11 = 3 * DT * DCACYT + r2_5 * r3 + (k_ryr / v) * DT * r2_5 * r3
            f_A12 = -3 * DT * DCACYT
            coe_spark[0][0] = f_A11
            coe_spark[0][1] = f_A12
            g_A11 = 3 * DT * DCAF + r2_5 * r3
            g_A12 = -3 * DT * DCAF
            coe_caf[0][0] = g_A11
            coe_caf[0][1] = g_A12
        elif j == NR - 1:
            r1 = spark_grid[NR - 2]
            r2 = spark_grid[NR - 1]
            r1_5 = (r1 + r2) / 2.0
            f_Akk_1 = -3 * DT * DCACYT * r1_5 ** 2
            f_Akk = 3 * DT * DCACYT * r1_5 ** 2 + (r2 ** 3 - r1_5 ** 3) * (r2 - r1)
            coe_spark[j][j-1] = f_Akk_1
            coe_spark[j][j] = f_Akk
            g_Akk_1 = -3 * DT * DCAF * r1_5 ** 2
            g_Akk = 3 * DT * DCAF * r1_5 ** 2 + (r2 ** 3 - r1_5 ** 3) * (r2 - r1)
            coe_caf[j][j - 1] = g_Akk_1
            coe_caf[j][j] = g_Akk
        else:
            r1 = spark_grid[j - 1]
            r2 = spark_grid[j]
            r3 = spark_grid[j + 1]
            r2_5 = (r2 + r3) / 2.0
            r1_5 = (r1 + r2) / 2.0
            f_A1 = -3 * DT * DCACYT * r1_5 ** 2 * (r3 - r2)
            f_A2 = 3 * DT * DCACYT * (r2_5 ** 2 * (r2 - r1) + r1_5 ** 2 * (r3 - r2)) + (r2_5 ** 3 - r1_5 ** 3) * (r3 - r2) * (r2 - r1)
            f_A3 = -3 * DT * DCACYT * r2_5 ** 2 * (r2 - r1)
            # print(f_A2)
            coe_spark[j][j-1] = f_A1
            coe_spark[j][j] = f_A2
            coe_spark[j][j + 1] = f_A3
            g_A1 = -3 * DT * DCAF * r1_5 ** 2 * (r3 - r2)
            g_A2 = 3 * DT * DCAF * (r2_5 ** 2 * (r2 - r1) + r1_5 ** 2 * (r3 - r2)) + (r2_5 ** 3 - r1_5 ** 3) * (r3 - r2) * (r2 - r1)
            g_A3 = -3 * DT * DCAF * r2_5 ** 2 * (r2 - r1)
            coe_caf[j][j-1] = g_A1
            coe_caf[j][j] = g_A2
            coe_caf[j][j + 1] = g_A3
    return coe_spark, coe_caf


# 产生常数矩阵，即文档中的B矩阵
def generate_constant_matrix(c_ca, c_caf, c_cam, c_trc, c_srm, c_slm, spark_grid, k_ryr):
    const_spark = np.full((NR), 0, dtype = float)
    const_caf = np.full((NR), 0, dtype = float)
    J_dye = - KPCAF * c_ca * (BCAF - c_caf) + KMCAF * c_caf
    J_buffers = cal_j_buffers(c_ca, c_cam, c_trc, c_srm, c_slm)

    for j in range(0, NR):
        if j == 0:
            r2 = spark_grid[0]
            r3 = spark_grid[1]
            r2_5 = (r2 + r3) / 2.0
            v = (4 * PI * r2_5 ** 3) / 3.0
            f_B = r2_5 * r3 * ((J_dye[j] + J_buffers[j] + (k_ryr / v) * CCASTORE) * DT + c_ca[0]) #CCASTORE是[Ca2+]jSR吗？
            const_spark[0] = f_B
            g_B = r2_5 * r3 * ((- J_dye[j]) * DT + c_caf[0])
            const_caf[0] = g_B
        elif j == NR - 1:
            r1 = spark_grid[NR - 2]
            r2 = spark_grid[NR - 1]
            r1_5 = (r1 + r2) / 2.0

            f_B = (r2 ** 3 - r1_5 ** 3) * (r2 - r1) * ((J_dye[j] + J_buffers[j]) * DT + c_ca[j])
            const_spark[j] = f_B
            g_B = (r2 ** 3 - r1_5 ** 3) * (r2 - r1) * ((- J_dye[j]) * DT + c_caf[j])
            const_caf[j] = g_B
        else:
            r1 = spark_grid[j - 1]
            r2 = spark_grid[j]
            r3 = spark_grid[j + 1]
            r2_5 = (r2 + r3) / 2.0
            r1_5 = (r1 + r2) / 2.0
            f_B = (r2_5 ** 3 - r1_5 ** 3) * (r3 - r2) * (r2 - r1) * ((J_dye[j] + J_buffers[j]) * DT + c_ca[j])
            const_spark[j] = f_B
            g_B = (r2_5 ** 3 - r1_5 ** 3) * (r3 - r2) * (r2 - r1) * ((- J_dye[j]) * DT + c_caf[j])
            const_caf[j] = g_B

    return const_spark, const_caf


# 求解矩阵方程，得到Ca浓度、CaF浓度
def cal_f_g_concentration(coefficient_matrix, constant_matrix):
    concentration_arr = np.linalg.solve(coefficient_matrix, constant_matrix)
    return concentration_arr


# 计算缓冲物浓度
def cal_h_concentration(c_ca, c_cam, c_trc, c_srm, c_slm):
    j_cam, j_trc, j_srm, j_slm = cal_j_i(c_ca, c_cam, c_trc, c_srm, c_slm)

    c_cam = c_cam + (-j_cam) * DT
    c_trc = c_trc + (-j_trc) * DT
    c_srm = c_srm + (-j_srm) * DT
    c_slm = c_slm + (-j_slm) * DT

    return c_cam, c_trc, c_srm, c_slm


# 计算加权平均值
def cal_avg_concentration(spark_grid, c_matrix):
    total_c = 0.0
    total_v = (4 * PI * spark_grid[NR - 1] ** 3) / 3.0

    for j in range(0, NR):
        if j == 0:
            r3 = (spark_grid[0] + spark_grid[1]) / 2.0
            r1 = 0
        elif j == NR - 1:
            r3 = spark_grid[NR - 1]
            r1 = (spark_grid[NR - 1] + spark_grid[NR - 2]) / 2.0
        else:
            r3 = (spark_grid[j] + spark_grid[j + 1]) / 2.0
            r1 = (spark_grid[j - 1] + spark_grid[j]) / 2.0

        ctrl_v = (4 * PI / 3.0) * (r3 ** 3 - r1 ** 3)
        total_c = total_c + c_matrix[j] * ctrl_v

    avg_c = total_c / total_v
    return avg_c
