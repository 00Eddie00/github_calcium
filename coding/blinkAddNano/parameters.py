from numpy import pi

# 网格相关
NANO_POINT_COUNT = 6413
BLINK_POINT_COUNT = 6413

# 时间相关，单位：秒
RELEASE_TIME = 0.02  # 释放时间，即到此刻RYR通道关闭
STOP_TIME = 0.1  # 结束时间

# 出入流点属性值
B_INNER = 0  # 内部
B_INFLOW = 2  # 入流
B_OUTFLOW = 4  # 出流

# 入流参数
C_CA_JSR = 1.0  # JSR中的钙离子浓度
S = 2 * pi * 0.5 * 15
# K_RYR = 1.24368 * 10 ** 10 / S
K_RYR = 2 * 6.5 * 10 ** 7

# 出流参数
delta_r = 1.0
K_OUT = 7.854 * 10 ** 6
# K_FOUT = 2 * 10 ** 7 / delta_r  # 出流扩散系数
# C_CAF_OUT = 1 / 245
C_CA_OUT = 0.0001

# 时间单位s，长度单位nm，浓度mM
# 时间相关，单位：秒
DT = 2 * 10 ** -6  # DT 时间间隔

# 迭代次数
ITERATION_TIMES = 10

# cytosolic spark
INITIAL_C_CA = 0.0001  # 肌质中钙离子初始浓度
D_CA = 3.5 * 10 ** 8  # 自由钙离子的扩散系数

# 钙离子缓冲物的相应系数
# 荧光染料 GCaMP6f
K_PLUS_CAG = 27000  # 50000
K_MINUS_CAG = 17  # 34
TOTAL_CAG = 0.01
K_DIV_CAG = K_MINUS_CAG / K_PLUS_CAG

# 荧光染料 Fluo-3
K_PLUS_CAF = 80000
K_MINUS_CAF = 90
TOTAL_CAF = 0.05
K_DIV_CAF = K_MINUS_CAF / K_PLUS_CAF
D_CAF = 2 * 10 ** 7

# # calmodulin
# K_PLUS_CAM = 100000
# K_MINUS_CAM = 31
# TOTAL_CAM = 0.036
# K_DIV_CAM = K_MINUS_CAM / K_PLUS_CAM
#
# # troponin c
# K_PLUS_TRC = 125000
# K_MINUS_TRC = 250
# TOTAL_TRC = 0.07
# K_DIV_TRC = K_MINUS_TRC / K_PLUS_TRC
#
# # sr membrane
# # 在 KD 增大，K-不变的情况下，K+会变小
# K_PLUS_SRM = 115000
# K_MINUS_SRM = 100
# TOTAL_SRM = 0.047   # 改为 13
# K_DIV_SRM = K_MINUS_SRM / K_PLUS_SRM  # 改为 0.013
#
# # sl membrane
# K_PLUS_SLM = 115000
# K_MINUS_SLM = 1000
# TOTAL_SLM = 1.124   # 改为 165
# K_DIV_SLM = K_MINUS_SLM / K_PLUS_SLM  # 改为 1.1


# 缓冲物参数替换为 1998 Peskoff 论文 Table 1
# calmodulin
K_PLUS_CAM = 100000
K_MINUS_CAM = 38
TOTAL_CAM = 0.024
K_DIV_CAM = K_MINUS_CAM / K_PLUS_CAM

# troponin c
K_PLUS_TRC = 39000
K_MINUS_TRC = 20
TOTAL_TRC = 0.07
K_DIV_TRC = K_MINUS_TRC / K_PLUS_TRC

# sr membrane
# 在 KD 增大，K-不变的情况下，K+会变小
K_MINUS_SRM = 100
TOTAL_SRM = 13
K_DIV_SRM = 0.013
K_PLUS_SRM = K_MINUS_SRM / K_DIV_SRM

# sl membrane
K_MINUS_SLM = 1000
TOTAL_SLM = 165
K_DIV_SLM = 1.1
K_PLUS_SLM = K_MINUS_SLM / K_DIV_SLM
