from numpy import pi

# 网格相关
POINT_TOTALS = 19001

# 时间相关，单位：秒
DT = 2 * 10 ** -6  # DT 时间间隔
RELEASE_TIME = 0.02  # 释放时间，即到此刻RYR通道关闭
STOP_TIME = 0.1  # 结束时间

# 迭代次数
ITERATION_TIMES = 10

# cytosolic spark
INITIAL_C_CA = 0.0001  # 肌质中钙离子初始浓度
C_CA_JSR = 1.0  # JSR中的钙离子浓度
D_CA = 3.5 * 10 ** 8  # 自由钙离子的扩散系数
V = pi * 5 ** 2 * 0.5  # 终池靠近T小管膜的一侧有一个直径为10nm的RyR通道
K_RYR = 1.22522 * 10 ** 10 / V  # RYR通道处钙离子释放的扩散系数

# 钙离子缓冲物的相应系数
# 荧光染料 GCaMP6f-J
K_PLUS_CAF = 27000
# K_PLUS_CAF = 26899
K_MINUS_CAF = 17
TOTAL_CAF = 0.01
K_DIV_CAF = K_MINUS_CAF / K_PLUS_CAF  # KF

# calmodulin
K_PLUS_CAM = 100000
K_MINUS_CAM = 31
TOTAL_CAM = 0.036
K_DIV_CAM = K_MINUS_CAM / K_PLUS_CAM

# troponin c
K_PLUS_TRC = 125000
K_MINUS_TRC = 250
TOTAL_TRC = 0.07
K_DIV_TRC = K_MINUS_TRC / K_PLUS_TRC

# sr membrane
K_PLUS_SRM = 115000
K_MINUS_SRM = 100
TOTAL_SRM = 0.047
K_DIV_SRM = K_MINUS_SRM / K_PLUS_SRM

# sl membrane
K_PLUS_SLM = 115000
K_MINUS_SLM = 1000
TOTAL_SLM = 1.124
K_DIV_SLM = K_MINUS_SLM / K_PLUS_SLM
