from math import pi

PI = pi

# 网格信息
NR = 1001  # 点总数
DR = 10.0  # 点与点之间的间距，单位nm

# 时间相关，单位：秒
DT = 2 * 10 ** -6  # 时间间隔
RELEASE_TIME = 0.02  # 释放时间，即到此刻ryr通道关闭
STOP_TIME = 0.1  # 结束时间

# 循环控制
SAVE_INTERVAL = 100  # 保存数据间隔
ISTART = 0

# cytosolic spark
CCACYTREST = 1 * 10 ** -4
CCASTORE = 1.0
DCACYT = 3.5 * 10 ** 8
KRYR2 = 1.22522 * 10 ** 10

# Fluo-3 parameter
KPCAF = 80000
KMCAF = 90
BCAF = 0.05
KDCAF = KMCAF / KPCAF
DCAF = 2 * 10 ** 7

# Calmodulin parameter
KPCAM = 100000
# KMCAM = 38
# BCAM = 0.024
KMCAM = 31
BCAM = 0.036
KDCAM = KMCAM / KPCAM

# Troponin C parameter
# KPTRC = 39000
# KMTRC = 20
KPTRC = 125000
KMTRC = 250
BTRC = 0.07
KDTRC = KMTRC / KPTRC

# SR membrane parameter
KPSRM = 115000
KMSRM = 100
BSRM = 0.047
KDSRM = KMSRM / KPSRM

# SL membrane parameter
KPSLM = 115000
KMSLM = 1000
BSLM = 1.124
KDSLM = KMSLM / KPSLM
