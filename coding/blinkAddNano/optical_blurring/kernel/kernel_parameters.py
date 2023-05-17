from math import log

# 用于产生卷积核
# kernel_v1
xy_count = 517  # (卷积核矩阵在x、y维度上点个数-1)/2
z_count = 15  # (卷积核矩阵在z维度上点个数-1)/2
sigma_xy = 400 ** 2 / (8 * log(2))
sigma_z = 800 ** 2 / (8 * log(2))
# kernel_v2 z_count 不好用
# xy_count = 517  # (卷积核矩阵在x、y维度上点个数-1)/2
# z_count = 1033  # (卷积核矩阵在z维度上点个数-1)/2
# sigma_xy = 400 ** 2 / (8 * log(2))
# sigma_z = 800 ** 2 / (8 * log(2))

# kernel_v3 for nanospark
# xy_count = 697  # (卷积核矩阵在x、y维度上点个数-1)/2
# z_count = 697  # (卷积核矩阵在z维度上点个数-1)/2
# sigma_xy = 540 ** 2 / (8 * log(2))
# sigma_z = 540 ** 2 / (8 * log(2))
