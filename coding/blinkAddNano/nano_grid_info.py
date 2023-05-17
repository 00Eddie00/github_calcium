import numpy as np


# 纳米空间网格信息
# 各结点坐标
parent_path = "E:\\Code\\Python\\PycharmProjects\\blinkAddNano\\blinkgrid\\"
grid = np.loadtxt(parent_path + "gridt.dat", dtype=np.float64)
point_x = grid[:, 0]  # 点的横坐标
point_y = grid[:, 1]  # 点的纵坐标
nano_point_count = len(grid)
# 三角单元组成结点
triangle_matrix = np.loadtxt(parent_path + "nod.dat", dtype=np.int64) - 1
triangle_count = len(triangle_matrix)
# 三角单元邻居
triangle_neighbor = np.loadtxt(parent_path + "noe.dat", dtype=np.int64) - 1
# 各节点属性
point_och = np.loadtxt(parent_path + "npoch_nano.dat", dtype=np.int64)
# 其它信息
pre_parent_path = "E:\\Code\\Python\\PycharmProjects\\blinkAddNano\\nano_grid\\blink_mirrior\\"
boundary_type = np.load(pre_parent_path + "boundary_type.npy")
abc_matrix = np.load(pre_parent_path + "abc_matrix.npy")
nl_matrix = np.load(pre_parent_path + "nl_matrix.npy")
ctrl_triangle_list = np.load(pre_parent_path + "ctrl_triangle_list.npy")
ctrl_triangle_area = np.load(pre_parent_path + "ctrl_triangle_area.npy")
# out_boundary_length = np.load(pre_parent_path + "out_boundary_length.npy")
max_ctrl_triangle_count = 8
# total_triangle_area=282532.50412262074
