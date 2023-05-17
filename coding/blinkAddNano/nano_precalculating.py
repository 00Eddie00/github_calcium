import numpy as np
from math import sqrt
from parameters import B_INFLOW, B_OUTFLOW


# 计算三角形单元相关的数据：结点控制三角形单元的最大数量、每个结点控制三角形单元的总面积、所有三角形单元的总面积
def cal_triangle_area(grid, triangle_matrix):
    triangle_count = len(triangle_matrix)
    point_count = len(grid)
    point_x = grid[:, 0]  # 结点横坐标
    point_y = grid[:, 1]  # 结点纵坐标

    triangle_det = np.empty(triangle_count, float)  # 存储每个三角形单元的面积*2，即公式中的D，Determinant（行列式）
    ctrl_triangle_area = np.zeros(point_count, float)  # 存储每一结点控制三角形单元的总面积
    total_triangle_area = 0.0

    # 辅助变量
    cal_area_matrix = np.empty((3, 3))  # 用来存储一个3*3的矩阵，进而计算三角形单元的面积
    count = np.zeros(point_count)

    for triangle_id in range(0, triangle_count):  # 遍历所有三角形单元
        for vertex_id in range(0, 3):
            cal_area_matrix[vertex_id, 0] = 1.0
            cal_area_matrix[vertex_id, 1] = point_x[triangle_matrix[triangle_id, vertex_id]]
            cal_area_matrix[vertex_id, 2] = point_y[triangle_matrix[triangle_id, vertex_id]]
        triangle_det[triangle_id] = np.linalg.det(cal_area_matrix)  # 1.det()返回传入矩阵的行列式

        area = triangle_det[triangle_id] / 2.0  # 2.此三角形单元的面积
        for vertex_id in range(0, 3):
            ctrl_triangle_area[triangle_matrix[triangle_id, vertex_id]] += area  # 3.每一结点控制三角形单元的总面积
            count[triangle_matrix[triangle_id, vertex_id]] += 1  # 记录每个结点控制三角形单元的个数
        total_triangle_area = total_triangle_area + area  # 4.三角形单元的总面积

    max_ctrl_triangle_count = int(max(count))  # 4.结点控制三角形单元数量的最大值，max()返回传入列表中值最大的元素值

    return triangle_det, ctrl_triangle_area, total_triangle_area, max_ctrl_triangle_count


# 计算与每个三角形单元相关的 a,b,c,N,L
def cal_abc_nl(triangle_matrix, grid, triangle_det):
    triangle_count = len(triangle_matrix)
    point_x = grid[:, 0]  # 结点横坐标
    point_y = grid[:, 1]  # 结点纵坐标

    abc_matrix = np.empty((triangle_count, 3, 3))  # abc[三角形编号][顶点编号][a、b、c]，顶点编号说明此a、b、c从属的顶点
    nl_matrix = np.empty((triangle_count, 3, 3))  # abc[三角形编号][顶点编号][l、nx、ny]，顶点编号说明此顶点作为所求点时需要用到这组l、nx、ny

    for triangle_id in range(0, triangle_count):
        for vertex_id in range(0, 3):
            p2 = triangle_matrix[triangle_id, (vertex_id + 1) % 3]  # 通过取模运算确定点
            p3 = triangle_matrix[triangle_id, (vertex_id + 2) % 3]

            abc_matrix[triangle_id, vertex_id, 0] = (point_x[p2] * point_y[p3] - point_x[p3] * point_y[p2]) / \
                                                    triangle_det[triangle_id]
            abc_matrix[triangle_id, vertex_id, 1] = (point_y[p2] - point_y[p3]) / triangle_det[triangle_id]
            abc_matrix[triangle_id, vertex_id, 2] = (point_x[p3] - point_x[p2]) / triangle_det[triangle_id]

            nl_matrix[triangle_id, vertex_id, 0] = sqrt(
                (point_x[p2] - point_x[p3]) ** 2 + (point_y[p2] - point_y[p3]) ** 2)
            nl_matrix[triangle_id, vertex_id, 1] = (point_y[p3] - point_y[p2]) / nl_matrix[triangle_id, vertex_id, 0]
            nl_matrix[triangle_id, vertex_id, 2] = (point_x[p2] - point_x[p3]) / nl_matrix[triangle_id, vertex_id, 0]

    return abc_matrix, nl_matrix


# 寻找绕此点的三角形单元编号（即每个点控制的三角形单元清单）
def search_triangle(triangle_matrix, point_count, max_ctrl_triangle_count):
    triangle_count = len(triangle_matrix)
    ctrl_triangle_list = np.full([point_count, max_ctrl_triangle_count, 2], -1)  # 存放此点控制的三角形单元

    print("------------ 请稍等，正在搜寻各点控制的三角形单元 ------------")
    for point_id in range(0, point_count):  # 遍历每个结点，point_id 为结点编号
        print(point_id)
        ctrl_triangle_index = 0  # 表示当前结点第几个控制三角形单元
        for triangle_id in range(0, triangle_count):  # 遍历每个三角形单元，triangle_id 为三角形单元编号
            for vertex_id in range(0, 3):  # 遍历三角形的三个顶点，vertex_id 为结点作为三角形单元顶点的编号
                if point_id == triangle_matrix[triangle_id][vertex_id]:
                    ctrl_triangle_list[point_id][ctrl_triangle_index][0] = triangle_id
                    ctrl_triangle_list[point_id][ctrl_triangle_index][1] = vertex_id  # 说明此结点在三角形单元中为几号顶点
                    ctrl_triangle_index = ctrl_triangle_index + 1
                    break
    print("搜寻完毕")
    return ctrl_triangle_list


def det_boundary_type(triangle_matrix, grid, point_och):
    # boundary_type[triangle_id, 0] 指明三角形单元是否存在出入流边界，若无为0
    # boundary_type[triangle_id, 1] 为对应边界长度，若无为0
    triangle_count = len(triangle_matrix)
    boundary_type = np.zeros((triangle_count, 2))
    point_x = grid[:, 0]  # 结点横坐标
    point_y = grid[:, 1]  # 结点纵坐标

    for triangle_id in range(triangle_count):
        in_boundary = []
        out_boundary = []
        for vertex_id in range(0, 3):
            point_id = triangle_matrix[triangle_id, vertex_id]
            if point_och[point_id] == B_INFLOW:
                in_boundary.append(point_id)
            elif point_och[point_id] == B_OUTFLOW:
                out_boundary.append(point_id)

        if len(in_boundary) == 2:
            boundary_type[triangle_id, 0] = B_INFLOW

            boundary_type[triangle_id, 1] = np.sqrt(
                (point_x[in_boundary[0]] - point_x[in_boundary[1]]) ** 2 + (
                        point_y[in_boundary[0]] - point_y[in_boundary[1]]) ** 2)
        elif len(out_boundary) == 2:
            boundary_type[triangle_id, 0] = B_OUTFLOW
            boundary_type[triangle_id, 1] = np.sqrt(
                (point_x[out_boundary[0]] - point_x[out_boundary[1]]) ** 2 + (
                        point_y[out_boundary[0]] - point_y[out_boundary[1]]) ** 2)

    return boundary_type


def cal_out_boundary_length(grid):
    out_boundary = grid[84:205]
    out_boundary_point_count = len(out_boundary)
    out_boundary_length = np.empty(out_boundary_point_count)
    for i in range(out_boundary_point_count):
        if i != out_boundary_point_count - 1:
            out_boundary_length[i] = np.sqrt(
                (out_boundary[i, 0] - out_boundary[i + 1, 0]) ** 2 + (out_boundary[i, 1] - out_boundary[i + 1, 1]) ** 2)
        else:
            out_boundary_length[i] = np.sqrt(
                (out_boundary[i, 0] - out_boundary[0, 0]) ** 2 + (out_boundary[i, 1] - out_boundary[0, 1]) ** 2)
    return out_boundary_length


def precalculating():
    parent_path = "E:\\Code\\Python\\PycharmProjects\\blinkAddNano\\nano_grid_450nm\\"
    pre_parent_path = "E:\\Code\\Python\\PycharmProjects\\blinkAddNano\\nano_grid_450nm\\pre_info\\"
    grid = np.loadtxt(f"{parent_path}gridt.dat", dtype=np.float64)
    triangle_matrix = np.loadtxt(f"{parent_path}nod.dat", dtype=np.int64) - 1
    point_och = np.loadtxt(f"{parent_path}npoch.dat", dtype=np.int64)

    point_count = len(grid)
    triangle_det, ctrl_triangle_area, total_triangle_area, max_ctrl_triangle_count = cal_triangle_area(grid,
                                                                                                       triangle_matrix)
    abc_matrix, nl_matrix = cal_abc_nl(triangle_matrix, grid, triangle_det)
    boundary_type = det_boundary_type(triangle_matrix, grid, point_och)

    # out_boundary_length = cal_out_boundary_length(grid)

    ctrl_triangle_list = search_triangle(triangle_matrix, point_count, max_ctrl_triangle_count)

    print("ctrl_triangle_area SAVED")
    np.save(f"{pre_parent_path}ctrl_triangle_area", ctrl_triangle_area)

    print(f"total_triangle_area: {total_triangle_area}")

    print("abc_matrix SAVED")
    np.save(f"{pre_parent_path}abc_matrix", abc_matrix)

    print("nl_matrix SAVED")
    np.save(f"{pre_parent_path}nl_matrix", nl_matrix)

    print(f"max_ctrl_triangle_count: {max_ctrl_triangle_count}")

    print("ctrl_triangle_list SAVED")
    np.save(f"{pre_parent_path}ctrl_triangle_list", ctrl_triangle_list)

    print("boundary_type SAVED")
    np.save(f"{pre_parent_path}boundary_type", boundary_type)

    # print("out_boundary_length SAVED")
    # np.save("../out_boundary_length", out_boundary_length)


precalculating()
