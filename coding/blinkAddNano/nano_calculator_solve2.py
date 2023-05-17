from parameters import *
from nano_grid_info import *


def cal_j_dye_g(c_ca, c_cag):
    j_dye_g = -K_PLUS_CAG * c_ca * (TOTAL_CAG - c_cag) + K_MINUS_CAG * c_cag
    return j_dye_g


def cal_j_dye_f(c_ca, c_caf):
    j_dye_f = -K_PLUS_CAF * c_ca * (TOTAL_CAF - c_caf) + K_MINUS_CAF * c_caf
    return j_dye_f


def cal_j_buffer(c_ca, c_cam, c_trc, c_srm, c_slm):
    j_cam = -K_PLUS_CAM * c_ca * (TOTAL_CAM - c_cam) + K_MINUS_CAM * c_cam
    j_trc = -K_PLUS_TRC * c_ca * (TOTAL_TRC - c_trc) + K_MINUS_TRC * c_trc
    j_srm = -K_PLUS_SRM * c_ca * (TOTAL_SRM - c_srm) + K_MINUS_SRM * c_srm
    j_slm = -K_PLUS_SLM * c_ca * (TOTAL_SLM - c_slm) + K_MINUS_SLM * c_slm
    return j_cam, j_trc, j_srm, j_slm


def cal_ca_concentration(c_ca, c_cag, c_cam, c_trc, c_srm, c_slm, c_ca_jsr, k_ryr):
    j_dye_g = cal_j_dye_g(c_ca, c_cag)
    # j_dye_f = cal_j_dye_f(c_ca, c_caf)
    j_cam, j_trc, j_srm, j_slm = cal_j_buffer(c_ca, c_cam, c_trc, c_srm, c_slm)
    j_buffers = j_cam + j_trc + j_srm + j_slm
    # j_item = j_dye_g + j_dye_f + j_buffers
    j_item = j_dye_g + j_buffers

    new_c_cag = DT * (-j_dye_g) + c_cag
    new_c_cam = DT * (-j_cam) + c_cam
    new_c_trc = DT * (-j_trc) + c_trc
    new_c_srm = DT * (-j_srm) + c_srm
    new_c_slm = DT * (-j_slm) + c_slm

    # 构建系数矩阵与右端向量
    A = np.zeros(nano_point_count)  # 系数矩阵
    b = np.zeros(nano_point_count)  # 右端向量
    for ctrl_point_id in range(nano_point_count):
        s = ctrl_triangle_area[ctrl_point_id]
        # 所求点控制的每个三角形单元的相关计算
        for ctrl_triangle_index in range(0, max_ctrl_triangle_count):
            triangle_id = ctrl_triangle_list[ctrl_point_id, ctrl_triangle_index, 0]
            if triangle_id != -1:
                # 在此三角形单元中的顶点编号
                ctrl_vertex_id = ctrl_triangle_list[ctrl_point_id, ctrl_triangle_index, 1]
                non1_vertex_id = (ctrl_vertex_id + 1) % 3
                non2_vertex_id = (ctrl_vertex_id + 2) % 3
                # 顶点对应的点编号
                non1_point_id = triangle_matrix[triangle_id, non1_vertex_id]  # 另外两点之一
                non2_point_id = triangle_matrix[triangle_id, non2_vertex_id]  # 另外两点之二
                # 扩散项
                # 判断所求点在此单元的对边是否为边界
                if point_och[non1_point_id] == B_INNER or point_och[non2_point_id] == B_INNER:
                    non1_point_c_ca = c_ca[non1_point_id]
                    non2_point_c_ca = c_ca[non2_point_id]
                    b1 = abc_matrix[triangle_id, ctrl_vertex_id, 1]
                    b2 = abc_matrix[triangle_id, non1_vertex_id, 1]
                    b3 = abc_matrix[triangle_id, non2_vertex_id, 1]
                    c1 = abc_matrix[triangle_id, ctrl_vertex_id, 2]
                    c2 = abc_matrix[triangle_id, non1_vertex_id, 2]
                    c3 = abc_matrix[triangle_id, non2_vertex_id, 2]
                    l = nl_matrix[triangle_id, ctrl_vertex_id, 0]
                    nx = nl_matrix[triangle_id, ctrl_vertex_id, 1]
                    ny = nl_matrix[triangle_id, ctrl_vertex_id, 2]
                    A[ctrl_point_id] += D_CA * (b1 * nx + c1 * ny) * l * DT
                    b[ctrl_point_id] += D_CA * ((non1_point_c_ca * b2 + non2_point_c_ca * b3) * nx + (
                            non1_point_c_ca * c2 + non2_point_c_ca * c3) * ny) * l * DT

                # 获取此单元的边界信息
                boundary_type_status = boundary_type[triangle_id, 0]
                boundary_length = boundary_type[triangle_id, 1]  # 边界长度
                # 入流项
                # 判断此单元是否存在入流边界
                if boundary_type_status == B_INFLOW and k_ryr != 0:
                    fj = (c_ca[non1_point_id] + c_ca[non2_point_id]) / 3.0
                    A[ctrl_point_id] += -k_ryr * boundary_length / 3 * DT
                    b[ctrl_point_id] += k_ryr * (c_ca_jsr - fj) * boundary_length * DT
                # 出流项
                # 判断此单元是否存在出流边界
                elif boundary_type_status == B_OUTFLOW:
                    fk = (c_ca[non1_point_id] + c_ca[non2_point_id]) / 3.0
                    A[ctrl_point_id] += -K_OUT * boundary_length / 3 * DT
                    b[ctrl_point_id] += K_OUT * (C_CA_OUT - fk) * boundary_length * DT

            else:
                break
        # 最终此点A矩阵的对应行向量与b向量的对应的标量
        A[ctrl_point_id] = s - A[ctrl_point_id]
        b[ctrl_point_id] += s * (j_item[ctrl_point_id] * DT + c_ca[ctrl_point_id])

    # 解方程
    new_c_ca = b / A
    return new_c_ca, new_c_cag, new_c_cam, new_c_trc, new_c_srm, new_c_slm


def cal_caf_concentration(c_ca, c_caf):
    j_dye_f = cal_j_dye_f(c_ca, c_caf)

    # 构建系数矩阵与右端向量
    A = np.zeros(nano_point_count)  # 系数矩阵
    b = np.zeros(nano_point_count)  # 右端向量
    for ctrl_point_id in range(nano_point_count):
        s = ctrl_triangle_area[ctrl_point_id]
        # 所求点控制的每个三角形单元的相关计算
        for ctrl_triangle_index in range(0, max_ctrl_triangle_count):
            triangle_id = ctrl_triangle_list[ctrl_point_id, ctrl_triangle_index, 0]
            if triangle_id != -1:
                ctrl_vertex_id = ctrl_triangle_list[ctrl_point_id, ctrl_triangle_index, 1]
                non1_vertex_id = (ctrl_vertex_id + 1) % 3
                non2_vertex_id = (ctrl_vertex_id + 2) % 3

                non1_point_id = triangle_matrix[triangle_id, non1_vertex_id]  # 另外两点之一
                non2_point_id = triangle_matrix[triangle_id, non2_vertex_id]  # 另外两点之二

                # 扩散项
                # 判断所求点在此单元的对边是否为边界
                if point_och[non1_point_id] == B_INNER or point_och[non2_point_id] == B_INNER:
                    non1_point_c_caf = c_caf[non1_point_id]
                    non2_point_c_caf = c_caf[non2_point_id]
                    b1 = abc_matrix[triangle_id, ctrl_vertex_id, 1]
                    b2 = abc_matrix[triangle_id, non1_vertex_id, 1]
                    b3 = abc_matrix[triangle_id, non2_vertex_id, 1]
                    c1 = abc_matrix[triangle_id, ctrl_vertex_id, 2]
                    c2 = abc_matrix[triangle_id, non1_vertex_id, 2]
                    c3 = abc_matrix[triangle_id, non2_vertex_id, 2]
                    l = nl_matrix[triangle_id, ctrl_vertex_id, 0]
                    nx = nl_matrix[triangle_id, ctrl_vertex_id, 1]
                    ny = nl_matrix[triangle_id, ctrl_vertex_id, 2]
                    A[ctrl_point_id] += D_CAF * (b1 * nx + c1 * ny) * l * DT
                    b[ctrl_point_id] += D_CAF * ((non1_point_c_caf * b2 + non2_point_c_caf * b3) * nx + (
                            non1_point_c_caf * c2 + non2_point_c_caf * c3) * ny) * l * DT

                # 获取此单元的边界信息
                boundary_type_status = boundary_type[triangle_id, 0]
                boundary_length = boundary_type[triangle_id, 1]
                # 出流项
                # 判断此单元是否存在出流边界
                if boundary_type_status == B_OUTFLOW:
                    fk = (c_caf[non1_point_id] + c_caf[non2_point_id]) / 3.0
                    # A[ctrl_point_id] += -K_FOUT * boundary_length / 3 * DT
                    # b[ctrl_point_id] += K_FOUT * (C_CAF_OUT - fk) * boundary_length * DT
            else:
                break
        A[ctrl_point_id] = s - A[ctrl_point_id]
        b[ctrl_point_id] += s * (-j_dye_f[ctrl_point_id] * DT + c_caf[ctrl_point_id])
    # 解方程
    new_c_caf = b / A
    return new_c_caf


def cal_ryr_cca(c_ca):
    from math import sqrt
    c_ca_ryr = 0.0
    arc_len = 0.0
    for i in range(0, triangle_count):
        for j in range(0, 3):
            p1 = triangle_matrix[i, (j + 1) % 3]
            p2 = triangle_matrix[i, (j + 2) % 3]
            if (point_och[p1] == B_OUTFLOW) and (point_och[p2] == B_OUTFLOW):
                length = sqrt((point_x[p1] - point_x[p2]) ** 2 + (point_y[p1] - point_y[p2]) ** 2)
                c_ca_ryr = c_ca_ryr + length * (c_ca[p1] + c_ca[p2]) / 2
                arc_len = arc_len + length
    ryr_cca = c_ca_ryr / arc_len
    return ryr_cca


def nano_core(c_ca, c_cag, c_cam, c_trc, c_srm, c_slm, c_ca_jsr, k_ryr):
    new_c_ca, new_c_cag, new_c_cam, new_c_trc, new_c_srm, new_c_slm = \
        cal_ca_concentration(c_ca, c_cag, c_cam, c_trc, c_srm, c_slm, c_ca_jsr, k_ryr)
    # new_c_caf = cal_caf_concentration(c_ca, c_caf)
    ryr_cca = cal_ryr_cca(new_c_ca)
    # n 步的浓度
    # out_boundary_c_ca = np.sum(c_ca[84:205] * out_boundary_length) / np.sum(out_boundary_length)
    # out_boundary_c_caf = np.sum(c_caf[84:205] * out_boundary_length) / np.sum(out_boundary_length)

    return new_c_ca, new_c_cag, new_c_cam, new_c_trc, new_c_srm, new_c_slm, ryr_cca
