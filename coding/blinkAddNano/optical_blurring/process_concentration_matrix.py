import numpy as np

parent_path = "../nano_grid/"
pre_parent_path = f"{parent_path}/pre_info/"
relations = np.load(
    "E:\\Code\\Python\\PycharmProjects\\unstructured_nano_spark\\optical_blurring\\relation\\refined_relations.npy")
abc_matrix = np.load(pre_parent_path + "abc_matrix.npy")
triangle_matrix = np.loadtxt(parent_path + "4RYRnod.dat", dtype=np.int64) - 1


def process_concentration(original_concentration, c_val):
    concentration = np.full((601, 601, 16), c_val, dtype=float)
    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            relation = relations[i, j, 0]
            c_id = relations[i, j, 1]
            if relation == 1:
                concentration[i, j, :] = original_concentration[c_id]
            elif relation == 2 or relation == 3:
                vertex1 = triangle_matrix[c_id, 0]
                vertex2 = triangle_matrix[c_id, 1]
                vertex3 = triangle_matrix[c_id, 2]

                f1 = original_concentration[vertex1]
                f2 = original_concentration[vertex2]
                f3 = original_concentration[vertex3]

                a1 = abc_matrix[c_id, 0, 0]
                b1 = abc_matrix[c_id, 0, 1]
                c1 = abc_matrix[c_id, 0, 2]

                a2 = abc_matrix[c_id, 1, 0]
                b2 = abc_matrix[c_id, 1, 1]
                c2 = abc_matrix[c_id, 1, 2]

                a3 = abc_matrix[c_id, 2, 0]
                b3 = abc_matrix[c_id, 2, 1]
                c3 = abc_matrix[c_id, 2, 2]
                concentration[i, j, :] = \
                    f1 * (a1 + b1 * x + c1 * y) + f2 * (a2 + b2 * x + c2 * y) + f3 * (a3 + b3 * x + c3 * y)

    return concentration


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    data = np.loadtxt('E:\\Study\\Unstructured ResultSets\\ResultSet_22-0406-2027_\\NANO\\CaG\\CaG00010000.csv')
    processed_concentration = process_concentration(data, 0)
    # 准备数据
    ryr = processed_concentration[:, :, 15]
    x = np.linspace(-300, 300, 601)  # r轴坐标
    y = np.linspace(-300, 300, 601)  # r轴坐标
    # 绘制
    plt.contourf(x, y, ryr, levels=1000)
    plt.xlabel('x')  # 设置x轴标签
    plt.ylabel('y')  # 设置y轴标签
    plt.colorbar()
    # 保存并显示图表
    plt.show()
