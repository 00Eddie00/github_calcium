import numpy as np
from decimal import Decimal as De

parent_path = "../../blinkgrid/"
grid = np.loadtxt(parent_path + "gridt.dat", dtype=np.float64)
triangle_matrix = np.loadtxt(parent_path + "nod.dat", dtype=np.int64) - 1


def cal_cross_product(x, y, x1, y1, x2, y2):
    return (De(str(x1)) - De(str(x))) * (De(str(y2)) - De(str(y))) - (De(str(x2)) - De(str(x))) * (De(str(y1)) - De(str(y)))


# 判断某点与某三角之间的关系
def judge_relation(x, y, triangle_id):
    # 第一个值表示关系类型：1在顶点，2在边上，3在内部
    # 第二个值：如果此点是某三角形的顶点，则此值为顶点的点编号；如果不是，则与其有关系的三角形编号
    exclusive_relation = np.full(2, -1)

    vertex1 = triangle_matrix[triangle_id, 0]
    vertex2 = triangle_matrix[triangle_id, 1]
    vertex3 = triangle_matrix[triangle_id, 2]
    x1 = grid[vertex1, 0]
    y1 = grid[vertex1, 1]
    x2 = grid[vertex2, 0]
    y2 = grid[vertex2, 1]
    x3 = grid[vertex3, 0]
    y3 = grid[vertex3, 1]
    # 判断是否为三角形单元的某顶点
    if x == x1 and y == y1:
        exclusive_relation[0] = 1
        exclusive_relation[1] = vertex1
        return exclusive_relation
    if x == x2 and y == y2:
        exclusive_relation[0] = 1
        exclusive_relation[1] = vertex2
        return exclusive_relation
    if x == x3 and y == y3:
        exclusive_relation[0] = 1
        exclusive_relation[1] = vertex3
        return exclusive_relation
    # 计算叉积
    c1 = cal_cross_product(x, y, x1, y1, x2, y2)
    c2 = cal_cross_product(x, y, x2, y2, x3, y3)
    c3 = cal_cross_product(x, y, x3, y3, x1, y1)
    # 判断是否在内部
    if (c1 > 0 and c2 > 0 and c3 > 0) or (c1 < 0 and c2 < 0 and c3 < 0):
        exclusive_relation[0] = 3
        exclusive_relation[1] = triangle_id
        return exclusive_relation
    # 判断是否在边界上
    if (c1 == 0 and c2 * c3 > 0) or (c2 == 0 and c3 * c1 > 0) or (c3 == 0 and c1 * c2 > 0):
        exclusive_relation[0] = 2
        exclusive_relation[1] = triangle_id
        return exclusive_relation
    # 在外部
    return exclusive_relation

# 得到所有点所在的三角形
def get_relations():
    print("开始")
    relations = np.full((601, 601, 2), dtype=int,fill_value=-1)
    count = 1
    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            if x ** 2 + y ** 2 <= 90000:
                for triangle_id in range(len(triangle_matrix)):
                    exclusive_relation = judge_relation(x, y, triangle_id)
                    if exclusive_relation[0] != -1:
                        relations[i, j] = exclusive_relation
                        break
            print(count)
            count += 1

    return relations

# 对某些点进行特殊处理
def relations_refine():
    relations = np.load("./refined_relations.npy")
    refined_relations = np.copy(relations)

    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            r_square = x ** 2 + y ** 2
            if  r_square <= 90000 and relations[i, j, 0] == -1:
                r = np.sqrt(r_square)
                if r > 299:
                    refined_relations[i, j, 0] = 1
                    refined_relations[i, j, 1] = 84
                elif r < 299:
                    refined_relations[i, j, 0] = 1
                    if x == 15 and y == 15:
                        refined_relations[i, j, 1] = 0
                    elif x == 15 and y == -15:
                        refined_relations[i, j, 1] = 21
                    elif x == -15 and y == -15:
                        refined_relations[i, j, 1] = 42
                    elif x == -15 and y == 15:
                        refined_relations[i, j, 1] = 63

    np.save("./refined_relations", refined_relations)


if __name__ == "__main__":
    # np.save("./relations", get_relations())

    relations = np.load("./relations.npy")
    refined_relations = np.copy(relations)
    count = 0
    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            if x ** 2 + y ** 2 <= 90000 and relations[i, j, 0] == -1:
                r = np.sqrt(x ** 2 + y ** 2)
                print(f"i:{i},j:{j} --- ({x},{y}) --- {r} --- {r-300}")
                count += 1
    print(count)
