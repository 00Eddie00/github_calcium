import matplotlib.pyplot as plt
import numpy as np

def fluo_contourf():
    # 准备数据
    x = np.linspace(-300, 300, 601)  # r轴坐标
    y = np.linspace(-300, 300, 601)  # r轴坐标
    # 处理浓度数据
    original_data = np.load("./refined_relations.npy")[:,:,0]
    count1 = 0
    count2 = 0
    count3 = 0
    count0 = 0

    for i in range(601):
        for j in range(601):
            if original_data[i ,j] == 1:
                print(i,j)
                count1 += 1
            elif original_data[i ,j] == 2:
                count2 += 1
            elif original_data[i ,j] == 3:
                count3 += 1
            elif original_data[i ,j] == -1:
                count0 += 1
    print(count1)
    print(count2)
    print(count3)
    print(count0)

    # # 绘制
    # plt.contourf(x, y, original_data, levels=10)
    # plt.xlabel('x')  # 设置x轴标签
    # plt.ylabel('y')  # 设置y轴标签
    # plt.colorbar()
    # # 保存并显示图表
    # plt.show()

fluo_contourf()

