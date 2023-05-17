import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def plotter(dir_name, species):
    y = np.loadtxt(f"{dir_name}avg_c_{species}.csv")
    x = np.arange(len(y))

    # 画折线图
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, y, ls='-', lw=2)

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.title(f"{species}均值")  # 设置折线图标题

    plt.grid()  # 添加网格
    plt.savefig(f"{dir_name}avg_c_{species}.png")
    plt.show()


if __name__ == '__main__':
    my_dir_name = "../Result/TS50000_kRyR=311999711.0000832_08-13-14-41-34/"
    plotter(my_dir_name, "Ca")
    plotter(my_dir_name, "CaF")
