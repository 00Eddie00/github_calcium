import numpy as np
import matplotlib.pyplot as plt
from species_distribution.plotter.tool_specified_line_plotter import specify_the_line

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def plotter(caf):
    # y = (caf - caf[0]) / caf[0]
    cag = np.loadtxt("E:\Code\Python\PycharmProjects/nano_normal_spk\Result\TS50000_kRyR=311999711.0000832_10-13-14-56-34/new_cag_t_dis.csv")
    y = caf / caf[0]
    y1 = cag / cag[0]
    x = [i * 2 * 10 ** -6 * 100 * 1000 for i in range(len(y))]
    max_cag = max(y1)
    max_index = y1.index(max_cag)
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, y, ls='-', lw=2, label="only GcaMP6f")
    plt.plot(x, y1, ls='-', lw=2, label="GcaMP6f & fluo-3")
    plt.text(max_index, max_cag, (max_index, max_cag), color='r')
    plt.xlabel('t(ms)', fontsize=14)  # 设置x轴标签
    # plt.ylabel('concentration(mM) at z=7.5', fontsize=14)  # 设置y轴标签
    plt.legend()
    x_major_locator = plt.MultipleLocator(10)
    # 把x轴的刻度间隔设置为1，并存在变量里
    y_major_locator = plt.MultipleLocator(0.25)
    # 把y轴的刻度间隔设置为10，并存在变量里
    ax = plt.gca()
    # ax为两条坐标轴的实例
    ax.xaxis.set_major_locator(x_major_locator)
    # 把x轴的主刻度设置为1的倍数
    ax.yaxis.set_major_locator(y_major_locator)
    plt.title("F/F0")  # 设置折线图标题

    plt.grid()  # 添加网格
    plt.savefig(f"new_f_div_f0.png")
    plt.show()


if __name__ == "__main__":
    path = "../Result/TS70000_kRyR=311999711.0000832_09-23-11-38-22-new_grid/t_dis.csv"
    is_dir = False
    if is_dir:
        line = 24
        store_file_name = ""
        specify_the_line(path, line, store_file_name)
    else:
        plotter(np.loadtxt(path))
