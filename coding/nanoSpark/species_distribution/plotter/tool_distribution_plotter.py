import numpy as np
import matplotlib.pyplot as plt
from math import log10

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def plotter(dir_name, dis_type):
    y = np.loadtxt(f"{dir_name}.csv")
    if dis_type == '空间分布':
        x = [i - 300 for i in range(len(y))]
        x_label = 'x(nm)'
    elif dis_type == '时间分布':
        x = [i * 2 * 10 ** -6 * 100 * 1000 for i in range(len(y))]
        x_label = 't(ms)'
    else:
        return False

    ca_bc = [log10(y[:, 0][i])for i in range(len(y))]
    caf_psf = [log10(y[:, 1][i])for i in range(len(y))]
    ca_bc_psf = [log10(y[:, 2][i])for i in range(len(y))]

    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, ca_bc, ls='-', lw=2, label='[Ca]BC')
    # plt.plot(x, caf_psf, ls='-', lw=2, label='[CaF]w/PSF')
    # plt.plot(x, ca_bc_psf, ls='-', lw=2, label='[Ca]BCw/PSF')
    plt.legend()

    plt.xlabel(x_label, fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration(mM) at z=7.5', fontsize=14)  # 设置y轴标签
    plt.title(dis_type)  # 设置折线图标题

    plt.grid()  # 添加网格
    plt.savefig(f"{dir_name}_log10_ca_bc.png")
    plt.show()


# plotter("../species_distribution/x_dis", '空间分布')
plotter("../species_distribution/t_dis", '时间分布')
