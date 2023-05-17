import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def multi_plot_temporal_distribution(area, species, suffix, save=True):
    target_file = f"{area}_{species}_{suffix}.csv"
    # y = np.loadtxt(f"E:\\Study\\sr&sl参数调整\\情况九子情况1\\{target_file}")[:501]
    # y1 = np.loadtxt(f"E:\\Study\\sr&sl参数调整\\情况九子情况2\\{target_file}")[:501]
    y2 = np.loadtxt(f"E:\\Study\\sr&sl参数调整\\情况九子情况3\\{target_file}")[:501]
    y4 = np.loadtxt(f"E:\\Study\\sr&sl参数调整\\情况八子情况1\\{target_file}")[:501]

    x = np.arange(len(y2))

    # 画折线图
    plt.figure(dpi=300, figsize=(10, 10))
    # plt.plot(x, y, ls='-', lw=2, label='情况九子情况1', color='steelblue')
    # plt.plot(x, y1, ls='-', lw=2, label='情况九子情况2', color='orange')
    plt.plot(x, y2 / y2[0] - 1, ls='-', lw=2, label='情况九子情况3', color='green')
    plt.plot(x, y4 / y4[0] - 1, ls='-', lw=2, label='情况八子情况1', color='firebrick')

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.legend()
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.suptitle(f"{area}_{species}_{suffix}")  # 设置折线图标题
    # data_source = result_set[result_set.find("ResultSet_"):] + f"/{area}_{species}_{suffix}.csv"
    # plt.title(f"数据来源：{data_source}")

    plt.grid()  # 添加网格
    if save:
        plt.savefig(f"E:\\Study\\sr&sl参数调整\\{area}_{species}_{suffix}_情况九子情况3与情况八对比.png")
    plt.show()


def plot_temporal_distribution(area, species, suffix, result_set, save=True, ):
    target_file = f"{area}_{species}_{suffix}.csv"
    y = np.loadtxt(f"{result_set}\\{target_file}")
    x = np.arange(len(y))

    # 画折线图
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, y, ls='-', lw=2, label=f'{result_set}')

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.legend()
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.suptitle(f"{area}_{species}_{suffix}")  # 设置折线图标题
    data_source = result_set[result_set.find("ResultSet_"):] + f"/{area}_{species}_{suffix}.csv"
    plt.title(f"数据来源：{data_source}")

    plt.grid()  # 添加网格
    if save:
        plt.savefig(f"{result_set}\\{area}_{species}_{suffix}.png")
    plt.show()


if __name__ == '__main__':
    my_result_set = "E:\\Study\\sr&sl参数调整\\情况八子情况1"
    # plot_temporal_distribution(my_result_set, "CaG", "avg_c")
    # plot_temporal_distribution(my_result_set, "CaF", "avg_c")
    multi_plot_temporal_distribution("NANO", "CaG", "avg_c")
