from pathlib import Path
import os
import linecache
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def specify_the_line(path, line):
    path_object = Path(path)
    if not path_object.exists():
        print("此路径不存在！")
        return

    result = []
    for filename in os.listdir(path):
        data = float(linecache.getline(f"{path}/{filename}", line).strip('\n'))
        result.append(data)
    return np.array(result)



#  画折线图
def plotter(path, title):
    path_object = Path(path)
    if path_object.exists():
        parent = path_object.parent
    else:
        print("此路径不存在！")
        return

    y = np.loadtxt(path)
    x = np.arange(len(y))

    # 画折线图
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, y, ls='-', lw=2)

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.suptitle(title, fontsize=24)  # 设置折线图更上一级标题
    plt.title(f"数据来源：{parent}")  # 设置折线图标题

    plt.grid()  # 添加网格
    plt.savefig(f"{parent}/{title}.png")
    plt.show()


if __name__ == '__main__':
    my_path = "../Result/TS50000_kRyR=311999711.0000832_08-08-14-06-05/Ca"
    my_line = 24
    my_title = "ryr_var"
    plotter(specify_the_line(my_path, my_line), my_title)
