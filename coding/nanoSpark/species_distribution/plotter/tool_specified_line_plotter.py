from pathlib import Path
import os
import linecache
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def specify_the_line(path, line, store_file_name):
    path_object = Path(path)
    if path_object.exists():
        parent = path_object.parent
    else:
        print("此路径不存在！")
        return
    store_file_path = f"{parent}/{store_file_name}"
    with open(store_file_path, "w") as store_file:
        for filename in os.listdir(path):
            data = linecache.getline(f"{path}{filename}", line)
            store_file.write(data)

    print("success")
    return store_file_path


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


def specified_line_plotter(path, line, title):
    plotter(specify_the_line(path, line, f"{title}.csv"), title)


if __name__ == '__main__':
    my_path = "../../Result/TS70000_kRyR=311999711.0000832_09-23-11-38-22-new_grid/CaF/"
    my_line = 24
    my_title = "ryr_caf_var"
    my_store = specify_the_line(my_path, my_line, f"{my_title}.csv")
    plotter(my_store, my_title)
