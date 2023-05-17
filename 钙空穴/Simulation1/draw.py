import os
import linecache
import csv
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


# 读取平均浓度到指定文N 件内
def avg_ca(dir_name, data_name):  # DATA\\DATA_blink\\  fn
    # line = 1 # 中心处
    # line = 500 # 非中心非边界
    # line = 1001 #边界处
    line = 6414  # 指定行
    # line1 = 2
    if data_name == "fn":
        avg_file = "fn.csv"
        sub_dir_name = "Fn"
    elif data_name == "gn":
        avg_file = "gn.csv"
        sub_dir_name = "Gn"
    else:
        return False

    i = 0
    with open(dir_name + avg_file, "w") as file_object:
        for filename in os.listdir(dir_name + sub_dir_name):
            data = linecache.getline(dir_name + sub_dir_name + "\\" + filename, line)
            file_object.write(data)
            i = i + 1
            print(i)
    return True


#  画折线图
def plotter(dir_name, data_name, data_source):
    x1 = []
    y1 = []

    if data_name == "fn":
        avg_file1 = "fn.csv"
        y_locator = 0.0001
    elif data_name == "gn":
        avg_file1 = "gn.csv"
        y_locator = 0.0001
    else:
        return False

    # 将文件内容读入数组
    i = 0
    with open(dir_name + avg_file1, 'r') as ca_csv_file:
        plots = csv.reader(ca_csv_file, delimiter=',')
        for row in plots:
            x1.append(i)
            y1.append(float(row[0]))
            i = i + 1

    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x1, y1, ls='-', lw=2)

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.suptitle(data_name + '均值', fontsize=24)  # 设置折线图更上一级标题
    plt.title("数据来源：" + data_source)  # 设置折线图标题
    # 设置刻度间隔
    # x_major_locator = MultipleLocator(2500)  # 把x轴的刻度间隔设置为1000，并存在变量里
    # y_major_locator = MultipleLocator(y_locator)  # 把y轴的刻度间隔设置为0.02，并存在变量里
    # ax = plt.gca()  # ax为两条坐标轴的实例
    # ax.xaxis.set_major_locator(x_major_locator)  # 把x轴的主刻度设置为1000的倍数
    # ax.yaxis.set_major_locator(y_major_locator)  # 把y轴的主刻度设置为0.001的倍数
    # 设置刻度范围
    # plt.xlim(0, 50000)  # 把x轴的刻度范围设置为-500到15000，因为500不满一个刻度间隔，所以数字不会显示出来，但是能看到一点空白
    # plt.ylim(0, 0.00016)  # 把y轴的刻度范围设置为-0.0005到1，同理，-0.0005不会标出来，但是能看到一点空白

    # plt.legend(loc=0)

    plt.grid()  # 添加网格
    plt.savefig(dir_name + data_name + '均值' + '.png')
    plt.show()


def fn_gn(data_source):
    dir_name = "DATA\\" + data_source + "\\"
    avg_ca(dir_name, "fn")
    print("avg_cafn 执行结束")
    avg_ca(dir_name, "gn")
    print("avg_cagn 执行结束")
    plotter(dir_name, "fn", data_source)
    print("plotterfn 执行结束")
    plotter(dir_name, "gn", data_source)
    print("plottergn 执行结束")


if __name__ == '__main__':
    my_data_source = "DATA_blink"
    fn_gn(my_data_source)
