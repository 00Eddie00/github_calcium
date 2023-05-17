import os
import linecache
import csv
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


# 读取平均浓度到指定文件内
def avg_ca(dir_name, data_name):
    #line = 1 # 中心处
    #line = 500 # 非中心非边界
    #line = 1001 #边界处
    line = 1002  # 指定行
    if data_name == "Ca":
        avg_file = "Ca.csv"
        sub_dir_name = "Ca"
    elif data_name == "CaF":
        avg_file = "CaF.csv"
        sub_dir_name = "CaF"
    else:
        return False
    with open(dir_name + avg_file, "w") as file_object:
        for filename in os.listdir(dir_name + sub_dir_name):
            data = linecache.getline(dir_name + sub_dir_name + "\\" + filename, line)
            file_object.write(data)
    return True


#  画折线图
def plotter(dir_name, data_name, data_source):
    x = []
    y = []

    if data_name == "Ca":
        avg_file = "Ca.csv"
        y_locator = 0.0001
    elif data_name == "CaF":
        avg_file = "CaF.csv"
        y_locator = 0.0001
    else:
        return False

    # 将文件内容读入数组
    i = 0
    with open(dir_name + avg_file, 'r') as ca_csv_file:
        plots = csv.reader(ca_csv_file, delimiter=',')
        for row in plots:
            x.append(i * 100)
            y.append(float(row[0]))
            i = i + 1

    # 画折线图
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, y, ls='-', lw=2)

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.suptitle(data_name+'均值', fontsize=24)  # 设置折线图更上一级标题
    plt.title("数据来源：" + data_source)  # 设置折线图标题
    # 设置刻度间隔
    #x_major_locator = MultipleLocator(2500)  # 把x轴的刻度间隔设置为1000，并存在变量里
    #y_major_locator = MultipleLocator(y_locator)  # 把y轴的刻度间隔设置为0.02，并存在变量里
    #ax = plt.gca()  # ax为两条坐标轴的实例
    #ax.xaxis.set_major_locator(x_major_locator)  # 把x轴的主刻度设置为1000的倍数
    #ax.yaxis.set_major_locator(y_major_locator)  # 把y轴的主刻度设置为0.001的倍数
    # 设置刻度范围
    #plt.xlim(0, 50000)  # 把x轴的刻度范围设置为-500到15000，因为500不满一个刻度间隔，所以数字不会显示出来，但是能看到一点空白
    #plt.ylim(0, 0.00016)  # 把y轴的刻度范围设置为-0.0005到1，同理，-0.0005不会标出来，但是能看到一点空白

    plt.grid()  # 添加网格
    plt.savefig(dir_name + data_name + '均值' + '.png')
    plt.show()


def fn_gn(data_source):
    dir_name = "..\\Result\\" + data_source + "\\"
    avg_ca(dir_name, "Ca")
    avg_ca(dir_name, "CaF")
    plotter(dir_name, "Ca", data_source)
    plotter(dir_name, "CaF", data_source)


if __name__ == '__main__':
    my_data_source = "result_TotalSteps50000_07-30-16-40zhang"

    fn_gn(my_data_source)
