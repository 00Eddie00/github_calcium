import numpy as npp
from math import sqrt


def not_empty(s):
    return s and s.strip()


def main():
    '''
    输入：sidegridt2.dat -->  大圆，小圆上所有点的坐标，第1行是点数，小圆上70个点，大圆776个点
        sidenpoch2.dat -->  大圆，小圆上所有点的属性，多1行，多出来的是大圆上第一个入流点和最后一个入流点属性重复了
    :return:
    '''
    boun = npp.zeros(100, int)
    nboun = 0
    npoch = npp.zeros(10000, int)
    nboun = 2
    boun[0] = 70    # 小圆的点数
    boun[1] = 846   # 小圆 + 大圆
    sideface = open("sideface.dat", "w")
    sidegridt = open("sidegridt.dat", "r")
    sidenpoch = open("sidenpoch.dat", "r")
    sidegridts = npp.zeros((1000, 2), float)
    j = 0
    # for line in sidegridt.readlines():
    #     current_line = list(filter(not_empty, line.strip("\n").split(" ")))
    #     sidegridts[j][0] = current_line[0]
    #     sidegridts[j][1] = current_line[1]
    #     j = j + 1
    # nb = j
    #sideface.write(str(nb))
    d = 0
    sidenpochs = npp.empty(1000, float)
    for line in sidenpoch.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        sidenpochs[d] = current_line[0]
        d = d + 1
    nb=d  # 边界点个数
    sideface.write(str(nb)+ "\n")   # sideface 第一行写：统计npoch有多少行
    for i in range(0, nb):
        npoch[i] = sidenpochs[i]
    istart = 0
    for i in range(0, nboun):  # nboun=2
        iend = boun[i]  # boun[0]=70 boun[1]=846
        for j in range(istart, iend - 1):
            k = npoch[j]
            if k != npoch[j + 1]:    # 线段两端一个是入流点，一个是壁面点，那么这个线段上的点都是壁面点
                k = 1
            sideface.write(str(j+1) + " " + str(j + 2) + " " + str(k)+"\n")
        k = npoch[iend-1]
        if k != npoch[istart]:        # 最后一个点和起点的比较
            k = 1
        sideface.write(str(iend) + " " + str(istart+1) + " " + str(k)+"\n")
        istart = iend
    sideface.close()
    sidegridt.close()
    sidenpoch.close()
main()