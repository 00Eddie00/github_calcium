from math import cos, sin, asin

import numpy as npp


def xyarcbo():
    '''
    可以理解为是一个大圆，圆上分为不同的圆弧，每段圆弧都有不同的属性，或是壁面或是入流
    :return:
    xbo: 每段小圆弧的起始点的x坐标
    ybo: 每段小圆弧的起始点的y坐标
    arcbo: 每段圆弧的弧长
    npoch: 每段圆弧的属性
    '''
    xbo = npp.zeros(100, float)
    ybo = npp.zeros(100, float)
    arcbo = npp.zeros(100, float)
    i = 0
    npoch = npp.zeros(100, int)
    pai = 0.00
    rad = 0.00
    arc = 0.00
    pai = 2 * asin(1.0)
    rad = 225  # 半径
    arc = 0

    file = open("arc.dat","w")
    # 半径缩小，fSR管子尺寸不变，需修改成4个
    arcbo[0] = 15
    arcbo[1] = rad * (pai / 2) - 30
    arcbo[2] = 30
    arcbo[3] = rad * (pai / 2) - 30
    arcbo[4] = 30
    arcbo[5] = rad * (pai / 2) - 30
    arcbo[6] = 30
    arcbo[7] = rad * (pai / 2) - 30
    arcbo[8] = 15
    npoch[0] = 2
    npoch[1] = 1
    npoch[2] = 2
    npoch[3] = 1
    npoch[4] = 2
    npoch[5] = 1
    npoch[6] = 2
    npoch[7] = 1
    npoch[8] = 2
    for i in range(0, 10):  # 最后一个点和第一个点重合(300,0)，18个点, 17段弧长
        xbo[i] = rad * cos(arc / rad)  # arc/rad = 圆心角度数, xbo,ybo是每段弧长起始点坐标
        ybo[i] = rad * sin(arc / rad)
        arc = arc + arcbo[i]
        file.write(str(i)+" "+str(xbo[i])+" "+ str(ybo[i])+ "\n")
    file.close()

xyarcbo()