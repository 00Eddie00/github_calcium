from math import asin, sin, cos

import numpy as npp

global boungridts, bounpochs
boungridts = npp.zeros((1000, 2), float)
bounpochs = npp.zeros(1000, float)

pmax = 100000
emax = 200000
bmax = 50000
fmin = 10 * -8
fm2 = 10 * -8

global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
np = 0
nboun = 0
sidetype = 0
nbs = 0
bounarc = 0
px = npp.zeros(bmax, float)
py = npp.zeros(bmax, float)
pstep = npp.zeros(bmax, float)
bx = npp.zeros(bmax, float)
by = npp.zeros(bmax, float)
bs = npp.zeros(bmax, float)

global gx, gy, gcx, gcy, gcr, x, y, xce, yce
gx = npp.zeros(3, float)
gy = npp.zeros(3, float)
gcx = 0.00
gcy = 0.00
gcr = 0.00
x = 0.00
y = 0.00
xce = 0.00
yce = 0.00

global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
bgnp = 0
bgne = 0
bgpx = npp.zeros(pmax, float)
bgpy = npp.zeros(pmax, float)
bgstep = npp.zeros(pmax, float)
bgnod = npp.zeros((3, emax), int)
bgnoe = npp.zeros((3, emax), int)
bgcx = npp.zeros(emax, float)
bgcy = npp.zeros(emax, float)
bgcr = npp.zeros(emax, float)


def not_empty(s):
    return s and s.strip()


def main():
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    iloop = 0
    npoch = npp.zeros(2000, int)
    arc = 0.00
    x0 = 0.00
    y0 = 0.00
    x1 = 0.00
    y1 = 0.00
    pai = 0.00
    arcboun = 0.00
    xbo = npp.zeros(2000, float)
    ybo = npp.zeros(2000, float)
    arcbo = npp.zeros(2000, float)
    filename = ""

    sidetype = 2
    boungridt_file = open("boungridt2.dat", "w+")
    bounpoch_file = open("bounpoch2.dat", "w+")
    # sidetype=1 离散形式的边界
    if sidetype == 1:
        filename = 'boun3.dat'
        loadboun(filename)
        x0 = bx[0]
        y0 = by[0]
        x1 = bx[nboun]
        y1 = by[nboun]
        arc = bs[nboun]
        arcscatter(x0, y0, x1, y1, arc)
        savedata(1, boungridt_file, bounpoch_file)
    # sidetype=2 单段连续边界，函数形式，须改动
    elif sidetype == 2:
        pai = 2 * asin(1.0)  # Π值
        x0 = 1.0  # 该弧起点的坐标值
        y0 = 0
        x1 = 1.0  # 该弧终点的坐标值
        y1 = 0
        arc = 2 * pai
        arcscatter(x0, y0, x1, y1, arc)
        savedata(4, boungridt_file, bounpoch_file)
        # sidetype=3 分段连续边界，函数形式，须改动
    elif sidetype == 3:
        xbo, ybo, arcbo, npoch = xyarcbo()
        arcboun = 0
        for iloop in range(0, 17):
            x0 = xbo[iloop]
            y0 = ybo[iloop]
            x1 = xbo[iloop + 1]
            y1 = ybo[iloop + 1]
            arc = arcbo[iloop]
            arcscatter(x0, y0, x1, y1, arc)
            np = np - 1
            savendata(npoch[iloop], boungridt_file, bounpoch_file)
            arcboun = arcboun + arcbo[iloop]
    boungridt_file.close()
    bounpoch_file.close()


# ****************************************************
# 单段连续边界，弧长arc和直角坐标x，y的关系，需改动
def arc_to_xy(arc):  # arc = (arc + arc0 - arcxy) / 2 , 插入的平均步长
    pai = 0.00
    pai = 2 * asin(1.0)
    x = cos(-arc)  # cos(-Π) = -1
    y = sin(-arc)  # sin(-Π) = 0
    return x, y


# 分段连续边界，各段起止点、弧长
def xyarcbo():
    xbo = npp.zeros(100, float)
    ybo = npp.zeros(100, float)
    arcbo = npp.zeros(100, float)
    i = 0
    npoch = npp.zeros(100, int)
    pai = 0.00
    rad = 0.00
    arc = 0.00
    pai = 2 * asin(1.0)
    rad = 300
    arc = 0
    arcbo[0] = 15
    arcbo[1] = rad * (pai / 4) - 30
    arcbo[2] = 30
    arcbo[3] = rad * (pai / 4) - 30
    arcbo[4] = 30
    arcbo[5] = rad * (pai / 4) - 30
    arcbo[6] = 30
    arcbo[7] = rad * (pai / 4) - 30
    arcbo[8] = 30
    arcbo[9] = rad * (pai / 4) - 30
    arcbo[10] = 30
    arcbo[11] = rad * (pai / 4) - 30
    arcbo[12] = 30
    arcbo[13] = rad * (pai / 4) - 30
    arcbo[14] = 30
    arcbo[15] = rad * (pai / 4) - 30
    arcbo[16] = 15
    npoch[0] = 2
    npoch[1] = 1
    npoch[2] = 2
    npoch[3] = 1
    npoch[4] = 2
    npoch[5] = 1
    npoch[6] = 2
    npoch[7] = 1
    npoch[8] = 2
    npoch[9] = 1
    npoch[10] = 2
    npoch[11] = 1
    npoch[12] = 2
    npoch[13] = 1
    npoch[14] = 2
    npoch[15] = 1
    npoch[16] = 2
    for i in range(0, 18):
        xbo[i] = rad * cos(arc / rad)
        ybo[i] = rad * sin(arc / rad)
        arc = arc + arcbo[i]
    return xbo, ybo, arcbo, npoch


# 分段连续边界，弧长arc和直角坐标x，y的关系，需改动
def arcn_to_xy(arc):
    global arcboun
    rad = 0.00
    rad = 300
    x = rad * cos((arc + arcboun) / rad)
    y = rad * sin((arc + arcboun) / rad)
    return x, y


# ************************************************
def savedata(npoch, boungridt_file, bounpoch_file):  # npoch=4
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    i = 0
    for i in range(0, np - 1):
        boungridt_file.write(str(px[i]) + 5 * " " + str(py[i]) + "\n")
        bounpoch_file.write(str(npoch) + "\n")


# *******************************
def savendata(npoch, boungridt_file, bounpoch_file):
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    i = 0
    for i in range(0, np):
        boungridt_file.write(str(px[i]) + 5 * " " + str(py[i]) + "\n")
    if npoch == 1:
        np = np - 1
    else:
        np = np + 1
    for i in range(0, np):
        bounpoch_file.write(str(npoch) + "\n")


def arcscatter(x0, y0, x1, y1, arc):
    # pai = 2 * asin(1.0) x0 = 1.0 y0 = 0 x1 = 1.0 y1 = 0 arc = 2 * pai
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    ee = 0
    arcxy = 0.00
    arc0 = 0.00
    arcstep = 0.00
    background_grid()
    arcstep = asin(1.0) * 10
    arc0 = 0.0
    np = 1  # 点的个数
    nbs = 1
    px[np - 1] = x0  # x0 = 1.0 y0 = 0
    py[np - 1] = y0
    x = px[np - 1]
    y = py[np - 1]
    ee, arcxy = search_boun(x, y, ee, arcxy)  # 判断边界点是不是合理 ee：面积最小背景单元编号 arcxy：step，新插入点尺度参数
    pstep[np - 1] = 0  # np = 1，这条弧剩下的长度
    arc0 = arc0 + arcxy  # arc0 = arcxy ，下一个点的位置
    while True:
        if (arc - arc0) <= 0:  # arc <= arc0 , arc = 2*Π , 新点超出这条弧的区域延伸到下一条弧
            arc0 = arc
            x = x1  # x1 = 1.0 y1 = 0
            y = y1
            ee, arcxy = search_boun(x, y, ee, arcxy)  # 判断边界点是不是合理 ee：面积最小背景单元编号 arcxy：step新插入点尺度参数
            if (arcxy * 1.5) > (arc0 - pstep[np - 1]):  # pstep这条弧剩下的长度
                break
        elif (arc - arc0) < arcxy:  # 还在原来的弧上，但是该点离弧终点近
            arc0 = (arc + arc0 - arcxy) / 2
            transfer_xy(arc0)  # 将新点位置转换为xy坐标
            ee, arcxy = search_boun(x, y, ee, arcxy)  # 计算新点被包含的背景网格单元和新点步长
        else:
            transfer_xy(arc0)  # 既没有超出该弧区域，有没有离终点太近
            ee, arcxy = search_boun(x, y, ee, arcxy)
        if (arcxy * 1.5) <= (arc0 - pstep[np - 1]):  # 新点位置减去上一个点位置即为上一个点的步长
            arc0 = arc0 - arcxy * 0.5
        else:
            np = np + 1
            px[np - 1] = x
            py[np - 1] = y
            pstep[np - 1] = arc0
            arc0 = arc0 + arcxy  # 下一个点位置
    np = np + 1  # 最后加一个点作为终点
    px[np - 1] = x1
    py[np - 1] = y1
    pstep[np - 1] = arc


# *****************************************
def background_grid():
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    i = 0
    j = 0
    bgn = open("bgnod2.dat", "r+")
    bgno = open("bgnoe.dat", "r+")
    bgg = open("bggridt2.dat", "r+")
    bgs = open("bgstep.dat", "r+")
    bgncount = 0
    bggcount = 0
    for line in bgg.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        bgpx[bggcount] = current_line[0]
        bgpy[bggcount] = current_line[1]
        bggcount += 1
    bgnp = bggcount  # 背景网格点的数量
    bgscount = 0
    for line in bgs.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        bgstep[bgscount] = current_line[0]
        bgscount = bgscount + 1  # 不知道是什么
    for line in bgn.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        bgnod[0, bgncount] = current_line[0]
        bgnod[1, bgncount] = current_line[1]
        bgnod[2, bgncount] = current_line[2]
        bgncount = bgncount + 1
    bgne = bgncount  # 背景网格三角形数量
    bgn0count = 0
    for line in bgno.readlines():  # 每个背景网格三角形每条边相邻的三角形编号
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        bgnoe[0, bgn0count] = current_line[0]
        bgnoe[1, bgn0count] = current_line[1]
        bgnoe[2, bgn0count] = current_line[2]
        bgn0count += 1
    bgn.close()
    bgno.close()
    bgg.close()
    bgs.close()
    for i in range(0, bgne):
        for j in range(0, 3):
            gx[j] = bgpx[bgnod[j, i] - 1]
            gy[j] = bgpy[bgnod[j, i] - 1]
        circumcircle()  # 计算该三角形外接圆圆心坐标及其半径
        bgcx[i] = gcx
        bgcy[i] = gcy
        bgcr[i] = gcr


# ***********************************
def search_boun(xx, yy, result, step):  # xx=1.0, yy=0, result=0, step=0
    ''' 和delaundry函数有点区别'''
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    i = 0
    j = 0
    k = 0
    nbg = npp.zeros(3, int)
    arr = npp.zeros([3, 3], float)
    dcp = 0
    bulk = npp.zeros(3, float)
    v1 = 0
    v2 = 0
    v4 = 0
    vol = npp.zeros(3, float)
    vmin = 1.0
    emin = 0
    for i in range(bgne - 1, -1, -1):  # bgne：背景网格三角形数量
        nbg[0] = bgnoe[0, i]
        nbg[1] = bgnoe[1, i]
        nbg[2] = bgnoe[2, i]
        # if nbg[0] != 0 and nbg[1] != 0 and nbg[2] != 0:
        # continue
        gcx = bgcx[i]
        gcy = bgcy[i]
        gcr = bgcr[i]
        dcp = npp.sqrt((xx - gcx) ** 2 + (yy - gcy) ** 2)
        #  if((dcp-gcr)<gcr*fm2*1000)then
        for k in range(0, 3):
            for j in range(0, 3):
                arr[0, j] = 1
                arr[1, j] = bgpx[bgnod[j, i] - 1]
                arr[2, j] = bgpy[bgnod[j, i] - 1]
            arr[1, k] = xx
            arr[2, k] = yy
            bulk[k] = det(arr, 3)  # 计算行列式的值
        v1 = bulk[0] + bulk[1] + bulk[2]
        v2 = abs(bulk[0]) + abs(bulk[1]) + abs(bulk[2])
        v4 = (v2 - v1) / v1
        # if(v4<1e-4)then
        if vmin > v4:  # vmin=1.0，即v1>v2/2时
            vmin = v4
            emin = i + 1  # 找到面积最小背景单元编号，即包含该点的背景网格单元
            vol[0] = bulk[0]
            vol[1] = bulk[1]
            vol[2] = bulk[2]
    if emin == 0:  # 没找到
        print('error')
        return result, step
    v1 = vol[0] + vol[1] + vol[2]
    v2 = abs(vol[0]) + abs(vol[1]) + abs(vol[2])
    v4 = v1 + v2
    step = 0
    for k in range(0, 3):  # 计算新插入点的尺度参数
        vol[k] = (vol[k] + abs(vol[k])) / v4
        step = step + bgstep[bgnod[k, emin - 1] - 1] * vol[k]  # bgstep取自文件
    result = emin
    return result, step


# *****************************************
def circumcircle():
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    i = 0
    arr = npp.zeros([2, 2], float)
    ver = npp.zeros(2, float)
    verx = npp.zeros(2, float)
    x1 = 0.0
    x2 = 0.0
    x3 = 0.0
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    for i in range(0, 2):
        arr[i, 0] = gx[i] - gx[2]
        arr[i, 1] = gy[i] - gy[2]
        ver[i] = 0.5 * (arr[i, 0] ** 2 + arr[i, 1] ** 2)
    axeqb(arr, ver, 2, verx)  # 计算圆心和半径，解那个方程用的
    gcx = verx[0] + gx[2]
    gcy = verx[1] + gy[2]

    x1 = gx[0] - gcx
    x2 = gy[0] - gcy
    x3 = npp.sqrt(x1 * x1 + x2 * x2)
    r1 = x3
    x1 = gx[1] - gcx
    x2 = gy[1] - gcy
    x3 = npp.sqrt(x1 * x1 + x2 * x2)
    r2 = x3
    r3 = npp.sqrt(verx[0] ** 2 + verx[1] ** 2)
    gcr = (r1 + r2 + r3) / 3


# *************************************
def axeqb(arr1, verb1, num, verx1):
    arr = npp.zeros([num, num], float)
    verb = npp.zeros(num, float)
    verx = npp.zeros(num, float)
    i = 0
    j = 0
    k = 0
    real_1 = 0.0
    real_2 = 0.0
    for i in range(0, num):
        for j in range(0, num):
            arr[i, j] = arr1[i, j]
        verb[i] = verb1[i]

    for i in range(0, num - 1):
        real_1 = arr[i, i]
        for j in range(i + 1, num):
            if abs(arr[j, i]) > abs(real_1):
                for k in range(i, num):
                    real_1 = arr[i, k]
                    arr[i, k] = arr[j, k]
                    arr[j, k] = real_1
                real_1 = verb[i]
                verb[i] = verb[j]
                verb[j] = real_1
            real_1 = arr[i, i]
        real_1 = arr[i, i]
        arr[i, i] = 1.0
        if abs(real_1) < (10 ** -6):
            print("error1:3 point in 1 line")
            print(i, real_1)
            k = input()
            return
        for j in range(i + 1, num):
            arr[i, j] = arr[i, j] / real_1
        verb[i] = verb[i] / real_1
        for j in range(i + 1, num):
            real_2 = arr[j, i]
            arr[j, i] = 0.0
            for k in range(i + 1, num):
                arr[j, k] = arr[j, k] - real_2 * arr[i, k]
            verb[j] = verb[j] - real_2 * verb[i]
    if abs(arr[num - 1, num - 1]) < 10 ** -6:
        print("error2:3 point in 1 line")
        print(arr[num, num])
        k = input()
        return
    verx[num - 1] = verb[num - 1] / arr[num - 1, num - 1]
    for i in range(num - 2, -1, -1):
        verx[i] = verb[i]
        for j in range(num - 1, i, -1):
            verx[i] = verx[i] - arr[i, j] * verx[j]

    for i in range(0, num):
        verx1[i] = verx[i]


# ****************************************
def det(matri, num):  # ?????
    '''I: num 数组的元素个数'''
    matrix = npp.zeros([num, num], float)
    i = 0
    j = 0
    k = 0
    sign_1 = 0
    det_matrix = 0.0
    real_1 = 0.0
    real_2 = 0.0
    for i in range(0, num):
        for j in range(0, num):
            matrix[j, i] = matri[i, j]
    det_matrix = 1
    sign_1 = 1
    for i in range(0, num - 1):
        real_1 = matrix[i, i]
        for j in range(i + 1, num):
            if abs(matrix[j, i]) > abs(real_1):
                for k in range(i, num):
                    real_1 = matrix[i, k]
                    matrix[i, k] = matrix[j, k]
                    matrix[j, k] = real_1
                real_1 = matrix[i, i]
                sign_1 = sign_1 * (-1)
        real_1 = matrix[i, i]
        if abs(real_1) < 10 ** (-15):
            return 0.0
        matrix[i, i] = 1.0
        det_matrix = det_matrix * real_1
        for j in range(i + 1, num):
            real_2 = matrix[j, i] / real_1
            matrix[j, i] = 0.0
            for k in range(i + 1, num):
                matrix[j, k] = matrix[j, k] - real_2 * matrix[i, k]
    det_matrix = det_matrix * matrix[num - 1, num - 1] * sign_1
    return det_matrix


# *********************************
def loadboun(filename):
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    i = 0
    fil = open(filename, "r")
    for line in fil.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        if i > 0:
            bx[i] = current_line[0]
            by[i] = current_line[1]
            bs[i] = current_line[2]
        else:
            nboun = current_line[0]
        i = i + 1
    fil.close()


# ***********************************
def transfer_xy(arc):  # arc = (arc + arc0 - arcxy) / 2 , 插入的平均步长
    global boungridts, bounpochs
    global np, nboun, sidetype, nbs, bounarc, px, py, pstep, bx, by, bs
    global gx, gy, gcx, gcy, gcr, x, y, xce, yce
    global bgnp, bgne, bgstep, bgpx, bgpy, bgcx, bgcy, bgcr, bgnod, bgnoe
    global arcboun
    ratio = 0.00
    i = 0
    n = 0
    # sidetype=1 离散形式的边界
    if sidetype == 1:
        if bs[nbs - 1] >= arc:
            n = 1
            for i in range(nbs - 1, -1, -1):
                if bs[i] < arc:
                    n = i + 1
                    break
        elif bs[nbs - 1] < arc:
            n = nboun
            for i in range(nbs - 1, nboun):
                if bs[i] >= arc:
                    n = i
                    break
        nbs = n
        if n == 1 or n == nboun:
            x = bx[n - 1]
            y = by[n - 1]
        else:
            ratio = (arc - bs[n - 1]) / (bs[n] - bs[n - 1])
            x = bx[n - 1] + ratio * (bx[n] - bx[n - 1])
            y = by[n - 1] + ratio * (by[n] - by[n - 1])
        # sidetype=2 连续边界，函数形式
    elif sidetype == 2:
        x, y = arc_to_xy(arc)
    elif sidetype == 3:
        x, y = arcn_to_xy(arc)


main()
