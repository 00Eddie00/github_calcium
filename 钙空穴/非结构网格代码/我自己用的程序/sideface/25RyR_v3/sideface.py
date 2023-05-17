import numpy as npp

def not_empty(s):
    return s and s.strip()


def main():
    '''
    输入：sidegridt2.dat -->  大圆，小圆上所有点的坐标，第1行是点数，小圆上45个点，大圆337个点
        sidenpoch2.dat -->  大圆，小圆上所有点的属性，多1行，多出来的是大圆上第一个入流点和最后一个入流点属性重复了
    :return:
    '''
    boun = npp.zeros(100, int)
    nboun = 26
    boun[0] = 11
    boun[1] = 23
    boun[2] = 35
    boun[3] = 47
    boun[4] = 56
    boun[5] = 69
    boun[6] = 81
    boun[7] = 93
    boun[8] = 105
    boun[9] = 119
    boun[10] = 131
    boun[11] = 143
    boun[12] = 155
    boun[13] = 167
    boun[14] = 180
    boun[15] = 192
    boun[16] = 204
    boun[17] = 216
    boun[18] = 228
    boun[19] = 240
    boun[20] = 250
    boun[21] = 264
    boun[22] = 278
    boun[23] = 291
    boun[24] = 300
    boun[25] = 773

    sideface = open("sideface.dat", "w")
    sidegridt = open("sidegridt.dat", "r")
    sidenpoch = open("sidenpoch.dat", "r")

    d = 0
    sidenpochs = npp.empty(1000, int)
    for line in sidenpoch.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        sidenpochs[d] = current_line[0]
        d = d + 1
    # nb = d
    # sideface.write(str(nb)+ "\n")   # sideface 第一行写：统计npoch有多少行
    istart = 0
    for i in range(0, nboun):
        iend = boun[i]
        for j in range(istart, iend - 1):
            k = sidenpochs[j]
            if k != sidenpochs[j + 1]:    # 线段两端一个是入流点，一个是壁面点，那么这个线段上的点都是壁面点
                k = 1
            sideface.write(str(j+1) + " " + str(j + 2) + " " + str(k)+"\n")
        k = sidenpochs[iend-1]
        if k != sidenpochs[istart]:        # 最后一个点和起点的比较
            k = 1
        sideface.write(str(iend) + " " + str(istart+1) + " " + str(k)+"\n")
        istart = iend
    sideface.close()
    sidegridt.close()
    sidenpoch.close()
main()