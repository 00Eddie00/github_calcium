import numpy as npp
import math


def not_empty(s):
    return s and s.strip()


# 这个文件的作用是为每一个边界点定义一个属性，内部点是0，流进的边界点是1

def aa():
    npoch = npp.zeros(50000, dtype="int")
    sideface = npp.zeros(3, dtype="int")
    si_file = open("4RYRsidegridt.dat", "r")
    sicount = 0
    for line in si_file.readlines():
        sicount = sicount + 1
    nbp = sicount
    si_file.close()
    gri_file = open("4RYRgridt.dat", "r")
    gricount = 0
    for line in gri_file.readlines():
        gricount = gricount + 1
    np = gricount
    gri_file.close()
    for i in range(0, nbp):
        npoch[i] = 1

    for i in range(nbp, np):
        npoch[i] = 0

    sideface_file = open("4RYRsideface.dat", "r")
    sidefacecount = 0
    for line in sideface_file.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        sideface[0] = current_line[0]
        sideface[1] = current_line[1]
        sideface[2] = current_line[2]
        if sideface[2] != 1:
            npoch[sideface[0]-1] = sideface[2]
            npoch[sideface[1]-1] = sideface[2]
        sidefacecount=sidefacecount+1
    ns=sidefacecount
    sideface_file.close()
    npoch_file = open("4RYRnpoch.dat", "w")
    for i in range(0, np):
        npoch_file.write(str(npoch[i])+"\n")
    npoch_file.close()
aa()