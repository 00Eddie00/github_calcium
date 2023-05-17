#-*- coding: UTF-8 -*-
import turtle as t

fp = open("bggrid.dat", 'rb')
data = fp.readlines()
x = []
y = []
for num in data:
    x.append(float(num.split()[0]))
    y.append(float(num.split()[1]))
fp=open("nod.dat", 'rb')
data = fp.readlines()
p_1 = []
p_2 = []
p_3 = []
for num in data:
    p_1.append(int(num.split()[0]))
    p_2.append(int(num.split()[1]))
    p_3.append(int(num.split()[2]))
fp.close()
t.screensize(800,800)
t.setup(width=0.9, height=0.9)
t.tracer(False)
for i in range(len(p_2)):
    t.penup()
    t.goto(x[p_1[i]-1], y[p_1[i]-1])
    t.pendown()
    t.goto(x[p_2[i]-1], y[p_2[i]-1])
    t.goto(x[p_3[i]-1], y[p_3[i]-1])
    t.goto(x[p_1[i]-1], y[p_1[i]-1])
t.hideturtle()
ts = t.getscreen()
ts.getcanvas().postscript(file="jsr-ryr-center-move.eps")
# t.done()