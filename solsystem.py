import turtle
from datetime import datetime

from PIL import Image


def setup(p, N, Nt, Ke, Pe, dt):
    scale = 7 * 10 ** 8
    w, h = 800, 800
    window = turtle.Screen()
    window.setup(w, h)
    window.bgcolor("black")
    window.colormode(255)
    Te = Ke + Pe
    orbit = turtle.Turtle()
    orbit.hideturtle()
    orbit.penup()
    orbit.speed(0)
    turtle.tracer(0, 0)
    write = turtle.Turtle()
    Ncolor = [[255, 204, 51], [219, 206, 202], [165, 124, 27], [79, 76, 176], [247, 234, 198], [188, 39, 50]]
    for j in range(Nt):
        write.color('white')
        style = ('Courier', 10, 'italic')
        write.penup()
        write.setposition(0, -200)
        write.pendown()
        write.write(
            'Sluttid:' + str(Nt * dt) + '\nTid:' + str(j * dt) + '\nAntal tidssteg:' + str(Nt) + '\nTidssteg:' + str(
                j) + '\nStorlek på tidsteg:' + str(dt) + '\nKinetisk energi:' + str(
                Ke[j]) + '\nPotentiell energi:' + str(Pe[j]) + '\nTotal energi:' + str(Te[j]), font=style,
            align='center')
        write.hideturtle()
        for i in range(N):
            orbit.setposition((p[i, 0, j]['pos'] / scale, p[i, 1, j]['pos'] / scale))
            orbit.color(Ncolor[i][0], Ncolor[i][1], Ncolor[i][2])
            orbit.pendown()
            orbit.setposition((p[i, 0, j + 1]['pos'] / scale, p[i, 1, j + 1]['pos'] / scale))
            orbit.penup()
        write.clear()
    date = (datetime.now()).strftime("%d%b%Y-%H%M%S")
    fileName = 'plot-' + date

    screen = turtle.Screen()
    screen.tracer(False)
    screen.tracer(True)
    canvas = screen.getcanvas()
    canvas.postscript(file=fileName + '.eps', width=w, height=h)
    img = Image.open(fileName + '.eps')
    img.save(fileName + '.jpg')
    turtle.clearscreen()
