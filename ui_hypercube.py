import re
import os
import math
import turtle
import argparse
import subprocess
from tkinter import *
from functools import lru_cache
from turtle import Turtle, Screen
from typing import Any, Dict, Tuple, List


@lru_cache(maxsize=None)
def coords(v: int, dim: int)-> Tuple[int, int]:
    x,y = 0, -1 * (v & 1)
    for i in range(1, dim):
        if (v & (1<<i)):
            a =  90 * (dim-i)/(dim-1)
            x += math.sin(a) * 100
            y -= math.cos(a) * 100
    return x,y

def print_results(t: turtle.Turtle, output: Dict[str, Any]):
    t.penup()
    t.setheading(90)
    t.forward(150)
    t.setheading(180)
    t.forward(25)
    t.write(output["data"], font=('Courier',20,'bold'), align='center')
    t.hideturtle()


def draw_cube(t: turtle.Turtle, args: Dict[str,int])-> None:
    dim = args['n_dimensions']
    t.penup(); t.shape('circle')
    t.shapesize(0.05); t.speed(0); t.setpos(-10,-10)
    cs = {i: coords(i*2,dim) for i in range(1<<dim)}

    for i,v in sorted(cs.items()):
        t.color("black"); dotsize = 3
        if i == args['start' ]: t.color("blue");  dotsize = 10
        if i == args['finish']: t.color("red"); dotsize = 10
        t.goto(v); t.dot(dotsize)

    for i in range(1<<dim):
        for s in range(dim):
            t.goto(cs[i]); t.pendown()
            t.goto(cs[i^(1<<s)]); t.penup()

@lru_cache(maxsize=None)
def get_output(es: Entry)-> Dict[str, Any]:
    keys = ['n_dimensions', 'start', 'finish', 'n_threads']
    args = {keys[i]: int(e[1].get()) for i,e in enumerate(es)}
    args['n_dimensions'] += 1
    t = subprocess.Popen(f"./scripts/run_hypercube {args['start']} {args['finish']} {args['n_threads']} {args['n_dimensions']} ",
                        shell=True, stdout=subprocess.PIPE).stdout.read()
    d = {"data": t.strip().decode(), "nums": list(map(int, re.findall(rb'\d+', t)))}

    return args, d

def hypercube(es: Entry)-> None:
    root = Tk()

    args, output = get_output(tuple(es))

    screen = Screen()
    screen.clear()
    t = Turtle()
    print_results(t, output)
    draw_cube(t, args)
    t.hideturtle()

@lru_cache(maxsize=None)
def makeform(root, fields):
    entries = list()
    l = Label(root, width=30, text="ANT Inc. Interdimensional Taxi Services", anchor='w')
    l.pack(side=TOP)
    for i,f in enumerate(fields):
        row = Frame(root)
        l = Label(row, width=15, text=f, anchor='w')
        e = Entry(row)
        l.pack(side=LEFT)
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        e.pack(side=RIGHT, expand=YES, fill=X)
        entries.append((f,e))
    return entries

if __name__ == '__main__':
    root = Tk()
    es = makeform(root, ('Which dimension?', 'Entry node index?', 'Exit node index?', 'How many threads?'))
    b = Button(root, text='End', command=root.quit)
    b.pack(side=BOTTOM, padx=5, pady=5)
    c = Button(root, text='Compute', command=(lambda es=es : hypercube(es)))
    c.pack(side=BOTTOM, padx=5, pady=5)
    root.mainloop()

