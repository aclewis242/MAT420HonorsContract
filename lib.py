import numpy as np
import matplotlib.colors as clrs
from math import *
from PIL import ImageColor, Image

def rk4old(d1f, d2f, xs, y, y1):
    # To use: d1f is 1st der., d2f is 2nd. y is initial value of f, y1 is initial value of d1f
    ys = np.array([y])
    y1s = np.array([y1])
    h = xs[1] - xs[0]
    for xn in xs[1:]:
        yn = ys[-1]
        y1n = y1s[-1]
        k11 = h*d1f(xn, yn, y1n)
        k12 = h*d2f(xn, yn, y1n)
        k21 = h*d1f(xn+h/2, yn+k11/2, y1n+k12/2)
        k22 = h*d2f(xn+h/2, yn+k11/2, y1n+k12/2)
        k31 = h*d1f(xn+h/2, yn+k21/2, y1n+k22/2)
        k32 = h*d2f(xn+h/2, yn+k21/2, y1n+k22/2)
        k41 = h*d1f(xn+h, yn+k31, y1n+k32)
        k42 = h*d2f(xn+h, yn+k31, y1n+k32)
        ys = np.append(ys, yn + (k11 + 2*k21 + 2*k31 + k41)/6)
        y1s = np.append(y1s, y1n + (k12 + 2*k22 + 2*k32 + k42)/6)
    return ys

def rk4(vars: np.ndarray, funcs: np.ndarray, ts: np.ndarray, refl):
    h = ts[1] - ts[0]
    while ts.any():
        t = ts[0]
        ts = np.delete(ts, 0)
        vns = np.array([])
        for v in vars[-1]:
            vns = np.append(vns, v)
        fns = np.array([])
        for f in funcs:
            fns = np.append(fns, f)
        k1s = np.array([])
        for i in range(funcs.size):
            k1s = np.append(k1s, h*funcs[i](t, vns))
        k2s = np.array([])
        for i in range(funcs.size):
            k2s = np.append(k2s, h*funcs[i](t+h/2, np.add(vns, k1s/2)))
        k3s = np.array([])
        for i in range(funcs.size):
            k3s = np.append(k3s, h*funcs[i](t+h/2, np.add(vns, k2s/2)))
        k4s = np.array([])
        for i in range(funcs.size):
            k4s = np.append(k4s, h*funcs[i](t+h, np.add(vns, k3s)))
        k12 = np.add(k1s, 2*k2s)
        k34 = np.add(2*k3s, k4s)
        ks = np.add(k12, k34)/6
        vars = np.append(vars, [np.add(vns, ks)], 0)
        if refl and sqrt(vars[-1,0]**2 + vars[-1,1]**2) >= refl:
            return vars, ts, refl
    return vars, ts, 0

def intg(vars: np.ndarray, funcs: np.ndarray, ts: np.ndarray, refl):
    vars, ts, refl = rk4(vars, funcs, ts, refl)
    while refl:
        v = [vars[-1,2], vars[-1,3]]
        l = [vars[-1,0], vars[-1,1]]
        l /= np.linalg.norm(l)
        [vars[-1,2], vars[-1,3]] = np.subtract(v, 2*np.dot(l, v)*l)
        vars, ts, refl = rk4(vars, funcs, ts, refl)
    return vars

def makegrad(pts: np.ndarray, func, color: str='black'):
    color_hex = clrs.CSS4_COLORS[color]
    pts_grad = np.array([])
    minN = -9e9999
    maxN = -minN
    maxX = 0
    maxY = 0
    for p in pts:
        gradVal = func(p.x,p.y)
        pts_grad = np.append(pts_grad, gradVal)
        if gradVal > minN: minN = gradVal
        if gradVal < maxN: maxN = gradVal
        if abs(p.x) > maxX: maxX = abs(p.x)
        if abs(p.y) > maxY: maxY = abs(p.y)
    pts_grad = castHex(255*np.subtract(pts_grad, minN)/(maxN - minN))
    dim = int(sqrt(np.size(pts_grad)))
    img = Image.new('RGBA', (dim,dim), color='#ffffff00')
    for i in range(pts.size):
        p = pts[i]
        imgX = int((p.x/maxX)*(dim-1)/2 + (dim-1)/2)
        imgY = int((p.y/maxY)*(dim-1)/2 + (dim-1)/2)
        img.putpixel((imgX,imgY), ImageColor.getrgb(color_hex+pts_grad[i]))
    return img

def castHex(arr: np.ndarray):
    return [fillPlaces(hex(int(v)).replace('0x', '')) for v in arr]

def fillPlaces(s: str):
    if len(s) == 1:
        return '0' + s
    return s