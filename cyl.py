import numpy as np
import math

R = 8.3145 # Ideal gas constant
T = 273.15 # Standard temperature

class Cylinder:
    w = None # Rotational velocity of cylinder
    R = None # Radius of cylinder
    d0 = None # Characteristic density of fluid
    mc = None # Material constant (for index of refraction)

    def __init__(self, w: float=10, R: float=5, d0: float=1, mc: float=0.5):
        self.w = w
        self.R = R
        self.d0 = d0
        self.mc = mc
    
    def rho(self, x, y):
        r = math.sqrt(x**2 + y**2)
        return self.d0*math.exp(((self.w*r)**2)/(2*R*T))
    
    def n(self, x, y):
        return self.mc*self.rho(x, y) + 1
    
    def dn(self, x, y):
        return (self.w**2/(R*T))*self.mc*self.rho(x, y)
    
    def rk4old(self, d1f, d2f, xs, y, y1):
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

    def rk4(self, vars: np.ndarray, funcs: np.ndarray, ts: np.ndarray):
        h = ts[1] - ts[0]
        for t in ts[1:]:
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
            print(vars)
        return vars