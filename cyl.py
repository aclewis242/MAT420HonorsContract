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
    
    def n_simple(self, x, y):
        if x < 2.4:
            return 1.2
        elif x >= 2.4 and x <= 2.8:
            return x**2/4.8
        else:
            return 1.63
    
    def dn_simple(self, x, y):
        if x >= 2.4 and x <= 2.8:
            return x/2.4
        else:
            return 0
    
    def n_const(self, x, y):
        return 1.3
    
    def dn_const(self, x, y):
        return 0