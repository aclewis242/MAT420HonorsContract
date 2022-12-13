from math import *

R = 8.3145 # Ideal gas constant
T = 273.15 # Standard temperature

# This class serves to define all the various equations for refractive index, not just those for the cylinder (despite the name).

class Cylinder:
    w = None # Rotational velocity of cylinder
    R = None # Radius of cylinder
    d0 = None # Characteristic density of fluid
    mc = None # Material constant (for index of refraction)

    def __init__(self, w: float=22, R: float=5, d0: float=1, mc: float=1):
        self.w = w
        self.R = R
        self.d0 = d0
        self.mc = mc
    
    def rho(self, x, y): # Density of fluid within rotating cylinder
        r = sqrt(x**2 + y**2)
        if r <= self.R: return self.d0*exp(((self.w*r)**2)/(2*R*T))
        else: return 0
    
    def n(self, x, y): # Refractive index for rotating cylinder
        return self.mc*self.rho(x,y) + 1
    
    def dnx(self, x, y): # Derivative of refractive index with respect to x
        return (self.w**2/(R*T))*x*self.mc*self.rho(x,y)
    
    def dny(self, x, y):
        return (self.w**2/(R*T))*y*self.mc*self.rho(x,y) # The above, with respect to y
    
    def n_simple(self, x, y): # Refractive index mimicking basic medium-to-medium (no gradient) refraction
        if x < 2.4: return 1.2
        elif x >= 2.4 and x <= 2.8: return x**2/4.8
        else: return 1.63
    
    def dn_simple(self, x, y):
        if x >= 2.4 and x <= 2.8: return x/2.4
        else: return 0
    
    def n_const(self, x, y): # Constant refractive index (for stationary cylinder)
        return self.mc*self.d0 + 1
    
    def dn_const(self, x, y):
        return 0