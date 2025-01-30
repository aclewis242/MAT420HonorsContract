from math import *

R = 8.3145 # Ideal gas constant (J/K*mol)
T = 273.15 # Standard temperature

# This class serves to define all the various equations for refractive index, not just those for the cylinder (despite the name).

class Cylinder:
    w = None # Rotational velocity of cylinder
    R = None # Radius of cylinder
    d0 = None # Characteristic density of fluid
    mc = None # Material constant (for index of refraction)

    def __init__(self, w: float=22, R: float=5, d0: float=1, mc: float=1):
        '''
        Initialises the cylinder.

        ### Parameters
        - `w`: Rotational velocity (rad/sec)
        - `R`: Radius of cylinder (m)
        - `d0`: Characteristic density of fluid, i.e. its density when stationary
        - `mc`: Material constant relating the fluid's density with its index of refraction

        *(The units of `d0` and `mc` are non-specific; they must simply agree with one another.)*
        '''
        self.w = w
        self.R = R
        self.d0 = d0
        self.mc = mc
    
    def rho(self, x, y):
        '''
        Returns the density of the fluid within the rotating cylinder at the given point.
        '''
        r = sqrt(x**2 + y**2)
        if r <= self.R: return self.d0*exp(((self.w*r)**2)/(2*R*T))
        else: return 0
    
    def n(self, x, y):
        '''
        Returns the refractive index of the fluid within the cylinder at the given point.
        '''
        return self.mc*self.rho(x,y) + 1
    
    def dnx(self, x, y):
        '''
        Returns the derivative of the fluid's refractive index with respect to x at the given point.
        '''
        return (self.w**2/(R*T))*x*self.mc*self.rho(x,y)
    
    def dny(self, x, y):
        '''
        Returns the derivative of the fluid's refractive index with respect to y at the given point.
        '''
        return (self.w**2/(R*T))*y*self.mc*self.rho(x,y)
    
    def n_simple(self, x, y):
        '''
        Produces the refractive index used to mimic basic medium-to-medium (no gradient) refraction.
        '''
        if x < 2.4: return 1.2
        elif x >= 2.4 and x <= 2.8: return x**2/4.8
        else: return 1.63
    
    def dn_simple(self, x, y):
        '''
        Produces the derivative of the refractive index used to mimic basic medium-to-medium (no gradient) refraction.
        '''
        if x >= 2.4 and x <= 2.8: return x/2.4
        else: return 0
    
    def n_const(self, x, y):
        '''
        Returns a constant refractive index, used in the stationary cylinder case.
        '''
        return self.mc*self.d0 + 1
    
    def dn_const(self, x, y):
        '''
        Returns the derivative of the constant refractive index in the stationary cylinder case, i.e. 0.
        '''
        return 0