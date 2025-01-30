from math import *

class Point:
    '''
    A simple class describing a point in space, able to store a color at that point. Used almost exclusively to create the color gradient for the cylinder.
    '''
    x = None
    y = None
    z = None
    color = None

    def __init__(self, x, y, z=0, color=None):
        self.x = x
        self.y = y
        self.z = z
        self.color = color
    
    def r(self):
        '''
        The point's distance from the vertical axis of the cylinder
        '''
        return sqrt(self.x**2 + self.y**2)