import numpy as np
from math import *

class Point:
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
        return sqrt(self.x**2 + self.y**2)