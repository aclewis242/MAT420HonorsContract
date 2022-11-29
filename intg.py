import math

class intg:
    n = None
    def __init__(self, n):
        self.n = n
    def __truediv__(self, n2):
        return abs(n2.n/self.n) + 1