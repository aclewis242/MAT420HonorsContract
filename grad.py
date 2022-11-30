import numpy as np
import matplotlib.colors as clrs

def makegrad(pts: np.ndarray, func, color: str='green'):
    color_hex = clrs.CSS4_COLORS[color]
    pts_grad = np.array([])
    minN = -9e9999
    maxN = -minN
    for p in pts:
        gradVal = func(p[0],p[1])
        pts_grad = np.append(pts_grad, gradVal)
        if gradVal > minN:
            minN = gradVal
        if gradVal < maxN:
            maxN = gradVal
    pts_grad = hex(int(255*np.subtract(pts_grad, minN)/(maxN - minN))).replace('0x', '')
    return np.core.defchararray.add(color) # unfinished!