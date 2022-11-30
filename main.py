import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import cyl
from lib import *

if __name__ == "__main__":
    c = cyl.Cylinder()
    def dx(t, vs): return vs[2]
    def dy(t, vs): return vs[3]
    def dxp(t, vs, n=c.n, dn=c.dn):
        x = vs[0]
        y = vs[1]
        xp = vs[2]
        yp = vs[3]
        return -dn(x,y)*xp**2/n(x,y) - 2*dn(x,y)*xp*yp/n(x,y) + dn(x,y)*yp**2/n(x,y)
    def dyp(t, vs, n=c.n, dn=c.dn):
        x = vs[0]
        y = vs[1]
        xp = vs[2]
        yp = vs[3]
        return dn(x,y)*xp**2/n(x,y) - 2*dn(x,y)*xp*yp/n(x,y) - dn(x,y)*yp**2/n(x,y)

    x0 = 0; y0 = 2
    xp0 = -1; yp0 = 0

    ### ROTATING CYLINDER ###
    ts = np.linspace(0, 100, 5000)
    vs = intg(np.array([[x0, y0, xp0, yp0]]), np.array([dx, dy, dxp, dyp]), ts, c.R)
    plt.plot(vs[:,0], vs[:,1])
    plt.gca().add_artist(ptch.Circle((0, 0), c.R, fill=False))
    plt.xlim([-c.R-1, c.R+1])
    plt.ylim([-c.R-1, c.R+1])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.gca().set_aspect('equal')
    plt.scatter(x0, y0, label='Start')
    plt.scatter(vs[-1,0], vs[-1,1], label='Finish', zorder=10)
    plt.title("Trajectory of light within rotating cylinder")
    plt.legend()
    plt.savefig("traj_rot.png", dpi=150)
    plt.show()

    ### STATIONARY CYLINDER ###
    ts = np.linspace(0, 100, 5000)
    vs = intg(np.array([[x0, y0, xp0, yp0]]), np.array([dx, dy, dxp, dyp]), ts, c.R)