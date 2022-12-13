import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import cyl
from point import Point as pt
from lib import *
from numpy.linalg import norm

# The control file - run this to produce output

if __name__ == "__main__":
    c = cyl.Cylinder()
    # Intended structure of 'vs': [var 1, var 2, ..., var 1 prime, var 2 prime, ...]
    def dx(t, vs): return vs[2]
    def dy(t, vs): return vs[3]
    def dxp(t, vs, n=c.n, dnx=c.dnx, dny=c.dny):
        x = vs[0]
        y = vs[1]
        xp = vs[2]
        yp = vs[3]
        return -dnx(x,y)*xp**2/n(x,y) - 2*dny(x,y)*xp*yp/n(x,y) + dnx(x,y)*yp**2/n(x,y)
    def dyp(t, vs, n=c.n, dnx=c.dnx, dny=c.dny):
        x = vs[0]
        y = vs[1]
        xp = vs[2]
        yp = vs[3]
        return dny(x,y)*xp**2/n(x,y) - 2*dnx(x,y)*xp*yp/n(x,y) - dny(x,y)*yp**2/n(x,y)

    x0 = 0; y0 = 2
    xp0 = -1; yp0 = 1

    ts_base = np.linspace(0, 50, 5000)
    ### ROTATING CYLINDER ###
    ts = ts_base
    h = 5e-2
    pts_base = np.mgrid[-c.R:c.R+h:h, -c.R:c.R+h:h].reshape(2,-1).T
    pts = np.array([])
    for p in pts_base:
        if norm(p) < c.R: pts = np.append(pts, pt(p[0], p[1]))
    img = makegrad(pts, c.n)
    vs = intg(np.array([[x0, y0, xp0, yp0]]), np.array([dx, dy, dxp, dyp]), ts, c.R)
    plt.imshow(img, extent=[-c.R, c.R, -c.R, c.R])
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
    ts = ts_base
    n = c.n_const
    dn = c.dn_const
    def dxpc(t, vs): return dxp(t, vs, n=n, dnx=dn, dny=dn)
    def dypc(t, vs): return dyp(t, vs, n=n, dnx=dn, dny=dn)
    vs = intg(np.array([[x0, y0, xp0, yp0]]), np.array([dx, dy, dxpc, dypc]), ts, c.R)
    plt.plot(vs[:,0], vs[:,1])
    plt.gca().add_artist(ptch.Circle((0, 0), c.R, fill=False))
    plt.xlim([-c.R-1, c.R+1])
    plt.ylim([-c.R-1, c.R+1])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.gca().set_aspect('equal')
    plt.scatter(x0, y0, label='Start')
    plt.scatter(vs[-1,0], vs[-1,1], label='Finish', zorder=10)
    plt.title("Trajectory of light within stationary cylinder")
    plt.legend()
    plt.savefig("traj_stat.png", dpi=150)
    plt.show()

    ### SIMPLE CASE ###
    ts = np.linspace(0, 5, 100)
    n = c.n_simple
    dn = c.dn_simple
    def dxps(t, vs): return dxp(t, vs, n=n, dnx=dn, dny=dn)
    def dyps(t, vs): return dyp(t, vs, n=n, dnx=dn, dny=dn)
    vs = intg(np.array([[0, 2, 1, -1]]), np.array([dx, dy, dxps, dyps]), ts, 0)
    xs, ys = vs[:,0], vs[:,1]
    plt.plot(xs, ys)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.gca().set_aspect('equal')
    plt.scatter(0, 2, label='Start')
    print(f"arctan: {np.arctan((vs[-1,1] - vs[-2,1])/(vs[-1,0] - vs[-2,0]))}")
    plt.scatter(vs[-1,0], vs[-1,1], label='Finish', zorder=10)
    plt.title("Trajectory of light within simple medium")
    plt.legend()
    plt.savefig("traj_simp.png", dpi=150)
    plt.show()