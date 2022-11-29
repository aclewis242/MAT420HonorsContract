import numpy as np
import matplotlib.pyplot as plt
import cyl
import math

c = cyl.Cylinder()
rs = np.linspace(0, 5, 100)
def dyp(t, vs): return -math.sin(t)
def dy(t, vs): return vs[1]
vs = c.rk4(np.array([[0, 1]]), np.array([dy, dyp]), rs)
vas = np.zeros_like(rs)
for ri in range(rs.size): vas[ri] = math.sin(rs[ri])
plt.plot(rs, vs[:,0])
plt.plot(rs, vas)
plt.show()