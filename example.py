"""
--- optipy ---

A simple optical raytracer

"""

import math
import numpy as np
import logging
from lib import *
from refractiveIndex.refractiveIndex import *
catalog = RefractiveIndex()
import logging

logging.basicConfig(level=logging.WARNING)

if __name__ == "__main__":
    # Define a fan of mathematical rays with two different wavelengths
    rayfan = Rays.RayFan(origin=Point(0,0,1),open_angle=np.radians(10),n_phi=100,n_theta=100,wavelength=650)
    rayfan2 = Rays.RayFan(origin=Point(0,0,1),open_angle=np.radians(10),n_phi=100,n_theta=100,wavelength=400)
    # Append them together
    rayfan.append(rayfan2.rays)

    # Pick two material
    bk7 = catalog.getMaterial('glass', 'BK7', 'SCHOTT')
    air = catalog.getMaterial('main', 'H2', 'Peck')
    # Define two objects
    sp = Objects.Sphere(Point(0,10,0),3,pupil=5,material=bk7)
    #sp2 = Objects.Sphere(Point(0,12,0),3,pupil=5,material=air)
    sp2 = Objects.Plane(Point(0,12,0),5,material=air)
    sp3 = Objects.Plane(Point(0,20,0),5,image=True)

    # Define the scene and trace the rays
    sc = Render.Scene(rayfan.rays,[sp,sp2,sp3])
    sc.render()
    sc.plot()
    #sc.objects[-1].plot_spots()
