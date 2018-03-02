import numpy as np
import matplotlib.pylab as plt
import math
import Wvl2RGB
from Rays import *
from Ipsh import *

class Sphere:
    """
    A sphere object.
    """

    def __init__(self, origin, radius, pupil=None, material=None):
        self.origin = origin
        #self.origin.y = self.origin.y + radius if radius > 0 else self.origin.y - radius
        self.origin.y = self.origin.y + radius #if radius > 0 else self.origin.y - radius
        self.radius = radius
        self.pupil = pupil
        self.material = material
        if self.pupil is None:
            self.pupil = 2*self.radius
        elif self.pupil/2. > np.abs(self.radius):
            self.pupil = 2*self.radius

    def intersects(self, ray, depth):
        """
        If `ray` intersects sphere, return the distance at which it does;
        otherwise, `None`.
        """

        sphere_to_ray = ray.origin[depth] - self.origin
        b = 2 * ray.direction[depth] * sphere_to_ray
        c = sphere_to_ray ** 2 - self.radius ** 2
        discriminant = b ** 2 - 4 * c

        if discriminant >= 0:
            if self.radius > 0:
                dist = (-b - math.sqrt(discriminant)) / 2
            else:
                dist = (-b + math.sqrt(discriminant)) / 2
            if dist > 0:
                intersection_point = ray.point_at_dist(dist,depth)
                return dist

    def surface_norm(self, pt):
        """
        Return the surface normal to the sphere at `pt`.
        """
        out = (pt - self.origin).normalize() if self.radius > 0 else (self.origin-pt).normalize()
        return out

    def plot(self,ax=None):
        # plot only semicircle
        max_angle = np.arcsin(self.pupil/2./np.float(np.abs(self.radius)))
        theta = np.linspace(-max_angle,max_angle)
        theta = theta + np.pi if self.radius > 0 else theta
        y = np.abs(self.radius)*np.cos(theta)+self.origin.y
        z = np.abs(self.radius)*np.sin(theta)+self.origin.z
        return ax.plot(y,z,'b')[0]

class Plane:

    def __init__(self, origin, diameter, material='air', normal=Vector(0,-1,0), image=False):
        """
        Optical Object: Plane with circular diameter
        """
        self.origin = origin
        self.diameter = diameter
        self.radius = diameter/2.
        self.material = material
        if self.material == 'air':
            from refractiveIndex.refractiveIndex import *
            catalog = RefractiveIndex()
            self.material = catalog.getMaterial('main', 'H2', 'Peck')
        self.normal = normal
        self.image = image
        if self.image:
            self.spots = []
            self.wavelengths = []

    def intersects(self, ray, depth):
        """
        If `ray` intersects plane, return the distance at which it does;
        otherwise, `None`.
        """
        ndotu = self.normal * ray.direction[depth]
        if abs(ndotu) < 1.e-5:
            # No intersection or line is within plane
            return None

        if self.image:

        w = ray.origin[depth] - self.origin
        si = - (self.normal * w) / ndotu
        Psi = w + si * ray.direction[depth] + self.origin
        if np.sqrt(Psi.x**2+Psi.z**2) > self.radius:
            return None
        return (ray.origin[depth]-Psi).norm()

    def surface_norm(self, pt):
        """
        Return the surface normal to the plane at Point pt
        """
        return self.normal

    def _projection(self,point,projection=Vector(0,1,0)):
        costheta = (self.normal*projection)/(self.normal.norm()*projection.norm())
        axis = self.normal.unitcross(projection)

    def plot(self,ax=None):
        x = self.radius
        projection = Vector(0,1,0)
        return ax.plot([self.origin.y,self.origin.y],[-self.radius,self.radius],'b')[0]

