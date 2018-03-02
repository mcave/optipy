import numpy as np
import matplotlib.pylab as plt
import math
import Wvl2RGB

class Vector:
    """
    A generic 3-element vector. All of the methods should be self-explanatory.
    """

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def from_numpy(self, array):
        return self.__init__(array[0],array[1],array[2])

    def norm(self):
        return math.sqrt(sum(num * num for num in self))

    def normalize(self):
        return self / self.norm()

    def reflect(self, other):
        other = other.normalize()
        return self - 2 * (self * other) * other

    def refract(self, other, r):
        """
        Refration rule: -> r=n1/n2, other = normalized normal vector
        """
        other = other.normalize()
        c = - np.dot(other,self)
        out = r*self + (r*c-np.sqrt(1-r**2*(1-c**2)))*other
        return out

    def cross(self,vector):
        return Vector.from_numpy(np.cross(self.to_numpy(),vector.to_numpy()))

    def unitcross(self,vector):
        return self.cross(vector).normalize()

    def __str__(self):
        return "Vector({}, {}, {})".format(*self)

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        if isinstance(other, Vector):
            return self.x * other.x + self.y * other.y + self.z * other.z;
        else:
            return Vector(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return Vector(self.x / other, self.y / other, self.z / other)

    def __div__(self, other):
        return Vector(self.x / other, self.y / other, self.z / other)

    def __pow__(self, exp):
        if exp != 2:
            raise ValueError("Exponent can only be two")
        else:
            return self * self

    def __iter__(self):
        yield self.x
        yield self.y
        yield self.z

    def spherical(self):
        theta = np.arccos(self.z/self.norm())
        phi = np.arctan2(self.y,self.x)
        return self.norm(),theta,phi

    def to_numpy(self):
        return np.array([self.x,self.y,self.z])

Point = Vector

class Ray:
    """
    A mathematical ray.
    """

    def __init__(self, origin, direction, wavelength=566, default_legth=0):
        if isinstance(origin, list):
            self.origin = origin
            self.direction = direction
        else:
            self.origin = []
            self.origin.append(origin)
            self.direction = []
            self.direction.append(direction)
        self.depth = len(self.origin)
        self.wavelength = wavelength
        self.end = []
        self.end.append(self.origin[0])

    def point_at_dist(self, dist, depth):
        return self.origin[depth] + self.direction[depth] * dist

    def refract(self, norm, r):
        new_direction = self.direction[self.depth-1].refract(norm,r)
        self.direction.append(new_direction.normalize())
        self.depth += 1

    def plot(self,ax=None):
        for jdepth in xrange(self.depth-1):
            line = ax.plot([self.origin[jdepth].y,self.origin[jdepth+1].y],\
                    [self.origin[jdepth].z,self.origin[jdepth+1].z],\
                    color=Wvl2RGB.wavelength_to_rgb(self.wavelength))
        return line[0]

class RayFan:
    """
    A fan of mathematical rays
    """

    def __init__(self, open_angle=np.radians(10), origin=Point(0,0,0), direction=Vector(0,1,0),\
              wavelength=566,n_theta=5, n_phi=5):
        self.open_angle = open_angle
        self.origin = origin
        self.direction = direction.normalize()
        self.rays = []
        _,theta0,phi0 = self.direction.spherical()
        # Theta and Phi follw spherical coordinates
        # n_theta and n_phi have to be odd to have a central ray
        n_theta=n_theta+1 if n_theta%2 == 0 else n_theta
        n_phi=n_phi+1 if n_phi%2 == 0 else n_phi
        dtheta = 0 if n_theta == 1 else self.open_angle/(n_theta-1)
        dphi = 0 if n_phi == 1 else self.open_angle/(n_phi-1)
        for jtheta in xrange(-(n_theta-1)/2,(n_theta-1)/2+1):
            theta = theta0-dtheta*jtheta
            for jphi in xrange(-(n_phi-1)/2,(n_phi-1)/2+1):
                phi = phi0-dphi*jphi
                # Spherical -> Cartesian
                x = np.sin(theta)*np.cos(phi)
                y = np.sin(theta)*np.sin(phi)
                z = np.cos(theta)
                direction = Vector(x,y,z)
                self.rays.append(Ray(self.origin,direction,wavelength=wavelength))

    def append(self,rays):
        for ray in rays:
            self.rays.append(ray)

    def point_at_dist(self, dist):
        return self.origin + self.direction * dist

    def plot(self):
        for ray in self.rays:
            ray.plot()
