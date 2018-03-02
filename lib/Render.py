import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Polygon
import logging
import math
import copy
from Objects import *
from Ipsh import *

from refractiveIndex.refractiveIndex import *
catalog = RefractiveIndex()
air = catalog.getMaterial('main', 'H2', 'Peck')

class Scene:
    """
    The scene that gets rendered. Contains information like the camera
    position, the different objects present, rays, etc..
    """

    def __init__(self, rays=None, objects=None, camera=None):
        self.rays = rays
        self.objects = objects

    def render(self):
        jray = 0
        for ray in self.rays:
            self.rays[jray] = self._trace_ray(copy.copy(ray))
            jray += 1

    def _trace_ray(self, ray, depth=0):

        intersection = self._get_intersection(ray,depth)
        logging.debug('Intersection: '+str(intersection)+' Depth: '+str(depth))
        if intersection is None:
            intersection_pt = ray.point_at_dist(10,depth)
            ray.origin.append(intersection_pt)
            ray.depth += 1
            return ray

        obj, dist = intersection
        intersection_pt = ray.point_at_dist(dist,depth)
        ray.origin.append(intersection_pt)
        # Check if object is a Camera
        if obj.image:
            # End propagation
            ray.depth += 1
            return ray
        # Ray starts from vacuum
        if depth == 0:
            n1 = 1.
        else:
            n1 = self.objects[depth-1].material.getRefractiveIndex(ray.wavelength)
        n2 = obj.material.getRefractiveIndex(ray.wavelength)
        ray.refract(obj.surface_norm(intersection_pt),n1/n2)
        return self._trace_ray(ray,depth+1)

    def _get_intersection(self, ray, depth):
        """
        If ray intersects any of `self.objects`, return `obj, dist` (the object
        itself, and the distance to it). Otherwise, return `None`.
        """
        intersection = None
        if depth >= len(self.objects):
            return intersection
        dist = self.objects[depth].intersects(ray, depth)
        if dist is not None and \
            (intersection is None or dist < intersection[1]):
            intersection = self.objects[depth], dist
        return intersection

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111,aspect='equal')
        defautl_wvl = 0.
        legend_lbl = []
        legend_line = []
        for ray in self.rays:
            if np.abs(ray.direction[0].x) > 1.e-10:
                continue
            line = ray.plot(ax=ax)
            if ray.wavelength != defautl_wvl:
                defautl_wvl = ray.wavelength
                line.set_label('%d nm'%ray.wavelength)
            ax.legend()
        objects_lines = []
        for obj in self.objects:
            line = obj.plot(ax=ax)
            objects_lines.append(line)
        for jobj,obj in enumerate(self.objects):
            if obj.material.refractiveIndex.coefficients != \
                       air.refractiveIndex.coefficients:
                if jobj < len(self.objects)-1:
                    index1 = np.argsort(objects_lines[jobj].get_ydata())
                    index2 = np.argsort(objects_lines[jobj+1].get_ydata())[::-1]
                    x = np.append(objects_lines[jobj].get_xdata()[index1],\
                                  objects_lines[jobj+1].get_xdata()[index2])
                    y = np.append(objects_lines[jobj].get_ydata()[index1],\
                                  objects_lines[jobj+1].get_ydata()[index2])
                else:
                    x = np.append(objects_lines[jobj].get_xdata(),[objects_lines[jobj].get_xdata()+500,\
                                                        objects_lines[jobj].get_xdata()+500])
                    y = objects_lines[jobj].get_ydata()
                    y = np.append(y,[y.max(),y.min()])
                verts = list(zip(x, y))
                poly = Polygon(verts,facecolor=(0,0,1,0.2),edgecolor=(0,0,1,1),ls='dashed')
                ax.add_patch(poly)
        plt.show()

