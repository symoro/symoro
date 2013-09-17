__author__ = 'Izzat'

from itertools import product
from numpy import sin, cos, pi, array
from OpenGL.GL import *

def create_cyl_array(radius, length, n_segments, centered=True):
    vertices, indices, normals = [], [], []
    z = -length / 2. if centered else 0
    # Side of the cylinder
    for i in range(n_segments):
        alpha = 2 * pi * i / n_segments
        xn, yn = cos(alpha), sin(alpha)
        x, y = radius*xn, radius*yn
        vertices += ( x, y, z, x, y, z + length )
        normals += (xn, yn, 0, xn, yn, 0)
    num_trian = 2 * n_segments
    for i in range(num_trian):
        indices += (i, (i + 1) % num_trian, (i + 2) % num_trian)
    # Disks of the cylinder
    for j in range(2):
        for i in range(n_segments):
            alpha = 2 * pi * i / n_segments
            x, y = radius*cos(alpha), radius*sin(alpha)
            vertices += ( x, y, z )
            normals += (0, 0, -1 + j*2)
            indices += (j*n_segments + num_trian + i, j*n_segments + num_trian + (i + 1) % n_segments, j*n_segments + n_segments + num_trian)
        vertices += (0, 0, z)
        normals += (0, 0, -1 + j*2)
        z += length
    return vertices, indices, normals

def create_sphere_array(radius, n_lat, n_long):
    vertices, indices, normals = [0., 0., radius, 0., 0., -radius], [], [0.,0.,1.,0.,0.,-1.]
    for i in range(1, n_lat):
        theta = i*pi/n_lat
        for j in range(n_long):
            phi = j*2*pi/n_long
            xn, yn, zn = cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)
            x, y, z = radius * xn, radius * yn, radius * zn
            normals += xn, yn, zn
            vertices += (x, y, z)
    # top and bottom indices
    for i, start in enumerate([2, len(vertices)/3 - n_long]):
        for j in range(n_long):
            indices += (i, j + start, (j + 1) % n_long + start)
    # in-between
    for i in range(n_lat - 2):
        start_i = 2 + i * n_long
        for j in range(n_long):
            indices += (j + start_i, j + start_i + n_long, (j + 1) % n_long + start_i + n_long)
            indices += (j + start_i, (j + 1) % n_long + start_i, (j + 1) % n_long + start_i + n_long)
    return vertices, indices, normals

def create_arrow_array(length, thick_rod_scale=0.04, thick_hat_scale=0.125, length_hat_scale=1.3):
    vertices, indices, normals = [], [], []
    # Arrow rod
    for i, j in product([-1, 1], repeat=2):
        vertices += (i*thick_rod_scale*length,i*j*thick_rod_scale*length,0)
        vertices += (i*thick_rod_scale*length,i*j*thick_rod_scale*length,length)
        normals += 2 * (0.707*i, 0.707*i*j, 0)
    # Arrow head
    for i, j in product([-1, 1], repeat=2):
        vertices += (i*thick_hat_scale*length,i*j*thick_hat_scale*length,length)
        normals += (0.707*i, 0.707*i*j, 0)
    vertices += (0, 0, length_hat_scale * length)
    normals += (0, 0, 1)
    for i in range(8):
        indices += (i, (i + 1) % 8, (i + 2) % 8)
    for i in range(4):
        indices += (i + 8, (i + 1) % 4 + 8, 12)
    return vertices, indices, normals

def create_box_array(length=3., width=1.):
    dim = [width/2., width/2., length/2.]
    vertices, normals = [], []
    for l in range(3):
        for i, j, k in product([-1, 1], repeat=3):
            coord = [0, 0, 0]
            coord[l] = i
            normals += coord
            coord[l] = i*dim[l]
            coord[(l+1)%3] = j*dim[(l+1)%3]
            coord[(l+2)%3] = j*k*dim[(l+2)%3]
            vertices += coord
    return vertices, normals

class Primitives:

    @classmethod
    def box_array(cls, length):
        """ Returns:
            vertices, normals
        """
        return create_box_array(length, length/3.)

    @classmethod
    def cyl_array(cls, length):
        """ Returns:
            vertices, indices, normals
        """
        return create_cyl_array(length/6., length, 8)

    @classmethod
    def rod_array(cls, length):
        """ Returns:
            vertices, indices, normals
        """
        return create_cyl_array(length/36., 1., 8, centered=False)

    @classmethod
    def sph_array(cls, length):
        """ Returns:
            vertices, indices, normals
        """
        return create_sphere_array(length/5., 8, 8)

    @classmethod
    def arr_array(cls, length):
        """ Returns:
            vertices, indices, normals
        """
        return create_arrow_array(length)

class PrimDrawer:
    def __init__(self, scale=1.):
        self.length = 3.0
        self.radius = 0.5
        self.box_vertices, self.box_normals = create_box_array(self.length, self.length/3.)
        self.cyl_vertices, self.cyl_indices, self.cyl_normals = create_cyl_array(self.radius, self.length, 16)
        self.rod_vertices, self.rod_indices, self.rod_normals = create_cyl_array(self.radius/6., 1., 8, centered=False)
        self.sph_vertices, self.sph_indices, self.sph_normals = create_sphere_array(self.length/6., 16, 16)
        self.arr_vertices, self.arr_indices, self.arr_normals = create_arrow_array(self.length)
        # self.box_vertices, self.box_normals = create_box_array(3 * scale, scale)
        # self.cyl_vertices, self.cyl_indices, self.cyl_normals = create_cyl_array(scale / 2., 3 * scale, 16)
        # self.rod_vertices, self.rod_indices, self.rod_normals = create_cyl_array(scale / 16., 1., 8, centered=False)
        # self.sph_vertices, self.sph_indices, self.sph_normals = create_sphere_array(scale / 2., 16, 16)
        # self.arr_vertices, self.arr_indices, self.arr_normals = create_arrow_array(3 * scale)

    def Box(self):
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glColor3f(1.0, 0.6, 0.)
        glVertexPointer(3, GL_FLOAT, 0, self.box_vertices)
        glNormalPointer(GL_FLOAT, 0, self.box_normals)
        glDrawArrays(GL_QUADS,0,24)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)

    def Cylinder(self):
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glColor3f(1., 1., 0.)
        glVertexPointer(3, GL_FLOAT, 0, self.cyl_vertices)
        glNormalPointer(GL_FLOAT, 0, self.cyl_normals)
        glDrawElements(GL_TRIANGLES, len(self.cyl_indices), GL_UNSIGNED_INT, self.cyl_indices)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)

    def Rod(self, length):
        glPushMatrix()
        glMultMatrixf(array([[1.,0.,0.,0.], [0.,1.,0.,0.], [0.,0.,length,0.],[0.,0.,0.,1.]]))
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glColor3f(0.8, 0.51, 0.25)
        glVertexPointer(3, GL_FLOAT, 0, self.rod_vertices)
        glNormalPointer(GL_FLOAT, 0, self.rod_normals)
        glDrawElements(GL_TRIANGLES, len(self.rod_indices), GL_UNSIGNED_INT, self.rod_indices)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)
        glPopMatrix()

    def Arrow(self):
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glVertexPointer(3, GL_FLOAT, 0, self.arr_vertices)
        glNormalPointer(GL_FLOAT, 0, self.arr_normals)
        glDrawElements(GL_TRIANGLES, len(self.arr_indices), GL_UNSIGNED_INT, self.arr_indices)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)

    def Sphere(self):
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glColor3f(1.0, 0., 1.)
        glVertexPointer(3, GL_FLOAT, 0, self.sph_vertices)
        glNormalPointer(GL_FLOAT, 0, self.sph_normals)
        glDrawElements(GL_TRIANGLES, len(self.sph_indices), GL_UNSIGNED_INT, self.sph_indices)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)
