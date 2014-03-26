#!/usr/bin/env python
# -*- coding: utf-8 -*-


from itertools import product

from numpy import sin, cos, pi


def create_cyl_array(radius, length, n_segments, centered=True):
    vertices, indices, normals = [], [], []
    z = -length / 2. if centered else 0
    # Side of the cylinder
    for i in range(n_segments):
        alpha = 2 * pi * i / n_segments
        xn, yn = cos(alpha), sin(alpha)
        x, y = radius*xn, radius*yn
        vertices += (x, y, z, x, y, z + length)
        normals += (xn, yn, 0, xn, yn, 0)
    num_trian = 2 * n_segments
    for i in range(num_trian):
        indices += (i, (i + 1) % num_trian, (i + 2) % num_trian)
    # Disks of the cylinder
    for j in range(2):
        for i in range(n_segments):
            alpha = 2 * pi * i / n_segments
            x, y = radius*cos(alpha), radius*sin(alpha)
            vertices += (x, y, z)
            normals += (0, 0, -1 + j*2)
            indices += (j*n_segments + num_trian + i,
                        j*n_segments + num_trian + (i + 1) % n_segments,
                        j*n_segments + n_segments + num_trian)
        vertices += (0, 0, z)
        normals += (0, 0, -1 + j*2)
        z += length
    return vertices, indices, normals


def create_sphere_array(radius, n_lat, n_long):
    vertices, indices = [0., 0., radius, 0., 0., -radius], []
    normals = [0., 0., 1., 0., 0., -1.]
    for i in range(1, n_lat):
        theta = i*pi/n_lat
        for j in range(n_long):
            phi = j*2.*pi/n_long
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
            indices += (j + start_i, j + start_i + n_long,
                        (j + 1) % n_long + start_i + n_long)
            indices += (j + start_i, (j + 1) % n_long + start_i,
                        (j + 1) % n_long + start_i + n_long)
    return vertices, indices, normals


def create_arrow_array(length, thick_rod=0.04, thick_hat=0.125, length_hat=1.3):
    vertices, indices, normals = [], [], []
    # Arrow rod
    for i, j in product([-1, 1], repeat=2):
        vertices += (i*thick_rod*length, i*j*thick_rod*length, 0)
        vertices += (i*thick_rod*length, i*j*thick_rod*length, length)
        normals += 2 * (0.707*i, 0.707*i*j, 0)
    # Arrow head
    for i, j in product([-1, 1], repeat=2):
        vertices += (i*thick_hat*length, i*j*thick_hat*length, length)
        normals += (0.707*i, 0.707*i*j, 0)
    vertices += (0, 0, length_hat * length)
    normals += (0, 0, 1)
    for i in range(8):
        indices += (i, (i + 1) % 8, (i + 2) % 8)
    for i in range(4):
        indices += (i + 8, (i + 1) % 4 + 8, 12)
    return vertices, indices, normals


def create_box_array(length=3., width=1.):
    dim = [0.5 * width, 0.5 * width, 0.5 * length]
    vertices, normals = [], []
    for l in range(3):
        for i, j, k in product([-1, 1], repeat=3):
            coord = [0, 0, 0]
            coord[l] = i
            normals += coord
            coord[l] = i*dim[l]
            coord[(l+1) % 3] = j*dim[(l+1) % 3]
            coord[(l+2) % 3] = j*k*dim[(l+2) % 3]
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


