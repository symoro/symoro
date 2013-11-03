#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Izzat'

import OpenGL.GL as gl
from numpy import degrees, identity, array

from primitives import Primitives


class Frame(object):
    def __init__(self, index=0, length=0.5, T=identity(4), show_frame=True):
        self.children = []
        self.T = T
        self.show_frame = show_frame
        self.index = index
        self.length = length
        self.arr_vertices, self.arr_indices, self.arr_normals = \
            Primitives.arr_array(self.length)
        # TODO: Change to variable

    def draw_frame(self):
        if self.show_frame:
            gl.glPushMatrix()
            gl.glColor3f(1, 0, 0)
            self.draw_arrow()
            gl.glRotatef(90, 0, 1, 0)
            gl.glColor3f(0, 1, 0)
            self.draw_arrow()
            gl.glPopMatrix()

    def draw_arrow(self):
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.arr_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.arr_normals)
        gl.glDrawElements(gl.GL_TRIANGLES, len(self.arr_indices),
                          gl.GL_UNSIGNED_INT, self.arr_indices)

    def set_show_frame(self, show=True):
        self.show_frame = show

    def add_child(self, child):
        self.children.append(child)

    def draw_frames(self):
        gl.glPushMatrix()
        gl.glMultMatrixf(self.T)
        self.draw_frame()
        for child in self.children:
            child.draw_frames()
        gl.glPopMatrix()

    def draw(self):
        gl.glPushMatrix()
        gl.glMultMatrixf(self.T)
        for child in self.children:
            child.draw()
        gl.glPopMatrix()

    def __str__(self):
        return '{0} Children: {1}'.format(self.index, self.children)


class JointObject(Frame):

    def __init__(self, index, length,
                 theta=0., r=0., alpha=0., d=0., gamma=0., b=0.):
        super(JointObject, self).__init__(index, length)
        self.theta = theta
        self.r = r
        self.alpha = alpha
        self.d = d
        self.gamma = gamma
        self.b = b
        self.shift = 0.
        self.length = length
        self.rod_vertices, self.rod_indices, self.rod_normals = \
            Primitives.rod_array(self.length)

    def draw_rod(self, length):
        gl.glPushMatrix()
        gl.glMultMatrixf(array([[1., 0., 0., 0.], [0., 1., 0., 0.],
                                [0., 0., length, 0.], [0., 0., 0., 1.]]))
        gl.glColor3f(0.8, 0.51, 0.25)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.rod_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.rod_normals)
        gl.glDrawElements(gl.GL_TRIANGLES, len(self.rod_indices),
                          gl.GL_UNSIGNED_INT, self.rod_indices)
        gl.glPopMatrix()

    def set_shift(self, shift):
        self.shift = shift

    def draw_frames(self):
        gl.glPushMatrix()
        gl.glRotatef(degrees(self.gamma), 0, 0, 1)
        gl.glTranslatef(0, 0, self.b)
        gl.glRotatef(degrees(self.alpha), 1, 0, 0)
        gl.glTranslatef(self.d, 0, 0)
        gl.glRotatef(degrees(self.theta), 0, 0, 1)
        gl.glTranslatef(0, 0, self.r)
        self.draw_frame()
        for child in self.children:
            child.draw_frames()
        gl.glPopMatrix()

    def draw(self):
        gl.glPushMatrix()
        if self.b:
            self.draw_rod(self.b)
            gl.glTranslatef(0, 0, self.b)
        gl.glRotatef(degrees(self.gamma), 0, 0, 1)
        if self.d:
            gl.glPushMatrix()
            gl.glRotatef(90, 0, 1, 0)
            self.draw_rod(self.d)
            gl.glPopMatrix()
            gl.glTranslatef(self.d, 0, 0)
        gl.glRotatef(degrees(self.alpha), 1, 0, 0)
        if self.r:
            self.draw_rod(self.r)
            gl.glTranslatef(0, 0, self.r)
        gl.glRotatef(degrees(self.theta), 0, 0, 1)
        if self.shift:
            gl.glPushMatrix()
            self.draw_rod(self.shift)
            gl.glTranslatef(0, 0, self.shift)
            self.draw_joint()
            gl.glPopMatrix()
        else:
            self.draw_joint()
        for child in self.children:
            child.draw()
        gl.glPopMatrix()


class RevoluteJoint(JointObject):

    def __init__(self, *args):
        super(RevoluteJoint, self).__init__(*args)
        self.cyl_vertices, self.cyl_indices, self.cyl_normals = \
            Primitives.cyl_array(self.length)
        self.q_init = self.theta

    def draw_joint(self):
        gl.glColor3f(1., 1., 0.)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.cyl_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.cyl_normals)
        gl.glDrawElements(gl.GL_TRIANGLES, len(self.cyl_indices),
                          gl.GL_UNSIGNED_INT, self.cyl_indices)

    @property
    def q(self):
        return self.theta

    @q.setter
    def q(self, theta):
        self.theta = theta


class PrismaticJoint(JointObject):

    def __init__(self, *args):
        super(PrismaticJoint, self).__init__(*args)
        self.box_vertices, self.box_normals = Primitives.box_array(self.length)
        self.q_init = self.r

    def draw_joint(self):
        gl.glColor3f(1., 0.6, 0.)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.box_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.box_normals)
        gl.glDrawArrays(gl.GL_QUADS, 0, 24)

    @property
    def q(self):
        return self.r

    @q.setter
    def q(self, r):
        self.r = r


class FixedJoint(JointObject):

    def __init__(self, *args):
        super(FixedJoint, self).__init__(*args)
        self.sph_vertices, self.sph_indices, self.sph_normals = \
            Primitives.sph_array(self.length)

    def draw_joint(self):
        gl.glColor3f(1., 0., 1.)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.sph_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.sph_normals)
        gl.glDrawElements(gl.GL_TRIANGLES, len(self.sph_indices),
                          gl.GL_UNSIGNED_INT, self.sph_indices)
