#!/usr/bin/env python
# -*- coding: utf-8 -*-


import OpenGL.GL as gL

from numpy import degrees, identity, array

from primitives import Primitives


class Frame(object):
    def __init__(self, index=0, T=identity(4), show_frame=True):
        self.children = []
        self.T = T
        self.show_frame = show_frame
        self.index = index

    def draw_frame(self):
        if self.show_frame:
            gL.glPushMatrix()
            gL.glColor3f(1, 0, 0)
            self.draw_arrow()
            gL.glRotatef(90, 0, 1, 0)
            gL.glColor3f(0, 1, 0)
            self.draw_arrow()
            gL.glPopMatrix()

    def draw_arrow(self):
        gL.glVertexPointer(3, gL.GL_FLOAT, 0, self.arr_vertices)
        gL.glNormalPointer(gL.GL_FLOAT, 0, self.arr_normals)
        gL.glDrawElements(gL.GL_TRIANGLES, len(self.arr_indices),
                          gL.GL_UNSIGNED_INT, self.arr_indices)

    def set_show_frame(self, show=True):
        self.show_frame = show

    def add_child(self, child):
        self.children.append(child)

    def draw_frames(self):
        gL.glPushMatrix()
        gL.glMultMatrixf(self.T)
        self.draw_frame()
        for child in self.children:
            child.draw_frames()
        gL.glPopMatrix()

    def draw(self):
        gL.glPushMatrix()
        gL.glMultMatrixf(self.T)
        for child in self.children:
            child.draw()
        gL.glPopMatrix()

    def __str__(self):
        return '{0} Children: {1}'.format(self.index, self.children)

    def set_length(self, new_length):
        self.arr_vertices, self.arr_indices, self.arr_normals = \
            Primitives.arr_array(new_length)


class JointObject(Frame):

    def __init__(self, index, theta=0., r=0., alpha=0., d=0., gamma=0., b=0.):
        super(JointObject, self).__init__(index)
        self.theta = theta
        self.r = r
        self.alpha = alpha
        self.d = d
        self.gamma = gamma
        self.b = b
        self.shift = 0.
        self.init_length = 0.

    def draw_rod(self, length):
        gL.glPushMatrix()
        gL.glMultMatrixf(array([[1., 0., 0., 0.], [0., 1., 0., 0.],
                                [0., 0., length, 0.], [0., 0., 0., 1.]]))
        gL.glColor3f(0.8, 0.51, 0.25)
        gL.glVertexPointer(3, gL.GL_FLOAT, 0, self.rod_vertices)
        gL.glNormalPointer(gL.GL_FLOAT, 0, self.rod_normals)
        gL.glDrawElements(gL.GL_TRIANGLES, len(self.rod_indices),
                          gL.GL_UNSIGNED_INT, self.rod_indices)
        gL.glPopMatrix()

    def draw_frames(self):
        gL.glPushMatrix()
        gL.glRotatef(degrees(self.gamma), 0, 0, 1)
        gL.glTranslatef(0, 0, self.b)
        gL.glRotatef(degrees(self.alpha), 1, 0, 0)
        gL.glTranslatef(self.d, 0, 0)
        gL.glRotatef(degrees(self.theta), 0, 0, 1)
        gL.glTranslatef(0, 0, self.r)
        self.draw_frame()
        for child in self.children:
            child.draw_frames()
        gL.glPopMatrix()

    def draw(self):
        gL.glPushMatrix()
        if self.b:
            self.draw_rod(self.b)
            gL.glTranslatef(0, 0, self.b)
        gL.glRotatef(degrees(self.gamma), 0, 0, 1)
        if self.d:
            gL.glPushMatrix()
            gL.glRotatef(90, 0, 1, 0)
            self.draw_rod(self.d)
            gL.glPopMatrix()
            gL.glTranslatef(self.d, 0, 0)
        gL.glRotatef(degrees(self.alpha), 1, 0, 0)
        if self.r:
            self.draw_rod(self.r)
            gL.glTranslatef(0, 0, self.r)
        gL.glRotatef(degrees(self.theta), 0, 0, 1)
        if self.shift:
            gL.glPushMatrix()
            shift = self.shift*self.length
            self.draw_rod(shift)
            gL.glTranslatef(0, 0, shift)
            self.draw_joint()
            gL.glPopMatrix()
        else:
            self.draw_joint()
        for child in self.children:
            child.draw()
        gL.glPopMatrix()

    def set_length(self, new_length):
        if not self.init_length:
            self.rod_vertices, self.rod_indices, self.rod_normals = \
                Primitives.rod_array(new_length)
            self.init_length = new_length
        self.length = new_length
        super(JointObject, self).set_length(new_length)


class RevoluteJoint(JointObject):

    def __init__(self, *args):
        super(RevoluteJoint, self).__init__(*args)
        self.q_init = self.theta

    def draw_joint(self):
        gL.glColor3f(1., 1., 0.)
        gL.glVertexPointer(3, gL.GL_FLOAT, 0, self.cyl_vertices)
        gL.glNormalPointer(gL.GL_FLOAT, 0, self.cyl_normals)
        gL.glDrawElements(gL.GL_TRIANGLES, len(self.cyl_indices),
                          gL.GL_UNSIGNED_INT, self.cyl_indices)

    @property
    def q(self):
        return self.theta

    @q.setter
    def q(self, theta):
        self.theta = theta

    def set_length(self, new_length):
        self.cyl_vertices, self.cyl_indices, self.cyl_normals = \
            Primitives.cyl_array(new_length)
        super(RevoluteJoint, self).set_length(new_length)


class PrismaticJoint(JointObject):

    def __init__(self, *args):
        super(PrismaticJoint, self).__init__(*args)
        self.q_init = self.r

    def draw_joint(self):
        gL.glColor3f(1., 0.6, 0.)
        gL.glVertexPointer(3, gL.GL_FLOAT, 0, self.box_vertices)
        gL.glNormalPointer(gL.GL_FLOAT, 0, self.box_normals)
        gL.glDrawArrays(gL.GL_QUADS, 0, 24)

    @property
    def q(self):
        return self.r

    @q.setter
    def q(self, r):
        self.r = r

    def set_length(self, new_length):
        self.box_vertices, self.box_normals = Primitives.box_array(new_length)
        super(PrismaticJoint, self).set_length(new_length)

    def draw(self):
        gL.glPushMatrix()
        if self.b:
            self.draw_rod(self.b)
            gL.glTranslatef(0, 0, self.b)
        gL.glRotatef(degrees(self.gamma), 0, 0, 1)
        if self.d:
            gL.glPushMatrix()
            gL.glRotatef(90, 0, 1, 0)
            self.draw_rod(self.d)
            gL.glPopMatrix()
            gL.glTranslatef(self.d, 0, 0)
        gL.glRotatef(degrees(self.alpha), 1, 0, 0)
        if self.shift:
            gL.glPushMatrix()
            shift = self.shift*self.length
            self.draw_rod(shift)
            gL.glTranslatef(0, 0, shift)
            self.draw_joint()
            gL.glPopMatrix()
        else:
            self.draw_joint()
        if self.r:
            self.draw_rod(self.r)
            gL.glTranslatef(0, 0, self.r)
        gL.glRotatef(degrees(self.theta), 0, 0, 1)
        for child in self.children:
            child.draw()
        gL.glPopMatrix()


class FixedJoint(JointObject):

    def __init__(self, *args):
        super(FixedJoint, self).__init__(*args)

    def draw_joint(self):
        gL.glColor3f(1., 0., 1.)
        gL.glVertexPointer(3, gL.GL_FLOAT, 0, self.sph_vertices)
        gL.glNormalPointer(gL.GL_FLOAT, 0, self.sph_normals)
        gL.glDrawElements(gL.GL_TRIANGLES, len(self.sph_indices),
                          gL.GL_UNSIGNED_INT, self.sph_indices)

    def set_length(self, new_length):
        self.sph_vertices, self.sph_indices, self.sph_normals = \
            Primitives.sph_array(new_length)
        super(FixedJoint, self).set_length(new_length)


