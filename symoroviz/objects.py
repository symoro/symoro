# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


import OpenGL.GL as gl

from numpy import degrees, identity, array

from primitives import Primitives


class Frame(object):
    def __init__(self, index=0, T=identity(4), show_frame=True):
        self.children = []
        self.T = T
        self.show_frame = show_frame
        self.index = index

    def __str__(self):
        return '(Frame {0} Children: {1})'.format(self.index, self.children)

    def __repr__(self):
        return self.__str__()

    def draw_frame(self):
        if self.show_frame:
            gl.glPushMatrix()
            # z-axis (joint axis) - blue
            gl.glColor3f(0, 0, 1)
            self.draw_arrow()
            # x-axis - red
            gl.glRotatef(90, 0, 1, 0)
            gl.glColor3f(1, 0, 0)
            self.draw_arrow()
            gl.glPopMatrix()

    def draw_arrow(self):
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.arr_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.arr_normals)
        gl.glDrawElements(
            gl.GL_TRIANGLES, len(self.arr_indices),
            gl.GL_UNSIGNED_INT, self.arr_indices
        )

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
        self.has_base = False
        self.has_end = False

    def draw_rod(self, length):
        gl.glPushMatrix()
        gl.glMultMatrixf(array([
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., length, 0.],
            [0., 0., 0., 1.]
        ]))
        gl.glColor3f(0.8, 0.51, 0.25)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.rod_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.rod_normals)
        gl.glDrawElements(
            gl.GL_TRIANGLES, len(self.rod_indices),
            gl.GL_UNSIGNED_INT, self.rod_indices
        )
        gl.glPopMatrix()

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
            shift = self.shift*self.length
            self.draw_rod(shift)
            gl.glTranslatef(0, 0, shift)
            self.draw_joint()
            self.draw_base()
            self.draw_end()
            gl.glPopMatrix()
        else:
            self.draw_joint()
            self.draw_base()
            self.draw_end()
        for child in self.children:
            child.draw()
        gl.glPopMatrix()

    def draw_base(self):
        if self.has_base:
            gl.glColor3f(0.0, 0.0, 0.0)
            gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.sph_vertices)
            gl.glNormalPointer(gl.GL_FLOAT, 0, self.sph_normals)
            gl.glDrawElements(
                gl.GL_TRIANGLES, len(self.sph_indices),
                gl.GL_UNSIGNED_INT, self.sph_indices,
            )

    def draw_end(self):
        if self.has_end:
            gl.glColor3f(0.2, 0.7, 0.0)
            gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.sph_vertices)
            gl.glNormalPointer(gl.GL_FLOAT, 0, self.sph_normals)
            gl.glDrawElements(
                gl.GL_TRIANGLES, len(self.sph_indices),
                gl.GL_UNSIGNED_INT, self.sph_indices,
            )

    def set_length(self, new_length):
        if not self.init_length:
            self.rod_vertices, self.rod_indices, self.rod_normals = \
                Primitives.rod_array(new_length)
            self.init_length = new_length
        self.length = new_length
        if self.has_base or self.has_end:
            self.sph_vertices, self.sph_indices, self.sph_normals = \
                Primitives.sph_array(1.5 * new_length)
        super(JointObject, self).set_length(new_length)


class RevoluteJoint(JointObject):
    def __init__(self, *args):
        super(RevoluteJoint, self).__init__(*args)
        self.q_init = self.theta

    def __str__(self):
        return '(RevoluteJoint {0} children: {1})'.format(
            self.index, self.children
        )

    def __repr__(self):
        return self.__str__()

    def draw_joint(self):
        gl.glColor3f(1., 1., 0.)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.cyl_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.cyl_normals)
        gl.glDrawElements(
            gl.GL_TRIANGLES, len(self.cyl_indices),
            gl.GL_UNSIGNED_INT, self.cyl_indices
        )

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

    def __str__(self):
        return '(PrismaticJoint {0} children: {1})'.format(
            self.index, self.children
        )

    def __repr__(self):
        return self.__str__()

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

    def set_length(self, new_length):
        self.box_vertices, self.box_normals = \
            Primitives.box_array(new_length)
        super(PrismaticJoint, self).set_length(new_length)

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
        if self.shift:
            gl.glPushMatrix()
            shift = self.shift*self.length
            self.draw_rod(shift)
            gl.glTranslatef(0, 0, shift)
            self.draw_joint()
            self.draw_base()
            self.draw_end()
            gl.glPopMatrix()
        else:
            self.draw_joint()
            self.draw_base()
            self.draw_end()
        if self.r:
            self.draw_rod(self.r)
            gl.glTranslatef(0, 0, self.r)
        gl.glRotatef(degrees(self.theta), 0, 0, 1)
        for child in self.children:
            child.draw()
        gl.glPopMatrix()


class FixedJoint(JointObject):
    def __init__(self, *args):
        super(FixedJoint, self).__init__(*args)

    def __str__(self):
        return '(FixedJoint {0} children: {1})'.format(
            self.index, self.children
        )

    def __repr__(self):
        return self.__str__()

    def draw_joint(self):
        gl.glColor3f(1., 0., 1.)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.sph_vertices)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.sph_normals)
        gl.glDrawElements(
            gl.GL_TRIANGLES, len(self.sph_indices),
            gl.GL_UNSIGNED_INT, self.sph_indices
        )

    def set_length(self, new_length):
        self.sph_vertices, self.sph_indices, self.sph_normals = \
            Primitives.sph_array(new_length)
        super(FixedJoint, self).set_length(new_length)
