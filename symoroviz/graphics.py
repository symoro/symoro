#!/usr/bin/env python
# -*- coding: utf-8 -*-


from math import atan2

import wx
import wx.lib.agw.floatspin as FS
from wx.glcanvas import GLCanvas

import OpenGL.GL as gl
import OpenGL.GLU as glu

from numpy import sin, cos, radians, pi, inf, nan

from sympy import Expr

from pysymoro.invgeom import loop_solve
from pysymoro.geometry import dgm
from symoroutils import samplerobots
from symoroutils import symbolmgr
from symoroutils import tools

from objects import Frame, RevoluteJoint, FixedJoint, PrismaticJoint

#TODO: Fullscreen camera rotation bug
#TODO: X-, Z-axis
#TODO: Random button

class myGLCanvas(GLCanvas):
    def __init__(self, parent, robo, params, size=(600, 600)):
        super(myGLCanvas, self).__init__(parent, size=size)
        self.Bind(wx.EVT_PAINT, self.OnPaintAll)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_RIGHT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_RIGHT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.size = self.GetClientSize()
        self.hor_angle = self.ver_angle = 0
        self.cen_x = self.cen_y = self.cen_z = 0
        self.robo = robo
        self.q_sym = self.robo.q_vec
        self.q_pas_sym = self.robo.q_passive
        self.q_act_sym = self.robo.q_active
        self.pars_num = params
        self.init = 0
        self.distance = 5.
        self.fov = 40.
        self.jnt_objs = []
        self.construct_hierarchy()
        self.dgms = {}
        self.l_solver = None

    def OnEraseBackground(self, event):
        # Do nothing, to avoid flashing on MSW.
        # This is needed
        pass

    def assign_mono_scale(self):
        """Sets the coefficients used to draw objects

        This function calculates coefficients which are used
        to draw the objects (Joints, links, end-effectors)
        It computes the minimum and maximum r or d different from 0.
        Then uses those sizes to determine the reference
        numbers which are used all over the class.
        """
        minv = inf
        for jnt in self.jnt_objs[1:]:
            dist = max(abs(jnt.r), abs(jnt.d))
            if dist < minv and dist != 0:
                minv = dist
        if minv == inf:
            minv = 1.
        self.length = 0.4 * minv
        #print self.length
        for jnt in self.jnt_objs:
            if isinstance(jnt, PrismaticJoint):
                jnt.r = 3.5 * self.length
            jnt.set_length(self.length)

    def add_items_to_frame(self, frame, index, jnt_hier):
        children = jnt_hier[index]
        for child_i in children:
            params = [child_i]
            for par in ['theta', 'r', 'alpha', 'd', 'gamma', 'b']:
                val = self.robo.get_val(child_i, par)
                try:
                    if isinstance(val, Expr):
                        params.append(float(val.subs(self.pars_num)))
                    else:
                        params.append(float(val))
                except:
                    if val in self.q_sym:
                        params.append(0.)

            if self.robo.sigma[child_i] == 0:
                child_frame = RevoluteJoint(*params)
            elif self.robo.sigma[child_i] == 1:
                child_frame = PrismaticJoint(*params)
            else:
                child_frame = FixedJoint(*params)
            self.jnt_objs.append(child_frame)
            frame.add_child(child_frame)
            self.add_items_to_frame(child_frame, child_i, jnt_hier)

    def get_joints_dictionary(self):
        result = {}
        for jnt_i in range(self.robo.NF):
            result[jnt_i] = [i for i, child in enumerate(self.robo.ant)
                             if child == jnt_i]
        return result

    def construct_hierarchy(self):
        jnt_hier = self.get_joints_dictionary()
        self.base = Frame(0)
        self.jnt_objs.append(self.base)
        self.add_items_to_frame(self.base, 0, jnt_hier)
        self.jnt_objs.sort(key=lambda jnt_obj: jnt_obj.index)
        self.jnt_dict = {}
        for i, jnt in enumerate(self.jnt_objs):
            sym = self.robo.get_q(i)
            if sym != 0:
                self.jnt_dict[sym] = jnt
        self.assign_mono_scale()
        self.representation(expanded=True)

    def OnSize(self, event):
        size = self.size = self.GetClientSize()
        if self.GetContext():
            self.SetCurrent()
            gl.glViewport(0, 0, size.width, size.height)
            #gl.glMatrixMode(gl.GL_PROJECTION)
            #gl.glLoadIdentity()
            #gl.gluPerspective(40.0, size.width/size.height, 1.0, 100.0)
            #gl.glMatrixMode(gl.GL_MODELVIEW)
        event.Skip()

    def OnMouseDown(self, evt):
        self.CaptureMouse()
        self.lastx, self.lasty = evt.GetPosition()

    def OnMouseUp(self, _):
        if self.HasCapture():
            self.ReleaseMouse()

    def OnMouseMotion(self, evt):
        if evt.Dragging():
            x, y = evt.GetPosition()
            if evt.LeftIsDown():
                dx, dy = x - self.lastx, y - self.lasty
                self.lastx, self.lasty = x, y
                coef = self.distance/self.length*sin(radians(self.fov/2.))
                self.hor_angle += dx*coef/self.size.width
                self.ver_angle += dy*coef/self.size.height
                self.ver_angle = max(min(pi/2, self.ver_angle), -pi/2)
            elif evt.RightIsDown():
                dy = y - self.lasty
                self.lasty = y
                self.distance *= 1 + 2 * float(dy) / self.size.height
            self.CameraTransformation()
            self.Refresh(False)

    def dgm_for_frame(self, i):
        if i not in self.dgms:
            jnt = self.jnt_objs[i]
            if i > 0 and jnt.r == 0 and jnt.d == 0 and jnt.b == 0:
                self.dgms[i] = self.dgm_for_frame(self.robo.ant[i])
            else:
                symo = symbolmgr.SymbolManager(sydi=self.pars_num)
                T = dgm(self.robo, symo, 0, i, fast_form=True, trig_subs=True)
                self.dgms[i] = symo.gen_func('dgm_generated', T, self.q_sym)
        return self.dgms[i]

    def find_solution(self, qs_act, qs_pas):
        slns = self.l_solver(qs_act)
        min_error = 99999999
        min_sln = None
        for sln in slns:
            if nan in sln:
                continue
            error = 0
            for j, (qp, sigma) in enumerate(qs_pas):
                if sigma == 0:
                    error += atan2(sin(qp-sln[j]), cos(qp-sln[j]))**2
                else:
                    error += (qp - sln[j])**2
            if error < min_error:
                min_sln = sln
                min_error = error
        if min_sln is not None:
            for i, sym in enumerate(self.q_pas_sym):
                self.jnt_dict[sym].q = min_sln[i]

    def solve(self):
        if self.robo.structure != tools.CLOSED_LOOP:
            return
        if self.l_solver is None:
            self.generate_loop_fcn()
        qs_act = []
        qs_pas = []
        for sym in self.q_sym:
            if sym in self.q_act_sym:
                qs_act.append(self.jnt_dict[sym].q)
            elif sym in self.q_pas_sym:
                i = self.jnt_dict[sym].index
                if i < self.robo.NL:
                    qs_pas.append((self.jnt_dict[sym].q, self.robo.sigma[i]))

        self.find_solution(qs_act, qs_pas)

    def generate_loop_fcn(self):
        symo = symbolmgr.SymbolManager(sydi=self.pars_num)
        loop_solve(self.robo, symo)
        self.l_solver = symo.gen_func('IGM_gen', self.q_pas_sym,
                                      self.q_act_sym, multival=True)

    def centralize_to_frame(self, index):
        q_vec = [self.jnt_dict[sym].q for sym in self.q_sym]
        T = self.dgm_for_frame(index)(q_vec)
        self.cen_x, self.cen_y, self.cen_z = T[0][3], T[1][3], T[2][3]
        self.CameraTransformation()
        self.Refresh(False)

    def OnPaintAll(self, event):
        self.SetCurrent()
        if not self.init:
            self.InitGL()
            self.init = 1
        self.OnDraw()
        event.Skip()

    def CameraTransformation(self):
        gl.glLoadIdentity()
        glu.gluLookAt(
            self.cen_x-self.distance*cos(self.ver_angle)*sin(self.hor_angle),
            self.cen_y-self.distance*cos(self.ver_angle)*cos(self.hor_angle),
            self.cen_z+self.distance*sin(self.ver_angle),
            self.cen_x, self.cen_y, self.cen_z,
            0.0, 0.0, 1.0)

    # def SetCameraForLabels(self):
    #     gl.glLoadIdentity()
    #     glu.gluLookAt(self.center_x,
    #                   self.center_y - self.distance, self.center_z,
    #                   self.center_x, self.center_y, self.center_z,
    #                   0.0, 0.0, 1.0)

    # def direction(self, jnt):
    #     """Returns 1 if the direction of previous r was positive otherwise 0
    #     It is used to determine the direction of shifts in expanded view.
    #     """
    #     while jnt:
    #         if jnt.r != 0:
    #             return jnt.r/abs(jnt.r)
    #         jnt = jnt.antc
    #     return 1

    @staticmethod
    def get_child_obj(jnt_obj):
        """ Returns a child frame which has the common perpendicular with
            the current frame.
        """
        for child in jnt_obj.children:
            if child.b == 0:
                return child
        if jnt_obj.children:
            return jnt_obj.children[0]
        return

    def representation(self, expanded=True):
        """ This method is used to shift joints that have the same origin
            in a way that is close to the scientific representation.
            Only joints (not frames) are shifted.
        """
        if expanded:
            for jnt in self.jnt_objs[1:]:
                i = jnt.index
                if i > self.robo.nl or self.robo.sigma[i] == 2:
                    # if cut or fixed don't shift
                    continue
                child = self.get_child_obj(jnt)
                if child is None and jnt.r == 0:
                    # The last frame (like Rx-90 6-th frame)
                    jnt.shift = 1
                elif child is not None:
                    if child.alpha != 0:
                        if self.robo.sigma[i] == 1 and child.d == 0:
                            jnt.shift = 1
                        elif jnt.r != 0 and child.d == 0 or i == 1:
                            jnt.shift = -1
                    elif child.d == 0 and child.r == 0:
                        s_child = self.get_child_obj(child)
                        if not s_child or s_child.d != 0:
                            jnt.shift = 1
                        else:
                            # jnt.shift = -2*self.length
                            # shift = -0.7*self.length
                            # jnt.shift = -(4*jnt.r + self.length)/6.
                            # child.shift = -(2*jnt.r - self.length)/6.
                            jnt.shift = 1
                            child.shift = -1

        else:
            for obj in self.jnt_objs[1:]:
                obj.shift = 0.
        self.OnDraw()

    def show_frames(self, lst):
        for jnt in self.jnt_objs:
            jnt.set_show_frame(False)
        for index in lst:
            self.jnt_objs[index].set_show_frame(True)
        self.OnDraw()

    def change_lengths(self, new_length):
        for jnt in self.jnt_objs:
            jnt.set_length(new_length)
        self.length = new_length
        self.OnDraw()

    def OnDraw(self):
        if not self.init:
            return
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
        self.base.draw()
        gl.glClear(gl.GL_DEPTH_BUFFER_BIT)
        self.base.draw_frames()
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
        gl.glFlush()
        self.SwapBuffers()

    def InitGL(self):
        # set viewing projection
        mat_specular = (1.0, 1.0, 1.0, 1.0)
        light_position = (.3, .3, 0.5, 0.0)
        light_position1 = (-0.3, 0.3, 0.5, 0.0)
        diffuseMaterial = (1., 1., 1., 1.0)
        ambientMaterial = (0.5, .5, .5, 1.0)

        gl.glClearColor(1.0, 1.0, 1.0, 1.0)
        gl.glShadeModel(gl.GL_SMOOTH)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_AMBIENT, ambientMaterial)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_DIFFUSE, diffuseMaterial)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_SPECULAR, mat_specular)
        gl.glMaterialf(gl.GL_FRONT, gl.GL_SHININESS, 25.0)
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, light_position)
        gl.glLightfv(gl.GL_LIGHT1, gl.GL_POSITION, light_position1)
        gl.glEnable(gl.GL_LIGHTING)
        gl.glEnable(gl.GL_LIGHT0)
        gl.glEnable(gl.GL_LIGHT1)

        gl.glColorMaterial(gl.GL_FRONT, gl.GL_DIFFUSE)
        gl.glEnable(gl.GL_COLOR_MATERIAL)

        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        glu.gluPerspective(40.0, 1., 0.2, 200.0)

        gl.glMatrixMode(gl.GL_MODELVIEW)
        self.CameraTransformation()


class MainWindow(wx.Frame):
    def __init__(self, prefix, robo, params=None, parent=None, identifier=-1):
        super(MainWindow, self).__init__(
            parent, identifier, prefix + ': Robot representation',
            style=wx.DEFAULT_FRAME_STYLE | wx.FULL_REPAINT_ON_RESIZE
        )
        self.robo = robo
        self.params_dict = params if params is not None else {}

        self.solve_loops = False
        self.canvas = myGLCanvas(self, robo, self.params_dict, size=(600, 600))
        self.p = wx.Panel(self)
        self.init_ui()
        self.update_spin_controls()

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer.Add(self.canvas, 1, wx.EXPAND)
        self.sizer.Add(self.p, 0, wx.EXPAND)
        self.SetSizerAndFit(self.sizer)

        self.Show()

    def init_ui(self):
        top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        gridControl = wx.GridBagSizer(hgap=10, vgap=10)
        cb = wx.CheckBox(self.p, label="Exploded View")
        cb.SetValue(True)
        gridControl.Add(cb, pos=(0, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        cb.Bind(wx.EVT_CHECKBOX, self.OnChangeRepresentation)

        self.tButton = wx.ToggleButton(self.p, label="All Frames")
        self.tButton.SetValue(True)
        self.tButton.Bind(wx.EVT_TOGGLEBUTTON, self.OnShowAllFrames)
        gridControl.Add(self.tButton, pos=(1, 0), flag=wx.ALIGN_CENTER)

        btnReset = wx.Button(self.p, label="Reset All")
        btnReset.Bind(wx.EVT_BUTTON, self.OnResetJoints)
        gridControl.Add(btnReset, pos=(3, 0), flag=wx.ALIGN_CENTER)

        btnRandom = wx.Button(self.p, label="Random")
        btnRandom.Bind(wx.EVT_BUTTON, self.OnFindRandom)
        gridControl.Add(btnRandom, pos=(4, 0), flag=wx.ALIGN_CENTER)

        self.spin_ctrls = {}
        gridJnts = wx.GridBagSizer(hgap=10, vgap=10)
        p_index = 0
        for sym in self.canvas.q_sym:
            if sym == 0:
                continue
            jnt = self.canvas.jnt_dict[sym]
            if isinstance(jnt, RevoluteJoint):
                label = 'th'
            else:
                label = 'r'
            gridJnts.Add(wx.StaticText(self.p, label=label+str(jnt.index)),
                         pos=(p_index, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            s = FS.FloatSpin(self.p, size=(70, -1), id=jnt.index,
                             increment=0.05, min_val=-10., max_val=10.)
            s.Bind(FS.EVT_FLOATSPIN, self.OnSetJointVar)
            s.SetDigits(2)
            # if sym in self.canvas.q_pas_sym:
            #     s.Enable(False)
            self.spin_ctrls[sym] = s
            gridJnts.Add(s, pos=(p_index, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            p_index += 1

        if self.robo.structure == tools.CLOSED_LOOP:
            choise_list = ['Break Loops', 'Make Loops']
            self.radioBox = wx.RadioBox(self.p, choices=choise_list,
                                        style=wx.RA_SPECIFY_ROWS)
            self.radioBox.Bind(wx.EVT_RADIOBOX, self.OnSelectLoops)
            gridControl.Add(self.radioBox, pos=(5, 0), flag=wx.ALIGN_CENTER)

        choices = []
        for jnt in self.canvas.jnt_objs:
            choices.append("Frame " + str(jnt.index))

        self.drag_pos = None
        self.check_list = wx.CheckListBox(self.p, choices=choices)
        self.check_list.SetChecked(range(len(choices)))
        self.check_list.Bind(wx.EVT_CHECKLISTBOX, self.CheckFrames)
        self.check_list.Bind(wx.EVT_LISTBOX, self.SelectFrames)
        gridControl.Add(self.check_list, pos=(2, 0), flag=wx.ALIGN_CENTER)

        q_box = wx.StaticBoxSizer(wx.StaticBox(self.p, label='Joint variables'))
        q_box.Add(gridJnts, 0, wx.ALL, 10)

        ver_sizer = wx.BoxSizer(wx.VERTICAL)
        ver_sizer.Add(q_box)
        lbl_length = wx.StaticText(self.p, label='Joint size')
        self.jnt_slider = wx.Slider(self.p, minValue=1, maxValue=100)
        self.jnt_slider.SetValue(100*self.canvas.length)
        self.jnt_slider.Bind(wx.EVT_SCROLL, self.OnSliderChanged)
        ver_sizer.Add(lbl_length, 0, wx.TOP | wx.ALIGN_CENTER_HORIZONTAL, 15)
        ver_sizer.Add(self.jnt_slider, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        top_sizer.Add(ver_sizer, 0, wx.ALL, 10)
        top_sizer.AddSpacer(10)
        top_sizer.Add(gridControl, 0, wx.ALL, 10)

        # button1 = wx.Button(self.panel, label="TEXT 1")
        # button2 = wx.Button(self.panel, label="TEXT 2")
        # check1 = wx.CheckBox(self.panel, label="Show Axes")
        # self.panel.Bind(wx.EVT_CHECKBOX, None)
        #
        # sizer = wx.BoxSizer(wx.VERTICAL)
        # sizer.Add(button1, flag=wx.BOTTOM, border=5)
        # sizer.Add(button2, flag=wx.BOTTOM, border=5)
        # sizer.Add(check1)
        #
        # border = wx.BoxSizer()
        # border.Add(sizer, flag=wx.ALL | wx.EXPAND, border=5)

        self.p.SetSizerAndFit(top_sizer)

    def OnChangeRepresentation(self, evt):
        self.canvas.representation(evt.EventObject.GetValue())

    def OnShowWorldFrame(self, evt):
        self.canvas.base.set_show_frame(evt.EventObject.GetValue())
        self.canvas.OnDraw()

    def OnShowAllFrames(self, _):
        """Shows or hides all the frames (Toggle button event handler)
        """
        if self.tButton.Value:
            indices = range(len(self.canvas.jnt_objs))
        else:
            indices = []
        self.canvas.show_frames(indices)
        self.check_list.SetChecked(indices)

    def update_spin_controls(self):
        for ctrl in self.spin_ctrls.values():
            ctrl.SetValue(self.canvas.jnt_objs[ctrl.Id].q)

    def OnSelectLoops(self, _):
        for q in self.canvas.q_pas_sym:
            self.spin_ctrls[q].Enable(not self.radioBox.Selection)
        if self.radioBox.Selection == 1:
            self.solve_loops = True
            self.canvas.solve()
            self.update_spin_controls()
            self.canvas.OnDraw()
        else:
            self.solve_loops = False

    def OnResetJoints(self, _):
        for ctrl in self.spin_ctrls.values():
            jnt_obj = self.canvas.jnt_objs[ctrl.Id]
            jnt_obj.q = jnt_obj.q_init
        if self.solve_loops:
            self.canvas.solve()
        self.update_spin_controls()
        self.canvas.OnDraw()

    def OnFindRandom(self, evt):
        pass

    def OnSetJointVar(self, evt):
        """Sets joint values from the spin-controls
        """
        jnt_index = evt.GetId()
        jnt = self.canvas.jnt_objs[jnt_index]
        value = evt.EventObject.GetValue()
        if self.robo.sigma[jnt_index] == 0:
            value = atan2(sin(value), cos(value))
            evt.EventObject.SetValue(value)
        jnt.q = value
        if self.solve_loops:
            self.canvas.solve()
            self.update_spin_controls()
        self.canvas.OnDraw()

    def CheckFrames(self, evt):
        self.canvas.show_frames(evt.EventObject.GetChecked())
        self.tButton.Value = False
        evt.EventObject.DeselectAll()

    def SelectFrames(self, evt):
        selections = evt.EventObject.GetSelections()
        if selections:
            self.canvas.centralize_to_frame(selections[0])

    def OnSliderChanged(self, _):
        self.canvas.change_lengths(self.jnt_slider.Value/100.)


if __name__ == '__main__':
    app = wx.PySimpleApp()
    from pysymoro import robot
    robo = samplerobots.rx90()
    robo.d[3] = 1.
    robo.r[4] = 1.
    frame = MainWindow(prefix='', robo=robo)
    import profile
    profile.run('frame.canvas.OnPaintAll(wx.CommandEvent())', sort='cumtime')
    app.MainLoop()
    app.Destroy()
