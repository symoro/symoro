# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module creates the dialog box for differrent kinematic model
parameters.
"""
from sympy import zeros

from symoroutils import symbolmgr
from symoroutils.tools import ZERO, FAIL
from pysymoro.geometry import compute_rot_trans, Z_AXIS
from pysymoro import kinematics
from symoroutils.paramsinit import ParamsInit
import wx


def velocities(robo):
    symo = symbolmgr.SymbolManager(None)
    symo.file_open(robo, 'vel')
    symo.write_params_table(robo, 'Link velocities')
    antRj, antPj = compute_rot_trans(robo, symo)
    w = ParamsInit.init_w(robo)
    v = ParamsInit.init_v(robo)
    for j in xrange(1, robo.NL):
        jRant = antRj[j].T
        qdj = Z_AXIS * robo.qdot[j]
        kinematics._omega_ij(robo, j, jRant, w, qdj)
        symo.mat_replace(w[j], 'W', j, forced=True)
        kinematics._v_j(robo, j, antPj, jRant, v, w, qdj)
        symo.mat_replace(v[j], 'V', j, forced=True)
    symo.file_close()
    return symo


def accelerations(robo):
    symo = symbolmgr.SymbolManager(None)
    symo.file_open(robo, 'acc')
    symo.write_params_table(robo, 'Link accelerations')
    antRj, antPj = compute_rot_trans(robo, symo)
    kinematics.compute_vel_acc(robo, symo, antRj, antPj, forced=True)
    symo.file_close()
    return symo


#Similar to compute_vel_acc.
def jdot_qdot(robo):
    symo = symbolmgr.SymbolManager(None)
    symo.file_open(robo, 'jpqp')
    symo.write_params_table(robo, 'JdotQdot')
    antRj, antPj = compute_rot_trans(robo, symo)
    w = ParamsInit.init_w(robo)
    wdot, vdot = ParamsInit.init_wv_dot(robo, gravity=False)
    U = ParamsInit.init_u(robo)
    for j in xrange(1, robo.NL):
        jRant = antRj[j].T
        qdj = Z_AXIS * robo.qdot[j]
        qddj = Z_AXIS * ZERO
        wi, w[j] = kinematics._omega_ij(robo, j, jRant, w, qdj)
        symo.mat_replace(w[j], 'W', j)
        symo.mat_replace(wi, 'WI', j)
        kinematics._omega_dot_j(robo, j, jRant, w, wi, wdot, qdj, qddj)
        symo.mat_replace(wdot[j], 'WPJ', j, forced=True)
        kinematics._v_dot_j(robo, symo, j, jRant, antPj, w,
                            wi, wdot, U, vdot, qdj, qddj)
        symo.mat_replace(vdot[j], 'VPJ', j, forced=True)
    symo.file_close()
    return symo


def jacobian(robo, n, i, j):
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'jac')
    title = "Jacobian matrix for frame {}\n"
    title += "Projection frame {}, intermediate frame {}"
    symo.write_params_table(robo, title.format(n, i, j))
    kinematics._jac(robo, symo, n, i, j, forced=True)
    symo.file_close()
    return symo


def jacobian_determinant(robo, n, i, j, rows, cols):
    symo = symbolmgr.SymbolManager(None)
    J, L = kinematics._jac(robo, symo, n, i, j, trig_subs=False)
    J_reduced = zeros(len(rows), len(cols))
    for i, i_old in enumerate(rows):
        for j, j_old in enumerate(cols):
            J_reduced[i, j] = J[i_old, j_old]
    symo.file_open(robo, 'det')
    symo.write_params_table(robo, 'Jacobian determinant for frame %s' % n)
    symo.write_line(kinematics._jac_det(robo, symo, J=J_reduced))
    symo.file_close()
    return symo


def kinematic_constraints(robo):
    symo = symbolmgr.SymbolManager(None)
    res = kinematics.kin_loop_constr(robo, symo)
    if res == FAIL:
        return FAIL
    W_a, W_p, W_ac, W_pc, W_c, _ = res
    symo.file_open(robo, 'ckel')
    symo.write_params_table(robo, 'Constraint kinematic equations of loop',
                            equations=False)
    symo.write_line('Active joint variables')
    symo.write_line([robo.get_q(i) for i in robo.indx_active])
    symo.write_line()
    symo.write_line('Passive joints variables')
    symo.write_line([robo.get_q(i) for i in robo.indx_passive])
    symo.write_line()
    symo.write_line('Cut joints variables')
    symo.write_line([robo.get_q(i) for i in robo.indx_cut])
    symo.write_line()
    symo.mat_replace(W_a, 'WA', forced=True)
    symo.mat_replace(W_p, 'WP', forced=True)
    symo.mat_replace(W_ac, 'WPA', forced=True)
    symo.mat_replace(W_pc, 'WPC', forced=True)
    symo.mat_replace(W_c, 'WC', forced=True)
    symo.file_close()
    return symo


class DialogJacobian(wx.Dialog):
    """
    Creates the dialog box to specify parameters for the calculation
    of the Jacobian matrix.
    """
    def __init__(self, prefix, robo, parent=None):
        super(DialogJacobian, self).__init__(parent, style=wx.SYSTEM_MENU |
                              wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN)
        self.robo = robo
        self.init_ui()
        self.SetTitle(prefix + ": Jacobian matrix (jac)")

    def init_ui(self):
        sizer = wx.BoxSizer(wx.VERTICAL)
        #title
        label_main = wx.StaticText(self, label="Calculation of i Jr j")
        sizer.Add(label_main, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 25)
        #input
        choices = [str(i) for i in range(self.robo.NF)]
        label = wx.StaticText(self, label='Frame number ( r )')
        sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 2)
        self.cmb_frame = wx.ComboBox(self, size=(50, -1),
                                     choices=choices, style=wx.CB_READONLY)
        self.cmb_frame.SetSelection(len(choices)-1)
        self.cmb_frame.Bind(wx.EVT_COMBOBOX, self.OnFrameChanged)
        sizer.Add(self.cmb_frame, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 10)
        chain = self.robo.chain(int(self.cmb_frame.Value))
        choices = [str(i) for i in reversed(chain + [0])]
        label = wx.StaticText(self, label='Projection frame ( i )')
        sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 2)
        self.cmb_proj = wx.ComboBox(self, size=(50, -1),
                                    choices=choices, style=wx.CB_READONLY)
        self.cmb_proj.SetSelection(0)
        sizer.Add(self.cmb_proj, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 10)
        label = wx.StaticText(self, label='Intermediate frame ( j )')
        sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 2)
        self.cmb_inter = wx.ComboBox(self, size=(50, -1),
                                     choices=choices, style=wx.CB_READONLY)
        self.cmb_inter.SetSelection(0)
        sizer.Add(self.cmb_inter, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 10)
        hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        hor_sizer.Add(ok_btn, 0, wx.ALL, 25)
        hor_sizer.Add(cancel_btn, 0, wx.ALL, 25)
        sizer.Add(hor_sizer)
        self.SetSizerAndFit(sizer)

    def OnFrameChanged(self, _):
        chain = self.robo.chain(int(self.cmb_frame.Value))
        choices = [str(i) for i in reversed(chain + [0])]
        self.cmb_inter.SetItems(choices)
        self.cmb_inter.SetSelection(0)
        self.cmb_proj.SetItems(choices)
        self.cmb_proj.SetSelection(0)

    def OnOK(self, _):
        n = int(self.cmb_frame.Value)
        i = int(self.cmb_proj.Value)
        j = int(self.cmb_inter.Value)
        self.symo = jacobian(self.robo, n, i, j)
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)


class DialogDeterminant(wx.Dialog):
    """
    Creates the dialog box to specify parameters for the
    calculation of the determinant of a Jacobian.
    """
    def __init__(self, prefix, robo, parent=None):
        st = wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN
        super(DialogDeterminant, self).__init__(parent, style=st)
        self.robo = robo
        self.init_ui()
        self.SetTitle(prefix + ": Determinant of a jacobian matrix (det)")

    def init_ui(self):
        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        #title
        mainSizer.AddSpacer(10)
        grid = wx.GridBagSizer(5, 15)
        choices = [str(i) for i in range(self.robo.NF)]
        label = wx.StaticText(self, label='Frame number ( r )')
        grid.Add(label, pos=(0, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL)
        self.cmb_frame = wx.ComboBox(self, size=(50, -1),
                                     choices=choices, style=wx.CB_READONLY)
        self.cmb_frame.SetSelection(len(choices)-1)
        self.cmb_frame.Bind(wx.EVT_COMBOBOX, self.OnFrameChanged)
        grid.Add(self.cmb_frame, pos=(1, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, border=15)
        chain = self.robo.chain(int(self.cmb_frame.Value))
        choices = [str(i) for i in reversed(chain + [0])]
        label = wx.StaticText(self, label='Projection frame ( j )')
        grid.Add(label, pos=(2, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL)
        self.cmb_proj = wx.ComboBox(self, size=(50, -1),
                                    choices=choices, style=wx.CB_READONLY)
        self.cmb_proj.SetSelection(len(choices)-1)
        self.cmb_proj.SetSelection(0)
        grid.Add(self.cmb_proj, pos=(3, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, border=15)
        label = wx.StaticText(self, label='Intermediate frame ( i )')
        grid.Add(label, pos=(4, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL)
        self.cmb_inter = wx.ComboBox(self, size=(50, -1),
                                     choices=choices, style=wx.CB_READONLY)
        self.cmb_inter.SetSelection(0)
        grid.Add(self.cmb_inter, pos=(5, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, border=15)
        label_main = wx.StaticText(self, label="Definition of sub-matrix")
        grid.Add(label_main, pos=(6, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL)
        label_main = wx.StaticText(self,
                     label="(Select rows and columns to be deleted)")
        grid.Add(label_main, pos=(7, 0), span=(1, 2),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, border=5)
        #input
        label_row = wx.StaticText(self, label="Rows:")
        label_col = wx.StaticText(self, label="Columns:")
        self.box_row = wx.CheckListBox(self, size=(100, 100),
                       choices=['VX', 'VY', 'VZ', 'WX', 'WY', 'WZ'])
        self.box_col = wx.CheckListBox(self, size=(100, 100),
                       choices=[str(i) for i in self.robo.q_vec])
        grid.Add(label_row, pos=(8, 0))
        grid.Add(label_col, pos=(8, 1))
        grid.Add(self.box_row, pos=(9, 0))
        grid.Add(self.box_col, pos=(9, 1))
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        grid.Add(ok_btn, pos=(10, 0),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, border=5)
        grid.Add(cancel_btn, pos=(10, 1),
                 flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, border=5)
        mainSizer.Add(grid, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 15)
        self.SetSizerAndFit(mainSizer)
        self.OnFrameChanged(None)

    def OnFrameChanged(self, _):
        chain = self.robo.chain(int(self.cmb_frame.Value))
        choices = [str(i) for i in reversed(chain + [0])]
        self.cmb_inter.SetItems(choices)
        self.cmb_inter.SetSelection(0)
        self.cmb_proj.SetItems(choices)
        self.cmb_proj.SetSelection(0)
        choices_list = self.robo.get_q_chain(int(self.cmb_frame.Value))
        self.box_col.SetItems([str(i) for i in choices_list])

    def OnOK(self, _):
        rows = [i for i in range(len(self.box_row.Items))
                if not self.box_row.IsChecked(i)]
        cols = [i for i in range(len(self.box_col.Items))
                if not self.box_col.IsChecked(i)]
        n = int(self.cmb_frame.Value)
        i = int(self.cmb_proj.Value)
        j = int(self.cmb_inter.Value)
        self.symo = jacobian_determinant(self.robo, n, i, j, rows, cols)
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)


