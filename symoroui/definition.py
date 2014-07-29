# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module creates the dialog box to define a new robot and
visualisation parameters.
"""


import wx
import wx.lib.scrolledpanel as scrolled

from sympy import Expr, Symbol

from symoroutils import tools


class DialogDefinition(wx.Dialog):
    """Creates the dialog box to define a new robot."""
    def __init__(
        self, prefix, name, nl, nj, structure, is_floating,
        is_mobile, parent=None
    ):
        super(DialogDefinition, self).__init__(
            parent,
            style=wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN
        )
        self.init_ui(name, nl, nj, is_floating, is_mobile, structure)
        self.SetTitle(prefix + ": New robot definition")

    def init_ui(self, name, nl, nj, is_floating, is_mobile, structure):
        szr_topmost = wx.BoxSizer(wx.VERTICAL)
        # title
        szr_topmost.Add(
            wx.StaticText(self, label="Robot definition"), 0,
            wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 25
        )
        # grid
        szr_grd = wx.GridBagSizer(15, 15)
        szr_grd.Add(
            wx.StaticText(self, label='Name of the robot:'),
            pos=(0, 0), flag=wx.BOTTOM | wx.ALIGN_LEFT, border=2
        )
        szr_grd.Add(
            wx.TextCtrl(self, size=(92, -1),
            name='name', value=name), pos=(0, 1)
        )
        szr_grd.Add(
            wx.StaticText(self, label='Number of moving links:'),
            pos=(1, 0), flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT,
            border=2
        )
        self.spin_links = wx.SpinCtrl(
            self, size=(92, -1), value=str(nl), min=1
        )
        self.spin_links.Bind(wx.EVT_SPINCTRL, self.OnSpinNL)
        szr_grd.Add(self.spin_links, pos=(1, 1))
        szr_grd.Add(
            wx.StaticText(self, label='Number of joints:'), pos=(2, 0),
            flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2
        )
        self.spin_joints = wx.SpinCtrl(
            self, size=(92, -1), value=str(nj), min=0
        )
        szr_grd.Add(self.spin_joints, pos=(2, 1))
        szr_grd.Add(
            wx.StaticText(self, label='Type of structure'), pos=(3, 0),
            flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2
        )
        self.cmb_structure = wx.ComboBox(
            self, size=(92, -1), name='structure', style=wx.CB_READONLY,
            choices=[tools.SIMPLE, tools.TREE, tools.CLOSED_LOOP],
            value=structure
        )
        szr_grd.Add(self.cmb_structure, pos=(3, 1))
        self.cmb_structure.Bind(wx.EVT_COMBOBOX, self.OnTypeChanged)
        self.OnTypeChanged(None)
        szr_topmost.Add(
            szr_grd, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 15
        )
        # radio box for base type
        self.rbx_base_type = wx.RadioBox(
            self, label='Base Type', style=wx.RA_SPECIFY_COLS,
            majorDimension=3, name='rbx_base',
            choices=['Fixed', 'Floating', 'Mobile']
        )
        self.rbx_base_type.Bind(wx.EVT_RADIOBOX, self.OnBaseType)
        self.chk_keep_geo = wx.CheckBox(
            self, label=' Keep geometric parameters'
        )
        self.chk_keep_geo.Value = True
        self.chk_keep_dyn = wx.CheckBox(
            self, label=' Keep dynamic parameters'
        )
        self.chk_keep_dyn.Value = True
        self.chk_keep_joint = wx.CheckBox(
            self, label=' Keep joint parameters'
        )
        self.chk_keep_joint.Value = True
        self.chk_keep_base = wx.CheckBox(
            self, label=' Keep base parameters'
        )
        self.chk_keep_base.Value = True
        szr_topmost.Add(
            self.rbx_base_type, 0,
            wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER, 15
        )
        szr_topmost.Add(
            self.chk_keep_geo, 0,
            wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15
        )
        szr_topmost.Add(
            self.chk_keep_dyn, 0,
            wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15
        )
        szr_topmost.Add(
            self.chk_keep_joint, 0,
            wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15
        )
        szr_topmost.Add(
            self.chk_keep_base, 0,
            wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15
        )
        szr_horizontal = wx.BoxSizer(wx.HORIZONTAL)
        btn_ok = wx.Button(self, wx.ID_OK, "OK")
        btn_ok.Bind(wx.EVT_BUTTON, self.OnOK)
        btn_cancel = wx.Button(self, wx.ID_CANCEL, "Cancel")
        btn_cancel.Bind(wx.EVT_BUTTON, self.OnCancel)
        szr_horizontal.Add(btn_ok, 0, wx.ALL, 25)
        szr_horizontal.Add(btn_cancel, 0, wx.ALL, 25)
        szr_topmost.Add(szr_horizontal)
        self.SetSizerAndFit(szr_topmost)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def OnTypeChanged(self, _):
        if self.cmb_structure.GetSelection() == 2:
            self.spin_joints.Enable(True)
        else:
            self.spin_joints.Enable(False)
            self.OnSpinNL(None)

    def OnBaseType(self, _):
        idx = self.rbx_base_type.GetSelection()
        value = True if idx == 0 else False
        self.chk_keep_base.Value = value

    def OnSpinNL(self, _):
        self.spin_joints.SetRange(int(self.spin_links.Value), 100)
        if self.cmb_structure.GetSelection() != 2:
            self.spin_joints.Value = self.spin_links.Value

    def get_values(self):
        name = self.FindWindowByName('name').Value
        nl = int(self.spin_links.Value)
        nj = int(self.spin_joints.Value)
        base_type_idx = self.rbx_base_type.GetSelection()
        is_floating = True if base_type_idx == 1 else False
        is_mobile = True if base_type_idx == 2 else False
        params = {
            'name': name,
            'num_links': nl,
            'num_joints': nj,
            'num_frames': (2 * nj) -nl,
            'structure': self.cmb_structure.Value,
            'is_floating': is_floating,
            'is_mobile': is_mobile,
            'keep_geo': self.chk_keep_geo.Value,
            'keep_dyn': self.chk_keep_dyn.Value,
            'keep_joint': self.chk_keep_joint.Value,
            'keep_base': self.chk_keep_base.Value
        }
        return params


class DialogVisualisation(wx.Dialog):
    """Creates the dialog box to define the visualisation parameters."""
    def __init__(self, prefix, robo, par_dict, parent=None):
        super(DialogVisualisation, self).__init__(parent)
        self.robo = robo
        self.par_dict = par_dict
        self.SetTitle(prefix + ": Enter numerical values")
        self.construct_sym()
        self.init_ui()

    def has_syms(self):
        return len(self.syms) > 0

    def construct_sym(self):
        params = ['gamma', 'b', 'alpha', 'd', 'theta', 'r']
        q_vec = self.robo.q_vec
        self.syms = set()
        for par in params:
            for i in range(1, self.robo.NF):
                val = self.robo.get_val(i, par)
                if val in q_vec:
                    continue
                if isinstance(val, Expr):
                    for at in val.atoms(Symbol):
                        self.syms.add(at)

    def init_ui(self):
        panel = scrolled.ScrolledPanel(self, -1)
        szr_vertical = wx.BoxSizer(wx.VERTICAL)
        self.widgets = {}
        for par in self.syms:
            if par in self.par_dict:
                val = str(self.par_dict[par])
            else:
                val = 1.
            szr_horizontal = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(
                panel, label=str(par),
                size=(60, -1), style=wx.ALIGN_RIGHT
            )
            szr_horizontal.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            szr_horizontal.AddSpacer(5)
            txt_box = wx.TextCtrl(panel, size=(120, -1), value=str(val))
            szr_horizontal.Add(txt_box)
            szr_vertical.Add(
                szr_horizontal, 0,
                wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5
            )
            self.widgets[str(par)] = txt_box
        panel.SetSizer(szr_vertical)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        szr_topmost = wx.BoxSizer(wx.VERTICAL)
        szr_topmost.Add(panel, 1, wx.ALL | wx.EXPAND, 2)
        szr_horizontal = wx.BoxSizer(wx.HORIZONTAL)
        btn_ok = wx.Button(self, wx.ID_OK, "OK")
        btn_ok.Bind(wx.EVT_BUTTON, self.OnOK)
        btn_cancel = wx.Button(self, wx.ID_CANCEL, "Cancel")
        btn_cancel.Bind(wx.EVT_BUTTON, self.OnCancel)
        szr_horizontal.Add(btn_ok, 0, wx.ALL, 15)
        szr_horizontal.Add(btn_cancel, 0, wx.ALL, 15)
        szr_topmost.Add(szr_horizontal, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.SetSizer(szr_topmost)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def get_values(self):
        result = {}
        for par in self.syms:
            result[par] = self.widgets[str(par)].Value
        return result


