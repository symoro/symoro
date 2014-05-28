#!/usr/bin/env python
# -*- coding: utf-8 -*-


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
    def __init__(self, prefix, name, nl, nj, structure, is_mobile, parent=None):
        super(DialogDefinition, self).__init__(parent, style=wx.SYSTEM_MENU |
                                wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN)
        self.init_ui(name, nl, nj, is_mobile, structure)
        self.SetTitle(prefix + ": New robot definition")

    def init_ui(self, name, nl, nj, is_mobile, structure):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        #title
        label_main = wx.StaticText(self, label="Robot definition")
        main_sizer.Add(label_main, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 25)
        #grid
        grid = wx.GridBagSizer(15, 15)
        grid.Add(wx.StaticText(self, label='Name of the robot:'), pos=(0, 0),
                 flag=wx.BOTTOM | wx.ALIGN_LEFT, border=2)
        grid.Add(wx.TextCtrl(self, size=(92, -1), name='name', value=name),
                 pos=(0, 1))
        label = wx.StaticText(self, label='Number of moving links:')
        grid.Add(label, pos=(1, 0),
                 flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2)
        self.spin_links = wx.SpinCtrl(self, size=(92, -1),
                                      value=str(nl), min=1)
        self.spin_links.Bind(wx.EVT_SPINCTRL, self.OnSpinNL)
        grid.Add(self.spin_links, pos=(1, 1))
        label = wx.StaticText(self, label='Number of joints:')
        grid.Add(label, pos=(2, 0),
                 flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2)
        self.spin_joints = wx.SpinCtrl(self, size=(92, -1),
                                       value=str(nj), min=0)
        grid.Add(self.spin_joints, pos=(2, 1))
        grid.Add(wx.StaticText(self, label='Type of structure'), pos=(3, 0),
                 flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2)
        self.cmb_structure = wx.ComboBox(
            self, size=(92, -1), name='structure', style=wx.CB_READONLY,
            choices=[tools.SIMPLE, tools.TREE, tools.CLOSED_LOOP],
            value=structure
        )
        grid.Add(self.cmb_structure, pos=(3, 1))
        self.cmb_structure.Bind(wx.EVT_COMBOBOX, self.OnTypeChanged)
        self.OnTypeChanged(None)
        main_sizer.Add(grid, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 15)
        self.ch_is_mobile = wx.CheckBox(self, label=' Is mobile')
        self.ch_is_mobile.Value = is_mobile
        self.ch_keep_geo = wx.CheckBox(self, label=' Keep geometric parameters')
        self.ch_keep_geo.Value = True
        self.ch_keep_dyn = wx.CheckBox(self, label=' Keep dynamic parameters')
        self.ch_keep_dyn.Value = True
        self.ch_keep_base = wx.CheckBox(self, label=' Keep base parameters')
        self.ch_keep_base.Value = True
        main_sizer.Add(self.ch_is_mobile, 0,
                       wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, 15)
        main_sizer.Add(self.ch_keep_geo, 0,
                       wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15)
        main_sizer.Add(self.ch_keep_dyn, 0,
                       wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15)
        main_sizer.Add(self.ch_keep_base, 0,
                       wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15)
        hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        hor_sizer.Add(ok_btn, 0, wx.ALL, 25)
        hor_sizer.Add(cancel_btn, 0, wx.ALL, 25)
        main_sizer.Add(hor_sizer)
        self.SetSizerAndFit(main_sizer)

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

    def OnSpinNL(self, _):
        self.spin_joints.SetRange(int(self.spin_links.Value), 100)
        if self.cmb_structure.GetSelection() != 2:
            self.spin_joints.Value = self.spin_links.Value

    def get_values(self):
        name = self.FindWindowByName('name').Value
        nl = int(self.spin_links.Value)
        nj = int(self.spin_joints.Value)
        return {'init_pars': (name, nl, nj, 2*nj - nl,
                               self.ch_is_mobile.Value,
                               self.cmb_structure.Value),
                'keep_geo': self.ch_keep_geo.Value,
                'keep_dyn': self.ch_keep_dyn.Value,
                'keep_base': self.ch_keep_base.Value}


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
        params = self.robo.get_geom_head()[4:]
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
        p = scrolled.ScrolledPanel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.widgets = {}
        for par in self.syms:
            if par in self.par_dict:
                val = str(self.par_dict[par])
            else:
                val = 1.
            hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(p, label=str(par), size=(60, -1),
                                  style=wx.ALIGN_RIGHT)
            hor_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.AddSpacer(5)
            txt_box = wx.TextCtrl(p, size=(120, -1), value=str(val))
            hor_sizer.Add(txt_box)
            vbox.Add(hor_sizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
            self.widgets[str(par)] = txt_box
        p.SetSizer(vbox)
        p.SetAutoLayout(1)
        p.SetupScrolling()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(p, 1, wx.ALL | wx.EXPAND, 2)
        hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        hor_sizer.Add(ok_btn, 0, wx.ALL, 15)
        hor_sizer.Add(cancel_btn, 0, wx.ALL, 15)
        main_sizer.Add(hor_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.SetSizer(main_sizer)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def get_values(self):
        result = {}
        for par in self.syms:
            result[par] = self.widgets[str(par)].Value
        return result


