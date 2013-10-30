__author__ = 'Izzat'
import wx
import wx.lib.scrolledpanel as scrolled
from core.symoro import SIMPLE, TREE, CLOSED_LOOP
from sympy import Expr, Symbol
#TODO: PROG_NAME


class DialogConversion(wx.Dialog):
    def __init__(self, prefix, robo, par_dict, parent=None):
        super(DialogConversion, self).__init__(parent)
        self.robo = robo
        self.par_dict = par_dict
        self.SetTitle(prefix + ": Enter numerical values")
        self.construct_sym()
        self.InitUI()

    def HasSyms(self):
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

    def InitUI(self):
        p = scrolled.ScrolledPanel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.widgets = {}
        for par in self.syms:
            if par in self.par_dict:
                val = str(self.par_dict[par])
            else:
                val = 1.
            horSizer = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(p, label=str(par), size=(60, -1), style=wx.ALIGN_RIGHT)
            horSizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            horSizer.AddSpacer(5)
            textBox = wx.TextCtrl(p, size=(120, -1), value=str(val))
            horSizer.Add(textBox)
            vbox.Add(horSizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
            self.widgets[str(par)] = textBox

        p.SetSizer(vbox)
        p.SetAutoLayout(1)
        p.SetupScrolling()

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(p, 1, wx.ALL | wx.EXPAND, 2)
        horSizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, wx.ID_OK, "OK")
        okButton.Bind(wx.EVT_BUTTON, self.OnOK)
        cancelButton = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancelButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        horSizer.Add(okButton, 0, wx.ALL, 15)
        horSizer.Add(cancelButton, 0, wx.ALL, 15)
        mainSizer.Add(horSizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.SetSizer(mainSizer)
        #if n_syms == 0:
        #   self.EndModal(wx.ID_OK)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def GetValues(self):
        result = {}
        for par in self.syms:
            result[par] = self.widgets[str(par)].Value
        return result


class DialogDefinition(wx.Dialog):
    def __init__(self, prefix, name, NL, NJ, structure, is_mobile, parent=None):
        super(DialogDefinition, self).__init__(parent, style=wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN)
        self.InitUI(name, NL, NJ, is_mobile, structure)
        self.SetTitle(prefix + ": New robot definition")

    def InitUI(self, name, NL, NJ, is_mobile, structure):
        mainSizer = wx.BoxSizer(wx.VERTICAL)

        #title
        label_main = wx.StaticText(self, label="Robot definition")
        mainSizer.Add(label_main, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 25)

        #grid
        grid = wx.GridBagSizer(15, 15)
        grid.Add(wx.StaticText(self, label='Name of the robot:'), pos=(0, 0), flag=wx.BOTTOM | wx.ALIGN_LEFT, border=2)
        grid.Add(wx.TextCtrl(self, size=(92, -1), name='name', value=name), pos=(0, 1))

        label = wx.StaticText(self, label='Number of moving links:')
        grid.Add(label, pos=(1, 0), flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2)
        self.spin_links = wx.SpinCtrl(self, size=(92, -1), value=str(NL), min=1)
        self.spin_links.Bind(wx.EVT_SPINCTRL, self.OnSpinNL)
        grid.Add(self.spin_links, pos=(1, 1))

        label = wx.StaticText(self, label='Number of joints:')
        grid.Add(label, pos=(2, 0), flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2)
        self.spin_joints = wx.SpinCtrl(self, size=(92, -1), value=str(NJ), min=0)
        grid.Add(self.spin_joints, pos=(2, 1))

        grid.Add(wx.StaticText(self, label='Type of structure'), pos=(3, 0),
                 flag=wx.BOTTOM | wx.TOP | wx.ALIGN_LEFT, border=2)
        self.cmb_structure = wx.ComboBox(self, size=(92, -1), name='structure', style=wx.CB_READONLY,
                             choices=[SIMPLE, TREE, CLOSED_LOOP], value=structure)
        grid.Add(self.cmb_structure, pos=(3, 1))
        self.cmb_structure.Bind(wx.EVT_COMBOBOX, self.OnTypeChanged)
        self.OnTypeChanged(None)

        mainSizer.Add(grid, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 15)
        self.check_is_mobile = wx.CheckBox(self, label=' Is mobile')
        self.check_is_mobile.Value = is_mobile
        self.check_keep_geo = wx.CheckBox(self, label=' Keep geometric parameters')
        self.check_keep_geo.Value = True
        self.check_keep_dyn = wx.CheckBox(self, label=' Keep dynamic parameters')
        self.check_keep_dyn.Value = True
        self.check_keep_base = wx.CheckBox(self, label=' Keep base parameters')
        self.check_keep_base.Value = True
        mainSizer.Add(self.check_is_mobile, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, 15)
        mainSizer.Add(self.check_keep_geo, 0, wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15)
        mainSizer.Add(self.check_keep_dyn, 0, wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15)
        mainSizer.Add(self.check_keep_base, 0, wx.TOP | wx.LEFT | wx.ALIGN_LEFT, 15)
        horSizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, wx.ID_OK, "OK")
        okButton.Bind(wx.EVT_BUTTON, self.OnOK)
        cancelButton = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancelButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        horSizer.Add(okButton, 0, wx.ALL, 25)
        horSizer.Add(cancelButton, 0, wx.ALL, 25)
        mainSizer.Add(horSizer)

        self.SetSizerAndFit(mainSizer)

    def OnOK(self, e):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, e):
        self.EndModal(wx.ID_CANCEL)

    def OnTypeChanged(self, evt):
        if self.cmb_structure.GetSelection() == 2:
            self.spin_joints.Enable(True)
        else:
            self.spin_joints.Enable(False)
            self.OnSpinNL(None)

    def OnSpinNL(self, evt):
        self.spin_joints.SetRange(int(self.spin_links.Value), 100)
        if self.cmb_structure.GetSelection() != 2:
            self.spin_joints.Value = self.spin_links.Value

    def GetValues(self):
        name = self.FindWindowByName('name').Value
        NL = int(self.spin_links.Value)
        NJ = int(self.spin_joints.Value)
        return {'init_pars' : (name, NL, NJ, 2*NJ - NL, self.check_is_mobile.Value, self.cmb_structure.Value),
                'keep_geo' : self.check_keep_geo.Value, 'keep_dyn' : self.check_keep_dyn.Value,
                'keep_base' : self.check_keep_base.Value}

