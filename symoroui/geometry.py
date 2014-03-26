#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module creates the dialog boxes to specify parameters for
geometric model calculations.
"""


import wx


class DialogTrans(wx.Dialog):
    """Creates the dialog box for transformation matrix selection."""
    def __init__(self, prefix, nf, parent=None):
        super(DialogTrans, self).__init__(parent, style=wx.SYSTEM_MENU |
                           wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN)
        self.nf = nf
        self.init_ui()
        self.SetTitle(prefix + ": Transformation matrix (trm)")

    def init_ui(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        main_sizer.Add(hor_sizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 25)
        # Insert, Delete buttons and comboboxes
        grid = wx.GridBagSizer(hgap=40, vgap=11)
        lab_top = wx.StaticText(self, label='Original frame')
        lab_bottom = wx.StaticText(self, label='Destination frame')
        grid.Add(lab_top, pos=(0, 0), flag=wx.ALIGN_CENTER_HORIZONTAL)
        grid.Add(lab_bottom, pos=(2, 0), flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.cmb_i = wx.ComboBox(self, size=(80, -1),
             choices=[str(i) for i in range(self.nf)], style=wx.CB_READONLY)
        self.cmb_i.SetSelection(0)
        self.cmb_j = wx.ComboBox(self, size=(80, -1),
             choices=[str(i) for i in range(self.nf)], style=wx.CB_READONLY)
        self.cmb_j.SetSelection(self.nf - 1)
        grid.Add(self.cmb_i, pos=(1, 0),
                 flag=wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, border=40)
        grid.Add(self.cmb_j, pos=(3, 0),
                 flag=wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, border=10)
        insert_btn = wx.Button(self, wx.ID_ANY, "Insert ->")
        insert_btn.Bind(wx.EVT_BUTTON, self.OnInsert)
        delete_btn = wx.Button(self, wx.ID_ANY, "Delete <-")
        delete_btn.Bind(wx.EVT_BUTTON, self.OnDelete)
        grid.Add(insert_btn, pos=(1, 1))
        grid.Add(delete_btn, pos=(3, 1))
        hor_sizer.Add(grid)
        # Transformations label and list
        ver_sizer = wx.BoxSizer(wx.VERTICAL)
        hor_sizer.AddSpacer(40)
        hor_sizer.Add(ver_sizer)
        lab_trans = wx.StaticText(self, label='Transformations')
        ver_sizer.Add(lab_trans, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        ver_sizer.AddSpacer(10)
        self.listbox = wx.ListBox(self, size=(80, 120), style=wx.LB_SINGLE)
        self.result = set()
        ver_sizer.Add(self.listbox)
        self.check_short = wx.CheckBox(self, label='Trigonometric short form')
        main_sizer.Add(self.check_short, 0, wx.LEFT | wx.ALIGN_LEFT, 25)
        main_sizer.AddSpacer(15)
        # OK Cancel
        hor_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        main_sizer.Add(hor_sizer2, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 12)
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        hor_sizer2.Add(ok_btn)
        hor_sizer2.AddSpacer(22)
        hor_sizer2.Add(cancel_btn)
        self.SetSizerAndFit(main_sizer)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def OnInsert(self, _):
        trans = (int(self.cmb_i.Value), int(self.cmb_j.Value))
        if trans not in self.result:
            self.listbox.Insert(str(trans), len(self.result), trans)
            self.result.add(trans)

    def OnDelete(self, _):
        sel_index = self.listbox.GetSelection()
        if sel_index >= 0:
            trans = self.listbox.GetClientData(sel_index)
            self.result.remove(trans)
            self.listbox.Delete(sel_index)

    def GetValues(self):
        return self.result, self.check_short.Value


class DialogFast(wx.Dialog):
    """Creates the dialog box Fast Geometric model parameters."""
    def __init__(self, prefix, nf, parent=None):
        super(DialogFast, self).__init__(parent, style=wx.SYSTEM_MENU |
                          wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN)
        self.nf = nf
        self.init_ui()
        self.SetTitle(prefix + ": Fast geometric model (fgm)")

    def init_ui(self):
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        #title
        label_main = wx.StaticText(self, label="Calculation of iTj")
        #input
        grid = wx.GridBagSizer(hgap=25, vgap=5)
        lab_left = wx.StaticText(self, label='Frame i')
        lab_right = wx.StaticText(self, label='Frame j')
        grid.Add(lab_left, pos=(0, 0), flag=wx.ALIGN_CENTER_HORIZONTAL)
        grid.Add(lab_right, pos=(0, 1), flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.cmb_i = wx.ComboBox(self, size=(50, -1),
             choices=[str(i) for i in range(self.nf)], style=wx.CB_READONLY)
        self.cmb_i.SetSelection(0)
        self.cmb_j = wx.ComboBox(self, size=(50, -1),
             choices=[str(i) for i in range(self.nf)], style=wx.CB_READONLY)
        self.cmb_j.SetSelection(self.nf - 1)
        grid.Add(self.cmb_i, pos=(1, 0),
                 flag=wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, border=20)
        grid.Add(self.cmb_j, pos=(1, 1),
                 flag=wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, border=20)
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        grid.Add(cancel_btn, pos=(2, 1))
        grid.Add(ok_btn, pos=(2, 0))
        mainSizer.AddSpacer(30)
        mainSizer.Add(label_main, 0,
                      wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER_HORIZONTAL, 60)
        mainSizer.AddSpacer(30)
        mainSizer.Add(grid, flag=wx.ALIGN_CENTER)
        mainSizer.AddSpacer(20)
        self.SetSizerAndFit(mainSizer)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def GetValues(self):
        return int(self.cmb_i.Value), int(self.cmb_j.Value)


class DialogPaul(wx.Dialog):
    """Creates the dialog box to specify Paul method parameters."""
    def __init__(self, prefix, endeffs, EMPTY, parent=None):
        super(DialogPaul, self).__init__(parent, style=wx.SYSTEM_MENU |
                          wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN)
        self.endeffs = endeffs
        self.init_ui(EMPTY)
        self.SetTitle(prefix + ": IGM Paul Method (pau)")

    def init_ui(self, EMPTY):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        #title
        label_cmb = wx.StaticText(self, label="For frame :")
        main_sizer.Add(label_cmb, 0, wx.TOP | wx.ALIGN_CENTER_HORIZONTAL, 20)
        self.cmb = wx.ComboBox(self, size=(80, -1),
             choices=[str(i) for i in self.endeffs], style=wx.CB_READONLY)
        self.cmb.SetSelection(0)
        main_sizer.AddSpacer(5)
        main_sizer.Add(self.cmb, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 10)
        lbl = wx.StaticText(self, label="Components taken into account :")
        main_sizer.Add(lbl, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 20)
        #input
        grid = wx.GridBagSizer(hgap=15, vgap=15)
        names = ['S', 'N', 'A', 'P']
        for i, name in enumerate(names):
            check_box = wx.CheckBox(self, wx.ID_ANY,
                                    label='   ' + name, name=name)
            check_box.SetValue(True)
            check_box.Bind(wx.EVT_CHECKBOX, self.OnVectorChecked)
            grid.Add(check_box, pos=(0, i), flag=wx.ALIGN_CENTER_HORIZONTAL)
            for j in range(1, 4):
                w_name = name + str(j)
                cmb = wx.ComboBox(self,
                                  choices=[EMPTY, '-1', '0', '1', w_name],
                                  name=w_name, style=wx.CB_READONLY,
                                  size=(90, -1), id=(j-1)*4 + i)
                cmb.SetSelection(4)
                cmb.Bind(wx.EVT_COMBOBOX, self.OnComboBox)
                grid.Add(cmb, pos=(j, i))
            label = wx.StaticText(self,
                                  label=(' 1' if i == 3 else ' 0'), id=12 + i)
            grid.Add(label, pos=(4, i))
        main_sizer.Add(grid, 0, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, 35)
        main_sizer.AddSpacer(20)
        #buttons
        hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ok_btn = wx.Button(self, wx.ID_OK, "OK")
        ok_btn.Bind(wx.EVT_BUTTON, self.OnOK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.OnCancel)
        hor_sizer.Add(ok_btn, 0, wx.ALL, 15)
        hor_sizer.Add(cancel_btn, 0, wx.ALL, 15)
        main_sizer.Add(hor_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.SetSizerAndFit(main_sizer)

    def OnOK(self, _):
        self.EndModal(wx.ID_OK)

    def OnCancel(self, _):
        self.EndModal(wx.ID_CANCEL)

    def OnVectorChecked(self, evt):
        name = evt.EventObject.Name
        index = 4 if evt.EventObject.Value else 0
        for i in range(1, 4):
            cmb = self.FindWindowByName(name + str(i))
            cmb.SetSelection(index)

    def OnComboBox(self, evt):
        name = evt.EventObject.Name
        if evt.EventObject.GetSelection() != 4:
            check_box = self.FindWindowByName(name[0])
            check_box.SetValue(False)

    def get_values(self):
        lst = []
        for i in range(16):
            widget = self.FindWindowById(i)
            if isinstance(widget, wx.ComboBox):
                lst.append(widget.Value)
            else:
                lst.append(widget.LabelText)
        return lst, int(self.cmb.Value)


