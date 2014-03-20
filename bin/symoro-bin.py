#!/usr/bin/env python
# -*- coding: utf-8 -*-


__author__ = 'Izzat'


import os

import wx

from pysymoro.symoro import Robot, FAIL
from pysymoro import geometry, kinematics, dynamics, invgeom
from pysymoro.parfile import readpar, writepar
from symoroviz import graphics
from symoroui import ui_definition, ui_geometry, ui_kinematics


PROG_NAME = 'OpenSYMORO'


class MainFrame(wx.Frame):
    def __init__(self):
        title = PROG_NAME + " - SYmbolic MOdeling of RObots"
        size = wx.Size(-1, -1)
        style = wx.DEFAULT_FRAME_STYLE ^ wx.MAXIMIZE_BOX ^ wx.RESIZE_BORDER
        wx.Frame.__init__(self, None, title=title, size=size, style=style)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        self.create_mnu()
        self.robo = Robot.RX90()
        self.widgets = {}
        self.par_dict = {}

        self.statusbar = self.CreateStatusBar()
        self.p = wx.Panel(self)
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)

        self.create_ui()
        self.p.SetSizerAndFit(self.mainSizer)
        self.Fit()
        self.feed_data()

    def params_grid(self, grid, rows, elems, handler=None, start_i=0, size=60):
        for i, name in enumerate(elems):
            horBox = wx.BoxSizer(wx.HORIZONTAL)
            horBox.Add(wx.StaticText(self.p, label=name,
                                     size=(40, -1), style=wx.ALIGN_RIGHT),
                       0, wx.ALL | wx.ALIGN_RIGHT, 5)
            textBox = wx.TextCtrl(self.p, size=(size, -1), name=name, id=i%rows)
            self.widgets[name] = textBox
            textBox.Bind(wx.EVT_KILL_FOCUS, handler)
            horBox.Add(textBox, 0, wx.ALL | wx.ALIGN_LEFT, 1)
            grid.Add(horBox, pos=(i/rows, i % rows + start_i),
                     flag=wx.ALL, border=2)

    def create_ui(self):
        descr_sizer = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Robot Description'), wx.HORIZONTAL)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        descr_sizer.AddSpacer(3)
        descr_sizer.Add(sizer1, 0, wx.ALL | wx.EXPAND, 5)
        sizer2 = wx.BoxSizer(wx.VERTICAL)
        descr_sizer.Add(sizer2, 0, wx.ALL | wx.EXPAND, 5)
        self.mainSizer.Add(descr_sizer, 0, wx.ALL, 10)

        # Left Side of main window
        robot_type_sizer = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Robot Type'), wx.HORIZONTAL)
        grid = wx.GridBagSizer(12, 10)

        gen_info = [('Name of the robot:', 'name'),
                    ('Number of moving links:', 'NL'),
                    ('Number of joints:', 'NJ'), ('Number of frames:', 'NF'),
                    ('Type of structure:', 'type'), ('Is Mobile:', 'mobile'),
                    ('Number of closed loops:', 'loops')]

        for i, (lab, name) in enumerate(gen_info):
            label = wx.StaticText(self.p, label=lab)
            grid.Add(label, pos=(i, 0), flag=wx.LEFT, border=10)
            label = wx.StaticText(self.p, size=(125, -1), name=name)
            self.widgets[name] = label
            grid.Add(label, pos=(i, 1), flag=wx.LEFT | wx.RIGHT, border=10)

        robot_type_sizer.Add(grid, 0, wx.TOP | wx.BOTTOM | wx.EXPAND, 6)
        sizer1.Add(robot_type_sizer)

        ##### Gravity components
        sizer1.AddSpacer(8)
        sbs_gravity = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Gravity components'), wx.HORIZONTAL)
        for iden, name in enumerate(['GX', 'GY', 'GZ']):
            sbs_gravity.AddSpacer(5)
            sbs_gravity.Add(wx.StaticText(self.p, label=name), 0, wx.ALL, 4)
            text_box = wx.TextCtrl(self.p, name=name, size=(60, -1), id=iden)
            self.widgets[name] = text_box
            text_box.Bind(wx.EVT_KILL_FOCUS, self.OnBaseTwistChanged)
            sbs_gravity.Add(text_box, 0, wx.ALL | wx.ALIGN_LEFT, 2)
        sizer1.Add(sbs_gravity, 0, wx.ALL | wx.EXPAND, 0)

        ##### Location of the robot
        sizer1.AddSpacer(8)
        sbs_location = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Location of the robot'), wx.HORIZONTAL)
        sizer1.Add(sbs_location, 0, wx.ALL | wx.EXPAND, 0)
        loc_tbl = wx.GridBagSizer(1, 1)
        for i in range(4):
            ver_lbl = wx.StaticText(self.p, label='Z' + str(i + 1) + ':  ')
            hor_lbl = wx.StaticText(self.p, label='Z - ' + str(i + 1))
            botLab = wx.StaticText(self.p, label='  ' + str(0 if i < 3 else 1))
            loc_tbl.Add(ver_lbl, pos=(i + 1, 0), flag=wx.RIGHT, border=3)
            loc_tbl.Add(hor_lbl, pos=(0, i + 1),
                         flag=wx.ALIGN_CENTER_HORIZONTAL, border=3)
            loc_tbl.Add(botLab, pos=(4, i+1), flag=wx.ALIGN_LEFT, border=3)
            for j in range(3):
                index = j*4+i
                name = 'Z'+str(index)
                text_box = wx.TextCtrl(self.p, name=name,
                                      size=(60, -1), id=index)
                self.widgets[name] = text_box
                text_box.Bind(wx.EVT_KILL_FOCUS, self.OnZParamChanged)
                loc_tbl.Add(text_box, pos=(j + 1, i + 1),
                             flag=wx.ALIGN_LEFT, border=5)
        sbs_location.Add(loc_tbl, 0, wx.ALL | wx.EXPAND, 5)

        ##### Geometric Params
        sbs = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Geometric Params'), wx.HORIZONTAL)
        grid = wx.GridBagSizer(0, 5)
        ver_sizer = wx.BoxSizer(wx.VERTICAL)
        ver_sizer.Add(wx.StaticText(self.p, label='Frame'),
                    0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        combo_box = wx.ComboBox(self.p, style=wx.CB_READONLY,
                               size=(60, -1), name='frame')
        self.widgets['frame'] = combo_box
        combo_box.Bind(wx.EVT_COMBOBOX, self.OnFrameChanged)
        ver_sizer.Add(combo_box)
        for i, name in enumerate(self.robo.get_geom_head()[1:4]):
            hor_box = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(self.p, label=name, size=(40, -1),
                                  style=wx.ALIGN_RIGHT)
            self.widgets[name] = label
            hor_box.Add(label, 0, wx.ALL | wx.ALIGN_RIGHT, 5)
            combo_box = wx.ComboBox(self.p, style=wx.CB_READONLY,
                                   size=(60, -1), name=name)
            self.widgets[name] = combo_box
            combo_box.Bind(wx.EVT_COMBOBOX, self.OnGeoParamChanged)
            hor_box.Add(combo_box, 0, wx.ALL | wx.ALIGN_LEFT, 1)
            grid.Add(hor_box, pos=(i, 1), flag=wx.ALL, border=2)

        grid.Add(ver_sizer, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                 span=(3, 1), border=5)

        self.params_grid(grid, 2, self.robo.get_geom_head()[4:],
                                   self.OnGeoParamChanged, 2, 111)

        sbs.Add(grid)
        sizer2.Add(sbs, 0, wx.ALL | wx.EXPAND, 0)

        ##### Dynamic Params and external forces
        lbl = 'Dynamic Params and external forces'
        sbs = wx.StaticBoxSizer(wx.StaticBox(self.p, label=lbl), wx.HORIZONTAL)
        grid = wx.GridBagSizer(0, 0)
        ver_sizer = wx.BoxSizer(wx.VERTICAL)
        ver_sizer.Add(wx.StaticText(self.p, label='Link'),
                    0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        combo_box = wx.ComboBox(self.p, style=wx.CB_READONLY,
                               size=(60, -1), name='link')
        combo_box.Bind(wx.EVT_COMBOBOX, self.OnLinkChanged)
        self.widgets['link'] = combo_box
        ver_sizer.Add(combo_box)
        grid.Add(ver_sizer, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                 span=(5, 1), border=5)

        params = self.robo.get_dynam_head()[1:]
        params += self.robo.get_ext_dynam_head()[1:-3]
        self.params_grid(grid, 4, params, self.OnDynParamChanged, 1)

        sbs.Add(grid)
        sbs.AddSpacer(4)
        sizer2.AddSpacer(8)
        sizer2.Add(sbs, 0, wx.ALL | wx.EXPAND, 0)

        ##### Speed and acceleration of the base
        lbl = 'Speed and acceleration of the base'
        sbs = wx.StaticBoxSizer(wx.StaticBox(self.p, label=lbl), wx.HORIZONTAL)
        grid = wx.GridBagSizer(0, 0)
        sbs.Add(grid)
        hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        hor_sizer.Add(sbs)
        hor_sizer.AddSpacer(8)

        params = []
        for name in self.robo.get_base_vel_head()[1:-1]:
            for c in ['X', 'Y', 'Z']:
                params.append(name + c)

        self.params_grid(grid, 3, params, self.OnBaseTwistChanged, 0)

        ##### Joints velocity and acceleration
        lbl = 'Joint velocity and acceleration'
        sbs2 = wx.StaticBoxSizer(wx.StaticBox(self.p, label=lbl), wx.VERTICAL)
        sbs2.AddSpacer(5)
        grid = wx.GridBagSizer(5, 5)
        sbs2.Add(grid)
        combo_box = wx.ComboBox(self.p, size=(70, -1),
                               style=wx.CB_READONLY, name='joint')
        self.widgets['joint'] = combo_box
        combo_box.Bind(wx.EVT_COMBOBOX, self.OnJointChanged)
        grid.Add(combo_box, pos=(0, 1),
                 flag=wx.ALIGN_CENTER_HORIZONTAL, border=0)
        names = ['Joint'] + self.robo.get_ext_dynam_head()[-3:]
        for i, name in enumerate(names):
            label = wx.StaticText(self.p, label=name,
                                  size=(55, -1), style=wx.ALIGN_RIGHT)
            grid.Add(label, pos=(i, 0), flag=wx.TOP | wx.RIGHT, border=3)
            if i > 0:
                text_box = wx.TextCtrl(self.p, name=name, size=(70, -1))
                self.widgets[name] = text_box
                text_box.Bind(wx.EVT_KILL_FOCUS, self.OnSpeedChanged)
                grid.Add(text_box, pos=(i, 1))

        hor_sizer.Add(sbs2, 1, wx.ALL | wx.EXPAND, 0)
        sizer2.AddSpacer(8)
        sizer2.Add(hor_sizer, 1, wx.ALL | wx.EXPAND, 0)
        self.mainSizer.AddSpacer(10)
        #sizer2.Add(sbs, 0, wx.ALL | wx.EXPAND, 0)

    def Change(self, index, name, evtObject):
        prev_value = str(self.robo.get_val(index, name))
        if evtObject.Value != prev_value:
            if self.robo.put_val(index, name, evtObject.Value) == FAIL:
                message = "Unacceptable value '%s' has been input in %s%s" \
                          % (evtObject.Value, name, index)
                self.message_error(message)
                evtObject.Value = prev_value
            else:
                self.changed = True

    def OnGeoParamChanged(self, evt):
        frame_index = int(self.widgets['frame'].Value)
        self.Change(frame_index, evt.EventObject.Name, evt.EventObject)
        if evt.EventObject.Name == 'ant':
            self.widgets['type'].SetLabel(self.robo.structure)
        if evt.EventObject.Name == 'sigma':
            self.update_geo_params()

    def OnDynParamChanged(self, evt):
        link_index = int(self.widgets['link'].Value)
        self.Change(link_index, evt.EventObject.Name, evt.EventObject)
        # print type(self.robo.get_val(link_index, evt.EventObject.Name))

    def OnSpeedChanged(self, evt):
        joint_index = int(self.widgets['joint'].Value)
        self.Change(joint_index, evt.EventObject.Name, evt.EventObject)

    def OnBaseTwistChanged(self, evt):
        index = int(evt.EventObject.Id)
        name = evt.EventObject.Name[:-1]
        self.Change(index, name, evt.EventObject)

    def OnZParamChanged(self, evt):
        index = int(evt.EventObject.Id)
        self.Change(index, 'Z', evt.EventObject)

    def OnFrameChanged(self, evt):
        frame_index = int(evt.EventObject.Value)
        cmb = self.widgets['ant']
        cmb.SetItems([str(i) for i in range(frame_index)])
        self.update_geo_params()

    def OnLinkChanged(self, _):
        self.update_dyn_params()

    def OnJointChanged(self, _):
        self.update_vel_params()

    def update_params(self, index, pars):
        for par in pars:
            widget = self.widgets[par]
            widget.ChangeValue(str(self.robo.get_val(index, par)))

    def update_geo_params(self):
        index = int(self.widgets['frame'].Value)
        for par in self.robo.get_geom_head()[1:4]:
            self.widgets[par].SetValue(str(self.robo.get_val(index, par)))
        self.update_params(index, self.robo.get_geom_head()[4:])

    def update_dyn_params(self):
        pars = self.robo.get_dynam_head()[1:]
        # cut first and last 3 elements
        pars += self.robo.get_ext_dynam_head()[1:-3]

        index = int(self.widgets['link'].Value)
        self.update_params(index, pars)

    def update_vel_params(self):
        pars = self.robo.get_ext_dynam_head()[-3:]
        index = int(self.widgets['joint'].Value)
        self.update_params(index, pars)

    def update_base_twist_params(self):
        for name in self.robo.get_base_vel_head()[1:]:
            for i, c in enumerate(['X', 'Y', 'Z']):
                widget = self.widgets[name + c]
                widget.ChangeValue(str(self.robo.get_val(i, name)))

    def update_z_params(self):
        T = self.robo.Z
        for i in range(12):
            widget = self.widgets['Z' + str(i)]
            widget.ChangeValue(str(T[i]))

    def feed_data(self):
        # Robot Type
        names = [('name', self.robo.name), ('NF', self.robo.nf),
                 ('NL', self.robo.nl), ('NJ', self.robo.nj),
                 ('type', self.robo.structure),
                 ('mobile', self.robo.is_mobile),
                 ('loops', self.robo.nj-self.robo.nl)]
        for name, info in names:
            label = self.widgets[name]
            label.SetLabel(str(info))

        lsts = [('frame', [str(i) for i in range(1, self.robo.NF)]),
                ('link',  [str(i) for i in range(int(not self.robo.is_mobile),
                                                 self.robo.NL)]),
                ('joint', [str(i) for i in range(1, self.robo.NJ)]),
                ('ant', ['0']), ('sigma', ['0', '1', '2']), ('mu', ['0', '1'])]
        for name, lst in lsts:
            cmb = self.widgets[name]
            cmb.SetItems(lst)
            cmb.SetSelection(0)

        self.update_geo_params()
        self.update_dyn_params()
        self.update_vel_params()
        self.update_base_twist_params()
        self.update_z_params()

        self.changed = False
        self.par_dict = {}

    def create_mnu(self):
        mnu_bar = wx.MenuBar()

        ##### FILE
        file_mnu = wx.Menu()
        m_new = file_mnu.Append(wx.ID_NEW, '&New...')
        self.Bind(wx.EVT_MENU, self.OnNew, m_new)
        m_open = file_mnu.Append(wx.ID_OPEN, '&Open...')
        self.Bind(wx.EVT_MENU, self.OnOpen, m_open)
        m_save = file_mnu.Append(wx.ID_SAVE, '&Save...')
        self.Bind(wx.EVT_MENU, self.OnSave, m_save)
        m_save_as = file_mnu.Append(wx.ID_SAVEAS, '&Save As...')
        self.Bind(wx.EVT_MENU, self.OnSaveAs, m_save_as)
        file_mnu.Append(wx.ID_PREFERENCES, '&Preferences...')
        file_mnu.AppendSeparator()

        m_exit = file_mnu.Append(wx.ID_EXIT, "E&xit\tAlt-X",
                                 "Close window and exit program.")
        self.Bind(wx.EVT_MENU, self.OnClose, m_exit)
        mnu_bar.Append(file_mnu, "&File")

        ##### GEOMETRIC
        geo_mnu = wx.Menu()
        m_trans_matrix = wx.MenuItem(geo_mnu, wx.ID_ANY,
                                     "Transformation matrix...")
        self.Bind(wx.EVT_MENU, self.OnTransformationMatrix, m_trans_matrix)
        geo_mnu.AppendItem(m_trans_matrix)
        fast_dgm = wx.MenuItem(geo_mnu, wx.ID_ANY, "Fast geometric model...")
        self.Bind(wx.EVT_MENU, self.OnFastGeometricModel, fast_dgm)
        geo_mnu.AppendItem(fast_dgm)
        igm_paul = wx.MenuItem(geo_mnu, wx.ID_ANY, "I.G.M. Paul method...")
        self.Bind(wx.EVT_MENU, self.OnIGMPaul, igm_paul)
        geo_mnu.AppendItem(igm_paul)
        constr_geom = wx.MenuItem(geo_mnu, wx.ID_ANY,
                                  "Constraint geometric equations of loops")
        self.Bind(wx.EVT_MENU, self.OnConstraintGeoEq, constr_geom)
        geo_mnu.AppendItem(constr_geom)

        mnu_bar.Append(geo_mnu, "&Geometric")

        ##### KINEMATIC
        kin_mnu = wx.Menu()
        jac_matrix = wx.MenuItem(kin_mnu, wx.ID_ANY, "Jacobian matrix...")
        self.Bind(wx.EVT_MENU, self.OnJacobianMatrix, jac_matrix)
        kin_mnu.AppendItem(jac_matrix)
        determ = wx.MenuItem(kin_mnu, wx.ID_ANY, "Determinant of a Jacobian...")
        self.Bind(wx.EVT_MENU, self.OnDeterminant, determ)
        kin_mnu.AppendItem(determ)
        #TODO: add the dialog, ask for projection frame
        ckel = wx.MenuItem(kin_mnu, wx.ID_ANY, "Kinematic constraints")
        self.Bind(wx.EVT_MENU, self.OnCkel, ckel)
        kin_mnu.AppendItem(ckel)
        vels = wx.MenuItem(kin_mnu, wx.ID_ANY, "Velocities")
        self.Bind(wx.EVT_MENU, self.OnVelocities, vels)
        kin_mnu.AppendItem(vels)
        accel = wx.MenuItem(kin_mnu, wx.ID_ANY, "Accelerations")
        self.Bind(wx.EVT_MENU, self.OnAccelerations, accel)
        kin_mnu.AppendItem(accel)
        jpqp = wx.MenuItem(kin_mnu, wx.ID_ANY, "Jpqp")
        self.Bind(wx.EVT_MENU, self.OnJpqp, jpqp)
        kin_mnu.AppendItem(jpqp)

        mnu_bar.Append(kin_mnu, "&Kinematic")

        ##### DYNAMIC
        dyn_mnu = wx.Menu()
        dyn_submnu = wx.MenuItem(dyn_mnu, wx.ID_ANY, 'Inverse dynamic model')
        self.Bind(wx.EVT_MENU, self.OnInverseDynamic, dyn_submnu)
        dyn_mnu.AppendItem(dyn_submnu)
        dyn_submnu = wx.MenuItem(dyn_mnu, wx.ID_ANY, 'Inertia Matrix')
        self.Bind(wx.EVT_MENU, self.OnInertiaMatrix, dyn_submnu)
        dyn_mnu.AppendItem(dyn_submnu)
        dyn_submnu = wx.MenuItem(dyn_mnu, wx.ID_ANY,
                                 'Centrifugal, Coriolis & Gravity torques')
        self.Bind(wx.EVT_MENU, self.OnCentrCoriolGravTorq, dyn_submnu)
        dyn_mnu.AppendItem(dyn_submnu)
        dyn_submnu = wx.MenuItem(dyn_mnu, wx.ID_ANY, 'Direct Dynamic Model')
        self.Bind(wx.EVT_MENU, self.OnDirectDynamicModel, dyn_submnu)
        dyn_mnu.AppendItem(dyn_submnu)

        mnu_bar.Append(dyn_mnu, "&Dynamic")

        ##### IDENTIFICATION
        iden = wx.Menu()
        base_inert = wx.MenuItem(
            iden, wx.ID_ANY, 'Base inertial parameters (symbolic or numeric)')
        self.Bind(wx.EVT_MENU, self.OnBaseInertialParams, base_inert)
        iden.AppendItem(base_inert)
        dyn_iden_model = wx.MenuItem(iden, wx.ID_ANY,
                                     'Dynamic identification model')
        self.Bind(wx.EVT_MENU, self.OnDynIdentifModel, dyn_iden_model)
        iden.AppendItem(dyn_iden_model)

        mnu_bar.Append(iden, "&Identification")

        ##### OPTIMIZER
        # optMenu = wx.Menu()
        # jac_matrix = wx.MenuItem(optMenu, wx.ID_ANY, "Jacobian matrix...")
        # self.Bind(wx.EVT_MENU, self.OnJacobianMatrix, jac_matrix)
        # optMenu.AppendItem(jac_matrix)
        #
        # menuBar.Append(optMenu, "&Optimizer")

        ##### VISUALIZATION
        vis_mnu = wx.Menu()
        vis_menu = wx.MenuItem(vis_mnu, wx.ID_ANY, "Visualisation")
        self.Bind(wx.EVT_MENU, self.OnVisualisation, vis_menu)
        vis_mnu.AppendItem(vis_menu)

        mnu_bar.Append(vis_mnu, "&Visualization")

        self.SetMenuBar(mnu_bar)

    def OnNew(self, _):
        dialog = ui_definition.DialogDefinition(
            PROG_NAME, self.robo.name, self.robo.nl,
            self.robo.nj, self.robo.structure, self.robo.is_mobile)
        if dialog.ShowModal() == wx.ID_OK:
            result = dialog.get_values()
            new_robo = Robot(*result['init_pars'])
            if result['keep_geo']:
                nf = min(self.robo.NF, new_robo.NF)
                new_robo.ant[:nf] = self.robo.ant[:nf]
                new_robo.sigma[:nf] = self.robo.sigma[:nf]
                new_robo.mu[:nf] = self.robo.mu[:nf]
                new_robo.gamma[:nf] = self.robo.gamma[:nf]
                new_robo.alpha[:nf] = self.robo.alpha[:nf]
                new_robo.theta[:nf] = self.robo.theta[:nf]
                new_robo.b[:nf] = self.robo.b[:nf]
                new_robo.d[:nf] = self.robo.d[:nf]
                new_robo.r[:nf] = self.robo.r[:nf]
            if result['keep_dyn']:
                nl = min(self.robo.NL, new_robo.NL)
                new_robo.Nex[:nl] = self.robo.Nex[:nl]
                new_robo.Fex[:nl] = self.robo.Fex[:nl]
                new_robo.FS[:nl] = self.robo.FS[:nl]
                new_robo.IA[:nl] = self.robo.IA[:nl]
                new_robo.FV[:nl] = self.robo.FV[:nl]
                new_robo.MS[:nl] = self.robo.MS[:nl]
                new_robo.M[:nl] = self.robo.M[:nl]
                new_robo.J[:nl] = self.robo.J[:nl]
            if result['keep_base']:
                new_robo.Z = self.robo.Z
                new_robo.w0 = self.robo.w0
                new_robo.wdot0 = self.robo.wdot0
                new_robo.v0 = self.robo.v0
                new_robo.vdot0 = self.robo.vdot0
                new_robo.G = self.robo.G
            self.robo = new_robo
            directory = os.path.join('robots', self.robo.name)
            if not os.path.exists(directory):
                os.makedirs(directory)
            self.robo.directory = directory
            self.feed_data()
        dialog.Destroy()

    def message_error(self, message):
        wx.MessageDialog(None, message,
                         'Error', wx.OK | wx.ICON_ERROR).ShowModal()

    def message_warning(self, message):
        wx.MessageDialog(None, message,
                         'Error', wx.OK | wx.ICON_WARNING).ShowModal()

    def message_info(self, message):
        wx.MessageDialog(None, message, 'Information',
                         wx.OK | wx.ICON_INFORMATION).ShowModal()

    def model_success(self, model_name):
        msg = 'The model has been saved in %s\\%s_%s.txt' % \
              (self.robo.directory, self.robo.name, model_name)
        self.message_info(msg)

    def OnOpen(self, _):
        if self.changed:
            dialog_res = wx.MessageBox('Do you want to save changes?',
                                       'Please confirm',
                                       wx.ICON_QUESTION |
                                       wx.YES_NO | wx.CANCEL,
                                       self)
            if dialog_res == wx.CANCEL:
                return
            elif dialog_res == wx.YES:
                if self.OnSave(None) == FAIL:
                    return
        dialog = wx.FileDialog(self, message="Choose PAR file", style=wx.OPEN,
                               wildcard='*.par', defaultFile='*.par')
        if dialog.ShowModal() == wx.ID_OK:
            new_robo, flag = readpar(dialog.GetDirectory(),
                                     dialog.GetFilename()[:-4])
            if new_robo is None:
                self.message_error('File could not be read!')
            else:
                if flag == FAIL:
                    self.message_warning('While reading file an error occured.')
                self.robo = new_robo
                self.feed_data()

    def OnSave(self, _):
        writepar(self.robo)
        self.changed = False

    def OnSaveAs(self, _):
        dialog = wx.FileDialog(self, message="Save PAR file",
                               defaultFile=self.robo.name+'.par',
                               defaultDir=self.robo.directory,
                               wildcard='*.par')
        if dialog.ShowModal() == wx.ID_CANCEL:
            return FAIL

        self.robo.directory = dialog.GetDirectory()
        self.robo.name = dialog.GetFilename()[:-4]
        writepar(self.robo)
        self.widgets['name'].SetLabel(self.robo.name)
        self.changed = False

    def OnTransformationMatrix(self, _):
        dialog = ui_geometry.DialogTrans(PROG_NAME, self.robo.NF)
        if dialog.ShowModal() == wx.ID_OK:
            frames, trig_subs = dialog.GetValues()
            geometry.direct_geometric(self.robo, frames, trig_subs)
            self.model_success('trm')
        dialog.Destroy()

    def OnFastGeometricModel(self, _):
        dialog = ui_geometry.DialogFast(PROG_NAME, self.robo.NF)
        if dialog.ShowModal() == wx.ID_OK:
            i, j = dialog.GetValues()
            geometry.direct_geometric_fast(self.robo, i, j)
            self.model_success('fgm')
        dialog.Destroy()

    def OnIGMPaul(self, _):
        dialog = ui_geometry.DialogPaul(PROG_NAME, self.robo.endeffectors,
                                        str(invgeom.EMPTY))
        if dialog.ShowModal() == wx.ID_OK:
            lst_T, n = dialog.get_values()
            invgeom.igm_Paul(self.robo, lst_T, n)
            self.model_success('igm')
        dialog.Destroy()

    def OnConstraintGeoEq(self, _):
        pass

    def OnJacobianMatrix(self, _):
        dialog = ui_kinematics.DialogJacobian(PROG_NAME, self.robo)
        if dialog.ShowModal() == wx.ID_OK:
            n, i, j = dialog.get_values()
            kinematics.jacobian(self.robo, n, i, j)
            self.model_success('jac')
        dialog.Destroy()

    def OnDeterminant(self, _):
        dialog = ui_kinematics.DialogDeterminant(PROG_NAME, self.robo)
        if dialog.ShowModal() == wx.ID_OK:
            kinematics.jacobian_determinant(self.robo, *dialog.get_values())
            self.model_success('det')
        dialog.Destroy()

    def OnCkel(self, _):
        if kinematics.kinematic_constraints(self.robo) == FAIL:
            self.message_warning('There are no loops')
        else:
            self.model_success('ckel')

    def OnVelocities(self, _):
        kinematics.velocities(self.robo)
        self.model_success('vlct')

    def OnAccelerations(self, _):
        kinematics.accelerations(self.robo)
        self.model_success('aclr')

    def OnJpqp(self, _):
        kinematics.jdot_qdot(self.robo)
        self.model_success('jpqp')

    def OnInverseDynamic(self, _):
        dynamics.inverse_dynamic_NE(self.robo)
        self.model_success('idm')

    def OnInertiaMatrix(self, _):
        dynamics.inertia_matrix(self.robo)
        self.model_success('inm')

    def OnCentrCoriolGravTorq(self, _):
        dynamics.pseudo_force_NE(self.robo)
        self.model_success('ccg')

    def OnDirectDynamicModel(self, _):
        dynamics.direct_dynamic_NE(self.robo)
        self.model_success('ddm')

    def OnBaseInertialParams(self, _):
        dynamics.base_paremeters(self.robo)
        self.model_success('regp')

    def OnDynIdentifModel(self, _):
        dynamics.dynamic_identification_NE(self.robo)
        self.model_success('dim')

    def OnVisualisation(self, _):
        dialog = ui_definition.DialogConversion(PROG_NAME,
                                                self.robo, self.par_dict)
        if dialog.has_syms():
            if dialog.ShowModal() == wx.ID_OK:
                self.par_dict = dialog.get_values()
                graphics.MainWindow(PROG_NAME, self.robo, self.par_dict, self)
        else:
            graphics.MainWindow(PROG_NAME, self.robo, self.par_dict, self)

    def OnClose(self, _):
        if self.changed:
            result = wx.MessageBox(
                'Do you want to save changes?', 'Please confirm',
                wx.ICON_QUESTION | wx.YES_NO | wx.CANCEL, self)
            if result == wx.YES:
                if self.OnSave(_) == FAIL:
                    return
            elif result == wx.CANCEL:
                return
        self.Destroy()


def main():
    app = wx.App(redirect=False)
    frame = MainFrame()
    frame.Show()
    app.MainLoop()


if __name__ == "__main__":
    main()


