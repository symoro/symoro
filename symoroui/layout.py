#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module creates the main window user interface and draws the
interface on the screen for the SYMORO package.
"""


import os
from collections import OrderedDict

import wx

from pysymoro.robot import Robot
from pysymoro import geometry
from pysymoro import kinematics
from pysymoro import dynamics
from pysymoro import invgeom
from symoroutils import parfile
from symoroutils import filemgr
from symoroutils import samplerobots
from symoroutils import tools
from symoroui import definition as ui_definition
from symoroui import geometry as ui_geometry
from symoroui import kinematics as ui_kinematics
from symoroui import labels as ui_labels
from symoroviz import graphics


class MainFrame(wx.Frame):
    """This Frame contains the main window for SYMORO"""
    def __init__(self, *args, **kwargs):
        """Constructor : creates the UI and draws it on the screen."""
        wx.Frame.__init__(self, *args, **kwargs)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        # create status bar
        self.statusbar = self.CreateStatusBar()
        # create menu bar
        self.create_menu()
        # set default robot
        self.robo = samplerobots.planar2r()
        # object to store different ui elements
        self.widgets = {}
        # object to store parameter values got from dialog box input
        self.par_dict = {}
        # setup panel and sizer for content
        self.panel = wx.Panel(self)
        self.szr_topmost = wx.BoxSizer(wx.VERTICAL)
        self.create_ui()
        self.panel.SetSizerAndFit(self.szr_topmost)
        self.Fit()
        # update fields with data
        self.feed_data()
        # configure status bar
        self.statusbar.SetFieldsCount(number=2)
        self.statusbar.SetStatusWidths(widths=[-1, -1])
        self.statusbar.SetStatusText(text="Ready", number=0)
        self.statusbar.SetStatusText(
            text="Location of robot files is %s"
            % filemgr.get_base_path(), number = 1
        )

    def params_in_grid(self, szr_grd, elements, rows, cols, width=70):
        """Method to display a set of fields in a grid."""
        for idx, key in enumerate(elements):
            label = elements[key].label
            name = elements[key].name
            control = elements[key].control
            place = elements[key].place
            handler = getattr(self, elements[key].handler)
            field_id = int(elements[key].id)
            if control is 'cmb':
                ctrl = wx.ComboBox(
                    parent=self.panel, style=wx.CB_READONLY,
                    size=(width, -1), name=name
                )
                ctrl.Bind(wx.EVT_COMBOBOX, handler)
            elif control is 'lbl':
                ctrl = wx.StaticText(
                    parent=self.panel, size=(width, -1), name=name
                )
            else:
                ctrl = wx.TextCtrl(
                    parent=self.panel, size=(width, -1),
                    name=name, id=field_id
                )
                ctrl.Bind(wx.EVT_KILL_FOCUS, handler)
            self.widgets[name] = ctrl
            szr_ele = wx.BoxSizer(wx.HORIZONTAL)
            szr_ele.Add(
                wx.StaticText(
                    self.panel, label=label, style=wx.ALIGN_RIGHT
                ), proportion=0, flag=wx.ALL | wx.ALIGN_RIGHT, border=5
            )
            szr_ele.Add(
                ctrl, proportion=0,
                flag=wx.ALL | wx.ALIGN_LEFT, border=1
            )
            szr_grd.Add(
                szr_ele, pos=(place[0], place[1]),
                flag=wx.ALL | wx.ALIGN_RIGHT, border=2
            )

    def create_ui(self):
        """Method to create the contents of the user interface"""
        # main box - robot description box
        szr_robot_des = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['robot_des']
            ), wx.HORIZONTAL
        )
        szr_robot_des.AddSpacer(3)
        szr_left_col = wx.BoxSizer(wx.VERTICAL)
        szr_right_col = wx.BoxSizer(wx.VERTICAL)
        szr_robot_des.Add(szr_left_col, 0, wx.ALL | wx.EXPAND, 5)
        szr_robot_des.Add(szr_right_col, 0, wx.ALL | wx.EXPAND, 5)
        self.szr_topmost.Add(szr_robot_des, 0, wx.ALL, 10)
        # left col - robot type box
        szr_robot_type = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['robot_type']
            ), wx.HORIZONTAL
        )
        szr_grd_robot_type = wx.GridBagSizer(10, 12)
        for idx, key in enumerate(ui_labels.ROBOT_TYPE):
            label = ui_labels.ROBOT_TYPE[key].label
            name = ui_labels.ROBOT_TYPE[key].name
            self.widgets[name] = wx.StaticText(
                self.panel, size=(150, -1),
                name=ui_labels.ROBOT_TYPE[key].name
            )
            szr_grd_robot_type.Add(
                wx.StaticText(self.panel, label=label),
                pos=(idx, 0), flag=wx.LEFT, border=10
            )
            szr_grd_robot_type.Add(
                self.widgets[name], pos=(idx, 1),
                flag=wx.LEFT | wx.RIGHT, border=10
            )
        szr_robot_type.Add(
            szr_grd_robot_type, 0, wx.TOP | wx.BOTTOM | wx.EXPAND, 6
        )
        szr_left_col.Add(szr_robot_type)
        szr_left_col.AddSpacer(8)
        # left col - gravity components box
        szr_gravity = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['gravity']
            ), wx.HORIZONTAL
        )
        szr_grd_gravity = wx.GridBagSizer(5, 5)
        self.params_in_grid(
            szr_grd_gravity, elements=ui_labels.GRAVITY_CMPNTS,
            rows=1, cols=3, width=70,
        )
        szr_gravity.Add(szr_grd_gravity)
        szr_left_col.Add(szr_gravity, 0, wx.ALL | wx.EXPAND, 0)
        szr_left_col.AddSpacer(8)
        # left col - location of the robot box
        # NOTE: the columns of the T matrix are filled first
        szr_location = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['location']
            ), wx.HORIZONTAL
        )
        szr_grd_loc = wx.GridBagSizer(1, 1)
        for i in range(4):
            for j in range(3):
                idx = (j*4) + i
                name = 'Z'+str(idx)
                txt_z_element = wx.TextCtrl(
                    parent=self.panel, name=name,
                    id=idx, size=(60, -1)
                )
                self.widgets[name] = txt_z_element
                txt_z_element.Bind(
                    wx.EVT_KILL_FOCUS, self.OnZParamChanged
                )
                szr_grd_loc.Add(
                    txt_z_element, pos=(j + 1, i + 1),
                    flag=wx.ALIGN_LEFT, border=5
                )
            lbl_row = wx.StaticText(self.panel, label='Z'+str(i + 1))
            lbl_col = wx.StaticText(self.panel, label='Z'+str(i + 1))
            lbl_last_row = wx.StaticText(
                self.panel, label='  '+str(0 if i < 3 else 1)
            )
            szr_grd_loc.Add(
                lbl_row, pos=(i+1, 0), flag=wx.RIGHT, border=3
            )
            szr_grd_loc.Add(
                lbl_col, pos=(0, i+1),
                flag=wx.ALIGN_CENTER_HORIZONTAL, border=3
            )
            szr_grd_loc.Add(
                lbl_last_row, pos=(4, i+1),
                flag=wx.ALIGN_LEFT, border=3
            )
        szr_location.Add(szr_grd_loc, 0, wx.ALL | wx.EXPAND, 5)
        szr_left_col.Add(szr_location, 0, wx.ALL | wx.EXPAND, 0)
        # right col - geometric params box
        szr_geom_params = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['geom_params']
            ), wx.HORIZONTAL
        )
        szr_grd_geom = wx.GridBagSizer(0, 5)
        self.params_in_grid(
            szr_grd_geom, elements=ui_labels.GEOM_PARAMS,
            rows=2, cols=5, width=70
        )
        szr_geom_params.Add(szr_grd_geom)
        szr_right_col.Add(szr_geom_params, 0, wx.ALL | wx.EXPAND, 0)
        szr_right_col.AddSpacer(8)
        # right col - dynamic params and external forces box
        szr_dyn_params = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['dyn_params']
            ), wx.VERTICAL
        )
        cmb_link = wx.ComboBox(
            self.panel, style=wx.CB_READONLY, size=(100, -1),
            name=ui_labels.DYN_PARAMS['link'].name
        )
        cmb_link.Bind(
            wx.EVT_COMBOBOX,
            getattr(self, ui_labels.DYN_PARAMS['link'].handler)
        )
        self.widgets['link'] = cmb_link
        szr_link = wx.BoxSizer(wx.HORIZONTAL)
        szr_link.Add(
            wx.StaticText(
                self.panel, label=ui_labels.DYN_PARAMS['link'].label
            ), proportion=0,
            flag=wx.ALL | wx.ALIGN_LEFT, border=5
        )
        szr_link.AddSpacer((4,4))
        szr_link.Add(cmb_link, flag=wx.ALL | wx.ALIGN_RIGHT)
        szr_dyn_params.Add(szr_link, flag=wx.ALL | wx.ALIGN_CENTER)
        szr_grd_dyn = wx.GridBagSizer(0, 0)
        # add dynamic params to the grid
        elements = OrderedDict(ui_labels.DYN_PARAMS.items()[1:])
        self.params_in_grid(
            szr_grd_dyn, elements=elements, rows=4, cols=6, width=75
        )
        szr_dyn_params.Add(szr_grd_dyn)
        szr_dyn_params.AddSpacer(4)
        szr_right_col.Add(szr_dyn_params, 0, wx.ALL | wx.EXPAND, 0)
        szr_right_col.AddSpacer(8)
        # box sizer for the last row in right col
        szr_velacc = wx.BoxSizer(wx.HORIZONTAL)
        # right col - velocity and acceleration of the base box
        szr_base_velacc = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['base_vel_acc']
            ), wx.HORIZONTAL
        )
        szr_grd_base_velacc = wx.GridBagSizer(0, 0)
        self.params_in_grid(
            szr_grd_base_velacc, elements=ui_labels.BASE_VEL_ACC,
            rows=3, cols=4, width=60
        )
        szr_base_velacc.Add(szr_grd_base_velacc)
        szr_velacc.Add(szr_base_velacc)
        szr_velacc.AddSpacer(8)
        # right col - joint velocity and acceleration box
        szr_joint_velacc = wx.StaticBoxSizer(
            wx.StaticBox(
                self.panel, label=ui_labels.BOX_TITLES['joint_vel_acc']
            ), wx.HORIZONTAL
        )
        szr_grd_joint_velacc = wx.GridBagSizer(5, 5)
        self.params_in_grid(
            szr_grd_joint_velacc, elements=ui_labels.JOINT_VEL_ACC,
            rows=3, cols=1, width=75,
        )
        szr_joint_velacc.Add(szr_grd_joint_velacc,
            flag=wx.ALL | wx.ALIGN_CENTER, border=2
        )
        szr_velacc.Add(szr_joint_velacc, 1, wx.ALL | wx.EXPAND, 0)
        szr_right_col.Add(szr_velacc, 1, wx.ALL | wx.EXPAND, 0)
        self.szr_topmost.AddSpacer(10)

    def Change(self, index, name, event_object):
        prev_value = str(self.robo.get_val(index, name))
        if event_object.Value != prev_value:
            if self.robo.put_val(index, name, event_object.Value) == tools.FAIL:
                message = "Unacceptable value '%s' has been input in %s%s" \
                          % (event_object.Value, name, index)
                self.message_error(message)
                event_object.Value = prev_value
            else:
                self.changed = True

    def OnGeoParamChanged(self, event):
        frame_index = int(self.widgets['frame'].Value)
        self.Change(frame_index, event.EventObject.Name, event.EventObject)
        if event.EventObject.Name == 'ant':
            self.widgets['type'].SetLabel(self.robo.structure)
        if event.EventObject.Name == 'sigma':
            self.update_geo_params()

    def OnDynParamChanged(self, event):
        link_index = int(self.widgets['link'].Value)
        self.Change(link_index, event.EventObject.Name, event.EventObject)

    def OnSpeedChanged(self, event):
        joint_index = int(self.widgets['joint'].Value)
        self.Change(joint_index, event.EventObject.Name, event.EventObject)

    def OnBaseTwistChanged(self, event):
        index = int(event.EventObject.Id)
        name = event.EventObject.Name[:-1]
        self.Change(index, name, event.EventObject)

    def OnZParamChanged(self, event):
        index = int(event.EventObject.Id)
        self.Change(index, 'Z', event.EventObject)

    def OnFrameChanged(self, event):
        frame_index = int(event.EventObject.Value)
        cmb = self.widgets['ant']
        cmb.SetItems([str(i) for i in range(frame_index)])
        self.update_geo_params()

    def OnLinkChanged(self, event):
        self.update_dyn_params()

    def OnJointChanged(self, event):
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

    def create_menu(self):
        """Method to create the menu bar"""
        menu_bar = wx.MenuBar()
        # menu item - file
        file_menu = wx.Menu()
        m_new = wx.MenuItem(
            file_menu, wx.ID_NEW, ui_labels.FILE_MENU['m_new']
        )
        self.Bind(wx.EVT_MENU, self.OnNew, m_new)
        file_menu.AppendItem(m_new)
        m_open = wx.MenuItem(
            file_menu, wx.ID_OPEN, ui_labels.FILE_MENU['m_open']
        )
        self.Bind(wx.EVT_MENU, self.OnOpen, m_open)
        file_menu.AppendItem(m_open)
        m_save = wx.MenuItem(
            file_menu, wx.ID_SAVE, ui_labels.FILE_MENU['m_save']
        )
        self.Bind(wx.EVT_MENU, self.OnSave, m_save)
        file_menu.AppendItem(m_save)
        m_save_as = wx.MenuItem(
            file_menu, wx.ID_SAVEAS, ui_labels.FILE_MENU['m_save_as']
        )
        self.Bind(wx.EVT_MENU, self.OnSaveAs, m_save_as)
        file_menu.AppendItem(m_save_as)
        m_pref = wx.MenuItem(
            file_menu, wx.ID_ANY, ui_labels.FILE_MENU['m_pref']
        )
        file_menu.AppendItem(m_pref)
        file_menu.AppendSeparator()
        m_exit = wx.MenuItem(
            file_menu, wx.ID_EXIT, ui_labels.FILE_MENU['m_exit']
        )
        self.Bind(wx.EVT_MENU, self.OnClose, m_exit)
        file_menu.AppendItem(m_exit)
        menu_bar.Append(file_menu, ui_labels.MAIN_MENU['file_menu'])
        # menu item - geometric
        geom_menu = wx.Menu()
        m_trans_matrix = wx.MenuItem(
            geom_menu, wx.ID_ANY, ui_labels.GEOM_MENU['m_trans_matrix']
        )
        self.Bind(
            wx.EVT_MENU, self.OnTransformationMatrix, m_trans_matrix
        )
        geom_menu.AppendItem(m_trans_matrix)
        m_fast_dgm = wx.MenuItem(
            geom_menu, wx.ID_ANY, ui_labels.GEOM_MENU['m_fast_dgm']
        )
        self.Bind(wx.EVT_MENU, self.OnFastGeometricModel, m_fast_dgm)
        geom_menu.AppendItem(m_fast_dgm)
        m_igm_paul = wx.MenuItem(
            geom_menu, wx.ID_ANY, ui_labels.GEOM_MENU['m_igm_paul']
        )
        self.Bind(wx.EVT_MENU, self.OnIgmPaul, m_igm_paul)
        geom_menu.AppendItem(m_igm_paul)
        m_geom_constraint = wx.MenuItem(
            geom_menu, wx.ID_ANY, ui_labels.GEOM_MENU['m_geom_constraint']
        )
        self.Bind(wx.EVT_MENU, self.OnConstraintGeoEq, m_geom_constraint)
        geom_menu.AppendItem(m_geom_constraint)
        menu_bar.Append(geom_menu, ui_labels.MAIN_MENU['geom_menu'])
        # menu item - kinematic
        kin_menu = wx.Menu()
        m_jac_matrix = wx.MenuItem(
            kin_menu, wx.ID_ANY, ui_labels.KIN_MENU['m_jac_matrix']
        )
        self.Bind(wx.EVT_MENU, self.OnJacobianMatrix, m_jac_matrix)
        kin_menu.AppendItem(m_jac_matrix)
        m_determinant = wx.MenuItem(
            kin_menu, wx.ID_ANY, ui_labels.KIN_MENU['m_determinant']
        )
        self.Bind(wx.EVT_MENU, self.OnDeterminant, m_determinant)
        kin_menu.AppendItem(m_determinant)
        #TODO: add the dialog, ask for projection frame
        m_kin_constraint = wx.MenuItem(
            kin_menu, wx.ID_ANY, ui_labels.KIN_MENU['m_kin_constraint']
        )
        self.Bind(wx.EVT_MENU, self.OnCkel, m_kin_constraint)
        kin_menu.AppendItem(m_kin_constraint)
        m_vel = wx.MenuItem(
            kin_menu, wx.ID_ANY, ui_labels.KIN_MENU['m_vel']
        )
        self.Bind(wx.EVT_MENU, self.OnVelocities, m_vel)
        kin_menu.AppendItem(m_vel)
        m_acc = wx.MenuItem(
            kin_menu, wx.ID_ANY, ui_labels.KIN_MENU['m_acc']
        )
        self.Bind(wx.EVT_MENU, self.OnAccelerations, m_acc)
        kin_menu.AppendItem(m_acc)
        m_jpqp = wx.MenuItem(
            kin_menu, wx.ID_ANY, ui_labels.KIN_MENU['m_jpqp']
        )
        self.Bind(wx.EVT_MENU, self.OnJpqp, m_jpqp)
        kin_menu.AppendItem(m_jpqp)
        menu_bar.Append(kin_menu, ui_labels.MAIN_MENU['kin_menu'])
        # menu item - dynamic
        dyn_menu = wx.Menu()
        m_idym = wx.MenuItem(
            dyn_menu, wx.ID_ANY, ui_labels.DYN_MENU['m_idym']
        )
        self.Bind(wx.EVT_MENU, self.OnInverseDynamic, m_idym)
        dyn_menu.AppendItem(m_idym)
        m_inertia_matrix = wx.MenuItem(
            dyn_menu, wx.ID_ANY, ui_labels.DYN_MENU['m_inertia_matrix']
        )
        self.Bind(wx.EVT_MENU, self.OnInertiaMatrix, m_inertia_matrix)
        dyn_menu.AppendItem(m_inertia_matrix)
        m_h_term = wx.MenuItem(
            dyn_menu, wx.ID_ANY, ui_labels.DYN_MENU['m_h_term']
        )
        self.Bind(wx.EVT_MENU, self.OnCentrCoriolGravTorq, m_h_term)
        dyn_menu.AppendItem(m_h_term)
        m_ddym = wx.MenuItem(
            dyn_menu, wx.ID_ANY, ui_labels.DYN_MENU['m_ddym']
        )
        self.Bind(wx.EVT_MENU, self.OnDirectDynamicModel, m_ddym)
        dyn_menu.AppendItem(m_ddym)
        menu_bar.Append(dyn_menu, ui_labels.MAIN_MENU['dyn_menu'])
        # menu item - identification
        iden_menu = wx.Menu()
        m_base_inertial_params = wx.MenuItem(
            iden_menu, wx.ID_ANY,
            ui_labels.IDEN_MENU['m_base_inertial_params']
        )
        self.Bind(
            wx.EVT_MENU, self.OnBaseInertialParams, m_base_inertial_params
        )
        iden_menu.AppendItem(m_base_inertial_params)
        m_dyn_iden_model = wx.MenuItem(
            iden_menu, wx.ID_ANY, ui_labels.IDEN_MENU['m_dyn_iden_model']
        )
        self.Bind(wx.EVT_MENU, self.OnDynIdentifModel, m_dyn_iden_model)
        iden_menu.AppendItem(m_dyn_iden_model)
        m_energy_iden_model = wx.MenuItem(
            iden_menu, wx.ID_ANY,
            ui_labels.IDEN_MENU['m_energy_iden_model']
        )
        # TODO: uncomment the 3 lines below to add the event
        #self.Bind(
        #    wx.EVT_MENU, self.OnEnergyIdentifModel, m_energy_iden_model
        #)
        iden_menu.AppendItem(m_energy_iden_model)
        menu_bar.Append(iden_menu, ui_labels.MAIN_MENU['iden_menu'])
        # menu item - visualisation
        viz_menu = wx.Menu()
        m_viz = wx.MenuItem(
            viz_menu, wx.ID_ANY, ui_labels.VIZ_MENU['m_viz']
        )
        self.Bind(wx.EVT_MENU, self.OnVisualisation, m_viz)
        viz_menu.AppendItem(m_viz)
        menu_bar.Append(viz_menu, ui_labels.MAIN_MENU['viz_menu'])
        # set menu bar
        self.SetMenuBar(menu_bar)

    def OnNew(self, event):
        dialog = ui_definition.DialogDefinition(
            ui_labels.MAIN_WIN['prog_name'],
            self.robo.name, self.robo.nl,
            self.robo.nj, self.robo.structure, self.robo.is_mobile
        )
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
            self.robo.directory = filemgr.get_folder_path(self.robo.name)
            self.feed_data()
        dialog.Destroy()

    def message_error(self, message):
        wx.MessageDialog(
            None,
            message,
            'Error',
            wx.OK | wx.ICON_ERROR
        ).ShowModal()

    def message_warning(self, message):
        wx.MessageDialog(
            None,
            message,
            'Error',
            wx.OK | wx.ICON_WARNING
        ).ShowModal()

    def message_info(self, message):
        wx.MessageDialog(
            None,
            message,
            'Information',
            wx.OK | wx.ICON_INFORMATION
        ).ShowModal()

    def model_success(self, model_name):
        msg = 'The model has been saved in %s_%s.txt' % (
            os.path.join(self.robo.directory, self.robo.name), model_name
        )
        self.message_info(msg)

    def OnOpen(self, event):
        if self.changed:
            dialog_res = wx.MessageBox(
                'Do you want to save changes?',
                'Please confirm',
                wx.ICON_QUESTION | wx.YES_NO | wx.CANCEL,
                self
            )
            if dialog_res == wx.CANCEL:
                return
            elif dialog_res == wx.YES:
                if self.OnSave(None) == tools.FAIL:
                    return
        dialog = wx.FileDialog(
            self,
            message="Choose PAR file",
            style=wx.OPEN,
            wildcard='*.par',
            defaultFile='*.par'
        )
        if dialog.ShowModal() == wx.ID_OK:
            new_robo, flag = parfile.readpar(
                dialog.GetFilename()[:-4], dialog.GetPath()
            )
            if new_robo is None:
                self.message_error('File could not be read!')
            else:
                if flag == tools.FAIL:
                    self.message_warning(
                        "While reading file an error occured."
                    )
                self.robo = new_robo
                self.feed_data()

    def OnSave(self, event):
        parfile.writepar(self.robo)
        self.changed = False

    def OnSaveAs(self, event):
        dialog = wx.FileDialog(
            self, message="Save PAR file",
            defaultFile=self.robo.name+'.par',
            defaultDir=self.robo.directory,
            wildcard='*.par'
        )
        if dialog.ShowModal() == wx.ID_CANCEL:
            return tools.FAIL
        self.robo.directory = dialog.GetDirectory()
        self.robo.name = dialog.GetFilename()[:-4]
        parfile.writepar(self.robo)
        self.widgets['name'].SetLabel(self.robo.name)
        self.changed = False

    def OnTransformationMatrix(self, event):
        dialog = ui_geometry.DialogTrans(
            ui_labels.MAIN_WIN['prog_name'], self.robo.NF
        )
        if dialog.ShowModal() == wx.ID_OK:
            frames, trig_subs = dialog.GetValues()
            geometry.direct_geometric(self.robo, frames, trig_subs)
            self.model_success('trm')
        dialog.Destroy()

    def OnFastGeometricModel(self, event):
        dialog = ui_geometry.DialogFast(
            ui_labels.MAIN_WIN['prog_name'], self.robo.NF
        )
        if dialog.ShowModal() == wx.ID_OK:
            i, j = dialog.GetValues()
            geometry.direct_geometric_fast(self.robo, i, j)
            self.model_success('fgm')
        dialog.Destroy()

    def OnIgmPaul(self, event):
        dialog = ui_geometry.DialogPaul(
            ui_labels.MAIN_WIN['prog_name'],
            self.robo.endeffectors,
            str(invgeom.EMPTY)
        )
        if dialog.ShowModal() == wx.ID_OK:
            lst_T, n = dialog.get_values()
            invgeom.igm_Paul(self.robo, lst_T, n)
            self.model_success('igm')
        dialog.Destroy()

    def OnConstraintGeoEq(self, event):
        pass

    def OnJacobianMatrix(self, event):
        dialog = ui_kinematics.DialogJacobian(
            ui_labels.MAIN_WIN['prog_name'], self.robo
        )
        if dialog.ShowModal() == wx.ID_OK:
            n, i, j = dialog.get_values()
            kinematics.jacobian(self.robo, n, i, j)
            self.model_success('jac')
        dialog.Destroy()

    def OnDeterminant(self, event):
        dialog = ui_kinematics.DialogDeterminant(
            ui_labels.MAIN_WIN['prog_name'], self.robo
        )
        if dialog.ShowModal() == wx.ID_OK:
            kinematics.jacobian_determinant(self.robo, *dialog.get_values())
            self.model_success('det')
        dialog.Destroy()

    def OnCkel(self, event):
        if kinematics.kinematic_constraints(self.robo) == tools.FAIL:
            self.message_warning('There are no loops')
        else:
            self.model_success('ckel')

    def OnVelocities(self, event):
        kinematics.velocities(self.robo)
        self.model_success('vlct')

    def OnAccelerations(self, event):
        kinematics.accelerations(self.robo)
        self.model_success('aclr')

    def OnJpqp(self, event):
        kinematics.jdot_qdot(self.robo)
        self.model_success('jpqp')

    def OnInverseDynamic(self, event):
        dynamics.inverse_dynamic_NE(self.robo)
        self.model_success('idm')

    def OnInertiaMatrix(self, event):
        dynamics.inertia_matrix(self.robo)
        self.model_success('inm')

    def OnCentrCoriolGravTorq(self, event):
        dynamics.pseudo_force_NE(self.robo)
        self.model_success('ccg')

    def OnDirectDynamicModel(self, event):
        dynamics.direct_dynamic_NE(self.robo)
        self.model_success('ddm')

    def OnBaseInertialParams(self, event):
        dynamics.base_paremeters(self.robo)
        self.model_success('regp')

    def OnDynIdentifModel(self, event):
        dynamics.dynamic_identification_NE(self.robo)
        self.model_success('dim')

    def OnVisualisation(self, event):
        dialog = ui_definition.DialogVisualisation(
            ui_labels.MAIN_WIN['prog_name'], self.robo, self.par_dict
        )
        if dialog.has_syms():
            if dialog.ShowModal() == wx.ID_OK:
                self.par_dict = dialog.get_values()
                graphics.MainWindow(
                    ui_labels.MAIN_WIN['prog_name'], self.robo,
                    self.par_dict, self
                )
        else:
            graphics.MainWindow(
                ui_labels.MAIN_WIN['prog_name'], self.robo,
                self.par_dict, self
            )

    def OnClose(self, event):
        if self.changed:
            result = wx.MessageBox(
                'Do you want to save changes?', 'Please confirm',
                wx.ICON_QUESTION | wx.YES_NO | wx.CANCEL, self)
            if result == wx.YES:
                if self.OnSave(_) == tools.FAIL:
                    return
            elif result == wx.CANCEL:
                return
        self.Destroy()
        wx.GetApp().ExitMainLoop()


