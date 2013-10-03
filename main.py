__author__ = 'Izzat'
import wx
from Core.symoro import Robot, FAIL
from Core import geometry, kinematics, dynamics, invgeom
from Graphics import graphics
from UI import ui_definition, ui_geometry, ui_kinematics
from Core.parfile import readpar, writepar
import os

PROG_NAME = 'SYMORO-Python'


class MainFrame(wx.Frame):
    def __init__(self):
        mainTitle = PROG_NAME + ": SYmbolic MOdeling of RObots"
        style = wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER ^ wx.MAXIMIZE_BOX
        wx.Frame.__init__(self, None, title=mainTitle, size=(0, 0), style=style)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        self.CreateMenu()
        self.robo = Robot.RX90()
        self.widgets = {}
        self.par_dict = {}

        self.statusbar = self.CreateStatusBar()
        self.p = wx.Panel(self)
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)

        m_text = wx.StaticText(self.p, -1, "SYmbolic MOdelling of RObots")
        m_text.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text.SetSize(m_text.GetBestSize())
        self.mainSizer.Add(m_text, 0, wx.TOP | wx.ALIGN_CENTER_HORIZONTAL, 15)

        self.CreateUI()
        self.p.SetSizerAndFit(self.mainSizer)
        self.Fit()
        self.FeedData()

    def LabelsTextboxesInGrid(self, grid, rows, elems, handler=None, start_i=0):
        for i, name in enumerate(elems):
            horBox = wx.BoxSizer(wx.HORIZONTAL)
            horBox.Add(wx.StaticText(self.p, label=name,
                                     size=(40, -1), style=wx.ALIGN_RIGHT),
                       0, wx.ALL | wx.ALIGN_RIGHT, 5)
            textBox = wx.TextCtrl(self.p, size=(60, -1), name=name, id=i % rows)
            self.widgets[name] = textBox
            textBox.Bind(wx.EVT_KILL_FOCUS, handler)
            horBox.Add(textBox, 0, wx.ALL | wx.ALIGN_LEFT, 1)
            grid.Add(horBox, pos=(i/rows, i % rows + start_i),
                     flag=wx.ALL, border=2)

    def CreateUI(self):
        robotDescSizer = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Robot Description'), wx.HORIZONTAL)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        robotDescSizer.AddSpacer(3)
        robotDescSizer.Add(sizer1, 0, wx.ALL | wx.EXPAND, 5)
        sizer2 = wx.BoxSizer(wx.VERTICAL)
        robotDescSizer.Add(sizer2, 0, wx.ALL | wx.EXPAND, 5)
        self.mainSizer.Add(robotDescSizer, 0, wx.ALL, 10)

        # Left Side of main window
        robotTypeSizer = wx.StaticBoxSizer(
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

        robotTypeSizer.Add(grid, 0, wx.TOP | wx.BOTTOM | wx.EXPAND, 6)
        sizer1.Add(robotTypeSizer)

        ##### Gravity components
        sizer1.AddSpacer(8)
        sbsGravity = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Gravity components'), wx.HORIZONTAL)
        for iden, name in enumerate(['GX', 'GY', 'GZ']):
            sbsGravity.AddSpacer(5)
            sbsGravity.Add(wx.StaticText(self.p, label=name), 0, wx.ALL, 4)
            textBox = wx.TextCtrl(self.p, name=name, size=(60, -1), id=iden)
            self.widgets[name] = textBox
            textBox.Bind(wx.EVT_KILL_FOCUS, self.OnBaseTwistChanged)
            sbsGravity.Add(textBox, 0, wx.ALL | wx.ALIGN_LEFT, 2)
        sizer1.Add(sbsGravity, 0, wx.ALL | wx.EXPAND, 0)

        ##### Location of the robot
        sizer1.AddSpacer(8)
        sbsLocation = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Location of the robot'), wx.HORIZONTAL)
        sizer1.Add(sbsLocation, 0, wx.ALL | wx.EXPAND, 0)
        locTable = wx.GridBagSizer(1, 1)
        for i in range(4):
            verLab = wx.StaticText(self.p, label='Z' + str(i + 1) + ':  ')
            horLab = wx.StaticText(self.p, label='Z - ' + str(i + 1))
            botLab = wx.StaticText(self.p, label='  ' + str(0 if i < 3 else 1))
            locTable.Add(verLab, pos=(i + 1, 0), flag=wx.RIGHT, border=3)
            locTable.Add(horLab, pos=(0, i + 1),
                         flag=wx.ALIGN_CENTER_HORIZONTAL, border=3)
            locTable.Add(botLab, pos=(4, i+1), flag=wx.ALIGN_LEFT, border=3)
            for j in range(3):
                index = j*4+i
                name = 'Z'+str(index)
                textBox = wx.TextCtrl(self.p, name=name,
                                      size=(60, -1), id=index)
                self.widgets[name] = textBox
                textBox.Bind(wx.EVT_KILL_FOCUS, self.OnZParamChanged)
                locTable.Add(textBox, pos=(j + 1, i + 1),
                             flag=wx.ALIGN_LEFT, border=5)
        sbsLocation.Add(locTable, 0, wx.ALL | wx.EXPAND, 5)

        ##### Geometric Params
        sbs = wx.StaticBoxSizer(
            wx.StaticBox(self.p, label='Geometric Params'), wx.HORIZONTAL)
        grid = wx.GridBagSizer(0, 39)
        vertBox = wx.BoxSizer(wx.VERTICAL)
        vertBox.Add(wx.StaticText(self.p, label='Frame'),
                    0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        comboBox = wx.ComboBox(self.p, style=wx.CB_READONLY,
                               size=(60, -1), name='frame')
        self.widgets['frame'] = comboBox
        comboBox.Bind(wx.EVT_COMBOBOX, self.OnFrameChanged)
        vertBox.Add(comboBox)
        for i, name in enumerate(self.robo.get_geom_head()[1:4]):
            horBox = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(self.p, label=name, size=(40, -1),
                                  style=wx.ALIGN_RIGHT)
            self.widgets[name] = label
            horBox.Add(label, 0, wx.ALL | wx.ALIGN_RIGHT, 5)
            comboBox = wx.ComboBox(self.p, style=wx.CB_READONLY,
                                   size=(60, -1), name=name)
            self.widgets[name] = comboBox
            comboBox.Bind(wx.EVT_COMBOBOX, self.OnGeoParamChanged)
            horBox.Add(comboBox, 0, wx.ALL | wx.ALIGN_LEFT, 1)
            grid.Add(horBox, pos=(i, 1), flag=wx.ALL, border=2)

        grid.Add(vertBox, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                 span=(3, 1), border=5)

        self.LabelsTextboxesInGrid(grid, 2, self.robo.get_geom_head()[4:],
                                   self.OnGeoParamChanged, 2)

        sbs.Add(grid)
        sizer2.Add(sbs, 0, wx.ALL | wx.EXPAND, 0)

        ##### Dynamic Params and external forces
        lbl = 'Dynamic Params and external forces'
        sbs = wx.StaticBoxSizer(wx.StaticBox(self.p, label=lbl), wx.HORIZONTAL)
        grid = wx.GridBagSizer(0, 0)
        vertBox = wx.BoxSizer(wx.VERTICAL)
        vertBox.Add(wx.StaticText(self.p, label='Link'),
                    0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        comboBox = wx.ComboBox(self.p, style=wx.CB_READONLY,
                               size=(60, -1), name='link')
        comboBox.Bind(wx.EVT_COMBOBOX, self.OnLinkChanged)
        self.widgets['link'] = comboBox
        vertBox.Add(comboBox)
        grid.Add(vertBox, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                 span=(5, 1), border=5)

        params = self.robo.get_dynam_head()[1:]
        params += self.robo.get_ext_dynam_head()[1:-3]
        self.LabelsTextboxesInGrid(grid, 4, params, self.OnDynParamChanged, 1)

        sbs.Add(grid)
        sbs.AddSpacer(4)
        sizer2.AddSpacer(8)
        sizer2.Add(sbs, 0, wx.ALL | wx.EXPAND, 0)

        ##### Speed and acceleration of the base
        lbl = 'Speed and acceleration of the base'
        sbs = wx.StaticBoxSizer(wx.StaticBox(self.p, label=lbl), wx.HORIZONTAL)
        grid = wx.GridBagSizer(0, 0)
        sbs.Add(grid)
        horSizer = wx.BoxSizer(wx.HORIZONTAL)
        horSizer.Add(sbs)
        horSizer.AddSpacer(8)

        params = []
        for name in self.robo.get_base_vel_head()[1:-1]:
            for c in ['X', 'Y', 'Z']:
                params.append(name + c)

        self.LabelsTextboxesInGrid(grid, 3, params, self.OnBaseTwistChanged, 0)

        ##### Joints velocity and acceleration
        lbl = 'Joint velocity and acceleration'
        sbs2 = wx.StaticBoxSizer(wx.StaticBox(self.p, label=lbl), wx.VERTICAL)
        sbs2.AddSpacer(5)
        grid = wx.GridBagSizer(5, 5)
        sbs2.Add(grid)
        comboBox = wx.ComboBox(self.p, size=(70, -1),
                               style=wx.CB_READONLY, name='joint')
        self.widgets['joint'] = comboBox
        comboBox.Bind(wx.EVT_COMBOBOX, self.OnJointChanged)
        grid.Add(comboBox, pos=(0, 1),
                 flag=wx.ALIGN_CENTER_HORIZONTAL, border=0)
        names = ['Joint'] + self.robo.get_ext_dynam_head()[-3:]
        for i, name in enumerate(names):
            label = wx.StaticText(self.p, label=name,
                                  size=(55, -1), style=wx.ALIGN_RIGHT)
            grid.Add(label, pos=(i, 0), flag=wx.TOP | wx.RIGHT, border=3)
            if i > 0:
                textBox = wx.TextCtrl(self.p, name=name, size=(70, -1))
                self.widgets[name] = textBox
                textBox.Bind(wx.EVT_KILL_FOCUS, self.OnSpeedChanged)
                grid.Add(textBox, pos=(i, 1))

        horSizer.Add(sbs2, 1, wx.ALL | wx.EXPAND, 0)
        sizer2.AddSpacer(8)
        sizer2.Add(horSizer, 1, wx.ALL | wx.EXPAND, 0)
        self.mainSizer.AddSpacer(10)
        #sizer2.Add(sbs, 0, wx.ALL | wx.EXPAND, 0)

    def Change(self, index, name, evtObject):
        prev_value = str(self.robo.get_val(index, name))
        if evtObject.Value != prev_value:
            if self.robo.put_val(index, name, evtObject.Value) == FAIL:
                message = "Unacceptable value '%s' has been input in %s%s" \
                          % (evtObject.Value, name, index)
                self.MessageError(message)
                evtObject.Value = prev_value
            else:
                self.changed = True

    def OnGeoParamChanged(self, evt):
        frame_index = int(self.widgets['frame'].Value)
        self.Change(frame_index, evt.EventObject.Name, evt.EventObject)
        if evt.EventObject.Name == 'ant':
            self.widgets['type'].SetLabel(self.robo.structure)
        if evt.EventObject.Name == 'sigma':
            self.UpdateGeoParams()

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
        comboBox = self.widgets['ant']
        comboBox.SetItems([str(i) for i in range(frame_index)])
        self.UpdateGeoParams()

    def OnLinkChanged(self, _):
        self.UpdateDynParams()

    def OnJointChanged(self, _):
        self.UpdateVelParams()

    def UpdateParams(self, index, pars):
        for par in pars:
            widget = self.widgets[par]
            widget.ChangeValue(str(self.robo.get_val(index, par)))

    def UpdateGeoParams(self):
        index = int(self.widgets['frame'].Value)
        for par in self.robo.get_geom_head()[1:4]:
            self.widgets[par].SetValue(str(self.robo.get_val(index, par)))
        self.UpdateParams(index, self.robo.get_geom_head()[4:])

    def UpdateDynParams(self):
        pars = self.robo.get_dynam_head()[1:]
        # cut first and last 3 elements
        pars += self.robo.get_ext_dynam_head()[1:-3]

        index = int(self.widgets['link'].Value)
        self.UpdateParams(index, pars)

    def UpdateVelParams(self):
        pars = self.robo.get_ext_dynam_head()[-3:]
        index = int(self.widgets['joint'].Value)
        self.UpdateParams(index, pars)

    def UpdateBaseTwistParams(self):
        for name in self.robo.get_base_vel_head()[1:]:
            for i, c in enumerate(['X', 'Y', 'Z']):
                widget = self.widgets[name + c]
                widget.ChangeValue(str(self.robo.get_val(i, name)))

    def UpdateZParams(self):
        T = self.robo.Z
        for i in range(12):
            widget = self.widgets['Z' + str(i)]
            widget.ChangeValue(str(T[i]))

    def FeedData(self):
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
            comboBox = self.widgets[name]
            comboBox.SetItems(lst)
            comboBox.SetSelection(0)

        self.UpdateGeoParams()
        self.UpdateDynParams()
        self.UpdateVelParams()
        self.UpdateBaseTwistParams()
        self.UpdateZParams()

        self.changed = False
        self.par_dict = {}

    def CreateMenu(self):
        menuBar = wx.MenuBar()

        ##### FILE
        fileMenu = wx.Menu()
        m_new = fileMenu.Append(wx.ID_NEW, '&New...')
        self.Bind(wx.EVT_MENU, self.OnNew, m_new)
        m_open = fileMenu.Append(wx.ID_OPEN, '&Open...')
        self.Bind(wx.EVT_MENU, self.OnOpen, m_open)
        m_save = fileMenu.Append(wx.ID_SAVE, '&Save...')
        self.Bind(wx.EVT_MENU, self.OnSave, m_save)
        m_save_as = fileMenu.Append(wx.ID_SAVEAS, '&Save As...')
        self.Bind(wx.EVT_MENU, self.OnSaveAs, m_save_as)
        fileMenu.Append(wx.ID_PREFERENCES, '&Preferences...')
        fileMenu.AppendSeparator()

        m_exit = fileMenu.Append(wx.ID_EXIT, "E&xit\tAlt-X",
                                 "Close window and exit program.")
        self.Bind(wx.EVT_MENU, self.OnClose, m_exit)
        menuBar.Append(fileMenu, "&File")

        ##### GEOMETRIC
        geoMenu = wx.Menu()
        m_trans_matrix = wx.MenuItem(geoMenu, wx.ID_ANY,
                                     "Transformation matrix...")
        self.Bind(wx.EVT_MENU, self.OnTransformationMatrix, m_trans_matrix)
        geoMenu.AppendItem(m_trans_matrix)
        fast_dgm = wx.MenuItem(geoMenu, wx.ID_ANY, "Fast geometric model...")
        self.Bind(wx.EVT_MENU, self.OnFastGeometricModel, fast_dgm)
        geoMenu.AppendItem(fast_dgm)
        igm_paul = wx.MenuItem(geoMenu, wx.ID_ANY, "I.G.M. Paul method...")
        self.Bind(wx.EVT_MENU, self.OnIGMPaul, igm_paul)
        geoMenu.AppendItem(igm_paul)
        constr_geom = wx.MenuItem(geoMenu, wx.ID_ANY,
                                  "Constraint geometric equations of loops")
        self.Bind(wx.EVT_MENU, self.OnConstraintGeoEq, constr_geom)
        geoMenu.AppendItem(constr_geom)

        menuBar.Append(geoMenu, "&Geometric")

        ##### KINEMATIC
        kinMenu = wx.Menu()
        jac_matrix = wx.MenuItem(kinMenu, wx.ID_ANY, "Jacobian matrix...")
        self.Bind(wx.EVT_MENU, self.OnJacobianMatrix, jac_matrix)
        kinMenu.AppendItem(jac_matrix)
        determ = wx.MenuItem(kinMenu, wx.ID_ANY, "Determinant of a Jacobian...")
        self.Bind(wx.EVT_MENU, self.OnDeterminant, determ)
        kinMenu.AppendItem(determ)

        menuBar.Append(kinMenu, "&Kinematic")

        ##### DYNAMIC
        dynMenu = wx.Menu()
        dynSubMenu = wx.MenuItem(dynMenu, wx.ID_ANY, 'Inverse dynamic model')
        self.Bind(wx.EVT_MENU, self.OnInverseDynamic, dynSubMenu)
        dynMenu.AppendItem(dynSubMenu)
        dynSubMenu = wx.MenuItem(dynMenu, wx.ID_ANY, 'Inertia Matrix')
        self.Bind(wx.EVT_MENU, self.OnInertiaMatrix, dynSubMenu)
        dynMenu.AppendItem(dynSubMenu)
        dynSubMenu = wx.MenuItem(dynMenu, wx.ID_ANY,
                                 'Centrifugal, Coriolis & Gravity torques')
        self.Bind(wx.EVT_MENU, self.OnCentrCoriolGravTorq, dynSubMenu)
        dynMenu.AppendItem(dynSubMenu)
        dynSubMenu = wx.MenuItem(dynMenu, wx.ID_ANY, 'Direct Dynamic Model')
        self.Bind(wx.EVT_MENU, self.OnDirectDynamicModel, dynSubMenu)
        dynMenu.AppendItem(dynSubMenu)

        menuBar.Append(dynMenu, "&Dynamic")

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

        menuBar.Append(iden, "&Identification")

        ##### OPTIMIZER
        # optMenu = wx.Menu()
        # jac_matrix = wx.MenuItem(optMenu, wx.ID_ANY, "Jacobian matrix...")
        # self.Bind(wx.EVT_MENU, self.OnJacobianMatrix, jac_matrix)
        # optMenu.AppendItem(jac_matrix)
        #
        # menuBar.Append(optMenu, "&Optimizer")

        ##### VISUALIZATION
        visMenu = wx.Menu()
        vis_menu = wx.MenuItem(visMenu, wx.ID_ANY, "Visualisation")
        self.Bind(wx.EVT_MENU, self.OnVisualisation, vis_menu)
        visMenu.AppendItem(vis_menu)

        menuBar.Append(visMenu, "&Visualization")

        self.SetMenuBar(menuBar)

    def OnNew(self, _):
        dialog = ui_definition.DialogDefinition(
            PROG_NAME, self.robo.name, self.robo.nl,
            self.robo.nj, self.robo.structure, self.robo.is_mobile)
        if dialog.ShowModal() == wx.ID_OK:
            result = dialog.GetValues()
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
            directory = 'Robots\\' + self.robo.name + '\\'
            if not os.path.exists(directory):
                os.makedirs(directory)
            self.robo.directory = directory
            self.FeedData()
        dialog.Destroy()

    def MessageError(self, message):
        wx.MessageDialog(None, message,
                         'Error', wx.OK | wx.ICON_ERROR).ShowModal()

    def MessageWarning(self, message):
        wx.MessageDialog(None, message,
                         'Error', wx.OK | wx.ICON_WARNING).ShowModal()

    def MessageInfo(self, message):
        wx.MessageDialog(None, message, 'Information',
                         wx.OK | wx.ICON_INFORMATION).ShowModal()

    def model_success(self, model_name):
        msg = 'The model has been saved in models\\%s_%s.txt' % \
              (self.robo.name, model_name)
        self.MessageInfo(msg)

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
                self.MessageError('File could not be read!')
            else:
                if flag == FAIL:
                    self.MessageWarning('While reading file an error occured.')
                self.robo = new_robo
                self.FeedData()

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
            lst_T, n = dialog.GetValues()
            invgeom.igm_Paul(self.robo, lst_T, n)
            self.model_success('igm')
        dialog.Destroy()

    def OnConstraintGeoEq(self, _):
        pass

    def OnJacobianMatrix(self, _):
        dialog = ui_kinematics.DialogJacobian(PROG_NAME, self.robo)
        if dialog.ShowModal() == wx.ID_OK:
            n, i, j = dialog.GetValues()
            kinematics.jacobian(self.robo, n, i, j)
            self.model_success('jac')
        dialog.Destroy()

    def OnDeterminant(self, _):
        dialog = ui_kinematics.DialogDeterminant(PROG_NAME, self.robo)
        if dialog.ShowModal() == wx.ID_OK:
            kinematics.jacobian_determinant(self.robo, *dialog.GetValues())
            self.model_success('det')
        dialog.Destroy()

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
        if dialog.HasSyms():
            if dialog.ShowModal() == wx.ID_OK:
                self.par_dict = dialog.GetValues()
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

app = wx.App(redirect=False)
main = MainFrame()
main.Show()
app.MainLoop()