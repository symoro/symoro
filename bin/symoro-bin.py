#!/usr/bin/env python
# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This is the main executable script of SYMORO.
"""


import os
import sys

import wx

from symoroui import layout
from symoroui import labels as ui_labels


def main():
    app = wx.App(redirect=False)
    style = wx.DEFAULT_FRAME_STYLE ^ wx.MAXIMIZE_BOX ^ wx.RESIZE_BORDER
    frame = layout.MainFrame(
        parent=None,
        id=wx.ID_ANY,
        title=ui_labels.MAIN_WIN['window_title'],
        size=(-1, -1),
        style=style
    )
    frame.Show()
    app.MainLoop()


if __name__ == "__main__":
    main()


