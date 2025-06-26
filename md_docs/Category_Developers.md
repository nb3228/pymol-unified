# Category: Developers

## Extensions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Extensions&oldid=12070](https://pymolwiki.org/index.php?title=Extensions&oldid=12070)"


---

## PluginArchitecture

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes PyMOL's plugin architechture since version 1.5.0.6. 

## Contents

  * 1 API Entry Points
  * 2 __init_plugin__
    * 2.1 Example for a PyQt5 GUI (PyMOL 2.x)
    * 2.2 Example for a legacy Tkinter GUI (PyMOL 1.x)
  * 3 More than one file per plugin
  * 4 Config files
  * 5 PyMOL OS Fellowship Project 2011/2012
  * 6 See Also



## API Entry Points

A plugin is a Python module which uses PyMOL's API. The following entrypoints add functionality to PyMOL: 

  1. [`pymol.cmd.extend(func)`](/index.php/Extend "Extend") registers a function as a PyMOL command
  2. `pymol.plugins.addmenuitemqt(label, callback)` adds a plugin menu item, should be called inside `__init_plugin__`
  3. `MyPlugin.__init_plugin__(app)` is called during plugin initialization (if defined by the plugin)



## __init_plugin__

A plugin can implement an `__init_plugin__(app)` function (_previously`__init__`, deprecated_) which is called during plugin initialization. It is passed an object which is compatible with the legacy `PMGApp` instance and provides access to the Tkinter parent (`app.root`). 

### Example for a PyQt5 GUI (PyMOL 2.x)

_See also the[Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")_
    
    
    def __init_plugin__(app=None):
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('My Qt Plugin', myqtdialog)
    
    def myqtdialog():
        from pymol.Qt import QtWidgets
        QtWidgets.QMessageBox.information(None, 'My Plugin Title', 'Hello World')
    

### Example for a legacy Tkinter GUI (PyMOL 1.x)

It is important that all Tkinter objects are created with the `app.root` parent, otherwise legacy plugins won't work in PyMOL 2.0. 
    
    
    def __init_plugin__(app):
        app.menuBar.addmenuitem('Plugin', 'command',
            label='My Tk Plugin',
            command=lambda: mytkdialog(app.root))
    
    def mytkdialog(parent):
        try:
            import tkMessageBox  # Python 2
        except ImportError:
            import tkinter.messagebox as tkMessageBox  # Python 3
        tkMessageBox.showinfo(parent=parent, title='My Plugin Title', message='Hello World')
    

## More than one file per plugin

Plugins can be either a single Python file, or a directory with an `__init__.py` file. 

**Single file layout:**
    
    
    MyPlugin.py (defines "__init_plugin__(app=None)" function)
    

**Directory layout:**
    
    
    MyPlugin/
    ├── data
    │   ├── datasheet.txt
    │   └── image.png
    ├── submodule1.py
    ├── submodule2.py
    └── __init__.py (defines "__init_plugin__(app=None)" function)
    

They can be zipped (`.zip`) for distribution, the [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager") handles the unzipping during installation. 
    
    
    zip -r MyPlugin-1.0.zip MyPlugin
    

## Config files

PyMOL offers a basic API to store custom settings in `~/.pymolpluginsrc.py`. Only basic types (`str`, `int`, `float`, `list`, `tuple`, `dict`, etc.) are suppored (those who produce executable code with `repr()`). 

**Store:**
    
    
    pymol.plugins.pref_set(key, value)
    
    # needed if pymol.plugins.pref_get("instantsave") == False
    pymol.plugins.pref_save()
    

**Load:**
    
    
    value = pymol.plugins.pref_get(key, default=None)
    

**Example:**
    
    
    apbsbinary = pymol.plugins.pref_get("APBS_BINARY_LOCATION", "/usr/bin/apbs")
    

## PyMOL OS Fellowship Project 2011/2012

[Thomas Holder](/index.php/User:Speleo3 "User:Speleo3") was working on the new plugin system as part of his [2011-2012 Fellowship](http://pymol.org/fellowship). The improvements were incorporated into PyMOL 1.5.0.5. 

  * graphical [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * multiple files per plugin/module (see #More than one file per plugin)
  * user plugin directory
  * rarely used plugins can be disabled for better performance and enabled on demand
  * metadata support (version, author, description, citation, tags, ...)
  * install from online repositories
  * settings API for plugins (see #Config files)



## See Also

  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")
  * [Plugins](/index.php/Plugins "Plugins")



Retrieved from "[https://pymolwiki.org/index.php?title=PluginArchitecture&oldid=12829](https://pymolwiki.org/index.php?title=PluginArchitecture&oldid=12829)"


---

## Plugins Tutorial

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This tutorial is about writing your own plugin for PyMOL 2.x. 

See the [plugins](/index.php/Plugins "Plugins") page for how to install and use exisiting plugins. 

## Contents

  * 1 Writing Plugins: Learn By Example
    * 1.1 Plugin Files
    * 1.2 Registering your Plugin
    * 1.3 Creating a GUI
    * 1.4 Filling the GUI with Widgets
    * 1.5 Make Buttons do something
    * 1.6 Deploy the final plugin
    * 1.7 Full Source
  * 2 Extending Plugins to the Command Line
  * 3 See Also



## Writing Plugins: Learn By Example

This tutorial shows how to write a PyMOL plugin with PyQt. **The full source of the demo plugin is[available on github](https://github.com/Pymol-Scripts/pymol2-demo-plugin)**. 

The demo plugin adds a dialog to render images at a custom resolution. 

### Plugin Files

We will create a plugin which consists of [multiple files inside a directory](/index.php/PluginArchitecture#More_than_one_file_per_plugin "PluginArchitecture"): 
    
    
    pymol2-demo-plugin/
    ├── demowidget.ui
    └── __init__.py (defines "__init_plugin__(app=None)" function)
    

### Registering your Plugin

First you must add your plugin to the _Plugins_ menu. This is done in the `__init_plugin__` function of your plugin. A callback (here: `run_plugin_gui`) is added. 
    
    
    # file pymol2-demo-plugin/__init__.py
    def __init_plugin__(app=None):
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Demo "Render" Plugin', run_plugin_gui)
    

_Legacy note_ : The `__init_plugin__` function takes one argument, a reference to the main Tk app to [support legacy plugins witten with Tkinter](/index.php/PluginArchitecture#init_plugin "PluginArchitecture") (unused with PyQt plugins). 

### Creating a GUI

The callback may do arbitrary stuff. Here we're going to create a dialog and show it to the user. 
    
    
    # global reference to avoid garbage collection of our dialog
    dialog = None
    
    def run_plugin_gui():
        # pymol.Qt provides the PyQt5 interface, but may support PyQt4 and/or PySide as well
        from pymol.Qt import QtWidgets
    
        global dialog
    
        if dialog is None:
            # create a new (empty) Window
            dialog = QtWidgets.QDialog()
    
            # TODO: FILL DIALOG WITH WIDGETS HERE
    
        dialog.show()
    

### Filling the GUI with Widgets

[![](/images/8/81/Designer-demo-plugin.png)](/index.php/File:Designer-demo-plugin.png)

[](/index.php/File:Designer-demo-plugin.png "Enlarge")

Screenshot of Qt Designer

The most convenient and maintainable way to create user interfaces with Qt is by using the [Qt Designer](http://doc.qt.io/qt-5/qtdesigner-manual.html). Follow these steps to get started: 

  1. Create a new "Widget" form (New Form > templates/forms > Widget > Create)
  2. Drag widgets like "Push Button" or "Line Edit" into your form
  3. Name your widgets ("objectName" in the "Property Editor"), this name will become a Python variable in your code
  4. Save as a UI file (`*.ui`) inside the `pymol2-demo-plugin` directory



PyMOL provides a utility function to load a UI file into a parent widget: `pymol.Qt.utils.loadUi(filename, widget)`
    
    
        # filename of our UI file
        uifile = os.path.join(os.path.dirname(__file__), 'demowidget.ui')
    
        # load the UI file into our dialog
        from pymol.Qt.utils import loadUi
        form = loadUi(uifile, dialog)
    

### Make Buttons do something

We need to connect the [signals](http://doc.qt.io/qt-5/signalsandslots.html) of those widgets we created with the Qt Designer, like our buttons' `clicked` signal. Our example has a "Ray" button (_objectName=button_ray_) that we'll connect to the `run()` function. 
    
    
        def run():
            # get form data
            height = form.input_height.value()
    
            # some debugging feedback
            print('User entered height', height)
    
            # TODO: DO SOMETHING WITH FORM DATA
    
        form.button_ray.clicked.connect(run)
    

### Deploy the final plugin

The `pymol2-demo-plugin` directory can be zipped for deployment (see [PluginArchitecture#More than one file per plugin](/index.php/PluginArchitecture#More_than_one_file_per_plugin "PluginArchitecture")). 

### Full Source

The full source of the demo plugin is [available on github](https://github.com/Pymol-Scripts/pymol2-demo-plugin)**.**

## Extending Plugins to the Command Line

While not really applicable to our simple "Render" plugin (it only wraps existing `ray` and `png` commands), most plugins should also provide their functionality as new PyMOL commands using [cmd.extend()](/index.php/Extend "Extend"). 

## See Also

  * [PluginArchitecture](/index.php/PluginArchitecture "PluginArchitecture")
  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins](/index.php/Plugins "Plugins")
  * [cmd.extend()](/index.php/Extend "Extend")



Retrieved from "[https://pymolwiki.org/index.php?title=Plugins_Tutorial&oldid=13173](https://pymolwiki.org/index.php?title=Plugins_Tutorial&oldid=13173)"


---

