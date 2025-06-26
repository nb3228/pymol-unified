# Category: GUI

## Internal gui

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

The internal GUI is the openGL-based user interface (that is--not the TK-window). Here you can enable/disable proteins from view, color then, change their representation--a whole host of things. If you want a nice clean presentation mode, you can turn off the internal GUI using this setting. 

  * [![Internal GUI on](/images/9/98/Ig_on.png)](/index.php/File:Ig_on.png "Internal GUI on")

Internal GUI on 

  * [![Internal GUI off](/images/f/f3/Ig_off.png)](/index.php/File:Ig_off.png "Internal GUI off")

Internal GUI off 




# Syntax
    
    
    # turn the internal gui off
    set internal_gui, off
    
    # turn it back on
    set internal_gui, on
    

# See Also

[GUI Category](/index.php/Category:GUI "Category:GUI") , [Group](/index.php/Group "Group")

Retrieved from "[https://pymolwiki.org/index.php?title=Internal_gui&oldid=10532](https://pymolwiki.org/index.php?title=Internal_gui&oldid=10532)"


---

## Internal gui control size

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This sets how large the internal gui controls are. An internal GUI control is the button you click on to enable/disable a protein, select colors, etc. 

  * [![Set to default, 18.](/images/c/cd/InteranlGuiControlSize_Default.png)](/index.php/File:InteranlGuiControlSize_Default.png "Set to default, 18.")

Set to default, 18. 

  * [![Set to 30](/images/4/46/Igcs_df30.png)](/index.php/File:Igcs_df30.png "Set to 30")

Set to 30 




# Syntax
    
    
    # set the size to some positive integer
    set internal_gui_control_size, int
    
    # for example,
    set internal_gui_control_size, 30
    

# See Also

[Get](/index.php/Get "Get"), [Internal_gui](/index.php/Internal_gui "Internal gui")

Retrieved from "[https://pymolwiki.org/index.php?title=Internal_gui_control_size&oldid=6164](https://pymolwiki.org/index.php?title=Internal_gui_control_size&oldid=6164)"


---

## Internal gui mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This sets how PyMOL renders the internal GUI & controls. The default setting, 0. 

  * [![internal_gui_mode, set to 0](/images/0/05/Ig_mode0.png)](/index.php/File:Ig_mode0.png "internal_gui_mode, set to 0")

internal_gui_mode, set to 0 

  * [![internal_gui_mode, set to 1. Notice the coloring under the controls.](/images/6/65/Ig_mode1.png)](/index.php/File:Ig_mode1.png "internal_gui_mode, set to 1. Notice the coloring under the controls.")

internal_gui_mode, set to 1. Notice the coloring under the controls. 




# Syntax
    
    
    # set to 0 or 1
    set internal_gui_mode, int
    
    # for example, put it in mode 1
    set internal_gui_mode, 1
    

# See Also

[GUI Category](/index.php/Category:GUI "Category:GUI"), [Internal_gui](/index.php/Internal_gui "Internal gui")

Retrieved from "[https://pymolwiki.org/index.php?title=Internal_gui_mode&oldid=6165](https://pymolwiki.org/index.php?title=Internal_gui_mode&oldid=6165)"


---

## Internal gui width

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting determines with width of the openGL-based internal GUI in pixels. If you have objects with very long names, this might come in handy so you can see more of the names. 

In the following images, notice the size of the GUI controls on the right: 

  * [![internal_gui_width set to the default, 220.](/images/0/0d/Igwidth220.png)](/index.php/File:Igwidth220.png "internal_gui_width set to the default, 220.")

internal_gui_width set to the default, 220. 

  * [![internal_gui_width set to 50.](/images/3/3a/Igwidth50.png)](/index.php/File:Igwidth50.png "internal_gui_width set to 50.")

internal_gui_width set to 50. 

  * [![internal_gui_width set to the default, 500.](/images/a/a9/Igwidth500.png)](/index.php/File:Igwidth500.png "internal_gui_width set to the default, 500.")

internal_gui_width set to the default, 500. 




# Syntax
    
    
    # set to a positive integer width (in pixels)
    set internal_gui_width, int
    
    # for example; set the width to 446 pixels.
    set internal_gui_width, 446
    

# See Also

[Get](/index.php/Get "Get"), [Set](/index.php/Set "Set"), [GUI Category](/index.php/Category:GUI "Category:GUI"),[Internal_gui](/index.php/Internal_gui "Internal gui")

Retrieved from "[https://pymolwiki.org/index.php?title=Internal_gui_width&oldid=6166](https://pymolwiki.org/index.php?title=Internal_gui_width&oldid=6166)"


---

## Internal feedback

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

The Internal_feedback is the openGL-based command line area. Setting this to 0 removes the area; setting it to 1 makes one line visible (so you can see what you're typing); or you can set it to a higher value and see more lines internally (see images below). 

  * [![Internal_feedback set to 1](/images/f/f3/InternalFeedback1.png)](/index.php/File:InternalFeedback1.png "Internal_feedback set to 1")

Internal_feedback set to 1 

  * [![Internal_feedback set to 10.](/images/2/28/InternalFeedback10.png)](/index.php/File:InternalFeedback10.png "Internal_feedback set to 10.")

Internal_feedback set to 10. 




# Syntax
    
    
    # set the number of lines shown to 'int', where int is
    # an integer >= 0.
    set internal_feedback, int
    
    # For example, show 10 lines of internal feedback
    set internal_feedback, 10
    

# See Also

[GUI Category](/index.php/Category:GUI "Category:GUI"), [Internal_gui](/index.php/Internal_gui "Internal gui")

Retrieved from "[https://pymolwiki.org/index.php?title=Internal_feedback&oldid=6162](https://pymolwiki.org/index.php?title=Internal_feedback&oldid=6162)"


---

## Keyboard Shortcut Menu

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The Keyboard Shortcut Menu is a GUI panel that can create, modify, delete, and reset keyboard shortcuts. This functionality mimics the [Set_Key](/index.php/Set_Key "Set Key") function. 

_New in PyMOL 2.5_

## Usage

The Keyboard Shortcut Menu is located in the "Setting" menu heading. 

[![Keyboard-Shortcut-Menu-Location.png](/images/9/95/Keyboard-Shortcut-Menu-Location.png)](/index.php/File:Keyboard-Shortcut-Menu-Location.png)

From this table, shortcuts can be edited directly by clicking and typing in the new command. 

[![Shortcut-Menu-Edit.png](/images/3/3e/Shortcut-Menu-Edit.png)](/index.php/File:Shortcut-Menu-Edit.png)

New shortcuts can be created by clicking the "Create New" button. This will open a new dialog that will fill the first box with the key you press (e.g. CTRL-K). The desired command can then be typed in the "Command" box. Click the "Help" button for useful tips for making new shortcuts. 

[![Shortcut-Menu-Create-New.png](/images/7/72/Shortcut-Menu-Create-New.png)](/index.php/File:Shortcut-Menu-Create-New.png)

Selected commands can also be deleted or reset to their default values with the "Delete Selected" and "Reset Selected" buttons respectively. 

[![Shortcut-Menu-Delete-Reset.png](/images/0/0e/Shortcut-Menu-Delete-Reset.png)](/index.php/File:Shortcut-Menu-Delete-Reset.png)

The "Reset All" button will restore all commands to their default values and remove any new shortcuts that have been created. 

The "Save" button will save the current shortcut configuration to a file. This configuration will be automatically loaded in when you open PyMOL. 

## See Also

  * [Set_Key](/index.php/Set_Key "Set Key")



Retrieved from "[https://pymolwiki.org/index.php?title=Keyboard_Shortcut_Menu&oldid=13428](https://pymolwiki.org/index.php?title=Keyboard_Shortcut_Menu&oldid=13428)"


---

## Properties Inspector

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The Properties Inspector panel lets you inspect and modify object, object-state, atom, and atom-state properties and settings. It opens when clicking the "Properties" button in the upper right. 

_New in PyMOL 2.0_

## Screenshot

Screenshot of PyMOL with the properties panel open, after applying "A > preset > pretty" to PDB 1rx1: 

[![Properties-Inspector.png](/images/b/ba/Properties-Inspector.png)](/index.php/File:Properties-Inspector.png)

## Notes

  * The panel doesn't refresh automatically. After loading new objects, press the refresh button (next to the atom index field).



## See Also

  * [iterate](/index.php/Iterate "Iterate")
  * [get](/index.php/Get "Get"), [set](/index.php/Set "Set"), [unset](/index.php/Unset "Unset")
  * [get_title](/index.php/Get_Title "Get Title"), [set_title](/index.php/Set_title "Set title")



Retrieved from "[https://pymolwiki.org/index.php?title=Properties_Inspector&oldid=13345](https://pymolwiki.org/index.php?title=Properties_Inspector&oldid=13345)"


---

## Wrap output

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting controls the width of the output in the GUI. 

# Syntax
    
    
    # set it on or off
    set wrap_output, boolean
    
    # for example, turn on wrapping?
    set wrap_output, 1
    

# See Also

[Category:GUI](/index.php/Category:GUI "Category:GUI"), [Internal_gui](/index.php/Internal_gui "Internal gui")

Retrieved from "[https://pymolwiki.org/index.php?title=Wrap_output&oldid=6161](https://pymolwiki.org/index.php?title=Wrap_output&oldid=6161)"


---

