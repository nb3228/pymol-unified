# Category: Mouse

## Button

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**button** can be used to redefine what the mouse buttons do. 

### USAGE
    
    
    button <button>,<modifier>,<action>
    

### PYMOL API
    
    
    cmd.button( string button, string modifier, string action )
    

### NOTES

button: 

  * L
  * M
  * R
  * S



modifers: 

  * None
  * Shft
  * Ctrl
  * CtSh



actions: 

  * Rota
  * Move
  * MovZ
  * Clip
  * RotZ
  * ClpN
  * ClpF
  * lb
  * mb
  * rb
  * +lb
  * +lbX
  * -lbX
  * +mb
  * +rb
  * PkAt
  * PkBd
  * RotF
  * TorF
  * MovF
  * Orig
  * Cent



Retrieved from "[https://pymolwiki.org/index.php?title=Button&oldid=13578](https://pymolwiki.org/index.php?title=Button&oldid=13578)"


---

## Mouse modes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Add the code below to the file _mouse_modes.py_ and run it from within pymol. 

After you run it, click on the mouse mode indicator to cycle it at least one time in order to get the new bindings. 

Replace `['three_button_viewing']` by `['three_button_editing']` to edit the **3-Button Editing** mode. 

## mouse_modes.py
    
    
    from pymol.controlling import ring_dict,mode_name_dict,mode_dict
    
    # redefine the three_button_viewing mode
    mode_name_dict['three_button_viewing'] = 'My 3-But View'
    mode_dict['three_button_viewing'] =  [ ('l','none','rota'),
                          ('m','none','move'),
                          ('r','none','movz'),
                          ('l','shft','+Box'),
                          ('m','shft','-Box'),
                          ('r','shft','clip'),                 
                          ('l','ctrl','+/-'),
                          ('m','ctrl','pkat'),
                          ('r','ctrl','pk1'),                 
                          ('l','ctsh','Sele'),
                          ('m','ctsh','orig'),
                          ('r','ctsh','menu'),
                          ('w','none','slab'),
                          ('w','shft','movs'),
                          ('w','ctrl','mvsz'),
                          ('w','ctsh','movz'),
                          ('double_left','none','menu'),
                          ('double_middle','none','none'),
                          ('double_right','none', 'pk1'),
                          ('single_left','none','+/-'),
                          ('single_middle','none','cent'),
                          ('single_right','none', 'pkat'),
                          ]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mouse_modes&oldid=8361](https://pymolwiki.org/index.php?title=Mouse_modes&oldid=8361)"


---

## Mouse Settings

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Mouse_Settings&oldid=6520](https://pymolwiki.org/index.php?title=Mouse_Settings&oldid=6520)"


---

## Pickable

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

The pickable setting is used to determine whether clicking atoms with the mouse can add or remove selections. 

## Syntax

To turn off click-to-select: 
    
    
    set pickable, 0
    

To turn on click-to-select (default): 
    
    
    set pickable, 1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pickable&oldid=6522](https://pymolwiki.org/index.php?title=Pickable&oldid=6522)"


---

## Three Button Motions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This is a [Config_mouse](/index.php?title=Config_mouse&action=edit&redlink=1 "Config mouse \(page does not exist\)") option that deserved a little explanation. 

Retrieved from "[https://pymolwiki.org/index.php?title=Three_Button_Motions&oldid=6518](https://pymolwiki.org/index.php?title=Three_Button_Motions&oldid=6518)"


---

