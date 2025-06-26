# Category: View Module

## Clip

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**clip** alters the near and far clipping planes 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



### USAGE
    
    
    clip {near|far|move|slab|atoms}, distance [,selection [,state ]]
    

### EXAMPLES
    
    
    clip near, -5           # moves near plane away from you by 5 A
    clip far, 10            # moves far plane towards you by 10 A
    clip move, -5           # moves the slab away from you by 5 A
    clip slab, 20           # sets slab thickness to 20 A
    clip slab, 10, resi 11  # clip 10 A slab about residue 11
    clip atoms, 5, pept     # clip atoms in "pept" with a 5 A buffer
                            # about their current camera positions
    

## PYMOL API
    
    
    cmd.clip( string mode, float distance, string selection = None)
    

## SEE ALSO

  * [Zoom](/index.php/Zoom "Zoom")
  * [Reset](/index.php/Reset "Reset")



Retrieved from "[https://pymolwiki.org/index.php?title=Clip&oldid=7443](https://pymolwiki.org/index.php?title=Clip&oldid=7443)"


---

## Disable

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **disable** command toggles off the display of all currently visible representations of an object. It is the equivalent of deselecting the object in the list in the top panel of the [Internal GUI](/index.php/Internal_gui "Internal gui"). 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    disable [name]
    

name
    the name of an object or a named selection (default="all")

### PYMOL API
    
    
    cmd.disable( string name='all' )
    

### EXAMPLES
    
    
    disable my_object
    disable (my_object1 or my_object2)
    disable my_object*
    disable
    

### SEE ALSO

[Show](/index.php/Show "Show"), [Hide](/index.php/Hide "Hide"), [Enable](/index.php/Enable "Enable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates")

Retrieved from "[https://pymolwiki.org/index.php?title=Disable&oldid=11561](https://pymolwiki.org/index.php?title=Disable&oldid=11561)"


---

## Enable

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **enable** command toggles on the display of all currently visible representations of an object. It is the equivalent of selecting the object in the list at the top of the [Internal GUI](/index.php/Internal_gui "Internal gui"). 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    enable [name [, parents]]
    

name
    the name of an object or a named selection (default="all")

parents
    when set to 1, enables the named object and all parent objects of the named object (recursively), so it won't be hidden by having a parent object (or [group](/index.php/Group "Group")) that is [disabled](/index.php/Disable "Disable") (default=0)

### PYMOL API
    
    
    cmd.enable( string name='all', int parents=0 )
    

### EXAMPLES
    
    
    enable my_object
    enable (my_object1 or my_object2)
    enable my_object*
    enable
    enable my_object, parents=1
    

### SEE ALSO

[Show](/index.php/Show "Show"), [Hide](/index.php/Hide "Hide"), [Disable](/index.php/Disable "Disable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates")

Retrieved from "[https://pymolwiki.org/index.php?title=Enable&oldid=11562](https://pymolwiki.org/index.php?title=Enable&oldid=11562)"


---

## Get View

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_view** returns and optionally prints out the current view information in a format which can be embedded into a command script and can be used in subsequent calls to **set_view**. 

If a log file is currently open, get_view will not write the view matrix to the screen unless the "output" parameter is 2. 

This command is very useful for saving the orientation of a scene for later. Authors of molecular movies may find this command very powerful. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 API USAGE
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    get_view
    

### PYMOL API
    
    
    cmd.get_view(output=1, quiet=1)
       
    my_view = cmd.get_view()
    

output control: 

  * 0 = output matrix to screen
  * 1 = don't output matrix to screen
  * 2 = force output to screen even if log file is open
  * 3 = return formatted string instead of a list



### API USAGE
    
    
    cmd.get_view(0) # zero option suppresses output (LEGACY approach)
    cmd.get_view(quiet=1) # suppresses output using PyMOL's normal "quiet" parameter.
    

### NOTES

Contents of the view matrix 

  * 0 - 8 = column-major 3x3 matrix which rotates model axes to camera axes
  * 9 - 11 = origin of rotation relative to the camera in camera space
  * 12 - 14 = origin of rotation in model space
  * 15 = front plane distance from the camera
  * 16 = rear plane distance from the camera
  * 17 = orthoscopic flag (not implemented in older versions)



### SEE ALSO

[Set View](/index.php/Set_View "Set View"), [View](/index.php/View "View"), [Model_Space_and_Camera_Space](/index.php/Model_Space_and_Camera_Space "Model Space and Camera Space")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_View&oldid=13321](https://pymolwiki.org/index.php?title=Get_View&oldid=13321)"


---

## Hide

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**hide** conceals atom and bond representations for a certain selection or other graphical objects like distances. 

  * [![Some normal scene, notice the waters shown as spheres](/images/6/6b/Show1.png)](/index.php/File:Show1.png "Some normal scene, notice the waters shown as spheres")

Some normal scene, notice the waters shown as spheres 

  * [![hide spheres issues and the spheres are now hidden.](/images/e/ea/Show2.png)](/index.php/File:Show2.png "hide spheres issues and the spheres are now hidden.")

_hide spheres_ issues and the spheres are now hidden. 




The available representations are: 

  * lines
  * spheres
  * mesh
  * ribbon
  * cartoon
  * sticks
  * dots
  * surface
  * labels
  * nonbonded
  * nb_spheres



## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    hide representation [,object]
    hide representation [,(selection)]
    hide (selection)
    

### PYMOL API
    
    
    cmd.hide( string representation="", string selection="")
    

### EXAMPLES
    
    
    # hides all lines
    hide lines,all
    
    # hides all ribbons
    hide ribbon
    
    # hides all distances
    hide dashes
    
    # hides sticks in protA and all residues that aren't in the range of 40-65
    hide sticks, protA and not i. 40-65
    
    # hide hydrogen atoms
    hide (hydro)  # or hide (h.)
    

### SEE ALSO

[Show](/index.php/Show "Show"), [Enable](/index.php/Enable "Enable"), [Disable](/index.php/Disable "Disable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates")

Retrieved from "[https://pymolwiki.org/index.php?title=Hide&oldid=11564](https://pymolwiki.org/index.php?title=Hide&oldid=11564)"


---

## Move

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**move** translates the camera along one of the three primary view axes. 

## Contents

  * 1 USAGE
    * 1.1 EXAMPLES
  * 2 PYMOL API
  * 3 SEE ALSO



### USAGE
    
    
    move axis,distance
    

#### EXAMPLES
    
    
    move x,3
    move y,-1
    

### PYMOL API
    
    
    cmd.move( string axis, float distance )
    

### SEE ALSO

[turn](/index.php/Turn "Turn"), [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [zoom](/index.php/Zoom "Zoom"), [center](/index.php/Center "Center"), [clip](/index.php/Clip "Clip"), [Translate](/index.php/Translate "Translate")

Retrieved from "[https://pymolwiki.org/index.php?title=Move&oldid=11634](https://pymolwiki.org/index.php?title=Move&oldid=11634)"


---

## Orient

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**orient** aligns the principal components of the atoms in the selection with the XYZ axes. The function is similar to the orient command in X-PLOR. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 EXAMPLES
  * 5 SEE ALSO



### USAGE
    
    
    orient object-or-selection [, state]
    orient (selection)
    

### PYMOL API
    
    
    cmd.orient( string object-or-selection [, state = 0] )
    

### NOTES
    
    
      state = 0 (default) use all coordinate states
      state = -1 use only coordinates for the current state
      state > 0  use coordinates for a specific state
    

### EXAMPLES

For models with NCS symmetry, orient will align the model with the symmetry axis centered along the viewport's z axis. For example, 
    
    
      fetch 1hiw, async=0
      as cartoon
      remove (!chain A,B,C)
      orient
      util.cbc
    

will produce the first image below. However, if there is a larger symmetry, e.g. two trimers, this will not work. In the above example, leaving out `remove (!chain A,B,C)` from the script results in the second image below. 

  * [![One trimer from 1hiw after "orient" command.](/images/8/8d/1hiw_orient.png)](/index.php/File:1hiw_orient.png "One trimer from 1hiw after "orient" command.")

One trimer from 1hiw after "orient" command. 

  * [![Both trimers from 1hiw after "orient" command.](/images/8/86/1hiw_orient2.png)](/index.php/File:1hiw_orient2.png "Both trimers from 1hiw after "orient" command.")

Both trimers from 1hiw after "orient" command. 




### SEE ALSO

[Zoom](/index.php/Zoom "Zoom"), [Origin](/index.php/Origin "Origin"), [Reset](/index.php/Reset "Reset")

Retrieved from "[https://pymolwiki.org/index.php?title=Orient&oldid=8464](https://pymolwiki.org/index.php?title=Orient&oldid=8464)"


---

## Origin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**origin** sets the center of rotation about a selection. If an object name is specified, it can be used to set the center of rotation for the object's [TTT matrix](/index.php/Object_Matrix "Object Matrix"). 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO
  * 5 Example



### USAGE
    
    
    origin selection [, object [,position, [, state]]]
    origin (selection)
    origin position=[1.0,2.0,3.0]
    

### PYMOL API
    
    
    cmd.origin( string object-or-selection )
    

### NOTES

  * state = 0 (default) use all coordinate states
  * state = -1 use only coordinates for the current state
  * state > 0 use coordinates for a specific state



### SEE ALSO

[Zoom](/index.php/Zoom "Zoom"), [Orient](/index.php/Orient "Orient"), [Reset](/index.php/Reset "Reset")

  


### Example

This example puts the camera 'inside' a protein near the ligand and turns the camera 360 degrees about that location. 
    
    
    load $TUT/1hpv.pdb
    
    set field_of_view, 45
    
    zoom organic
    orient organic
    show stick, polymer
    show surface
    set surface_color, white
    set transparency, 0.5
    clip slab, 30
    set line_width, 5
    
    cmd.move("z",-cmd.get_view()[11])
    
    python
    def spin():
       from pymol import cmd
       for x in range(360):
          cmd.move("z",cmd.get_view()[11])
          cmd.turn("y",1)
          cmd.move("z",-cmd.get_view()[11])
          cmd.refresh()
    
    threading.Thread(target=spin).start()
    python end
    

Retrieved from "[https://pymolwiki.org/index.php?title=Origin&oldid=9541](https://pymolwiki.org/index.php?title=Origin&oldid=9541)"


---

## Rebuild

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rebuild** forces PyMOL to recreate geometric objects in case any of them have gone out of sync. 

### USAGE
    
    
    rebuild [selection [, representation ]]
    

### PYMOL API
    
    
    cmd.rebuild(string selection = 'all', string representation = 'everything')
    

### SEE ALSO

[refresh](/index.php/Refresh "Refresh")

Retrieved from "[https://pymolwiki.org/index.php?title=Rebuild&oldid=7600](https://pymolwiki.org/index.php?title=Rebuild&oldid=7600)"


---

## Refresh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**refresh** causes the scene to be refresh as soon as it is safe to do so. 

### USAGE

refresh 

### PYMOL API
    
    
    cmd.refresh()
    

### SEE ALSO

[rebuild ](/index.php/Rebuild "Rebuild")

Retrieved from "[https://pymolwiki.org/index.php?title=Refresh&oldid=7590](https://pymolwiki.org/index.php?title=Refresh&oldid=7590)"


---

## Reset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**reset** restores the rotation matrix to identity, sets the origin to the center of mass (approx.) and zooms the window and clipping planes to cover all objects. 

### USAGE
    
    
    reset
    

### PYMOL API
    
    
    cmd.reset()
    

## SEE ALSO

[Clip](/index.php/Clip "Clip") [Zoom](/index.php/Zoom "Zoom")

Retrieved from "[https://pymolwiki.org/index.php?title=Reset&oldid=8323](https://pymolwiki.org/index.php?title=Reset&oldid=8323)"


---

## Rock

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rock** toggles Y axis rocking. 

### USAGE
    
    
    rock
    

### PYMOL API
    
    
    cmd.rock()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Rock&oldid=7587](https://pymolwiki.org/index.php?title=Rock&oldid=7587)"


---

## Set View

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_view** sets viewing information for the current scene, including the rotation matrix, position, origin of rotation, clipping planes, and the orthoscopic flag. 

This command is extremely useful for making movies. One may set up the scene to be rendered, then save the exact orientation, with respect to the camera, of the scene using, the [Get_View](/index.php/Get_View "Get View") command. The output from the [Get_View](/index.php/Get_View "Get View") command may then be used by the **set_view** command to restore the orientation of the scene. 

## Contents

  * 1 USAGE
  * 2 EXAMPLE
  * 3 PYMOL API
    * 3.1 Example
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    set_view (...)  where ... is 18 floating point numbers
    

### EXAMPLE

This works in a `.pml` script: 
    
    
       set_view (\
           0.999876618,   -0.000452542,   -0.015699286,\
           0.000446742,    0.999999821,   -0.000372844,\
           0.015699454,    0.000365782,    0.999876678,\
           0.000000000,    0.000000000, -150.258514404,\
           11.842411041,   20.648729324,    8.775371552,\
           118.464958191,  182.052062988,    0.000000000 )
    

### PYMOL API
    
    
    cmd.set_view(string-or-sequence view)
    

#### Example

This works in a `.py` script: 
    
    
    cmd.set_view((
        0.999876618,   -0.000452542,   -0.015699286,
        0.000446742,    0.999999821,   -0.000372844,
        0.015699454,    0.000365782,    0.999876678,
        0.000000000,    0.000000000, -150.258514404,
        11.842411041,   20.648729324,    8.775371552,
        118.464958191,  182.052062988,    0.000000000))
    

The result of [get_view](/index.php/Get_view "Get view") is valid input for **set_view** : 
    
    
    myview = cmd.get_view()
    
    cmd.set_view(myview)
    

### NOTES

Contents of the view matrix 

  * 0 - 8: column-major 3x3 matrix which rotates model axes to camera axes
  * 9 - 11: origin of rotation relative to the camera in camera space
  * 12 - 14: origin of rotation in model space
  * 15: front plane distance from the camera
  * 16: rear plane distance from the camera
  * 17: orthoscopic flag (not implemented in older versions)



### SEE ALSO

[Get View](/index.php/Get_View "Get View"), [View](/index.php/View "View")

Retrieved from "[https://pymolwiki.org/index.php?title=Set_View&oldid=13269](https://pymolwiki.org/index.php?title=Set_View&oldid=13269)"


---

## Show

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Show** displays atom and bond representations for certain selections. 

## Contents

  * 1 Details
    * 1.1 USAGE
    * 1.2 PYMOL API
    * 1.3 EXAMPLES
      * 1.3.1 Example
      * 1.3.2 Example
      * 1.3.3 Example
      * 1.3.4 Example
    * 1.4 NOTES
    * 1.5 SEE ALSO



# Details

  * [![Example of a shown surface.](/images/e/e7/Ray2.png)](/index.php/File:Ray2.png "Example of a shown surface.")

Example of a shown surface. 

  * [![Ball and sticks shown.](/images/7/76/Ray_trace_gain2.png)](/index.php/File:Ray_trace_gain2.png "Ball and sticks shown.")

Ball and sticks shown. 

  * [![A cartoon inside a mesh shown.](/images/6/6b/Mesh_w05.png)](/index.php/File:Mesh_w05.png "A cartoon inside a mesh shown.")

A cartoon inside a mesh shown. 




The **Show** command, is one of the most often used commands in PyMOL. For example, you can _show_ certain atoms as _[Lines](/index.php/Lines "Lines")_ or _[Sticks](/index.php/Sticks "Sticks")_ or _[Cartoons](/index.php/Cartoon "Cartoon")_ or any of the following representations: 

  * [lines](/index.php/Lines "Lines")
  * [spheres](/index.php/Spheres "Spheres")
  * [mesh](/index.php/Mesh "Mesh")
  * [ribbon](/index.php/Ribbon "Ribbon")
  * [cartoon](/index.php/Cartoon "Cartoon")
  * [sticks](/index.php/Sticks "Sticks")
  * [dots](/index.php/Dots "Dots")
  * [surface](/index.php/Surface "Surface")
  * [label](/index.php/Label "Label")
  * [extent](/index.php?title=Extent&action=edit&redlink=1 "Extent \(page does not exist\)")
  * [nonbonded](/index.php?title=Nonbonded&action=edit&redlink=1 "Nonbonded \(page does not exist\)")
  * [nb_spheres](/index.php/Nb_spheres "Nb spheres")
  * [slice](/index.php/Slice "Slice")
  * [cell](/index.php/Cell "Cell")



## USAGE
    
    
    show
    show reprentation [,object]
    show reprentation [,(selection)]
    show (selection)
    

## PYMOL API
    
    
    cmd.show( string representation="", string selection="" )
    

## EXAMPLES

#### Example
    
    
    # show the backbone using lines.
    show lines,(name ca or name c or name n)
    

#### Example
    
    
    # show the ribbon representation for all objects
    show ribbon
    

#### Example
    
    
    # show all hetero atoms as spheres
    show spheres, het
    

#### Example
    
    
    # show only polar hydrogens
    hide everything, ele h
    show lines, ele h and neighbor (ele n+o)
    # hide nonpolar hydrogens
    hide (h. and (e. c extend 1))
    

Note: 

The above code hides all representations of nonpolar hydrogens, including surface representations, resulting in broken surface representations. It might be better to remove the nonpolar hydrogens instead: 
    
    
    # show only polar hydrogens
    hide everything, ele h
    show lines, ele h and neighbor (ele n+o)
    # remove nonpolar hydrogens
    remove (h. and (e. c extend 1))
    

## NOTES

**selection** can be an object name. 

**show** alone will turn on lines and nonbonded for all bonds. 

**show cell** will show the triclinic unit cell, provided it's described in the PDB. For cell packing, you need a script. See [The Script Library](/index.php/Category:Script_Library "Category:Script Library") or [Robert Campbell's Site](http://adelie.biochem.queensu.ca/~rlc/work/pymol/)

## SEE ALSO

[Hide](/index.php/Hide "Hide"), [Enable](/index.php/Enable "Enable"), [Disable](/index.php/Disable "Disable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates"), [Lines](/index.php/Lines "Lines"), [Cartoon](/index.php/Cartoon "Cartoon"), [Spheres](/index.php/Spheres "Spheres"), [Mesh](/index.php/Mesh "Mesh"), [Ribbon](/index.php/Ribbon "Ribbon"), [Sticks](/index.php/Sticks "Sticks"), [Dots](/index.php/Dots "Dots"), [Surface](/index.php/Surface "Surface"), [Label](/index.php/Label "Label"), [Extent](/index.php?title=Extent&action=edit&redlink=1 "Extent \(page does not exist\)"), [Slice](/index.php/Slice "Slice"), [Cell](/index.php/Cell "Cell"), [Select](/index.php/Select "Select"), [Show_as](/index.php/Show_as "Show as")

Retrieved from "[https://pymolwiki.org/index.php?title=Show&oldid=11744](https://pymolwiki.org/index.php?title=Show&oldid=11744)"


---

## Show as

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

show_as turns on and off atom and bond [representations](/index.php/Category:Representations "Category:Representations"). 

## Contents

  * 1 Details
    * 1.1 USAGE
    * 1.2 ARGUMENTS
    * 1.3 EXAMPLES
    * 1.4 PYMOL API
    * 1.5 NOTES
    * 1.6 SEE ALSO
    * 1.7 References



# Details

The available representations are the usual: 

(History: Python 2.6 came out with the keyword [as](/index.php/As "As") which is also a PyMOL keyword. So, we had to change the PyMOL-named keyword to show_as.) 

  * [lines](/index.php/Lines "Lines")
  * [spheres](/index.php/Spheres "Spheres")
  * [mesh](/index.php/Mesh "Mesh")
  * [ribbon](/index.php/Ribbon "Ribbon")
  * [cartoon](/index.php/Cartoon "Cartoon")
  * [sticks](/index.php/Sticks "Sticks")
  * [dots](/index.php/Dots "Dots")
  * [surface](/index.php/Surface "Surface")
  * [labels](/index.php?title=Labels&action=edit&redlink=1 "Labels \(page does not exist\)")
  * [extent](/index.php?title=Extent&action=edit&redlink=1 "Extent \(page does not exist\)")
  * [nonbonded](/index.php?title=Nonbonded&action=edit&redlink=1 "Nonbonded \(page does not exist\)")
  * [nb_spheres](/index.php/Nb_spheres "Nb spheres")
  * [slice](/index.php/Slice "Slice")



## USAGE
    
    
    show_as representation [, selection ]
    

## ARGUMENTS

  * **representation** = lines, spheres, mesh, ribbon, cartoon, sticks, dots, surface, labels, extent, nonbonded, nb_spheres, slice, extent, slice, dashes, angles, dihedrals, cgo, cell, callback, everything
  * **selection** = string {default: all}



## EXAMPLES
    
    
     
    # show the backbone as lines
    show_as lines, name ca or name c or name n
    
    # show everything as a ribbon
    show_as ribbon
    

## PYMOL API
    
    
     
    cmd.show_as(string representation, string selection)
    

## NOTES

  * **selection** can be an object name
  * **as** alone will turn on lines and nonbonded and hide everything else.



## SEE ALSO

[show](/index.php/Show "Show"), [hide](/index.php/Hide "Hide"), [enable](/index.php/Enable "Enable"), [disable](/index.php/Disable "Disable"), [All PyMOL Representations](/index.php/Category:Representations "Category:Representations")

## References

  * PyMOL Source code



Retrieved from "[https://pymolwiki.org/index.php?title=Show_as&oldid=7571](https://pymolwiki.org/index.php?title=Show_as&oldid=7571)"


---

## Turn

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**turn** rotates the camera about one of the three primary axes, centered at the origin. Also consider [rotate](/index.php/Rotate "Rotate"). 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



### USAGE
    
    
    turn axis, angle
    

### EXAMPLES
    
    
    turn x,90
    turn y,45
    

### PYMOL API
    
    
    cmd.turn( string axis, float angle )
    

### SEE ALSO

[move](/index.php/Move "Move"), [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [zoom](/index.php/Zoom "Zoom"), [center](/index.php/Center "Center"), [clip](/index.php/Clip "Clip")

Retrieved from "[https://pymolwiki.org/index.php?title=Turn&oldid=7556](https://pymolwiki.org/index.php?title=Turn&oldid=7556)"


---

## View

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**view** makes it possible to save and restore viewpoints on a given scene within a single session. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 FUNCTION KEY PRESETS
  * 4 EXAMPLES
  * 5 SEE ALSO



### USAGE
    
    
    view key[,action]
    view *
    

key can be any string action should be 'store' or 'recall' (default: 'recall') 

### PYMOL API
    
    
    cmd.view(string key,string action)
    

### FUNCTION KEY PRESETS

Views F1 through F12 are automatically bound to function keys provided that "set_key" has not been used to redefine the behaviour of the respective key, and that a "scene" has not been defined for that key. 

### EXAMPLES
    
    
    view 0,store
    view 0
    

### SEE ALSO

[Scene](/index.php/Scene "Scene"), [Set View](/index.php/Set_View "Set View"), [Get View](/index.php/Get_View "Get View")

Retrieved from "[https://pymolwiki.org/index.php?title=View&oldid=11633](https://pymolwiki.org/index.php?title=View&oldid=11633)"


---

## Zoom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**zoom** scales and translates the window and the origin to cover the atom selection. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    zoom [ selection [,buffer [, state [, complete ]]]]
    

### EXAMPLES
    
    
    # auto zoom depending on what's loaded in PyMOL
    zoom
    
    #
    zoom complete=1
    
    # zoom on just chain A
    zoom (chain A)
    
    # zoom on residue 142
    zoom 142/
    
    # zoom consistenly 20 Ang from each object at the center
    center prot1
    zoom center, 20
    
    # prot1 and prot2 will have the same exact zoom factor
    center prot2
    zoom center, 20
    

### PYMOL API
    
    
     
    cmd.zoom( string selection, float buffer=0.0,
              int state=0, int complete=0 )
    

### NOTES
    
    
    state = 0 (default) use all coordinate states
    state = -1 use only coordinates for the current state
    state > 0  use coordinates for a specific state
    
    
    
    complete = 0 or 1:
    

Normally the zoom command tries to guess an optimal zoom level for visualization, balancing closeness against occasional clipping of atoms out of the field of view. You can change this behavior by setting the complete option to 1, which will guarantee that the atom positions for the entire selection will fit in the field of an orthoscopic view. To absolutely prevent clipping, you may also need to add a buffer (typically 2 A) to account for the perspective transformation and for graphical representations which extend beyond the atom coordinates. 

### SEE ALSO

  * [Origin](/index.php/Origin "Origin")
  * [Orient](/index.php/Orient "Orient")
  * [Center](/index.php/Center "Center")



Retrieved from "[https://pymolwiki.org/index.php?title=Zoom&oldid=7545](https://pymolwiki.org/index.php?title=Zoom&oldid=7545)"


---

