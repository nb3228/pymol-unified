# Category: Advanced Issues

## Model Space and Camera Space

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Objects are defined in a Cartesian coordinate system, i.e. as xyz coordinates for each atom (or surface point or corner or...). This coordinate system is model space. To make objects visible, PyMOL places a "camera" such that it sees most of the molecule. The picture that is shown on the screen is the image that is "taken" by the camera. Initially, after loading an object, the camera is placed such that it looks on the object parallel to the model's z axis. The model's x and y axes are horizontal and vertical, respectively. So, initially, the model axes correspond to the physical directions of the screen: x and y are horizontal and vertical, z is perpendicular to the screen. 

In normal operation, to change the view, the camera is moved and the view of all objects is changed simultaneously; that is, the objects are not moved relative to one other. After changing the view, the x, y and z axes of model space no longer correspond to the directions of the screen. The commands [turn](/index.php/Turn "Turn"), [move](/index.php/Move "Move"), [center](/index.php/Center "Center"), [zoom](/index.php/Zoom "Zoom"), and [orient](/index.php/Orient "Orient"), as well as mouse action in viewing mode, can be used to move the camera in object space. 

For easy description of object movement relative to the camera, an additional coordinate system is defined: Camera space. The camera is situated in the origin of camera space, x and y are horizontal respectively vertical and z is the viewing direction. Camera space always corresponds to the physical directions of the screen--x and y are horizontal and vertical, and z is perpendicular to the screen. The relation between camera space and model space is described in the view matrix (see [get_view](/index.php/Get_view "Get view")). 

Changing the coordinates in model space is only necessary when objects must be moved relative to each other or when the transformed coordinates are to be written into a file. Model space coordinates can be changed with the commands [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), and [transform_selection](/index.php/Transform_selection "Transform selection"), or with the mouse in editing mode. 

  


## Space and File saving

When model coordinates are [saved](/index.php/Save "Save") (e.g. as a pdb file), PyMOL uses model space coordinates. In contrast, when a specific view is saved as [PovRay](/index.php/PovRay "PovRay") or [vrml](/index.php/Vrml "Vrml"), camera space coordinates are used (unless [geometry_export_mode](/index.php/Geometry_export_mode "Geometry export mode") is turned on, in which case no camera or lighting information is included in the output file). 

When image files in model space coordinates are needed, first choose colors and representations as needed, then align camera space to model space: 
    
    
    set_view (\
         1,    0,    0,\
         0,    1,    0,\
         0,    0,    1,\
         0,    0,    0,\
         0,    0,    0,\
         0,    300,  1 )
    

After setting the view, the molecule will probably be out of sight but will be saved to the file nonetheless. This can be used e.g. to receive a list of surface points in model space. 

## See Also

[turn](/index.php/Turn "Turn"), [move](/index.php/Move "Move"), [center](/index.php/Center "Center"), [zoom](/index.php/Zoom "Zoom"), [orient](/index.php/Orient "Orient"), [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [Transform_selection](/index.php/Transform_selection "Transform selection"), [Get_View](/index.php/Get_View "Get View")

Retrieved from "[https://pymolwiki.org/index.php?title=Model_Space_and_Camera_Space&oldid=11529](https://pymolwiki.org/index.php?title=Model_Space_and_Camera_Space&oldid=11529)"


---

## Modeling and Editing Structures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Modeling in PyMOL
    * 1.1 Saving with transformed coordinates
    * 1.2 Translate or rotate individual objects
    * 1.3 Moving one segment relative to the rest
    * 1.4 Split states to objects
    * 1.5 Altering secondary structures
    * 1.6 Altering van der Waals radii
    * 1.7 Altering atom coordinates
    * 1.8 Deleting bonds
    * 1.9 Converting D- to L- amino acids
    * 1.10 Adding disulfide bonds
    * 1.11 Adding hydrogen bonds
    * 1.12 Protonating ligands
    * 1.13 Superposition of two molecules
    * 1.14 Manual superposition of two molecules
    * 1.15 Adding and using your own fragments



# Modeling in PyMOL

## Saving with transformed coordinates

Here is a simple script that saves the molecule with coordinates from the current orientation. (invoke it with 'run save_transformed.py' and type the new save_transformed.py command thereafter). 
    
    
    # Adds the command save_transformed
    # Usage: save_transformed object, file
    def save_transformed(object,file):
        m = cmd.get_view(0)
        ttt = [m[0], m[1], m[2], 0.0,
               m[3], m[4], m[5], 0.0,
               m[6], m[7], m[8], 0.0,
               0.0,   0.0,  0.0, 1.0]
        cmd.transform_object(object,ttt,transpose=1)
        cmd.save(file,object)
    
    cmd.extend('save_transformed',save_transformed)
    

  


## Translate or rotate individual objects

There is a "[translate](/index.php/Translate "Translate")" function similar to "[rotate](/index.php/Rotate "Rotate")", feel free to use them in the following forms: 
    
    
    translate vector,object-name,state
       (vector needs to be something like [x,y,z])
    
       translate [1,0,0],pept
    
    rotate axis,angle,object-name,state
       (axis can be either the letter x,y,z or a 3D vector [x,y,z])
       rotate x,90,pept
       rotate [1,1,1],10,pept
    

## Moving one segment relative to the rest

This means moving two parts of one object into different directions. The easiest way to do this is to split the objects and then use the rotate command. 
    
    
    load 1FJ1.pdb
    
    # split PDB file
    
    create anti=(chain F) 
    create fab=(chain A,B)
    
    # delete original object
    delete 1FJ1
    
    # color objects
    color green,fab
    color pink,anti
    
    # color interface
    select inter = (byres ((fab within 5 of anti)\
       or (anti within 5 of fab)))
    
    color yellow,inter
    
    # splay apart
    orient
    origin fab
    rotate y,60,fab
    origin anti
    rotate y,-60, anti
    
    # zoom interface region
    zoom inter
    show sph,inter
    disable inter
    

## [Split states](/index.php/Split_States "Split States") to objects

There is also a new command in the 0.95 series: 
    
    
       split_states object-name
    

which will spread a PDB "biological unit" (or any multi-state object -- including SD files) over a series of independent objects. This makes it possible to interact with such objects more naturally than with "all_states = 1". 

  


## Altering secondary structures

Examples: 
    
    
     alter A/10:34/, ss='H'
     alter A/35:40/, ss='L'
     alter A/41:60/, ss='S'
    

## Altering van der Waals radii

Example: 
    
    
    alter (elem Fe),vdw=1.8
    rebuild
    

(The value for Fe is wrecked in PyMOL at the moment, so running the above line might be a good idea). 

  


## Altering atom coordinates

Example: 
    
    
    alter_state 1,(pdb1cse),x=x-10.0
    

The latter section can contain formulae involving at least the xyz coordinates, lots of constants and the (+-*/) operators. 

  


## Deleting bonds

Select the bond using Ctrl-right-click, then either 
    
    
    unbond pk1,pk2
    

or hit Ctrl-D. 

  


## Converting D- to L- amino acids

The inversion function was changed in version 0.95 to take advantage of multiple picked atoms. To invert a center, Ctrl-middle-click to pick the center atom as pk1 and two stationary atoms as pk2 and pk3. Then type Ctrl-E to invert. 

  


## Adding disulfide bonds

You can use the [bond](/index.php/Bond "Bond") command to attach them: 
    
    
    bond 24/sg,26/sg
    bond 56/sg,99/sg
    unpick
    

(unpick will hide the bond baton which gets displayed.) Additionally, the residue names can be changed for bonded cysteines: 
    
    
    alter cys/,name='CYX'
    

or for specific residues 
    
    
    alter 24+26+56+99/,name='CYX'
    

  


## Adding hydrogen bonds

    _See[Displaying biochemical properties](/index.php/Displaying_Biochemical_Properties#Hydrogen_bonds_and_Polar_Contacts "Displaying Biochemical Properties")._

## Protonating ligands

If your ligands come in with valid valencies and formal charges, PyMOL's h_add command can protonate ligands. (NOTE that there is a minor technical hiccup with SD-files which are loaded by default as immutable "discrete" objects.) Suffice it to say that in order to make changes to the chemical structure, an object must be loaded with the "discrete" flag set to zero. Unfortunately, much of the molecular editing stuff remains to be documented. Here's an example sequence, but I'm not sure it will help to much...as indicated in the manual, this is immature functionality with some major gaps. Attach in particular is very limited... 
    
    
    # show valences
    set valence=0.05
    
    # load cysteine fragment
    fragment cys
    
    # remove hydrogens
    remove (hydro)
    
    # edit gamma S
    edit cys////sg
    
    # add hydrogen
    attach H,1,1
    
    # add planer, trivalent nitrogen onto C terminus
    edit cys////C
    attach N,3,3
    
    # edit that nitrogen
    edit (elem N and neighbor cys////C)
    
    # attach a tetrahedral methyl (note random position)
    attach C,4,4
    
    # here's an example of adding a whole residue from the library
    edit cys////N
    editor.attach_amino_acid("pk1","ace")
    
    # now restore missing hydrogens (note that the names are off...)
    h_add
    

  


## Superposition of two molecules

Using pair_fit requires that you specify a set of paired atoms in each structure. Fortunately, you no longer have to specify each pair separately, so long as the ordering is the same in each selection (almost always true). 
    
    
    pair_fit ( trna10 and resid 10:15 and name P ), ( ref4 and resid 10:15 and name P )
    

Another example: 
    
    
    pair_fit prot1///11-26/CA, prot2///34-49/CA
    

would superimpose prot1 on prot2 using C-alphas from residues 11-26 in prot1 and 34-49 in prot2. 

## Manual superposition of two molecules

You can also align to structures using mouse rotation/translation. For this, you need to protect those molecules you don't want to move with (action menu -> movement -> protect) in the selection menu. 

Protect one object, deprotect the other, grab the deprotected object and move with Shift-Mouse. Don't forget to switch to Mouse Editing mode. 

  


## Adding and using your own fragments

Pymol has some build-in fragments (amino acids and simple functional groups). You can add your own fragments, eg. sugars, in this way: 

Create the molecule you want to use as a fragment. Save it as a .pkl file in <pymol_path>/data/chempy/fragments. 

How to use the fragment: 

Pick the atom (ctrl-middle) where you want to add the fragment. This will usually be a hydrogen atom (which will be removed). Then use the command: 
    
    
    editor.attach_fragment('pk1','my_fragment_name',11,0)
    

where my_fragment_name is the name of the pkl-file (w/o .pkl extension) and 11 is the number of the connecting (hydrogen) atom in the fragment. To determine this number, press '[L]abel' -> 'atom identifiers' -> 'index' and choose the hydrogen atom you want. 

If you want a menu item for your fragment, you can probably put it in <pymol_path>/modules/pmg_tk/skins/normal/__init__.py, but I haven't tried this. 

Retrieved from "[https://pymolwiki.org/index.php?title=Modeling_and_Editing_Structures&oldid=12249](https://pymolwiki.org/index.php?title=Modeling_and_Editing_Structures&oldid=12249)"


---

## Object Matrix

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Every object has an object matrix associated with it. The object matrix tells PyMOL whether the object coordinates are to be transformed before the object is displayed. The matrix is not a standard matrix in mathematical sense, it is something PyMOL-specific, also called TTT matrix: It is 4X4, with the upper left 3x3 forming a rotation matrix, the fourth column and row representing pre-rotation and post-rotation translation vectors respectively, and the 16th element always being 1.0. TTT stands for translation - transformation - translation. 

Initially, the object matrix is the identity matrix, so no transformation takes place. It is changed during structure alignment (e.g. align and super). By hand, it can be changed with the commands rotate, translate, and origin. 

It can be copied to other objects ([matrix_copy](/index.php/Matrix_copy "Matrix copy")) and reset to identity ([matrix_reset](/index.php/Matrix_reset "Matrix reset")). With [get_object_matrix](/index.php/Get_object_matrix "Get object matrix"), the current matrix can be read out. 

## See also

[rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [origin](/index.php/Origin "Origin"), [get_object_matrix](/index.php/Get_object_matrix "Get object matrix"), [matrix_copy](/index.php/Matrix_copy "Matrix copy"), [matrix_reset](/index.php/Matrix_reset "Matrix reset"), [align](/index.php/Align "Align"), [super](/index.php/Super "Super")

Retrieved from "[https://pymolwiki.org/index.php?title=Object_Matrix&oldid=8447](https://pymolwiki.org/index.php?title=Object_Matrix&oldid=8447)"


---

## Category:Image Manipulation

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Pages in category "Image Manipulation"

The following 5 pages are in this category, out of 5 total. 

### A

  * [Antialias](/index.php/Antialias "Antialias")



### D

  * [Display CCP4 Maps](/index.php/Display_CCP4_Maps "Display CCP4 Maps")



### P

  * [Publication Quality Images](/index.php/Publication_Quality_Images "Publication Quality Images")



### R

  * [Ray opaque background](/index.php/Ray_opaque_background "Ray opaque background")
  * [Ray orthoscopic](/index.php/Ray_orthoscopic "Ray orthoscopic")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Image_Manipulation&oldid=2946](https://pymolwiki.org/index.php?title=Category:Image_Manipulation&oldid=2946)"


---

## Category:Modeling and Editing Structures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Modeling in PyMOL
    * 1.1 Saving with transformed coordinates
    * 1.2 Translate or rotate individual objects
    * 1.3 Moving one segment relative to the rest
    * 1.4 Split states to objects
    * 1.5 Altering secondary structures
    * 1.6 Altering van der Waals radii
    * 1.7 Altering atom coordinates
    * 1.8 Deleting bonds
    * 1.9 Converting D- to L- amino acids
    * 1.10 Adding disulfide bonds
    * 1.11 Adding hydrogen bonds
    * 1.12 Protonating ligands
    * 1.13 Superposition of two molecules
    * 1.14 Manual superposition of two molecules
    * 1.15 Adding and using your own fragments



# Modeling in PyMOL

## Saving with transformed coordinates

Here is a simple script that saves the molecule with coordinates from the current orientation. (invoke it with 'run save_transformed.py' and type the new save_transformed.py command thereafter). 
    
    
    # Adds the command save_transformed
    # Usage: save_transformed object, file
    def save_transformed(object,file):
        m = cmd.get_view(0)
        ttt = [m[0], m[1], m[2], 0.0,
               m[3], m[4], m[5], 0.0,
               m[6], m[7], m[8], 0.0,
               0.0,   0.0,  0.0, 1.0]
        cmd.transform_object(object,ttt)
        cmd.save(file,object)
    
    cmd.extend('save_transformed',save_transformed)
    

  


## Translate or rotate individual objects

There is a "translate" function similar to "rotate", the docs for these don't exist yet, because the implementation isn't finished. However, feel free to use them in the following forms: 
    
    
    translate vector,object-name,state
       (vector needs to be something like [x,y,z])
    
       translate [1,0,0],pept
    
    rotate axis,angle,object-name,state
       (axis can be either the letter x,y,z or a 3D vector [x,y,z])
       rotate x,90,pept
       rotate [1,1,1],10,pept
    

## Moving one segment relative to the rest

This means moving two parts of one object into different directions. The easiest way to do this is to split the objects and then use the rotate command. 
    
    
    load 1FJ1.pdb
    
    # split PDB file
    
    create anti=(chain F) 
    create fab=(chain A,B)
    
    # delete original object
    delete 1FJ1
    
    # color objects
    color green,fab
    color pink,anti
    
    # color interface
    select inter = (byres ((fab within 5 of anti)\
       or (anti within 5 of fab)))
    
    color yellow,inter
    
    # splay apart
    orient
    origin fab
    rotate y,60,fab
    origin anti
    rotate y,-60, anti
    
    # zoom interface region
    zoom inter
    show sph,inter
    disable inter
    

## Split states to objects

There is also a new command in the 0.95 series: 
    
    
       split_states object-name
    

which will spread a PDB "biological unit" (or any multi-state object -- including SD files) over a series of independent objects. This makes it possible to interact with such objects more naturally than with "all_states = 1". 

  


## Altering secondary structures

Examples: 
    
    
     alter A/10:34/, ss='H'
     alter A/35:40/, ss='L'
     alter A/41:60/, ss='S'
    

## Altering van der Waals radii

Example: 
    
    
    alter (elem Fe),vdw=1.8
    rebuild
    

(The value for Fe is wrecked in PyMOL at the moment, so running the above line might be a good idea). 

  


## Altering atom coordinates

Example: 
    
    
    alter_state 1,(pdb1cse),x=x-10.0
    

The latter section can contain formulae involving at least the xyz coordinates, lots of constants and the (+-*/) operators. 

  


## Deleting bonds

Select the bond using Ctrl-right-click, then either 
    
    
    unbond pk1,pk2
    

or hit Ctrl-D. 

  


## Converting D- to L- amino acids

The inversion function was changed in version 0.95 to take advantage of multiple picked atoms. To invert a center, Ctrl-middle-click to pick the center atom as pk1 and two stationary atoms as pk2 and pk3. Then type Ctrl-E to invert. 

  


## Adding disulfide bonds

You can use the [Cmd bond](/index.php/Cmd_bond "Cmd bond") command to attach them: 
    
    
    bond 24/sg,26/sg
    bond 56/sg,99/sg
    unpick
    

(unpick will hide the bond baton which gets displayed.) Additionally, the residue names can be changed for bonded cysteines: 
    
    
    alter cys/,name='CYX'
    

or for specific residues 
    
    
    alter 24+26+56+99/,name='CYX'
    

  


## Adding hydrogen bonds

See 'displaying biochemical properties'. 

  


## Protonating ligands

If your ligands come in with valid valencies and formal charges, PyMOL's h_add command can protonate ligands. (NOTE that there is a minor technical hiccup with SD-files which are loaded by default as immutable "discrete" objects.) Suffice it to say that in order to make changes to the chemical structure, an object must be loaded with the "discrete" flag set to zero. Unfortunately, much of the molecular editing stuff remains to be documented. Here's an example sequence, but I'm not sure it will help to much...as indicated in the manual, this is immature functionality with some major gaps. Attach in particular is very limited... 
    
    
    # show valences
    set valence=0.05
    
    # load cysteine fragment
    fragment cys
    
    # remove hydrogens
    remove (hydro)
    
    # edit gamma S
    edit cys////sg
    
    # add hydrogen
    attach H,1,1
    
    # add planer, trivalent nitrogen onto C terminus
    edit cys////C
    attach N,3,3
    
    # edit that nitrogen
    edit (elem N and neighbor cys////C)
    
    # attach a tetrahedral methyl (note random position)
    attach C,4,4
    
    # here's an example of adding a whole residue from the library
    edit cys////N
    editor.attach_amino_acid("pk1","ace")
    
    # now restore missing hydrogens (note that the names are off...)
    h_add
    

  


## Superposition of two molecules

Using pair_fit requires that you specify a set of paired atoms in each structure. Fortunately, you no longer have to specify each pair separately, so long as the ordering is the same in each selection (almost always true). 
    
    
    pair_fit ( trna10 and resid 10:15 and name P ), ( ref4 and resid 10:15 and name P )
    

Another example: 
    
    
    pair_fit prot1///11-26/CA, prot2///34-49/CA
    

would superimpose prot1 on prot2 using C-alphas from residues 11-26 in prot1 and 34-49 in prot2. 

## Manual superposition of two molecules

You can also align to structures using mouse rotation/translation. For this, you need to protect those molecules you don't want to move with (action menu -> movement -> protect) in the selection menu. 

Protect one object, deprotect the other, grab the deprotected object and move with Shift-Mouse. Don't forget to switch to Mouse Editing mode. 

  


## Adding and using your own fragments

Pymol has some build-in fragments (amino acids and simple functional groups). You can add your own fragments, eg. sugars, in this way: 

Create the molecule you want to use as a fragment. Save it as a .pkl file in <pymol_path>/data/chempy/fragments. 

How to use the fragment: 

Pick the atom (ctrl-middle) where you want to add the fragment. This will usually be a hydrogen atom (which will be removed). Then use the command: 
    
    
    editor.attach_fragment('pk1','my_fragment_name',11,0)
    

where my_fragment_name is the name of the pkl-file (w/o .pkl extension) and 11 is the number of the connecting (hydrogen) atom in the fragment. To determine this number, press '[L]abel' -> 'atom identifiers' -> 'index' and choose the hydrogen atom you want. 

If you want a hotkey for your fragment, you can probably put it in <pymol_path>/modules/pmg_tk/skins/normal/__init__.py, but I haven't tried this. 

## Pages in category "Modeling and Editing Structures"

The following 6 pages are in this category, out of 6 total. 

### C

  * [CavitOmiX](/index.php/CavitOmiX "CavitOmiX")



### E

  * [Editing atoms](/index.php/Editing_atoms "Editing atoms")



### H

  * [Homology Modeling](/index.php/Homology_Modeling "Homology Modeling")



### M

  * [Molecular Sculpting](/index.php/Molecular_Sculpting "Molecular Sculpting")
  * [Morph](/index.php/Morph "Morph")



### R

  * [Rigimol.morph](/index.php/Rigimol.morph "Rigimol.morph")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Modeling_and_Editing_Structures&oldid=3224](https://pymolwiki.org/index.php?title=Category:Modeling_and_Editing_Structures&oldid=3224)"


---

