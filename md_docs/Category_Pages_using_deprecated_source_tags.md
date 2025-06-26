# Category: Pages using deprecated source tags

## AAindex

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/aaindex.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/aaindex.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.aaindex](http://pymol.org/psicoredirect.php?psico.aaindex)  
  
## Contents

  * 1 Introduction
  * 2 Python Example
  * 3 PyMOL Example
  * 4 See Also



## Introduction

[![](/images/4/44/AAindexExample.png)](/index.php/File:AAindexExample.png)

[](/index.php/File:AAindexExample.png "Enlarge")

Hydrophobicity coloring with KYTJ820101

AAindex is a database of numerical indices representing various physicochemical and biochemical properties of amino acids and pairs of amino acids. See <http://www.genome.jp/aaindex/>

This script is a python parser for the AAindex flat files which will be downloaded from <ftp://ftp.genome.jp/pub/db/community/aaindex/>

The script provides two PyMOL commands (but can also be used without PyMOL). 

  * aaindex2b: Loads numerical indices from aaindex1 as b-factors into your structure
  * pmf: Potential of Mean Force (aaindex3)



## Python Example

Consider the script is called aaindex.py, it is placed somewhere in your PYTHONPATH and the aaindex flatfiles are found in the current directory. 
    
    
    import aaindex
    aaindex.init(path='.')
    
    aaindex.grep('volume')
    x = aaindex.get('KRIW790103')
    print(x)
    print(x.get('A'))
    
    aaindex.init(index='2')
    aaindex.grep('blosum')
    x = aaindex.get('HENS920102')
    print(x.get('A', 'K'))
    

## PyMOL Example

Solvent-accessible surface coloring by [Hydropathy index (Kyte-Doolittle, 1982)](https://www.genome.jp/dbget-bin/www_bget?aaindex:KYTJ820101)
    
    
    run aaindex.py
    aaindex2b KYTJ820101
    spectrum b, white yellow forest
    set surface_solvent
    show surface
    

## See Also

  * Protscale from [rTools](http://www.rubor.de/pymol_extensions_de.html) does a similar job in coloring by amino acid properties



Retrieved from "[https://pymolwiki.org/index.php?title=AAindex&oldid=13884](https://pymolwiki.org/index.php?title=AAindex&oldid=13884)"


---

## Annotate v

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/annotate_v.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/annotate_v.py)  
Author(s)  | [Xin Yu](/index.php/User:XinYu "User:XinYu")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Annotation Schemes
  * 3 Dependencies and Limitations
  * 4 How to use
  * 5 Contact



## Description

_Version: 1.0_

A simple-to-use script for annotation of VH or VL sequence of an antibody. It creates a Pymol object for each FR or CDR region. It utilizes the REST API interface of _Abnum_ (<http://www.bioinf.org.uk/abs/abnum/>) from Dr Andrew Martin's group at UCL 

## Annotation Schemes

Currently supports **Kabat, Chothia, Contact, IMGT** 1

1Definitions for Kabat, Chothia, Contact, and IIMGT are same as listed in the table (<http://www.bioinf.org.uk/abs/info.html#kabatnum>), except that for IMGT, H-CDR2 is defined as **H51-H57** in this script, as opposed to of **H51-H56** in the table. This slight revision generates result that matches that from _IMGT website_ (<http://www.imgt.org/>) 

## Dependencies and Limitations

1\. Import **request** module 

2\. Relies on internet connection to _Abnum_

3\. Incomplete VH or VL sequence might not be annotated 

## How to use

1\. Create a selection object for the V region, (such as **VH** and **VL** shown here). Only 1 V-region (VH OR VL) can be annotated each time. Alternatively, simply select the residues into the default **< sele>** group. 

  


[![](/images/8/82/Annotate_v_pic1.png)](/index.php/File:Annotate_v_pic1.png)

[](/index.php/File:Annotate_v_pic1.png "Enlarge")

Create a selection group for input seq, or use the default sele

2\. Copy the entire script and create a `.py` file. Run script. 

3\. The general syntax is: 
    
    
    annotate_v("selection_group_name", "scheme")
    

**selection_group_name:** name of the selection group where you have the input sequence. Must be single-letter amino acid coded. Either VH or VL but not both 

**scheme:** currently supports Kabat, Chothia, Contact, IMGT. Must be lowercase 

For example: 
    
    
    annotate_v("VH", "kabat") #input sequence from VH selection group, annotate using "kabat"
    annotate_v("sele", "imgt") #input sequence from the default sele group, annotate using "imgt". You can also try "contact", "chothia"
    

  
4\. In the output window, the script will print the FR and CDR regions (if found). It also automatically creates a selection group for each of the FR and CDR region. 

[![](/images/6/63/Annotate_v_pic2.png)](/index.php/File:Annotate_v_pic2.png)

[](/index.php/File:Annotate_v_pic2.png "Enlarge")

example output from command lines

[![](/images/7/72/Annotate_v_pic3.png)](/index.php/File:Annotate_v_pic3.png)

[](/index.php/File:Annotate_v_pic3.png "Enlarge")

new FR and CDR selection groups are created

  


## Contact

For bugs & questions, emai @ xinyu18018@gmail.com. Enjoy! 

Retrieved from "[https://pymolwiki.org/index.php?title=Annotate_v&oldid=13886](https://pymolwiki.org/index.php?title=Annotate_v&oldid=13886)"


---

## Average b

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 average_b.py
    * 1.1 usage
    * 1.2 author
    * 1.3 code



# average_b.py

  * calculate the average B-factor of a selection.



## usage

  * copy code to text file and save as average_b.py
  * type "run average_b.py" in PyMOL
  * then type "average_b (selection)"



## author

Gregor Hagelueken 

## code
    
    
    from pymol import cmd,stored
    def average_b(selection):
    	stored.tempfactor = 0
    	stored.atomnumber = 0
    	cmd.iterate(selection, "stored.tempfactor = stored.tempfactor + b")
    	cmd.iterate(selection, "stored.atomnumber = stored.atomnumber + 1")
    	print("Your selection: %s" % selection)
    	print("sum of B factors: %s" % stored.tempfactor)
    	print("number of atoms: %s" % stored.atomnumber)
    	averagetempfactor = stored.tempfactor / stored.atomnumber
    	print("average B of '%s': %s" % (selection, averagetempfactor))
    cmd.extend("average_b", average_b)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Average_b&oldid=13555](https://pymolwiki.org/index.php?title=Average_b&oldid=13555)"


---

## BbPlane

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/bbPlane.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/bbPlane.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Blaine Bell](/index.php?title=User:Bell&action=edit&redlink=1 "User:Bell \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO plane between the backbone atoms of two neighboring residues. This is to show the planarity of the atoms. The image style this is meant to represent can be found many places, like "Introduction to Protein Structure" by Branden and Tooze (2nd ed. pp. 8). 

  * [![Close up of planar atoms](/images/8/8d/BbPlane3.png)](/index.php/File:BbPlane3.png "Close up of planar atoms")

Close up of planar atoms 

  * [![A few more](/images/8/89/BbPlane1.png)](/index.php/File:BbPlane1.png "A few more")

A few more 

  * [![Global view](/images/b/b8/BbPlane2.png)](/index.php/File:BbPlane2.png "Global view")

Global view 




# Examples
    
    
    # download the source and save as bbPlane.py
    run bbPlane.py
    fetch 1cll
    # make planes for residues 4-9
    bbPlane i. 4-10
    

Retrieved from "[https://pymolwiki.org/index.php?title=BbPlane&oldid=13888](https://pymolwiki.org/index.php?title=BbPlane&oldid=13888)"


---

## Cell color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Related settings
  * 4 See Also



## Overview

set **cell_color** allows one to color [cells](/index.php/Cell "Cell") independently from the rest of the representations. 

## Syntax
    
    
    set cell_color, theColor
    

where _theColor_ can be : 

  * any usual colors (blue, yellow, grey50,...)
  * number-coded colors (1:black, 2:blue, 3:greenish, ...)
  * special code -1 to revert to original chameleon setting (_set cartoon_color,-1_)



## Related settings

[sphere_color](/index.php/Sphere_color "Sphere color") [cartoon_color](/index.php/Cartoon_color "Cartoon color")

## See Also

[Color](/index.php/Color "Color")

Retrieved from "[https://pymolwiki.org/index.php?title=Cell_color&oldid=13573](https://pymolwiki.org/index.php?title=Cell_color&oldid=13573)"


---

## Center of mass

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/center_of_mass.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/center_of_mass.py)  
Author(s)  | [Henschel](/index.php?title=User:Henschel&action=edit&redlink=1 "User:Henschel \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Examples
  * 4 See Also



### Description

This script calculates the true center-of-mass (COM) or the center-of-geometry (COG) for a given selection and returns the x, y, z values in the form of a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") (rather than a CGO sphere). The benefit of using a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") is that it can be selected and used in calculations. In addition, this script also iteratively goes through all states of a selection if more than one state exists and appends the corresponding COM/COG values as states into the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). The script itself is quite simple yet robust enough to be applied in different settings. As well, the calculation of the COM/COG is handled independently from the formation of the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") and can be called as an independent function where applicable. 

### Usage
    
    
    com selection [,state=None [,mass=None [,object=None]]]
    

  


### Examples
    
    
    import center_of_mass
    fetch 1c3y, finish=1, multiplex=0
    
    com 1c3y, state=1
    #Create a pseudoatom representing the 1c3y COG and store it as "1c3y_COM"
    #The "1c3y_COM" object will contain 1 state only
    
    com 1c3y, state=1, object=COG
    #Create a pseudoatom representing the 1c3y COG and store it as "COG"
    #The "COG" object will contain 1 state only
    
    com 1c3y, state=1, object=COM, mass=1
    #Create a pseudoatom representing the 1c3y COM and store it as "COM"
    #The "COM" object will contain 1 state only
    
    com 1c3y, object=COM, mass=1
    #Create a single pseudoatom containing the COM for each state found in 1c3y and store it as "COM"
    #The "COM" object will contain MULTIPLE states!
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [COM](/index.php/COM "COM")
  * [get_extent](/index.php/Get_extent "Get extent")
  * [get_position](/index.php/Get_position "Get position")



Retrieved from "[https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892](https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892)"


---

## Centroid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/centroid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/centroid.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
    * 3.1 See Also



## Overview

Centroid is a small script that returns the value of the geometric center (or centroid) of your selection. It also can translate the object of that selection to the origin. 

## Syntax
    
    
    centroid (selection=PyMOLSelection), [center=boolean]
    

## Examples
    
    
    # get the centroid of the polymer
    import centroid
    fetch 4ins, async=0
    centroid polymer
    
    # move some 'ligand' to the origin
    centroid ligand, center=1
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [Center_Of_Mass](/index.php/Center_Of_Mass "Center Of Mass")



Retrieved from "[https://pymolwiki.org/index.php?title=Centroid&oldid=13893](https://pymolwiki.org/index.php?title=Centroid&oldid=13893)"


---

## Cgo grid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/cgo_grid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_grid.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![cgo_grid creates flowing mesh objects](/images/2/27/Cgo_grid.gif)](/index.php/File:Cgo_grid.gif "cgo_grid creates flowing mesh objects")

## Contents

  * 1 About cgo_grid
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 SEE ALSO



## About cgo_grid

**cgo_grid** will generate a flowing mesh object using the points provided or the current view. By default is will generate a flowing membrane. The shape is affected substantially by the arguments!  


  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage

**cgo_grid** has many arguments, but not all need to be set necessarily (see arguments or examples). 
    
    
    cgo_grid [ pos1 [, pos2 [, pos3 [, length_x [, length_z [, npoints_x [, npoints_z
    [, nwaves_x [, nwaves_z [, offset_x [, offset_z [, gain_x [, gain_z
    [, thickness [, color [, nstates [, startframe [, endframe
    [, mode [, view [, name [, quiet ]]]]]]]]]]]]]]]]]]]]]]
    

|   
---|---  
  
  


## Arguments
    
    
        pos1 = single atom selection (='pk1') or list of 3 floats {default: [0,0,0]}
    
        pos2 = single atom selection (='pk2') or list of 3 floats {default: [1,0,0]}
    
        pos3 = single atom selection (='pk3') or list of 3 floats {default: [0,0,1]}
    
        --> the plane is defined by pos1 (origin) and vectors to pos2 and pos3, respectively
    
        length_x = <float>: length of membrane {default: 30}
        length_z = <float>: length of membrane {default: ''} # same as length_x
    
        npoints_x = <int>: number of points(lines) along x-direction
                    {default: ''} #will be set to give a ~1 unit grid
        npoints_z = <int>: number of points(lines) along z-direction
                    {default: ''} #will be set to give a ~1 unit grid
                    {minimum: 1 # automatic}
    
        nwaves_x =   <float>: number of complete sin waves along object x-axis
                     {default: 2}
        nwaves_z =  <float>: number of complete sin waves along object z-axis
                    {default: ''} # same as nwaves_x
                    define separately to adjust number of waves in each direction
    
        offset_x = <float> phase delay (in degrees) of sin wave in x-axis
                 can be set to affect shape and starting amplitude {default: 0}
        offset_z = <float> phase delay (in degrees) of sin wave in z-axis
                 can be set to affect shape and starting amplitude
                 {default: ''} # same as  offset_x
        offset_x and offset_z can be used together to phase
        otherwise identical objects
    
        gain_x = <float>: multiplication factor for y-amplitude for x-direction
                 {default: 1}
        gain_z = <float>: multiplication factor for y-amplitude for z-direction
                 {default: ''} #=gain_x
    
        thickness = <float>: line thickness {default: 2}
    
        color = color name <string> (e.g. 'skyblue') OR
                rgb-value list of 3 floats (e.g. [1.0,1.0,1.0]) OR
                {default: ''} // opposite of background
                input illegal values for random coloring
    
        nstates =  <int>: number of states; {default: 60}
                   this setting will define how many states
                   the object will have (per wave) and how fluent and fast the
                   animation will be.
                   Higher values will promote 'fluent' transitions,
                   but decrease flow speed.
                       Note: Frame animation cycles thought the states one at a time
                       and needs to be set accordingly. Can also be used to phase
                       otherwise identical objects.
                   Set to 1 for static object {automatic minimum}
    
        startframe: specify starting frame <int> or set (='') to use current frame
                    set to 'append' to extend movie from the last frame {default: 1}
          endframe: specify end frame <int> or set (='') to use last frame
                    if 'append' is used for startframe,
                    endframe becomes the number of frames to be appended instead
                    {default: 1}
                    Note: if start- and endframe are the same, movie animation will
                    be skipped, the object will be loaded and can be used afterwards
    
        mode: defines positioning {default: 0}:
        0: pos1 is center
        1: pos1 is corner
    
        view {default: 0}:
        '0': off/ uses provided points to create CGO
        '1': overrides atom selections and uses current orienatation for positioning
             - pos1 = origin/center
             - pos2 = origin +1 in camera y
             - pos3 = origin +1 in camera z
    
        name: <string> name of cgo object {default: ''} / automatic
    
        quiet: <boolean> toggles output
    

| Concept sketch:  
[![cgo_grid concept sketch](/images/d/d1/Cgo_grid.png)](/index.php/File:Cgo_grid.png "cgo_grid concept sketch")  
---|---  
  
  


## Examples

The behaviour or shape of the cgo_grid object are substantially influenced by the arguments 
    
    
    delete all
    set_view (\
         0.263772517,   -0.113038681,    0.957937598,\
        -0.040853567,    0.990910411,    0.128179103,\
        -0.963716805,   -0.072944686,    0.256756991,\
         0.000000000,    0.000000000, -131.816467285,\
         0.000000000,    0.000000000,    0.000000000,\
       -50.008331299,  353.641235352,  -20.000000000 )
    
    #membrane
    cgo_grid color=blue
    
    #swimming worm, random color
    cgo_grid \
    pos1=[0,-5,0], pos2=[1,-5,1], pos3=[0,-5,1],\
    length_x=15,\
    npoints_z=1,\
    gain_x=2,\
    gain_z=0,\
    thickness=20,\
    color=3,\
    mode=1,\
    name="worm"
    
    #Moving Ladder
    cgo_grid \
    length_x=15,\
    pos1=[0,10,0], pos2=[0,10,1], pos3=[0,9,0],\
    npoints_x=2, npoints_z=30,\
    name="ladder"
    
    #Roof
    cgo_grid \
    nstates=1,\
    npoints_x=15,\
    npoints_z=15,\
    gain_x=20,\
    gain_z=20,\
    nwaves_x=0.5,\
    thickness=5,\
    color=cyan,\
    name="roof"
    
    #Boxes
    cgo_grid \
    pos1=[0,-10,0], pos2=[1,-10,0], pos3=[0,-10,1],\
    nstates=1,\
    npoints_x=50,\
    npoints_z=50,\
    nwaves_x=0,\
    color=[0.00 , 0.53 , 0.22],\
    thickness=5,\
    name="bottom"
    
    cgo_grid \
    nstates=1,\
    npoints_x=2,\
    npoints_z=2,\
    nwaves_x=0,\
    color=gray60,\
    thickness=10,\
    name="top"
    
    cgo_grid \
    pos1=[-15,-10,15], pos2=[-14,-10,15], pos3=[-15,-9,15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=-2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="front"
    
    cgo_grid \
    pos1=[-15,-10,-15], pos2=[-14,-10,-15], pos3=[-15,-9,-15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="back"
    
    set ray_trace_frames, 0
    set movie_loop,1
    mplay
    
    # play around with the ARGUMENTS! :-)
    

|   
---|---  
  
  


## SEE ALSO

  * [CgoCircle](/index.php/CgoCircle "CgoCircle")



Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938](https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938)"


---

## Color cbcobj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [ Oscar Conchillo-Solé](/index.php?title=User:OscarConchillo-Sol%C3%A9&action=edit&redlink=1 "User:OscarConchillo-Solé \(page does not exist\)")  
License  |   
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 Code



# Overview

color_cbcobj 

Color all chains in different objects with a different color 

color_cbcobj(selection='(all)', first_color=0, quiet=0, _self=cmd): 

uses colors and order as in _color_cycle defined here: <https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/util.py>

quiet=0 means verbose (not quiet) 

# Usage
    
    
    color_cbcobj selection, [first_color=0..n, [quiet=0/1 ]]
    

where: 

  * **selection** can be any existing or newly-defined selection
  * **first_color** first color tu use from "_color_cycle" (first is 0)
  * **quiet** (default 0) print colored chains and models quiet=0 means verbose (not quiet)



# Examples
    
    
    PyMOL>color_cbcobj all and elem C,2,1
    

colors all chains in all objects differently, only carbon atoms, starting by color 154 (lightmagenta), does not show what has been colored 
    
    
    PyMOL>color_cbcobj all and elem C,2
     util.cbcobj: color 154 and model 6Y6C_BC', (chain B)
     Executive: Colored 670 atoms.
     util.cbcobj: color 6 and model 6Y6C_BC', (chain C)
     Executive: Colored 1131 atoms.
     util.cbcobj: color 9 and model 6Y6C_AD', (chain A)
     Executive: Colored 664 atoms.
     util.cbcobj: color 29 and model 6Y6C_AD', (chain D)
     Executive: Colored 1135 atoms.
    

# Code

Copy the following text and save it as color_cbcobj.py 
    
    
    from pymol import cmd, CmdException
    
    _color_cycle = [
        26   , # /* carbon */
        5    , # /* cyan */
        154  , # /* lightmagenta */
        6    , # /* yellow */
        9    , # /* salmon */
        29   , # /* hydrogen */
        11   , # /* slate */
        13   , # /* orange */
        10   , # /* lime */
        5262 , # /* deepteal */
        12   , # /* hotpink */
        36   , # /* yelloworange */
        5271 , # /* violetpurple */
        124  , # /* grey70 */
        17   , # /* marine */
        18   , # /* olive */
        5270 , # /* smudge */
        20   , # /* teal */
        5272 , # /* dirtyviolet */
        52   , # /* wheat */
        5258 , # /* deepsalmon */
        5274 , # /* lightpink */
        5257 , # /* aquamarine */
        5256 , # /* paleyellow */
        15   , # /* limegreen */
        5277 , # /* skyblue */
        5279 , # /* warmpink */
        5276 , # /* limon */
        53   , # /* violet */
        5278 , # /* bluewhite */
        5275 , # /* greencyan */
        5269 , # /* sand */
        22   , # /* forest */
        5266 , # /* lightteal */
        5280 , # /* darksalmon */
        5267 , # /* splitpea */
        5268 , # /* raspberry */
        104  , # /* grey50 */
        23   , # /* deepblue */
        51   , # /* brown */
        ]
    
    _color_cycle_len = len(_color_cycle)
    
    def color_cbcobj(selection='(all)', first_color=0, quiet=0, _self=cmd):
        '''
    DESCRIPTION
    
        Color all chains in different objects a different color
    
    SEE ALSO
    
        util.cbc, https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/util.py
        split_chains.py, https://pymolwiki.org/index.php/Split_chains
        '''
        quiet = int(quiet)
        count = int(first_color)
        models = cmd.get_object_list('(' + selection + ')')
        for model in models:
            for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
                if len(chain.split()) != 1:
                    chain = '"' + chain + '"'
                color = _color_cycle[count % _color_cycle_len]
                if not quiet:
                    print(" util.cbcobj: color %s and model %s', (chain %s)" % (color, model, chain))
                #_self.color(color, "(chain %s and (%s))" % (chain, model), quiet=quiet)
                _self.color(color, "(chain %s and (%s) and (%s))" % (chain, model, selection), quiet)
                count+=1
    
    cmd.extend('color_cbcobj', color_cbcobj)
    
    #color_cbcobj all and elem C,1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_cbcobj&oldid=13496](https://pymolwiki.org/index.php?title=Color_cbcobj&oldid=13496)"


---

## ColorByRMSD

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/colorbyrmsd.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorbyrmsd.py)  
Author(s)  | [Shivender Shandilya](/index.php/User:Shiven "User:Shiven"), [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Arguments
  * 4 Examples



## Introduction

This script allows you to color two structures by Root Mean Square Deviation (RMSD). The distances between aligned C-alpha atom pairs are stored as B-factors of these residues, which are colored by a color spectrum, with blue specifying the minimum pairwise RMSD and red indicating the maximum. Unaligned residues are colored gray. 

## Usage
    
    
    colorbyrmsd mobile, target [, doAlign [, doPretty [, guide [, method ]]]]
    

## Arguments

  * **mobile** = string: atom selection for mobile atoms
  * **target** = string: atom selection for target atoms
  * **doAlign** = 0 or 1: Superpose selections before calculating distances {default: 1}
  * **doPretty** = 0 or 1: Show nice representation and colors {default: 1}
  * **guide** = 0 or 1: Only use C-alpha atoms {default: 1}
  * **method** = align or super: Method to match atoms {default: super}



## Examples
    
    
    # example #1
    colorbyrmsd 1cbs, 1hmt, doAlign=1, doPretty=1
    # example #2
    colorbyrmsd 1eaz, 1fao, doAlign=1, doPretty=1
    

  * [![1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.](/images/5/57/ColorByRMSD_1cbs_1hmt.png)](/index.php/File:ColorByRMSD_1cbs_1hmt.png "1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.")

1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white. 

  * [![1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.](/images/b/b6/ColorByRMSD_1eaz_1fao.png)](/index.php/File:ColorByRMSD_1eaz_1fao.png "1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.")

1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white. 




Retrieved from "[https://pymolwiki.org/index.php?title=ColorByRMSD&oldid=13898](https://pymolwiki.org/index.php?title=ColorByRMSD&oldid=13898)"


---

## Declare Command

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## The problem

Current PyMOL approach to new plugin commands is outdated. 

## The proposal

Introduce a new system based on modern Python with type checking and implicit type coercion. 

## What works right now?

There's some support on PyMOL open-source, but not on Incentive. However it isn't working for all possible cases. 

  

    
    
    @cmd.declare_command
    def new_command(
        my_var: int | float,
        dirname: Path = '.',
        # Optional values aren't currently supported.
        # Tuples are.
        nullable_point: Optional[Tuple[int, int, int]] = None,
        extended_calculation: bool = True,
        old_style: Any = "Support anything currently not supported",
    ):
        """
        A cool docstring.
        """
        pass
    

Retrieved from "[https://pymolwiki.org/index.php?title=Declare_Command&oldid=13883](https://pymolwiki.org/index.php?title=Declare_Command&oldid=13883)"


---

## Delete states

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**delete_states** removes states from a multi-state object like a trajectory. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    delete_states name, 1 2 3  # delete states 1, 2, and 3
    delete_states name, 1-3 9-13  # delete states 1 through 3 and 9 through 13
    

name = name of object or name expression (wildcard supported) 

### PYMOL API
    
    
    delete_states(name: str, states: str) -> None
    
    
    
    cmd.delete_states(string name = object-name, string states = states string)
    

### EXAMPLES
    
    
       delete_state 1nmr, 1-5     # delete states 1 to 5 from 1nmr
       delete_state *, 1-3 10-40  # deletes states 1 to 3 and 10 to 40 from all applicable objects
    

Note that special care needs to be taken to ensure that the object or selection name does not contain any quotes when passed as an argument. 

Note This function currently only applies to non-discrete multistate molecular objects. 

### SEE ALSO

[Remove](/index.php/Remove "Remove") [Delete](/index.php/Delete "Delete")

Retrieved from "[https://pymolwiki.org/index.php?title=Delete_states&oldid=13577](https://pymolwiki.org/index.php?title=Delete_states&oldid=13577)"


---

## Distance

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Distance
  * 2 Usage
  * 3 PYMOL API
    * 3.1 EXAMPLES
  * 4 See Also



# Distance

**distance** creates a new distance object between two selections. It will display all distances within the cutoff. Distance is also used to make hydrogen bonds. Calling distance without arguments will show distances between selections (pk1) and (pk1), which can be set in editing mode or using the PkAt mouse action (usually CTRL-middle-click).  
  
Note: For interactive use, the **measurement** [wizard](/index.php/Wizard "Wizard") (from the PyMOL menu) makes measuring distances easier than using the distance command.  
If you want to measure distance and avoid creating a distance object, use [Get Distance](/index.php/Get_Distance "Get Distance") or [Distancetoatom](/index.php/Distancetoatom "Distancetoatom").  
|  [![Example of distance. See example #1.](/images/e/ea/Dist_ex1.png)](/index.php/File:Dist_ex1.png "Example of distance. See example #1.") |  [![Example of distance. See example #2.](/images/4/4e/Dist_ex2.png)](/index.php/File:Dist_ex2.png "Example of distance. See example #2.") |   
---|---|---|---  
  
# Usage
    
    
    distance [ name [, selection1 [, selection2 [, cutoff [, mode ]]]]]
    

**name**

    

    string: name of the distance object to create

**selection1**

    

    string: first atom selection

**selection2**

    

    string: second atom selection

**cutoff**

    

    float: longest distance to show

**mode**

    

    0: all interatomic distances
    1: only bond distances
    2: only show polar contact distances
    3: like mode=0, but use [distance_exclusion](/index.php?title=Distance_exclusion&action=edit&redlink=1 "Distance exclusion \(page does not exist\)") setting
    4: distance between centroids (_new in 1.8.2_)
    5: pi-pi and pi-cation interactions
    6: pi-pi interactions
    7: pi-cation interactions
    8: like mode=3, but cutoff is the ratio between distance and sum of VDW radii

    

    modes 5 to 8 are available only in incentive PyMOL.

# PYMOL API
    
    
    cmd.distance( string name, string selection1, string selection2,
                  float cutoff, int mode )
       # returns the average distance between all atoms/frames
    

### EXAMPLES

  * Get and show the distance from residue 10's alpha carbon to residue 40's alpha carbon in 1ESR:


    
    
    # make the first residue 0.
    zero_residues 1esr, 0
    distance i. 10 and n . CA, i. 40 and n. CA
    

  * Get and show the distance from residue 10's alpha carbon to residue 35-42's alpha carbon in 1ESR:


    
    
    # make the first residue 0.
    zero_residues 1esr, 0
    distance i. 10 and n . CA, i. 35-42 and n. CA
    

  * This neat example shows how to create distance measurements from an atom in a molecule to all other atoms in the molecule (since PyMol supports wildcards).


    
    
     cmd.distance("(/mol1///1/C)","(/mol1///2/C*)")
    

or written without the PyMolScript code, 
    
    
    dist /mol1///1/C, /mol1///2/C*
    

  * Create multiple distance objects


    
    
    for at1 in cmd.index("resi 10"): \
       for at2 in cmd.index("resi 11"): \
           cmd.distance(None, "%s`%d"%at1, "%s`%d"%at2)
    
    
    
    distance (selection1), (selection2)
    
    # example
    dist i. 158 and n. CA, i. 160 and n. CA 
    
    distance mydist, 14/CA, 29/CA
    distance hbonds, all, all, 3.2, mode=2
    

# See Also

  * [Get Distance](/index.php/Get_Distance "Get Distance") # Avoid creating an object
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom") # Automated script for distances to a point
  * [Measure_Distance](/index.php/Measure_Distance "Measure Distance") # basic script
  * [Lines](/index.php/Lines "Lines")



Retrieved from "[https://pymolwiki.org/index.php?title=Distance&oldid=13867](https://pymolwiki.org/index.php?title=Distance&oldid=13867)"


---

## Distancetoatom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/distancetoatom.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/distancetoatom.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar") and [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![get distance to atoms within cutoff](/images/c/c2/Distancetoatom_0.png)](/index.php/File:Distancetoatom_0.png "get distance to atoms within cutoff")

## Contents

  * 1 About distancetoatom
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
  * 4 SEE ALSO



## About distancetoatom

**distancetoatom** prints all distances between a specified atom, coordinate or group selection center and all atoms within cutoff distance that are part of the selection.  
All coordinates and distances can be saved in a csv-style text file report and/or can be stored in a (custom) atom property, if defined. 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage
    
    
        distancetoatom [ origin [, cutoff [, filename [, selection [, state [, property_name [, coordinates [, decimals [, sort [, quiet ]]]]]]]]]]
    

|   
---|---  
  
### Arguments

Arguments for distancetoatom   
---  
Keyword  | Default  | What it does   
origin  | pk1  | "origin" defines the coordinates for the origin and can be: 1\. a list with coordinates [x,y,z]  
2\. a single atom selection string  
3\. a multi-atom selection string (center will be used)   
cutoff  | 10  | "cutoff" defines the maximum distance, atoms further away are not considered   
filename  | None  | "filename" is the name of a report file that will be created (optional). set to e.g. 'report.txt' to create a report. This file is CSV style and can be imported into EXCEL.  
(omit or set to "", None, 0 or False to disable)   
selection  | all  | only atoms within "selection" (and cutoff distance) will be used for calculation   
state  | 0  | "state" is the state that will be used for calculations use 0 (omit), to automatically use the current state   
property_name  | p.dist  | "property_name" defines a (custom) property in which the calculated distance will be stored. This can be used e.g. for labeling or spectrum coloring etc. (cf. examples)  
NB! older PyMOL versions only support b or q and the distance will overwrite this property if set.   
coordinates  | 0  | "coordinates" toggles whether, besides distance, the atom coordinates will be included in the report.   
decimals  | 3  | "decimals" is the number of decimal places calculated. Note that PDB files do not support a higher resolution than 3.   
sort  | 1  | "sort" defines how the output will be sorted: 1: ascending (default)  
0: no sorting (by selection)  
-1: descending   
quiet  | 1  | toggle verbosity   
  
  


## Examples
    
    
    # make function available to PyMOL
    import distancetoatom
    
    ##############################
    # Example 1: print all distances
    frag LYS
    edit name CA
    distancetoatom
    
    ##############################
    # Example 2: print all distances to file
    fetch 1cll, async=0
    distancetoatom origin=/1cll//A/CA`150/CA, filename=distances.txt, selection=elem O, property_name=b
    # the file is CSV-format and can be imported, e.g. to EXCEL
    
    # format
    hide everything
    select byres (resi 150 expand 5 and not resn hoh) 
    show_as sticks, sele
    util.cbaw sele
    show_as spheres, resi 150
    zoom sele
    
    ##########
    # Label by stored distance (property_name)
    label sele and elem O, "%.2f"%b
    
    
    ##############################
    # Example 3: color an object by distance to a coordinate
    
    fetch 1cll, async=0
    distancetoatom [0,0,0], cutoff=1000, property_name=b
    # the distance is stored in the b-factor in this case
    # newer PyMOL versions support custom properties, e.g. "p.dist"
    spectrum b, rainbow_rev
    

|  [![Example 2](/images/5/5c/Distancetoatom_2.png)](/index.php/File:Distancetoatom_2.png "Example 2")  
  
[![Example 3](/images/3/3a/Distancetoatom_1.png)](/index.php/File:Distancetoatom_1.png "Example 3")  
---|---  
  
  


## SEE ALSO

  * [Distance](/index.php/Distance "Distance")
  * [Get Distance](/index.php/Get_Distance "Get Distance")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")



Retrieved from "[https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939](https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939)"


---

## Draw Protein Dimensions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/Draw_Protein_Dimensions.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/Draw_Protein_Dimensions.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw the dimensions of a protein based on an Inertia Axis Aligned Bounding Box (IABB). 

  
The idea behind this script is to calculate an approximate minimal bounding box (MBB) to extract the cell dimensions of a protein. To calculate the MBB is not trivial and usually the Axis Aligned Bounding Box (AABB) does not show up the real dimensions of the protein. This script calculates the inertia tensor of the object, extract the eigenvectors and use them to rotate the molecule (using as rotation matrix the transpose of the eigenvector matrix). The result is a molecule oriented with the inertia axis aligned with the cartesian axis. A new Bounding Box is calculated, which is called Inertia Axis Aligned Bounding Box (IABB), with a volume always lower than the AABB volume, and in many cases may correspond with the MBB. 

As always with these type of things, you have to use at your own risk. I did not try all the possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado (at) gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run path-to-the-script/Draw_Protein_Dimensions.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_Protein_Dimensions _selection_**

This will draw the cell dimensions of your selection based on a IABB box. It also generates the IABB box and the inertia axis, you just need to do "show cgo" to display them. 

You could also try: 

**draw_BB _selection_**

This will draw the AABB and IABB boxes with their cell dimensions and show up their volumes, you can compare them. 

  


# Examples
    
    
    # download the source and save as Draw_Protein_Dimensions.py
    run Draw_Protein_Dimensions.py
    fetch 2vak
    # calculate the dimensions of the full ASU
    draw_Protein_Dimensions 2vak
    

  * [![Dimensions of 2vak ASU](/images/7/77/2vak_ASU.png)](/index.php/File:2vak_ASU.png "Dimensions of 2vak ASU")

Dimensions of 2vak ASU 



    
    
    # you can extract only one chain and calculates the dimensions
    draw_Protein_Dimensions obj01
    

  * [![Dimensions of one protomer \(chain A - 2vak\)](/images/e/ed/2vak_A.png)](/index.php/File:2vak_A.png "Dimensions of one protomer \(chain A - 2vak\)")

Dimensions of one protomer (chain A - 2vak) 



    
    
    # You can also draw the Bounding boxes (AABB and IABB) to compare them.
    fetch 2vak
    draw_BB 2vak
    

  


  * [![Axis-aligned bounding box](/images/4/46/2vak_AABB.png)](/index.php/File:2vak_AABB.png "Axis-aligned bounding box")

Axis-aligned bounding box 

  * [![Inertia-axis-aligned bounding box](/images/4/41/2vak_IABB.png)](/index.php/File:2vak_IABB.png "Inertia-axis-aligned bounding box")

Inertia-axis-aligned bounding box 




Retrieved from "[https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902](https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902)"


---

## DrawGridBox

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/drawgridbox.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/drawgridbox.py)  
Author(s)  | [Cunliang Geng](/index.php?title=User:Clgeng&action=edit&redlink=1 "User:Clgeng \(page does not exist\)")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the [drawgridbox](/index.php?title=Drawgridbox&action=edit&redlink=1 "Drawgridbox \(page does not exist\)") command, which draw grid boxes around a selection. 

## Usage
    
    
     drawgridbox [ selection,  [nx, [ny, [nz, [padding, [lw, [r, [g, [b]]]]]]]]]
    

  


## Prameters

  * **selection** = str: atom selection {default: all}
  * **nx** = int: number of grids on axis X {default: 10}
  * **ny** = int: number of grids on axis Y {default: 10}
  * **nz** = int: number of grids on axis Z {default: 10}
  * **padding** = float: default 0.0
  * **lw** = float: line width {default: 2.0}
  * **r** = float: red color component, valid range is [0.0, 1.0] {default 1.0}
  * **g** = float: green color component, valid range is [0.0, 1.0] {default 1.0}
  * **b** = float: blue color component, valid range is [0.0, 1.0] {default 1.0}



  


## Examples
    
    
    load drawgridbox.py
    fetch 1cbh
    show surface
    drawgridbox 1cbh, nx=5, ny=5, nz=5,  lw=1, g=0, b=0
    

[![Drawgridbox.png](/images/5/5a/Drawgridbox.png)](/index.php/File:Drawgridbox.png)

  

    
    
    drawgridbox 1cbh, nx=5, ny=5, nz=5, padding=5, lw=1, r=0, g=0
    

[![Drawgridbox padding.png](/images/3/32/Drawgridbox_padding.png)](/index.php/File:Drawgridbox_padding.png)

Retrieved from "[https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904](https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904)"


---

## Dssr block

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [script/dssr_block.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/script/dssr_block.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![Trna.png](/images/d/dd/Trna.png)](/index.php/File:Trna.png)

[](/index.php/File:Trna.png "Enlarge")

This script adds the dssr_block command, which is a simple wrapper for the **DSSR** program and can create "block" shaped cartoons for nucleic acid bases and base pairs. 

_Requires the**x3dna-dssr** executable, obtainable from <http://x3dna.org>_

Instead of installing **x3dna-dssr** and the plugin locally, you can also use the web server at <http://skmatic.x3dna.org/>

See also the blog post by DSSR author Xiang-Jun Lu: <http://x3dna.org/highlights/dssr-base-blocks-in-pymol-interactively>

## Contents

  * 1 Recommended Setup
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 See Also



## Recommended Setup

Place the **x3dna-dssr** executable into **$PATH**. Examples: 

  * Linux/Mac: **/usr/bin/x3dna-dssr**
  * Windows: **C:\Program Files\PyMOL\PyMOL\x3dna-dssr.exe**



The script can be [run](/index.php/Run "Run"), imported as a Python module, or installed with the [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager"). 

## Usage
    
    
    dssr_block [ selection [, state [, block_file [, block_depth
           [, block_color [, name [, exe ]]]]]]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **state** = int: object state (0 for all states) {default: -1, current state}
  * **block_file** = face|edge|wc|equal|minor|gray: Corresponds to the **\--block-file** option (see [DSSR manual](http://x3dna.org/)). Values can be combined, e.g. "wc-minor". {default: face}
  * **block_depth** = float: thickness of rectangular blocks {default: 0.5}
  * **block_color** = str: Corresponds to the **\--block-color** option (new in DSSR v1.5.2) {default: }
  * **name** = str: name of new CGO object {default: dssr_block##}
  * **exe** = str: path to "x3dna-dssr" executable {default: x3dna-dssr}



## Examples

Combining DSSR block representation with regular PyMOL cartoons: 
    
    
    fetch 1ehz, async=0
    as cartoon
    set cartoon_ladder_radius, 0.1
    set cartoon_ladder_color, gray
    set cartoon_nucleic_acid_mode, 1
    dssr_block
    

Joined base-pair blocks (block_file=wc): 
    
    
    fetch 1ehz, async=0
    dssr_block block_file=wc
    

Multi-state Example: 
    
    
    fetch 2n2d, async=0
    dssr_block 2n2d, 0
    set all_states
    

Custom coloring: 
    
    
    fetch 1msy, async=0
    dssr_block block_color=N red | minor 0.9 | major yellow
    

## See Also

  * [3DNA](/index.php/3DNA "3DNA") (old, obsoleted by x3dna-dssr)
  * [Overview of nucleic acid cartoons](/index.php/Overview_of_nucleic_acid_cartoons "Overview of nucleic acid cartoons")
  * <http://skmatic.x3dna.org/>



Retrieved from "[https://pymolwiki.org/index.php?title=Dssr_block&oldid=13940](https://pymolwiki.org/index.php?title=Dssr_block&oldid=13940)"


---

## Example Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Multiple Object Manipulation Scripts
    * 1.1 Dis/Enable Objects
    * 1.2 Comments in Scripts
      * 1.2.1 Sources



## Multiple Object Manipulation Scripts

### Dis/Enable Objects

If someone has many (possibly hundreds) of objects -- say distance objects -- one can turn them all on or off by using Python scripts within PyMol. The following script will extend the commands "sd" and "hd" for "show distance" and "hide distance," respectively. The first script, 
    
    
    from pymol import cmd 
    num_dist = 100 
            
    def show_dist(): 
        """ show all of my distance objects """ 
        for i in range(num_dist): 
            cmd.enable('_dist%s'%i) 
            
    def hide_dist(): 
        """ hide all of my distance objects """ 
        for i in range(num_dist): 
            cmd.disable('_dist%s'%i) 
            
    cmd.extend('sd',show_dist) 
    cmd.extend('hd',hide_dist)
    

works on 100 objects. We can extend the idea with more elegant scripting to work w/o forcing us to keep track of the number of objects: 
    
    
    def show_dist(): 
        dists = [name for name in cmd.get_names() if cmd.get_type(name) == 'object:distance'] 
        for name in dists: cmd.enable(name)
    

Now, just running "sd" would show all the distance objects. 

  


### Comments in Scripts

The hash-mark, "#" is the Python comment. However, when scripting in Python for PyMol note that the hash character (#) works to include comments on separate lines, but not at the end of a line that contains commands. So you can do 
    
    
    # Create separate dimer
    create dimer,(chain A,B)
    

but not 
    
    
    create dimer,(chain A,B)  # Create separate dimer
    

instead use: 
    
    
    create dimer,(chain A,B);  # Create separate dimer
    

#### Sources

Taken from the PyMol Users list. Python source by Michael Lerner. 

Retrieved from "[https://pymolwiki.org/index.php?title=Example_Scripts&oldid=13570](https://pymolwiki.org/index.php?title=Example_Scripts&oldid=13570)"


---

## Fab

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**fab** builds peptide entities from sequence. The sequence must be specified in one-letter code. Several fragments will be created if the sequence contains spaces (chain breaks). To specify chain-id and residue number, an advanced syntax pattern like _chain/resi/_ can be put before each sequence fragment. 

Similar functionality is also provided by the graphical [Builder](/index.php/Builder "Builder") in "Protein" mode. 

_New in PyMOL version 1.2_

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLE
  * 4 SEE ALSO



## USAGE
    
    
    fab input [, name [, mode [, resi [, chain [, segi [, state [, dir [, hydro [, ss [, async ]]]]]]]]]]
    

## ARGUMENTS

  * input = string: Amino acid sequence in one-letter code


  * name = string: name of object to create or modify {default: obj??}


  * mode = string: Only supported mode is "peptide" {default: peptide}


  * resi = integer: Residue number to start numbering from {default: 1}


  * chain, segi = string: Chain id and segment id to assign


  * state = integer: {default: -1}


  * dir = 0/1: 0=append to N-terminus, 1=append to C-terminus {default: 1}


  * hydro = 0/1: With or without hydrogens (BROKEN! Use [auto_remove_hydrogens](/index.php/Settings "Settings")}


  * ss = int: Secondary structure 1=alpha helix, 2=antiparallel beta, 3=parallel beta, 4=flat {default: 0}



## EXAMPLE

The two examples produce the same result: 
    
    
    # use resi and chain attributes
    fab KVRISAEL, myprot1, resi=10, chain=B, ss=2
    
    # use advanced syntax
    fab B/10/ KVRISAEL, myprot2, ss=2
    

## SEE ALSO

  * [Builder](/index.php/Builder "Builder")
  * [Peptide Sequence](/index.php/Peptide_Sequence "Peptide Sequence")
  * [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure")
  * [[1]](https://github.com/zigeuner/robert_campbell_pymol_scripts/tree/master/work_pymol/build_seq.py) by Robert Campbell



Retrieved from "[https://pymolwiki.org/index.php?title=Fab&oldid=13844](https://pymolwiki.org/index.php?title=Fab&oldid=13844)"


---

## Findseq

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/findseq.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findseq.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview & Motivation

Anyone ever give you a protein and then say, find the sequence "FLVEW"? Well, this script will find any string or regular expression in a given object and return that selection for you. Here's an example, 
    
    
    reinitialize
    import findseq
    # fetch two sugar-binding PDB
    fetch 1tvn, async=0
    # Now, find FATEW in 1tvn, similarly
    findseq FATEW, 1tvn
    # lower-case works, too
    findseq fatew, 1tvn
    # how about a regular expression?
    findseq F.*W, 1tvn
     
    # Find the regular expression:
    #  ..H[TA]LVWH
    # in the few proteins loaded.
    # I then showed them as sticks and colored them to highlight matched AAs
    for x in cmd.get_names(): findseq.findseq("..H[TA]LVWH", x, "sele_"+x, firstOnly=1)
    

[![](/images/c/c3/SeqFinder.png)](/index.php/File:SeqFinder.png)

[](/index.php/File:SeqFinder.png "Enlarge")

Red residues were those matching the regular expression '..H[TA]LVWH'.

# Usage

I built this to be rather flexible. You call it as: 
    
    
    findseq needle, haystack[, selName[, het[, firstOnly ]]]
    

where the options are: 

    

    **needle** the sequence of amino acids to find. Should be a string of one letter amino acid abbreviations. Can also be a string-style regular expression (eg. FW.*QQ).
    **haystack** the PyMOL object or selection in which to search
    **selName** the name of the returned selection. If you leave this blank, it'll be _foundSeqXYZ_ where XYZ is some random integer (eg. foundSeq1435); if you supply _sele_ then the usual PyMOL _(sele)_ is used; and, finally, if it's anything else, then that will be used verbatim. Defaults to _foundSeqXYZ_ so as not to overwrite any selections you might have in _sele_.
    **het** 0/1 -- if 0 then heteroatoms are not considered; if 1 then they are; defaults to 0
    **firstOnly** 0/1 -- if 0 then all matches are selected and returned; if 1 then only the first is returned

# See Also

select_pepseq@[Psico](/index.php/Psico "Psico")

Retrieved from "[https://pymolwiki.org/index.php?title=Findseq&oldid=13907](https://pymolwiki.org/index.php?title=Findseq&oldid=13907)"


---

## FindSurfaceCharge

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/findSurfaceCharge.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findSurfaceCharge.py)  
Author(s)  | [Teddy Warner](/index.php/User:TeddyWarner "User:TeddyWarner")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Drawing upon the [findSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") script, the findSurfaceCharge script will identify and output a list of all charged residues on the surface of a selectionand calculates the ionization state of a protein at a given pH. The charge can be calculated for either a folded or denatured protein. This function is intended to be used to give buffer conditions for mass spectrometry. 

# Usage
    
    
       findSurfaceCharge [selection, [pH, [folded ,[cutoff]]]]
    

# Arguments

  * **selection** = str: The object or selection for which to find exposed residues {default: all}
  * **pH** = float: The pH to calculate the surface charge at {default: 7.0}
  * **folded** = bool: Whether the program should calculate the charge of a folded (True) or denatured (False) protein.
  * **cutoff** = float: The cutoff in square Angstroms that defines exposed or not. Those atoms with > cutoff Ang^2 exposed will be considered _exposed_ {default: 2.5 Ang^2}



# Examples

  * [![Result of 4FIX.pdb at pH 7.0.](/images/2/2a/SurfaceCharge.PNG)](/index.php/File:SurfaceCharge.PNG "Result of 4FIX.pdb at pH 7.0.")

Result of 4FIX.pdb at pH 7.0. 



    
    
    run findSurfaceResiduesListCharged.py
    fetch 4FIX
    
    findSurfaceResiduesListCharged
    
    # see how pH changes the protein surface charge:
    findSurfaceCharge("4fix",7.0)
        Exposed charged residues: 
            ERREDRKEE...
        The expected surface charge of 4fix at pH 7.0 is: +3.24
    
    findSurfaceCharge("4fix",7.0,False)
        Charged residues: 
            HHHHHHRHERREDRKEE...
        The expected charge of denatured 4fix at pH 7.0 is: +0.13
    
    findSurfaceCharge("4fix",10.0)
        Charged residues: ...
        The expected surface charge of 4fix at pH 10.0 is: -3.86
    

# See Also

  * [findSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [get_area](/index.php/Get_Area "Get Area")



Retrieved from "[https://pymolwiki.org/index.php?title=FindSurfaceCharge&oldid=13951](https://pymolwiki.org/index.php?title=FindSurfaceCharge&oldid=13951)"


---

## FindSurfaceResidues

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/findSurfaceResidues.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findSurfaceResidues.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The findSurfaceResidues script will select (and color if requested) surface residues and atoms on an object or selection. See the options below. 

Each time, the script will create two new selections called, **exposed_res_XYZ** and **exposed_atm_XYZ** where _XYZ_ is some random number. This is done so that no other selections/objects are overwritten. 

# Usage
    
    
    findSurfaceResidues [ selection=all [, cutoff=2.5 [, doShow=0 ]]]
    

# Arguments

  * **selection** = str: The object or selection for which to find exposed residues {default: all}
  * **cutoff** = float: The cutoff in square Angstroms that defines exposed or not. Those atoms with > cutoff Ang^2 exposed will be considered _exposed_ {default: 2.5 Ang^2}
  * **doShow** = 0/1: Change the visualization to highlight the exposed residues vs interior {default: 0}



# Examples

  * [![Result of $TUT/1hpv.pdb at 2.5 Ang cutoff.](/images/f/fe/FindExRes.png)](/index.php/File:FindExRes.png "Result of $TUT/1hpv.pdb at 2.5 Ang cutoff.")

Result of $TUT/1hpv.pdb at 2.5 Ang cutoff. 

  * [![Example coloring of surface residues](/images/8/83/Surface_residues_ex.png)](/index.php/File:Surface_residues_ex.png "Example coloring of surface residues")

Example coloring of surface residues 



    
    
    run findSurfaceResidues.py
    fetch 1hpv, async=0
    
    findSurfaceResidues
    
    # now show the exposed
    findSurfaceResidues doShow=1
    
    # watch how the visualization changes:
    findSurfaceResidues doShow=1, cutoff=0.5
    findSurfaceResidues doShow=1, cutoff=1.0
    findSurfaceResidues doShow=1, cutoff=1.5
    findSurfaceResidues doShow=1, cutoff=2.0
    findSurfaceResidues doShow=1, cutoff=2.5
    findSurfaceResidues doShow=1, cutoff=3.0
    

# See Also

  * [get_sasa_relative](/index.php/Get_sasa_relative "Get sasa relative")
  * [get_area](/index.php/Get_Area "Get Area")



Retrieved from "[https://pymolwiki.org/index.php?title=FindSurfaceResidues&oldid=13953](https://pymolwiki.org/index.php?title=FindSurfaceResidues&oldid=13953)"


---

## Format bonds

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/format_bonds.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/format_bonds.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The script **format_bonds** will automatically format bonds in amino acids. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Notes
  * 4 SEE ALSO



## Usage
    
    
    format_bonds [ selection [, bonds ]]
    

## Examples

  * [![format_bonds bonds=1](/images/1/10/PHE_valence_0.png)](/index.php/File:PHE_valence_0.png "format_bonds bonds=1")

format_bonds bonds=1 

  * [![format_bonds bonds=2](/images/9/9f/PHE_valence_1_mode_1.png)](/index.php/File:PHE_valence_1_mode_1.png "format_bonds bonds=2")

format_bonds bonds=2 

  * [![format_bonds](/images/c/ce/PHE_delocalized.png)](/index.php/File:PHE_delocalized.png "format_bonds")

format_bonds 



    
    
    import format_bonds
    
    frag PHE
    format_bonds
    
    format_bonds bonds=2
    

  


## Notes

  * Remember to correctly configure plugin import (see: [Git intro](/index.php/Git_intro "Git intro"))
  * **format_bonds** will introduce delocalized bonds by default or when _bonds_ is larger than 2.  

  * Setting _bonds=1_ will simply disable valence display (globally!)
  * The _selection_ argument is 'all' by default and can be used to restrict editing to selected residues.  

  * Note that **format_bonds** will also format acidic residues, the C-terminus, arginine and nitro groups.
  * Tip: press the TAB key after entering format_bonds to get argument suggestions



  


# SEE ALSO

[Git intro](/index.php/Git_intro "Git intro"), [Valence](/index.php/Valence "Valence"), [Modeling and Editing Structures](/index.php/Modeling_and_Editing_Structures "Modeling and Editing Structures")

Retrieved from "[https://pymolwiki.org/index.php?title=Format_bonds&oldid=13941](https://pymolwiki.org/index.php?title=Format_bonds&oldid=13941)"


---

## Forster distance calculator

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/forster_distance_calculator.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/forster_distance_calculator.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Spectre input
  * 3 Getting spectres
  * 4 Input values
  * 5 How the script works
  * 6 How to run the script



## Introduction

This script can be handsome, if one is working with Förster resonance energy transfer in proteins. The script calculates the [Förster distance](http://en.wikipedia.org/wiki/F%C3%B6rster_resonance_energy_transfer): R0, from two dyes excitation/emission spectres.   
This script is very handsome, if one want's to pick two dyes from different companies, and the Förster distance is not provided in a table. 

This script does no calculation of proteins in pymol, but is made as a "hack/shortcut" for a python script.   
We use the python part of pymol to do the calculations, so a student would not need to install python at home, but simply pymol.   


## Spectre input

**There should be provided the path to four datafiles:**   
Excitation and Emission spectre for the Donor dye: D_Exi="path.txt" D_Emi="path.txt   
Excitation and Emission spectre for the Acceptor dye: A_Exi="path.txt" A_Emi="path.txt   


Each of the files shall be a two column file. Separated by space. Numbers are "." dot noted and not by "," comma. One can for example prepare the files by search-and-replace in a small text editor.   
The first column is the wavelength in nanometers "nm" or centimetres "cm". Second column is the arbitrary units of excitation/emission. **For example:**
    
    
    Absorption
    Wavelength "Alexa Fluor 488 H2O (AB)" 
    
    300.00 0.1949100000 
    301.00 0.1991200000 
    302.00 0.2045100000 
    303.00 0.2078800000 
    ....
    

## Getting spectres

For example provides the company [ATTO-TEC](http://www.atto-tec.com/index.php?id=128&L=1&language=en) spectre for their dyes in excel format [Spectre(excel-file)](https://www.atto-tec.com/index.php?id=113&no_cache=1&L=1). These can easily be copied from excel to a flat two column file. 

Some of the most cited dyes are the Alexa Fluor dyes from the company [invitrogen](http://www.invitrogen.com/site/us/en/home/References/Molecular-Probes-The-Handbook/Fluorophores-and-Their-Amine-Reactive-Derivatives/Alexa-Fluor-Dyes-Spanning-the-Visible-and-Infrared-Spectrum.html). Invitrogen does not provide spectres for their dyes in dataformat, but as flat picture files. 

Luckily, a group on University of Arizona have traced several spectre of dyes from the literature with a graphics program. They have made these spectre easily public at <http://www.spectra.arizona.edu/>. I highly recommend this homepage. With these Spectra, the script can calculate the Forster Distance for different dyes from different companies.   


Download one spectrum at the time by "deselecting" in the right side of the graph window. Then get the datafile with the button in the lower right corner. 

## Input values

The calculation needs input for: 

  * **Qd** : Fluorescence quantum yield of the donor in the absence of the acceptor. 

    This is normally provided from the manufacture, and this information is normally also available in the [database](http://www.spectra.arizona.edu/).
    It is recognized in the literature, that the value for **Qd** , changes on behalf on which position its located on the protein.
    Changes in 5-10 % is normal, but much larger values can be expected if the dye make hydrophobic interactions with the protein. This is dye-nature dependant.
    This value can only be determined precisely with experiments for the protein at hand.
    But anyway, for a "starting" selection for a FRET pair, the manufacture **Qd** should be "ok".
  * **Kappa2 = 2.0/3.0** : Dipole orientation factor. 

    In the literature, I often see this equal 2/3=0,667 This corresponds for free diffusing dye.
    For a dye attached to a protein, and restricted in one direction, the value is equal Kappa2=0.476
    But since we take the 6th root of the value, the difference ends up only being 5% in relative difference for the R0 value.
    Kappa2=2/3 can be a valid fast approximation, for the purpose of ordering dyes. But be careful and check the literature for discussions.
  * **n = 1.33** : 

    water=1.33, protein=1.4, n(2M GuHCl)=1.375
    Again, we take the sixth root, so the differences are "not so important".



## How the script works

  * The script integrate the Donor emission spectre, to divide with it. 

    This is done by simple numerical integration by the [Rectangle method](http://en.wikipedia.org/wiki/Rectangle_method)
  * All datapoints for the Acceptor excitation spectre are scaled, so that **e**(molar extinction coefficient) fits the wavelength where it is determined.
  * For the multiplication of the spectre, the script test that both data points exist before multiplication. 

    This can be troublesome, if there is 1-2 nm between each data point. Or if they are "un-even". Like Donor has 100.2 100.4 100.4 and Acceptor has 100.1 100.3 100.5
    But most spectre have much better resolution. Mostly 0.2 nm, and "even spaced"
    So this should work quite good as well.



The output is a text file, with all the datapoints. 

  1. Column: Donor: The input Emission wavelength.
  2. Column: Donor: The input Emission, arbitrary units of excitation.
  3. Column: Donor: The input Emission, arbitrary units of excitation, divided by the integrated area.
  4. Column: Acceptor: The input excitation wavelength.
  5. Column: Acceptor: The input excitation **e**(molar extinction coefficient).
  6. Column: Acceptor: The input excitation **e**(molar extinction coefficient), scaled correctly.
  7. Column: The calculated overlap for this datapoint.
  8. Column: The summed values for calculated overlap, until this datapoint.



  
It also automatically generates a [gnuplot](http://www.gnuplot.info/) .plt file, for super fast making graphs of the spectres. In linux, gnuplot is normally part of the installation. Just write: gnuplot FILENAME.plt If you are in windows, [download gnuplot](http://sourceforge.net/projects/gnuplot/files/) and open the .plt file. 

Gnuplot makes three graphs in .eps format. They can be converted to pdf with the program: epstopdf or eps2pdf. They are for example part of LaTeX: C:\Program Files (x86)\MiKTeX 2.9\miktex\bin or you can [download it here](http://tinyurl.com/eps2pdf). The format of .eps is choses, since gnuplot then allows for math symbols in the graphs. 

  * [![All graphs plotted together, but not scaled correctly.](/images/7/77/1-ALEXA488Emi-ALEXA633Exi-overlap-all-spectre.png)](/index.php/File:1-ALEXA488Emi-ALEXA633Exi-overlap-all-spectre.png "All graphs plotted together, but not scaled correctly.")

All graphs plotted together, but not scaled correctly. 

  * [![The Donor emission \(y1\) and Acceptor excitation \(y2\) scaled correctly.](/images/7/7d/2-ALEXA488Emi-ALEXA633Exi-overlap-normalized-spectre.png)](/index.php/File:2-ALEXA488Emi-ALEXA633Exi-overlap-normalized-spectre.png "The Donor emission \(y1\) and Acceptor excitation \(y2\) scaled correctly.")

The Donor emission (y1) and Acceptor excitation (y2) scaled correctly. 

  * [![The Rectangle integral datapoint \(y1\) and the summed integration \(y2\).](/images/0/09/3-ALEXA488Emi-ALEXA633Exi-overlap-integral.png)](/index.php/File:3-ALEXA488Emi-ALEXA633Exi-overlap-integral.png "The Rectangle integral datapoint \(y1\) and the summed integration \(y2\).")

The Rectangle integral datapoint (y1) and the summed integration (y2). 




## How to run the script

Make a pymol .pml file like this. Write in the required information, and execute/run the script with pymol. Then open the .plt file afterwards with gnuplot. 
    
    
    ## Change to your directory
    cd /homes/YOU/Documents/Atto-dyes/Spectre/ALEXA488-ALEXA633
    import forster_distance_calculator
    forster D_Exi=ALEXA488Exi.txt, D_Emi=ALEXA488Emi.txt, A_Exi=ALEXA633Exi.txt, A_Emi=ALEXA633Emi.txt, A_e_Max_Y=159000, A_e_Max_X=621, Qd=0.92
    

Retrieved from "[https://pymolwiki.org/index.php?title=Forster_distance_calculator&oldid=13908](https://pymolwiki.org/index.php?title=Forster_distance_calculator&oldid=13908)"


---

## Get colors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/get_colors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/get_colors.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  * **get_colors** contains two functions that can be useful when working with colors 
    * _get_colors_ : returns all available PyMOL colors
    * _get_random_color_ : returns a random available PyMOL color


  * Note that basic colors can be accessed manually without this script from the PyMOL menu under **Setting** \--> **Colors...**
  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



[![1LSD colored randomly by residue](/images/1/13/1LSD_random_colors.png)](/index.php/File:1LSD_random_colors.png "1LSD colored randomly by residue")

  


## Usage
    
    
    get_colors [ selection [, quiet ]]
    get_random_color [ selection [, quiet ]]
    

## Examples
    
    
    # basic example
    get_colors # basic colors
    get colors all # larger range with intermediates
    
    
    
    #get disco funky
    import get_colors
    from get_colors import get_colors
    from get_colors import get_random_color
    
    cmd.delete('all')
    cmd.fetch('1LSD', async=0) # :P
    cmd.hide('everything')
    cmd.show_as('sticks','not hetatm')
    cmd.orient()
    
    python # start a python block
    from pymol import stored
    stored.atom_list=[]
    cmd.iterate('name CA','stored.atom_list.append([model, resi])')
    resi_list=["model %s and resi %s"%(value[0],value[1]) for value in stored.atom_list]
    for resi in resi_list: cmd.color(get_random_color(),resi)
    python end # end python block
    cmd.ray()
    

|   
---|---  
  
  


## SEE ALSO

  * [Color](/index.php/Color "Color")
  * [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_colors&oldid=13942](https://pymolwiki.org/index.php?title=Get_colors&oldid=13942)"


---

## Get Model

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_model** returns a model object. 

## Contents

  * 1 PYMOL API
  * 2 USAGE
  * 3 NOTES
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.get_model(string "selection", integer "state" )
    

### USAGE
    
    
    cmd.get_model("chain A")
    

### NOTES

It can be useful to loop through all the atoms of a selection (rather than using the iterate command) 
    
    
    atoms = cmd.get_model("chain A")
    for at in atoms.atom:
        print("ATOM DEFINITION: "+at.model+" "\
                                 +at.chain+" "\
                                 +at.resn+" "\
                                 +str(at.resi)+" "\
                                 +at.name+" "\
                                 +str(at.index)+" "\
                                 +"%.2f " % (at.b)\
                                 +str(at.coord[0])+" "\
                                 +str(at.coord[1])+" "\
                                 +str(at.coord[2]))
    

### SEE ALSO

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Model&oldid=13576](https://pymolwiki.org/index.php?title=Get_Model&oldid=13576)"


---

## Get Title

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_title** retrieves a text string to the state of a particular object which will be displayed when the state is active. This is useful for printing the names of objects (given a state). 

## USAGE
    
    
    get_title object,state
    

## PYMOL API
    
    
     cmd.get_title(string object, int state)
    

## Example

Print out all the object names of the ensemble of states loaded in: 
    
    
    for x in range(numStates):
      print cmd.get_title("objName", x)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Title&oldid=13528](https://pymolwiki.org/index.php?title=Get_Title&oldid=13528)"


---

## Group

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The Group command creates or updates a "group" object. The grouped objects are collected underneath a **+** sign in the object tree (see images) in the Pymol Internal Gui. 

Group is tremendously helpful with multi-state or multi-structure sessions. Wildcards work great, for example: 
    
    
    # put all of objState into the group "ensemble".
    group ensemble, objState*
    

  * [![Three EF-Hand proteins loaded into PyMOL](/images/d/d2/Group_off.png)](/index.php/File:Group_off.png "Three EF-Hand proteins loaded into PyMOL")

Three EF-Hand proteins loaded into PyMOL 

  * [![Applied the group command to the proteins via: "group efHand, *"](/images/0/00/Group_on1.png)](/index.php/File:Group_on1.png "Applied the group command to the proteins via: "group efHand, *"")

Applied the group command to the proteins via: "group efHand, *" 

  * [![The plus was clicked and expanded to show the hierarchy of objects.](/images/7/7f/Group_on2.png)](/index.php/File:Group_on2.png "The plus was clicked and expanded to show the hierarchy of objects.")

The plus was clicked and expanded to show the hierarchy of objects. 




## Contents

  * 1 Usage
  * 2 Examples
    * 2.1 Creating, opening and closing
    * 2.2 More advanced usage of groups and naming
  * 3 Notes
  * 4 See Also



## Usage
    
    
    group name, members, action
    

Actions: 

  * **add** \- Add member to group
  * **remove** \- Remove members from group
  * **open** \- Open the group in the panel so objects can be dragged in
  * **close** \- Close the group in the panel so nothing can be dragged in
  * **toggle** \- Switch between open or close based on current state
  * **auto** \- (Deprecated)
  * **ungroup** \- (Deprecated) use the ungroup command instead
  * **empty** \- Move members to top level but do not delete groups
  * **purge** \- Delete members but do not delete groups
  * **excise** \- Delete groups but do not delete members
  * **raise** \- (Incentive 3.1+ only) Move the specified group to the top level, relevant for groups within groups



## Examples

### Creating, opening and closing
    
    
    group efHand, 1cll 1ggz 1sra
    
    # allow addition and removal from the group
    # If a group is open, objects can be added to or removed from 
    # it by right-click+drag from the control panel
    group efHand, open
    # disallow addition/removal from the group
    group efHand, close
    

### More advanced usage of groups and naming
    
    
    # names with dots are treated special
    
    set group_auto_mode, 2
    
    # load the example protein
    
    load $TUT/1hpv.pdb, 1hpv.other
    
    # create the new entry called ".protein" in group 1hpv
    
    extract 1hpv.protein, 1hpv.other and polymer
    
    # create ".ligand in the 1hpv group
    
    extract 1hpv.ligand, 1hpv.other and organic
    
    # supports wildcards
    
    show sticks, *.ligand
    
    hide lines, *.protein
    
    show surface, *.protein within 6 of *.ligand
    
    show lines, byres *.protein within 4 of *.ligand
    
    set two_sided_lighting
    
    set transparency, 0.5
    
    set surface_color, white
    
    # Also, to lexicographically sort the names in the control panel:
    
    order *, yes
    

## Notes

Group objects can usually be used as arguments to commands. It can be processed as a group or as a selection, in which case all the atoms from all objects in the group will be used. 

  


## See Also

[ungroup](/index.php?title=Ungroup&action=edit&redlink=1 "Ungroup \(page does not exist\)") [order](/index.php/Order "Order") [select](/index.php/Select "Select")

Retrieved from "[https://pymolwiki.org/index.php?title=Group&oldid=13873](https://pymolwiki.org/index.php?title=Group&oldid=13873)"


---

## Inertia tensor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/inertia_tensor.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/inertia_tensor.py)  
Author(s)  | [Mateusz Maciejewski](http://www.mattmaciejewski.com)  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![It.png](/images/d/dd/It.png)](/index.php/File:It.png)

This script will draw the eigenvectors of the inertia tensor of the selected part of the molecule. 

## Usage
    
    
    tensor selection, name, model
    

## Examples

To draw the inertia tensor of the first model in PDB file 4B2R do: 
    
    
    fetch 4b2r, async=0
    hide lines
    show cartoon
    run inertia_tensor.py
    tensor 4b2r & n. n+ca+c & i. 569-623, tensor_dom11, 1
    

## Reference

A figure generated using this script can be seen in the following reference: 

  * _Estimation of interdomain flexibility of N-terminus of factor H using residual dipolar couplings_. Maciejewski M, Tjandra N, Barlow PN. **Biochemistry**. 50(38) 2011, p. 8138-49, fig. 1 [doi:10.1021/bi200575b](http://dx.doi.org/10.1021/bi200575b)



Retrieved from "[https://pymolwiki.org/index.php?title=Inertia_tensor&oldid=13943](https://pymolwiki.org/index.php?title=Inertia_tensor&oldid=13943)"


---

## Label size

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Sets how large the labels are rendered. You can use positive numbers 2, 3, 4, etc for point sizes, or negative numbers for Angstrom-based sizes. Default is 14 points. Labels in Angstrom-size scale with the distance from the front plane, labels in point-size don't. This automatic scaling works with the [draw](/index.php/Draw "Draw") command but not with the [ray](/index.php/Ray "Ray") command. 

## Syntax
    
    
    # set the label size to 10pt
    set label_size, 10
    
    # set the label size to 1.5 Ang. -- large!
    set label_size, -1.5
    

Retrieved from "[https://pymolwiki.org/index.php?title=Label_size&oldid=13561](https://pymolwiki.org/index.php?title=Label_size&oldid=13561)"


---

## Main Page

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

hosted by [![SBGridlogo2.jpg](/images/e/ec/SBGridlogo2.jpg)](https://sbgrid.org/)  
---  
Welcome to the PyMOL Wiki!  The community-run support site for the [PyMOL](http://pymol.org) molecular viewer.   
---  
To request a new account, email SBGrid at: accounts (@) sbgrid dot org   
Quick Links  **[Tutorials](/index.php/Category:Tutorials "Category:Tutorials")** | **[Table of Contents](/index.php/TOPTOC "TOPTOC")** | **[Commands](/index.php/Category:Commands "Category:Commands")**  
---|---|---  
**[Script Library](/index.php/Category:Script_Library "Category:Script Library")** | **[Plugins](/index.php/Category:Plugins "Category:Plugins")** | **[FAQ](/index.php/Category:FAQ "Category:FAQ")**  
**[Gallery](/index.php/Gallery "Gallery")** | **[Covers](/index.php/Covers "Covers")** | **[PyMOL Cheat Sheet](/index.php/CheatSheet "CheatSheet")** (_[PDF](/images/7/77/PymolRef.pdf "PymolRef.pdf")_)  | **[Getting Help](/index.php/PyMOL_mailing_list "PyMOL mailing list")**  
News & Updates  | New Setup  | [PyMOL-open-source-windows-setup v3.1](https://github.com/kullik01/pymol-open-source-windows-setup/releases/tag/v3.1.0) has been released on January 20, 2025. More information under [Windows Install](/index.php/Windows_Install "Windows Install").   
---|---  
New Plugin  | [PySSA](/index.php/PySSA "PySSA") aims to combine PyMOL and [ColabFold](https://github.com/sokrypton/ColabFold) to enable the prediction and analysis of 3D protein structures for the scientific end-user. [v1.0 has been released](https://github.com/urban233/PySSA/releases/tag/v1.0.1) on July 10, 2024.   
Official Release  | [PyMOL v3.0 has been released](https://pymol.org) on March 12, 2024.   
New Plugin  | [CavitOmiX](/index.php/CavitOmiX "CavitOmiX") calculate [Catalophore™ cavities](https://innophore.com), predict protein structures with [OpenFold by NVIDIA-BioNeMo](https://www.nvidia.com/en-us/gpu-cloud/bionemo), [ESMFold](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/) and retrieve [Alphafold](https://www.deepmind.com/research/highlighted-research/alphafold) models   
Official Release  | [PyMOL v2.5 has been released](https://pymol.org) on May 10, 2021.   
Python 3  | New [Python 3 compatibility guide](/index.php/2to3 "2to3") for scripts and plugins   
POSF  | [New PyMOL fellowship announced for 2022-2023](https://pymol.org/fellowship)  
Tutorial  | [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial") updated for PyQt5   
New Plugin  | [PICv](/index.php/PICv "PICv") is a new plugin for clustering protein-protein interactions and visualization with available data from PDBe   
Selection keywords  | New [polymer.protein and polymer.nucleic](/index.php/Selection_Algebra "Selection Algebra") selection keywords. Thanks everyone who participated in the [poll](https://goo.gl/forms/r0Ck03VTytZQxN4A2)!   
Plugin Update  | [MOLE 2.5](/index.php/MOLE_2.0:_advanced_approach_for_analysis_of_biomacromolecular_channels "MOLE 2.0: advanced approach for analysis of biomacromolecular channels") is an updated version of channel analysis software in PyMOL   
New Script  | [dssr_block](/index.php/Dssr_block "Dssr block") is a wrapper for DSSR (3dna) and creates block-shaped nucleic acid cartoons   
Older News  | See [Older News](/index.php/Older_News "Older News").   
Did you know...  | 

### [Color h](/index.php/Color_h "Color h")

= Overview = This script colors the selection passed in based on the hydrophobicity scale as defined by: 

    

    Source: <http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html>
    Amino acid scale: Normalized consensus hydrophobicity scale
    Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
    Reference: J. Mol. Biol. 179:125-142 (1984)
Or check out hydrophobicity coloring, with rTools from Kristian Rother. <http://www.rubor.de/pymol_extensions_de.html>

# The Code

<source lang="python">

  1. color_h
  2. \-------


  1. PyMOL command to color protein molecules according to the Eisenberg hydrophobicity scale


  1.   2. Source: <http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html>
  3. Amino acid scale: Normalized consensus hydrophobicity scale
  4. Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
  5. Reference: J. Mol. Biol. 179:125-142 (1984)
  6.   7. Amino acid scale values:
  8.   9. Ala: 0.620
  10. Arg: -2.530
  11. Asn: -0.780
  12. Asp: -0.900
  13. Cys: 0.290
  14. Gln: -0.850
  15. Glu: -0.740
  16. Gly: 0.480
  17. His: -0.400
  18. Ile: 1.380
  19. Leu: 1.060
  20. Lys: -1.500
  21. [..→](/index.php/Color_h "Color h")

  
---  
|  [![](/images/1/13/070803_h.a.steinberg_molec_cell.jpg)](/index.php/File:070803_h.a.steinberg_molec_cell.jpg) [](/index.php/File:070803_h.a.steinberg_molec_cell.jpg "Enlarge")A Random PyMOL-generated Cover. See [Covers](/index.php/Covers "Covers").   
  
  
Retrieved from "[https://pymolwiki.org/index.php?title=Main_Page&oldid=13837](https://pymolwiki.org/index.php?title=Main_Page&oldid=13837)"


---

## Modevectors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/modevectors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/modevectors.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
### DESCRIPTION

[![](/images/3/34/Modevectors.png)](/index.php/File:Modevectors.png)

[](/index.php/File:Modevectors.png "Enlarge")

Example showing modevectors in action. (See Examples below).

**modevectors.py** is a PyMol script that was originally written to visualize results obtained from Normal Mode Analysis (NMA) by drawing arrows or vectors from a starting structure to a final structure. However, this script can be used to visualize the direction of motion between two specified states (e.g. targeted MD, comparing open and closed structures, etc). The strength of this script is that the arrows are highly customizable with numerous options to choose from (see script below). It is important to note that the two states MUST be identical except for the X, Y, and Z coordinates. That is, you cannot draw vectors between homologous structures. The default settings sets the background to white and displays the arrows along with the first object frame (in cartoon representation). 

  
Note: Ray tracing these CGO arrows appears to be a memory intensive process so please save a session before ray tracing and use "pymol -qc" where necessary. If anybody can come up with a method to improve this please e-mail me (see address in script). 

Update: A new way of drawing the arrows has been implemented for PyMOL v.1.1 and has been updated in this code (which automatically detects the version that you are using). The new method uses the built in cone primitive rather than concentric circles decreasing in size. Thus, the new method is much less memory intensive and is a significant improvement over the old method. Additionally, you may want to turn off [ray shadow](/index.php/Ray_shadow "Ray shadow") depending on the [scene](/index.php/Scene "Scene") that you are viewing as shadows can be confusing to others when seen in publications. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    modevectors first_obj_frame, last_obj_frame [,first_state=1 [,last_state=1 [,outname=modevectors \
      [,head=1.0 [,tail=0.3 \[,head_length=1.5 [,headrgb=(1.0,1.0,1.0) [,tailrgb=(1.0,1.0,1.0) [,cutoff=4.0 
      [,skip=0 [,cut=0.5 [,atom=CA [,stat=show [,factor=1.0 [,notail=0]]]]]]]]]]]]]]
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "modevectors" (which can be changed through the options). 

  


### EXAMPLES
    
    
    modevectors 1E3M, 1W7A
    modevectors 1E3M, 1W7A, outname="arrows"
    modevectors 1E3M, 1W7A, atom="P"
    

Copy/paste the following code to see an example of modevectors. This uses a multistate protein and the arrows are connected between the first and last states. 
    
    
    import modevectors
    # fetch the PDBs from pdb.org
    fetch 1c3y, finish=1, multiplex=0, async=0
    # separate the first and last states of the NMR ensemble to individual objects
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
    # run the modevectors code
    modevectors 1c3y_0001, 1c3y_0023
    # just setup a nice representation
    as cartoon, 1c3y_0001 or 1c3y_0023
    show cgo, modevectors
    color marine, 1c3y_0001
    color purpleblue, 1c3y_0023
    

The following set of examples will illustrate the power of each optional argument. Each example should be compared to the default figure in the table below. 

[![](/images/3/35/Mv_default.png)](/index.php/File:Mv_default.png) [](/index.php/File:Mv_default.png "Enlarge")Fig.1 - Default Settings | [![](/images/b/b0/Mv_fig2.png)](/index.php/File:Mv_fig2.png) [](/index.php/File:Mv_fig2.png "Enlarge")Fig.2 - Arrow direction drawn from CA to CA | [![](/images/4/4f/Mv_fig3.png)](/index.php/File:Mv_fig3.png) [](/index.php/File:Mv_fig3.png "Enlarge")Fig.3 - Arrow head radius  
---|---|---  
[![](/images/4/4f/Mv_fig4.png)](/index.php/File:Mv_fig4.png) [](/index.php/File:Mv_fig4.png "Enlarge")Fig.4 - Arrow tail radius | [![](/images/f/f1/Mv_fig5.png)](/index.php/File:Mv_fig5.png) [](/index.php/File:Mv_fig5.png "Enlarge")Fig.5 - Arrow head length | [![](/images/e/e4/Mv_fig6.png)](/index.php/File:Mv_fig6.png) [](/index.php/File:Mv_fig6.png "Enlarge")Fig.6 - Arrow head RGB  
[![](/images/f/fb/Mv_fig7.png)](/index.php/File:Mv_fig7.png) [](/index.php/File:Mv_fig7.png "Enlarge")Fig.7 - Arrow tail RGB | [![](/images/b/bc/Mv_fig8.png)](/index.php/File:Mv_fig8.png) [](/index.php/File:Mv_fig8.png "Enlarge")Fig.8 - Arrow length cutoff | [![](/images/7/7f/Mv_fig9.png)](/index.php/File:Mv_fig9.png) [](/index.php/File:Mv_fig9.png "Enlarge")Fig.9 - Skip every other arrow  
[![](/images/d/d7/Mv_fig10.png)](/index.php/File:Mv_fig10.png) [](/index.php/File:Mv_fig10.png "Enlarge")Fig.10 - Subtract value from vector length | [![](/images/4/42/Mv_fig11.png)](/index.php/File:Mv_fig11.png) [](/index.php/File:Mv_fig11.png "Enlarge")Fig.11 - Draw arrows from CB atoms | [![](/images/b/b1/Mv_fig12.png)](/index.php/File:Mv_fig12.png) [](/index.php/File:Mv_fig12.png "Enlarge")Fig.12 - Shorten by 50%  
| [![](/images/a/a3/Mv_fig13.png)](/index.php/File:Mv_fig13.png) [](/index.php/File:Mv_fig13.png "Enlarge")Fig.13 - Multiple options being used (see final example below)  
      
    
    reinitialize
    import modevectors
    
    fetch 1c3y, async=0
    
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
     
    #This is the default setting (Fig.1)
    modevectors 1c3y_0001, 1c3y_0023
     
    #This is shows that the direction of the arrows is drawn from the 1c3y_0001 towards 1c3y_0023 (Fig.2)
    show cartoon, 1c3y_0023
    color red, 1c3y_0023
     
    #This controls the base radius of the cone/arrow head in angstroms (Fig.3)
    modevectors 1c3y_0001, 1c3y_0023, head=2.5
     
    #This controls the radius of the cylinders/arrow tails in angstroms (Fig.4)
    modevectors 1c3y_0001, 1c3y_0023, tail=0.75
     
    #This controls the length of the cone/arrow head in angstroms (Fig.5)
    #Note that this option does NOT increase the vector length but simply changes the tail length
    modevectors 1c3y_0001, 1c3y_0023, head_length=3.0
     
    #This controls the colour of the cone/arrow head using RGB values (Fig.6)
    modevectors 1c3y_0001, 1c3y_0023, headrgb=(1.0,0.0,0.0)
     
    #This controls the colour of the cylinder/arrow tails using RGB values (Fig.7)
    modevectors 1c3y_0001, 1c3y_0023, tailrgb=(1.0,0.0,0.0)
     
    #This controls the which vectors to show based on a specific cutoff value in angstroms.  Vector lengths that are less 
    #than the cutoff value will not be displayed (Fig.8)
    modevectors 1c3y_0001, 1c3y_0023, cutoff=30.0
     
    #This controls how many vectors to skip (integer value) and is useful when there are too many vectors showing.  
    #Skip=1 will show every other arrow (Fig.9)
    modevectors 1c3y_0001, 1c3y_0023, skip=1
     
    #This controls how much to cut from each vector (in angstroms).  Note that arrows will point in the opposite directions 
    #when too much is cutoff (resulting in a negative vector length) (Fig.10) and should be used with caution!
    modevectors 1c3y_0001, 1c3y_0023, cut=15.0
     
    #This controls which atom to draw a vector from (Fig.11).  Note that this is case-sensitive and is really only useful 
    #when atom=CA or when atom=P (for DNA backbones)
    modevectors 1c3y_0001, 1c3y_0023, atom=CB
     
    #This controls how much to multiple the length of each vector by (percentage increase/decrease) (Fig.12)
    #This example halves the length of each vector (50%)
    modevectors 1c3y_0001, 1c3y_0023, factor=0.5
     
    #This hides the statistics which count atoms skipped, atoms counted (number of arrows showing), and number of atoms 
    #that did not meet the cutoff and are not shown
    modevectors 1c3y_0001, 1c3y_0023, stat=hide
     
    #This example shows multiple options being used (Fig.13)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, tail=1.0, head_length=2.0, headrgb=(1.0,0.0,0.0), tailrgb=(0.0,0.0,1.0),
    cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8
     
    #Finally, this example hides all arrow tails and only uses arrow heads via the notail option(No Figure)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8, notail=1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Modevectors&oldid=13914](https://pymolwiki.org/index.php?title=Modevectors&oldid=13914)"


---

## Movie color fade

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/movie_color_fade.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/movie_color_fade.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![movie_color_fade example](/images/0/05/Movie_color_fade_example_1.gif)](/index.php/File:Movie_color_fade_example_1.gif "movie_color_fade example")

**movie_color_fade** is like [movie_fade](/index.php/Movie_fade "Movie fade"), but will fade colors in a movie. Simply specify the arguments (see below). 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Examples
  * 3 SEE ALSO



## Usage
    
    
     movie_color_fade [ startframe [, startcolor [, endframe [, endcolor [, selection ]]]]]
     help movie_color_fade
    

|   
---|---  
  
### Arguments

  * startframe, endframe = beginning and end movie frame for fading 
    * if omitted, the current or last frame will be respectively set automatically
  * startcolor, endcolor = coloring at start and end
  * selection: target selection



## Examples

This example yields the movie to the top right 
    
    
    # make function available to PyMOL
    import movie_color_fade
    
    # object prep
    bg black
    delete all
    fetch 1hpv, async=0
    as cartoon
    orient
    color yellow, chain A
    color skyblue, chain B
    # movie prep
    mset 1x120
    
    # color fading
    movie_color_fade 1, yellow, 60, skyblue, chain A
    movie_color_fade 60, skyblue, 120, yellow, chain A
    movie_color_fade 1, skyblue, 60, yellow, chain B
    movie_color_fade 60, yellow, 120, skyblue, chain B
    

|   
---|---  
  
  


## SEE ALSO

  * [Movie_fade](/index.php/Movie_fade "Movie fade")
  * [mappend](/index.php/Mappend "Mappend")
  * [Color](/index.php/Color "Color")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_color_fade&oldid=13945](https://pymolwiki.org/index.php?title=Movie_color_fade&oldid=13945)"


---

## Movie fade

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/movie_fade.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/movie_fade.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will help fade in and out settings in a movie. Just specify the setting, it's initial value at an initial frame and it's ending value and frame. 

## Usage
    
    
    movie_fade setting, startFrame, startVal, endFrame, endVal [, selection ]
    

## Examples

To fade in sticks from fully transparent to fully opaque across 60 frames do: 
    
    
    mset 1x60
    movie_fade stick_transparency, 1, 1., 60, 0.
    

More complex example which involves camera motion: 
    
    
    fetch 1rx1, async=0
    as cartoon
    show surface
    mset 1x80
    movie.roll
    movie_fade transparency,  1, 0., 40, 1.
    movie_fade transparency, 41, 1., 80, 0.
    

**Example file:**

[![](/images/9/98/Movie_fade_example.gif)](/index.php/File:Movie_fade_example.gif)

[](/index.php/File:Movie_fade_example.gif "Enlarge")

Result
    
    
    import movie_fade
    
    fetch 1hpv, async=0
    orient
    
    #format
    bg black
    show_as cartoon
    color marine, ss s
    color red, ss h
    color white, ss l+""
    show sticks, resn 478
    util.cbao resn 478
    set surface_color, white, 1hpv
    show surface
    
    #movie
    mset 1x120
    movie_fade transparency,1, 0, 60, 1, 1hpv
    movie_fade transparency,61, 1, 120, 0, 1hpv
    

## See Also

  * [mdo](/index.php/Mdo "Mdo")
  * [mappend](/index.php/Mappend "Mappend")
  * [set](/index.php/Set "Set")
  * [movie_color_fade](/index.php/Movie_color_fade "Movie color fade")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_fade&oldid=13946](https://pymolwiki.org/index.php?title=Movie_fade&oldid=13946)"


---

## Outline

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [outline.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/outline.py)  
Author(s)  | [Jarrett Johnson](/index.php/User:JarrettJohnson "User:JarrettJohnson")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
  
## Contents

  * 1 Outline
    * 1.1 Relevant Settings
    * 1.2 Example
    * 1.3 Automation



# Outline

As of April 2023, the plugin referred to here will generate outlines given a selection and a list representations to outline with a specified color. 

(Pre-2023 using `ImageMagick`): 

Previously, user could employ `ImageMagick` to perform outlines by compositing ray-traced image. See information below. 

## Relevant Settings

  * [ray_opaque_background](/index.php/Ray_opaque_background "Ray opaque background")
  * [ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")
  * [ray_trace_color](/index.php?title=Ray_trace_color&action=edit&redlink=1 "Ray trace color \(page does not exist\)")
  * [ray_trace_gain](/index.php/Ray_trace_gain "Ray trace gain")
  * [surface_quality](/index.php/Surface_quality "Surface quality")



## Example

  * [![Structure of Cytochrome C \(PDB 3cyt\) shown in surface representation.](/images/2/20/Outline_base.png)](/index.php/File:Outline_base.png "Structure of Cytochrome C \(PDB 3cyt\) shown in surface representation.")

Structure of Cytochrome C (PDB 3cyt) shown in surface representation. 

  * [![Outline of residues within 2Å of Lys85 on chain I. Notice the image background is transparent.](/images/d/d9/Outline_overlay.png)](/index.php/File:Outline_overlay.png "Outline of residues within 2Å of Lys85 on chain I. Notice the image background is transparent.")

Outline of residues within 2Å of Lys85 on chain I. Notice the image background is transparent. 

  * [![Composite of the previous 2 images.](/images/a/a8/Outline_composite.png)](/index.php/File:Outline_composite.png "Composite of the previous 2 images.")

Composite of the previous 2 images. 

  * [![The same composite image, cleaned up.](/images/6/6b/Outline_cleaned.png)](/index.php/File:Outline_cleaned.png "The same composite image, cleaned up.")

The same composite image, cleaned up. 




  
The following script was used to generate the first two images from PyMOL. 
    
    
    bg_color white
    
    fetch 3cyt, async=0
    
    as surface
    color marine
    
    select outline, br. chain I and resi 85 around 2
    
    set_view (\
         0.061975956,   -0.950684488,    0.303902954,\
         0.703773856,    0.257531703,    0.662103057,\
        -0.707715809,    0.172843516,    0.685028315,\
         0.000000000,    0.000000000, -152.244812012,\
        25.274658203,    8.288025856,    9.110867500,\
        51.143974304,  253.345642090,  -20.000000000 )
    
    png base.png, ray=1
    
    hide everything
    as surface, outline
    set ray_trace_mode, 2
    set ray_trace_color, yellow
    set ray_opaque_background, 0
    
    png overlay.png, ray=1
    

In an image editing program, position overlay.png in a layer exactly overlapping base.png. The extra outline lines can also (optionally) be cleaned up by carefully erasing them using the graphics program's eraser tool. 

## Automation

If you need to have a fully automated solution, you can use ImageMagick or GraphicMagick to do this. If you use ImageMagick in a Unix box, you can type this command on the terminal: 
    
    
    composite -gravity center overlay.png base.png composite.png
    

See [this page](http://www.imagemagick.org/script/composite.php) on ImageMagick for more information. If you need to use it with GraphicMagick or in Windows or Mac, search for the relevant guidances. 

If you need to erase the lines inside the border, you can use the following command to do it. The `color` variable is the desired final output color of the border; change the `width` variable to change the border's thickness. 
    
    
    infile="overlay.png"
    outfile="overlay_cleaned.png"
    color="black"
    width=3
    w2=`convert $infile -format "%[fx:w-2]" info:`
    h2=`convert $infile -format "%[fx:h-2]" info:`
    
    convert $infile -background white -flatten -fill black +opaque white -bordercolor none -border 2 \
    -fill none -draw "matte 2,2 floodfill matte $w2,2 floodfill matte $w2,$h2 floodfill matte 2,$h2 floodfill" \
    -fill white +opaque none -fill black -opaque none \
    -alpha off -morphology edge octagon:$width \
    -channel rgba -fill none +opaque white -fill $color -opaque white -shave 2x2 $outfile
    

See [this thread](http://www.imagemagick.org/discourse-server/viewtopic.php?f=1&t=25970&p=113686#p113686) in ImageMagick forum for more detail. 

Retrieved from "[https://pymolwiki.org/index.php?title=Outline&oldid=13552](https://pymolwiki.org/index.php?title=Outline&oldid=13552)"


---

## Peptide Sequence

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Building a Peptide Sequence From Hand

There are more than one method in PyMOL for building a peptide sequence from hand. 

  * Simply hold down the **alt** (**option** on Mac) key and type in the one-letter code for the sequence.
  * Use the [fab](/index.php/Fab "Fab") command
  * Use the [Builder](/index.php/Builder "Builder") in "Protein" mode
  * You can write a script like the following which will build the amino acid sequence "DCAHWLGELVWCT".


    
    
    for aa in "DCAHWLGELVWCT": cmd._alt(string.lower(aa))
    

## Other Sources

Robert Campbell has notified us of [his two scripts](https://github.com/zigeuner/robert_campbell_pymol_scripts/tree/master/work_pymol) to solve the problem. You can apparently specify phi/psi angles, too. Look for **build_seq.py** and **build_seq_phi_psi.py**. 

Also check out [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure") , which you can build repeating units of different types of secondary structure (a-helix,b-sheets, loops,etc) (Dan Kulp) 

In addition, [Rotamer Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle") can set the sidechains to different rotamers or specific side-chain angle sets. (Dan Kulp) 

Retrieved from "[https://pymolwiki.org/index.php?title=Peptide_Sequence&oldid=13843](https://pymolwiki.org/index.php?title=Peptide_Sequence&oldid=13843)"


---

## Quickdisplays

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/quickdisplays.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/quickdisplays.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The module **quickdisplays** is a package containing several functions for quick standard displays:  


List of functions   
---  
Function  | What it does  | Usage   
disp_list  | prints a list of functions  | disp_list   
disp_ss  | secondary structure  | disp_ss [ selection [, colors [, only ]]]   
disp_ball_stick  | balls and sticks  | disp_ball_stick [ selection [, hydrogens [, only ]]]   
disp_mesh  | surface mesh  | disp_mesh [ selection [, color_m [, hydrogens [, only [, limits ]]]]]   
disp_surf  | surface  | disp_surf [ selection [, color_s [, transparency [, hydrogens [, solvent [, ramp_above [, only [, limits ]]]]]]]]   
disp_putty  | putty b-factor sausage  | disp_putty [ selection [, limits [, only ]]]   
  
  * Note that each function has an **individual help description** , callable by entering, e.g.:


    
    
    help disp_surf
    

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



  


Examples   
---  
[![Quickdisplays disp surf.png](/images/3/3d/Quickdisplays_disp_surf.png)](/index.php/File:Quickdisplays_disp_surf.png)**disp_surf**  
 _surface_ | [![Quickdisplays disp ss.png](/images/a/a8/Quickdisplays_disp_ss.png)](/index.php/File:Quickdisplays_disp_ss.png)**disp_ss**  
 _secondary structure_ | [![Quickdisplays disp putty.png](/images/5/51/Quickdisplays_disp_putty.png)](/index.php/File:Quickdisplays_disp_putty.png)**disp_putty**  
 _putty b-factor sausage_ |  The _combined displays_ example was produced like this: 
    
    
    import quickdisplays
    fetch 1hpv, async=0
    
    disp_putty all, limits=10
    disp_surf color_s=lightblue, transparency=0.8
    disp_ss chain B
    disp_ball_stick hetatm
    util.cbao hetatm
    set mesh_mode, 1 # make sure hetams are meshed
    set mesh_cutoff, 4.5
    disp_mesh resn 478, color_m=green
      
  
[![Quickdisplays disp mesh.png](/images/0/07/Quickdisplays_disp_mesh.png)](/index.php/File:Quickdisplays_disp_mesh.png)**disp_mesh**  
 _mesh_ | [![Quickdisplays disp ball stick.png](/images/4/4b/Quickdisplays_disp_ball_stick.png)](/index.php/File:Quickdisplays_disp_ball_stick.png)**disp_ball_stick**  
 _balls and sticks_ | [![Quickdisplays combo.png](/images/f/fc/Quickdisplays_combo.png)](/index.php/File:Quickdisplays_combo.png)_combined displays_  
see example   
  
## Some notes on the arguments

  * **selection** can be used to restrict the display to certain objects or selections, the default is 'all'
  * **only** can be _True_ or _False_ and will toggle whether the display will be added (_show_) or replace others (_show_as_)
  * **disp_mesh** and **disp_surf** support the color **putty** , which will color the object/selection by b-factor
  * **limits** defines the b-factor color range: 
    * a list entry will define absolute values **[min;max]**
    * a float value will calculate the corresponding percentiles **±value%** , _default=5_
  * setting **hydrogens** will add hydrogen to the model; if set to _=False_ , hydrogens are removed, if omitted the state of hydogens will be 'as is'
  * **colors** in **disp_ss** is flexible: 
    * set e.g. three colors, e.g. 'red green blue' for sheets, helices and loops, respectively (cf. example)
    * alternatively certain 'util.' functions can be used (cf. example)
    * setting _False_ once or for individual colors will suppress coloring



  


## Applied example

Enter the following lines one-by-one and observe how the display will be affected 
    
    
    # import function to PyMOL
    import quickdisplays
    
    # prep
    fetch 1hpv, async=0
    bg black
    orient
    disp_list # list of functions
    
    help disp_ss # check out help
    disp_ss
    disp_ss only=T
    disp_ss colors=smudge orange grey
    disp_ss colors=chartreuse False False # False suppresses coloring
    disp_ss colors=util.chainbow # You can use util.cbc, util.rainbow, ...
    
    help disp_ball_stick # check out help
    disp_ball_stick hetatm
    disp_ball_stick all, only=T
    util.cbaw
    disp_ball_stick hydrogens=1
    disp_ball_stick hydrogens=-1
    disp_stick_ball # will this work?
    
    help disp_mesh # check out help
    disp_mesh
    disp_mesh color_m=green, only=False
    disp_mesh color_m=green, only=T
    disp_ball_stick hetatm
    set mesh_skip, 2 # omits lines in mesh
    disp_mesh color_m=default, hydrogens=1, only=T
    disp_mesh color_m=putty, hydrogens=0, limits=20
    disp_mesh color_m=putty, limits=30
    disp_mesh color_m=putty, limits=[15,50]
    
    help disp_putty # check out help
    disp_putty chain A, only=False, limits=[15,50] #only default is True
    disp_putty all, limits=0
    disp_putty all, limits=10
    
    help disp_surf # check out help
    disp_surf color_s=white, solvent=1
    disp_surf color_s=white# solvent=0
    disp_surf color_s=lightblue, transparency=0.8
    

  


## SEE ALSO

[Git intro](/index.php/Git_intro "Git intro"), [Displaying Biochemical Properties](/index.php/Displaying_Biochemical_Properties "Displaying Biochemical Properties")

Retrieved from "[https://pymolwiki.org/index.php?title=Quickdisplays&oldid=13947](https://pymolwiki.org/index.php?title=Quickdisplays&oldid=13947)"


---

## Ray texture

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 See Also



## Overview

**ray_texture** sets what texture PyMol applies to the surface of your object during ray tracing. Options now consist of 

  1. None
  2. Matte 1
  3. Matte 2
  4. Swirl 1
  5. Swirl 2
  6. Fiber



## Usage
    
    
    set ray_texture, X      # where X = 0 ... 5
    

## Examples

The thumbnails are most likely too small to see the detail. Please click on each image to see a close up. 

  * [![Ray Texture 0 \(None\)](/images/6/68/Ray_texture_0.png)](/index.php/File:Ray_texture_0.png "Ray Texture 0 \(None\)")

Ray Texture 0 (None) 

  * [![Ray Texture 1 \(Matte 1\)](/images/d/d7/Ray_texture_1.png)](/index.php/File:Ray_texture_1.png "Ray Texture 1 \(Matte 1\)")

Ray Texture 1 (Matte 1) 

  * [![Ray Texture 2 \(Matte 2\)](/images/1/11/Ray_texture_2.png)](/index.php/File:Ray_texture_2.png "Ray Texture 2 \(Matte 2\)")

Ray Texture 2 (Matte 2) 

  * [![Ray Texture 3 \(Swirl1\)](/images/7/79/Ray_texture_3.png)](/index.php/File:Ray_texture_3.png "Ray Texture 3 \(Swirl1\)")

Ray Texture 3 (Swirl1) 

  * [![Ray Texture 4 \(Swirl 2\)](/images/7/79/Ray_texture_4.png)](/index.php/File:Ray_texture_4.png "Ray Texture 4 \(Swirl 2\)")

Ray Texture 4 (Swirl 2) 

  * [![Ray Texture 5 \(Fiber\)](/images/b/bb/Ray_texture_5.png)](/index.php/File:Ray_texture_5.png "Ray Texture 5 \(Fiber\)")

Ray Texture 5 (Fiber) 




# See Also

  * [ray_texture_settings](/index.php/Ray_texture_settings "Ray texture settings")



Retrieved from "[https://pymolwiki.org/index.php?title=Ray_texture&oldid=13562](https://pymolwiki.org/index.php?title=Ray_texture&oldid=13562)"


---

## RotationAxis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/draw_rotation_axis.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/draw_rotation_axis.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO cylinder representing a rotation axis for a given transformation. It's very useful for drawing the axes of rotational symmetry in an oligomeric assembly. 

The idea is to align two molecules/domains/chains/selections (using cmd.super) and extract the trasformation (TTT) matrix (T). The direction of the rotation axis, and a point are obtained from T and are used to create a cgo object representing the axis. The script generates two measures: one in the graphical screen (the rotation axis value, and the norm of the translation vector along the rotation axis) and some basic information in the command-line (the transformation matrix, the rotation angle, distance between centroids, and some pml lines that you can use to reproduce the axis...) 

As always with these type of things, you have to use it at your own risk. I did not try all possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado@gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run _path-to-the-script_ /draw_rotation_axis.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_axis('selection1', 'selection2')**

Please, pay attention to the apostrophes around the selections, you MUST use them. Also works with chains: 

**draw_axis('chain A', 'chain B')**

or objects: 

**draw_axis('obj01', 'obj02')**

Also, you can play a bit with the length, width and color of the axis. 

**draw_axis('selection1', 'selection2', scale_factor, width, r1, g1, b1, r2, g2, b2)**

**scale_factor** = to control the length of the axis, the default is 20 

**width** = to control the width of the axis. Default is 0.6 

**r1, g1, b1** = first RGB color code. Default is 1, 1, 1 

**r2, g2, b2** = second RGB color code to create the gradient. Default is 1, 0, 0. 

To create a single color axis, just made r1,g1,b1=r2,g2,b2 

# Examples
    
    
    # download the source and save as draw_rotation_axis.py
    run draw_rotation_axis.py
    fetch 2vak
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/d/df/2vak4.png)](/index.php/File:2vak4.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/3/37/2vak2.png)](/index.php/File:2vak2.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Another example, calculate the rotation axis of an homodimer
    run draw_rotation_axis.py
    fetch 3rfg
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/5/55/3rfg1.png)](/index.php/File:3rfg1.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/d/d6/3rfg4.png)](/index.php/File:3rfg4.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Clearly, each of the domains are made up with motifs with internal symmetry
    draw_axis('resi 200-236 and chain A', 'resi 238-274 and chain A', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # Also, you can create first the selections and use them to calculate the axis
    sele selection1, resi 200-236 and chain A
    sele selection2, resi 238-274 and chain A
    draw_axis('selection1', 'selection2', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # will produce the same result.
    

  * [![Internal rotation axis](/images/6/60/3rfg2.png)](/index.php/File:3rfg2.png "Internal rotation axis")

Internal rotation axis 

  * [![Some basic information is printed in the screen](/images/0/03/3rfg3.png)](/index.php/File:3rfg3.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 




Retrieved from "[https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923](https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923)"


---

## Spectrumany

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/spectrumany.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrumany.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.viewing](http://pymol.org/psicoredirect.php?psico.viewing)  
  
[![](/images/8/8b/SpectrumanyExample.png)](/index.php/File:SpectrumanyExample.png)

[](/index.php/File:SpectrumanyExample.png "Enlarge")

Coloring with a gradient of different green shades

This script works similar to the [spectrum](/index.php/Spectrum "Spectrum") command, but instead of predefined palettes, any color sequence can be used. 

The color sequence is given by a space separated list of colors, so palette "red_white_blue" is the same as color sequence "red white blue". 

# Example
    
    
    fetch 2x19
    
    # these two produce the same result
    spectrum count, red_white_blue, chain B
    spectrumany count, red white blue, chain B
    
    # gradient of different green colors
    spectrumany count, smudge palegreen limegreen limon green forest, chain B
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spectrumany&oldid=13930](https://pymolwiki.org/index.php?title=Spectrumany&oldid=13930)"


---

## Spectrumbar

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/spectrumbar.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrumbar.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
    * 1.1 USAGE
    * 1.2 EXAMPLES
  * 2 Script



## Description

Often times one may color their structure based upon some predefined color spectrum. However, this alone is not helpful when dealing with publication images. The purpose of this script is to allow the user to draw their own custom spectrum bar to accompany their structural images. 

Further work will be done to improve the functionality of this script. Please feel free to contact the author for function requests. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "spectrumbar" (which can be changed through the options) which is a cylinder colored according to the user-specified colors. 

Note: It is strongly recommended to turn off the specular reflections before ray tracing 

### EXAMPLES

[![](/images/b/ba/Sb_red.png)](/index.php/File:Sb_red.png) [](/index.php/File:Sb_red.png "Enlarge")Fig.1 - spectrumbar red | [![](/images/b/ba/Sb_red.png)](/index.php/File:Sb_red.png) [](/index.php/File:Sb_red.png "Enlarge")Fig.2 - spectrumbar 1,0,0 | [![](/images/f/f5/Sb_red_orange_yellow.png)](/index.php/File:Sb_red_orange_yellow.png) [](/index.php/File:Sb_red_orange_yellow.png "Enlarge")Fig.3 - spectrumbar red, orange, yellow  
---|---|---  
[![](/images/0/00/Sb_default.png)](/index.php/File:Sb_default.png) [](/index.php/File:Sb_default.png "Enlarge")Fig.4 - spectrumbar red,orange,yellow,green,blue,purple | [![](/images/0/00/Sb_default.png)](/index.php/File:Sb_default.png) [](/index.php/File:Sb_default.png "Enlarge")Fig.5 - spectrumbar 1,0,0, orange, yellow, 0,1,0, 0,0,1, purple | [![](/images/5/58/Sb_radius.png)](/index.php/File:Sb_radius.png) [](/index.php/File:Sb_radius.png "Enlarge")Fig.6 - spectrum red,radius=2  
[![](/images/6/63/Sb_long_red_orange_long_yellow.png)](/index.php/File:Sb_long_red_orange_long_yellow.png) [](/index.php/File:Sb_long_red_orange_long_yellow.png "Enlarge")Fig 7 - spectrumbar red, red, red, orange, yellow, yellow | [![](/images/b/ba/Sb_rounded.png)](/index.php/File:Sb_rounded.png) [](/index.php/File:Sb_rounded.png "Enlarge")Fig 8 - spectrumbar red, radius=2, ends=rounded  
  
## Script

load the script using the [run](/index.php/Run "Run") command 
    
    
    from pymol.cgo import *
    from math import *
    from pymol import cmd
    from re import *
    
    def spectrumbar (*args, **kwargs):
    
        """
        Author Sean M. Law
        University of Michigan
        seanlaw_(at)_umich_dot_edu
    
        USAGE
    
        While in PyMOL
    
        run spectrumbar.py
    
        spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)
    
        Parameter     Preset         Type     Description
        RGB_Colors    [1.0,1.0,1.0]  N/A      RGB colors can be specified as a
                                              triplet RGB value or as PyMOL
                                              internal color name (i.e. red)
        radius        1.0            float    Radius of cylindrical spectrum bar
        name          spectrumbar    string   CGO object name for spectrum bar
        head          (0.0,0.0,0.0)  float    Starting coordinate for spectrum bar
        tail          (10.0,0.0,0.0) float    Ending coordinate for spectrum bar
        length        10.0           float    Length of spectrum bar
        ends          square         string   For rounded ends use ends=rounded
    
        Examples:
    
        spectrumbar red, green, blue
        spectrumbar 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0
    
        The above two examples produce the same spectrumbar!
    
        spectrumbar radius=5.0
        spectrumbar length=20.0
    
        """
    
        rgb=[1.0, 1.0, 1.0]
        name="spectrumbar"
        radius=1.0
        ends="square"
        x1=0
        y1=0
        z1=0
        x2=10
        y2=0
        z2=0
        num=re.compile('[0-9]')
        abc=re.compile('[a-z]')
    
        for key in kwargs:
            if (key == "radius"):
                radius = float(kwargs["radius"])
            elif (key == "name"):
                name=kwargs["name"]
            elif (key == "head"):
                head=kwargs["head"]
                head=head.strip('" []()')
                x1,y1,z1=map(float,head.split(','))
            elif (key == "tail"):
                tail=kwargs["tail"]
                tail=tail.strip('" []()')
                x2,y2,z2=map(float,tail.split(','))
            elif (key == "length"):
                if (abc.search(kwargs["length"])):
                    print "Error: The length must be a value"
                    return
                else:
                    x2=float(kwargs["length"]);
            elif (key == "ends"):
                ends=kwargs["ends"]
            elif (key != "_self"):
                print "Ignoring unknown option \""+key+"\""
            else:
                continue
    
        args=list(args)
        if (len(args)>=1):
            rgb=[]
        while (len(args)>=1):
            if (num.search(args[0]) and abc.search(args[0])):
                if (str(cmd.get_color_tuple(args[0])) != "None"):
                    rgb.extend(cmd.get_color_tuple(args.pop(0)))
                else:
                    return
            elif (num.search(args[0])):
                rgb.extend([float(args.pop(0))])
            elif (abc.search(args[0])):
                if (str(cmd.get_color_tuple(args[0])) != "None"):
                    rgb.extend(cmd.get_color_tuple(args.pop(0)))
                else:
                    return
            else:
                print "Error: Unrecognized color format \""+args[0]+"\""
                return
        
        if (len(rgb) % 3):
            print "Error: Missing RGB value"
            print "Please double check RGB values"
            return
    
        dx=x2-x1
        dy=y2-y1
        dz=z2-z1
        if (len(rgb) == 3):
            rgb.extend([rgb[0]])
            rgb.extend([rgb[1]])
            rgb.extend([rgb[2]])
        t=1.0/(len(rgb)/3.0-1)
        c=len(rgb)/3-1
        s=0
        bar=[]
        
        while (s < c):
            if (len(rgb) >0):
                r=rgb.pop(0)
                g=rgb.pop(0)
                b=rgb.pop(0)
            if (s == 0 and ends == "rounded"):
                bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz, radius])
            bar.extend([CYLINDER])
            bar.extend([x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz])
            bar.extend([x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz])
            bar.extend([radius, float(r), float(g), float(b)])
            if (len(rgb) >= 3):
                bar.extend([float(rgb[0]), float(rgb[1]), float(rgb[2])])
                r=rgb[0]
                g=rgb[1]
                b=rgb[2]
            else:
                bar.extend([float(r), float(g), float(b)])
            if (s == c-1 and ends == "rounded"):
                bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz, radius])
            s=s+1
    
        cmd.delete(name)
        cmd.load_cgo(bar,name)
        
    
        return
    cmd.extend("spectrumbar",spectrumbar)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spectrumbar&oldid=13949](https://pymolwiki.org/index.php?title=Spectrumbar&oldid=13949)"


---

## Sphere quality

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



## Overview

**sphere_quality** controls the rendering quality of sphere objects. This setting only affects sphere rendering when not using shaders. 

## Syntax
    
    
    # the default value is 1
    set sphere_quality, <integer>
    

Larger values of <integer> result in higher quality sphere rendering. Values >1 may result in poor [performance](/index.php/Category:Performance "Category:Performance") during real-time rotation or translation. 

**Note** : Selecting values larger than 2 with **stick_ball** = 1 (enabled) causes PyMol to crash in the Windows version. 

## Examples

Open the images to see rendering details. 

  * [![sphere_quality 0](/images/0/05/Sphere_quality_0.png)](/index.php/File:Sphere_quality_0.png "sphere_quality 0")

sphere_quality 0 

  * [![sphere_quality 1 \(default\)](/images/6/6e/Sphere_quality_1.png)](/index.php/File:Sphere_quality_1.png "sphere_quality 1 \(default\)")

sphere_quality 1 (default) 

  * [![sphere_quality 2](/images/2/2f/Sphere_quality_2.png)](/index.php/File:Sphere_quality_2.png "sphere_quality 2")

sphere_quality 2 




# See Also

[cgo_sphere_quality](/index.php?title=Cgo_sphere_quality&action=edit&redlink=1 "Cgo sphere quality \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_quality&oldid=13559](https://pymolwiki.org/index.php?title=Sphere_quality&oldid=13559)"


---

## Spheres

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/8/8a/Spheres_ex.png)](/index.php/File:Spheres_ex.png)

[](/index.php/File:Spheres_ex.png "Enlarge")

Normal Sphere Representation Example

## Representation

To enable the **spheres** representation do the following for any selection SEL, 
    
    
    show spheres, SEL
    

## Adjusting Sphere Size
    
    
    alter selection, vdw=number
    

### Examples

Shrink the size of all Iron atoms: 
    
    
    alter elem fe, vdw=1.0
    rebuild
    

Dramatically enlarge all spheres in an object 
    
    
    alter object, vdw=4.0
    rebuild
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spheres&oldid=13560](https://pymolwiki.org/index.php?title=Spheres&oldid=13560)"


---

## Sync

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

sync is an API-only function which waits until all current commands have been executed before returning. A timeout can be used to insure that this command eventually returns. 

## PYMOL API
    
    
    cmd.sync(timeout: float = 1.0, poll: float = 0.05)
    

## SEE ALSO

[frame](/index.php/Frame "Frame")

Retrieved from "[https://pymolwiki.org/index.php?title=Sync&oldid=13551](https://pymolwiki.org/index.php?title=Sync&oldid=13551)"


---

## Timeline Custom Programs

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Introduction

[Programs](/index.php?title=Programs&action=edit&redlink=1 "Programs \(page does not exist\)") are animations that occur over a function of time. Before PyMOL3, there were several included programs users could choose from: nutate, rock, roll, etc... For the Timeline, the user can create their own animations as a function of relative time (t). 

All information on this page is below subject to change during PyMOL3 beta (and beyond). 

## Class structure

In order to create a custom program, the only requirements currently is to create a class that contains a method `animate` that contains a single float parameter (t). This parameter represents the relative time within the range of [0, 1] that describes any change of animation or property over time. Programs are represented as clips (rectangle shapes) on a track and must not overlap with any other clip or keyframe on the same subtrack. 
    
    
    class MyCustomProgram:
        def animate(t: float) -> None:
            # ... logic here based on t
    

and to register this program to PyMOL: 
    
    
    timeline_register_program(program_name: str, program_type: object)
    

See [Timeline_Python_API](/index.php/Timeline_Python_API "Timeline Python API") for more information about this function 

  


## Example
    
    
    class RailLook:
        """
        Example PyMOL Timeline Program
        RailLook will move an object along a curve object while orienting
        itself toward a specified direction given by a target object.
        """
        def __init__(self, track) -> None:
            self.this_object = "_Camera"
            self.curve_object = "Curve01"
            self.target_object = "1obyA"
    
        def animate(self, t: float) -> None:
            cmd.move_on_curve(self.this_object, self.curve_object, t)
            cmd.look_at(self.target_object, self.this_object)
    
    cmd.timeline_register_program(program_name="rail_look", program_type=RailLook)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Timeline_Custom_Programs&oldid=13568](https://pymolwiki.org/index.php?title=Timeline_Custom_Programs&oldid=13568)"


---

## Timeline Python API

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Concepts & Classes
    * 1.1 Keyframe
    * 1.2 Clip
    * 1.3 Program
    * 1.4 Track
    * 1.5 Composition
    * 1.6 Sequence Composition
  * 2 Timeline Python API
    * 2.1 General comments of the API
    * 2.2 Explanation of terms
    * 2.3 API



## Concepts & Classes

**Note:** This is a incentive-only 3.0+ feature 

### Keyframe

Represents a state of an entity at a given time. Here, state does not necessarily pertain to the PyMOL concept of state but any property (transform, PyMOL state, setting, etc…) for a given entity (PyMOL Object or Camera). 

### Clip

Represents a collection of keyframes that describes an animation over time. 

### Program

Represents an animation of an object as a function of t. Where t is the interpolation factor from [0, 1] (0 represents the beginning of the animation and 1 represents the end). Visually, they are represented as clips with a keyframe count of two. PyMOL supplies several builtin programs including: nutate, rock, roll, state loop, etc… See Creating Custom Programs below to add complex, custom animations in the timeline. 

### Track

Represents a set of related keyframes and clips. These keyframes and clips are interpolated to and through each other within the same track and not outside. For each animated property, there shall be one Timeline track. 

### Composition

Represents a set of tracks. Allows for multiple objects and/or properties to be animated simultaneously. 

### Sequence Composition

A special composition that only stores clips of “normal” compositions. This is useful to combine multiple compositions into a single animation and blend between them seamlessly. Currently, the blending effect is not yet implemented. 

## Timeline Python API

### General comments of the API

All timeline functions force named-arguments for all parameters of each function. Thus, you must specify the name of the parameter followed by the argument value you want to pass. Examples are shown below. Also, uncommon from the rest of most of the other PyMOL API, builtin types are not commonly passed nor returned. Parameter and return types are specified for each function and a full description of their meaning is explained above. Examples of movies using the Python API can be found at: www.pymol.org. 

### Explanation of terms

  * _**Arguments**_ : a sequence of named parameters. For the timeline API, all arguments are named; positional arguments are disallowed. No arguments needed are shown by an empty sequence []. Defaulted arguments are shown by a proceeding assignment with value (ex: myParam = “myValue”). All non-defaulted parameters in this sequence must be provided. Sequence representation is for the purpose of documentation only. A sequence is not expected as input for the API.
  * _**Return**_ : a sequence of returned value(s). No returned values are depicted by empty sequence []. Sequence representation is for the purpose of documentation only. A sequence is not expected as input for the API.
  * _**Description**_ : explanation of the use of the API function call.
  * _**Example**_ : brief example of how to use the API.
  * _**Query**_ : As these functions mutate the timeline state, developers may often need to retrieve the current state of a timeline composition or element. These properties can be fetched from a variable that is synchronized with the Timeline backend and doesn’t require a call via `cmd`.
  * _**See Also**_ : Related timeline API calls. Many of these examples could be paired with the presented function.



### API

timeline_add_composition

_Signature_
    
    
    def timeline_add_composition() -> Composition
    

_Description_

    Appends a new composition to the timeline.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    

_See Also_

    

timeline_delete_composition

* * *

timeline_delete_composition

_Signature_
    
    
    def timeline_delete_composition(comp: Composition) -> None
    

_Description_

    Deletes the composition from the timeline. The composition and all track, clip, and keyframe instances managed under the deleted composition are thus invalidated and should not be used in the API again.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    cmd.timeline_delete_composition(comp=myComp)
    # myComp now invalidated
    

_See Also_

    

_timeline_add_composition_

* * *

timeline_add_sequence_composition

_Signature_
    
    
    def timeline_add_sequence_composition() -> Composition
    

_Description_

    Appends a new sequence composition to the Timeline.

_Example_
    
    
    mySeqComp = cmd.timeline_add_sequence_composition()
    

_See Also_

    

_timeline_add_composition_

* * *

timeline_set_composition_duration

_Signature_
    
    
    def timeline_set_composition_duration(comp: Composition, duration: int) -> None:
    

_Description_

    Sets the duration of comp to duration number of frames.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    cmd.timeline_set_composition_duration(comp=myComp, duration=120)
    # myComp is now 2 seconds (assuming myComp is set to at 60 fps (default))
    

_Query_

    Duration can be queried by duration.

_Example_
    
    
    cmd.timeline_set_composition_duration(comp=myComp, duration=240)
    print(myComp.duration) # prints “240”
    

_See Also_

    

timeline_add_composition

* * *

timeline_scrub_to

_Signature_
    
    
    def timeline_scrub_to(comp: Composition, frame_num: int)
    

_Description_

    Sets the current frame number of the composition to the given frame number. The argument frame_num may be outside the bounds of the composition’s duration but will not be reached during the animation.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    cmd.timeline_scrub_to(comp=myComp, frame_num=50)
    

_See Also_

    

timeline_set_composition_duration

* * *

timeline_add_object_track

_Signature_
    
    
    def timeline_add_object_track(comp: Composition, obj_name: str, track_type: TrackType=TrackType.TRANSFORM, track_type_detail: str="") -> Track
    

_Description_

    Appends a new object track onto comp. obj_name must be a valid name of an object currently managed in PyMOL.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    foo_track = cmd.timeline_add_object_track(comp=myComp, obj_name="foo")
    

_See Also_

    

timeline_delete_track

* * *

timeline_delete_track

_Signature_
    
    
    def timeline_delete_track(track: Track) -> None:
    

  
_Description_

    Deletes track from its composition. The track and all clip and keyframe instances managed under the deleted track are thus invalidated and should not be used in the API again.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    foo_track = cmd.timeline_add_object_track(comp=myComp, obj_name="foo")
    

_See Also_

    

timeline_add_object_track

* * *

timeline_store_scene_at

_Signature_
    
    
    def timeline_store_scene_at(comp: Composition, frame_num: int, scene_name: str) -> Keyframe
    

_Description_

    Appends a new scene keyframe onto comp at frame frame_num with the scene of name scene_name.

_Example_
    
    
    scene_kf = cmd.timeline_store_scene_at(comp=myComp, frame_num=42, scene_name="foo")
    cmd.timeline_delete_keyframe(keyframe=scene_kf)

_See Also_

    

timeline_add_keyframe

* * *

timeline_add_keyframe

_Signature_
    
    
    def timeline_add_keyframe(track: Track, frame_num: int, obj_state: int = -1) -> Keyframe
    

_Description_

    Appends a new object keyframe onto track at frame frame_num.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    

_See Also_

    

timeline_add_keyframe

* * *

timeline_delete_keyframe

_Signature_
    
    
    def timeline_delete_keyframe(keyframe: Keyframe) -> None
    

_Description_

    Deletes a keyframe

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_delete_keyframe(keyframe=bar_kf)
    

* * *

timeline_set_keyframe_interpolation

_Signature_
    
    
    def timeline_set_keyframe_interpolation(keyframe: Keyframe, interp_mode, KeyframeInterpolationMode) -> None
    

_Description_

    Describes how the properties between keyframe and the following keyframe are interpolated. By default, virtually all keyframes interpolate linearly (each time step between two keyframe changes in property by the same amount). Values of _KeyframeInterpolationMode_ include: 

    _CONSTANT_ ,
    _LINEAR_ ,
    _EASEIN_ ,
    _EASEOUT_ ,
    _EASEINOUT_

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_set_keyframe_interpolation(keyframe=bar_kf, interp_mode=KeyframeInterpolationMode.LINEAR)
    

* * *

timeline_move_keyframe_to

_Signature_
    
    
    def timeline_move_keyframe_to(keyframe: Track, frame_num: int) -> None
    

_Description_

    Moves keyframe to frame frame_num.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_move_keyframe_to(keyframe=bar_kf, frame_num=24)
    

_See Also_

    

timeline_move_clip_to

* * *

timeline_move_clip_to

_Signature_
    
    
    def timeline_move_clip_to(clip: Clip, frame_num: int)
    

_Description_

    Moves clip to frame frame_num.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_delete_keyframe(keyframe=bar_kf)
    

_See Also_

    

timeline_move_keyframe_to

* * *

timeline_set_clip_duration

_Signature_
    
    
    def timeline_set_clip_duration(clip: Clip, duration: int) -> None
    

_Description_

    Sets the duration of clip to duration frames.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_delete_keyframe(keyframe=bar_kf)
    

_See Also_

    

timeline_set_clip_speed, timeline_move_clip_to

* * *

timeline_set_clip_speed

_Signature_
    
    
    def timeline_set_clip_speed(clip: Clip, duration: int) -> None
    

_Description_

    Sets the duration of clip to duration frames.

_Example_
    
    
    bar_kf = cmd.timeline_set_clip_speed(track=foo_track, frame_num=42)
    

_See Also_

    

timeline_set_clip_duration, timeline_move_clip_to

* * *

timeline_register_program

_Signature_
    
    
    def timeline_register_program(program_name: str, program_type: object) -> None
    

_Description_

    Registers the program with program_type type and program_name name to an internal program database. The types are registered so that instances of it can be generated and represented as animation clips on the timeline. After a program is registered, you can store a program using timeline_store_program_at using program_name (case-sensitive) to identify the program.

_Example_
    
    
    class Drift:
        def animate(t: float) -> None:
            # ...
    
    cmd.timeline_register_program(program_name="drift", program_type=Drift)
    

_See Also_

    

timeline_store_program_at

* * *

timeline_store_program_at

_Signature_
    
    
    def timeline_store_program_at(track: Track, frame_num: int, duration: int, program: Union[object, str]) -> Clip
    

_Description_

    Adds an animation clip described by program onto track at frame frame_num for duration frames.

_Example_
    
    
    cmd.timeline_register_program(program_name="drift", program_type=Drift)
    cam_track = myComp.get_main_camera_track()
    cmd.timeline_store_program_at(track=cam_track, frame_num=60, duration=360, program="drift")
    

_See Also_

    

timeline_register_program

_Notes_

    To customize the properties of a program, instantiate a program, edit its properties, and provide it to the function. Example:
    
    
    roll_prg= pymol.timeline_programs.Roll(cam_track)
    roll_prg.num_loops = 2
    clip = cmd.timeline_store_program_at(track=cam_track, frame_num=60, duration=360, program=roll_prg)
    

    Since the program itself is agnostic to the duration of the clip that it’s stored in, in order to change the duration of the program, change the clip’s duration instead via timeline_set_clip_duration.

* * *

timeline_produce (NOT YET SUPPORTED)

_Signature_
    
    
    def timeline_produce(comp: Composition, render_mode: str = "draw", first_frame: int = None, last_frame: int = None, encoder: str = None, width: int = None, height: int = None, quality: int=100, quiet: bool=False)
    

[comp: Composition, render_mode: str = "draw", first_frame: int = None, last_frame: int = None, encoder: str = None, width: int = None, height: int = None, quality: int=100, quiet=False] 

_Description_

    Exports a comp to a playable movie (file format determined by encoder).

_Example_
    
    
    cmd.movie.timeline_produce(comp=comp)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Timeline_Python_API&oldid=13874](https://pymolwiki.org/index.php?title=Timeline_Python_API&oldid=13874)"


---

## ToGroup

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/togroup.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/togroup.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Examples
  * 3 The Code
  * 4 See Also



# Overview

toGroup will convert a multistate object into a group of single-state objects. _Be warned, by default it deletes your original object (since it's extracting a copy)._

PyMOL does a great job at handling multistate objects and grouping them together. One thing that I found myself doing over and over again was 

  * loading a multistate object (say a PDBQT file with 100 ligand poses)
  * splitting that object into all 100 states, with some given prefix
  * then grouping them into their own group
  * and then finally removing the original.



This became tedious, so I automated that with this script. 

# Examples
    
    
    # A multistate object (20 NMR states)
    fetch 1nmr
    
    # Create the group called,  "nmrEnsemble"
    # from '1nmr' and name all the new states state1,
    # state2, state3, etc.
    toGroup nmrEnsemble, 1nmr, prefix=state
    

# The Code
    
    
    import pymol
    from pymol import cmd
    
    def toGroup(groupName,sel,prefix="",delOrig=True):
        """
        DESCRIPTION
        toGroup will take a multistate object and extract it
        to a group with N objects all in state #1.  It essentially
        performs the following:
    
        split_states myObj, prefix=somePrefix
        group newGroup, somePrefix*
        delete myObj
    
        PARAMETERS:
        
        groupName (string)
            The name of the group to create
    
        sel (string)
            The name of the selection/object from which
            to make the group
    
        prefix (string)
            The prefix of the names of each of split states.
            For example, if your prefix is ''obj'' and is in
            states 1 through 100 then the states will be labeled
            obj1, obj2, obj3, ..., obj100.
    
        delOrig (string/boolean)
            If true then delete the original selection, otherwise not.
    
        RETURN
    
        Nothing, it makes a new group.
    
        """
        if prefix=="":
            prefix="grouped"
    
        cmd.split_states(sel, prefix=prefix)
        cmd.group(groupName,prefix+"*")
        
        if delOrig:
            cmd.delete(sel)
    
    cmd.extend("toGroup", toGroup)
    

# See Also

[group](/index.php/Group "Group"), [saveGroup](/index.php/SaveGroup "SaveGroup"), [select](/index.php/Select "Select"),, [split_states](/index.php/Split_states "Split states"), [delete](/index.php/Delete "Delete"), [extend](/index.php/Extend "Extend"). 

Retrieved from "[https://pymolwiki.org/index.php?title=ToGroup&oldid=13950](https://pymolwiki.org/index.php?title=ToGroup&oldid=13950)"


---

## Valence

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Examples for the settings: valence and valence_mode
    * 2.1 Syntax
  * 3 The valence command
  * 4 Editing bonds
  * 5 Automatic editing of bonds
  * 6 SEE ALSO



## Overview

Turning on the **valence** setting will enable the display of double bonds.  
Toggling **valence_mode** alters the positioning of double bonds (for representation as [Lines](/index.php/Lines "Lines"))  
**valence_size** alters the distance of double bonds.  
Note that bonds can be edited to be delocalized using [Unbond](/index.php/Unbond "Unbond") and [Bond](/index.php/Bond "Bond").  
There is also a command called **valence**.  


  


## Examples for the settings: valence and valence_mode

  * [![set valence, 0 #\(no double bonds\)](/images/1/10/PHE_valence_0.png)](/index.php/File:PHE_valence_0.png "set valence, 0 #\(no double bonds\)")

set valence, 0  
#(no double bonds) 

  * [![set valence, 1 set valence_mode, 1 #bonds inside](/images/9/9f/PHE_valence_1_mode_1.png)](/index.php/File:PHE_valence_1_mode_1.png "set valence, 1 set valence_mode, 1 #bonds inside")

set valence, 1  
set valence_mode, 1  
#bonds inside 

  * [![set valence, 1 set valence_mode, 0 #bonds centered](/images/8/8e/PHE_valence_1_mode_0.png)](/index.php/File:PHE_valence_1_mode_0.png "set valence, 1 set valence_mode, 0 #bonds centered")

set valence, 1  
set valence_mode, 0  
#bonds centered 

  * [![set valence, 1 #delocalized bonds #see section "Editing bonds" below or try command "valence guess, all"](/images/c/ce/PHE_delocalized.png)](/index.php/File:PHE_delocalized.png "set valence, 1 #delocalized bonds #see section "Editing bonds" below or try command "valence guess, all"")

set valence, 1  
#delocalized bonds  
#see section "Editing bonds" below or try command "valence guess, all" 




  
**valence_size** alters the distance of double bonds, but behaves slightly different depending on valence_mode  


valence_size  valence_size with valence_mode 1  
inside | valence_size with valence_mode 0   
centered   
---|---  
[![Valence size mode1.gif](/images/b/b1/Valence_size_mode1.gif)](/index.php/File:Valence_size_mode1.gif) | [![Valence size mode0.gif](/images/d/d9/Valence_size_mode0.gif)](/index.php/File:Valence_size_mode0.gif)  
  
### Syntax
    
    
    set valence, 0 # off
    set valence, 1 # on
    
    set valence_mode, 0 # centered
    set valence_mode, 1 # inside
    
    set valence_size, 0.1 # default: 0.06 # range 0 - ~0.5
    

## The valence command

The **valence** command automatically formats existing bonds and can even guess the bonds for standard amino acids.  

    
    
    #USAGE:
    valence order, selection1 [, selection2 [, source [, target_state [, source_state [, reset [, quiet ]]]]]]
    order can be either: 1, 2, 3, 4, aromatic, copy, guess
    
    #make PyMOL guess/autoformat bonds in proteins
    valence guess, all
    

  


## Editing bonds
    
    
    # In editing mode: select the bond using Ctrl-right-click, then enter:
    unbond pk1,pk2
    bond pk1,pk2,4
    # 1: single bond, 2: double bond, 3:triple bond, 4:delocalized
    

  


## Automatic editing of bonds

Try using the **valence** command first.  
Secondly, [Format_bonds](/index.php/Format_bonds "Format bonds") is a script that automatically formats valence in all amino acids and has additional options. 

  


## SEE ALSO

[Bond](/index.php/Bond "Bond"), [Unbond](/index.php/Unbond "Unbond")

Retrieved from "[https://pymolwiki.org/index.php?title=Valence&oldid=13567](https://pymolwiki.org/index.php?title=Valence&oldid=13567)"


---

## Windows Install

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes how to install PyMOL on Microsoft Windows. 

## Contents

  * 1 Incentive PyMOL
  * 2 Open-Source PyMOL
    * 2.1 One-click Installer
    * 2.2 Pre-compiled Open-Source PyMOL
    * 2.3 Building from source
    * 2.4 Extend PyMOL with additional scripts
  * 3 See Also



## Incentive PyMOL

[Schrödinger](http://www.schrodinger.com) provides an installer to paying sponsors (EXE for PyMOL 2.0, MSI for previous version). The bundle also includes ready-to-use [APBS](/index.php/APBS "APBS"), [RigiMOL](/index.php/Morph "Morph"), an MPEG encoder for movie export, and a small molecule energy minimization engine. 

Download: <https://pymol.org/#download>

## Open-Source PyMOL

Open-Source PyMOL is available [free of charge](https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE). It also allows sponsors to create highly customized PyMOL installations which might not be possible with the MSI installer. 

### One-click Installer

A convenient and easy-to-use one-click installer is available that runs out-of-the-box without the need for an existing Python installation. You can download it [here](https://github.com/kullik01/pymol-open-source-windows-setup/releases/download/v3.1.0/PyMOL_Open_source_v3.1.0a0_WINx64_setup.exe) or go to the [GitHub repository](https://github.com/kullik01/pymol-open-source-windows-setup) and then to [Releases](https://github.com/kullik01/pymol-open-source-windows-setup/releases/tag/v3.1.0). It features an alternate splash screen and a slightly updated menu stylesheet. 

### Pre-compiled Open-Source PyMOL

Pre-compiled Open-Source PyMOL is available free from [Christoph Gohlke of the Laboratory for Fluorescence Dynamics, University of California, Irvine](https://www.cgohlke.com/). 

  1. Install the latest version of Python 3 for Windows (e.g., by going to <http://www.python.org/downloads/> and choosing the x64 EXE installer). Use the standard options, which should mean that the installation directory is most likely C:\Users\<Your Username>\AppData\Local\Programs\Python\Python38). Make sure the option to add environment variables is selected or add the folder of python.exe to system PATH.
  2. Install [the current Microsoft Visual C++ Redistributable for Visual Studio 2015, 2017 and 2019](https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads). Otherwise the installed PyMOL binary may fail to run (without any error message!). 
     1. Have a look in "Add/Remove" programs if it's installed already
  3. Download the [appropriate wheel files](https://github.com/cgohlke/pymol-open-source-wheels/), along with all requirement wheel files ([Numpy+MKL](https://github.com/cgohlke/numpy-mkl-wheels/)) into a single file directory, e.g., `C:\Users\<Your Username>\Downloads`
  4. Example of filenames 2023-01-12 
     1. pymol_launcher-2.5-cp311-cp311-win_amd64.whl
     2. pymol-2.6.0a0-cp311-cp311-win_amd64.whl
     3. numpy-1.22.4+mkl-cp311-cp311-win_amd64.whl



Navigate to the installation directory in a CMD window (Not PowerShell!) 
    
    
    cd C:\Users\<Your Username>\Downloads
    

(or where ever you put the files) and begin the installation using the command: ("%CD%" only works in CMD) 
    
    
    python -m pip pmw
    python -m pip install --no-index --find-links="%CD%" pymol_launcher-2.5-cp311-cp311-win_amd64.whl
    where.exe pymol
    pymol
    

PyMOL.exe should now be in C:\Users\<Your Username>\AppData\Local\Programs\Python 

To use the newer single-window Qt interface, also install the optional PyQt5 dependency for your Python installation: 
    
    
    python -m pip install pyqt5
    

To update PyMOL update the files in the PyMOL install directory and run: 
    
    
    pip install --upgrade --no-deps pymol.whl
    

where `pymol.whl` is replaced by the PyMOL wheel file name (not the launcher, the launcher should not require updating). 

### Building from source

This [GitHub repository](https://github.com/urban233/pymol-open-source-windows-build) provides the tools to build PyMOL for Windows from source. It automatically downloads all necessary dependencies via the vcpkg package manager. You can either build the complete wheel file or just the _cmd C/C++ extension module using the provided automation script. A pre-built wheel file is available for download directly from the [Releases](https://github.com/urban233/pymol-open-source-windows-build/releases/tag/v3.1.0a0) tab. 

### Extend PyMOL with additional scripts

If you now want to extend the capabilities of PyMOL, and take advantage of all the available plugins+scripts "out there", then do the following.   


  1. First install "numpy" as an available module to Python. [Select appropriate installer from here](https://github.com/cgohlke/numpy-mkl-wheels/)
  2. Download the script/plugin collection [ Pymol-script-repo](/index.php/Git "Git") from [a .zip file from here](https://github.com/Pymol-Scripts/Pymol-script-repo/zipball/master)


    
    
    git clone <https://github.com/Pymol-Scripts/Pymol-script-repo>
    

  1. Unpack it to here: **C:\Python27\Lib\site-packages\pymol\pymol_path\Pymol-script-repo** Double check that the folder name is correct and the same.



Open "Notepad" and write. 
    
    
    # Add paths to sys.path so PyMOL can find modules and scripts
    import sys, os
    pymol_git = os.path.abspath(os.path.join(os.environ['PYMOL_PATH'], 'Pymol-script-repo'))
    os.environ['PYMOL_GIT_MOD'] = os.path.join(pymol_git,'modules')
    sys.path.append(pymol_git)
    sys.path.append(os.environ['PYMOL_GIT_MOD'])
    
    # Make setting changes to Plugin Manager
    import pymol.plugins
    pymol.plugins.preferences = {'instantsave': False, 'verbose': False}
    pymol.plugins.autoload = {'apbs_tools': False}
    pymol.plugins.set_startup_path([os.path.join(pymol_git, 'plugins'), os.path.join(sys.prefix, 'Lib', 'site-packages', 'pmg_tk', 'startup')])
    pymol.plugins.preferences = {'instantsave': True, 'verbose': False}
    

**Then "File- >Save as->All files-> C:\Python27\Lib\site-packages\pymol\pymol_path\run_on_startup.py**

Now start pymol, and enjoy all the plugins available from the menu. 

**PyMOL** shortcut  
Make a **pymol** directory in your homepath. **mkdir %HOMEPATH%\pymol** Then make sure, PyMOL starts here, when you open the shortcut.  
Make a shortcut to the .cmd file, and modify it.   
Target: C:\python27\PyMOL\pymol.cmd   
Start in: %HOMEPATH%\pymol 

## See Also

  * [pymolrc](/index.php/Pymolrc "Pymolrc")
  * [Linux Install](/index.php/Linux_Install "Linux Install")
  * [MAC Install](/index.php/MAC_Install "MAC Install")



Retrieved from "[https://pymolwiki.org/index.php?title=Windows_Install&oldid=13838](https://pymolwiki.org/index.php?title=Windows_Install&oldid=13838)"


---

