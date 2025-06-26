# Category: Representations

## Ball and Stick

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Ball and Stick
  * 2 Usage
    * 2.1 Builtin
      * 2.1.1 GUI
      * 2.1.2 Command Line/API
    * 2.2 Hand Made



## Ball and Stick

The **ball and stick** representation is very often used to display macromolecules. PyMOL allows the user the ability to turn on this representation for certain selections, or roll their own hand-made versions of the command (see below). 

  * [![Example of Bal and Stick Repr.](/images/9/98/Ball_stick.png)](/index.php/File:Ball_stick.png "Example of Bal and Stick Repr.")

Example of Bal and Stick Repr. 




## Usage

### Builtin

#### GUI

Using the mouse, for the desired object or selection click **A > Preset > ball and stick**. 

#### Command Line/API
    
    
    # turn on the representation for everything in PyMOL
    preset.ball_and_stick(selection='all', mode=1)
    

or... (keeping in mind that the [single-word selectors](/index.php/Single-word_selectors "Single-word selectors") page includes more possible selections to substitute `r. STI` with) 
    
    
    # turn on representation for one ligand (in this example STI or imatinib) in PyMOL
    preset.ball_and_stick(selection='r. STI', mode=1)
    

Parameter _mode_ can be 1 or 2. With _mode=1_ you'll have... 
    
    
    set_bond stick_color, white, selection, selection
    set_bond stick_radius, 0.14, selection, selection
    set sphere_scale, 0.25, selection
    show sticks, selection
    show spheres, selection
    

If _mode=2_ than you'll get... 
    
    
    set_bond stick_color, white, selection, selection
    set_bond stick_radius, -0.14, selection, selection
    set stick_ball, 1
    set stick_ball_ratio, -1.0
    set stick_ball_color, atomic
    show sticks, selection
    

### Hand Made

Balls with sticks really look nice. You can even create your own style of this, with control over **[sphere_scale](/index.php/Sphere_scale "Sphere scale")** and **[stick_radius](/index.php/Stick_radius "Stick radius")**. 
    
    
    hide lines
    show sticks
    show spheres
    set stick_radius, 0.1, (all)
    set sphere_scale, 0.25, (all)
    

Also OK: 
    
    
    show sticks
    set valence, on
    set stick_ball, on
    set stick_ball_ratio, 3
    set stick_radius, 0.12
    

You can adjust the two numbers above to your taste. 

* * *

Ed. As of 0.98bXX there is a GUI-enable Ball & Stick representation available to users. [Tree](/index.php/User:Tree "User:Tree")

Retrieved from "[https://pymolwiki.org/index.php?title=Ball_and_Stick&oldid=11961](https://pymolwiki.org/index.php?title=Ball_and_Stick&oldid=11961)"


---

## Cartoon

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Cartoon Command
    * 1.1 DESCRIPTION
    * 1.2 USAGE
    * 1.3 PYMOL API
    * 1.4 EXAMPLES
    * 1.5 NOTES
  * 2 Adjusting width of cartoon
    * 2.1 Forcing Exact Boundaries in Coloring Secondary Structures
  * 3 Sausage Representation
  * 4 Black and White Representation
  * 5 CA (Alpha Carbon) Trace
  * 6 Various Transparency Levels
  * 7 Nucleic Acid Representation
    * 7.1 Other Nucleic Acids & Cartoon Settings
  * 8 Breaking Up a Visualization with Cartoon Skip
  * 9 See Also



## Cartoon Command

### DESCRIPTION

**cartoon** changes the default cartoon for a set of atoms. 

### USAGE
    
    
     cartoon type, (selection)
    

**type:**

  * skip (-1)
  * automatic (0) (use [ss](/index.php/Dss "Dss") property)
  * loop (1)
  * rectangle (2)
  * oval (3)
  * tube (4)
  * arrow (5)
  * dumbbell (6) (see also [cartoon_fancy_helices](/index.php/Cartoon_fancy_helices "Cartoon fancy helices"))
  * putty (7) (b-factor scaling)
  * dash (8) _(new in PyMOL 1.8.2)_


  * [![cartoon_skip; yes, blank.](/images/2/29/Cartoon_skip.png)](/index.php/File:Cartoon_skip.png "cartoon_skip; yes, blank.")

cartoon_skip; yes, blank. 

  * [![cartoon_automatic](/images/5/5f/Cartoon_automatic.png)](/index.php/File:Cartoon_automatic.png "cartoon_automatic")

cartoon_automatic 

  * [![cartoon_loop](/images/f/f0/Cartoon_loop.png)](/index.php/File:Cartoon_loop.png "cartoon_loop")

cartoon_loop 

  * [![cartoon_rectangle](/images/9/9c/Cartoon_rectangle.png)](/index.php/File:Cartoon_rectangle.png "cartoon_rectangle")

cartoon_rectangle 

  * [![cartoon_oval](/images/8/81/Cartoon_oval.png)](/index.php/File:Cartoon_oval.png "cartoon_oval")

cartoon_oval 

  * [![cartoon_tube](/images/8/85/Cartoon_tube.png)](/index.php/File:Cartoon_tube.png "cartoon_tube")

cartoon_tube 

  * [![cartoon_arrow](/images/1/16/Cartoon_arrow.png)](/index.php/File:Cartoon_arrow.png "cartoon_arrow")

cartoon_arrow 

  * [![cartoon_dumbbell](/images/c/ca/Cartoon_dumbbell.png)](/index.php/File:Cartoon_dumbbell.png "cartoon_dumbbell")

cartoon_dumbbell 




### PYMOL API
    
    
    cmd.cartoon(string type, string selection )
    

### EXAMPLES
    
    
    cartoon rectangle,(chain A)
    cartoon skip,(resi 145:156)
    

### NOTES

the "automatic" mode utilizes ribbons according to the information in the PDB HELIX and SHEET records. 

## Adjusting width of cartoon

Try varying the following. 

For β-strands: 
    
    
    cartoon_rect_length
    cartoon_rect_width
    

For α-helices: 
    
    
    cartoon_oval_length
    cartoon_oval_width
    

For loops: 
    
    
    cartoon_loop_radius
    

For nucleic acid backbones which resemble 'loops'; however, are not classified as such by PyMOL (see more about nucleic acid representation settings at bottom of page): 
    
    
    cartoon_tube_radius,0.8
    

For "fancy" α-helices: 
    
    
    cartoon_dumbbell_length
    cartoon_dumbbell_width
    cartoon_dumbbell_radius  (radius of cylinder at edge of helix ribbon)
    

  
In each case "length" refers to what some might call the width and "width" refers to what some might call the thickness. 

  


[![](/images/9/9f/Cartoon_ex.png)](/index.php/File:Cartoon_ex.png)

[](/index.php/File:Cartoon_ex.png "Enlarge")

Cartoon Representation Example

  


### Forcing Exact Boundaries in Coloring Secondary Structures

To force PyMOL to respect secondary structural elements color-wise (PyMOL smooths out colors near color chagnes for a prettier image) use the following PyMOL command: ` set cartoon_discrete_colors, on `

[![](/images/c/c7/Cartoon_discrete_color0.png)](/index.php/File:Cartoon_discrete_color0.png)

[](/index.php/File:Cartoon_discrete_color0.png "Enlarge")

Discrete Coloring Off

[![](/images/2/2d/Cartoon_discrete_color1.png)](/index.php/File:Cartoon_discrete_color1.png)

[](/index.php/File:Cartoon_discrete_color1.png "Enlarge")

Discrete Coloring On

## Sausage Representation

The familiar sausage representation in PyMOL is called, "putty". To enable the putty/sausage view simply do, 
    
    
    show cartoon
    cartoon putty
    unset cartoon_smooth_loops
    unset cartoon_flat_sheets
    

As of v 0.98 or so, there's a Putty option. Use this. 

[![](/images/4/40/B_factor_putty.png)](/index.php/File:B_factor_putty.png)

[](/index.php/File:B_factor_putty.png "Enlarge")

Example of B-factor Putty

## Black and White Representation

**UPDATE** : This method is essentially obseleted by the new setting **set ray_trace_mode,2**. More information on this at [Ray](/index.php/Ray "Ray"). For those who want a nifty black and white representation of their protein try the following: 

  1. Ray trace your protein of choice in a cartoon representation use a BLACK background
  2. Save the image
  3. Load the image in GIMP. 

[![](/images/d/d3/Bw1.jpeg)](/index.php/File:Bw1.jpeg)

[](/index.php/File:Bw1.jpeg "Enlarge")

Black BG Ribbon

  4. Desaturate or Grayscale the image. 

[![](/images/d/dc/Bw2.jpeg)](/index.php/File:Bw2.jpeg)

[](/index.php/File:Bw2.jpeg "Enlarge")

Grayscale

  5. Run the filter: Filter->Edge-Detect->Edge. 

[![](/images/2/20/Bw3.jpeg)](/index.php/File:Bw3.jpeg)

[](/index.php/File:Bw3.jpeg "Enlarge")

Edge Detect

  6. Select: Layers->Color->Invert. 

[![](/images/c/c1/Bw4.jpeg)](/index.php/File:Bw4.jpeg)

[](/index.php/File:Bw4.jpeg "Enlarge")

Invert Color

  7. Different methods of edge detection will give you different results. In the last example, I used Laplace Edge-Detect, then painted an all white layer beneath the current layer to achieve the results. 

[![](/images/7/7c/Bw5.jpeg)](/index.php/File:Bw5.jpeg)

[](/index.php/File:Bw5.jpeg "Enlarge")

Comments




  
I'm sure there are other ways to do this. If you want to include it in a publication make sure you ray traced it large enough. For that, see [Creating Publication Quality Images](/index.php/Publication_Quality_Images "Publication Quality Images"). 

## CA (Alpha Carbon) Trace

If you have a structure with just a alpha carbon trace, you can get a cartoon by 
    
    
    set cartoon_trace,1
    show cartoon
    

If your structure is more than just the CA backbone, the cartoon representation will look incorrect, so use it just with CA trace. 

## Various Transparency Levels

[![](/images/4/43/Cartoon_multi_transp.png)](/index.php/File:Cartoon_multi_transp.png)

[](/index.php/File:Cartoon_multi_transp.png "Enlarge")

Example of Cartoon Multi-level Transparency. The near cartoon has transparency setting **0.2** , the segment in the BG **0.5**.

One can make different cartoon selections have different transparency values, in PyMOL. The trick here is to use "create" or "extract" instead of "select". Create makes new objects that can have independent settings while leaving the original object intact, whereas extract removes the specified atoms from the original object when creating the new object. 
    
    
    load mol_obj.pdb
    
    # transfer a piece of the molecule into a new object
    extract new_obj, chain A
    
    # adjust trasparency for the new object
    set cartoon_transparency, 0.5, new_obj
    

## Nucleic Acid Representation

[![](/images/c/cb/Nucleic1.png)](/index.php/File:Nucleic1.png)

[](/index.php/File:Nucleic1.png "Enlarge")

Showing Nucleic Acids

To control radius of nucleic acids default backbone cartoon: 
    
    
    set cartoon_tube_radius,0.8   #0.5 seems close to the default setting
    

To show nucleic acids in a nicer format do: 
    
    
     set cartoon_ring_mode,1
     show cartoon
    

  


### Other Nucleic Acids & Cartoon Settings

Here are some things to try: 
    
    
    set cartoon_ring_mode, 1   # (or 2 or 3)
    set cartoon_ring_finder, 1 # (or 2 or 3 or 4)
    set cartoon_nucleic_acid_mode, 0   # (or 1 or 2 or 3 or 4)
    
    set cartoon_side_chain_helper
    rebuild
    
    set cartoon_ring_transparency, 0.5
    
    set cartoon_ladder_mode, 0 # or 1
    
    set cartoon_ladder_color, color-name
    set cartoon_nucleic_acid_color, color-name
    
    cartoon oval
    set cartoon_oval_width, 0.8
    
    cartoon rect
    
    cartoon dumbbell
    set cartoon_dumbbell_width, 0.4
    set cartoon_dumbbell_radius, 0.4
    

[Overview of nucleic acid cartoons](/index.php/Overview_of_nucleic_acid_cartoons "Overview of nucleic acid cartoons")

[Examples of nucleic acid cartoons](/index.php/Examples_of_nucleic_acid_cartoons "Examples of nucleic acid cartoons")

## Breaking Up a Visualization with Cartoon Skip
    
    
    fetch 1cll, async=0
    
    create a_copy, 1cll
    
    # show them as cartoons
    
    color blue, 1cll
    
    set cartoon_highlight_color, green, 1cll
    
    color marine, a_copy
    
    set cartoon_highlight_color, magenta, a_copy
    
    # select the region of interest
    
    select to_hide, i. 0-79
    
    show_as cartoon
    
    cartoon skip, to_hide and 1cll
    
    orient
    
    # hide/skip this section in the original
    
    cartoon skip, i. 81-200 and a_copy
    
    ray
    

  


## See Also

[Displaying_Biochemical_Properties](/index.php/Displaying_Biochemical_Properties "Displaying Biochemical Properties")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon&oldid=12708](https://pymolwiki.org/index.php?title=Cartoon&oldid=12708)"


---

## Dots

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/d/dd/Dots_ex.png)](/index.php/File:Dots_ex.png)

[](/index.php/File:Dots_ex.png "Enlarge")

Dots Representation Example

## Contents

  * 1 Overview
    * 1.1 Dots Settings
      * 1.1.1 dot_radius (dot weight)
      * 1.1.2 dot_width (dot weight)
      * 1.1.3 dot_density (dot density)
    * 1.2 Possible Other Settings
    * 1.3 Unused Settings
    * 1.4 See Also



# Overview

A simple PyMol representation. Use ` `

`
    
    
    show dots, SELECTION
    

```

`` where SELECTION is a valid selection or previously defined selection name. 

  


## Dots Settings

The differences in dot_radius and dot_width arise from the differences between the rasterization of OpenGL and ray tracing. PyMOL line and dot widths have the name meaning as their underlying OpenGL primitives (pixel widths), and they dominate over line and dot radii settings such that you get WYSIWYG raytracing output which is faithful to the OpenGL view in terms of cross-sectional area. If you set the radii values explicitly (Angstroms), then the ray tracer will be affected, but not the OpenGL viewer. 

### dot_radius (dot weight)

Controls the dot radius when ray tracing (that is, it is apparent only after the scene is ray traced). This option can take any value, values from 0 to whatever (4 is huge!) If dot_radius is set to zero, then it is subordinated to dot_width; if dot_radius is a non-zero value, then it will take priority. 

To check/change dots_radius: 
    
    
    set dot_radius, VALUE
    # equivalent to
    set dot_radius=VALUE
    
    # to get the current setting, use
    get dot_radius
    

  * The default value is: 0



### dot_width (dot weight)

This setting defines the size of the dots, in (screen) pixels. It is subordinate to dot_radius, that is dot_radius takes precedence over dot_width unless dot_radius=0. 
    
    
    set dot_width, VALUE
    # equivalent to
    set dot_width=VALUE
    
    # to get the current setting, use
    get dot_width
    

  * The default value is: 2



### dot_density (dot density)

Determines the density of dots in the dot representation; the existing dot densities arise from successive cycles of isocosohedral interpolation.. It takes only integer values from 0 to 4. 
    
    
    set dot_density, VALUE
    # equivalent to
    set dot_density=VALUE
    
    # to get the current setting, use
    get dot_density
    

  * The default value is: 2



## Possible Other Settings

  * dot_solvent, default is off. This setting controls whether the dots are representing the vdW surface; when the setting is on, it's the solvent accessible surface area.
  * dot_color, default is default
  * cgo_dot_width, default is 2.00000
  * cgo_dot_radius, default is -1.00000
  * dot_normals, default is on
  * dot_lighting, default is on



## Unused Settings

  * dot_mode, default is 0



  


## See Also

[Cmd isodot](/index.php/Cmd_isodot "Cmd isodot")

Retrieved from "[https://pymolwiki.org/index.php?title=Dots&oldid=6572](https://pymolwiki.org/index.php?title=Dots&oldid=6572)"


---

## Ellipsoids

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

The **Ellipsoids** representation draws the anisotropic (ANISOU fields in the PDB record) thermal ellipsoids from high resolution structures. It is the newest addition to the PyMOL [Category:Representations](/index.php/Category:Representations "Category:Representations"). 

  * [![Showing ellipsoids](/images/d/d6/Ell1.png)](/index.php/File:Ell1.png "Showing ellipsoids")

Showing ellipsoids 

  * [![Showing ellipsoids and main chain as sticks.](/images/c/c3/Ell2.png)](/index.php/File:Ell2.png "Showing ellipsoids and main chain as sticks.")

Showing ellipsoids and main chain as sticks. 




## Syntax
    
    
    show ellipsoids, SEL
    

## Example
    
    
    fetch 1tjx
    show ellipsoids, 1tjx
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoids&oldid=4580](https://pymolwiki.org/index.php?title=Ellipsoids&oldid=4580)"


---

## Examples of nucleic acid cartoons

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Default settings
  * 2 Cartoon ring mode
    * 2.1 Settings
    * 2.2 Examples
  * 3 Cartoon ring finder
    * 3.1 Settings
    * 3.2 Examples
  * 4 Cartoon ladder mode
    * 4.1 Settings
    * 4.2 Examples
  * 5 Cartoon nucleic acid mode
    * 5.1 Settings
    * 5.2 Examples



## Default settings

The defaults give a phosphate backbone with single sticks passing across the full width of the base plane. 
    
    
    set cartoon_nucleic_acid_mode, 0 # backbone follows phosphates; actually Pymol itself uses setting '4' as default
    set cartoon_ladder_mode, 1 # sticks from backbone into nucleotide
    set cartoon_ring_mode, 0 # no nucleotide rings
    set cartoon_ring_finder, 1 # ribose and base rings (not displayed since ring mode 0)
    

[![](/images/b/bd/DNA-default-ring0-ladder1-na0-finder1.png)](/index.php/File:DNA-default-ring0-ladder1-na0-finder1.png) [](/index.php/File:DNA-default-ring0-ladder1-na0-finder1.png "Enlarge")default view  
---  
  
## Cartoon ring mode

### Settings
    
    
    set cartoon_ring_mode, value

value | effect   
---|---  
0 | stick from backbone atom to N1 of purines and N3 of pyrimidines   
1 | simple plane for ribose and base rings covering area between ring bonds   
2 | simple plane for ribose and base rings covering area inside sticks (slightly smaller than mode 1)   
3 | plane bounded by sticks for ribose and base rings   
4 | large sphere of ring diameter at centre of ribose and each base ring   
5 | small sphere of 1/10 diameter at centre of ribose and each base ring   
  
### Examples

[![](/images/b/bb/DNA-ring0-ladder1-na0-finder1.png)](/index.php/File:DNA-ring0-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring0-ladder1-na0-finder1.png "Enlarge")cartoon_ring_mode,0 | [![](/images/b/b1/DNA-ring1-ladder1-na0-finder1.png)](/index.php/File:DNA-ring1-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring1-ladder1-na0-finder1.png "Enlarge")cartoon_ring_mode,1 | [![](/images/5/5c/DNA-ring2-ladder1-na0-finder1.png)](/index.php/File:DNA-ring2-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring2-ladder1-na0-finder1.png "Enlarge")cartoon_ring_mode,2 | [![](/images/d/de/DNA-ring3-ladder1-na0-finder1.png)](/index.php/File:DNA-ring3-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder1.png "Enlarge")cartoon_ring_mode,3 | [![](/images/8/88/DNA-ring4-ladder1-na0-finder1.png)](/index.php/File:DNA-ring4-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring4-ladder1-na0-finder1.png "Enlarge")cartoon_ring_mode,4 | [![](/images/b/bb/DNA-ring5-ladder1-na0-finder1.png)](/index.php/File:DNA-ring5-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring5-ladder1-na0-finder1.png "Enlarge")cartoon_ring_mode,5  
---|---|---|---|---|---  
all with defaults:  _cartoon_ladder_mode,1 cartoon_nucleic_acid_mode,0 cartoon_ring_finder,1_

## Cartoon ring finder

### Settings
    
    
    set cartoon_ring_finder, value

value | effect   
---|---  
0 | no rings or sticks joining them   
1 | both ribose and base ring   
2 | only base ring(s), stick connects directly from phosphate to ring   
3 | very similar to ring finder 1, slight effect on transparency = distinct behaviour?   
4 | very similar to ring finder 1: finds ribose and base of nucleotides, and aromatic side chains of proteins   
5 | sticks visible but rings invisible   
  
### Examples

[![](/images/d/d3/DNA-ring3-ladder1-na0-finder0.png)](/index.php/File:DNA-ring3-ladder1-na0-finder0.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder0.png "Enlarge")cartoon_ring_finder,0 | [![](/images/d/de/DNA-ring3-ladder1-na0-finder1.png)](/index.php/File:DNA-ring3-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder1.png "Enlarge")cartoon_ring_finder,1 | [![](/images/e/eb/DNA-ring3-ladder1-na0-finder2.png)](/index.php/File:DNA-ring3-ladder1-na0-finder2.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder2.png "Enlarge")cartoon_ring_finder,2 | [![](/images/8/89/DNA-ring3-ladder1-na0-finder3.png)](/index.php/File:DNA-ring3-ladder1-na0-finder3.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder3.png "Enlarge")cartoon_ring_finder,3 | [![](/images/3/36/DNA-ring3-ladder1-na0-finder4.png)](/index.php/File:DNA-ring3-ladder1-na0-finder4.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder4.png "Enlarge")cartoon_ring_finder,4 | [![](/images/0/09/DNA-ring3-ladder1-na0-finder5.png)](/index.php/File:DNA-ring3-ladder1-na0-finder5.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder5.png "Enlarge")cartoon_ring_finder,5  
---|---|---|---|---|---  
  
all with:  _cartoon_ladder_mode,1 cartoon_ring_mode,3 cartoon_nucleic_acid_mode,0_

## Cartoon ladder mode

### Settings
    
    
    set cartoon_ladder_mode, value

value | effect   
---|---  
0 | no sticks shown   
1 | sticks show   
  
### Examples

[![](/images/2/22/DNA-ring3-ladder0-na0-finder1.png)](/index.php/File:DNA-ring3-ladder0-na0-finder1.png) [](/index.php/File:DNA-ring3-ladder0-na0-finder1.png "Enlarge")cartoon_ladder_mode,0 | [![](/images/d/de/DNA-ring3-ladder1-na0-finder1.png)](/index.php/File:DNA-ring3-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder1.png "Enlarge")cartoon_ladder_mode,1  
---|---  
  
all with:  _cartoon_ring_mode,3 cartoon_nucleic_acid_mode,0 cartoon_ring_finder,1_

Note that the visibility of the ladder sticks depends on ring mode >0, ring finder >0, nucleic acid mode = 0 

## Cartoon nucleic acid mode

### Settings
    
    
    set cartoon_nucleic_acid_mode, value

value | effect   
---|---  
0 | smooth backbone passing through phosphorus atoms, backbone terminates at last phosphorus on either end of chain   
1 | smooth backbone passing through ribose C3' atoms, backbone terminates at last C3' on either end of chain   
2 | smooth backbone passing through phosphorus atoms, backbone terminates at last phosphorus on 5' end and O3' on 3' end (note takes O3' colour at terminus in default colouring)   
3 | appears same as mode 0?   
4 | appears same as mode 2? Seems to be what Pymol uses when it first opens nucleic acid containing file because any other settings change ends and colors.   
  
### Examples

[![](/images/d/de/DNA-ring3-ladder1-na0-finder1.png)](/index.php/File:DNA-ring3-ladder1-na0-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na0-finder1.png "Enlarge")cartoon_nucleic_acid_mode,0 | [![](/images/3/33/DNA-ring3-ladder1-na1-finder1.png)](/index.php/File:DNA-ring3-ladder1-na1-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na1-finder1.png "Enlarge")cartoon_nucleic_acid_mode,1 | [![](/images/4/4f/DNA-ring3-ladder1-na2-finder1.png)](/index.php/File:DNA-ring3-ladder1-na2-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na2-finder1.png "Enlarge")cartoon_nucleic_acid_mode,2 | [![](/images/c/c4/DNA-ring3-ladder1-na3-finder1.png)](/index.php/File:DNA-ring3-ladder1-na3-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na3-finder1.png "Enlarge")cartoon_nucleic_acid_mode,3 | [![](/images/5/53/DNA-ring3-ladder1-na4-finder1.png)](/index.php/File:DNA-ring3-ladder1-na4-finder1.png) [](/index.php/File:DNA-ring3-ladder1-na4-finder1.png "Enlarge")cartoon_nucleic_acid_mode,4  
---|---|---|---|---  
  
all with:  _cartoon_ladder_mode,0 cartoon_ring_mode,3 cartoon_ring_finder,1_

Retrieved from "[https://pymolwiki.org/index.php?title=Examples_of_nucleic_acid_cartoons&oldid=7997](https://pymolwiki.org/index.php?title=Examples_of_nucleic_acid_cartoons&oldid=7997)"


---

## Light

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 User Hint
    * 3.1 The Code
  * 4 See Also



# Overview

Lighting is important for high-quality shots. PyMOL supports of up to 10 virtual lights. You can turn the lights on/off and also position them where you want (behind the camera). 

Light2..Light9 control where the other lights go. 

  * [![Normal lighting](/images/7/76/L1.png)](/index.php/File:L1.png "Normal lighting")

Normal lighting 

  * [![Light1 moved a bit. Notice how the shadows have changed.](/images/9/93/L2.png)](/index.php/File:L2.png "Light1 moved a bit. Notice how the shadows have changed.")

Light1 moved a bit. Notice how the shadows have changed. 




# Syntax
    
    
    # set the light to some position.  The 'position'
    # must be a vector specifying the XYZ location
    # to put the light.
    set light, position
    
    # for example
    set light, [ -0.55, -0.70, 0.15 ]
    

# User Hint

  * [![Moving lights make for a cool effect; you can make this look better by smoothing out the distances and ranges.](/images/1/1b/Ll.gif)](/index.php/File:Ll.gif "Moving lights make for a cool effect; you can make this look better by smoothing out the distances and ranges.")

Moving lights make for a cool effect; you can make this look better by smoothing out the distances and ranges. 




: 

One neat trick, for rendering a "sunset" on your protein is to turn off all the lights, then render the scene as you move the light across the scene. The shadows move across the protein based on the light position and it looks like the sun is setting. 

  


## The Code

Here's the code for the animated GIF shown above. 
    
    
    python
    cmd.set("light", ll)
    for x in range(10):
            l = [ -0.4 + 2*float(x/10.), -0.4, -1  ]
            print l
            cmd.set("light", l)
            cmd.ray()
            cmd.png("lll" + str(x) + ".png" )
    for x in range(10):
            l[0] -= 2*float(x/10.)
            print l
            cmd.set("light", l)
            cmd.ray()
            cmd.png("llll" + str(x) + ".png" )
    python end
    

Then, in the shell do 
    
    
    convert lll* llll* light_movie.gif
    

# See Also

[Ray](/index.php/Ray "Ray")

Retrieved from "[https://pymolwiki.org/index.php?title=Light&oldid=6152](https://pymolwiki.org/index.php?title=Light&oldid=6152)"


---

## Lines

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Lines** is name of the basic representation for atoms and bonds in PyMOL. **Lines** is a very simple representation, where each atom bond is displayed as a single colored line, and each atom is displayed as the intersection of any two or more non-terminal bonds. 

## Contents

  * 1 Usage
    * 1.1 Examples
      * 1.1.1 Example: Displaying dashed lines between two atoms
    * 1.2 See Also



# Usage
    
    
    # show everything as lines
    show lines
    
    # only show residues 50-80 as lines
    show lines, i.50-80
    

## Examples

### Example: Displaying dashed lines between two atoms

The following commands will create a dashed line between two atoms. 
    
    
    # first, create two named selections
    select a, ///A/501/02
    select b, ///B/229/N
    # calculate & show the distance from selection a to selection b.
    distance d, a, b
    # hide just the distance labels; the 
    # dashed bars should still be shown
    hide labels, d
    

Technically, the object _d_ is a labelled distance, only the label is hidden. When ray-tracing the image, the dashes come out a bit fat. You can slim them with 
    
    
    set dash_gap, 0.5
    set dash_radius, 0.1
    

before the 'ray' command. 

[![](/images/d/d4/Lines_ex.png)](/index.php/File:Lines_ex.png)

[](/index.php/File:Lines_ex.png "Enlarge")

Lines Representation Example

## See Also

Please read about other representations in the **[Representation Category](/index.php/Category:Representations "Category:Representations")**.   
[Measure_Distance](/index.php/Measure_Distance "Measure Distance")   
[Distance](/index.php/Distance "Distance")   


Retrieved from "[https://pymolwiki.org/index.php?title=Lines&oldid=7527](https://pymolwiki.org/index.php?title=Lines&oldid=7527)"


---

## Mesh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Example
    * 3.1 Settings
    * 3.2 Notes
    * 3.3 See Also



# Overview

PyMOL has a web-like, Mesh representation, as shown in the example images below. 

Using the internal GUI, you may enable mesh by clicking the box with the letter **S** in it and selecting **mesh**. 

# Syntax
    
    
    # show the selection, SEL, in mesh
    show mesh, SEL
    
    # using the 'as' keyword
    as mesh, SEL
    

# Example
    
    
    select B, i. 40-110
    show mesh, B
    

  


  * [![Mesh Representation Example](/images/3/32/Mesh_ex.png)](/index.php/File:Mesh_ex.png "Mesh Representation Example")

Mesh Representation Example 

  * [![Mesh in ray_trace_mode,0](/images/4/48/Mesh_rm0.png)](/index.php/File:Mesh_rm0.png "Mesh in ray_trace_mode,0")

Mesh in ray_trace_mode,0 

  * [![Mesh in ray_trace_mode,3; burnt and probably not what you're looking for.](/images/5/5a/Mesh_rm3.png)](/index.php/File:Mesh_rm3.png "Mesh in ray_trace_mode,3; burnt and probably not what you're looking for.")

Mesh in ray_trace_mode,3; burnt and probably not what you're looking for. 

  * [![Mesh with width of 0.5, half of normal](/images/6/6b/Mesh_w05.png)](/index.php/File:Mesh_w05.png "Mesh with width of 0.5, half of normal")

Mesh with width of 0.5, half of normal 

  * [![Using mesh is a great way to show enclosed pockes.](/images/4/49/5ABP.png)](/index.php/File:5ABP.png "Using mesh is a great way to show enclosed pockes.")

Using mesh is a great way to show enclosed pockes. 




## Settings

  * [cavity_cull](/index.php/Cavity_cull "Cavity cull")
  * [mesh_carve_cutoff](/index.php?title=Mesh_carve_cutoff&action=edit&redlink=1 "Mesh carve cutoff \(page does not exist\)")
  * [mesh_cutoff](/index.php?title=Mesh_cutoff&action=edit&redlink=1 "Mesh cutoff \(page does not exist\)")
  * [mesh_quality](/index.php/Mesh_quality "Mesh quality")
  * [mesh_carve_selection](/index.php?title=Mesh_carve_selection&action=edit&redlink=1 "Mesh carve selection \(page does not exist\)")
  * [mesh_grid_max](/index.php?title=Mesh_grid_max&action=edit&redlink=1 "Mesh grid max \(page does not exist\)")
  * [mesh_radius](/index.php?title=Mesh_radius&action=edit&redlink=1 "Mesh radius \(page does not exist\)")
  * [mesh_carve_state](/index.php?title=Mesh_carve_state&action=edit&redlink=1 "Mesh carve state \(page does not exist\)")
  * [mesh_lighting](/index.php/Mesh_lighting "Mesh lighting")
  * [mesh_skip](/index.php?title=Mesh_skip&action=edit&redlink=1 "Mesh skip \(page does not exist\)")
  * [mesh_clear_cutoff](/index.php?title=Mesh_clear_cutoff&action=edit&redlink=1 "Mesh clear cutoff \(page does not exist\)")
  * [mesh_mode](/index.php/Mesh_mode "Mesh mode")
  * [mesh_solvent](/index.php?title=Mesh_solvent&action=edit&redlink=1 "Mesh solvent \(page does not exist\)")
  * [mesh_clear_selection](/index.php?title=Mesh_clear_selection&action=edit&redlink=1 "Mesh clear selection \(page does not exist\)")
  * [mesh_negative_color](/index.php/Mesh_negative_color "Mesh negative color")
  * [mesh_type](/index.php/Mesh_type "Mesh type")
  * [mesh_clear_state](/index.php?title=Mesh_clear_state&action=edit&redlink=1 "Mesh clear state \(page does not exist\)")
  * [mesh_negative_visible](/index.php/Mesh_negative_visible "Mesh negative visible")
  * [mesh_width](/index.php/Mesh_width "Mesh width")
  * [mesh_color](/index.php/Mesh_color "Mesh color")
  * [mesh_normals](/index.php?title=Mesh_normals&action=edit&redlink=1 "Mesh normals \(page does not exist\)")



## Notes

  * Mesh doesn't ray trace well in _set [ray_trace_mode]], 3_. Try setting this to _0_ if you mesh looks wonky.



## See Also

[Isomesh](/index.php/Isomesh "Isomesh") [Surface](/index.php/Surface "Surface")

Retrieved from "[https://pymolwiki.org/index.php?title=Mesh&oldid=8388](https://pymolwiki.org/index.php?title=Mesh&oldid=8388)"


---

## Ribbon

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Ribbon CA (Alpha Carbon) Trace

If you have a structure with just a alpha carbon trace, you can get a cartoon by 
    
    
    set ribbon_trace_atoms,1
    show ribbon
    

If your structure is more than just the CA backbone, the cartoon representation will look incorrect, so use it just with CA trace. 

Retrieved from "[https://pymolwiki.org/index.php?title=Ribbon&oldid=4419](https://pymolwiki.org/index.php?title=Ribbon&oldid=4419)"


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

## Sticks

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
    * 1.1 Settings
      * 1.1.1 Example Settings
      * 1.1.2 Color Sticks
        * 1.1.2.1 Sticks Radius (Sticks Weight)
        * 1.1.2.2 Sticks Transparency



# Overview

A simple PyMol representation where bonds are drawn as sticks. Use 
    
    
    # using the show command, for some SELECTION
    show sticks, SELECTION
    
    # using the as command
    as sticks, SELECTION
    

where SELECTION is a valid selection or previously defined selection name. 

  * [![Example Sticks Representation](/images/7/77/Sticks_ex.png)](/index.php/File:Sticks_ex.png "Example Sticks Representation")

Example Sticks Representation 




## Settings

  * [stick_ball](/index.php/Stick_ball "Stick ball")
  * [stick_nub](/index.php/Stick_nub "Stick nub")
  * [stick_transparency](/index.php/Stick_transparency "Stick transparency")
  * [stick_ball_ratio](/index.php/Stick_ball_ratio "Stick ball ratio")
  * [stick_overlap](/index.php?title=Stick_overlap&action=edit&redlink=1 "Stick overlap \(page does not exist\)")
  * [stick_valence_scale](/index.php?title=Stick_valence_scale&action=edit&redlink=1 "Stick valence scale \(page does not exist\)")
  * [stick_color](/index.php/Stick_color "Stick color")
  * [stick_quality](/index.php?title=Stick_quality&action=edit&redlink=1 "Stick quality \(page does not exist\)")
  * [stick_fixed_radius](/index.php?title=Stick_fixed_radius&action=edit&redlink=1 "Stick fixed radius \(page does not exist\)")
  * [stick_radius](/index.php/Stick_radius "Stick radius")
  * [set_bond](/index.php/Set_bond "Set bond")



### Example Settings

### Color Sticks

Use [set_bond](/index.php/Set_bond "Set bond") to set stick-bond settings, like color: 
    
    
    set_bond 1foo and i. XYZ, color red
    

#### Sticks Radius (Sticks Weight)

To change the radius for sticks, enter the following: 
    
    
    set stick_radius, VALUE
    

where 
    
    
    0.0<=VALUE<=1.0
    

  * **1.0** is 100% or full radius
  * **0.0** is 0% or invisible -- so use at least 0.1 or greater
  * The default value is: ~0.3



#### Sticks Transparency

To enable transparency for sticks, enter the following: 
    
    
    set stick_transparency, VALUE
    

where _0.0 <=VALUE<=1.0_

  * **1.0** is 100% transparent -- so invisible
  * **0.0** is 0% transparent -- so opaque


    
    
    set stick_transparency, 0.45
    

Retrieved from "[https://pymolwiki.org/index.php?title=Sticks&oldid=6996](https://pymolwiki.org/index.php?title=Sticks&oldid=6996)"


---

## Surface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The surface representation of a protein, in PyMol, shows the ["Connolly" surface](http://en.wikipedia.org/wiki/Connolly_surface) or the surface that would be traced out by the **surfaces** of waters in contact with the protein at all possible positions. 

[![](/images/7/79/Surface_ex.png)](/index.php/File:Surface_ex.png)

[](/index.php/File:Surface_ex.png "Enlarge")

Surface Representation Example

  


## Contents

  * 1 Enabling
  * 2 Settings
    * 2.1 Examples
      * 2.1.1 Transparency
      * 2.1.2 Quality
      * 2.1.3 Probe Radius
  * 3 Tips
    * 3.1 Exporting Surface/Mesh Coordinates to File
      * 3.1.1 Older PyMOL Versions
      * 3.1.2 Newer PyMOL Versions
    * 3.2 Representation-independent Color Control
    * 3.3 Displaying a protein as surface with a ligand as sticks
    * 3.4 Calculating a partial surface
    * 3.5 Displaying surface inside a molecule
    * 3.6 Creating a Closed Surface
    * 3.7 Smooth surface quick (blob)
    * 3.8 Smooth surface accurate (blob)
    * 3.9 Huge Surfaces
  * 4 Performance



## Enabling

To enable the surface representation do 
    
    
    show surface, SEL
    

for any proper selection SEL. 

## Settings

  * [cavity_cull](/index.php/Cavity_cull "Cavity cull")
  * [surface_best](/index.php?title=Surface_best&action=edit&redlink=1 "Surface best \(page does not exist\)")
  * [surface_negative_color](/index.php/Surface_negative_color "Surface negative color")
  * [surface_carve_cutoff](/index.php/Surface_carve_cutoff "Surface carve cutoff")
  * [surface_negative_visible](/index.php/Surface_negative_visible "Surface negative visible")
  * [surface_carve_normal_cutoff](/index.php/Surface_carve_normal_cutoff "Surface carve normal cutoff")
  * [surface_normal](/index.php?title=Surface_normal&action=edit&redlink=1 "Surface normal \(page does not exist\)")
  * [surface_carve_selection](/index.php/Surface_carve_selection "Surface carve selection")
  * [surface_optimize_subsets](/index.php?title=Surface_optimize_subsets&action=edit&redlink=1 "Surface optimize subsets \(page does not exist\)")
  * [surface_carve_state](/index.php?title=Surface_carve_state&action=edit&redlink=1 "Surface carve state \(page does not exist\)")
  * [surface_poor](/index.php?title=Surface_poor&action=edit&redlink=1 "Surface poor \(page does not exist\)")
  * [surface_circumscribe](/index.php?title=Surface_circumscribe&action=edit&redlink=1 "Surface circumscribe \(page does not exist\)")
  * [surface_proximity](/index.php/Surface_proximity "Surface proximity")
  * [surface_clear_cutoff](/index.php?title=Surface_clear_cutoff&action=edit&redlink=1 "Surface clear cutoff \(page does not exist\)")
  * [surface_quality](/index.php/Surface_quality "Surface quality")
  * [surface_clear_selection](/index.php?title=Surface_clear_selection&action=edit&redlink=1 "Surface clear selection \(page does not exist\)")
  * [surface_ramp_above_mode](/index.php/Surface_ramp_above_mode "Surface ramp above mode")
  * [surface_clear_state](/index.php?title=Surface_clear_state&action=edit&redlink=1 "Surface clear state \(page does not exist\)")
  * [surface_solvent](/index.php/Surface_solvent "Surface solvent")
  * [surface_color](/index.php/Surface_color "Surface color")
  * [surface_trim_cutoff](/index.php?title=Surface_trim_cutoff&action=edit&redlink=1 "Surface trim cutoff \(page does not exist\)")
  * [surface_debug](/index.php?title=Surface_debug&action=edit&redlink=1 "Surface debug \(page does not exist\)")
  * [surface_trim_factor](/index.php?title=Surface_trim_factor&action=edit&redlink=1 "Surface trim factor \(page does not exist\)")
  * [surface_miserable](/index.php?title=Surface_miserable&action=edit&redlink=1 "Surface miserable \(page does not exist\)")
  * [surface_type](/index.php/Surface_type "Surface type")
  * [surface_mode](/index.php/Surface_mode "Surface mode")



### Examples

#### Transparency

To adjust the transparency of surfaces try: 
    
    
    set transparency, 0.5
    

Where 1.0 will be an invisible and 0.0 a completely solid surface. 

#### Quality

To smooth your surface representation try: 
    
    
    set surface_quality, 1
    

or higher if you wish, though it will take longer and might look odd. 

#### Probe Radius

To change the probe radius other than default 1.4 Å, you need to change the solvent radius, say, 1.6 Å: 
    
    
    set solvent_radius, 1.6
    

If the surface does not change correspondingly, use: 
    
    
    rebuild
    

## Tips

### Exporting Surface/Mesh Coordinates to File

PyMOL can export its coordinates as WRL wireframe model files for VRML input. 

#### Older PyMOL Versions
    
    
    # export the coordinates to povray
    open("surface.inc","w").write(cmd.get_povray()[1])
    

#### Newer PyMOL Versions
    
    
    # export the coordinates to .wrl file
    save myscene.wrl
    

or 
    
    
    # export the coordinates to .obj file. Only surface representation can be saved as .obj.
    # NOTE: the coordinates are saved in camera coordinate system.
    save myscene.obj
    

### Representation-independent Color Control

To color the surface representation a different color than the underlying cartoon or ligand representations, simply duplicate the object, show only the surface in the duplicate, and show only the cartoon and/or ligands in the original object. 

Or use the [surface_color](/index.php/Surface_color "Surface color") setting that is available. 

  


[![](/images/5/51/Representation_independent_color_control.jpg)](/index.php/File:Representation_independent_color_control.jpg)

[](/index.php/File:Representation_independent_color_control.jpg "Enlarge")

Representation-independent Color Control Example

### Displaying a protein as surface with a ligand as sticks

An easy way to do this is to create separate objects for each type of display. 

1 Load your protein 

2 Select the ligand 

3 Create a separate object for the ligand 

4 Remove ligand atoms from the protein 

5 Display both objects separately 

Example: 
    
    
    load prot.ent,protein
    select ligand,resn FAD
    create lig_sticks,ligand
    remove ligand
    show sticks,lig_sticks
    show surface,protein
    

  
Even easier is to: 

1 Load the protein 

2 S (Show) > organic > stick 

3 S (Show) > surface 

### Calculating a partial surface

There is, until now, an undocumented way to calculate a surface for only a part of an object without creating a new one: 
    
    
    flag ignore, not A/49-63/, set
    delete indicate
    show surface
    

If the surface was already computed, then you'll also need to issue the command: 
    
    
    rebuild
    

See [Get_Area](/index.php/Get_Area "Get Area") for more information on surface area calculations. 

### Displaying surface inside a molecule

As far as I can tell, setting ambient to zero alone doesn't quite do the job, since some triangles still get lit by the light source. The best combination I can find is: 
    
    
    set ambient=0
    set direct=0.7
    set reflect=0.0
    set backface_cull=0
    

Which gives no shadows and only a few artifacts. 

As an alternative, you might just consider showing the inside of the surface directly...that will create less visual artifacts, and so long as ambient and direct are sufficiently low, it will look reasonable in "ray". 
    
    
    util.ray_shadows("heavy")
    set two_sided_lighting=1
    set backface_cull=0
    

### Creating a Closed Surface

  * [![Example OPEN Surface](/images/6/6a/Surface_open.png)](/index.php/File:Surface_open.png "Example OPEN Surface")

Example OPEN Surface 

  * [![Example CLOSED Surface](/images/f/f5/Surface_closed.png)](/index.php/File:Surface_closed.png "Example CLOSED Surface")

Example CLOSED Surface 




To create what I'll call a **closed surface** (see images), you need to first make your atom selections, then create a new object for that selection then show the surface for that object. Here's an example. 
    
    
     sel A, id 1-100
     create B, A
     show surface, B
    

### Smooth surface quick (blob)

To get a quick blob type surface (not as accurate): 
    
    
    set solvent_radius, 4   
    alter all, vdw=4 
    sort
    set surface_quality, 1
    

### Smooth surface accurate (blob)

To get an accurate blob type surface: 
    
    
    set surface_quality, 1
    alter all, b=50
    alter all, q=1
    set gaussian_resolution,5
    map_new mapA, gaussian, 1, sele or pdb, 6
    isosurface surfA, mapA
    

**Notes:** Set gaussian resolution is variable with a larger number causing a more smooth surface (4 is medium and 8 is very smooth). The temperature factor field (b) has at least as much impact as the resolution setting, so increasing b factors is the more computationally efficient way of increasing the blur effect. If you are displaying more then one surface in a .pse file you must create a new map for each one (if you have three you will create mapA for the first, mapB for the second and mapC for the third), then you apply an isosurface to each map (isosurface surfA, mapA - isosurface surfB, mapB - isosurface surfC, mapC). 

### Huge Surfaces

If your protein or complex is too large to render ([ray](/index.php/Ray "Ray") runs out of RAM, for example) then check out the [tip for huge surfaces](/index.php/Huge_surfaces "Huge surfaces"). 

## Performance

To optimize performance and responsiveness, PyMOL tends to defer compute-intensive tasks until their results are actually needed. Thus, 
    
    
    cmd.show("surface")
    

doesn't actually show a surface, it only sets the surface visibility flag on the atoms present (for future reference). An actual surface won't be computed until PyMOL is asked to refresh or render the display. When running a script, you can force an update by calling: 
    
    
    cmd.refresh()
    

after cmd.show. 

Retrieved from "[https://pymolwiki.org/index.php?title=Surface&oldid=12463](https://pymolwiki.org/index.php?title=Surface&oldid=12463)"


---

## Surface type

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



# Overview

**Surface_type** sets the type (dots, filled triangles, empty triangles) of the surface that PyMOL will render. 

# Syntax
    
    
    set surface_type, typeNumber  # typeNumber = [1..3]
    

# Examples
    
    
    set surface_type, 2
    

  * [![surface_type, 0 # usual](/images/1/18/St0.png)](/index.php/File:St0.png "surface_type, 0 # usual")

surface_type, 0 # usual 

  * [![surface_type, 1 # dots](/images/2/2f/St1.png)](/index.php/File:St1.png "surface_type, 1 # dots")

surface_type, 1 # dots 

  * [![surface_type, 2 # triangles](/images/1/19/St2.png)](/index.php/File:St2.png "surface_type, 2 # triangles")

surface_type, 2 # triangles 

  * [![surface_type, 3 # is this the same as surface_type 0?](/images/f/f2/St3.png)](/index.php/File:St3.png "surface_type, 3 # is this the same as surface_type 0?")

surface_type, 3 # is this the same as surface_type 0? 




# See Also

[surface](/index.php/Surface "Surface"), [as](/index.php/As "As")

Retrieved from "[https://pymolwiki.org/index.php?title=Surface_type&oldid=5049](https://pymolwiki.org/index.php?title=Surface_type&oldid=5049)"


---

## Volume

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/b/b1/1okyVol.png)](/index.php/File:1okyVol.png)

[](/index.php/File:1okyVol.png "Enlarge")

Volume visualization of electron density for PDB 1oky

[![](/images/5/5d/1okyVolPanel.png)](/index.php/File:1okyVolPanel.png)

[](/index.php/File:1okyVolPanel.png "Enlarge")

Volume panel for the 1oky volume example. It has the iso-levels on the x-axis and the opacity (1.0 - transparency) on the y-axis.

  
Volume creates a new volume object from a map object. The data (3D scalar fields) are shown as a true 3D object using coloring and transparencies defined by the user to illustrate the data values. This technique supports single and multiple isosurfaces. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Screencasts
  * 5 Changes with PyMOL Version
  * 6 Ray Tracing
  * 7 Known Limitations
  * 8 See Also



## Usage
    
    
    volume name, map [, ramp [, selection [, buffer [, state [, carve ]]]]]
    

## Arguments

  * name = the name for the new volume object.
  * map = the name of the map object to use for computing the volume.
  * ramp = str: named color ramp {default: }
  * selection = an atom selection about which to display the mesh with an additional "buffer" (if provided).
  * carve = a radius about each atom in the selection for which to include density. If "carve" is not provided, then the whole brick is displayed.



## Example
    
    
    fetch 1oky, type=2fofc, async=0
    volume 1okyVol, 1oky_2fofc
    

## Screencasts

  * [Silent demo movie](http://www.youtube.com/watch?v=tuAo_8-_HIc) showing the basics of loading and using a volume in PyMOL. There are more capabilities, but this is the basic functionality.



## Changes with PyMOL Version

  * 1.4.0: first version with volume support
  * 1.7.2: 
    * pre-integrated volume rendering (volume_mode=1) as Incentive-PyMOL-only feature.
    * scripting support with custom color ramp ([volume_color](/index.php/Volume_color "Volume color"), [volume_ramp_new](/index.php?title=Volume_ramp_new&action=edit&redlink=1 "Volume ramp new \(page does not exist\)"))
    * improved volume panel, panel can be opened from the object menu ("C > panel")
    * lots of bugs fixed



## Ray Tracing

**There is no actual ray tracing support**. The volume rendering is implemented exclusively with OpenGL shaders. The recommended way to render a high resolution image is to use the [draw](/index.php/Draw "Draw") command. Example: 
    
    
    # render high resolution image on screen
    draw 4000, 3000, antialias=2
    png highres.png
    

Ray trace specific features like shadows or outlines ([ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")) are not available with this approach. If such features are needed, the [ray_volume](/index.php?title=Ray_volume&action=edit&redlink=1 "Ray volume \(page does not exist\)") setting activates a hybrid solution, which blends the OpenGL rendered volume with the ray traced non-volume objects. Image sizes other than the current window size are not possible. 
    
    
    # compose on-screen volume with ray traced image
    set ray_volume
    ray
    png composed.png
    

Neither of these two solutions work with [headless (batch) mode](/index.php/Launching_PyMOL#Running_PyMOL_in_batch_mode "Launching PyMOL"). 

## Known Limitations

  * No real ray-tracing support yet
  * Multiple volume objects don't blend properly



## See Also

  * <http://pymol.org/volume>
  * <http://pymol.org/d/media:volumevisualization>
  * [volume_color](/index.php/Volume_color "Volume color")
  * [volume_ramp_new](/index.php?title=Volume_ramp_new&action=edit&redlink=1 "Volume ramp new \(page does not exist\)")
  * [map_new](/index.php/Map_new "Map new")
  * [isomesh](/index.php/Isomesh "Isomesh")
  * [isosurface](/index.php/Isosurface "Isosurface")



Retrieved from "[https://pymolwiki.org/index.php?title=Volume&oldid=12536](https://pymolwiki.org/index.php?title=Volume&oldid=12536)"


---

