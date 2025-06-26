# Category: Nucleic Acids

## 3DNA

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_**Update:** See also the [dssr_block](/index.php/Dssr_block "Dssr block") wrapper script for **x3dna-dssr**._

3DNA[[1]](http://x3dna.org/) provides a ruby script (**blocview**) which produces Calladine-Drew style block representations based on the standard reference frame for nucleic acids. The produced files are in raster3D format, and are also processed by molscript. Therefore to produce the following representations you will need to install: 

\- 3DNA [[2]](http://x3dna.org/)

\- Raster 3D [[3]](http://skuld.bmsc.washington.edu/raster3d/raster3d.html)

\- Molscript [[4]](http://www.avatar.se/molscript/download.html)

Once you have the previous software installed, and the path to their binaries correctly configured, you can get the block representation for any nucleic acid. The next example shows how to do it for tRNA (PDB:ID 1ehz) 

## Contents

  * 1 Invoking
  * 2 Example 1
  * 3 Syntax
  * 4 Example 2
  * 5 More



## Invoking
    
    
     
    pdb_get.py 1ehz
    blocview -o 1ehz.pdb
    

The pdb_get python script comes from [[5]](http://pldserver1.biochem.queensu.ca/~rlc/work/scripts/pdb_get.py)Dr. Robert Campbell's website. 

Once blocview is run it will generate various r3d files which can be combined with the original pdb files to produce the following image: 

## Example 1

  * [![Using tb.pdb and t2.r3d](/images/d/dd/Trna.png)](/index.php/File:Trna.png "Using tb.pdb and t2.r3d")

Using tb.pdb and t2.r3d 




To generate the previous image the following script was used: trna.pml 

## Syntax
    
    
    load tb.pdb
    load t2.r3d
    bg_color white
    hide lines
    zoom *,-5
    set ray_trace_mode, 3
    set ray_trace_fog, 0
    set ray_shadows, 0
    set orthoscopic, 1
    set antialias, 5
    set valence, 1
    util.cba(29)
    color grey, (elem C)
    cartoon arrow
    set cartoon_ladder_mode, 0
    set cartoon_rect_width, 0.2
    set cartoon_rect_length, 0.5
    show cartoon
    set stick_radius, 0.14
    show sticks
    set_view (\
         0.680474579,   -0.153203458,   -0.716576934,\
         0.658882320,   -0.300013036,    0.689829707,\
        -0.320666909,   -0.941552401,   -0.103208199,\
        -0.000084338,    0.000022471, -263.030426025,\
        57.723434448,   45.338260651,   20.895099640,\
       218.710235596,  307.348541260,    1.000000000 )
    set cartoon_color, green, resn G
    set cartoon_color, yellow, resn C
    set cartoon_color, red, resn A
    set cartoon_color, cyan, resn U
    ray 1024,768
    png trna.png
    quit
    

## Example 2

An automatically generated pymol ray traced image can also be obtained by running **blocview** from 3DNA v. 2.0 directly. Using the following command: 
    
    
    blocview -o -t=100 1ehz.pdb
    

And the result is: 

  * [![X3dna r3d pymol.png](/images/9/97/X3dna_r3d_pymol.png)](/index.php/File:X3dna_r3d_pymol.png)




  


## More

For more examples of figures obtained using 3DNA and pymol follow the next link: [[[6]](http://mesguerra.org/render/render.htm)] 

Retrieved from "[https://pymolwiki.org/index.php?title=3DNA&oldid=12273](https://pymolwiki.org/index.php?title=3DNA&oldid=12273)"


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

## Cartoon nucleic acid color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



## Overview

This setting determines the color of the sugar bb, when in nucleic acid mode. 

## Syntax
    
    
    # set the sugar bb to red
    set cartoon_nucleic_acid_color, red
    

## Examples

  * [![Color Red](/images/6/69/Cna_color2.png)](/index.php/File:Cna_color2.png "Color Red")

Color Red 

  * [![Color Blue](/images/e/e1/Cna_color1.png)](/index.php/File:Cna_color1.png "Color Blue")

Color Blue 




## See Also

[Cartoon#Nucleic_Acid_Representation](/index.php/Cartoon#Nucleic_Acid_Representation "Cartoon")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_nucleic_acid_color&oldid=5252](https://pymolwiki.org/index.php?title=Cartoon_nucleic_acid_color&oldid=5252)"


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

## Cartoon ring mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



## Overview

PyMOL provides a variety of representations for illustrating nucleic acid structures. This setting, **ring_mode** , affects how PyMol draws the bases in a nucleic acid representation. 

## Syntax
    
    
    # try other natural numbers
    set cartoon_ring_mode, 3
    

Currently, there are 7 modes, numbered **0** through **6**. See the illustrations below. 

## Examples

  * [![Set to 0](/images/0/06/Crm_0.png)](/index.php/File:Crm_0.png "Set to 0")

Set to 0 

  * [![Set to 1](/images/a/ac/Crm_1.png)](/index.php/File:Crm_1.png "Set to 1")

Set to 1 

  * [![Set to 2](/images/e/e0/Crm_2.png)](/index.php/File:Crm_2.png "Set to 2")

Set to 2 

  * [![Set to 3](/images/c/ce/Crm_3.png)](/index.php/File:Crm_3.png "Set to 3")

Set to 3 

  * [![Set to 4](/images/6/65/Crm_4.png)](/index.php/File:Crm_4.png "Set to 4")

Set to 4 

  * [![Set to 6](/images/4/47/Crm_6.png)](/index.php/File:Crm_6.png "Set to 6")

Set to 6 




## See Also

[Cartoon#Nucleic_Acid_Representation](/index.php/Cartoon#Nucleic_Acid_Representation "Cartoon")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_ring_mode&oldid=4294](https://pymolwiki.org/index.php?title=Cartoon_ring_mode&oldid=4294)"


---

## Cartoon ring transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



## Overview

Sets the level of transparency in nucleic acid representation. A value of 0.0 is completely opaque; 1.0 is completely transparent. 

  


## Syntax
    
    
    set cartoon_ring_transparency, 0.5  # valid values reals in the range of 0.0 - 1.0
    

## Examples

  * [![Set to 1.0](/images/2/24/Crt_05.png)](/index.php/File:Crt_05.png "Set to 1.0")

Set to 1.0 

  * [![Set to 0.5](/images/6/65/Crt_1.png)](/index.php/File:Crt_1.png "Set to 0.5")

Set to 0.5 




## See Also

[Nucleic Acid Representations](/index.php/Cartoon#Nucleic_Acid_Representation "Cartoon") and the [Nucleic Acids Category](/index.php/Category:Nucleic_Acids "Category:Nucleic Acids")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_ring_transparency&oldid=6296](https://pymolwiki.org/index.php?title=Cartoon_ring_transparency&oldid=6296)"


---

## Caver

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/loschmidt/caver-pymol-plugin/archive/master.zip>  
Author(s)  | CAVER Team   
License  | [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)  
<http://www.caver.cz>  
  
  
The Caver3 plugin is the successor of [Caver2](/index.php/Caver2 "Caver2"). CAVER is a software tool for analysis and visualization of tunnels and channels in protein structures. 

## Installation

  1. Install [Java](https://java.com/en/download/)
  2. Install [caver-pymol-plugin-master.zip](https://github.com/loschmidt/caver-pymol-plugin/archive/master.zip) using the PyMOL [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



[![PluginManagerCaverPlugin301 fixed.png](/images/b/b2/PluginManagerCaverPlugin301_fixed.png)](/index.php/File:PluginManagerCaverPlugin301_fixed.png)

## See Also

  * [Mole2](/index.php/Mole2 "Mole2")
  * [Caver2](/index.php/Caver2 "Caver2") (old version)



Retrieved from "[https://pymolwiki.org/index.php?title=Caver3&oldid=13180](https://pymolwiki.org/index.php?title=Caver3&oldid=13180)"


---

## Displaying Biochemical Properties

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Selecting secondary structures
    * 1.1 Manually Assigning Secondary Structure
    * 1.2 See Also
  * 2 Coloring
    * 2.1 Color by atom type from a script
    * 2.2 Assign color by B-factor
    * 2.3 Representation-independent color control
  * 3 Bonds
    * 3.1 Displaying double bonds
    * 3.2 Hydrogen bonds and Polar Contacts
      * 3.2.1 Hydrogen bonds between specific atoms
      * 3.2.2 Hydrogen bonds where find->polar contacts doesn't do what you need
      * 3.2.3 Details
  * 4 Calculating dihedral angles
  * 5 Cavities
  * 6 Surface-Related
    * 6.1 Surface Area
    * 6.2 Polar surface area
    * 6.3 Display solvent accessible surface
    * 6.4 Contact Potential
    * 6.5 Electric Field Lines
    * 6.6 Residues with functional groups
  * 7 Backbones
    * 7.1 Displaying the C-Alpha trace of proteins
    * 7.2 Displaying the Amino Acid Backbone
    * 7.3 Displaying the Phosphate backbone of nucleic acids
      * 7.3.1 Native Nucleic Acid Rendering in PyMol
  * 8 Align proteins with CA fit



## Selecting secondary structures

A few examples: 
    
    
    select helix, (ss h)
    select sheet, (ss s)
    select loop, (ss l+'')
    

### Manually Assigning Secondary Structure

You can manually assign secondary stuctures to your protein by 
    
    
    alter 96-103/, ss='S'
    alter 96-103/, ss='H'
    alter 96-103/, ss='L'
    

to set residues 96-103 to beta Strand, alpha Helix, and Loop respectively. 

### See Also

[DSSP](/index.php/DSSP "DSSP") [Dss](/index.php/Dss "Dss") [Caver](/index.php/Caver "Caver")

[FAQ](/index.php/Category:FAQ "Category:FAQ") [Displaying Biochemical Properties](/index.php/Category:Objects_and_Selections "Category:Objects and Selections")

## Coloring

See also [Category:Coloring](/index.php/Category:Coloring "Category:Coloring"). 

### Color by atom type from a script

See [Color](/index.php/Color "Color") for this. 

### Assign color by B-factor

See section [Color](/index.php/Color "Color") for this. 

### Representation-independent color control

See section [Surface#Representation-independent Color Control](/index.php/Surface#Representation-independent_Color_Control "Surface"). 

## Bonds

PyMOL can deduce bonds from the PDB structure file, even if the **CONECT** records are missing. In fact, PyMOL guesses bonding connectivity based on proximity, based on the empirical observation that two atoms of a given radius will not be generally closer than a certain distance unless they are bonded. 

### Displaying double bonds

  * [![Image showing double bonds in PyMOL. Double bonds are supported in lines and sticks.](/images/d/d2/DoubleBonds.png)](/index.php/File:DoubleBonds.png "Image showing double bonds in PyMOL. Double bonds are supported in lines and sticks.")

Image showing double bonds in PyMOL. Double bonds are supported in [lines](/index.php/Lines "Lines") and [sticks](/index.php/Sticks "Sticks"). 




You can go into the [lines](/index.php/Lines "Lines") mode and turning on the valence display: 
    
    
    hide
    as lines
    set valence, 0.1
    

A higher value for valence spreads things out more. See the [Valence](/index.php/Valence "Valence") page for more information on options for changing the appearance of the valence lines. 

### Hydrogen bonds and Polar Contacts

[![](/images/a/ad/Polar_contacts_small.png)](/index.php/File:Polar_contacts_small.png)

[](/index.php/File:Polar_contacts_small.png "Enlarge")

Polar Contacts in PyMol

Using the actions [A] button for an object or selection you can display Hydrogen bonds and Polar Contacts. [A]->find->polar contacts-><select from menu>

The command behind the menus is the **dist** ance command called with the additional argument mode=2. 

Parameters that control the the identification of H-bonds are defined as 
    
    
    set h_bond_cutoff_center, 3.6
    

with ideal geometry and 
    
    
    set h_bond_cutoff_edge, 3.2
    

with minimally acceptable geometry. 

These settings can be changed *before* running the detection process (dist command mode=2 or via the menus). 

Note that the hydrogen bond geometric criteria used in PyMOL was designed to emulate that used by [DSSP](http://swift.cmbi.kun.nl/gv/dssp/). 

#### Hydrogen bonds between specific atoms
    
    
    dist name, sele1, sele2, mode=2
    

  


#### Hydrogen bonds where find->polar contacts doesn't do what you need

You can show H-bonds between two objects using atom selections so long as hydrogens are present in both molecules. If you don't have hydrogens, you can use [h_add](/index.php/H_add "H add") on the proteins, or provide ligands with valence information and then use h_add. 

Two examples are below. For clarity, they draw dashes between the heavy atoms and hide the hydrogens. 
    
    
    # EXAMPLE 1: Show hydrogen bonds between protein 
    # and docked ligands (which must have hydrogens)
    
    load target.pdb,prot
    load docked_ligs.sdf,lig
    
    # add hydrogens to protein
    
    h_add prot
    
    select don, (elem n,o and (neighbor hydro))
    select acc, (elem o or (elem n and not (neighbor hydro)))
    dist HBA, (lig and acc),(prot and don), 3.2
    dist HBD, (lig and don),(prot and acc), 3.2
    delete don
    delete acc
    hide (hydro)
    
    hide labels,HBA
    hide labels,HBD
    
    
    
    # EXAMPLE 2
    # Show hydrogen bonds between two proteins
    
    load prot1.pdb
    load prot2.pdb
    
    h_add prot1
    h_add prot2
    
    select don, (elem n,o and (neighbor hydro))
    select acc, (elem o or (elem n and not (neighbor hydro)))
    dist HBA, (prot1 and acc),(prot2 and don), 3.2
    dist HBD, (prot1 and don),(prot2 and acc), 3.2
    delete don
    delete acc
    hide (hydro)
    
    hide labels,HBA
    hide labels,HBD
    
    # NOTE: that you could also use this approach between two
    # non-overlapping selections within a single object.
    

The "polar contacts" mentioned above are probably better at finding hydrogen bonds than these scripts. "Polar contacts" check geometry as well as distance. 

#### Details

Generally speaking, PyMOL does not have sufficient information to rigorously determine hydrogen bonds, since typical PDB file are ambiguous with respect to charge states, bonds, bond valences, and tautomers. As it stands, all of those things are guessed heuristically. Rigorously determining the location of lone pair electrons and proton coordinates from raw PDB files is a nontrival problem especially when arbitrary small molecule structures are present. In addition, PyMOL would also need to consider the implied coordinate error due to overall structure resolution and local temperature factors before rigorously asserting than any specific hydrogen bond does or does not exist. 

Furthermore, our hydrogen bond detection machinery was originally developed for purposes of mimicking Kabsch and Sander's DSSP secondary structure assignment algorithm (Biopolymers 22, 2577, 1983) which is based on a rather generous notion of hydrogen bonding (see Kabsch Figure 1). 

Although this approximate capability can be accessed via the distance command using mode=2, the criteria applied by our implementation may be based on heavy-atom coordinates (only) and does not necessarily correspond to anything rigorous or published. So the bottom line is that PyMOL merely offers up putative polar contacts and leaves it to the user to determine whether or not the interactions present are in fact hydrogen bonds, salt bridges, polar interactions, or merely artifacts of incorrect assignments (i.e. two carbonyls hydrogen bonding because they're being treated like hydroxyls). 

With respect to the h_bond_* settings, the angle in question for h_bond_cutoff_* and h_bond_max_angle is ADH, assuming H exists. If H does not exist, then PyMOL will guess a hypothetical coordinate which may not actually be valid (in plane, etc.). Tthe hydrogen must also lie within a cone of space with its origin on A (along B->A) and angular width h_bond_cone. Since h_bond_cone in 180 by default, the present behavior is to simply reject any hydrogen bond where the hydrogen would lie behind the plane defined by the acceptor atom (A) in relation to its bonded atom(s) B (if any). In other words, if B is only one atom (e.g. C=O vs. C-O-C), then by default, HAB cannot be less then 90 degrees. 

The two h_bond_power_* settings are merely fitting parameters which enable PyMOL to reproduce a curve shape reflecting Kabsch Figure 1. The endpoints of the effective cutoff curve is a function of the two h_bond_cutoff_* setting. 

## Calculating dihedral angles

The get_dihedral function requires four single-atom selections to work: 
    
    
    get_dihedral prot1///9/C, prot1///10/N, prot1///10/CA, prot1///10/C
    

## Cavities

See [Surfaces_and_Voids](/index.php/Surfaces_and_Voids "Surfaces and Voids"). Also [Caver](/index.php/Caver "Caver") and [CASTp](/index.php/CASTp "CASTp"). 

## Surface-Related

### Surface Area

To calculate the surface area of a selection, see [Get_Area](/index.php/Get_Area "Get Area"). 

### Polar surface area

For a solvent accessible PSA approximation: 
    
    
    set dot_density, 3
    remove hydro
    remove solvent
    show dots
    set dot_solvent, on
    get_area elem N+O
    get_area elem C+S
    get_area all
    

For molecular PSA approximation 
    
    
    set dot_density, 3
    remove hydro
    remove solvent
    set dot_solvent, off
    get_area elem N+O
    get_area elem C+S
    get_area all
    

Showing dots isn't mandatory, but it's a good idea to confirm that you're getting the value for the atom dot surface you think you're using. Please realize that the resulting numbers are only approximate, reflecting the sum of partial surface areas for all the dots you see. To increase accuracy, set dot_density to 4, but be prepared to wait... 

### Display solvent accessible surface

Using the surface display mode, PyMOL doesn't show the solvent accessible surface, rather it shows the solvent/protein contact surface. The solvent accessible surface area is usually defined as the surface traced out by the center of a water sphere, having a radius of about 1.4 angstroms, rolled over the protein atoms. The contact surface is the surface traced out by the vdw surfaces of the water atoms when in contact with the protein. 

PyMOL can show solvent accessible surfaces using the dot or sphere representations: 

for dots: 
    
    
    show dots
    set dot_mode,1
    set dot_density,3
    

for spheres: 
    
    
    alter all,vdw=vdw+1.4
    show spheres
    

Once the Van der Waals radii for the selection have been altered, the surface representation will also be "probe-inflated" to show a pseudo solvent accessible surface, as detailed above. 

for surfaces: 
    
    
    alter all,vdw=vdw+1.4
    show surface
    

[![](/images/d/d0/Solvent-accessible_surface.jpg)](/index.php/File:Solvent-accessible_surface.jpg)

[](/index.php/File:Solvent-accessible_surface.jpg "Enlarge")

Solvent-Accessible Surface Example

  


  


  


  


  


  


  


  


  
Note that to display both the molecular surface and the solvent-accessible surface, the object must be duplicated, as is done for [Surface#Representation-independent Color Control](/index.php/Surface#Representation-independent_Color_Control "Surface"). This also applies if the spheres representation is to be used to display "real" atoms. 

### Contact Potential

See [Protein_contact_potential](/index.php/Protein_contact_potential "Protein contact potential") and [APBS](/index.php/APBS "APBS"). 

[![Prot contact pot.png](/images/6/64/Prot_contact_pot.png)](/index.php/File:Prot_contact_pot.png)

[](/index.php/File:Prot_contact_pot.png "Enlarge")

### Electric Field Lines

[![](/images/6/66/Elines.png)](/index.php/File:Elines.png)

[](/index.php/File:Elines.png "Enlarge")

PyMOL and APBS used to show electronic field lines.

To produce an image with electric field lines, first run APBS. Then, input the following: 
    
    
    gradient my_grad, pymol-generated
    ramp_new my_grad_ramp, pymol-generated
    color my_grad_ramp, my_grad
    

### Residues with functional groups

[![](/images/c/c4/1igt_cys_lys_asp_glu_colored.png)](/index.php/File:1igt_cys_lys_asp_glu_colored.png)

[](/index.php/File:1igt_cys_lys_asp_glu_colored.png "Enlarge")

Whole residues colored (Cys: yellow, Lys: blue, Asp and Glu: red)

Poor man's solution: Display protein as surface, colorize all Lys (-NH2), Asp and Glu (-COOH) and Cys (-SH): 
    
    
    remove resn hoh    # remove water
    h_add              # add hydrogens
    
    as surface
    color grey90
    
    color slate, resn lys       # lysines in light blue
    color paleyellow, resn cys  # cysteines in light yellow
    color tv_red, (resn asp or(resn glu))  # aspartic and glutamic acid in light red
    

[![](/images/3/33/1igt_functional_groups_colored.png)](/index.php/File:1igt_functional_groups_colored.png)

[](/index.php/File:1igt_functional_groups_colored.png "Enlarge")

Only central atoms of functional groups colored (Cys: S, Lys: NH2, Asp and Glu: CO2)

Not-_so_ -poor-man's solution: In order to have the functional groups better localized, only the central atoms can be colored: 

  * the S atom of cystein,
  * the N and H atoms of the free amine of lysine (may be displayed with three H atoms at all three possible positions)
  * the C and two O atoms of free carboxylic groups in aspartic and glutamic acid



In this way, they are better visible through the surface compared to only one colored atom, both amines and carboxylic groups consist of three colored atoms each. 
    
    
    remove resn hoh    # remove water
    h_add              # add hydrogens
    
    as surface
    color grey90
    
    select sulf_cys, (resn cys and (elem S))      # get the sulfur atom of cystein residues
    color yellow, sulf_cys
    
    select nitro_lys, (resn lys and name NZ)              # get the nitrogens of free amines ("NZ" in PDB file)
    select hydro_lys, (elem H and (neighbor nitro_lys))   # get the neighboring H atoms 
    select amine_lys, (nitro_lys or hydro_lys)
    color tv_blue, amine_lys
    
    
    select oxy_asp, (resn asp and (name OD1 or name OD2))             # get the two oxygens of -COOH  ("OD1", "OD2")
    select carb_asp, (resn asp and (elem C and (neighbor oxy_asp)))   # get the connecting C atom
    select oxy_glu, (resn glu and (name OE1 or name OE2))             # oxygens "OE1" and "OE2" in PDB file
    select carb_glu, (resn glu and (elem c and (neighbor oxy_glu)))
    select carboxy, (carb_asp or oxy_asp or carb_glu or oxy_glu)
    color tv_red, carboxy
    

By displaying the protein as non-transparent surface, only the functional groups (colored atoms) at the surface are visible. The visualization of those groups can be pronounced by displaying the corresponding atoms as spheres, e.g. "as spheres, carboxy + amine_lys + sulf_cys", in this way it might become more clear how accessible they are. 

When displaying the protein as cartoon, the functional groups can be shown as spheres, and the whole residues cys, lys, asp and glu as sticks connected to the backbone, with the atoms of the functional groups again as spheres. However, then also the not accessible residues inside the protein are visible. 

## Backbones

### Displaying the C-Alpha trace of proteins
    
    
    hide
    show ribbon
    set ribbon_sampling,1
    

And if your model only contains CA atoms, you'll also need to issue: 
    
    
    set ribbon_trace,1
    

### Displaying the Amino Acid Backbone

The easiest way to see the backbone of the protein is to do 
    
    
    hide all
    show ribbon
    

If you don't like the ribbon representation, you can also do something like 
    
    
    hide all
    show sticks, name C+O+N+CA
    

You can replace **sticks** in the above by other representations like **spheres** or **lines**. 

### Displaying the Phosphate backbone of nucleic acids

#### Native Nucleic Acid Rendering in PyMol

PyMol now better supports viewing nucleic acid structure. [Nuccyl](/index.php/Nuccyl "Nuccyl") still seems to be the reigning champ for image quality, but see PyMol's native [Cartoon](/index.php/Cartoon "Cartoon") command. For more information on representing nucleic acids, please see the [Nucleic Acids](/index.php/Category:Nucleic_Acids "Category:Nucleic Acids") Category. 

[![Cnam 0.png](/images/0/0c/Cnam_0.png)](/index.php/File:Cnam_0.png)

[](/index.php/File:Cnam_0.png "Enlarge")

Should you ever want to show the phosphate trace of a nucleic acid molecule: 
    
    
    def p_trace(selection="(all)"):
        s = str(selection)
        cmd.hide('lines',"("+s+")")
        cmd.hide('spheres',"("+s+")")
        cmd.hide('sticks',"("+s+")")
        cmd.hide('ribbon',"("+s+")")
        cmd.show('cartoon',"("+s+")")
        cmd.set('cartoon_sampling',1,"("+s+")")
        cmd.set('cartoon_tube_radius',0.5,"("+s+")")
    cmd.extend('p_trace',p_trace)
    

and then: 
    
    
    p_trace (selection)
    

## Align proteins with CA fit

If two proteins have significant homology, you can use the [Align](/index.php/Align "Align") command: 
    
    
    align prot1////ca,prot2
    

which will perform a sequence alignment of prot1 against prot2, and then an optimizing fit using the CA positions. I'm not sure if the help text for align got into 0.82, but the next version will definitely have it. 

Retrieved from "[https://pymolwiki.org/index.php?title=Displaying_Biochemical_Properties&oldid=11686](https://pymolwiki.org/index.php?title=Displaying_Biochemical_Properties&oldid=11686)"


---

## Symmetry

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

There are many ways that symmetry can be important/useful/beautiful to look at in macromolecular structures. A typically application is the reconstruction of a symmetric oligomer from a few subunits. PyMOL has tools that can help with this type of analysis or depiction. 

  


## Expanding crystallographic symmetry

Structures determined by X-ray crystallography are typically deposited as files containing the coordinates for one asymmetric unit (ASU). Knowledge of the symmetry operators that describe how the ASUs are arranged relative to each other allows the arrangement of the crystal lattice to be recreated. PyMOL can read this symmetry information from the input coordinate file and recreate the neigbouring copies of the ASU using symmexp. 

PyMOL's built in symmetry expansion functionality is available as A->generate->symmetry mates for an object or as the [symexp](/index.php/Symexp "Symexp") command. 

The [SuperSym](/index.php/SuperSym "SuperSym") plugin has additional unit cell and symmetry axis tools. 

## Displaying symmetry axes

Often it is necessary to be able to find and draw symmetry axes. There are some contributed scripts that help do this including [Symmetry Axis](/index.php/Symmetry_Axis "Symmetry Axis"), for which you need to know the coordinates and direction of the axis you would like to draw, and [RotationAxis](/index.php/RotationAxis "RotationAxis"), which does an excellent job of displaying symmetry relationships between selections. 

Retrieved from "[https://pymolwiki.org/index.php?title=Symmetry&oldid=13341](https://pymolwiki.org/index.php?title=Symmetry&oldid=13341)"


---

## MAC Install

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes how to install PyMOL on Mac OS X. 

## Contents

  * 1 Incentive PyMOL
    * 1.1 Launching from Command Line
    * 1.2 X11 Hybrid
    * 1.3 Stereo on Second Monitor
  * 2 Open-Source PyMOL
    * 2.1 Package managers
    * 2.2 Install from Source
    * 2.3 Install APBS with Fink
    * 2.4 Stereo issues
  * 3 See Also



## Incentive PyMOL

[Schrödinger](http://www.schrodinger.com) provides pre-compiled PyMOL to paying sponsors. The bundle also includes ready-to-use [APBS](/index.php/APBS "APBS"), [RigiMOL](/index.php/Morph "Morph"), an MPEG encoder for movie export, and a small molecule energy minimization engine. 

Download: <https://pymol.org/>

Installation: Drag **PyMOL.app** on the **/Applications** shortcut. (In principle, you could drag it into any Finder window and run it from there, it doesn’t have to live in /Applications). 

Uninstallation: Move **/Applications/PyMOL.app** to Trash 

### Launching from Command Line

The unix executable resides at **/Applications/PyMOL.app/Contents/MacOS/PyMOL**

### X11 Hybrid

_Applies to PyMOL 1.x, not to PyMOL 2.x_

MacPyMOL can optionally run with the same two-window GUI which PyMOL uses on Windows and Linux. This GUI has some additional features, like the [Plugin](/index.php/Plugins "Plugins") menu and the [Builder](/index.php/Builder "Builder"). 

Requires [XQuartz](http://xquartz.macosforge.org/). 

There are two ways to launch the X11 interface: 

  1. Rename or copy/duplicate **/Applications/MacPyMOL.app** to **/Applications/MacPyMOLX11Hybrid.app** or to **/Applications/PyMOLX11Hybrid.app**
  2. Launch the unix executable with the **-m** flag: **/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL -m**



### Stereo on Second Monitor

The trick to getting MacPyMOL to work in stereo on the second monitor is to force it to initially open on that display by providing an appropriate "-X #" (and perhaps -Y #) option on launch. That way the OpenGL context will be created with stereo support. 
    
    
    ./MacPyMOL.app/Contents/MacOS/MacPyMOL -X -1000
    ./MacPyMOL.app/Contents/MacOS/MacPyMOL -X -1000 -Y 100
    

**Source:** [Warren DeLano; PyMOL Users Archive](https://sourceforge.net/p/pymol/mailman/message/11671952/)

## Open-Source PyMOL

### Package managers

Open-Source PyMOL is available [free of charge](https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE) and may be readily installed via the [Homebrew](http://brew.sh/) (recommended), [MacPorts](https://www.macports.org/), or [Fink](http://www.finkproject.org/) package managers. 
    
    
    # Homebrew (recommended)
    brew install brewsci/bio/pymol
    
    # Fink
    fink install pymol-py27
    
    # MacPorts
    sudo port install pymol
    

You may need to make sure that the dependencies are installed with the required flags, e.g. for MacPorts: 
    
    
    # MacPorts
    sudo port install tcl -corefoundation
    sudo port install tk -quartz
    

If PyMOL complains that it wasn't able to find X11, try starting xquartz first, then run pymol from the console. 

### Install from Source

If you want the latest PyMOL code (warning: might include experimental changes), then follow the [Linux installation instructions](/index.php/Linux_Install#Install_from_source "Linux Install"). You will need an environment like Fink, MacPorts or Homebrew to install the dependencies. Make sure you use the appropriate python interpreter (e.g. **/sw/bin/python2.7** when using Fink). 

To run PyMOL with a native PyQt library (linked against macOS OpenGL framework, not against XQuartz), it needs to be built with the `--osx-frameworks` option: 
    
    
    python setup.py --osx-frameworks install
    

### Install APBS with Fink

To use the electrostatics plugin, you will need [APBS](http://apbs.sourceforge.net/) and its dependencies. These are also available as Fink packages, and include [APBS](http://pdb.finkproject.org/pdb/package.php/apbs), [maloc](http://pdb.finkproject.org/pdb/package.php/maloc) and [pdb2pqr](http://pdb.finkproject.org/pdb/package.php/pdb2pqr). If you have multiple processors available, you might wish to install the [MPI version of APBS](http://pdb.finkproject.org/pdb/package.php/apbs-mpi-openmpi). 

Issuing the command 
    
    
    fink install apbs
    

will install apbs and its required dependencies for you. The fink pymol package is already preconfigured to do the right thing to use apbs as a plugin. 

### Stereo issues

Some older Macs seem to crash with stereo graphics. If this happens to you, a workaround is to launch PyMOL explicitly in Mono mode with `pymol -M`. You can also set up an alias in your ~/.profile: 
    
    
    alias pymol='pymol -M'
    

## See Also

  * [pymolrc](/index.php/Pymolrc "Pymolrc")
  * [Linux Install](/index.php/Linux_Install "Linux Install")
  * [Windows Install](/index.php/Windows_Install "Windows Install")
  * [FreeMOL installation for MPEG movie export](/index.php/MovieSchool_6#Exporting_your_Movie "MovieSchool 6")
  * [Bill Scott’s](/index.php/User:Wgscott "User:Wgscott") [MacOSX-specific .pymolrc file](/index.php/MacOSX-specific_.pymolrc_file "MacOSX-specific .pymolrc file") and his crystallographic software [wiki](http://xanana.ucsc.edu/~wgscott/xtal/wiki/index.php/Main_Page) and [website](http://chemistry.ucsc.edu/~wgscott/xtal/), including instructions on [how to install precompiled binary packages using fink](http://xanana.ucsc.edu/~wgscott/xtal/wiki/index.php/Getting_your_fink_installation_to_use_packages_that_I_have_pre-compiled).
  * [Launching_PyMOL#MacOS_X](/index.php/Launching_PyMOL#MacOS_X "Launching PyMOL")



Retrieved from "[https://pymolwiki.org/index.php?title=MAC_Install&oldid=12800](https://pymolwiki.org/index.php?title=MAC_Install&oldid=12800)"


---

## Nuccyl

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

See [Nuccyl Home](http://www.biosci.ki.se/groups/ljo/software/nuccyl.html)

## Overview

**nuccyl** is a [Perl](http://www.perl.org/) program that allows [PyMOL](http://www.pymol.org/), a powerful open-source molecular modeling system, to display atomic models of nucleic acids in a highly simplified representation. By depicting the nucleic acid backbone as a cartoon and the bases of nucleotides as cylinders, similarly to the excellent programs [Ribbons](http://sgce.cbse.uab.edu/ribbons/) and [Drawna](http://www-ibmc.u-strasbg.fr/upr9002/westhof/download.html), nuccyl and PyMOL allow to quickly grasp the overall folding of an RNA or DNA molecule 

[![](/images/6/63/Phe_trna_small.jpg)](/index.php/File:Phe_trna_small.jpg)

[](/index.php/File:Phe_trna_small.jpg "Enlarge")

Example Nuccyl Image

[![](/images/5/5a/45s_rna_small.jpg)](/index.php/File:45s_rna_small.jpg)

[](/index.php/File:45s_rna_small.jpg "Enlarge")

Example Nuccyl II

[![](/images/b/be/50s_small.jpg)](/index.php/File:50s_small.jpg)

[](/index.php/File:50s_small.jpg "Enlarge")

Example Nuccyl III

[![](/images/7/7a/P4p6_small.jpg)](/index.php/File:P4p6_small.jpg)

[](/index.php/File:P4p6_small.jpg "Enlarge")

Example Nuccyl IV

nuccyl produces base cylinder coordinates and command files to display them with PyMOL. The base pair analysis programs [RNAView](http://beta-ndb.rutgers.edu/services/download/index.html#rnaview) and [3DNA](http://rutchem.rutgers.edu/%7exiangjun/3DNA/) can also be optionally used to obtain starting input files for nuccyl. 

[RNAView](http://beta-ndb.rutgers.edu/services/download/index.html#rnaview) can be downloaded and installed as [described](http://ndbserver.rutgers.edu/services/help/rnaview-readme.html); a [publication](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12824344&dopt=Abstract) describing the program is also available. Information and download links for [3DNA](http://rutchem.rutgers.edu/%7exiangjun/3DNA/) may be found at the program's [home page](http://rutchem.rutgers.edu/%7exiangjun/3DNA/); a [paper](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12930962&dopt=Abstract) has also been published. 

* * *

The four example images Copyright © 2000-2005 Luca Jovine. Do **NOT** reuse them without contacting him or reading his release license for the images! 

Retrieved from "[https://pymolwiki.org/index.php?title=Nuccyl&oldid=5397](https://pymolwiki.org/index.php?title=Nuccyl&oldid=5397)"


---

## Cartoon nucleic acid mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Examples
  * 3 Syntax
  * 4 See Also



## Overview

This setting modifies how PyMol draws the sugar bb of nucleic acids. Some differences between the settings are subtle, however PyMol gives you the ability to get the representation you want. Click on each images for a closer look. 

## Examples

  * [![Set to 0](/images/0/0c/Cnam_0.png)](/index.php/File:Cnam_0.png "Set to 0")

Set to 0 

  * [![Set to 1](/images/7/7d/Cnam_1.png)](/index.php/File:Cnam_1.png "Set to 1")

Set to 1 

  * [![Set to 2](/images/7/73/Cnam_2.png)](/index.php/File:Cnam_2.png "Set to 2")

Set to 2 




## Syntax
    
    
    set cartoon_nucleic_acid_mode, 1  # try different natural numbers
    

## See Also

[Cartoon#Nucleic_Acid_Representation](/index.php/Cartoon#Nucleic_Acid_Representation "Cartoon")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_nucleic_acid_mode&oldid=4242](https://pymolwiki.org/index.php?title=Cartoon_nucleic_acid_mode&oldid=4242)"


---

## S2S

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# S2S

Efficient RNA sequence manipulations (like multiple alignments) need to be constrained by rules of RNA structure folding. The structural knowledge has increased dramatically in the last years with the accumulation of several large RNA structures like those of the bacterial ribosome subunits. However, no tool in the RNA community provides an easy way to link and integrate progress made at the sequence level using the available three-dimensional information. S2S proposes a framework in which an user can easily display, manipulate, and interconnect heterogeneous RNA data like multiple sequence alignments, secondary and tertiary structures. S2S has been implemented with the Java language and has been developed and tested under UNIX systems like Linux and MacOSX. (From: <http://bioinformatics.org/S2S/>) 

Check out the [S2S Homepage](http://bioinformatics.org/S2S/). 

Retrieved from "[https://pymolwiki.org/index.php?title=S2S&oldid=4394](https://pymolwiki.org/index.php?title=S2S&oldid=4394)"


---

