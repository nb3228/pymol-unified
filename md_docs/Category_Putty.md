# Category: Putty

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

## Cartoon putty scale max

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

cartoon_putty_scale_max sets the max scale that PyMOL will use for rendering cartoon putty. 

# Usage
    
    
    # set a max scale
    set cartoon_putty_scale_max, 3
    
    # no limits on scaling
    set cartoon_putty_scale_max, -1
    

  


# See Also

  * [putty](/index.php/Cartoon#Sausage_Representation "Cartoon")
  * [Cartoon_putty_scale_min](/index.php/Cartoon_putty_scale_min "Cartoon putty scale min")
  * [cartoon_putty_transform](/index.php/Cartoon_putty_transform "Cartoon putty transform")



Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_putty_scale_max&oldid=8241](https://pymolwiki.org/index.php?title=Cartoon_putty_scale_max&oldid=8241)"


---

## Cartoon putty scale min

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

cartoon_putty_scale_min sets the minimum scale that PyMOL will use for rendering cartoon putty. 

# Usage
    
    
    # set a minimum scale
    set cartoon_putty_scale_min, 3
    
    # no limits on scaling
    set cartoon_putty_scale_min, -1
    

  


# See Also

  * [putty](/index.php/Cartoon#Sausage_Representation "Cartoon")
  * [Cartoon_putty_scale_max](/index.php/Cartoon_putty_scale_max "Cartoon putty scale max")
  * [cartoon_putty_transform](/index.php/Cartoon_putty_transform "Cartoon putty transform")



Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_putty_scale_min&oldid=8240](https://pymolwiki.org/index.php?title=Cartoon_putty_scale_min&oldid=8240)"


---

## Cartoon putty transform

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting determines how PyMOL renders the transformation of original values into putty related settings. 

  


# Usage
    
    
    # normalized nonlinear scaling
    set cartoon_putty_transform, 0
    # relative nonlinear scaling
    set cartoon_putty_transform, 1
    # scaled nonlinear scaling
    set cartoon_putty_transform, 2
    # absolute nonlinear scaling
    set cartoon_putty_transform, 3
    
    # normalized linear scaling
    set cartoon_putty_transform, 4
    # relative linear scaling
    set cartoon_putty_transform, 5
    # scaled linear scaling
    set cartoon_putty_transform, 6
    # absolute linear scaling from the B factor
    set cartoon_putty_transform, 7
    
    # implied RMS scaling
    set cartoon_putty_transform, 8
    

  


# See Also

  * [putty](/index.php/Cartoon#Sausage_Representation "Cartoon")
  * [Cartoon_putty_scale_max](/index.php/Cartoon_putty_scale_max "Cartoon putty scale max")
  * [Cartoon_putty_scale_min](/index.php/Cartoon_putty_scale_min "Cartoon putty scale min")



Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_putty_transform&oldid=8239](https://pymolwiki.org/index.php?title=Cartoon_putty_transform&oldid=8239)"


---

