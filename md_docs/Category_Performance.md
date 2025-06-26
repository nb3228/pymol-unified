# Category: Performance

## Async builds

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

A _beta_ setting. This setting controls how PyMOL sends out the geometry tasks for certain commands. Turning this on, should [speed](/index.php/Category:Performance "Category:Performance") up PyMOL in general. 

## Syntax
    
    
    set async_builds, 1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Async_builds&oldid=5511](https://pymolwiki.org/index.php?title=Async_builds&oldid=5511)"


---

## Cache

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Cache manages storage of precomputed results, such as molecular surfaces. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Notes
  * 4 API



# Usage
    
    
    cache action [, scenes [, state ]]
    

where, the arguments are: 

**action**

    

    string: enable, disable, read_only, clear, or optimize

**scenes**

    

    string: a space-separated list of scene names (default: _)_

**state**

    

    integer: state index (default: -1)

# Examples
    
    
    cache enable
    cache optimize
    cache optimize, F1 F2 F5
    

# Notes

"cache optimize" will iterate through the list of scenes provided (or all defined scenes), compute any missing surfaces, and store them in the cache for later reuse. 

# API
    
    
    cmd.cache(string action, string scenes, int state, int quiet)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cache&oldid=9091](https://pymolwiki.org/index.php?title=Cache&oldid=9091)"


---

## Cartoon sampling

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting changes the number of segments that makes up a given length of cartoon. The default value is ten; larger values (more segments) provide a smoother curve while requiring greater amounts of memory and slowing [performance](/index.php/Category:Performance "Category:Performance"). 

## Settings
    
    
    set cartoon_sampling, 20   # increases the number of cartoon segments
    

## Examples

Open the images to actually see the details! 

  * [![cartoon sampling 5](/images/1/1c/Cartoon_05.png)](/index.php/File:Cartoon_05.png "cartoon sampling 5")

cartoon sampling 5 

  * [![cartoon sampling 10 \(default\)](/images/b/b0/Cartoon_10.png)](/index.php/File:Cartoon_10.png "cartoon sampling 10 \(default\)")

cartoon sampling 10 (default) 

  * [![cartoon sampling 20](/images/4/48/Cartoon_20.png)](/index.php/File:Cartoon_20.png "cartoon sampling 20")

cartoon sampling 20 

  * [![cartoon sampling 30](/images/e/ec/Cartoon_30.png)](/index.php/File:Cartoon_30.png "cartoon sampling 30")

cartoon sampling 30 




Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_sampling&oldid=5503](https://pymolwiki.org/index.php?title=Cartoon_sampling&oldid=5503)"


---

## Defer builds mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

The defer_builds_mode setting for improved [performance](/index.php/Category:Performance "Category:Performance") with long trajectories. This now makes it possible to work with files containing thousands of states, and to render impossibly long movies piecewise. This setting, as shown below, stops PyMOL caching the geometry of trajectory data in RAM. 

# Usage
    
    
    # improve PyMOL performance for many-state objects and long movies.
    set defer_builds_mode, 3
    

# See Also

  * [Load_Traj](/index.php/Load_Traj "Load Traj")
  * [async_builds](/index.php/Async_builds "Async builds")



Retrieved from "[https://pymolwiki.org/index.php?title=Defer_builds_mode&oldid=5801](https://pymolwiki.org/index.php?title=Defer_builds_mode&oldid=5801)"


---

## Editing atoms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Altering atom types/elements
  * 2 Altering atom coordinates
  * 3 Translate or rotate individual objects
  * 4 Get coordinates from Python



## Altering atom types/elements

Example: 

To transform a carbon into an oxygen, pick an atom by a _select_ command or CTRL+middle. 
    
    
    alter pk1,elem='O'
    alter pk1,name='O2'
    

Note that the _name_ field should contain something more descriptive than just the element symbol. For reasons unknown (?), util.cbag() and similar methods will not recognize the new element. 

## Altering atom coordinates

Example: 
    
    
    alter_state 1,(pdb1cse),x=x-10.0
    

The latter section can contain formulae involving at least the xyz coordinates, lots of constants and the (+-*/) operators. 

  


## Translate or rotate individual objects

There is a "translate" function similar to "rotate", the docs for these don't exist yet, because the implementation isn't finished. However, feel free to use them in the following forms: 
    
    
    translate vector,object-name,state
       vector needs to be something like [x,y,z]
    
       translate [1,0,0],pept
    
    
    
    rotate axis,angle,object-name,state
       axis can be either the letter x,y,z or a 3D vector [x,y,z]
    
       rotate x,90,pept
       rotate [1,1,1],10,pept
    

## Get coordinates from Python

  1. The actual C-langauge arrays aren't exposed, but there are at least three different ways you can modify coordinates from within Python: You can get a python object which contains the molecular information, modify the coordinates in that object, load the modified molecule into PyMOL, update the modified coordinates to the original model, and then delete the modified object. (link to example)
  2. Another approach is the "alter_state" function, which can perform the same transformation in a single PyMOL command statement:
         
         alter_state 1,pept,(x,y)=(-y,x)
         

Likewise sub-selections can be transformed as well:
         
         alter_state 1,(pept and name ca),(x,y,z)=(x+5,y,z)
         

  3. A third approach is to use alter_state with the global "stored" object: Example in a Python script)



Approaches 2 gives the best [performance](/index.php/Category:Performance "Category:Performance"), approach 3 gives more flexibility, and approach 1 gives you a reusable and fully modifiable Python object. 

Retrieved from "[https://pymolwiki.org/index.php?title=Editing_atoms&oldid=10793](https://pymolwiki.org/index.php?title=Editing_atoms&oldid=10793)"


---

## Hash max

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Example
  * 4 See Also



## Overview

hash_max sets how much memory PyMOL uses when ray tracing. Higher values will enable PyMOL to ray trace more quickly, provided you can secure the necessary memory. Thus for large scenes, it can be useful to increase this value. Simpler scenes probably don't need it. If **hash_max** is set too high, then PyMOL can (and will) crash when it attempts to use more memory than available. Likewise, **hash_max** can be used to limit the use of memory to avoid crashes at the expense of increased computing time. 

## Syntax
    
    
    set hash_max, int
    

where **int** is a positive integer. The default value is 100. 

## Example
    
    
    set hash_max, 200
    

## See Also

[ray](/index.php/Ray "Ray"), [png](/index.php/Png "Png"), [draw](/index.php/Draw "Draw")

Retrieved from "[https://pymolwiki.org/index.php?title=Hash_max&oldid=11525](https://pymolwiki.org/index.php?title=Hash_max&oldid=11525)"


---

## Huge surfaces

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL can render **very large** (huge) surfaces of proteins, using a few tricks. First off, you should know that because of the size, you do not get the level of accuracy that you would for a smaller molecule. However, if you need to represent a molecule (or object) with a few million atoms, this may be your option of choice. 

# Example

The output from the example below is shown below. The PDB 1AON has nearly 60,000 atoms. This isn't the largest PDB, but it's a good example. 

  * [![Large surface.png](/images/f/f5/Large_surface.png)](/index.php/File:Large_surface.png)

  * [![Large surface approx.png](/images/7/7d/Large_surface_approx.png)](/index.php/File:Large_surface_approx.png)




## Strategy 1: Low resolution Gaussian Map
    
    
    # load a whopping big PDB
    fetch 1aon, struct, async=0
    
    # create a color spectrum over the object
    spectrum count, selection=struct
    
    # === now create a pseudo-fcalc map (a 3D volumetric scalar field) ===
    # set the B-factors nice and high for smoothness
    set gaussian_b_floor, 40
    
    # ~10 A map resolution
    set gaussian_resolution, 10
    
    # ~10 A map spacing with a 10 A surrounding buffer
    #  (you may need to vary this)
    map_new map, gaussian, 10, struct, 10
    
    # create a surface from the map
    isosurface surf, map, 1.0
    
    # now color the map based on the underlying protein
    ramp_new ramp, struct, [0,10], [atomic, atomic]
    
    color ramp, surf
    
    # hide the ramp
    disable ramp
    

## Strategy 2: CA-only model with oversized radii
    
    
    # load a whopping big PDB
    fetch 1aon, struct, async=0
    
    # reduce to CA atoms only
    remove struct & not guide
    
    # create a color spectrum over the object
    spectrum count, selection=struct
    
    # increase VDW radii to compensate volume of missing non-CA atoms
    alter struct, vdw=4
    
    # increase VDW radius of solvent
    set solvent_radius, 4, struct
    
    # show molecular surface
    as surface
    

Retrieved from "[https://pymolwiki.org/index.php?title=Huge_surfaces&oldid=12740](https://pymolwiki.org/index.php?title=Huge_surfaces&oldid=12740)"


---

## Load traj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**load_traj** loads a trajectory as "states" into an already loaded molecular object. 

Since version 1.0, PyMOL uses the [Molfile Plugin](http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/) backend, which supports a variety of trajectory file formats. Older versions only supported the ascii AMBER format (".trj" file extension). 

Loading a large trajectory may take up a lot of RAM, unless the [defer_builds_mode](/index.php/Defer_builds_mode "Defer builds mode") is set to 3. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 Notes
  * 5 See Also



## Usage
    
    
    load_traj filename [,object [,state [,format [,interval [,average ]
                       [,start [,stop [,max [,selection [,image [,shift 
                       [, plugin ]
                       ]]]]]]]]]
    

## Arguments

  * **filename** = str: trajectory file path
  * **object** = str: name of the molecular object where the trajectory should be appended as states {default: guess from filename}
  * **state** = int: object state where to start appending states. To discard the currently loaded coordinates, use _state=1_. To append new states, use _state=0_ {default: 0}
  * **format** = str: specify file type instead of guessing from file extension (only affects AMBER .trj format, use "plugin" argument for Molfile Plugin types) {default: }



## Examples
    
    
    # topology from PDB file, trajectory from DCD file
    load      sampletrajectory.pdb
    load_traj sampletrajectory.dcd
    
    # gromacs trajectory, using "mytraj" as object name
    load      sampletrajectory.gro, mytraj
    load_traj sampletrajectory.xtc, mytraj
    
    # desmond trajectory
    load      sample-out.cms, mytraj
    load_traj sample_trj/clickme.dtr, mytraj
    
    # playing through states, memory optimized (but eventually slower)
    set defer_builds_mode, 3
    mplay
    

## Notes

  * PyMOL does not know how to wrap the truncated octahedron used by Amber You will need to use the [cpptraj](http://ambermd.org/tutorials/analysis/tutorial0/index.htm) program first to do this.
  * The average option is not a running average. To perform this type of average, use the [smooth](/index.php/Smooth "Smooth") command after loading the trajectory file.
  * For quickly viewing Trajectories as a movie, use the [mset](/index.php/Mset "Mset") command to map each state to a movie frame.



useful notes from the email list:   
<http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg04272.html>   
<http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10266.html>

in one line convert dcd and psf to pdb : 

catdcd -o all.pdb -otype pdb -s autopsf.psf -stype psf out.dcd 

## See Also

[Load](/index.php/Load "Load"), [defer_builds_mode](/index.php/Defer_builds_mode "Defer builds mode")

Retrieved from "[https://pymolwiki.org/index.php?title=Load_traj&oldid=12888](https://pymolwiki.org/index.php?title=Load_traj&oldid=12888)"


---

## Mesh quality

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Altering the **mesh_quality** setting adjusts how precisely the mesh represents the surface. The higher the value the more curved and precise the surface is defined. The default is good compromise for quality and [speed](/index.php/Category:Performance "Category:Performance") (mesh calculation time). 

## Syntax
    
    
    # The default is 2.
    set mesh_quality, <integer>
    

## Example

  * [![mesh_quality 0 \(not smooth\)](/images/8/8e/Mesh_quality_0.jpg)](/index.php/File:Mesh_quality_0.jpg "mesh_quality 0 \(not smooth\)")

mesh_quality 0 (not smooth) 

  * [![mesh_quality 1](/images/3/35/Mesh_quality_1.jpg)](/index.php/File:Mesh_quality_1.jpg "mesh_quality 1")

mesh_quality 1 

  * [![mesh_quality 2 \(default\)](/images/6/65/Min_mesh_spacing_default.jpg)](/index.php/File:Min_mesh_spacing_default.jpg "mesh_quality 2 \(default\)")

mesh_quality 2 (default) 




Retrieved from "[https://pymolwiki.org/index.php?title=Mesh_quality&oldid=5510](https://pymolwiki.org/index.php?title=Mesh_quality&oldid=5510)"


---

## Orthoscopic

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
    * 1.1 Performance Notes
  * 2 Examples
  * 3 Syntax
  * 4 See Also



## Overview

**orthoscopic** turns on and off the perspective handling. Turning **on** orthoscopic rendering removes the effects of perspective. To adjust perspective rendering see, the [Field of View](/index.php/Field_Of_View "Field Of View") command. 

### [Performance](/index.php/Category:Performance "Category:Performance") Notes

Turning **on** orthoscopic rendering improves the [rendering speed](/index.php/Ray "Ray") by a factor of 2x-4x! So, if you have large images, or are rendering a movie that does not need perspective in the images, try turning off perspective by setting orthoscopic to **1**. 

## Examples

  * [![orthoscopic on \(perspective off\)](/images/e/ea/Orthoon2.png)](/index.php/File:Orthoon2.png "orthoscopic on \(perspective off\)")

orthoscopic on (perspective off) 

  * [![orthoscopic off \(perspective on\)](/images/f/f8/Orthooff2.png)](/index.php/File:Orthooff2.png "orthoscopic off \(perspective on\)")

orthoscopic off (perspective on) 

  * For many images, the difference is hardly visible.



  * [![orthoscopic on \(perspective off\)](/images/4/46/Orthoon.png)](/index.php/File:Orthoon.png "orthoscopic on \(perspective off\)")

orthoscopic on (perspective off) 

  * [![orthoscopic off \(perspective on\)](/images/3/3d/Orthooff.png)](/index.php/File:Orthooff.png "orthoscopic off \(perspective on\)")

orthoscopic off (perspective on) 

  * Especially when straight features parallel to the z axis are shown, the effect may be large.




  


## Syntax
    
    
    set orthoscopic, on
    set orthoscopic, off
    

## See Also

[Field of View](/index.php/Field_Of_View "Field Of View"), [Perspective Settings](/index.php/Ray#Perspective_Example_Images "Ray")

Retrieved from "[https://pymolwiki.org/index.php?title=Orthoscopic&oldid=8694](https://pymolwiki.org/index.php?title=Orthoscopic&oldid=8694)"


---

## Seq view

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

The **seq_view** command turns the sequence viewer on and off. The sequence viewer is a very handy tool, you can use it for example to select residues very easily. When the sequence viewer is turned on, you can select individual residues or multiple residues by selecting residues, using the mouse, on the sequence viewer. 

_Note:_ The sequence viewer, seems to decrease PyMOL's [performance](/index.php/Category:Performance "Category:Performance") when very many proteins are loaded into PyMOL. Therefore, if you have many proteins loaded, turn off the sequence viewer to increase performance if you don't need the viewer. You can remedy this, by turning on [texture_fonts](/index.php/Texture_fonts "Texture fonts"). Simply do, 
    
    
    set texture_fonts, 1
    

## Syntax
    
    
    # to turn the sequence viewer on.
    set seq_view, 1
    
    # to turn the sequence viewer off.
    set seq_view, 0
    

## Related settings

  * [seq_view_color](/index.php/Seq_view_color "Seq view color")
  * [seq_view_discrete_by_state](/index.php/Seq_view_discrete_by_state "Seq view discrete by state")
  * [seq_view_fill_char](/index.php/Seq_view_fill_char "Seq view fill char")
  * [seq_view_fill_color](/index.php?title=Seq_view_fill_color&action=edit&redlink=1 "Seq view fill color \(page does not exist\)")
  * [seq_view_format](/index.php/Seq_view_format "Seq view format")
  * [seq_view_label_color](/index.php/Seq_view_label_color "Seq view label color")
  * [seq_view_label_mode](/index.php/Seq_view_label_mode "Seq view label mode")
  * [seq_view_label_spacing](/index.php/Seq_view_label_spacing "Seq view label spacing")
  * [seq_view_label_start](/index.php?title=Seq_view_label_start&action=edit&redlink=1 "Seq view label start \(page does not exist\)")
  * [seq_view_location](/index.php/Seq_view_location "Seq view location")
  * [seq_view_overlay](/index.php/Seq_view_overlay "Seq view overlay")
  * [seq_view_unaligned_color](/index.php/Seq_view_unaligned_color "Seq view unaligned color")
  * [seq_view_unaligned_mode](/index.php/Seq_view_unaligned_mode "Seq view unaligned mode")



Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view&oldid=6507](https://pymolwiki.org/index.php?title=Seq_view&oldid=6507)"


---

## Png

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**png** writes a png format image file of the current image to disk. 

## Contents

  * 1 Usage
  * 2 Example
  * 3 PyMOL API
  * 4 Comments
    * 4.1 Transparent Backgrounds
    * 4.2 DPI Setting



## Usage
    
    
    png filename[, width[, height[, dpi[, ray[, quiet]]]]]
    

  * **filename** = string: file path to be written
  * **width** = integer or string: width in pixels (integer or string without units), inches (in), or centimeters (cm). If unit suffix is given, `dpi` argument is required as well. If only one of `width` or `height` is given, the aspect ratio of the viewport is preserved. {default: 0 (current)}
  * **height** = integer or string: height (see width) {default: 0 (current)}
  * **dpi** = float: dots-per-inch {default -1.0 (unspecified)}
  * **ray** = 0 or 1: should ray be run first {default: 0 (no)}
  * **quiet** = 0 or 1: if 1, logged output is suppressed. {default: 0}



## Example
    
    
    png ~/Desktop/test.png, width=10cm, dpi=300, ray=1
    

## PyMOL API
    
    
    cmd.png(string filename, int width=0, int height=0, float dpi=-1, int ray=0, int quiet=0)
    

## Comments

#### Transparent Backgrounds

Use the `[ray_opaque_background](/index.php/Ray_opaque_background "Ray opaque background")` setting to output images with transparent backgrounds. 
    
    
    set ray_opaque_background, 0
    

This can be useful for presentations, images that are placed on top of a background of nonuniform color (e.g. gradients), and images that overlap text or other images. 

#### DPI Setting

Use the DPI option to have PyMol set the DPI of your image. Executing the command 
    
    
    png /tmp/ex.png, width=1200, height=1200, dpi=300, ray=1
    

will ouput a four-inch square image at 300dpi. Leaving off the **dpi** parameter would yield a 1200x1200 image at your system's default pixel density (e.g. 72 or 96 dpi). This saves the intermediate step of having to use GIMP/PhotoShop/etc to rescale your photos for publication. 

Retrieved from "[https://pymolwiki.org/index.php?title=Png&oldid=11924](https://pymolwiki.org/index.php?title=Png&oldid=11924)"


---

## Ray

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**ray** creates a ray-traced image of the current frame. 

[![](/images/b/bc/Gslike.png)](/index.php/File:Gslike.png)

[](/index.php/File:Gslike.png "Enlarge")

Varying settings to play with rendering options

## Contents

  * 1 Details
    * 1.1 Usage
    * 1.2 PyMol API
    * 1.3 Settings
      * 1.3.1 Modes
      * 1.3.2 Perspective
        * 1.3.2.1 Perspective Example Images
          * 1.3.2.1.1 Notes
      * 1.3.3 Renderer
    * 1.4 Performance
      * 1.4.1 Memory
    * 1.5 Examples
      * 1.5.1 Simple
      * 1.5.2 Specify Image Size
      * 1.5.3 Specify Renderer
      * 1.5.4 High Quality B&W Rendering
      * 1.5.5 High Quality Color
      * 1.5.6 Ray Tracing Stereo Images
    * 1.6 See also
    * 1.7 User comments



# Details

This command is used to make high-resolution photos fit for publication and formal movies. Please note, the **ray** command can take some time (up to several minutes, depending on image complexity and size). 

For those who are making movies with PyMOL, **Ray** is one of the most commonly used last steps before stitching the frames together to compile the movie. **Ray** has many powerful features such as setting the size of the image -- and it still works even if the [Viewport](/index.php/Viewport "Viewport") or screen is smaller than the requested output file dimensions. 

[![](/images/b/b9/No_ray_trace.png)](/index.php/File:No_ray_trace.png)

[](/index.php/File:No_ray_trace.png "Enlarge")

Image, not ray traced.

[![](/images/3/30/Ray_traced.png)](/index.php/File:Ray_traced.png)

[](/index.php/File:Ray_traced.png "Enlarge")

Image, ray traced.

## Usage
    
    
    ray [width,height [,renderer [,angle [,shift ]]]
    

**angle** and **shift** can be used to generate matched stereo pairs 

**width** and **height** can be set to any non-negative integer. If both are set to zero than the current window size is used and is equivalent to just using **ray** with no arguments. If one is set to zero (or missing) while the other is a positive integer, then the argument set to zero (or missing) will be scaled to preserve the current aspect ratio. 

## PyMol API
    
    
    cmd.ray(int width,int height,int renderer=-1,float shift=0)
    

## Settings

### Modes

Setting the **[Ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")** variable in PyMOL changes the way PyMOL's internal renderer represents proteins in the final output. New modes were recently added to give the user more options of molecular representation. New modes are: normal rendering, but with a black outline (nice for presentations); black and white only; quantized color with black outline (also, very nice for presentations; more of a _cartoony_ appearance). 

**Note:** Mode 3, the quantized color one, sort of **burns** the background if you're using this setting. This will make a pure white background somewhat "offwhite"; thus, a poster would look poor because you could see the border for the image. If you'll be using this mode, try the [ray_opaque_background](/index.php/Ray_opaque_background "Ray opaque background") setting. 
    
    
    # normal color
    set ray_trace_mode, 0
    
    # normal color + black outline
    set ray_trace_mode,  1
    
    # black outline only
    set ray_trace_mode,  2
    
    # quantized color + black outline
    set ray_trace_mode,  3
    
    set ray_trace_mode, 1 # (or 2 or 3; best with "bg_color white;set antialias,2")
    # These two new modes -- 2 and 3 -- are cool cartoon looking modes.
    
    # change the color of the outline to a named color, or a hex-code
    set ray_trace_color, magenta
    set ray_trace_color, 0x0033ff
    

Here are the example images for the new modes 

  * [![set ray_trace_mode,1](/images/a/a8/Ray_mode_1_ex.png)](/index.php/File:Ray_mode_1_ex.png "set ray_trace_mode,1")

set ray_trace_mode,1 

  * [![set ray_trace_mode,2](/images/9/95/Ray_mode_2_ex.png)](/index.php/File:Ray_mode_2_ex.png "set ray_trace_mode,2")

set ray_trace_mode,2 

  * [![set ray_trace_mode,3](/images/5/51/Ray_mode_3_ex.png)](/index.php/File:Ray_mode_3_ex.png "set ray_trace_mode,3")

set ray_trace_mode,3 




### Perspective

#### Perspective Example Images

  * [![Example with Perspective Turned Off](/images/3/31/No_persp.png)](/index.php/File:No_persp.png "Example with Perspective Turned Off")

Example with Perspective Turned Off 

  * [![Example with Perspective Turned On](/images/c/cc/Persp1.png)](/index.php/File:Persp1.png "Example with Perspective Turned On")

Example with Perspective Turned On 

  * [![Example with Perspective Turned On and Field of View Set High \(70\).](/images/8/88/Persp2.png)](/index.php/File:Persp2.png "Example with Perspective Turned On and Field of View Set High \(70\).")

Example with Perspective Turned On and Field of View Set High (70). 




##### Notes

PyMol 0.97 and prior used **orthoscopic** rendering -- that is, no perspective. Upon the arrival of 0.98 and later, we get perspective based rendering at a cost of a 4x decrease in render speed. If you want perspective 
    
    
    set orthoscopic, off
    

Otherwise 
    
    
    set orthoscopic, on
    

To magnify the effect of perspective on the scene, 
    
    
    set field_of_view, X
    

where 50<X<70\. Default is 20. 50-70 gives a very strong perspective effect. Nb. the field of view is in Y, not X as one would expect. 

  


### Renderer

**renderer = -1** is default (use value in ray_default_renderer) 

**renderer = 0** uses PyMOL's internal renderer 

**renderer = 1** uses PovRay's renderer. This is Unix-only and you must have "povray" in your path. It utilizes two temporary files: "tmp_pymol.pov" and "tmp_pymol.png". Also works on Mac via Povray37UnofficialMacCmd but it needs to be in your path as "povray". 

## Performance

  * The ray performance depends on distance between camera and molecule.



If the distance is big rendering takes much time. If the distance is too small distant parts of molecule dissolve. 

  * [![Too close to molecule](/images/7/70/Close_ray.png)](/index.php/File:Close_ray.png "Too close to molecule")

Too close to molecule 

  * [![Normal distance](/images/1/1c/Middle_ray.png)](/index.php/File:Middle_ray.png "Normal distance")

Normal distance 



  * Tip: If you have a rather complicated scene that is zoomed into only a part of the molecule, you can speed up the ray tracing by hiding everything else outside of a certain range of the zoomed-on point. For example, if I have a large molecule and I'm looking only at the 30-atom ligand bound to it, then I can do something like the following:


    
    
    # setup your complex scene
    ...
    
    # zoom on the hetero atom (ligand and not water) within 5 Angstroms
    select hh, het and not resn HOH
    zoom hh, 5
    
    # turn on depth cueing
    set depth_cue, 1
    
    # now, select stuff to hide; we select everything that is 
    # farther than 8 Ang from our main selection
    select th, (all) and not ( (all) within 8 of hh) )
    
    hide everything, th
    
    # any additional commands you want
    ...
    
    ray
    

As an example of the efficacy of this method, I ray traced a rather complex scene with all the atoms visible here's the output of ray: 
    
    
    PyMOL>ray
     Ray: render time: 24.50 sec. = 146.9 frames/hour (941.88 sec. accum.).
    

and here is the result when I soft-clipped everything else using the above method: 
    
    
    PyMOL>ray
     Ray: render time: 47.93 sec. = 75.1 frames/hour (989.80 sec. accum.).
    

The two images in the following gallery show the results of the rendering. 

  * [![normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI.](/images/1/1a/Ray_method_off.png)](/index.php/File:Ray_method_off.png "normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI.")

normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI. 

  * [![manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left.](/images/7/7d/Ray_method_on.png)](/index.php/File:Ray_method_on.png "manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left.")

manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left. 




### Memory

If memory is an issue for you in PyMOL, try executing your rendering from a script rather than a PyMOL session file. An unfortunate unavoidable consequence of the fact that we use Python's portable, platform-independent "pickle" machinery for PyMOL session files. Packing or unpacking a Session or Scene file thus requires that there be two simultanous copies of the information to reside in RAM simultaneously: one native and a second in Python itself. 

So when memory is a limiting factor, scripts are recommended over sessions. 

## Examples

### Simple
    
    
    # ray trace the current scene using the default size of the viewport
    ray
    

### Specify Image Size
    
    
    # ray trace the current scene, but scaled to 1024x768 pixels
    ray 1024,768
    

### Specify Renderer
    
    
    # ray trace with an external renderer.
    ray renderer=0
    

### High Quality B&W Rendering

[![](/images/6/6c/1l9l.png)](/index.php/File:1l9l.png)

[](/index.php/File:1l9l.png "Enlarge")

Black and White (ray_trace_mode,2); click to see full image
    
    
    # Black and White Script
    load /tmp/3fib.pdb;
    show cartoon;
    set ray_trace_mode, 2;  # black and white cartoon
    bg_color white;
    set antialias, 2;
    ray 600,600
    png /tmp/1l9l.png
    

### High Quality Color

[![](/images/7/7e/1l9l_2.png)](/index.php/File:1l9l_2.png)

[](/index.php/File:1l9l_2.png "Enlarge")

Color mode (ray_trace_mode,3); click to see full image
    
    
    # Color Script
    load /tmp/thy_model/1l9l.pdb;
    hide lines;
    show cartoon;
    set ray_trace_mode, 3; # color
    bg_color white;
    set antialias, 2;
    remove resn HOH
    remove resn HET
    ray 600,600
    png /tmp/1l9l.png
    

### Ray Tracing Stereo Images

    _See[Stereo_Ray](/index.php/Stereo_Ray "Stereo Ray")_

## See also

  1. "help faster" for optimization tips with the builtin renderer. "help povray" for how to use PovRay instead of PyMOL's built-in ray-tracing engine. For high-quality photos, please also see the [Antialias](/index.php/Antialias "Antialias") command. [Ray shadows](/index.php/Ray_shadows "Ray shadows") for controlling shadows.
  2. See also [Ray Tracing](/index.php/Ray_Tracing "Ray Tracing").
  3. [Desaturation Tutorial](http://www.gimp.org/tutorials/Color2BW) \-- A good resource for making nice B&W images from color images (desaturation).
  4. [Ray Trace Gain](/index.php/Ray_Trace_Gain "Ray Trace Gain")



## User comments

How do I ray trace a publication-ready (~300dpi) image using PyMol?
    This answer is in the [Advanced Issues](/index.php/Category:Advanced_Issues "Category:Advanced Issues") (Image Manipulation Section).

Retrieved from "[https://pymolwiki.org/index.php?title=Ray&oldid=11725](https://pymolwiki.org/index.php?title=Ray&oldid=11725)"


---

## Read Molstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**read_molstr** reads an MDL MOL format file as a string 

### PYMOL API ONLY
    
    
    cmd.read_molstr( string molstr, string name, int state=0, \
      int finish=1, int discrete=1 )
    

### NOTES

  * **state** is a 1-based state index for the object, or 0 to append.


  * **finish** is a flag (0 or 1) which can be set to zero to improve [performance](/index.php/Category:Performance "Category:Performance") when loading large numbers of objects, but you must call **finish_object** when you are done.


  * **discrete** is a flag (0 or 1) which tells PyMOL that there will be no overlapping atoms in the file being loaded. **discrete** objects save memory but can not be edited.



Retrieved from "[https://pymolwiki.org/index.php?title=Read_Molstr&oldid=7602](https://pymolwiki.org/index.php?title=Read_Molstr&oldid=7602)"


---

## Read Pdbstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**read_pdbstr** in an API-only function which reads a pdb file from a Python string. This feature can be used to load or update structures into PyMOL without involving any temporary files. 

### PYMOL API ONLY
    
    
    cmd.read_pdbstr( string pdb-content, string object name 
       [ ,int state [ ,int finish [ ,int discrete ] ] ] )
    

### NOTES

**state** is a 1-based state index for the object. 

**finish** is a flag (0 or 1) which can be set to zero to improve [performance](/index.php/Category:Performance "Category:Performance") when loading large numbers of objects, but you must call "finish_object" when you are done. 

**discrete** is a flag (0 or 1) which tells PyMOL that there will be no overlapping atoms in the PDB files being loaded. **discrete** objects save memory but can not be edited. 

Retrieved from "[https://pymolwiki.org/index.php?title=Read_Pdbstr&oldid=7601](https://pymolwiki.org/index.php?title=Read_Pdbstr&oldid=7601)"


---

## Ribbon sampling

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting changes the number of segments that makes up a given length of ribbon. The default value is ten; larger values (more segments) provide a smoother curve while requiring greater amounts of memory and slowing [performance](/index.php/Category:Performance "Category:Performance"). 

## Settings
    
    
    set ribbon_sampling, 20   # increases the number of ribbon segments
    

## Examples

Open the images to actually see the details! 

  * [![ribbon sampling 5](/images/b/b4/Ribbon_05.png)](/index.php/File:Ribbon_05.png "ribbon sampling 5")

ribbon sampling 5 

  * [![ribbon sampling 10 \(default\)](/images/1/1a/Ribbon_10.png)](/index.php/File:Ribbon_10.png "ribbon sampling 10 \(default\)")

ribbon sampling 10 (default) 

  * [![ribbon sampling 20](/images/4/41/Ribbon_20.png)](/index.php/File:Ribbon_20.png "ribbon sampling 20")

ribbon sampling 20 

  * [![ribbon sampling 30](/images/0/04/Ribbon_30.png)](/index.php/File:Ribbon_30.png "ribbon sampling 30")

ribbon sampling 30 




Retrieved from "[https://pymolwiki.org/index.php?title=Ribbon_sampling&oldid=5502](https://pymolwiki.org/index.php?title=Ribbon_sampling&oldid=5502)"


---

## Ribbon smooth

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This PyMOL setting determines smoothly how PyMOL draws ribbons. Turning _on_ the setting smooths out the ribbons and make the image look cleaner. The setting is by default **off** , and yet, when turned **on** shows absolutely no decreased [performance](/index.php/Category:Performance "Category:Performance"). 

## Usage
    
    
    #show ribbons, first
    hide; show ribbon
    
    # turn on smoothing
    set ribbon_smooth, 1
    
    # turn off smoothing
    set ribbon_smooth, 0
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ribbon_smooth&oldid=5504](https://pymolwiki.org/index.php?title=Ribbon_smooth&oldid=5504)"


---

## Smooth

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Smooth performs a window average of coordinate states. 

# Usage
    
    
    smooth [ selection [, passes [, window [, first [, last [, ends]]]]]]
    
    # for example, smooth and object with 30 passes and
    # a window size of 100.
    smooth myObj, 30, 100
    

  * ends = 0 or 1: controls whether or not the end states are also smoothed using a weighted asymmetric window



To see smooth in context see this [example](http://www.pymolwiki.org/index.php/Protect#Example)

# NOTES

  * This type of averaging is often used to suppress high-frequency vibrations in a molecular dynamics trajectory.
  * This function is not memory efficient. For reasons of flexibility, it uses two additional copies of every atomic coordinate for the calculation. If you are memory-constrained in visualizing MD trajectories, then you may want to use an external tool such as ptraj to perform smoothing before loading coordinates into PyMOL.



# See Also

[Load_Traj](/index.php/Load_Traj "Load Traj"), [protect](/index.php/Protect "Protect"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Smooth&oldid=12096](https://pymolwiki.org/index.php?title=Smooth&oldid=12096)"


---

## Specular

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This controls whether or not specular reflections are shown unrendered, in the GUI. Turning this off only makes a slight improvement in speed. 

  * [![Specular set to 0](/images/6/64/Spec0.png)](/index.php/File:Spec0.png "Specular set to 0")

Specular set to 0 

  * [![Specular set to 1](/images/2/2c/Spec01.png)](/index.php/File:Spec01.png "Specular set to 1")

Specular set to 1 

  * [![Specular set to 5](/images/e/e1/Spec5.png)](/index.php/File:Spec5.png "Specular set to 5")

Specular set to 5 

  * [![Specular set to 10](/images/0/0a/Spec10.png)](/index.php/File:Spec10.png "Specular set to 10")

Specular set to 10 




Notice in the images above that the setting is controlling the amount of directly reflected light and not the shininess of the reflection. 

# Syntax
    
    
    # turn the setting on or off
    set specular, boolean
    
    # turn it off
    set specular, off
    
    # turn it on
    set specular, on
    

# See Also

[Pages on Specular Reflections](/index.php/Category:Specular_Reflections "Category:Specular Reflections")

Retrieved from "[https://pymolwiki.org/index.php?title=Specular&oldid=7156](https://pymolwiki.org/index.php?title=Specular&oldid=7156)"


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

## Texture fonts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

Not sure exactly what this does. The reason this page exists is to let the users know that this setting can affect rendering performance in PyMOL. 

To test how this affects your PyMOL. Turn on the frames per second (FPS) meter. Load a complex scene turn on the sequence viewer ([Seq_View](/index.php?title=Seq_View&action=edit&redlink=1 "Seq View \(page does not exist\)")) and then rotate the molecule. While doing this look at the maximum frame rate. Then turn on Texture_fonts and rotate the scene while watching the frame rate. Here are my example numbers: 

    

    **texture_fonts off** max = 142 FPS
    **texture_fonts on** max = 363 FPS

That's a huge difference. 

# Syntax
    
    
    # turn on texture fonts
    set texture_fonts, 1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Texture_fonts&oldid=6901](https://pymolwiki.org/index.php?title=Texture_fonts&oldid=6901)"


---

