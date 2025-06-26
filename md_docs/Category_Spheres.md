# Category: Spheres

## Sphere color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

set sphere_color allows to color spheres independently from the rest of the representations. 

## Syntax
    
    
    set sphere_color, theColor
    

where _theColor_ can be : 

  * any usual colors (blue, yellow, grey50,...)
  * number-coded colors (1:black, 2:blue, 3:greenish, ...)
  * special code -1 to revert to original chameleon setting (_set sphere_color,-1_)



## Related settings

[cartoon_color](/index.php/Cartoon_color "Cartoon color")

Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_color&oldid=6300](https://pymolwiki.org/index.php?title=Sphere_color&oldid=6300)"


---

## Sphere mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

The **sphere_mode** setting controls how PyMOL draws spheres in the OpenGL interface. It does not affect rendering, using [ray](/index.php/Ray "Ray"). 

Sphere mode set to **5** , which looks beautiful, requires have [OpenGL Shaders Installed](/index.php/Spheres#Enabling_Shaders "Spheres"). 

## Syntax
    
    
    set sphere_mode, VAL
    

where 

VAL = [1..5]

## Example

  * [![Sphere mode unset](/images/9/97/Sphere_mode_off.png)](/index.php/File:Sphere_mode_off.png "Sphere mode unset")

Sphere mode unset 

  * [![Sphere mode 1](/images/2/22/Sphere_mode_1.png)](/index.php/File:Sphere_mode_1.png "Sphere mode 1")

Sphere mode 1 

  * [![Sphere mode 2](/images/4/47/Sphere_mode_2.png)](/index.php/File:Sphere_mode_2.png "Sphere mode 2")

Sphere mode 2 

  * [![Sphere mode 3](/images/9/91/Sphere_mode_3.png)](/index.php/File:Sphere_mode_3.png "Sphere mode 3")

Sphere mode 3 

  * [![Sphere mode 4](/images/4/42/Sphere_mode_4.png)](/index.php/File:Sphere_mode_4.png "Sphere mode 4")

Sphere mode 4 

  * [![Sphere mode 5](/images/3/33/Sphere_mode_5.png)](/index.php/File:Sphere_mode_5.png "Sphere mode 5")

Sphere mode 5 




Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_mode&oldid=6301](https://pymolwiki.org/index.php?title=Sphere_mode&oldid=6301)"


---

## Sphere point max size

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_point_max_size&oldid=6302](https://pymolwiki.org/index.php?title=Sphere_point_max_size&oldid=6302)"


---

## Sphere point size

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_point_size&oldid=6303](https://pymolwiki.org/index.php?title=Sphere_point_size&oldid=6303)"


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

## Sphere scale

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 Related settings



## Overview

This setting affects the apparent radius of spheres in the sphere representation. Default scale is set to 1.0. 

## Syntax
    
    
    # set the sphere scale to size for selection.
    cmd.set ("sphere_scale", size=1, selection='', state=0, updates=1, log=0, quiet=1)
    # generally it's simpler to use the console form
    set sphere_scale, size, selection
    # you can print the current value for sphere_scale
    get sphere_scale
    

_size_ can be any floating point number, _selection_ is the name of a selection. 

## Examples

Using 0.25 gives a nice balls&sticks representation with both lines and spheres turned on. 
    
    
    set_bond stick_color, white, (all), (all)
    set_bond stick_radius, 0.14, (all), (all)
    set sphere_scale, 0.25, (all)
    show sticks
    show spheres
    

**set sphere_scale** by itself will revert to default. Here you'll get a simple VDW rapresentation. 
    
    
    set_bond stick_color, blue, (all), (all)
    set_bond stick_radius, 0.3, (all), (all)
    set sphere_transparency, 0.3
    set sphere_scale
    

## Related settings

  * [sphere_color](/index.php/Sphere_color "Sphere color")
  * [sphere_mode](/index.php/Sphere_mode "Sphere mode")
  * [sphere_point_max_size](/index.php/Sphere_point_max_size "Sphere point max size")
  * [sphere_point_size](/index.php/Sphere_point_size "Sphere point size")
  * [sphere_quality](/index.php/Sphere_quality "Sphere quality")
  * [sphere_solvent](/index.php/Sphere_solvent "Sphere solvent")
  * [sphere_transparency](/index.php/Sphere_transparency "Sphere transparency")



Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_scale&oldid=8294](https://pymolwiki.org/index.php?title=Sphere_scale&oldid=8294)"


---

## Sphere solvent

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_solvent&oldid=6305](https://pymolwiki.org/index.php?title=Sphere_solvent&oldid=6305)"


---

## Sphere transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

set sphere_transparency is used to adjust the transparency of spheres! 
    
    
    set sphere_transparency, 0.5
    

or 
    
    
    set sphere_transparency=0.5, selection
    

Where 1.0 is invisible and 0.0 completely solid 

# Examples

  * [![Example of transparent spheres](/images/c/c2/Sphere_transparency_ex1.png)](/index.php/File:Sphere_transparency_ex1.png "Example of transparent spheres")

Example of transparent spheres 

  * [![Example of transparent spheres](/images/6/65/Sphere_transparency_ex3.png)](/index.php/File:Sphere_transparency_ex3.png "Example of transparent spheres")

Example of transparent spheres 




These images were made through the following script: 
    
    
    fetch 1ifr; color wheat; hide; show spheres;
    color red, i. 506-509;
    color marine, byres all within 5 of color red;
    set sphere_transparency,0.5,color marine;
    

Retrieved from "[https://pymolwiki.org/index.php?title=Sphere_transparency&oldid=6306](https://pymolwiki.org/index.php?title=Sphere_transparency&oldid=6306)"


---

