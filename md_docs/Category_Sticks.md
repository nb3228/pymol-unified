# Category: Sticks

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

## Stick ball

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

**This setting is deprecated in v1.5 and later as it is always enabled.**

The setting "stick_ball" controls whether bonded atoms are shown simply as joined sticks (set stick_ball, off) or as traditional "ball-and-stick" representation (set stick_ball, on). Note that simply setting stick_ball on will result in balls with the same radius as the sticks and so will appear only slightly different (the joins will be smoother). 

## Settings
    
    
    set stick_ball, on   # displays atoms as balls joined by sticks
    set stick_ball, off  # displays only connected sticks
    
    set stick_ball_ratio, 1.7 # change the radius of the balls
    

## Examples

Open the images to actually see the details! 

  * [![stick_ball, off](/images/c/cc/Stick_ball_off.png)](/index.php/File:Stick_ball_off.png "stick_ball, off")

stick_ball, off 

  * [![stick_ball, on](/images/0/00/Stick_ball_on.png)](/index.php/File:Stick_ball_on.png "stick_ball, on")

stick_ball, on 

  * [![stick_ball "on" with stick_ball_ratio at 1.5](/images/4/4d/Stick_ball_ratio_1.5.png)](/index.php/File:Stick_ball_ratio_1.5.png "stick_ball "on" with stick_ball_ratio at 1.5")

stick_ball "on" with stick_ball_ratio at 1.5 




Retrieved from "[https://pymolwiki.org/index.php?title=Stick_ball&oldid=10909](https://pymolwiki.org/index.php?title=Stick_ball&oldid=10909)"


---

## Stick ball ratio

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Settings
  * 3 Related
  * 4 Examples



## Overview

The setting "stick_ball_ratio" controls the relative ratio between the radius of sticks connecting bonded atoms and the radius of the atom spheres. Note that simply setting "stick_ball, on" will result in balls with the same radius as the sticks and so will appear only slightly different (the joins will be smoother). Changing the stick_ball_ratio without setting "stick_ball, on" will -- obviously -- have no apparent effect. 

## Settings
    
    
    set stick_ball_ratio, 1.5
    

## Related
    
    
    set stick_ball, on   # displays atoms as balls joined by sticks
    set stick_ball, off  # displays only connected sticks
    

## Examples

Open the images to actually see the details! 

  * [![stick_ball, off](/images/c/cc/Stick_ball_off.png)](/index.php/File:Stick_ball_off.png "stick_ball, off")

stick_ball, off 

  * [![stick_ball, on](/images/0/00/Stick_ball_on.png)](/index.php/File:Stick_ball_on.png "stick_ball, on")

stick_ball, on 

  * [![stick_ball "on" with stick_ball_ratio at 1.5](/images/4/4d/Stick_ball_ratio_1.5.png)](/index.php/File:Stick_ball_ratio_1.5.png "stick_ball "on" with stick_ball_ratio at 1.5")

stick_ball "on" with stick_ball_ratio at 1.5 




Retrieved from "[https://pymolwiki.org/index.php?title=Stick_ball_ratio&oldid=5270](https://pymolwiki.org/index.php?title=Stick_ball_ratio&oldid=5270)"


---

## Stick color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Example
  * 4 See Also



## Overview

_Stick_color_ controls the color of sticks for a selection or object. For example, one could color all atoms by their element as transparent spheres, but for a nice underlying structure, have all the sticks shown in marine. 

## Syntax
    
    
    set stick_color, color
    

## Example

As noted in the above text -- transparent spheres, where sticks are all one color: 

  * [![Example of color sticks one color](/images/3/33/Stick_color.png)](/index.php/File:Stick_color.png "Example of color sticks one color")

Example of color sticks one color 




## See Also

[Color](/index.php/Color "Color")

Retrieved from "[https://pymolwiki.org/index.php?title=Stick_color&oldid=5272](https://pymolwiki.org/index.php?title=Stick_color&oldid=5272)"


---

## Stick nub

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

The stick_nub setting controls the height of the cones at the end of the sticks in stick representation. This only affects the simple representation: this setting doesn't have any effect for ray tracing. 

  * [![Default stick_nub setting.](/images/2/28/Sticknub_default.png)](/index.php/File:Sticknub_default.png "Default stick_nub setting.")

Default stick_nub setting. 

  * [![Stubby nubs--stick nub set to 0.](/images/5/55/Sticknub_0.png)](/index.php/File:Sticknub_0.png "Stubby nubs--stick nub set to 0.")

Stubby nubs--stick nub set to 0. 

  * [![Nub spears: stick_nub=2.0](/images/f/f4/Sticknub_2.png)](/index.php/File:Sticknub_2.png "Nub spears: stick_nub=2.0")

Nub spears: stick_nub=2.0 




## Syntax
    
    
    set stick_nub, ''float''
    

Where _float_ is a floating point number. The default value is 0.7. 

## Example
    
    
    set stick_nub, 0
    

Retrieved from "[https://pymolwiki.org/index.php?title=Stick_nub&oldid=5281](https://pymolwiki.org/index.php?title=Stick_nub&oldid=5281)"


---

## Stick radius

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 If the above commands do not work
  * 4 Related settings



## Overview

This setting affects the radius of sticks in the sticks representation. Default scale is set to 0.25. 

In newer versions of PyMOL, one may set the Stick_radius on a per-bond basis. So, you can set for example, the radius of only selected bonds if you want. This is done through the [Set_bond](/index.php/Set_bond "Set bond") command. 

  * [![stick_radius set to 0.05](/images/1/18/Stick_rad_0.05.png)](/index.php/File:Stick_rad_0.05.png "stick_radius set to 0.05")

stick_radius set to 0.05 

  * [![stick_radius set to the default 0.25](/images/d/de/Stick_radius_default.png)](/index.php/File:Stick_radius_default.png "stick_radius set to the default 0.25")

stick_radius set to the default 0.25 

  * [![stick_radius set to 0.85](/images/4/41/Stick_rad_0.85.png)](/index.php/File:Stick_rad_0.85.png "stick_radius set to 0.85")

stick_radius set to 0.85 




## Syntax
    
    
    set_bond stick_radius, ''size'', selection
    

where, 

  * _size_ can be any float number. Using 0.25 (default value) is usually appropriate for most representations, although 0.15 migh be preferred for comparing closely related structures, e.g., conformers.



_Note:_
    
    
    set_bond stick_radius
    

by itself will revert to 1.00. 

## If the above commands do not work

You can do something like below 
    
    
    To set on the entire object
    
    set stick_radius=0.12
    
    OR
    
    create myObj, <selection>
    
    Ex : create myObj, hetatm
    
    set stick_radius,0.2,myObj
    

## Related settings

  * [Set_bond](/index.php/Set_bond "Set bond")
  * [sphere_scale](/index.php/Sphere_scale "Sphere scale")
  * [stick_ball_ratio](/index.php/Stick_ball_ratio "Stick ball ratio")



Retrieved from "[https://pymolwiki.org/index.php?title=Stick_radius&oldid=7204](https://pymolwiki.org/index.php?title=Stick_radius&oldid=7204)"


---

## Stick transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 Note



## Overview

The setting "stick_transparency" allows one to set the degree of transparency for stick objects, independent of all other objects. Allowable values range from 0 (fully opaque) to 1 (fully transparent, i.e. invisible). 

## Usage
    
    
    set stick_transparency, F, objectname
    set_bond stick_transparency, F, selection
    

where **F** is a floating point number in the range _[0.0 - 1.0]_ , and **objectname** is the object to apply the changes to. If ObjectName is omitted, then the transparency of all stick representations will be changed. 

To apply this setting to a selection, use "set_bond" syntax 
    
    
    set_bond stick_transparency, 0.7, */n+c+ca+o
    

For the value of _F_ =1.0 sticks will be completely transparent/invisible and for _F_ =0.0 sticks are solid/opaque. 

## Examples
    
    
    set stick_transparency, 0.50   # Makes sticks 50-percent transparent
    

Open the images to actually see the details! 

  * [![sticks with no transparency](/images/b/bf/Stick_trans_zero.png)](/index.php/File:Stick_trans_zero.png "sticks with no transparency")

sticks with no transparency 

  * [![sticks with 0.5 transparency](/images/9/99/Stick_trans_50.png)](/index.php/File:Stick_trans_50.png "sticks with 0.5 transparency")

sticks with 0.5 transparency 

  * [![sticks with 0.9 transparency](/images/e/e0/Stick_trans_90.png)](/index.php/File:Stick_trans_90.png "sticks with 0.9 transparency")

sticks with 0.9 transparency 




## Note

Stick transparency works best with "unilayer" transparency (Setting menu > transparency > Unilayer) rather than "multilayer", which leads to odd artifacts where the sticks join. 

  * [![sticks with 0.5 transparency and multilayer transparency](/images/a/a0/Stick_trans_50-multi.png)](/index.php/File:Stick_trans_50-multi.png "sticks with 0.5 transparency and multilayer transparency")

sticks with 0.5 transparency and multilayer transparency 




Retrieved from "[https://pymolwiki.org/index.php?title=Stick_transparency&oldid=12792](https://pymolwiki.org/index.php?title=Stick_transparency&oldid=12792)"


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

