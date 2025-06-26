# Category: Specular Reflections

## Spec direct

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting controls how PyMOL renders the direct specular reflection from a representation. 

This image shows an increasing Spec_direct setting from 0 to 1.0 

[![](/images/7/76/Spec_direct.gif)](/index.php/File:Spec_direct.gif)

[](/index.php/File:Spec_direct.gif "Enlarge")

Increasing the spec_direct setting from 0 to 1.0.

# Syntax
    
    
    # set spec direct to real-valued number
    set spec_direct, float
    
    # turn up the spec_direct
    set spec_direct, 0.75
    

# See Also

[Specular Reflections Category](/index.php/Category:Specular_Reflections "Category:Specular Reflections")

Retrieved from "[https://pymolwiki.org/index.php?title=Spec_direct&oldid=6199](https://pymolwiki.org/index.php?title=Spec_direct&oldid=6199)"


---

## Spec direct power

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This controls the direct specular power. Essentially, this is how strongly surfaces that are perpendicular to the ray drawn from the camera are rendered. The features seems to be inversely proportional to the value. Typical values here range from [0, 200+]. 

This image shows the spec_direct_power setting increasing **from** 0 **to** 100\. 

[![](/images/2/27/Spec_dir_power.gif)](/index.php/File:Spec_dir_power.gif)

[](/index.php/File:Spec_dir_power.gif "Enlarge")

spec_direct_power going from 0 to 100.

# Syntax
    
    
    # set it to some value
    set spec_direct_power, float
    
    # for example turn it up
    set spec_direct_power, 200
    

# See Also

[Specular Reflections Category](/index.php/Category:Specular_Reflections "Category:Specular Reflections")

Retrieved from "[https://pymolwiki.org/index.php?title=Spec_direct_power&oldid=6205](https://pymolwiki.org/index.php?title=Spec_direct_power&oldid=6205)"


---

## Spec power

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting controls the power of specular reflections. 

# Syntax
    
    
    # set the power
    set spec_power, float
    
    # turn up the reflections
    set spec_power, 200
    

# See Also

[Spec_reflect](/index.php/Spec_reflect "Spec reflect")

[Category:Specular_Reflections](/index.php/Category:Specular_Reflections "Category:Specular Reflections")

Retrieved from "[https://pymolwiki.org/index.php?title=Spec_power&oldid=6207](https://pymolwiki.org/index.php?title=Spec_power&oldid=6207)"


---

## Spec reflect

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



## Overview

This setting changes how PyMol handles specular reflections. Essentially, this setting, when combined with [spec_power](/index.php/Spec_power "Spec power") adjusts how sharp or diffuse the reflection is, giving you shiny looking surface or duller surfaces. 

## Syntax
    
    
    set spec_reflect, 2
    set spec_power, 1500
    

## Examples

  * [![Specular Reflection Example. Settings are spec_reflect=2.0; spec_power=1500](/images/b/b0/Spec.png)](/index.php/File:Spec.png "Specular Reflection Example. Settings are spec_reflect=2.0; spec_power=1500")

Specular Reflection Example. Settings are spec_reflect=2.0; spec_power=1500 

  * [![Specular Reflection Example. Settings are spec_reflect=0.0; spec_power=1500](/images/3/3d/Spec1.png)](/index.php/File:Spec1.png "Specular Reflection Example. Settings are spec_reflect=0.0; spec_power=1500")

Specular Reflection Example. Settings are spec_reflect=0.0; spec_power=1500 




## See Also

  * [spec_power](/index.php/Spec_power "Spec power")
  * [About Specular Reflections](http://en.wikipedia.org/wiki/Specular_reflection)



Retrieved from "[https://pymolwiki.org/index.php?title=Spec_reflect&oldid=8897](https://pymolwiki.org/index.php?title=Spec_reflect&oldid=8897)"


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

## Specular intensity

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This controls how PyMOL displays specular reflections in the GUI (unrendered). This is analogous to [Spec_power](/index.php/Spec_power "Spec power") for ray traced images. 

Values range form [0, 200+]. 

# Syntax
    
    
    # set the specular intensity to some real valued number
    set specular_intensity, float
    
    # for example, set it to 200
    set specular_intensity, 200.0
    

# See Also

[Spec_power](/index.php/Spec_power "Spec power")

Retrieved from "[https://pymolwiki.org/index.php?title=Specular_intensity&oldid=6209](https://pymolwiki.org/index.php?title=Specular_intensity&oldid=6209)"


---

