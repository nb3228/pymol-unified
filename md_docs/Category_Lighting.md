# Category: Lighting

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

## Light count

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 PyMol API
  * 4 Example
  * 5 See Also



## Overview

set light_count defines the number of light sources.  
Setting to 0 or 1 removes any directional light source resulting in no shadows. Light_count does affect [ambient](/index.php/Ambient "Ambient") lighting. 

## Syntax
    
    
    set light_count, <integer>                          #default 2
    

## PyMol API
    
    
    cmd.set(light_count,int sources)
    

  


## Example

  * [![light_count 0](/images/f/ff/Light_count_0.png)](/index.php/File:Light_count_0.png "light_count 0")

light_count 0 

  * [![light_count 2 \(default\)](/images/b/b2/Light_count_2.png)](/index.php/File:Light_count_2.png "light_count 2 \(default\)")

light_count 2 (default) 

  * [![light_count 10](/images/3/35/Light_count_10.png)](/index.php/File:Light_count_10.png "light_count 10")

light_count 10 




# See Also

[Light](/index.php/Light "Light")

Retrieved from "[https://pymolwiki.org/index.php?title=Light_count&oldid=12005](https://pymolwiki.org/index.php?title=Light_count&oldid=12005)"


---

