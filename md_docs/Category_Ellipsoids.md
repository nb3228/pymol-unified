# Category: Ellipsoids

## Ellipsoid color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting controls the color of ellipsoids. 

# Syntax
    
    
    # proper settings are colors
    set ellipsoid_color, color
    
    # for example set it to blue
    set ellipsoid_color, blue
    

# See Also

[Color](/index.php/Color "Color"), [The Coloring Category](/index.php/Category:Coloring "Category:Coloring")

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoid_color&oldid=6190](https://pymolwiki.org/index.php?title=Ellipsoid_color&oldid=6190)"


---

## Ellipsoid quality

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting determines how much refinement PyMOL uses in rendering ellipsoids. Other [representations](/index.php/Category:Representations "Category:Representations") have similar settings. The most often used is [Surface_quality](/index.php/Surface_quality "Surface quality"). 

This set to 0 is a rough approximation. Higher values, like 1, 2, 3 make truer representations. 

# Syntax
    
    
    # set the quality to some positive integer
    set ellipsoid_quality, int
    
    # for example, turn up the quality for rendering
    set ellipsoid_quailty, 3
    

# See Also

[Ray](/index.php/Ray "Ray"), [Surface_Quality](/index.php?title=Surface_Quality&action=edit&redlink=1 "Surface Quality \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoid_quality&oldid=8362](https://pymolwiki.org/index.php?title=Ellipsoid_quality&oldid=8362)"


---

## Ellipsoid scale

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This sets the relative scale of drawing [Ellipsoids](/index.php/Ellipsoids "Ellipsoids"). Larger numbers make larger ellipsoids. 

# Syntax
    
    
    # set to some real valued size
    set ellipsoid_scale, float
    
    # for example, magnify the ellipses
    # by drawing them 1.4x their normal size
    set ellipsoid_scale, 1.4
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoid_scale&oldid=6184](https://pymolwiki.org/index.php?title=Ellipsoid_scale&oldid=6184)"


---

## Ellipsoid transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This controls the level of [Transparency](/index.php/Transparency "Transparency") with which PyMOL renders ellipsiods. Setting this to 0 means full opacity (no transparency); setting this to 1.0 means fully transparent. 

# Syntax
    
    
    # set to a float in the range [0,1]
    set ellipsoid_transparency, float
    
    # for example, set it to 0.35
    set ellipsoid_transparency, 0.35
    

# See Also

[Category:Representations](/index.php/Category:Representations "Category:Representations"), [Transparency](/index.php/Transparency "Transparency"), [Ray](/index.php/Ray "Ray")

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoid_transparency&oldid=6183](https://pymolwiki.org/index.php?title=Ellipsoid_transparency&oldid=6183)"


---

