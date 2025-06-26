# Category: Lines

## Line color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

Line_color controls per atom/bond coloring in objects or selections. 

[![](/images/b/ba/Lw2.png)](/index.php/File:Lw2.png)

[](/index.php/File:Lw2.png "Enlarge")

Line_color changed using set_bond. (Also, line_width was changed using set_bond)

# Syntax
    
    
    # set the color of the lines to colorName
    set line_color, colorName
    
    # set per atom/bond line colors (see examples)
    set_bond line_color, colorName, objSel
    
    # example
    fetch 1te1
    as lines
    orient
    # draw all lines red
    set line_color, red
    # draw just chain A, blue
    set_bond line_color, marine, blue
    # color the lysines magnesium!
    set_bond line_color, magnesium, resn lys
    

# See Also

  * [Set](/index.php/Set "Set")
  * [Set_bond](/index.php/Set_bond "Set bond")
  * [Lines](/index.php/Lines "Lines")
  * [Line_width](/index.php/Line_width "Line width")



Retrieved from "[https://pymolwiki.org/index.php?title=Line_color&oldid=7012](https://pymolwiki.org/index.php?title=Line_color&oldid=7012)"


---

## Line smooth

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This is a display setting. 

The **line_smooth setting** determines whether lines are **[antialias](/index.php/Antialias "Antialias")** ed. The default is to use antialiased lines. 

## Syntax

To turn off antialiasing: 
    
    
    set line_smooth, 0
    

To turn on antialising (default): 
    
    
    set line_smooth, 1
    

From python: 
    
    
    cmd.set('line_smooth', '0')
    

or 
    
    
    cmd.set('line_smooth', '1')
    

Retrieved from "[https://pymolwiki.org/index.php?title=Line_smooth&oldid=5284](https://pymolwiki.org/index.php?title=Line_smooth&oldid=5284)"


---

