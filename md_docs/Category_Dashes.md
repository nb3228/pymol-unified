# Category: Dashes

## Dash color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

Sets the color of dashes. The default is a yellow. You can use any color PyMOL knows about. 

  * [![Dash_color set to 'red'](/images/f/f3/DashColorRed.png)](/index.php/File:DashColorRed.png "Dash_color set to 'red'")

Dash_color set to 'red' 

  * [![Dash_color set to 'marine'](/images/c/c9/DashColorMarine.png)](/index.php/File:DashColorMarine.png "Dash_color set to 'marine'")

Dash_color set to 'marine' 

  * [![PyMOL accepts hex notation for colors: set dash_color, 0x6644FF used here.](/images/6/60/DashColor64f.png)](/index.php/File:DashColor64f.png "PyMOL accepts hex notation for colors: set dash_color, 0x6644FF used here.")

PyMOL accepts hex notation for colors: set dash_color, 0x6644FF used here. 




# Syntax
    
    
    # set dash_color to 'color'--any valid color
    set dash_color, color
    
    # for example
    set dash_color, red
    set dash_color, marine
    set dash_color, 0xFFFFFF
    
    # for example
    

# See Also

[Color](/index.php/Color "Color"), [Label](/index.php/Label "Label"), [Dash_gap](/index.php/Dash_gap "Dash gap"), [Dash_width](/index.php/Dash_width "Dash width")

Retrieved from "[https://pymolwiki.org/index.php?title=Dash_color&oldid=6112](https://pymolwiki.org/index.php?title=Dash_color&oldid=6112)"


---

## Dash gap

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

Sets the distance (gap) between individual dashes. 

  * [![dash_gap set to 0.](/images/d/d8/DashGap0.png)](/index.php/File:DashGap0.png "dash_gap set to 0.")

dash_gap set to 0. 

  * [![dash_gap set to 3.0--very big gaps](/images/8/85/DashGap3.png)](/index.php/File:DashGap3.png "dash_gap set to 3.0--very big gaps")

dash_gap set to 3.0--very big gaps 




  


# Syntax
    
    
    # set dash_gap to 0.20
    set dash_gap, 0.20
    
    # make a line; set the dash_gap to 0
    set dash_gap, 0
    
    # make large gaps
    set dash_gap, 
    
    # for example
    

# See Also

[Dash_width](/index.php/Dash_width "Dash width"), [Label](/index.php/Label "Label")

Retrieved from "[https://pymolwiki.org/index.php?title=Dash_gap&oldid=6107](https://pymolwiki.org/index.php?title=Dash_gap&oldid=6107)"


---

## Dash length

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See also



## Overview

Sets the length of dashes in labels. 

## Syntax
    
    
    set dash_length, 0.3500 # make longer dashes
    

## Examples

  * [![Dash Length set to 0.1500](/images/8/88/Dashes3500.png)](/index.php/File:Dashes3500.png "Dash Length set to 0.1500")

Dash Length set to 0.1500 

  * [![Dash Length set to 0.3500](/images/2/2f/Dl3500.png)](/index.php/File:Dl3500.png "Dash Length set to 0.3500")

Dash Length set to 0.3500 




## See also

[Dash Radius](/index.php/Dash_Radius "Dash Radius"), [Dash Width](/index.php?title=Dash_Width&action=edit&redlink=1 "Dash Width \(page does not exist\)"), [Dash Gap](/index.php/Dash_Gap "Dash Gap"), [Dash Round Ends](/index.php/Dash_Round_Ends "Dash Round Ends")

Retrieved from "[https://pymolwiki.org/index.php?title=Dash_length&oldid=12087](https://pymolwiki.org/index.php?title=Dash_length&oldid=12087)"


---

## Dash round ends

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting determines if PyMOL draws the ends of the gaps as flat or rounded. The default is set to 'off'. 

# Syntax
    
    
    # set dash_round_ends to 'on' or 'off'.
    set dash_round_ends, boolean
    
    # for example round the dash ends
    set dash_round_ends, on
    
    # set flat ends on
    set dash_round_ends, off
    

# See Also

[Label](/index.php/Label "Label")

Retrieved from "[https://pymolwiki.org/index.php?title=Dash_round_ends&oldid=6113](https://pymolwiki.org/index.php?title=Dash_round_ends&oldid=6113)"


---

## Dash width

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

PyMOL has great control over its representation of dashes. You can use the default Dash_width, or 

  * [![dash_width increasing from 1 to 9.](/images/5/5c/DashWidthEx.jpg)](/index.php/File:DashWidthEx.jpg "dash_width increasing from 1 to 9.")

dash_width increasing from 1 to 9. 




# Syntax
    
    
    # set dash_width to some positive value, default value 2.5
    set dash_width, float
    
    # for example,
    set dash_width, 4
    

# See Also

[Dash Radius](/index.php/Dash_Radius "Dash Radius"), [Dash Gap](/index.php/Dash_Gap "Dash Gap"), [Dash Length](/index.php/Dash_Length "Dash Length"), [Dash Round Ends](/index.php/Dash_Round_Ends "Dash Round Ends")

Retrieved from "[https://pymolwiki.org/index.php?title=Dash_width&oldid=12035](https://pymolwiki.org/index.php?title=Dash_width&oldid=12035)"


---

