# Category: Seqeuence Viewer

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

## Seq view alignment

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

# Syntax

# See Also

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_alignment&oldid=6182](https://pymolwiki.org/index.php?title=Seq_view_alignment&oldid=6182)"


---

## Seq view color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

This setting controls the color that the sequences in the sequence viewer are displayed. 

# Syntax
    
    
    # set the color to 'color'
    # valid colors are any PyMOL-known colors
    set seq_view_color, color
    
    # set the seq view color to marine
    set seq_view_color, marine
    

# See Also

[Color](/index.php/Color "Color")

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_color&oldid=6226](https://pymolwiki.org/index.php?title=Seq_view_color&oldid=6226)"


---

## Seq view discrete by state

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

# Syntax

# See Also

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_discrete_by_state&oldid=6181](https://pymolwiki.org/index.php?title=Seq_view_discrete_by_state&oldid=6181)"


---

## Seq view fill char

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

This set the character used for filler. The default is '-'. 

# Syntax
    
    
    # set the fill character to some character, c
    set seq_view_fill_char, c
    
    # for example, set the fill char to X
    set seq_view_fill_char, X
    

# See Also

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_fill_char&oldid=6279](https://pymolwiki.org/index.php?title=Seq_view_fill_char&oldid=6279)"


---

## Seq view format

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

### OVERVIEW

The **seq_view_format** setting controls how PyMOL displays the sequence viewer. (To turn on the sequence viewer, type, **set seq_view, on**.) The available formats are currently: 

  * 0 = Display residues as single letter amino acid names
  * 1 = Display residues as triple letter amino acid names
  * 2 = Displays all atoms in each residue based on their atom name
  * 3 = Displays each peptide chain



### SYNTAX
    
    
    # usage
    set seq_view_format, number
    
    # default
    set seq_view_format, 0
    
    # triple letter amino acids
    set seq_view_format, 1
    

### SEE ALSO

[Seq_view](/index.php/Seq_view "Seq view")

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_format&oldid=6293](https://pymolwiki.org/index.php?title=Seq_view_format&oldid=6293)"


---

## Seq view gap mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The seq_view_gap_mode setting controls if gap indicators are displayed in the [sequence viewer](/index.php/Category:Seqeuence_Viewer "Category:Seqeuence Viewer"). 

_New in PyMOL 2.3_

## Values

  * 0: no gap indicator display
  * 1: number of dashes equals number of missing residues (based on residue numbers) {default}
  * 2: one dash per gap (independent of size)



## Example
    
    
    fetch 2xwu, type=pdb, async=0
    set seq_view_gap_mode, 1
    set seq_view
    

Scroll sequence viewer to chain B residue 152, it should display 3 dashes. 

## See Also

  * [seq_view](/index.php/Seq_view "Seq view")



Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_gap_mode&oldid=12863](https://pymolwiki.org/index.php?title=Seq_view_gap_mode&oldid=12863)"


---

## Seq view label color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

### OVERVIEW

The **Seq_view_label_color** setting controls the setting of the color that the information about the sequence in displayed in. 

### USAGE
    
    
    # make the information cyan
    set seq_view_label_color, cyan
    

### SEE ALSO

[Seq_view](/index.php/Seq_view "Seq view")

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_label_color&oldid=6225](https://pymolwiki.org/index.php?title=Seq_view_label_color&oldid=6225)"


---

## Seq view label mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

### OVERVIEW

The **seq_view_label_mode** setting controls the display of the sequence viewer in PyMOL. Currently, the settings are: 

  * 0 = Just the object name and chosen sequence format
  * 1 = The object name and chosen sequence format along with the residue numbering and chain name
  * 2 = Residue numbering and chain name
  * 3 = Just the residue sequence (no numbers, object names, etc.)



### USAGE
    
    
    # turn off everything but the sequence display
    set seq_view_label_mode, 3
    

### SEE ALSO

[Seq_view](/index.php/Seq_view "Seq view")

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_label_mode&oldid=6171](https://pymolwiki.org/index.php?title=Seq_view_label_mode&oldid=6171)"


---

## Seq view label spacing

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

This controls how the sequence indices are shown in the sequence viewer. Setting this to 1 means every index will be shown, setting this to 5, means that every 5th residue index will be shown in the sequence viewer. 

  * [![The spacing here is 5.](/images/0/06/Svl5.png)](/index.php/File:Svl5.png "The spacing here is 5.")

The spacing here is 5. 

  * [![The spacing here is 10.](/images/c/ca/Svl10.png)](/index.php/File:Svl10.png "The spacing here is 10.")

The spacing here is 10. 




Note: Don't forget to turn on the sequence viewer, **set seq_view**. 

# Syntax
    
    
    # set to any positive integer
    set seq_view_label_spacing, int
    
    # for example, set it to 5
    set seq_view_label_spacing, 5
    

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_label_spacing&oldid=6178](https://pymolwiki.org/index.php?title=Seq_view_label_spacing&oldid=6178)"


---

## Seq view location

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

seq_view_location affects the location of the sequence viewer. 

## Syntax

set seq_view_location, _integer_

set _integer_ to 0 to get the default (the viewer is on top of the internal window) 

set _integer_ to 1 or higher to locate the sequence viewer just above the command line. 

  


## See Also

  * [seq_view](/index.php/Seq_view "Seq view")
  * [seq_view_color](/index.php/Seq_view_color "Seq view color")
  * [seq_view_discrete_by_state](/index.php/Seq_view_discrete_by_state "Seq view discrete by state")
  * [seq_view_fill_char](/index.php/Seq_view_fill_char "Seq view fill char")
  * [seq_view_fill_color](/index.php?title=Seq_view_fill_color&action=edit&redlink=1 "Seq view fill color \(page does not exist\)")
  * [seq_view_format](/index.php/Seq_view_format "Seq view format")
  * [seq_view_label_color](/index.php/Seq_view_label_color "Seq view label color")
  * [seq_view_label_mode](/index.php/Seq_view_label_mode "Seq view label mode")
  * [seq_view_label_spacing](/index.php/Seq_view_label_spacing "Seq view label spacing")
  * [seq_view_label_start](/index.php?title=Seq_view_label_start&action=edit&redlink=1 "Seq view label start \(page does not exist\)")
  * [seq_view_overlay](/index.php/Seq_view_overlay "Seq view overlay")
  * [seq_view_unaligned_color](/index.php/Seq_view_unaligned_color "Seq view unaligned color")
  * [seq_view_unaligned_mode](/index.php/Seq_view_unaligned_mode "Seq view unaligned mode")



Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_location&oldid=6172](https://pymolwiki.org/index.php?title=Seq_view_location&oldid=6172)"


---

## Seq view overlay

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

### OVERVIEW

The **seq_view_overaly** setting controls whether or not the background of the sequence viewer is opaque or not, allowing access to the underlying molecules. 

  


### USAGE
    
    
    # make the sequence viewer background transparent
    set seq_view_overlay, 1
    
    # make it opaque again
    set seq_view_overlay, 0
    

### SEE ALSO

[Seq_view](/index.php/Seq_view "Seq view")

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_overlay&oldid=6173](https://pymolwiki.org/index.php?title=Seq_view_overlay&oldid=6173)"


---

## Seq view unaligned color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

Sets the colors of unaligned residues in the sequence viewer. 

# Syntax
    
    
    # set the color
    set seq_view_unaligned_color, color
    
    # for example
    set seq_view_unaligned_color, cyan
    

# See Also

[Color](/index.php/Color "Color")

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_unaligned_color&oldid=6227](https://pymolwiki.org/index.php?title=Seq_view_unaligned_color&oldid=6227)"


---

## Seq view unaligned mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overivew

# Syntax

# See Also

Retrieved from "[https://pymolwiki.org/index.php?title=Seq_view_unaligned_mode&oldid=6180](https://pymolwiki.org/index.php?title=Seq_view_unaligned_mode&oldid=6180)"


---

