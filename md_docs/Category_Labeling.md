# Category: Labeling

## Label

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![Pl.png](/images/b/b9/Pl.png)](/index.php/File:Pl.png)

The Label command controls how PyMOL draws text labels for PyMOL objects. 

## Contents

  * 1 Details
    * 1.1 Built-in Object Properties
  * 2 Syntax
  * 3 Settings
  * 4 Examples
  * 5 User Comments
    * 5.1 Labels Using ID Numbers
    * 5.2 Labels Using One Letter Abbreviations
    * 5.3 Labels and defer_builds_mode
  * 6 See Also



# Details

Labeling is important so there are many options for your fine tuning needs. You can change the [label size](/index.php/Label_size "Label size"), [label color](/index.php/Label_color "Label color"), positioning, [font](/index.php/Label_font_id "Label font id"), the [label outline color](/index.php/Label_outline_color "Label outline color") that masks the font and much, much more. 

You can have PyMOL label atoms by properties or arbitrary strings as you want; you can even use Unicode fonts for special symbols like,  α , β , ± , \mathrm {\AA}  {\displaystyle \alpha ,\beta ,\pm ,{\textrm {\mathrm {\AA} }}} ![{\\displaystyle \\alpha ,\\beta ,\\pm ,{\\textrm {\\mathrm {\\AA} }}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/07cd06e9fdbee0a625881a9a0ee07dadcde15a66), etc. 

The following gallery shows some examples of how extensible the Label command is. 

  * [![Simple label](/images/7/78/Label_pre.png)](/index.php/File:Label_pre.png "Simple label")

Simple label 

  * [![Example showing usage of Unicode fonts for special characters, see label_font_id.](/images/e/e0/New_fonts.jpeg)](/index.php/File:New_fonts.jpeg "Example showing usage of Unicode fonts for special characters, see label_font_id.")

Example showing usage of Unicode fonts for special characters, see [label_font_id](/index.php/Label_font_id "Label font id"). 

  * [![Another example with Unicode fonts](/images/2/2f/Font_ex.png)](/index.php/File:Font_ex.png "Another example with Unicode fonts")

Another example with Unicode fonts 

  * [![Example label](/images/c/c2/Label_ex.png)](/index.php/File:Label_ex.png "Example label")

Example label 

  * [![Label shadows turned off](/images/0/05/Ls0.png)](/index.php/File:Ls0.png "Label shadows turned off")

Label shadows turned off 

  * [![Label shadows turned on](/images/8/89/Ls2.png)](/index.php/File:Ls2.png "Label shadows turned on")

Label shadows turned on 




## Built-in Object Properties

Aside from arbitrary string labels, like "This is the catalytic residue" for an atom/residue you can also use the following built-in molecular properties: 

  * **name** , the atom name
  * **resn** , the residue name
  * **resi** , the residue number/identifier
  * **chain** , the chain name
  * **q** , charge
  * **b** , the occupancy/b-factor
  * **segi** , the segment identifier
  * **type** _(ATOM,HETATM)_ , the type of atom
  * **formal_charge** , the formal charge
  * **partial_charge** , the partial charge
  * **numeric_type** , the numeric type
  * **text_type** , the text type



You can use one of these properties as: 
    
    
    # simple example: label residue 22's atoms with their names
    label i. 22, name
    
    # Label residue #44's alpha carbon with it's residue name, residue number and B-factor.
    label n. CA and i. 44, "(%s, %s, %s)" % (resn, resi, b)
    

See the syntax and examples below for more info. 

# Syntax

To use the label command follow this syntax: 
    
    
    # labeling syntax
    label [ selection[, expression]]
    

where **selection** is some object/selection you want to label and **expression** is some string (or set of strings) which PyMOL is to use to label the given selection. 

We have plenty of examples. See the examples below. 

# Settings

Here are all the label settings and their general effect. For each label setting, see the respective web page for more details. 

**[label_angle_digits](/index.php/Label_angle_digits "Label angle digits")**

    

    sets the number of decimals in angle label.

**[label_distance_digits](/index.php/Label_distance_digits "Label distance digits")**

    

    sets the number of decimals in distance label.

**[label_shadow_mode](/index.php/Label_shadow_mode "Label shadow mode")**

    

    sets whether or not PyMOL will ray trace shadows for your label text. Eg: 
    
    
    set label_shadow_mode, 2
    

**[label_color](/index.php/Label_color "Label color")**

    

    sets the color of the label text. Note that you can have labels of different colors for different objects or selections. Some examples:
    
    
    # per-object:
    set label_color, color-name, object-name  #eg, set label-color, magenta, /protein
    
    # per-atom:
    set label_color, color-name, selection    #eg, set label-color, marine, /protein/A/A/23/CA
    
    # another example
    fragment arg
    label all, name
    set label_color, yellow, arg
    set label_color, red, elem c
    

**[label_font_id](/index.php/Label_font_id "Label font id")**

    

    sets the font to render your label. There are 12 different fonts from 5—16. Numbers 15 and 16 are special for unicode. Eg: 
    
    
    set label_font_id, 12
    

. See the [label_font_id](/index.php/Label_font_id "Label font id") page for explicit examples on how to use unicode characters in PyMOL labels.

**[label_size](/index.php/Label_size "Label size")**

    

    sets the size of the text. You can use positive numbers 2, 3, 4, etc for point sizes, or negative numbers for Angstrom-based sizes. Default is 14 points. Labels in Angstrom-size scale with the distance from the front plane, labels in point-size don't. Eg: 
    
    
    set label_size, -2  #results in a size of 2 Angstroms
    

**[label_digits](/index.php?title=Label_digits&action=edit&redlink=1 "Label digits \(page does not exist\)")**

    

    sets the number of decimals in label. It affects all digits only if label_distance_digits or label_dihedral_digits or label_angle_digits are set to -1.

**[label_outline_color](/index.php/Label_outline_color "Label outline color")**

    

    each label is outlined (so you can do white-on-white labels, for example). This options sets the color of the label outline. Eg. 
    
    
    set label_outline_color, orange
    

**[label_dihedral_digits](/index.php/Label_dihedral_digits "Label dihedral digits")**

    

    sets the number of decimals in dihedral label.

**[label_position](/index.php/Label_position "Label position")**

    

    sets any offset from the original X,Y,Z coordinates for the label. If you like to use the mouse, you can enter [edit_mode](/index.php?title=Edit_mode&action=edit&redlink=1 "Edit mode \(page does not exist\)") and **ctrl-left_click** to drag labels around; **ctrl-shift-left_click** will let you move the labels in the z-direction. **"Save labels"-workaround** If you want to save the position of your labels, the best way might be to create a new object and move the atoms in this object. Since the labels are positioned from the atom positions this is an indirect way of moving the labels and being able to save them.

# Examples
    
    
    #1.
    # make a very simple label on the 14th alpha carbon.
    label n. CA and i. 14, "This is carbon 14."
    
    #2.
    # make a fake scene label; use this to label entire scenes, not just atoms/bonds.
    pseudoatom foo
    label foo, "Once upon a time..."
    
    #3.
    # make a huge label
    set label_size, -5
    pseudoatom foo
    label foo, "This is large text"
    
    #4. Partial Charge
    label (chain A),chain
    label (n;ca),"%s-%s" % (resn,resi)
    label (resi 200),"%1.3f" % partial_charge
    
    
    #5. The gallery image above Label_ex.png was created with this code
    #   and finally, some labels were moved around in '''edit_mode'''.
    label (resi 200),"%1.3f" % b
    set label_font_id, 10
    set label_size, 10
    
    #6. This example shows how to label a selection with the 
    #   XYZ coordinates of the atoms 
    from pymol import stored
    stored.pos = []
    # select the carbon atoms in my hetero atoms to label
    select nn, het and e. C
    # get the XYZ coordinates and put them into stored.pos
    # insert at the front because pop() will read the array in reverse
    iterate_state 1, (nn), stored.pos.insert(0,(x,y,z))
    # label all N atoms.  You need the pop() function or else
    # PyMOL will complain b/c you didn't provide enough coords.
    label nn, ("%5.5s, %5.5s, %5.5s") %  stored.pos.pop()
    

# User Comments

## Labels Using ID Numbers

The following commnent, 
    
    
    label SELECTION, " %s" % ID
    

labels the SELECTION with atom ID numbers. 

You can make more complicated selections/lables such as 
    
    
    label SELECTION, " %s:%s %s" % (resi, resn, name)
    

which will give you something like "GLU:139 CG" 

## Labels Using One Letter Abbreviations

  * First, Add this to your $HOME/.pymolrc file:


    
    
    # start $HOME/.pymolrc modification
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}
    # end modification
    

  * Second, instead of:


    
    
    label n. ca, resn
    

use: 
    
    
    label n. ca, one_letter[resn]
    

or: ( to get something like D85) 
    
    
    label n. ca, "%s%s" %(one_letter[resn],resi)
    

## Labels and defer_builds_mode

If You have a weak video card, You might want to set 
    
    
    set defer_builds_mode, 5
    

It helps a lot but breaks labels rendering. You can use 
    
    
    set defer_builds_mode, 4
    

instead. 

# See Also

[Category:Labeling](/index.php/Category:Labeling "Category:Labeling")

All the settings posted above. 

Retrieved from "[https://pymolwiki.org/index.php?title=Label&oldid=12514](https://pymolwiki.org/index.php?title=Label&oldid=12514)"


---

## Talk:Label

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Fixes
  * 2 New Page Content
    * 2.1 New Page Overview
  * 3 Old Page
    * 3.1 DESCRIPTION
    * 3.2 USAGE
    * 3.3 SETTINGS
      * 3.3.1 FONT
        * 3.3.1.1 UTF8 Fonts
        * 3.3.1.2 Unicode Fonts
      * 3.3.2 SIZE
      * 3.3.3 COLOR
      * 3.3.4 EXPRESSION
      * 3.3.5 POSITION
    * 3.4 EXAMPLES
      * 3.4.1 Partial Charge
      * 3.4.2 Example 2
      * 3.4.3 More Advanced
    * 3.5 Users Comments
      * 3.5.1 Labels Using ID Numbers
      * 3.5.2 Labels Using One Letter Abbreviations



## Fixes

  * Updates needed



# New Page Content

## New Page Overview

This is the content for the new labels page. 

  


# Old Page

[![PyMol Labels](/images/7/78/Label_pre.png)](/index.php/File:Label_pre.png "PyMol Labels")

### DESCRIPTION

**label** allows one to configure the appearance of text labels for PyMOL objects. It labels one or more atoms properties over a selection using the python evaluator with a separate name space for each atom. The symbols defined in the name space are: 

  * **name** , the atom name
  * **resn** , the residue name
  * **resi** , the residue number/identifier
  * **chain** , the chain name
  * **q** ,
  * **b** , the occupancy/b-factor
  * **segi** , the segment identifier
  * **type** _(ATOM,HETATM)_ , the type of atom
  * **formal_charge** , the formal charge
  * **partial_charge** , the partial charge
  * **numeric_type** , the numeric type
  * **text_type** , the text type



All strings in the expression must be explicitly quoted. This operation typically takes several seconds per thousand atoms altered. To clear labels, simply omit the expression or set it to _._

[Label](/index.php/Label "Label") is great for labeling atoms, residues and objects. For a scene label, see [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). 

### USAGE
    
    
    label (selection),expression
    

  


### SETTINGS

#### FONT

There are 10 different scalable fonts. 
    
    
    set label_font_id, number
    

where number is 5 through 14. 

##### UTF8 Fonts

[![](/images/e/e0/New_fonts.jpeg)](/index.php/File:New_fonts.jpeg)

[](/index.php/File:New_fonts.jpeg "Enlarge")

New fonts in PyMol. Notice the alpha and beta characters.

Newer versions support UTF8 fonts; use **label_font_id** from above to 15 or 16. The good news about the UTF8 fonts is that they support the alpha and beta characters. (See image.) 

Here's some example code for the image at right: 
    
    
    # roman
    set label_font_id, 15
    set label_shadow_mode, 3
    label 5/CA, "\316\261-Helix"
    label 10/CA, "\316\262-Sheet"
    
    # italic
    set label_font_id, 16
    
    # make bigger
    set label_size, 50
    

##### Unicode Fonts

[![](/images/2/2f/Font_ex.png)](/index.php/File:Font_ex.png)

[](/index.php/File:Font_ex.png "Enlarge")

Notice the Angstrom and superscript 2 characters. You can add other characters as well.

PyMOL gives you the flexibility to use encoded unicode fonts. This allows us to insert various symbols, like the symbol used for Angstrom. Here are the steps to insert a character from the unicode character set. 

  * Find the code for your character at [Unicode Charts](http://www.unicode.org/charts). The Angstrom character,  \mathrm {\AA}  {\displaystyle {\textrm {\mathrm {\AA} }}} ![{\\displaystyle {\\textrm {\\mathrm {\\AA} }}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/fe5b45baabc495a7ef582330f075a45a1d6bb5fc) is u"\u00c5" and  ± {\displaystyle \pm } ![{\\displaystyle \\pm }](https://wikimedia.org/api/rest_v1/media/math/render/svg/869e366caf596564de4de06cb0ba124056d4064b) is u"\u00b1".
  * Label the selection. For simple strings, just type the string in double quote, -- "like this" -- and append to the end of that .encode('utf-8') -- "like this".encode('utf-8'). A working example is shown here,


    
    
    # label residue 30 with "4.1 Ang^2 +/- 0.65 Ang^2; see the image at right
    label i. 30, "4.1" + u"\u00c5\u00b2  \u00b1 0.65 \u00c5\u00b2 ".encode('utf-8')
    

#### SIZE

The font size can be adjusted 
    
    
    set label_size, number
    

where number is the point size (or -number for Angstroms) 

#### COLOR

Set a label's color by 
    
    
    set label_color, color
    

where color is a valid PyMol color. 

If the coloring of the labels is not _exactly_ the same as you'd expect (say black turns out grey, or red turns out pink), then try the following settings: 
    
    
    unset depth_cue
    unset ray_label_specular
    

#### EXPRESSION

To set what the label reads (see above) 
    
    
    label selection, expression
    

For example 
    
    
     label all, name
     label resi 10, b
    

#### POSITION

To position labels 
    
    
    edit_mode
    

then ctrl-middle-click-and-drag to position the label in space. (On Windows systems this appears to be shift-left-click-and-drag, presumably because those mice lack a true middle button.) 

ctrl-shift-left-click-and-drag alters a label's z-plane. (Windows only? This may use the middle button, rather than shift-left, under *NIX / 3-button mice systems.) 

### EXAMPLES

#### Partial Charge
    
    
    label (chain A),chain
    label (n;ca),"%s-%s" % (resn,resi)
    label (resi 200),"%1.3f" % partial_charge
    

#### Example 2

The following image was created with 
    
    
    label (resi 200),"%1.3f" % b
    set label_font_id, 10
    set label_size, 10
    

and finally, some labels were moved around in **edit_mode**. 

[![](/images/c/c2/Label_ex.png)](/index.php/File:Label_ex.png)

[](/index.php/File:Label_ex.png "Enlarge")

Labels.

  


#### More Advanced

This example shows how to label a selection with the XYZ coordinates of the atoms 
    
    
    from pymol import stored
    stored.pos = []
    
    # select the carbon atoms in my hetero atoms to label
    select nn, het and e. C
    
    # get the XYZ coordinates and put htem into stored.pos
    iterate_state 1, (nn), stored.pos.append((x,y,z))
    
    # label all N atoms.  You need the pop() function or else
    # PyMOL will complain b/c you didn't provide enough coords.
    label nn, ("%5.5s, %5.5s, %5.5s") %  stored.pos.pop()
    

### Users Comments

#### Labels Using ID Numbers

The following commnent, 
    
    
    label SELECTION, " %s" % ID 
    

labels the SELECTION with atom ID numbers. 

You can make more complicated selections/lables such as 
    
    
    label SELECTION, " %s:%s %s" % (resi, resn, name)
    

which will give you something like "GLU:139 CG" 

#### Labels Using One Letter Abbreviations

  * First, Add this to your $HOME/.pymolrc file:


    
    
    # start $HOME/.pymolrc modification
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}
    # end modification
    

  * . Second, instead of:


    
    
    label n. ca, resn
    

use: 
    
    
    label n. ca, one_letter[resn]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Talk:Label&oldid=12387](https://pymolwiki.org/index.php?title=Talk:Label&oldid=12387)"


---

## Label angle digits

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

label_angle_digits sets the precision of angle labels 

# Usage
    
    
    # set angle labels to display 2 decimals places to the right of the period
    set label_angle_digits, 2
    

# See Also

  * [label_distance_digits](/index.php/Label_distance_digits "Label distance digits")
  * [label_dihedral_digits](/index.php/Label_dihedral_digits "Label dihedral digits")



Retrieved from "[https://pymolwiki.org/index.php?title=Label_angle_digits&oldid=8422](https://pymolwiki.org/index.php?title=Label_angle_digits&oldid=8422)"


---

## Label color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/b/b0/Lce.png)](/index.php/File:Lce.png)

[](/index.php/File:Lce.png "Enlarge")

Label colors and sized changed, per object.

  


## Contents

  * 1 Overview
  * 2 Syntax
  * 3 User Comments
  * 4 See Also



## Overview

Sets the color that PyMol uses to draws/renders labels. This can be set for all objects/selections or for one in particular. 

## Syntax
    
    
    # set object's color to colorName
    set label_color, colorName, object
    
    # example showing two different objects
    # each with their own coloring.
    pseudoatom foo
    label foo, "foo"
    pseudoatom another
    label another, "Another label"
    set label_color, green, foo
    set label_color, lightpink, another
    translate [0, -10, 0], object=another
    set label_size, -2
    zoom foo or another, 10
    

# User Comments

If the coloring of the labels is not _exactly_ the same as you'd expect (say black turns out grey, or red turns out pink), then try the following settings: 
    
    
    unset depth_cue
    unset ray_label_specular
    

# See Also

[Color_Values](/index.php/Color_Values "Color Values"), [Category:Coloring](/index.php/Category:Coloring "Category:Coloring")

Retrieved from "[https://pymolwiki.org/index.php?title=Label_color&oldid=6913](https://pymolwiki.org/index.php?title=Label_color&oldid=6913)"


---

## Label dihedral digits

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

label_dihedral_digits sets the precision of dihedral labels 

# Usage
    
    
    # set dihedral labels to display 2 decimals places to the right of the period
    set label_dihedral_digits, 2
    

# See Also

  * [label_angle_digits](/index.php/Label_angle_digits "Label angle digits")
  * [label_distancel_digits](/index.php?title=Label_distancel_digits&action=edit&redlink=1 "Label distancel digits \(page does not exist\)")



Retrieved from "[https://pymolwiki.org/index.php?title=Label_dihedral_digits&oldid=8423](https://pymolwiki.org/index.php?title=Label_dihedral_digits&oldid=8423)"


---

## Label distance digits

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

label_distance_digits sets the precision of distance labels 

# Usage
    
    
    # set distance labels to display 2 decimals places to the right of the period
    set label_distance_digits, 2
    

# See Also

  * [label_angle_digits](/index.php/Label_angle_digits "Label angle digits")
  * [label_dihedral_digits](/index.php/Label_dihedral_digits "Label dihedral digits")



Retrieved from "[https://pymolwiki.org/index.php?title=Label_distance_digits&oldid=8420](https://pymolwiki.org/index.php?title=Label_distance_digits&oldid=8420)"


---

## Label font id

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The label_font_id setting sets the typeface (font family) for atom [labels](/index.php/Label "Label"). 

## Contents

  * 1 Syntax
  * 2 Available Font Families
  * 3 Special Characters
  * 4 Example 1
  * 5 Example 2
  * 6 See Also



## Syntax
    
    
    # use "Serif Bold" font family
    set label_font_id, 10
    

## Available Font Families

_Font IDs 0-4 were fixed size GLUT fonts, their support was dropped in PyMOL 1.6._

Name | label_font_id   
---|---  
Sans | 5   
Sans Oblique | 6   
Sans Bold | 7   
Sans Bold Oblique | 8   
Serif | 9   
Serif Oblique | 17   
Serif Bold | 10   
Serif Bold Oblique | 18   
Mono | 11   
Mono Oblique | 12   
Mono Bold | 13   
Mono Bold Oblique | 14   
Gentium Roman | 15   
Gentium Italic | 16   
  
## Special Characters

Several [Unicode](http://www.unicode.org/) characters are supported. They can be entered with [unicode literals](https://docs.python.org/2/howto/unicode.html#unicode-literals-in-python-source-code) as 4-digit hexadecimal escape sequences. 

**Example characters** (find the code for your character at [Unicode Charts](http://www.unicode.org/charts)): 

Code | Character | Name   
---|---|---  
`u"\u03b1"` | α | Alpha   
`u"\u03b2"` | β | Beta   
`u"\u00c5"` | Å | Ångström   
`u"\u00b1"` | ± | plus/minus   
`u"\u00b2"` | ² | superscript 2   
  
## Example 1

[![](/images/e/e0/New_fonts.jpeg)](/index.php/File:New_fonts.jpeg)

[](/index.php/File:New_fonts.jpeg "Enlarge")

The alpha and beta are Unicode characters.
    
    
    # if Python is configured with utf-8 default encoding (all incentive builds)
    label  5/CA, u"\u03b1-Helix"
    label 10/CA, u"\u03b2-Sheet"
    
    # if Python is configured differently, explicit utf-8 encoding is necessary
    label  5/CA, u"\u03b1-Helix".encode("utf-8")
    label 10/CA, u"\u03b2-Sheet".encode("utf-8")
    
    # italic
    set label_font_id, 16
    
    # make bigger
    set label_size, 50
    
    # cast shadows in ray tracing
    set label_shadow_mode, 3
    

## Example 2

[![](/images/2/2f/Font_ex.png)](/index.php/File:Font_ex.png)

[](/index.php/File:Font_ex.png "Enlarge")

Notice the Angstrom and superscript 2 characters. You can add other characters as well.
    
    
    # label residue 30 with "4.1 Ang^2 +/- 0.65 Ang^2
    label i. 30, "4.1" + u"\u00c5\u00b2  \u00b1 0.65 \u00c5\u00b2 "
    

## See Also

  * [label](/index.php/Label "Label")
  * [label_size](/index.php/Label_size "Label size")
  * [label_shadow_mode](/index.php/Label_shadow_mode "Label shadow mode")



Retrieved from "[https://pymolwiki.org/index.php?title=Label_font_id&oldid=12622](https://pymolwiki.org/index.php?title=Label_font_id&oldid=12622)"


---

## Label outline color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Sets which color PyMol uses to outlines the labels. This helps offset labels from their backgrounds. 

## Syntax
    
    
    set label_outline_color, red  # will draw red outline around the labels
    

Retrieved from "[https://pymolwiki.org/index.php?title=Label_outline_color&oldid=6222](https://pymolwiki.org/index.php?title=Label_outline_color&oldid=6222)"


---

## Label position

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

label_position sets the offset of the label relative to the atomic position. 

example: 
    
    
    set label_position,(3,2,1)
    

Offsets the labels 3 Å in X, 2 in Y and 1 in Z, relative to the [Viewport](/index.php/Viewport "Viewport"). Can be useful if [spheres](/index.php/Spheres "Spheres") needs [labeling](/index.php/Label "Label")

Retrieved from "[https://pymolwiki.org/index.php?title=Label_position&oldid=8315](https://pymolwiki.org/index.php?title=Label_position&oldid=8315)"


---

## Label shadow mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Sets whether or not PyMol renders shadows casted onto labels and whether or not the labels themselves cast shadows. 

Values 

  * 0, no shadows
  * 1, object shadows are projected onto the labels, but the labels don't cast shadows
  * 2, object shadows are projected onto the labels, AND the lables DO cast shadows
  * 3, object shadows are NOT projected onto the labels, BUT the lables DO cast shadows.



## Syntax
    
    
    set label_shadow_mode, 0  # turn off shadows projected onto lables, turn off shadows cast from labels
    set label_shadow_mode, 1  # turn on shadows projected onto labels, turn off shadows cast from labels
    set label_shadow_mode, 2  # turn off shadows projected onto lables, turn on shadows cast from labels
    set label_shadow_mode, 3  # turn on shadows projected onto lables, turn on shadows cast from labels
    

Retrieved from "[https://pymolwiki.org/index.php?title=Label_shadow_mode&oldid=4317](https://pymolwiki.org/index.php?title=Label_shadow_mode&oldid=4317)"


---

## Label size

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Sets how large the labels are rendered. You can use positive numbers 2, 3, 4, etc for point sizes, or negative numbers for Angstrom-based sizes. Default is 14 points. Labels in Angstrom-size scale with the distance from the front plane, labels in point-size don't. This automatic scaling works with the [draw](/index.php/Draw "Draw") command but not with the [ray](/index.php/Ray "Ray") command. 

## Syntax
    
    
    # set the label size to 10pt
    set label_size, 10
    
    # set the label size to 1.5 Ang. -- large!
    set label_size, -1.5
    

Retrieved from "[https://pymolwiki.org/index.php?title=Label_size&oldid=13561](https://pymolwiki.org/index.php?title=Label_size&oldid=13561)"


---

## Pseudoatom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

pseudoatom creates a molecular object with a pseudoatom or adds a pseudoatom to a molecular object if the specified object already exists. Default position is in the middle of the viewing window. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES FOR USE
  * 4 References



## USAGE
    
    
    pseudoatom object [, selection [, name [, resn [, resi [, chain
            [, segi [, elem [, vdw [, hetatm [, b [, q [, color [, label
            [, pos [, state [, mode [, quiet ]]]]]]]]]]]]]]]]]
    

## ARGUMENTS

  * object = string: name of object to create or modify


  * selection = string: optional atom selection. If given, calculate position (pos) as center of this selection (like [center of mass](/index.php/Center_Of_Mass "Center Of Mass") without atom masses).


  * name, resn, resi, chain, segi, elem, vdw, hetatm, b, q = string: [atom properties](/index.php/Property_Selectors "Property Selectors") of pseutoatom


  * color = string: [color](/index.php/Color "Color") of pseudoatom


  * label = string: place a [text label](/index.php/Label "Label") and hide all other representations


  * pos = 3-element tuple of floats: location in space (only if no selection given)


  * state = integer: [state](/index.php/State "State") to modify, 0 for current state, -1 for all states {default: 0}


  * mode = string: determines the vdw property if vdw is not given. Either as RMS distance (mode=rms, like [radius of gyration](/index.php/Radius_of_gyration "Radius of gyration") without atom masses) or as maximum distance (mode=extent) from pseudoatom to selection {default: rms}


    
    
    # create the pseudoatom
    pseudoatom tmpPoint
     ObjMol: created tmpPoint/PSDO/P/PSD`1/PS1
    # show it as a sphere.
    show spheres, tmpPoint
    
    # create another, with more options.
    pseudoatom tmpPoint2, resi=40, chain=ZZ, b=40, color=tv_blue, pos=[-10, 0, 10]
     ObjMol: created tmpPoint2/PSDO/ZZ/PSD`40/PS1
    

## EXAMPLES FOR USE

pseudoatom can be used for a wide variety of tasks where on must place an atom or a label in 3D space, e.g. as a placeholder for distance measurement or distance specifications. 
    
    
    # A pseudoatom as a placeholder for selections according to distance:
    load $TUT/1hpv.pdb
    pseudoatom tmp, pos=[10.0, 17.0, -3.0]
    show sticks, tmp expand 6
    delete tmp
    
    # A pseudoatom as placeholder for distance measurement: 
    # position it at the center of an aromatic ring.  Then 
    # calc the distance from another atom to the pseudoatom.
    load $TUT/1hpv.pdb
    pseudoatom pi_cent,b/53/cg+cz
    dist pi_cent////ps1, b/met`46/ce
    

You can use a pseudoatom to make a label for a scene title. Move the protein to the bottom of the window before the pseudoatom is created. Or move the label after creating it (Shift + Middle mouse button in editing mode). 
    
    
    fetch 1rq5
    pseudoatom forLabel
    label forLabel, "This Protein is a Carbohydrate Active Enzyme"
    set label_color, black
    # png ray=1
    

  * [![You can use a pseudoatom for a movable scene label.](/images/3/32/Patom.png)](/index.php/File:Patom.png "You can use a pseudoatom for a movable scene label.")

You can use a pseudoatom for a movable scene label. 

  * [![Example of using a pseudoatom for distance measurement.](/images/3/33/Pseu1.png)](/index.php/File:Pseu1.png "Example of using a pseudoatom for distance measurement.")

Example of using a pseudoatom for distance measurement. 




  


  


  


## References

PyMOL Mailing List 

Retrieved from "[https://pymolwiki.org/index.php?title=Pseudoatom&oldid=9090](https://pymolwiki.org/index.php?title=Pseudoatom&oldid=9090)"


---

