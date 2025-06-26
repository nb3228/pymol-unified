# Category: Publication Quality

## Advanced Coloring

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Coloring Molecules
    * 1.1 Basic Coloring
    * 1.2 Coloring secondary structures
    * 1.3 Coloring by atom type
    * 1.4 CMYK-safe Colors
    * 1.5 Coloring with 'chainbows' from a script
    * 1.6 Assign color by B-factor
  * 2 See Also
    * 2.1 Creating a Color bar
    * 2.2 Coloring insides and outsides of helices differently
    * 2.3 Coloring all objects differently
    * 2.4 List the color of atoms
    * 2.5 See also



# Coloring Molecules

## Basic Coloring

Any molecule in PyMOL can be assigned a color using the small rightmost buttons in the object list (in the upper right part of the main GUI window. The [Color](/index.php/Color "Color") command will do the same. 

PyMOL has a predefined set of colors that can be edited in the _Settings- >Colors_ menu. Alternatively, you can use the [Set_Color](/index.php/Set_Color "Set Color") command. 

  


## Coloring secondary structures

To assign helices, sheets and loops individual colors, do: 
    
    
    color red, ss h
    color yellow, ss s
    color green, ss l+''
    

When the colour bleeds from the ends of helices and sheets into loops, do: 
    
    
    set cartoon_discrete_colors, 1
    

Or activate _Cartoon - > Discrete Colors_ in the GUI menu. 

## Coloring by atom type

The util.cba* ("Color By Atom") commands color atoms according to type: oxygen in red, nitrogen in blue, hydrogen in white. Carbon will get a different colors, depending on the command: 

command  | carbon color   
---|---  
util.cba**g** | green   
util.cba**c** | cyan   
util.cba**m** | light magenta   
util.cba**y** | yellow   
util.cba**s** | salmon   
util.cba**w** | white/grey   
util.cba**b** | slate   
util.cba**o** | bright orange   
util.cba**p** | purple   
util.cba**k** | pink   
  
For instance: 
    
    
      util.cbay three
    

will color the object _three_ by atom type, with the carbon atoms in yellow. 

The util.cnc command will color all the atoms according to type, as in the util.cba* commands stated above, except for the C-atoms. 

For instance: 
    
    
      util.cnc three
    

will color the object _three_ by atom type, but leave the color of the C-atom unaltered. 

## CMYK-safe Colors

There are two distinct color spaces on computers: RGB (red-green-blue), which is for screens, and CMYK (cyan-magenta-yellow-black), which is for printing. Some RGB triplets do not have equivalents in CMYK space. As a result, a figure that looks great on a screen can come out with unpredictable colors when printed. 

Most applications do a good job with RGB-to-CMYK conversions for photos, but do not do such a good job with graphics that use pure primary colors. For example, reds are generally OK, but pure blues and greens do not translate very well. 

Here are some RGB values that are within the CMYK gamut (i.e. are "CMYK-safe"): 
    
    
    #optimized rgb values for cmyk output:
    set_color dblue= [0.05 , 0.19 , 0.57]
    set_color blue=  [0.02 , 0.50 , 0.72]
    set_color mblue= [0.5  , 0.7  , 0.9 ]
    set_color lblue= [0.86 , 1.00 , 1.00]
    
    set_color green= [0.00 , 0.53 , 0.22]
    set_color lgreen=[0.50 , 0.78 , 0.50]
    set_color yellow=[0.95 , 0.78 , 0.00]
    set_color orange=[1.00 , 0.40 , 0.0 ]
    
    # these are trivial
    set_color red=   [1.00 , 0.00 , 0.00]
    set_color mred=  [1.00 , 0.40 , 0.40]
    set_color lred=  [1.00 , 0.80 , 0.80]
    set_color vlred= [1.00 , 0.90 , 0.90]
    set_color white= [1.00 , 1.00 , 1.00]
    set_color vlgray=[0.95 , 0.95 , 0.95]
    set_color lgray= [0.90 , 0.90 , 0.90]
    set_color gray=  [0.70 , 0.70 , 0.70]
    set_color dgray= [0.50 , 0.50 , 0.50]
    set_color vdgray=[0.30 , 0.30 , 0.30]
    set_color black= [0.00 , 0.00 , 0.00]
    ##
    

Note that there are default atom colors such as "carbon", "nitrogen", "oxygen", "hydrogen", "sulfur", etc. which should also be redefined: 
    
    
    set_color carbon= [0.00 , 0.53 , 0.22]
    etc.
    

## Coloring with 'chainbows' from a script

The chainbow function can be invoked by: 
    
    
    util.chainbow("object-name")
    

  


## Assign color by B-factor

B-factor coloring can be done with the [spectrum](/index.php/Spectrum "Spectrum") command. Example: 
    
    
    spectrum b, blue_white_red, minimum=20, maximum=50
    as cartoon
    cartoon putty
    

# See Also

[Color](/index.php/Color "Color"), [Spectrum](/index.php/Spectrum "Spectrum")

## Creating a Color bar

To show a vertical/horizontal color bar indiacting the b-factor variation, use the script pseudobar.pml on the structure pseudobar.pdb, or do the following: 

  1. Create a pdb-file which contains CA positions only, whereas the numbers correspond to your wanted increments of colors. Be sure that CA's are separated by a contant value, say 5 Angstroem.
  2. Load this new pseudobar-pdb file into PyMOL, make bonds between increment 1 and increment 2 [increment 2 and increment 3 and so on...], define/assign a smooth color for each increment (copy colors definition from automatically created colors made by b-factor script) and show the b-factor bar as lines (or sticks).



Also, see the newly created [spectrumbar](/index.php/Spectrumbar "Spectrumbar") script! 

## Coloring insides and outsides of helices differently

The inside of helices can be adressed with: 
    
    
    set cartoon_highlight_color, red
    

  


## Coloring all objects differently

Is there a simple way to colour each object currently loaded, with a different colour? There is a script [color_obj.py](/index.php/Color_Objects "Color Objects") that does the job. 

USAGE 
    
    
           color_obj(rainbow=0)
    

This function colours each object currently in the PyMOL heirarchy with a different colour. Colours used are either the 22 named colours used by PyMOL (in which case the 23rd object, if it exists, gets the same colour as the first), or are the colours of the rainbow 

## List the color of atoms

To retrieve the color for all residues in a selection, you can iterate over it from the PyMOL command line 
    
    
    iterate all, print color
    

In Python, it looks like this: 
    
    
    import pymol
    pymol.color_list = []
    cmd.iterate('all', 'pymol.color_list.append(color)')
    print pymol.color_list
    

The colors listed will be in terms of Pymol indexing system, see [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices") for converting to names or rgb values. 

## See also

  * [Ramp_New](/index.php/Ramp_New "Ramp New")



Retrieved from "[https://pymolwiki.org/index.php?title=Advanced_Coloring&oldid=12310](https://pymolwiki.org/index.php?title=Advanced_Coloring&oldid=12310)"


---

## Antialias

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

The **antialias** setting in pymol controls the level of antialiasing that pymol does in the **[ray](/index.php/Ray "Ray")** -tracing routine. Higher numbers take longer to ray trace but provide much smoother and better looking images. 

## Examples

Notice the overall **quality** of the images improves as the antialiasing setting is increased. Click each image to see the high-resolution version; the small image doesn't show the differences in AA settings so well. 

  * [![AA set to 0](/images/e/e4/Aa0.png)](/index.php/File:Aa0.png "AA set to 0")

AA set to 0 

  * [![AA set to 1](/images/3/3a/Aa1.png)](/index.php/File:Aa1.png "AA set to 1")

AA set to 1 

  * [![AA set to 2](/images/e/e0/Aa2.png)](/index.php/File:Aa2.png "AA set to 2")

AA set to 2 




## Syntax
    
    
    set antialias,0  # low setting/off
    set antialias,2  # higher setting, better image quality
    

Retrieved from "[https://pymolwiki.org/index.php?title=Antialias&oldid=4275](https://pymolwiki.org/index.php?title=Antialias&oldid=4275)"


---

## CMYK

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Introduction
    * 1.1 Value
    * 1.2 CMYK vs. RGB
      * 1.2.1 Comment 1
      * 1.2.2 Response to Comment 1
      * 1.2.3 Response to Response to Comment 1
  * 2 Activating CMYK In PyMol



## Introduction

**CMYK** is a device-dependent colour space, i.e. the colour that is chosen to match to a particular RGB colour depends on the machine it is going to be printed on (and hence the inks that will be used in the printing). Therefore there is no one-to-one mapping of CMYK values to RGB values (even with a profile defining how the RGB values should be displayed). 

### Value

In the real world, scientists do not have time or resources to mess with color profiles and color-calibrated hardware. They simply want to go from screen to hardcopy knowing their colors will not be grossly distorted by their hardware, which is often consumer-grade, non-calibrated equipment. 

Though CMYK is indeed device-dependent, there are some color regions that, *in practice*, are more problematic than others when printing from an RGB source. If you stick to safer regions of RGB space that are pragmatically captured through PyMOL's "CMYK" space command, then you will get closer to a WYSIWYG experience even in the absense of calibration. 

On the other hand, blithely working in RGB space and then relying upon automatic RGB->CMYK color translations in the *absense* of color calibration for both display and printer almost always results in unacceptably poor color quality, and that is the practical real-world issue facing PyMOL users. 

### CMYK vs. RGB

There's an ongoing discussion about what's better for printing from PyMol: CMYK or getting your RGB to look just right. The merit of both are captured here: 

#### Comment 1

Since pymol's idea of CMYK is so limited (and from a printer's point of view probably shouldn't even exist) you are far better off getting the RGB image out of pymol to look just the way you want it to (assuming a calibrated display and approprate colour profile) and then use a full-featured converter to generate the appropriate CMYK image for the device it is going to be printed on. One suggestion for a converter is GIMP (soon to be renamed something else) because it does know how to use/respect ICC profiles. 

#### Response to Comment 1

Unless you have calibrated hardware for both display and printer, you are definitely _not far better off getting the RGB image out of pymol to look just the way you want it to_. Instead, you are far better off avoiding areas of RGB color space that are difficult or impossible to handle without professional-grade color hardware, and that is the sole task PyMOL's CMYK capability is designed to help you with. 

#### Response to Response to Comment 1

If somebody has gotten to the point where they are worrying about their own CMYK separations then telling them to just use a reduced gamut space to limit the conversion problems is less than helpful to my mind. They need to understand what is going on so they can make the decisions that are necessary if they are to get the best colour reproduction on the journal pages. 

## Activating CMYK In PyMol

Select 
    
    
    Display -> Color Space -> CMYK
    

Please read the above note about color spaces, too. 

Retrieved from "[https://pymolwiki.org/index.php?title=CMYK&oldid=8466](https://pymolwiki.org/index.php?title=CMYK&oldid=8466)"


---

## Color Values

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Background
  * 2 Spectral range colours
  * 3 Color by Element/Atom
  * 4 Simple named colours
  * 5 Chemical element colours
  * 6 Interactive menu colours
    * 6.1 reds
    * 6.2 greens
    * 6.3 blues
    * 6.4 yellows
    * 6.5 magentas
    * 6.6 cyans
    * 6.7 oranges
    * 6.8 tints
    * 6.9 grays
  * 7 Secondary structure colour schemes
  * 8 Carbon colour schemes



## Background

  * Naming 
    * All listed colours can be specified by name (eg color red, mymolecule)
    * American and English spellings of gray/grey and color/colour can be used
  * Sources 
    * Pymol colours are specified as triples of independent Red, Green and Blue contribution in range 0-1.0 in square brackets (e.g. [1.0,1.0,1.0] for white)
    * Colour specs from source files Color.c, menu.py and appearance.py of Pymol 0.99 beta29
    * Colour names can be defined/redefined using set_color command or interactively in Colors item within Setting menu



## Spectral range colours

Spectral ranges are available by name including numerical value (eg grey56, s532) 

  * Grays from 'gray00' to 'gray99' (white to black)
  * Spectrum from 's000' to 's999' (violet to red)
  * Original spectrum with extra blue and red at ends from 'o000' to 'o999'
  * Reversed offset spectrum from 'r000' to 'r999'
  * Complementary spectrum from 'c000' to 'c999'
  * Complementary spectrum separated by white from 'w000' to 'w999'



See Colour.c in source code for details about how spectra are calculated 

## Color by Element/Atom

[ See: Color by Atom](/index.php/Advanced_Coloring#Coloring_by_atom_type "Advanced Coloring")

## Simple named colours

| name | R | G | B | note   
---|---|---|---|---|---  
| aquamarine | 0.5 | 1.0 | 1.0 |   
| black | 0.0 | 0.0 | 0.0 |   
| blue | 0.0 | 0.0 | 1.0 |   
| bluewhite | 0.85 | 0.85 | 1.00 |   
| br0 | 0.1 | 0.1 | 1.0 |   
| br1 | 0.2 | 0.1 | 0.9 |   
| br2 | 0.3 | 0.1 | 0.8 |   
| br3 | 0.4 | 0.1 | 0.7 |   
| br4 | 0.5 | 0.1 | 0.6 |   
| br5 | 0.6 | 0.1 | 0.5 |   
| br6 | 0.7 | 0.1 | 0.4 |   
| br7 | 0.8 | 0.1 | 0.3 |   
| br8 | 0.9 | 0.1 | 0.2 |   
| br9 | 1.0 | 0.1 | 0.1 |   
| brightorange | 1.0 | 0.7 | 0.2 |   
| brown | 0.65 | 0.32 | 0.17 |   
| carbon | 0.2 | 1.0 | 0.2 |   
| chartreuse | 0.5 | 1.0 | 0.0 |   
| chocolate | 0.555 | 0.222 | 0.111 |   
| cyan | 0.0 | 1.0 | 1.0 |   
| darksalmon | 0.73 | 0.55 | 0.52 |   
| dash | 1.0 | 1.0 | 0.0 |   
| deepblue | 0.25 | 0.25 | 0.65 | deep   
| deepolive | 0.6 | 0.6 | 0.1 |   
| deeppurple | 0.6 | 0.1 | 0.6 |   
| deepsalmon | 1.0 | 0.5 | 0.5 | duplicated?   
| deepsalmon | 1.00 | 0.42 | 0.42 | duplicated?   
| deepteal | 0.1 | 0.6 | 0.6 |   
| density | 0.1 | 0.1 | 0.6 |   
| dirtyviolet | 0.70 | 0.50 | 0.50 |   
| firebrick | 0.698 | 0.13 | 0.13 |   
| forest | 0.2 | 0.6 | 0.2 |   
| gray | 0.5 | 0.5 | 0.5 | american spelling   
| green | 0.0 | 1.0 | 0.0 |   
| greencyan | 0.25 | 1.00 | 0.75 |   
| grey | 0.5 | 0.5 | 0.5 | english spelling   
| hotpink | 1.0 | 0.0 | 0.5 |   
| hydrogen | 0.9 | 0.9 | 0.9 |   
| lightblue | 0.75 | 0.75 | 1.0 |   
| lightmagenta | 1.0 | 0.2 | 0.8 |   
| lightorange | 1.0 | 0.8 | 0.5 |   
| lightpink | 1.00 | 0.75 | 0.87 |   
| lightteal | 0.4 | 0.7 | 0.7 |   
| lime | 0.5 | 1.0 | 0.5 |   
| limegreen | 0.0 | 1.0 | 0.5 |   
| limon | 0.75 | 1.00 | 0.25 |   
| magenta | 1.0 | 0.0 | 1.0 |   
| marine | 0.0 | 0.5 | 1.0 |   
| nitrogen | 0.2 | 0.2 | 1.0 |   
| olive | 0.77 | 0.70 | 0.00 |   
| orange | 1.0 | 0.5 | 0.0 |   
| oxygen | 1.0 | 0.3 | 0.3 |   
| palecyan | 0.8 | 1.0 | 1.0 |   
| palegreen | 0.65 | 0.9 | 0.65 |   
| paleyellow | 1.0 | 1.0 | 0.5 |   
| pink | 1.0 | 0.65 | 0.85 |   
| purple | 0.75 | 0.00 | 0.75 |   
| purpleblue | 0.5 | 0.0 | 1.0 | legacy name   
| raspberry | 0.70 | 0.30 | 0.40 |   
| red | 1.0 | 0.0 | 0.0 |   
| ruby | 0.6 | 0.2 | 0.2 |   
| salmon | 1.0 | 0.6 | 0.6 | was 0.5   
| sand | 0.72 | 0.55 | 0.30 |   
| skyblue | 0.20 | 0.50 | 0.80 |   
| slate | 0.5 | 0.5 | 1.0 |   
| smudge | 0.55 | 0.70 | 0.40 |   
| splitpea | 0.52 | 0.75 | 0.00 |   
| sulfur | 0.9 | 0.775 | 0.25 | far enough from yellow   
| teal | 0.00 | 0.75 | 0.75 |   
| tv_blue | 0.3 | 0.3 | 1.0 |   
| tv_green | 0.2 | 1.0 | 0.2 |   
| tv_orange | 1.0 | 0.55 | 0.15 |   
| tv_red | 1.0 | 0.2 | 0.2 |   
| tv_yellow | 1.0 | 1.0 | 0.2 |   
| violet | 1.0 | 0.5 | 1.0 |   
| violetpurple | 0.55 | 0.25 | 0.60 |   
| warmpink | 0.85 | 0.20 | 0.50 |   
| wheat | 0.99 | 0.82 | 0.65 |   
| white | 1.0 | 1.0 | 1.0 |   
| yellow | 1.0 | 1.0 | 0.0 |   
| yelloworange | 1.0 | 0.87 | 0.37 |   
  
## Chemical element colours

| name | R | G | B | note   
---|---|---|---|---|---  
| actinium | 0.439215686 | 0.670588235 | 0.980392157 |   
| aluminum | 0.749019608 | 0.650980392 | 0.650980392 |   
| americium | 0.329411765 | 0.360784314 | 0.949019608 |   
| antimony | 0.619607843 | 0.388235294 | 0.709803922 |   
| argon | 0.501960784 | 0.819607843 | 0.890196078 |   
| arsenic | 0.741176471 | 0.501960784 | 0.890196078 |   
| astatine | 0.458823529 | 0.309803922 | 0.270588235 |   
| barium | 0.000000000 | 0.788235294 | 0.000000000 |   
| berkelium | 0.541176471 | 0.309803922 | 0.890196078 |   
| beryllium | 0.760784314 | 1.000000000 | 0.000000000 |   
| bismuth | 0.619607843 | 0.309803922 | 0.709803922 |   
| bohrium | 0.878431373 | 0.000000000 | 0.219607843 |   
| boron | 1.000000000 | 0.709803922 | 0.709803922 |   
| bromine | 0.650980392 | 0.160784314 | 0.160784314 |   
| cadmium | 1.000000000 | 0.850980392 | 0.560784314 |   
| calcium | 0.239215686 | 1.000000000 | 0.000000000 |   
| californium | 0.631372549 | 0.211764706 | 0.831372549 |   
| carbon | 0.2 | 1.0 | 0.2 |   
| cerium | 1.000000000 | 1.000000000 | 0.780392157 |   
| cesium | 0.341176471 | 0.090196078 | 0.560784314 |   
| chlorine | 0.121568627 | 0.941176471 | 0.121568627 |   
| chromium | 0.541176471 | 0.600000000 | 0.780392157 |   
| cobalt | 0.941176471 | 0.564705882 | 0.627450980 |   
| copper | 0.784313725 | 0.501960784 | 0.200000000 |   
| curium | 0.470588235 | 0.360784314 | 0.890196078 |   
| deuterium | 0.9 | 0.9 | 0.9 |   
| dubnium | 0.819607843 | 0.000000000 | 0.309803922 |   
| dysprosium | 0.121568627 | 1.000000000 | 0.780392157 |   
| einsteinium | 0.701960784 | 0.121568627 | 0.831372549 |   
| erbium | 0.000000000 | 0.901960784 | 0.458823529 |   
| europium | 0.380392157 | 1.000000000 | 0.780392157 |   
| fermium | 0.701960784 | 0.121568627 | 0.729411765 |   
| fluorine | 0.701960784 | 1.000000000 | 1.000000000 |   
| francium | 0.258823529 | 0.000000000 | 0.400000000 |   
| gadolinium | 0.270588235 | 1.000000000 | 0.780392157 |   
| gallium | 0.760784314 | 0.560784314 | 0.560784314 |   
| germanium | 0.400000000 | 0.560784314 | 0.560784314 |   
| gold | 1.000000000 | 0.819607843 | 0.137254902 |   
| hafnium | 0.301960784 | 0.760784314 | 1.000000000 |   
| hassium | 0.901960784 | 0.000000000 | 0.180392157 |   
| helium | 0.850980392 | 1.000000000 | 1.000000000 |   
| holmium | 0.000000000 | 1.000000000 | 0.611764706 |   
| hydrogen | 0.9 | 0.9 | 0.9 |   
| indium | 0.650980392 | 0.458823529 | 0.450980392 |   
| iodine | 0.580392157 | 0.000000000 | 0.580392157 |   
| iridium | 0.090196078 | 0.329411765 | 0.529411765 |   
| iron | 0.878431373 | 0.400000000 | 0.200000000 |   
| krypton | 0.360784314 | 0.721568627 | 0.819607843 |   
| lanthanum | 0.439215686 | 0.831372549 | 1.000000000 |   
| lawrencium | 0.780392157 | 0.000000000 | 0.400000000 |   
| lead | 0.341176471 | 0.349019608 | 0.380392157 |   
| lithium | 0.800000000 | 0.501960784 | 1.000000000 |   
| lutetium | 0.000000000 | 0.670588235 | 0.141176471 |   
| magnesium | 0.541176471 | 1.000000000 | 0.000000000 |   
| manganese | 0.611764706 | 0.478431373 | 0.780392157 |   
| meitnerium | 0.921568627 | 0.000000000 | 0.149019608 |   
| mendelevium | 0.701960784 | 0.050980392 | 0.650980392 |   
| mercury | 0.721568627 | 0.721568627 | 0.815686275 |   
| molybdenum | 0.329411765 | 0.709803922 | 0.709803922 |   
| neodymium | 0.780392157 | 1.000000000 | 0.780392157 |   
| neon | 0.701960784 | 0.890196078 | 0.960784314 |   
| neptunium | 0.000000000 | 0.501960784 | 1.000000000 |   
| nickel | 0.313725490 | 0.815686275 | 0.313725490 |   
| niobium | 0.450980392 | 0.760784314 | 0.788235294 |   
| nitrogen | 0.2 | 0.2 | 1.0 |   
| nobelium | 0.741176471 | 0.050980392 | 0.529411765 |   
| osmium | 0.149019608 | 0.400000000 | 0.588235294 |   
| oxygen | 1.0 | 0.3 | 0.3 |   
| palladium | 0.000000000 | 0.411764706 | 0.521568627 |   
| phosphorus | 1.000000000 | 0.501960784 | 0.000000000 |   
| platinum | 0.815686275 | 0.815686275 | 0.878431373 |   
| plutonium | 0.000000000 | 0.419607843 | 1.000000000 |   
| polonium | 0.670588235 | 0.360784314 | 0.000000000 |   
| potassium | 0.560784314 | 0.250980392 | 0.831372549 |   
| praseodymium | 0.850980392 | 1.000000000 | 0.780392157 |   
| promethium | 0.639215686 | 1.000000000 | 0.780392157 |   
| protactinium | 0.000000000 | 0.631372549 | 1.000000000 |   
| radium | 0.000000000 | 0.490196078 | 0.000000000 |   
| radon | 0.258823529 | 0.509803922 | 0.588235294 |   
| rhenium | 0.149019608 | 0.490196078 | 0.670588235 |   
| rhodium | 0.039215686 | 0.490196078 | 0.549019608 |   
| rubidium | 0.439215686 | 0.180392157 | 0.690196078 |   
| ruthenium | 0.141176471 | 0.560784314 | 0.560784314 |   
| rutherfordium | 0.800000000 | 0.000000000 | 0.349019608 |   
| samarium | 0.560784314 | 1.000000000 | 0.780392157 |   
| scandium | 0.901960784 | 0.901960784 | 0.901960784 |   
| seaborgium | 0.850980392 | 0.000000000 | 0.270588235 |   
| selenium | 1.000000000 | 0.631372549 | 0.000000000 |   
| silicon | 0.941176471 | 0.784313725 | 0.627450980 |   
| silver | 0.752941176 | 0.752941176 | 0.752941176 |   
| sodium | 0.670588235 | 0.360784314 | 0.949019608 |   
| strontium | 0.000000000 | 1.000000000 | 0.000000000 |   
| sulfur | 0.9 | 0.775 | 0.25 |   
| tantalum | 0.301960784 | 0.650980392 | 1.000000000 |   
| technetium | 0.231372549 | 0.619607843 | 0.619607843 |   
| tellurium | 0.831372549 | 0.478431373 | 0.000000000 |   
| terbium | 0.188235294 | 1.000000000 | 0.780392157 |   
| thallium | 0.650980392 | 0.329411765 | 0.301960784 |   
| thorium | 0.000000000 | 0.729411765 | 1.000000000 |   
| thulium | 0.000000000 | 0.831372549 | 0.321568627 |   
| tin | 0.400000000 | 0.501960784 | 0.501960784 |   
| titanium | 0.749019608 | 0.760784314 | 0.780392157 |   
| tungsten | 0.129411765 | 0.580392157 | 0.839215686 |   
| uranium | 0.000000000 | 0.560784314 | 1.000000000 |   
| vanadium | 0.650980392 | 0.650980392 | 0.670588235 |   
| xenon | 0.258823529 | 0.619607843 | 0.690196078 |   
| ytterbium | 0.000000000 | 0.749019608 | 0.219607843 |   
| yttrium | 0.580392157 | 1.000000000 | 1.000000000 |   
| zinc | 0.490196078 | 0.501960784 | 0.690196078 |   
| zirconium | 0.580392157 | 0.878431373 | 0.878431373 |   
  
## Interactive menu colours

Accessible from Colour submenu of objects 

### reds

| name | R | G | B |   
---|---|---|---|---|---  
| red | 1.0 | 0.0 | 0.0 |   
| tv_red | 1.0 | 0.2 | 0.2 |   
| raspberry | 0.70 | 0.30 | 0.40 |   
| darksalmon | 0.73 | 0.55 | 0.52 |   
| salmon | 1.0 | 0.6 | 0.6 |   
| deepsalmon | 1.00 | 0.42 | 0.42 |   
| warmpink | 0.85 | 0.20 | 0.50 |   
| firebrick | 0.698 | 0.13 | 0.13 |   
| ruby | 0.6 | 0.2 | 0.2 |   
| chocolate | 0.555 | 0.222 | 0.111 |   
| brown | 0.65 | 0.32 | 0.17 |   
  
### greens

| name | R | G | B |   
---|---|---|---|---|---  
| green | 0.0 | 1.0 | 0.0 |   
| tv_green | 0.2 | 1.0 | 0.2 |   
| chartreuse | 0.5 | 1.0 | 0.0 |   
| splitpea | 0.52 | 0.75 | 0.00 |   
| smudge | 0.55 | 0.70 | 0.40 |   
| palegreen | 0.65 | 0.9 | 0.65 |   
| limegreen | 0.0 | 1.0 | 0.5 |   
| lime | 0.5 | 1.0 | 0.5 |   
| limon | 0.75 | 1.00 | 0.25 |   
| forest | 0.2 | 0.6 | 0.2 |   
  
### blues

| name | R | G | B |   
---|---|---|---|---|---  
| blue | 0.0 | 0.0 | 1.0 |   
| tv_blue | 0.3 | 0.3 | 1.0 |   
| marine | 0.0 | 0.5 | 1.0 |   
| slate | 0.5 | 0.5 | 1.0 |   
| lightblue | 0.75 | 0.75 | 1.0 |   
| skyblue | 0.20 | 0.50 | 0.80 |   
| purpleblue | 0.5 | 0.0 | 1.0 |   
| deepblue | 0.25 | 0.25 | 0.65 |   
| density | 0.1 | 0.1 | 0.6 |   
  
### yellows

| name | R | G | B |   
---|---|---|---|---|---  
| yellow | 1.0 | 1.0 | 0.0 |   
| tv_yellow | 1.0 | 1.0 | 0.2 |   
| paleyellow | 1.0 | 1.0 | 0.5 |   
| yelloworange | 1.0 | 0.87 | 0.37 |   
| limon | 0.75 | 1.00 | 0.25 |   
| wheat | 0.99 | 0.82 | 0.65 |   
| sand | 0.72 | 0.55 | 0.30 |   
  
### magentas

| name | R | G | B |   
---|---|---|---|---|---  
| magenta | 1.0 | 0.0 | 1.0 |   
| lightmagenta | 1.0 | 0.2 | 0.8 |   
| hotpink | 1.0 | 0.0 | 0.5 |   
| pink | 1.0 | 0.65 | 0.85 |   
| lightpink | 1.00 | 0.75 | 0.87 |   
| dirtyviolet | 0.70 | 0.50 | 0.50 |   
| violet | 1.0 | 0.5 | 1.0 |   
| violetpurple | 0.55 | 0.25 | 0.60 |   
| purple | 0.75 | 0.00 | 0.75 |   
| deeppurple | 0.6 | 0.1 | 0.6 |   
  
### cyans

| name | R | G | B |   
---|---|---|---|---|---  
| cyan | 0.0 | 1.0 | 1.0 |   
| palecyan | 0.8 | 1.0 | 1.0 |   
| aquamarine | 0.5 | 1.0 | 1.0 |   
| greencyan | 0.25 | 1.00 | 0.75 |   
| teal | 0.00 | 0.75 | 0.75 |   
| deepteal | 0.1 | 0.6 | 0.6 |   
| lightteal | 0.4 | 0.7 | 0.7 |   
  
### oranges

| name | R | G | B |   
---|---|---|---|---|---  
| orange | 1.0 | 0.5 | 0.0 |   
| tv_orange | 1.0 | 0.55 | 0.15 |   
| brightorange | 1.0 | 0.7 | 0.2 |   
| lightorange | 1.0 | 0.8 | 0.5 |   
| yelloworange | 1.0 | 0.87 | 0.37 |   
| olive | 0.77 | 0.70 | 0.00 |   
| deepolive | 0.6 | 0.6 | 0.1 |   
  
### tints

| name | R | G | B |   
---|---|---|---|---|---  
| wheat | 0.99 | 0.82 | 0.65 |   
| palegreen | 0.65 | 0.9 | 0.65 |   
| lightblue | 0.75 | 0.75 | 1.0 |   
| paleyellow | 1.0 | 1.0 | 0.5 |   
| lightpink | 1.00 | 0.75 | 0.87 |   
| palecyan | 0.8 | 1.0 | 1.0 |   
| lightorange | 1.0 | 0.8 | 0.5 |   
| bluewhite | 0.85 | 0.85 | 1.00 |   
  
### grays

| name | R | G | B |   
---|---|---|---|---|---  
| white | 1.0 | 1.0 | 1.0 |   
| grey90 | 0.9 | 0.9 | 0.9 |   
| grey80 | 0.8 | 0.8 | 0.8 |   
| grey70 | 0.7 | 0.7 | 0.7 |   
| grey60 | 0.6 | 0.6 | 0.6 |   
| grey50 | 0.5 | 0.5 | 0.5 |   
| grey40 | 0.4 | 0.4 | 0.4 |   
| grey30 | 0.3 | 0.3 | 0.3 |   
| grey20 | 0.2 | 0.2 | 0.2 |   
| grey10 | 0.1 | 0.1 | 0.1 |   
| black | 0.0 | 0.0 | 0.0 |   
  
## Secondary structure colour schemes

Accessible from 'by ss' section within Colour submenu of object 

helix | sheet | loop   
---|---|---  
red | yellow | green   
cyan | magenta | salmon   
cyan | red | magenta   
  
## Carbon colour schemes

  * Accessible from 'by element' section within Colour submenu of object
  * Only carbon colour is affected, other elements are set to their default element colour
  * Separated into five sets: 
    * main set

| name   
---|---  
| (carbon not changed)   
| carbon   
| cyan   
| lightmagenta   
| yellow   
| salmon   
| hydrogen   
| slate   
| orange   
  
  *     * set2

| name   
---|---  
| lime   
| deepteal   
| hotpink   
| yelloworange   
| violetpurple   
| grey70   
| marine   
| olive   
  
  *     * set3

| name   
---|---  
| smudge   
| teal   
| dirtyviolet   
| wheat   
| deepsalmon   
| lightpink   
| aquamarine   
| paleyellow   
  
  *     * set4

| name   
---|---  
| limegreen   
| skyblue   
| warmpink   
| limon   
| violet   
| bluewhite   
| greencyan   
| sand   
  
  *     * set 5

| name   
---|---  
| forest   
| lightteal   
| darksalmon   
| splitpea   
| raspberry   
| grey50   
| deepblue   
| brown   
  
Retrieved from "[https://pymolwiki.org/index.php?title=Color_Values&oldid=8384](https://pymolwiki.org/index.php?title=Color_Values&oldid=8384)"


---

## Draw

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/4/4a/Draw_ex.png)](/index.php/File:Draw_ex.png)

[](/index.php/File:Draw_ex.png "Enlarge")

Showing Draw in Action

**Draw** creates an oversized and antialiased OpenGL image using the current window. It's like [Ray](/index.php/Ray "Ray") but not ray traced. Also, as now with [Ray](/index.php/Ray "Ray") the oversized images are scaled and shown in the viewer window. As Draw doesn't ray trace the shadows of the scene, it is **far** faster than ray. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Match ray tracing appearance
  * 4 Notes
  * 5 See Also



## Usage
    
    
    draw [ width [, height [, antialias ]]]
    

## Examples
    
    
    draw 1600
    

will create an 1600-pixel wide image with an aspect ratio equal to that of the current screen. 
    
    
    draw 2000, 1500, 0
    

will create a 2000 by 1500 pixel image with antialiasing disabled 
    
    
    draw 600, 400, 2
    

will create a 600 by 500 pixel image with maximum (16X) antialiasing 

## Match ray tracing appearance

Since PyMOL 1.6, all "line"-type representations can be rendered as cylinders if [shaders](/index.php?title=Use_shaders&action=edit&redlink=1 "Use shaders \(page does not exist\)") are available and all ***_as_cylinders** settings are set. Example: 
    
    
    set line_as_cylinders
    set nonbonded_as_cylinders
    set ribbon_as_cylinders
    set mesh_as_cylinders
    set dash_as_cylinders
    set render_as_cylinders
    draw 3000, 2000
    

## Notes

This is a quick alternative to ray tracing with antialiasing enabled by default. 

Note that image size and antialiasing may be limited in some cases due to OpenGL hardware limits, such as screen size. For example, high-end ATI and NVidia cards max out at 4096 x 4096. 

## See Also

  * [ray](/index.php/Ray "Ray")
  * [png](/index.php/Png "Png")



Retrieved from "[https://pymolwiki.org/index.php?title=Draw&oldid=12430](https://pymolwiki.org/index.php?title=Draw&oldid=12430)"


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

## Publication Quality Images

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 High Resolution (DPI) Images
    * 1.1 Creating Resolution Specific Images
      * 1.1.1 GIMP 2.0
  * 2 Important Note
    * 2.1 See Also



## High Resolution (DPI) Images

### Creating Resolution Specific Images

We often get the question: 
    
    
    Q.  How do I create a X DPI image using PyMol's 'ray' command?
    

The answer is deceivingly simple. There are two steps in the process. First, setup your scene and use [ray](/index.php/Cmd_ray "Cmd ray") like this 
    
    
    ray 2400, 2400
    

From this you will get a large square image. You can now resize it to X inches where X=A*B, A and B factors of 2400 (eg: 8 inch photo @ 300DPI, 4 inch photo at 600DPI, 2 inch photo @ 1200DPI). The formula for creating images is: 
    
    
    rayVal = inches * DPI
    

where rayVal is the value you pass to "ray". Thus, an 8 inch square photo at 72, 100, and 300 DPI would be created by the following commands: 
    
    
    ray 576,576   # 8inch * 72dpi
    ray 800,800   # 8inch * 100dpi; or a 4inch * 200 DPI photo; or 1x800.
    ray 2400,2400 # 8inch * 300dpi; 6"x400dpi, etc...
    

  
The **new and better solution** is to use the [Png](/index.php/Png "Png") command with the **dpi** setting, such as 
    
    
    png fileName, dpi=300
    

If you don't have a new enough version of PyMOL that supports this, then the **old solution** is to use [PhotoShop](/index.php?title=PhotoShop&action=edit&redlink=1 "PhotoShop \(page does not exist\)") or [GIMP](/index.php?title=GIMP&action=edit&redlink=1 "GIMP \(page does not exist\)") to resize the photo to the appropriate height and width while making sure also to tell the [PhotoShop](/index.php?title=PhotoShop&action=edit&redlink=1 "PhotoShop \(page does not exist\)") or [GIMP](/index.php?title=GIMP&action=edit&redlink=1 "GIMP \(page does not exist\)") what resolution the image is. 

  


##### GIMP 2.0

In Gimp, load your image (it'll most likely be very large) then 

  * Select:


    
    
    Image ->
    Scale Image ->
    'Print Size & Scale Image Section' change 'New Width' to whatever width you decided on before.
    

  * Save Image (File->Save)



You should be set. 

  


## Important Note

For those who are hung up on the idea that DPI is important please read this interesting [explanation](http://www.rideau-info.com/photos/mythdpi.html). 

### See Also

  * [Png](/index.php/Png "Png")
  * [Ray](/index.php/Ray "Ray")
  * [Draw](/index.php/Draw "Draw")
  * [CMYK](/index.php/CMYK "CMYK")
  * the rest of the pages in this [Category](/index.php/Category:Publication_Quality "Category:Publication Quality")



Retrieved from "[https://pymolwiki.org/index.php?title=Publication_Quality_Images&oldid=6490](https://pymolwiki.org/index.php?title=Publication_Quality_Images&oldid=6490)"


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

## Ray opaque background

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting changes how PyMol treats the background. If this option is ON then the background is whatever you specify -- like black or white; however, if the setting is OFF, then the background will be treated as a transparent alpha channel. 

## Settings

Note: turning this setting **off** creates the transparent backgrounds. 
    
    
    # turn on transparent alpha channel
    set ray_opaque_background, off
    

## Examples

  * [![ray_opaque_background set to 0. The image from PyMOL is shown over a checkerboard.](/images/e/e7/Rob0.png)](/index.php/File:Rob0.png "ray_opaque_background set to 0. The image from PyMOL is shown over a checkerboard.")

ray_opaque_background set to 0. The image from PyMOL is shown over a checkerboard. 

  * [![ray_opaque_background set to 1. The image from PyMOL is shown over a checkerboard, however because the background is opaque, the checkerboard underneath does not show through.](/images/a/ab/Rob1.png)](/index.php/File:Rob1.png "ray_opaque_background set to 1. The image from PyMOL is shown over a checkerboard, however because the background is opaque, the checkerboard underneath does not show through.")

ray_opaque_background set to 1. The image from PyMOL is shown over a checkerboard, however because the background is opaque, the checkerboard underneath does not show through. 




Retrieved from "[https://pymolwiki.org/index.php?title=Ray_opaque_background&oldid=11598](https://pymolwiki.org/index.php?title=Ray_opaque_background&oldid=11598)"


---

## Ray orthoscopic

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting controls whether ray-traced images are rendered with or without perspective. Note that this can be in conflict with the setting "orthoscopic"; by default, images are rendered with the same orthoscopic setting as the viewport, unless "ray_orthoscopic" is deliberately set otherwise. 

## Settings
    
    
    set ray_orthoscopic, off   # render ray-traced images with perspective (2-4x SLOWER)
    set ray_orthoscopic, on    # render ray-traced images without perspective
    

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




Retrieved from "[https://pymolwiki.org/index.php?title=Ray_orthoscopic&oldid=8551](https://pymolwiki.org/index.php?title=Ray_orthoscopic&oldid=8551)"


---

## Space

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Space selects a color palette (or color space). Current choices are 'RGB', 'CMYK' (usually used in printing) and, 'pymol'. 

Whereas computer displays use the RGB color space, computer printers typically use the CMYK color space. The two spaces are non-equivalent, meaning that certain RGB colors cannot be expressed in the CMYK space and vice-versa. And as a result, molecular graphics images prepared using RGB often turn out poorly when converted to CMYK, with purplish blues or yellowish greens. "space cmyk" forces PyMOL to restrict its use of the RGB color space to subset that can be reliably converted to CMYK using common tools such as Adobe Photoshop. Thus, what you see on the screen is much closer to what you will get in print. 

Analog video systems as well as digital video compression codecs based on the YUV color space also have incompatibilities with RGB. Oversaturated colors usually cause the most problems. 

Although PyMOL lacks "space yuv", "space pymol" will help PyMOL avoid oversaturated colors can cause problems when exporting animations to video. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



# USAGE
    
    
    space space
    

**space**

    

    the color space from which to choose. Options are, 'rgb', 'cmyk' or 'pymol'. The default is 'RGB'.

# EXAMPLES
    
    
    space rgb
    space cmyk
    space pymol
    

  


# PYMOL API
    
    
    cmd.space(string space)
    

# SEE ALSO

[Color](/index.php/Color "Color"), [Coloring Pages](/index.php/Category:Coloring "Category:Coloring"), [Publication Quality Pages](/index.php/Category:Publication_Quality "Category:Publication Quality")

Retrieved from "[https://pymolwiki.org/index.php?title=Space&oldid=7565](https://pymolwiki.org/index.php?title=Space&oldid=7565)"


---

## Stereo Figures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The stereo command sets the current stereo mode. Stereo mode is a convenient way to "see" 3D with two images from slightly different angles. 

_There are corresponding**stereo** and [stereo_mode](/index.php/Stereo_mode "Stereo mode") settings which are controlled by the **stereo** command, so you don't need to set them directly._

## Usage
    
    
    stereo [ toggle ]
    

Valid values for the **toggle** argument are: on, swap, off, quadbuffer, crosseye, walleye, geowall, sidebyside, byrow, bycolumn, checkerboard, custom, anaglyph, dynamic, clonedynamic (see also [stereo_mode](/index.php/Stereo_mode "Stereo mode")) 

[![](/images/b/b8/Stereo_on.png)](/index.php/File:Stereo_on.png)

[](/index.php/File:Stereo_on.png "Enlarge")

Example of 1ESR shown in cross-eyed stereo

## Example
    
    
    fetch 1ESR, async=0
    as cartoon
    set cartoon_smooth_loops
    spectrum
    bg white
    stereo crosseye
    

## See Also

  * [stereo_mode](/index.php/Stereo_mode "Stereo mode")
  * [stereo_angle](/index.php/Stereo_angle "Stereo angle")
  * [stereo_shift](/index.php?title=Stereo_shift&action=edit&redlink=1 "Stereo shift \(page does not exist\)")
  * [stereo_ray](/index.php/Stereo_ray "Stereo ray")



Retrieved from "[https://pymolwiki.org/index.php?title=Stereo&oldid=10665](https://pymolwiki.org/index.php?title=Stereo&oldid=10665)"


---

## Transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
    * 3.1 Whole Surface
    * 3.2 Selected Surface Elements



## Overview

**Transparency** is used to adjust the transparency of **[Surfaces](/index.php/Surface "Surface")** and **[Slices](/index.php/Slice "Slice")**. (For other transparencies in PyMOL, see [Cartoon Transparency](/index.php/Cartoon_transparency "Cartoon transparency"), [Sphere Transparency](/index.php/Sphere_transparency "Sphere transparency"), and [Stick Transparency](/index.php/Stick_transparency "Stick transparency"). 

## Usage
    
    
    set transparency, F, selection
    

where **F** is a floating point number in the range _[0.0 - 1.0]_ , where **selection** is the selected surface to apply the change to (for examples, see below). 

For the value of _F_ , 1.0 will be an invisible and 0.0 a completely solid surface. 

  


## Examples

### Whole Surface

Change the transparency of the whole surface to 50%. 
    
    
    # show all surfaces with 50% transparency.
    set transparency, 0.5
    

  


  * [![Image showing partial surface transparencies](/images/2/2e/Surf065.png)](/index.php/File:Surf065.png "Image showing partial surface transparencies")

Image showing partial surface transparencies 

  * [![Image showing 100% of the surface with 65% transparency.](/images/d/d6/Surfall_065.png)](/index.php/File:Surfall_065.png "Image showing 100% of the surface with 65% transparency.")

Image showing 100% of the surface with 65% transparency. 




### Selected Surface Elements

Simple example showing how to do partial surface transparency. This allows different selections to have different transparencies on the same object (or also on different objects). 
    
    
    # load a random protein
    fetch 1rty
    
    # set the partial transparency for the selected residues
    set transparency, 0.65, i. 1-100
    

  * [![Different transparency settings on different objects.](/images/f/f5/Diffy_transp.png)](/index.php/File:Diffy_transp.png "Different transparency settings on different objects.")

Different transparency settings on different objects. 

  * [![Different transparency settings on the same object.](/images/f/f3/Diffy_transp2.png)](/index.php/File:Diffy_transp2.png "Different transparency settings on the same object.")

Different transparency settings on the same object. 




Retrieved from "[https://pymolwiki.org/index.php?title=Transparency&oldid=10608](https://pymolwiki.org/index.php?title=Transparency&oldid=10608)"


---

