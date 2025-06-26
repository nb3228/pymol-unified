# Category: Coloring

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

## Bg Color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Bg_Color sets the background color 

### USAGE
    
    
    bg_color [color]
    

### PYMOL API
    
    
    cmd.bg_color(string color="black")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Bg_Color&oldid=7437](https://pymolwiki.org/index.php?title=Bg_Color&oldid=7437)"


---

## Bg rgb

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting is used to set the color of the background. Any color of the spectrum can be rendered. White backgrounds are often desirable for publication images. 

## Syntax

set bg_rgb,[_float1_ ,_float2_ ,_float3_] 

each float must be between 0.0 and 1.0 

float1: red component 

float2: green component 

float3: blue component 

  


example: set bg_rgb,[1,1,1] 

will color the background in white. 

## See Also

  * [Bg Color](/index.php/Bg_Color "Bg Color"), [bg_gradient](/index.php/Bg_gradient "Bg gradient")



Retrieved from "[https://pymolwiki.org/index.php?title=Bg_rgb&oldid=9542](https://pymolwiki.org/index.php?title=Bg_rgb&oldid=9542)"


---

## Cartoon color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Related settings
  * 4 See Also



## Overview

set **cartoon_color** allows one to color cartoons independently from the rest of the representations. 

## Syntax
    
    
    set cartoon_color, theColor
    

where _theColor_ can be : 

  * any usual colors (blue, yellow, grey50,...)
  * number-coded colors (1:black, 2:blue, 3:greenish, ...)
  * special code -1 to revert to original chameleon setting (_set cartoon_color,-1_)



## Related settings

[sphere_color](/index.php/Sphere_color "Sphere color")

## See Also

[Color](/index.php/Color "Color")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_color&oldid=5516](https://pymolwiki.org/index.php?title=Cartoon_color&oldid=5516)"


---

## Cartoon discrete colors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting toggles whether adjacent colors on a cartoon share a discrete boundary or whether they blend into each other. 

## Settings
    
    
    set cartoon_discrete_colors, on   # sets color blending off
    set cartoon_discrete_colors, off  # allows the blending of adjacent colors
    

## Examples

Open the images to actually see the details! 

  * [![cartoon_discrete_colors ON](/images/a/a9/Cartoon_discrete_colors_on.png)](/index.php/File:Cartoon_discrete_colors_on.png "cartoon_discrete_colors ON")

cartoon_discrete_colors ON 

  * [![cartoon_discrete_colors OFF](/images/4/47/Cartoon_discrete_colors_off.png)](/index.php/File:Cartoon_discrete_colors_off.png "cartoon_discrete_colors OFF")

cartoon_discrete_colors OFF 




Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_discrete_colors&oldid=5251](https://pymolwiki.org/index.php?title=Cartoon_discrete_colors&oldid=5251)"


---

## Cartoon highlight color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting allows one to specify a contrasting color for the interior face of helices and the side faces of strands. 

## Settings
    
    
    set cartoon_highlight_color, grey50   # sets cartoon highlight color to middle-grey
    set cartoon_highlight_color, -1      # turns this feature off (default)
    

## Examples

Open the images to actually see the details! 

  * [![cartoon_highlight_colors ON](/images/a/a9/Cartoon_discrete_colors_on.png)](/index.php/File:Cartoon_discrete_colors_on.png "cartoon_highlight_colors ON")

cartoon_highlight_colors ON 

  * [![cartoon_highlight_color OFF](/images/2/24/Cartoon_highlight_color_off.png)](/index.php/File:Cartoon_highlight_color_off.png "cartoon_highlight_color OFF")

cartoon_highlight_color OFF 




Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_highlight_color&oldid=9551](https://pymolwiki.org/index.php?title=Cartoon_highlight_color&oldid=9551)"


---

## Cartoon nucleic acid color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



## Overview

This setting determines the color of the sugar bb, when in nucleic acid mode. 

## Syntax
    
    
    # set the sugar bb to red
    set cartoon_nucleic_acid_color, red
    

## Examples

  * [![Color Red](/images/6/69/Cna_color2.png)](/index.php/File:Cna_color2.png "Color Red")

Color Red 

  * [![Color Blue](/images/e/e1/Cna_color1.png)](/index.php/File:Cna_color1.png "Color Blue")

Color Blue 




## See Also

[Cartoon#Nucleic_Acid_Representation](/index.php/Cartoon#Nucleic_Acid_Representation "Cartoon")

Retrieved from "[https://pymolwiki.org/index.php?title=Cartoon_nucleic_acid_color&oldid=5252](https://pymolwiki.org/index.php?title=Cartoon_nucleic_acid_color&oldid=5252)"


---

## CBC

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Util.cbc](/index.php/Util.cbc "Util.cbc") stands for (the utilities library's) **color by chain**. This simply colors molecules by their chains, as the name suggests. 

The affects of CBC are shown in the following images. 

  * [![Molecule loaded](/images/e/ed/Cbc0.png)](/index.php/File:Cbc0.png "Molecule loaded")

Molecule loaded 

  * [![util.cbc command issued. Each chain gets its own color.](/images/4/47/Cbc1.png)](/index.php/File:Cbc1.png "util.cbc command issued. Each chain gets its own color.")

util.cbc command issued. Each chain gets its own color. 




## Usage
    
    
    # simple command
    util.cbc selection, first_color, quiet
    
    # api usage
    util.cbc(selection='(all)',first_color=7,quiet=1,legacy=0,_self=cmd)
    

where 

  * selection defaults to 'all',
  * first_color defaults to 7,
  * quiet defaults to 1



## Example
    
    
    # color everything by chain
    util.cbc
    

## See Also

  * [Chainbow](/index.php/Advanced_Coloring#Coloring_with_.27chainbows.27_from_a_script "Advanced Coloring")
  * [Util.rainbow](/index.php?title=Util.rainbow&action=edit&redlink=1 "Util.rainbow \(page does not exist\)")
  * [util.cba](/index.php/Advanced_Coloring#Coloring_by_atom_type "Advanced Coloring")



Retrieved from "[https://pymolwiki.org/index.php?title=CBC&oldid=7440](https://pymolwiki.org/index.php?title=CBC&oldid=7440)"


---

## Cell color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Related settings
  * 4 See Also



## Overview

set **cell_color** allows one to color [cells](/index.php/Cell "Cell") independently from the rest of the representations. 

## Syntax
    
    
    set cell_color, theColor
    

where _theColor_ can be : 

  * any usual colors (blue, yellow, grey50,...)
  * number-coded colors (1:black, 2:blue, 3:greenish, ...)
  * special code -1 to revert to original chameleon setting (_set cartoon_color,-1_)



## Related settings

[sphere_color](/index.php/Sphere_color "Sphere color") [cartoon_color](/index.php/Cartoon_color "Cartoon color")

## See Also

[Color](/index.php/Color "Color")

Retrieved from "[https://pymolwiki.org/index.php?title=Cell_color&oldid=13573](https://pymolwiki.org/index.php?title=Cell_color&oldid=13573)"


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

## Color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**color** sets the color of an object or an atom selection to a predefined, named color, an RGB hex color, or a color [ramp](/index.php/Ramp_new "Ramp new"). For an overview of predefined colors, see [Color Values](/index.php/Color_Values "Color Values"). For a script that enumerates all the colors see, [List_Colors](/index.php/List_Colors "List Colors"). If you want to define your own colors, see [Set_Color](/index.php/Set_Color "Set Color"). 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
    * 3.1 Color all carbons yellow
    * 3.2 Color by element, except carbons
    * 3.3 Use a custom RGB color
    * 3.4 Color by Spectrum Example
    * 3.5 B-Factors
      * 3.5.1 Reassigning B-Factors and Coloring
      * 3.5.2 Reassigning B-Factors and Coloring - from file
    * 3.6 Expanding to Surface
    * 3.7 Getting Atom Colors
    * 3.8 What colors does PyMOL currently have?
    * 3.9 Color States Individually
  * 4 See Also



## Usage
    
    
    color color-name
    color color-name, object-name
    color color-name, (selection)
    

## Arguments

  * color-name = str or int: [named color](/index.php/Color_Values "Color Values"), index of a named color, `0xRRGGBB` string, `0x40RRGGBB` integer, name of a [ramp](/index.php/Ramp_new "Ramp new") object, or special value `atomic`
  * object-name = str: Name of an object or name pattern (with Asterisk (`*`)) which matches multiple objects. Changes the [object color](/index.php?title=Set_object_color&action=edit&redlink=1 "Set object color \(page does not exist\)"), and if it is a molecular object, also the color of all atoms. {default: all}
  * selection = str: [atom selection](/index.php/Selection_Algebra "Selection Algebra"), this is any expression which is not a valid object pattern. E.g. an expression which starts with a parentheses. This will only color atoms, not objects.



## Examples

### Color all carbons yellow
    
    
    color yellow, (name C*)
    

### Color by element, except carbons
    
    
    color atomic, (not elem C)
    

### Use a custom RGB color
    
    
    color 0x006600
    

### Color by Spectrum Example

Color by spectrum is in the GUI menu but did you realize that the spectrum is not limited to a simple rainbow? 
    
    
    spectrum count, palette, object_name
    

For available palettes and more details see: [spectrum](/index.php/Spectrum "Spectrum")

### B-Factors

The command to color a molecule by B-Factors (B Factors) is: 
    
    
    spectrum b, selection=SEL
    

where **SEL** is a valid selection, for example, "protA and n. CA", for protein A's alpha carbons. 

For more details see: [spectrum](/index.php/Spectrum "Spectrum")

#### Reassigning B-Factors and Coloring

It is commonplace to replace the B-Factor column of a protein with some other biochemical property at that residue, observed from some calculation or experiment. PyMOL can easily reassign the B-Factors and color them, too. The following example will load a protein, set ALL it's B Factors to "0", read in a list of properties for each alpha carbon in the proteins, assign those new values as the B-Factor values and color by the new values. This example is possible because commands PyMOL does not recognize are passed to the Python interpreter --- a very powerful tool. 
    
    
    # load the protein
    cmd.load("protA.pdb")
    
    # open the file of new values (just 1 column of numbers, one for each alpha carbon)
    inFile = open("newBFactors", 'r')
    
    # create the global, stored array
    stored.newB = []
    
    # read the new B factors from file
    for line in inFile.readlines(): stored.newB.append( float(line) )
    
    # close the input file
    inFile.close()
    
    # clear out the old B Factors
    alter protA, b=0.0
    
    # update the B Factors with new properties
    alter protA and n. CA, b=stored.newB.pop(0)
    
    # color the protein based on the new B Factors of the alpha carbons
    cmd.spectrum("b", "protA and n. CA")
    

If you want to save the file with the new B Factor values for each alpha carbon, 
    
    
    cmd.save("protA_newBFactors.pdb", "protA")
    

or similar is all you need. 

A script (data2bfactor.py) that loads data into the B-factor (b) or occupancy (q) columns from an external file can be found in Robert Campbell's PyMOL script repository (<http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/>) 

#### Reassigning B-Factors and Coloring - from file
    
    
    reinitialize
    prot="1XYZ"
    cmd.fetch(prot,async=0)
    
    # Set b value to zero
    cmd.alter(prot,b=0.0)
    cmd.show_as("cartoon",prot)
    
    python
    inFile = open("phi_values.txt", 'r')
    val_list = []
    for line in inFile.readlines()[1:]:
        split = line.split()
        resn = split[0][0] 
        resi = split[0][1:-1]
        phi_ppm2 = float(split[1])
        phi_ppm2_err = float(split[3])
        R2_0 = float(split[4])
        R2_0_err = float(split[6])
        print "Resn=%s Resi=%s, phi_ppm2=%2.2f, phi_ppm2_err=%2.2f, R2_0=%2.2f, R2_0_err=%2.2f"%(resn,resi,phi_ppm2,phi_ppm2_err,R2_0,R2_0_err)
    
        val_list.append(phi_ppm2)
        cmd.alter("%s and resi %s and n. CA"%(prot,resi), "b=%s"%phi_ppm2)
    
    python end
    minval = min(val_list)
    print minval
    maxval = max(val_list)
    print maxval
    cmd.spectrum("b", "blue_white_red", "%s and n. CA"%prot, minimum=0, maximum=maxval)
    cmd.ramp_new("ramp_obj", prot, range=[0, minval, maxval], color="[blue, white, red ]")
    cmd.save("%s_newBFactors.pdb"%prot, "%s"%prot)
    

### Expanding to Surface

See [Expand_To_Surface](/index.php/Expand_To_Surface "Expand To Surface"). 

If you have run the above code and would like the colors to be shown in the [Surface](/index.php/Surface "Surface") representation, too, then you need to do the following: 
    
    
    # Assumes alpha carbons colored from above.
    create ca_obj, your-object-name and name ca 
    ramp_new ramp_obj, ca_obj, [0, 10], [-1, -1, 0]
    set surface_color, ramp_obj, your-object-name
    

Thanks to Warren, for this one. 

### Getting Atom Colors
    
    
    stored.list = []
    iterate all, stored.list.append(color) # use cmd.get_color_tuple(color) to convert color index to RGB values
    print stored.list
    

Or, you can [label](/index.php/Label "Label") each atom by it's color code: 
    
    
    label all, color
    

### What colors does PyMOL currently have?

basic colors can be manually accessed and edited from the PyMOL menu under **Setting** \--> **Colors...**  
At the Pymol prompt, you can use the command [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices") to get a list of Pymols named colors.  
[Get_colors](/index.php/Get_colors "Get colors") is a script that allows accessing colors as well.  


### Color States Individually
    
    
    fetch 1nmr
    set all_states
    
    # the object has 20 states, so we can set separate line colors
    # for each state as follows:
    for a in range(1,21): cmd.set("line_color","auto","1nmr",a)
    

Or, we can do it differently, 
    
    
    # start over,
    fetch 1nmr
    
    # break apart the object by state
    split_states 1nmr
    
    # delete the original
    dele 1nmr
    
    # and color by object (carbons only)
    util.color_objs("elem c")
    
    # (all atoms)
    util.color_objs("all")
    

## See Also

  * [Advanced Coloring](/index.php/Advanced_Coloring "Advanced Coloring")
  * [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices")
  * [spectrum](/index.php/Spectrum "Spectrum")
  * [Ramp_New](/index.php/Ramp_New "Ramp New")



Retrieved from "[https://pymolwiki.org/index.php?title=Color&oldid=13190](https://pymolwiki.org/index.php?title=Color&oldid=13190)"


---

## Color Objects

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | Gareth Stockwell   
License  |   
      
    
    #####################################################################
    #
    # Colour by object
    #
    #####################################################################
     
    from pymol import cmd
    
    def color_obj(rainbow=0):
     
            """
     
    AUTHOR 
     
            Gareth Stockwell
     
    USAGE
     
            color_obj(rainbow=0)
     
            This function colours each object currently in the PyMOL heirarchy
            with a different colour.  Colours used are either the 22 named
            colours used by PyMOL (in which case the 23rd object, if it exists,
            gets the same colour as the first), or are the colours of the rainbow
     
    SEE ALSO
    
            util.color_objs()
            """
     
            # Process arguments
            rainbow = int(rainbow)
     
            # Get names of all PyMOL objects
            obj_list = cmd.get_names('objects')
     
            if rainbow:
     
               print "\nColouring objects as rainbow\n"
     
               nobj = len(obj_list)
     
               # Create colours starting at blue(240) to red(0), using intervals
               # of 240/(nobj-1)
               for j in range(nobj):
                  hsv = (240-j*240/(nobj-1), 1, 1)
                  # Convert to RGB
                  rgb = hsv_to_rgb(hsv)
                  # Define the new colour
                  cmd.set_color("col" + str(j), rgb)
                  print obj_list[j], rgb
                  # Colour the object
                  cmd.color("col" + str(j), obj_list[j])
     
            else:
     
               print "\nColouring objects using PyMOL defined colours\n"
     
               # List of available colours
               colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
               'wheat', 'white', 'grey' ]
               ncolours = len(colours)
     
               # Loop over objects
               i = 0
               for obj in obj_list:
                  print "  ", obj, colours[i]
                  cmd.color(colours[i], obj)
                  i = i+1
                  if(i == ncolours):
                     i = 0
     
     
    # HSV to RGB routine taken from Robert L. Campbell's color_b.py script
    #   See http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
    # Original algorithm from: http://www.cs.rit.edu/~ncs/color/t_convert.html
    def hsv_to_rgb(hsv):
     
            h = float(hsv[0])
            s = float(hsv[1])
            v = float(hsv[2])
     
            if( s == 0 ) :
                    #achromatic (grey)
                    r = g = b = v
     
            else:
                    # sector 0 to 5
                    h = h/60.            
                    i = int(h)
                    f = h - i                       # factorial part of h
                    #print h,i,f
                    p = v * ( 1 - s )
                    q = v * ( 1 - s * f )
                    t = v * ( 1 - s * ( 1 - f ) )
     
                    if i == 0:
                            (r,g,b) = (v,t,p)
                    elif i == 1:
                            (r,g,b) = (q,v,p)
                    elif i == 2:
                            (r,g,b) = (p,v,t)
                    elif i == 3:
                            (r,g,b) = (p,q,v)
                    elif i == 4:
                            (r,g,b) = (t,p,v)
                    elif i == 5:
                            (r,g,b) = (v,p,q)
                    else:
                            (r,g,b) = (v,v,v)
                            print "error, i not equal 1-5"
     
            return [r,g,b]
     
     
     
    # Add color_obj to the PyMOL command list 
    cmd.extend("color_obj",color_obj)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_Objects&oldid=10269](https://pymolwiki.org/index.php?title=Color_Objects&oldid=10269)"


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

## Colorama

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/colorama.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/colorama.py)  
Author(s)  | [Gregor Hagelueken](/index.php?title=User:Gha&action=edit&redlink=1 "User:Gha \(page does not exist\)")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 COLORAMA
    * 1.1 Program features
    * 1.2 Screenshot
    * 1.3 Usage
      * 1.3.1 Single color
      * 1.3.2 Color gradients
    * 1.4 Contact



# COLORAMA

COLORAMA is a PyMOL plugin which allows to color objects using adjustable scale bars. 

## Program features

In RGB mode, each color R, G, B is represented by one scale bar which can be manually adjusted while the selected object is colored in real time. For convenience, it is as well possible to switch to the HSV color system. 

Additionally, a color gradient with user-defined start- and end-colors can be created for the selected molecule. 

## Screenshot

[![COLORAMA-screenshot.jpg](/images/9/95/COLORAMA-screenshot.jpg)](/index.php/File:COLORAMA-screenshot.jpg)

## Usage

This plugin is included in the project [ Pymol-script-repo](/index.php/Git_intro "Git intro"). 

Manual install the plugin by copying the code below into an empty text file (e.g. "colorama.py") located in the \Pymol\modules\pmg_tk\startup directory. After PyMOL has been started, the program can be launched from the PLUGINS menu. COLORAMA has not been tested with PyMOL versions older than 1.0. 

### Single color

  1. Type the name of a PyMOL object to be colored into the COLORAMA entry field.
  2. Push the "Set" button.
  3. The scales are adjusted to the current color which is additionally visualized in a field left to the scales.
  4. If one of the scales is moved, the color of the selected object will change in real-time.
  5. Pushing the RGB or HSV buttons on the left allows to switch between both color systems.



### Color gradients

  1. After an object has been selected, push the "G" button (gradient).
  2. Select the start color by pushing "C1" and adjusting it using the scales.
  3. Select the end color "C2" in the same way.
  4. To create the gradient, push "Set gradient".



A new object will be created which is called "dummy-OLD_OBJECT". The B-factor column of this object is overwritten and now contains the number of each residue. The original object is left unchanged. The gradient mode can be left by pushing "M" (monochrome). This part of the program uses a modified version of the [color_b script](https://github.com/zigeuner/robert_campbell_pymol_scripts/tree/master/work_pymol) by Robert L. Campbell & James Stroud. 

## Contact

Gregor Hagelueken, hagelueken'at'pc.uni-bonn.de 

Retrieved from "[https://pymolwiki.org/index.php?title=Colorama&oldid=13845](https://pymolwiki.org/index.php?title=Colorama&oldid=13845)"


---

## Colorblindfriendly

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/colorblindfriendly.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorblindfriendly.py)  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  


## Contents

  * 1 Introduction
  * 2 Colors
  * 3 Usage
    * 3.1 Import as a module
    * 3.2 Run the latest version of the script from Github
    * 3.3 Download the script and run locally
  * 4 Requirements



## Introduction

Certain colors are indistinguishable to people with the various forms of color blindness, and therefore are better not used in figures intended for public viewing. 

This script generates a palette of named colors for PyMOL that are unambiguous both to colorblind and non-colorblind individuals. 

The colors listed here are defined according to recommendations by Okabe and Ito previously available at [J*FLY (archived)](https://web.archive.org/web/20170608112249/http://jfly.iam.u-tokyo.ac.jp/color/#pallet). This website is a good reference to consult when making all kinds of figures, not just those made using PyMOL. 

## Colors

These are the 0-255 RGB values from the J*FLY page that are used in the script, with the defined color names and alternate names. 

| name | R | G | B | alternate names   
---|---|---|---|---|---  
| cb_black | 0 | 0 | 0 |   
| cb_orange | 230 | 159 | 0 |   
| cb_sky_blue | 86 | 180 | 233 | cb_skyblue, cb_light_blue, cb_lightblue   
| cb_bluish_green | 0 | 158 | 115 | cb_bluishgreen, cb_green   
| cb_yellow | 240 | 228 | 66 |   
| cb_blue | 0 | 114 | 178 |   
| cb_vermillion | 213 | 94 | 0 | cb_red, cb_red_orange, cb_redorange   
| cb_reddish_purple | 204 | 121 | 167 | cb_rose, cb_violet, cb_magenta   
  
## Usage

### Import as a module

[![](/images/c/ce/Colorblindfriendly_menu.png)](/index.php/File:Colorblindfriendly_menu.png)

[](/index.php/File:Colorblindfriendly_menu.png "Enlarge")

Screenshot of the `cb_colors` color menu in the OpenGL GUI in PyMOL 2.0.

After importing the module, call the `set_colors()` function to add the colors to PyMOL's color palette. Then, use these color names just like any other named color, using the `[color](/index.php/Color "Color")` command. 
    
    
    import colorblindfriendly as cbf
    
    # Add the new colors
    cbf.set_colors()
    color cb_red, myObject
    

The colors can also be made to replace the built-in colors (i.e. they are created both with and without the "`cb_`" prefix.). Do this by passing the `replace` keyword argument. 
    
    
    # Replace built-in colors with cbf ones
    cbf.set_colors(replace=True)
    color yellow, myOtherObject   # actually cb_yellow
    

One can also add an entry to the color menu in the right-side OpenGL GUI. So clicking on [C], there will now be a `cb_colors` menu item, which expands to give all the color blind-friendly colors, except black, which is available in the `grays` menu. 
    
    
    # Add a `cb_colors` menu to the OpenGL GUI (see screenshot)
    # This will also add the colors if `set_colors()` hasn't yet been run.
    cbf.add_menu()
    

### Run the latest version of the script from Github

In a PyMOL session, run at the command line: 
    
    
    run <https://github.com/Pymol-Scripts/Pymol-script-repo/raw/master/colorblindfriendly.py>
    

This will add all the colors as well as the OpenGL menu. 

  


### Download the script and run locally

Save the script from the link in the box at the top right to your computer. 

In a PyMOL session, run at the command line: 
    
    
    run /path/to/colorblindfriendly.py
    

This will add all the colors as well as the OpenGL menu. 

## Requirements

The `cb_colors` GUI menu (generated by the `add_menu()` function) requires PyMOL 1.6.0 and later. 

Retrieved from "[https://pymolwiki.org/index.php?title=Colorblindfriendly&oldid=13896](https://pymolwiki.org/index.php?title=Colorblindfriendly&oldid=13896)"


---

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

## Get color index

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_color_index** in combination with **get_color_tuple** will retrieve the RGB values for a single color. 
    
    
    PyMOL>print cmd.get_color_index('black')
    1
    PyMOL>print cmd.get_color_tuple(1)
    (0.0, 0.0, 0.0)

  


## See Also

  * [Get_Color_Indices](/index.php/Get_Color_Indices "Get Color Indices")
  * [Get_Color_Tuples](/index.php/Get_Color_Tuples "Get Color Tuples")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_color_index&oldid=12665](https://pymolwiki.org/index.php?title=Get_color_index&oldid=12665)"


---

## Get Color Indices

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_color_indices** in combination with **get_color_tuple** will retrieve the RGB values for colors. 
    
    
    print cmd.get_color_indices()
    

will retrieve the Pymol color names and corresponding internal color indices. The Pymol names can be used to designate color for objects, see [Color](/index.php/Color "Color"). To retrieve a single index for a specific color name, use **[get_color_index](/index.php/Get_color_index "Get color index")** instead. 
    
    
    print cmd.get_color_tuple(index-number)
    

will retrieve individual RGB components when _index-number_ is replaced with one of the color indices from above. 

The color index, an integer, gets returned when color is returned while employing [Iterate](/index.php/Iterate "Iterate"). You can thus use the **get_color_tuple** command above to convert that to RGB color values if you need to use the colors outside Pymol. 

Tangentially related is the fact you can name additional colors, 
    
    
    set_color color-name, [r,b,g]
    

will create a new color that will appear in the GUI list. From the open-source GUI you can use the "add" button in the color list viewer. In MacPyMOL, enter the new name into the MacPyMOL color editor window, set the RGBs, and then click Apply. See [Set Color](/index.php/Set_Color "Set Color") for more details and examples. The colors created will be added to the end of the list of Pymol's color indices that you can view the **get_color_indices()** command. 

### EXAMPLES
    
    
    # rename a color 
    cmd.set_color('myfavoritecolor', cmd.get_color_tuple(cmd.get_color_index('red')))
    

  


## See Also

  * [Get_color_index](/index.php/Get_color_index "Get color index")
  * [Get_Color_Tuples](/index.php/Get_Color_Tuples "Get Color Tuples")
  * [Iterate](/index.php/Iterate "Iterate")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_Color_Indices&oldid=12766](https://pymolwiki.org/index.php?title=Get_Color_Indices&oldid=12766)"


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

## Label outline color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Sets which color PyMol uses to outlines the labels. This helps offset labels from their backgrounds. 

## Syntax
    
    
    set label_outline_color, red  # will draw red outline around the labels
    

Retrieved from "[https://pymolwiki.org/index.php?title=Label_outline_color&oldid=6222](https://pymolwiki.org/index.php?title=Label_outline_color&oldid=6222)"


---

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

## List Colors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    # This is how to do it from the PyMOL command line or .pml script:
    #
    iterate all, print color
    
    #! /usr/bin/python
    #
    # and this in a Python script
    #
    import pymol
    pymol.color_list = []
    cmd.iterate('all', 'pymol.color_list.append(color)')
    print pymol.color_list
    

Retrieved from "[https://pymolwiki.org/index.php?title=List_Colors&oldid=6345](https://pymolwiki.org/index.php?title=List_Colors&oldid=6345)"


---

## Load new B-factors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [figshare](http://dx.doi.org/10.6084/m9.figshare.1176991)  
Author(s)  | [Pietro Gatti-Lafranconi](/index.php/User:PietroGattiLafranconi "User:PietroGattiLafranconi")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Limitations
  * 6 Examples
  * 7 The Code



## Introduction

Quick and simple script to replace B-factor values in a PDB structure. 

New B-factors are provided in an external .txt file, one per line. 

By default, the script will also redraw the selected molecule as cartoon putty and colour by B-factor 

  


## Usage
    
    
    loadBfacts mol, [startaa, [source, [visual Y/N]]]
    

  


## Required Arguments

  * **mol** = any object selection (within one single object though)



  


## Optional Arguments

  * **startaa** = number of amino acid the first value in 'source' refers to (default=1)
  * **source** = name of the file containing new B-factor values (default=newBfactors.txt)
  * **visual** = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)



## Limitations

For its very nature, this script is not suitable for complex cases in which only some B-factors are to be replaced, or on complex selections. In such cases, [data2bfactor.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/) is a better choice. 

B-factors have to be entered one per line, with no labels or amino acid ID. They will replace values of a continuous stretch of amino acids of equal length (or shorter) starting from the position provided with 'startaa'. 

'newBfactors.txt' looks like 
    
    
    2
    0
    0
    3
    7
    ...
    

  


## Examples
    
    
    PyMOL>loadbfacts 1LVM
    

[![example #1](/images/8/86/LoadBfacts1.png)](/index.php/File:LoadBfacts1.png "example #1")

  

    
    
    PyMOL>loadbfacts 1LVM and chain a, startaa=135
    

[![example #2](/images/2/2b/LoadBfacts02.png)](/index.php/File:LoadBfacts02.png "example #2")

  


## The Code
    
    
    from pymol import cmd, stored, math
    	
    def loadBfacts (mol,startaa=1,source="newBfactors.txt", visual="Y"):
    	"""
    	Replaces B-factors with a list of values contained in a plain txt file
    	
    	usage: loadBfacts mol, [startaa, [source, [visual]]]
     
    	mol = any object selection (within one single object though)
    	startaa = number of first amino acid in 'new B-factors' file (default=1)
    	source = name of the file containing new B-factor values (default=newBfactors.txt)
    	visual = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)
     
    	example: loadBfacts 1LVM and chain A
    	"""
    	obj=cmd.get_object_list(mol)[0]
    	cmd.alter(mol,"b=-1.0")
    	inFile = open(source, 'r')
    	counter=int(startaa)
    	bfacts=[]
    	for line in inFile.readlines():	
    		bfact=float(line)
    		bfacts.append(bfact)
    		cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%bfact)
    		counter=counter+1
    	if visual=="Y":
    		cmd.show_as("cartoon",mol)
    		cmd.cartoon("putty", mol)
    		cmd.set("cartoon_putty_scale_min", min(bfacts),obj)
    		cmd.set("cartoon_putty_scale_max", max(bfacts),obj)
    		cmd.set("cartoon_putty_transform", 0,obj)
    		cmd.set("cartoon_putty_radius", 0.2,obj)
    		cmd.spectrum("b","rainbow", "%s and n. CA " %mol)
    		cmd.ramp_new("count", obj, [min(bfacts), max(bfacts)], "rainbow")
    		cmd.recolor()
    
    cmd.extend("loadBfacts", loadBfacts);
    

Retrieved from "[https://pymolwiki.org/index.php?title=Load_new_B-factors&oldid=11711](https://pymolwiki.org/index.php?title=Load_new_B-factors&oldid=11711)"


---

## Mesh color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

set mesh_color simply sets the colour of the mesh 

## Syntax
    
    
    set mesh_color, <colour>            #default is what ever the protein already is
    

## Example
    
    
    set mesh_color, blue
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mesh_color&oldid=5244](https://pymolwiki.org/index.php?title=Mesh_color&oldid=5244)"


---

## Palette Colorbars

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The palette colorbars herein correspond to the palettes enumerated on the [Spectrum](/index.php/Spectrum#Notes "Spectrum") page. 

  


## Contributed Public Domain Colorbars

These colorbars are placed in the public domain and may be used in any publication freely. 

_blue_green_

    

    [![Colorbar-blue green.png](/images/5/5d/Colorbar-blue_green.png)](/index.php/File:Colorbar-blue_green.png)

_blue_magenta_

    

    [![Pymol-blue magenta.png](/images/9/96/Pymol-blue_magenta.png)](/index.php/File:Pymol-blue_magenta.png)

_blue_red_

    

    [![Pymol-blue red.png](/images/d/dc/Pymol-blue_red.png)](/index.php/File:Pymol-blue_red.png)

_blue_white_green_

    

    [![Pymol-blue white green.png](/images/8/88/Pymol-blue_white_green.png)](/index.php/File:Pymol-blue_white_green.png)

_blue_white_magenta_

    

    [![Pymol-blue white magenta.png](/images/6/67/Pymol-blue_white_magenta.png)](/index.php/File:Pymol-blue_white_magenta.png)

_rainbow_

    

    [![Pymol-rainbow.png](/images/8/88/Pymol-rainbow.png)](/index.php/File:Pymol-rainbow.png)

_red_blue_

    

    [![Pymol-red blue.png](/images/b/b1/Pymol-red_blue.png)](/index.php/File:Pymol-red_blue.png)

# See Also

  * [Color](/index.php/Color "Color")
  * [Spectrum](/index.php/Spectrum "Spectrum")



Retrieved from "[https://pymolwiki.org/index.php?title=Palette_Colorbars&oldid=9609](https://pymolwiki.org/index.php?title=Palette_Colorbars&oldid=9609)"


---

## Ramp New

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[ramp_new](/index.php/Ramp_new "Ramp new") creates a color ramp based on a map potential value or based on proximity to a molecular object. 

## Contents

  * 1 Usage
  * 2 Color and range interpretation
    * 2.1 Number of colors and values
    * 2.2 Supported special colors
  * 3 Examples
    * 3.1 Color surface by an APBS calculated electrostatic map
    * 3.2 Color binding pocket by proximity to ligand
    * 3.3 Color isosurface by closest atom
    * 3.4 More user provided examples
      * 3.4.1 Ramp + Distance Measure
      * 3.4.2 Surface Colored by Distance from a Point
      * 3.4.3 Coloring a Viral Capsid by Distance from Core
  * 4 See Also



## Usage
    
    
    ramp_new name, map_name [, range [, color [, state [, selection [,
           beyond [, within [, sigma [, zero ]]]]]]]]
    

  * **name** = string: name of the ramp object to create
  * **map_name** = string: name of the map (for potential) or molecular object (for proximity)
  * **range** = list: values corresponding to slots in the ramp {default: [-1, 0, 1]}
  * **color** = list or str: colors corresponding to slots in the ramp or a given palette name: afmhot, grayscale, object, rainbow, traditional, grayable, hot, ocean, and sludge {default: [red, white, blue]}
  * **state** = integer: state identifier {default: 1}
  * **selection** = selection: for automatic ranging (maps only) {default: }
  * with automatic ranging only: 
    * **beyond** = number: are we excluding values beyond a certain distance from the selection? {default: 2.0}
    * **within** = number: are we only including valuess within a certain distance from the selection? {default: 6.0}
    * **sigma** = number: how many standard deviations from the mean do we go? {default: 2.0}
    * **zero** = integer: do we force the central value to be zero? {default: 1}



## Color and range interpretation

### Number of colors and values

  * Obvious case: Same number of values and colors (_len(range) == len(color)_)
  * Color palettes: Only first and last value are used as the palette boundaries, with linear interpolation in between
  * Less colors than values: The last color will be repeated to fill up the color array
  * More colors than values: **This behavior is deprecated and will change in PyMOL 1.8.**. If **N** is the number of values, then: 
    * One extra color: Last color (_color[N]_) will be used for _values > range[N-1]_
    * Two extra colors: _color[N]_ will be used for _values < range[0]_ and _color[N+1]_ will be used for _values > range[N-1]_
    * More than two extra colors: Like with two extra colors, but _color[N+2:]_ will be ignored
    * **Recommended practice: Instead of providing out-of-bounds colors with the last two colors, put them in the beginning and end of the color list and repeat the first and last values in the range list.** Example with white out-of-bounds coloring: _range=[0, 0, 10, 10], color=[white, red, blue, white]_
    * **Planned change in PyMOL 1.8:** With two values and more than two colors, colors will be evenly spaced between the two values (like color palettes)



### Supported special colors

With proximity ramps, in addition to regular named colors (like "blue"), the following special colors are supported: 

  * atomic: color of closest atom in proximity object
  * default: alias for atomic
  * object: color of proximity object
  * <name of another ramp>: recursive ramp coloring



## Examples

### Color surface by an APBS calculated electrostatic map

Example map: <http://pymol.org/tmp/1ubq_apbs.dx>
    
    
    load 1ubq_apbs.dx, e_pot_map
    fetch 1ubq, async=0
    as surface
    
    ramp_new e_pot_color, e_pot_map, [-5, 0, 5], [red, white, blue]
    set surface_color, e_pot_color
    set surface_ramp_above_mode
    

### Color binding pocket by proximity to ligand

  * Ligand: blue
  * Protein within 4 Angstrom of ligand: red
  * Protein beyond 8 Angstrom of ligand: yellow


    
    
    fetch 1rx1, async=0
    extract ligand, organic
    color blue, ligand
    show surface
    
    ramp_new prox, ligand, [4, 8], [red, yellow]
    color prox, 1rx1
    

### Color isosurface by closest atom

_See also the[huge surfaces](/index.php/Huge_surfaces "Huge surfaces") example._

Color the isosurface within 2 Angstrom of the protein (without solvent) by atom colors. Everything beyond 2 Angstrom will be gray. 
    
    
    fetch 1rx1, map, async=0, type=2fofc
    isosurface surf, map
    
    fetch 1rx1, mol, async=0
    remove solvent
    
    ramp_new prox, mol, [0, 2, 2], [atomic, atomic, gray]
    color prox, surf
    

Updating the atom colors is possible, for example: 
    
    
    spectrum
    recolor
    

### More user provided examples

#### Ramp + Distance Measure

Using a ramp and a distance measure, we can color the surface by some property--here, I'll chose distance from some important atom in the receptor to the ligand atom. 

  * [![Ligand surface colored by distance from some given atom. The remainder of the protein is hidden to more clearly visualize the calculated distances and surface color](/images/9/9a/Surface_by_prop3.png)](/index.php/File:Surface_by_prop3.png "Ligand surface colored by distance from some given atom. The remainder of the protein is hidden to more clearly visualize the calculated distances and surface color")

Ligand surface colored by distance from some given atom. The remainder of the protein is hidden to more clearly visualize the calculated distances and surface color 

  * [![Surface by prop.png](/images/5/5e/Surface_by_prop.png)](/index.php/File:Surface_by_prop.png)




To reproduce the results shown here you must do the following: 

  * obtain a protein
  * calculate some property for some set of atoms (like distance from some central location) and stuff the values into the b-factor
  * create a new object from the atoms for which you just measured a property
  * create a new ramp from the object with ramp_new
  * set the surface color of the new object



  
Another possible application of the ramp_new command can be the representation of the ELF function [[1]](http://en.wikipedia.org/wiki/Electron_localization_function). This function can be calculated with the TopMod software [[2]](http://www.lct.jussieu.fr/suite64.html). 

  * Load the cube file containing the ELF function, e.g. H2O_elf.cube.
  * Create an isosurface with a contour level of 0.8.


    
    
    isosurface elf, H2O_elf, 0.8
    

  * Load the cube containing the basin information, e.g. H20_esyn.cube. Basically in this cube for each point in the first cube you have either one of the numbers from 1 to 5. More details on what these numbers mean can be found in the TopMod manual.
  * Create a new ramp.


    
    
    ramp_new ramp, H2O_esyn, [1, 2, 3, 5], [tv_orange, lightblue, palegreen, deeppurple]
    

  * Assign the color ramp to the ELF isosurface.


    
    
    set surface_color, ramp, elf
    

  * Rebuild if necessary.


    
    
    rebuild
    

  * [![Localization domains of H2O. The bounding isosurface is ELF=0.8](/images/4/44/H2O.png)](/index.php/File:H2O.png "Localization domains of H2O. The bounding isosurface is ELF=0.8")

Localization domains of H2O. The bounding isosurface is ELF=0.8 

  * [![Localization domains of BeCl2. The bounding isosurface is ELF=0.8](/images/c/cd/BeCl2.png)](/index.php/File:BeCl2.png "Localization domains of BeCl2. The bounding isosurface is ELF=0.8")

Localization domains of BeCl2. The bounding isosurface is ELF=0.8 




#### Surface Colored by Distance from a Point

See [Spectrum](/index.php/Spectrum "Spectrum") for another method that allows for more flexible coloring schemes, but needs more work to get there. 

  * [![Measure.png](/images/f/ff/Measure.png)](/index.php/File:Measure.png)




This example shows how to color a protein surface by its distance from a given point: 
    
    
    # fetch a friendly protein
    
    fetch 1hug, async=0
    
    # show it as a surface
    
    as surface
    
    # create a pseudoatom at the origin; we will
    # measure the distance from this point
    
    pseudoatom pOrig, pos=(0,0,0), label=origin
    
    # create a new color ramp, measuring the distance
    # from pOrig to 1hug, colored as rainbow
    
    ramp_new proximityRamp, pOrig, selection=1hug, range=[5,65], color=rainbow
    
    # set the surface color to the ramp coloring
    
    set surface_color, proximityRamp, 1hug
    
    # some older PyMOLs need this recoloring/rebuilding
    
    recolor; rebuild
    

#### Coloring a Viral Capsid by Distance from Core

  * [![Capsid by dist.png](/images/3/33/Capsid_by_dist.png)](/index.php/File:Capsid_by_dist.png)



    
    
    # create a pseudoatom at the origin-- we will
    # measure the distance from this point
    
    pseudoatom pOrig, pos=(0,0,0), label=origin
    
    # fetch and build the capsid
    
    fetch 2xpj, async=0, type=pdb1
    split_states 2xpj
    delete 2xpj
    
    # show all 60 subunits it as a surface
    # this will take a few minutes to calculate
    
    as surface
    
    # create a new color ramp, measuring the distance
    # from pOrig to 1hug, colored as rainbow
    
    ramp_new proximityRamp, pOrig, selection=(2xpj*), range=[110,160], color=rainbow
    
    # set the surface color to the ramp coloring
    
    set surface_color, proximityRamp, (2xpj*)
    
    # some older PyMOLs need this recoloring/rebuilding
    
    recolor
    

## See Also

[load](/index.php/Load "Load"), [color](/index.php/Color "Color"), [create](/index.php/Create "Create"), [slice](/index.php/Slice "Slice"), [gradient](/index.php/Gradient "Gradient"), [Expand_To_Surface](/index.php/Expand_To_Surface "Expand To Surface")

Retrieved from "[https://pymolwiki.org/index.php?title=Ramp_New&oldid=12061](https://pymolwiki.org/index.php?title=Ramp_New&oldid=12061)"


---

## Rendering plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/rendering_plugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/rendering_plugin.py)  
Author(s)  | [Michael G. Lerner](/index.php/User:Mglerner "User:Mglerner")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Description

Here is a small plugin to render images with a given DPI. 

The "Draw" button just draws the image without raytracing (a fast way to see that the height/width look good).  
The "Ray" button raytraces first, before saving. 

The functionality is also available as a script (see my .pymolrc [here](/index.php/User:Mglerner "User:Mglerner")). 

First the imperial version. The metric version follows. 

To install, save the script as e.g. rendering_plugin.py or rendering_plugin_metric.py and install via PyMOL's Plugin --> Manage Plugins --> Install menu. 

The plugins are available through the project, [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro). 

Retrieved from "[https://pymolwiki.org/index.php?title=Rendering_plugin&oldid=10439](https://pymolwiki.org/index.php?title=Rendering_plugin&oldid=10439)"


---

## Resicolor plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/resicolor_plugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/resicolor_plugin.py)  
Author(s)  | [Philaltist](/index.php/User:Philaltist "User:Philaltist")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Description

A small plugin to color proteins according to residue type. This functionality is also available as a script ([resicolor](/index.php/Resicolor "Resicolor")). 

## Example of use

If you follow the instructions for the [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro), it is available from the Plugin menu. 

Retrieved from "[https://pymolwiki.org/index.php?title=Resicolor_plugin&oldid=10440](https://pymolwiki.org/index.php?title=Resicolor_plugin&oldid=10440)"


---

## Resicolor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Description

Here is a small script to color proteins according to residue type. Call it from your .pymolrc or from within pymol: "run resicolor.py". Then, issuing "resicolor" will apply the scheme. Selections are supported: "resicolor chain A" will apply the scheme only to chain A. This functionality is also available as a plugin ([Resicolor_plugin](/index.php/Resicolor_plugin "Resicolor plugin")). 
    
    
    #script contributed by Philippe Garteiser; garteiserp@omrf.ouhsc.edu
    from pymol import cmd
    
    def resicolor(selection='all'):
        
        '''USAGE: resicolor <selection>
        colors all or the given selection with arbitrary
        coloring scheme.
        '''
        cmd.select ('calcium','resn ca or resn cal')
        cmd.select ('acid','resn asp or resn glu or resn cgu')
        cmd.select ('basic','resn arg or resn lys or resn his')
        cmd.select ('nonpolar','resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
        cmd.select ('polar','resn ser or resn thr or resn asn or resn gln or resn tyr')
        cmd.select ('cys','resn cys or resn cyx')
        cmd.select ('backbone','name ca or name n or name c or name o')
        cmd.select ('none')
    
        print selection
        code={'acid'    :  'red'    ,
              'basic'   :  'blue'   ,
              'nonpolar':  'orange' ,
              'polar'   :  'green'  ,
              'cys'     :  'yellow'}
        cmd.select ('none')
        for elem in code:
            line='color '+code[elem]+','+elem+'&'+selection
            print line
            cmd.do (line)
        word='color white,backbone &'+selection
        print word
        cmd.do (word)                  #Used to be in code, but looks like
                                       #dictionnaries are accessed at random
        cmd.hide ('everything','resn HOH')
    
    cmd.extend ('resicolor',resicolor)
    

## Download

[Resicolor-0.1.tar.zip‎](/images/f/ff/Resicolor-0.1.tar.zip "Resicolor-0.1.tar.zip")

Retrieved from "[https://pymolwiki.org/index.php?title=Resicolor&oldid=12198](https://pymolwiki.org/index.php?title=Resicolor&oldid=12198)"


---

## Ribbon color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 Advanced details



## Overview

Ribbon_color allows one to explicitly state the color to be applied to a ribbon object. 

## Syntax
    
    
    # set it to a color
    set ribbon_color, color
    
    # examples:
    # default auto-cycles the color for each new object
    set ribbon_color, marine
    
    # apply green to the ribbon only for object obj01
    set ribbon_color, green /obj01
    
    # apply orange to the ribbon only for chain B of object obj01; NB: not fully enabled until version 1.0
    set ribbon_color, orange /obj01//B
    

## Examples

  * [![ribbon_color green](/images/1/17/Ribbon_color_green.png)](/index.php/File:Ribbon_color_green.png "ribbon_color green")

ribbon_color green 

  * [![ribbon_color marine](/images/0/0d/Ribbon_color_marine.png)](/index.php/File:Ribbon_color_marine.png "ribbon_color marine")

ribbon_color marine 

  * [![two chains in same object with different colors](/images/f/f1/Ribbon_color_chains.png)](/index.php/File:Ribbon_color_chains.png "two chains in same object with different colors")

two chains in same object with different colors 




  


## Advanced details

Some notes on general [Set](/index.php/Set "Set") syntax (From PyMOL "help set"): 

DESCRIPTION 
    
    
      "set" changes one of the PyMOL state variables,
    
    

USAGE 
    
    
      set name, [,value [,object-or-selection [,state ]]]
    
      set name = value  # (DEPRECATED)
    
    

PYMOL API 
    
    
      cmd.set ( string name, string value=1,
                string selection=_, int state=0,_
                 int updates=1, quiet=1)
    
    

NOTES 
    
    
      The default behavior (with a blank selection) changes the global
      settings database.  If the selection is 'all', then the settings
      database in all individual objects will be changed.  Likewise, for
      a given object, if state is zero, then the object database will be
      modified.  Otherwise, the settings database for the indicated state
      within the object will be modified.
    
      If a selection is provided, then all objects in the selection will
      be affected.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ribbon_color&oldid=6230](https://pymolwiki.org/index.php?title=Ribbon_color&oldid=6230)"


---

## Selection Macros

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/a/ac/Right-click-macro.png)](/index.php/File:Right-click-macro.png)

The atom-right-click context menu title shows a selection macro

Selection Macros allow to represent a long [atom selection](/index.php/Selection_Algebra "Selection Algebra") phrase such as 
    
    
    PyMOL> select pept and segi lig and chain B and resi 142 and name CA
    

in a more compact form: 
    
    
    PyMOL> select /pept/lig/B/142/CA
    

## Contents

  * 1 Syntax and Semantics
  * 2 Details
    * 2.1 Macros Beginning With a Slash
    * 2.2 Macros Not Beginning With a Slash
    * 2.3 Omitting Fields in a Macro
  * 3 See Also



## Syntax and Semantics

An atom selection macro uses slashes to define fields corresponding to identifiers. The macro is used to select atoms using the boolean "and," that is, the selected atoms must have all the matching identifiers: 
    
    
    /object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
    

  * must contain **at least one slash**
  * no spaces allowed
  * **empty** fields are interpreted as **wildcards**
  * starting with a slash: 
    * **yes** : fields from the **right** can be omitted
    * **no** : fields from the **left** can be omitted



## Details

Macros come in two flavors: those that begin with a slash and those that don't. The presence or absence of a slash at the beginning of the macro determines how it is interpreted. If the macro begins with a slash, PyMOL expects to find the fields starting from the top of the hierarchy: the first field to the right of the slash is interpreted as an object-name; the second field as an identifier to segi; the third as an identifier to chain, and so on. It may take any of the following forms: 

### Macros Beginning With a Slash
    
    
      
       /object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
       /object-name/segi-identifier/chain-identifier/resi-identifier
       /object-name/segi-identifier/chain-identifier
       /object-name/segi-identifier
       /object-name
    
    
    EXAMPLES
       PyMOL> zoom /pept
       PyMOL> show spheres, /pept/lig/
       PyMOL> show cartoon, /pept/lig/A
       PyMOL> color pink, /pept/lig/A/10
       PyMOL> color yellow, /pept/lig/A/10/CA
    

### Macros Not Beginning With a Slash

If the macro does not begin with a slash, it is interpreted differently. In this case, PyMOL expects to find the fields ending with the bottom of the hierarchy. Macros that don't start with a slash may take the following forms: 
    
    
                                                 resi-identifier/name-identifier
                                chain-identifier/resi-identifier/name-identifier
                segi-identifier/chain-identifier/resi-identifier/name-identifier
    object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
    
    
    EXAMPLES
       PyMOL> zoom 10/CB
       PyMOL> show spheres, A/10-12/CB
       PyMOL> show cartoon, lig/B/6+8/C+O
       PyMOL> color pink, pept/enz/C/3/N
    

### Omitting Fields in a Macro

You can also omit fields between slashes. Omitted fields will be interpreted as wildcards, as in the following forms: 
    
    
       
       resi-identifier/
       resi-identifier/name-identifier
       chain-identifier//
       object-name//chain-identifier                
    
       
    EXAMPLES
       PyMOL> zoom 142/                  # Residue 142 fills the viewer. 
       
       PyMOL> show spheres, 156/CA       # The alpha carbon of residue 156
                                         # is shown as a sphere     
                                         
       PyMOL> show cartoon, A//          # Chain "A" is shown as a cartoon.  
       
       PyMOL> color pink, /pept//B       # Chain "B" in object "pept"
                                         # is colored pink.
    

## See Also

  * [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")



Retrieved from "[https://pymolwiki.org/index.php?title=Selection_Macros&oldid=12699](https://pymolwiki.org/index.php?title=Selection_Macros&oldid=12699)"


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

## Set Color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Set_Color defines a new color with color indices (0.0-1.0). Numbers between 0 an 255 can be used as well. (If at least one value is larger than 1, pymol will interpret all 3 values as between 0 and 255). If an existing color name is used, the old color will be overridden. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 See Also



### USAGE
    
    
    set_color name, [ red-float, green-float, blue-float ]
    set_color name = [ red-float, green-float, blue-float ]  #(DEPRECATED)
    

### PYMOL API
    
    
    cmd.set_color( string name, float-list rgb )
    

### EXAMPLES
    
    
    PyMOL>set_color red, [1,0.01,0.01]
     Color: "red" defined as [ 1.000, 0.010, 0.010 ].
    PyMOL>set_color khaki, [195,176,145]
     Color: "khaki" defined as [ 0.765, 0.690, 0.569 ].
    

These will be added to the end of the list of Pymol's color indices that you can view the [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices") command. 

## See Also

  * [Get_Color_Tuples](/index.php/Get_Color_Tuples "Get Color Tuples")
  * [set_object_color](/index.php?title=Set_object_color&action=edit&redlink=1 "Set object color \(page does not exist\)")



Retrieved from "[https://pymolwiki.org/index.php?title=Set_Color&oldid=13184](https://pymolwiki.org/index.php?title=Set_Color&oldid=13184)"


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

## Spectrum states

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/spectrum_states.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrum_states.py)  
Author(s)  | [Takanori Nakane](/index.php?title=User:TakanoriNakane&action=edit&redlink=1 "User:TakanoriNakane \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.viewing](http://pymol.org/psicoredirect.php?psico.viewing)  
  
**spectrum_states** colors each state in a multi-state object different. Colors are along a user defined palette, just like with the [spectrumany](/index.php/Spectrumany "Spectrumany") command. 

## Usage
    
    
    spectrum_states [ selection [, representations [, color_list ]]]
    

## Example
    
    
    import spectrum_states
    
    # fetch a morph from molmovdb.org
    load http://molmovdb.org/uploads/284066-6299/movie.pdb.gz
    
    # show ribbon in all states
    as ribbon
    set all_states
    
    # color states from dark gray to red
    spectrum_states *, ribbon, gray20 orange red
    

## See Also

  * [spectrumany](/index.php/Spectrumany "Spectrumany")
  * [split_states](/index.php/Split_states "Split states")
  * [color_obj](/index.php/Color_Objects "Color Objects")



Retrieved from "[https://pymolwiki.org/index.php?title=Spectrum_states&oldid=13929](https://pymolwiki.org/index.php?title=Spectrum_states&oldid=13929)"


---

## Spectrumany

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/spectrumany.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrumany.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.viewing](http://pymol.org/psicoredirect.php?psico.viewing)  
  
[![](/images/8/8b/SpectrumanyExample.png)](/index.php/File:SpectrumanyExample.png)

[](/index.php/File:SpectrumanyExample.png "Enlarge")

Coloring with a gradient of different green shades

This script works similar to the [spectrum](/index.php/Spectrum "Spectrum") command, but instead of predefined palettes, any color sequence can be used. 

The color sequence is given by a space separated list of colors, so palette "red_white_blue" is the same as color sequence "red white blue". 

# Example
    
    
    fetch 2x19
    
    # these two produce the same result
    spectrum count, red_white_blue, chain B
    spectrumany count, red white blue, chain B
    
    # gradient of different green colors
    spectrumany count, smudge palegreen limegreen limon green forest, chain B
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spectrumany&oldid=13930](https://pymolwiki.org/index.php?title=Spectrumany&oldid=13930)"


---

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

## Surface color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 See Also



# Overview

**surface_color** controls the color of surfaces as drawn in PyMOL. 

  * [![Usual surface coloring](/images/b/bc/Surface_color_ex1.png)](/index.php/File:Surface_color_ex1.png "Usual surface coloring")

Usual surface coloring 

  * [![Color the entire surface "marine".](/images/d/d3/Surface_color_ex2.png)](/index.php/File:Surface_color_ex2.png "Color the entire surface "marine".")

Color the entire surface "marine". 




# Syntax
    
    
    # color the surface
    set surface_color, (color), (selection)
    

# Examples
    
    
    # color the surface white
    set surface_color, white, *
    
    # return surface coloring to the default scheme
    set surface_color, default, *
    

# See Also

[Color](/index.php/Color "Color"), [Color_Values](/index.php/Color_Values "Color Values")

Retrieved from "[https://pymolwiki.org/index.php?title=Surface_color&oldid=12614](https://pymolwiki.org/index.php?title=Surface_color&oldid=12614)"


---

