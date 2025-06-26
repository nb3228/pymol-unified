# Category: Objects and Selections

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

## Displaying Biochemical Properties

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Selecting secondary structures
    * 1.1 Manually Assigning Secondary Structure
    * 1.2 See Also
  * 2 Coloring
    * 2.1 Color by atom type from a script
    * 2.2 Assign color by B-factor
    * 2.3 Representation-independent color control
  * 3 Bonds
    * 3.1 Displaying double bonds
    * 3.2 Hydrogen bonds and Polar Contacts
      * 3.2.1 Hydrogen bonds between specific atoms
      * 3.2.2 Hydrogen bonds where find->polar contacts doesn't do what you need
      * 3.2.3 Details
  * 4 Calculating dihedral angles
  * 5 Cavities
  * 6 Surface-Related
    * 6.1 Surface Area
    * 6.2 Polar surface area
    * 6.3 Display solvent accessible surface
    * 6.4 Contact Potential
    * 6.5 Electric Field Lines
    * 6.6 Residues with functional groups
  * 7 Backbones
    * 7.1 Displaying the C-Alpha trace of proteins
    * 7.2 Displaying the Amino Acid Backbone
    * 7.3 Displaying the Phosphate backbone of nucleic acids
      * 7.3.1 Native Nucleic Acid Rendering in PyMol
  * 8 Align proteins with CA fit



## Selecting secondary structures

A few examples: 
    
    
    select helix, (ss h)
    select sheet, (ss s)
    select loop, (ss l+'')
    

### Manually Assigning Secondary Structure

You can manually assign secondary stuctures to your protein by 
    
    
    alter 96-103/, ss='S'
    alter 96-103/, ss='H'
    alter 96-103/, ss='L'
    

to set residues 96-103 to beta Strand, alpha Helix, and Loop respectively. 

### See Also

[DSSP](/index.php/DSSP "DSSP") [Dss](/index.php/Dss "Dss") [Caver](/index.php/Caver "Caver")

[FAQ](/index.php/Category:FAQ "Category:FAQ") [Displaying Biochemical Properties](/index.php/Category:Objects_and_Selections "Category:Objects and Selections")

## Coloring

See also [Category:Coloring](/index.php/Category:Coloring "Category:Coloring"). 

### Color by atom type from a script

See [Color](/index.php/Color "Color") for this. 

### Assign color by B-factor

See section [Color](/index.php/Color "Color") for this. 

### Representation-independent color control

See section [Surface#Representation-independent Color Control](/index.php/Surface#Representation-independent_Color_Control "Surface"). 

## Bonds

PyMOL can deduce bonds from the PDB structure file, even if the **CONECT** records are missing. In fact, PyMOL guesses bonding connectivity based on proximity, based on the empirical observation that two atoms of a given radius will not be generally closer than a certain distance unless they are bonded. 

### Displaying double bonds

  * [![Image showing double bonds in PyMOL. Double bonds are supported in lines and sticks.](/images/d/d2/DoubleBonds.png)](/index.php/File:DoubleBonds.png "Image showing double bonds in PyMOL. Double bonds are supported in lines and sticks.")

Image showing double bonds in PyMOL. Double bonds are supported in [lines](/index.php/Lines "Lines") and [sticks](/index.php/Sticks "Sticks"). 




You can go into the [lines](/index.php/Lines "Lines") mode and turning on the valence display: 
    
    
    hide
    as lines
    set valence, 0.1
    

A higher value for valence spreads things out more. See the [Valence](/index.php/Valence "Valence") page for more information on options for changing the appearance of the valence lines. 

### Hydrogen bonds and Polar Contacts

[![](/images/a/ad/Polar_contacts_small.png)](/index.php/File:Polar_contacts_small.png)

[](/index.php/File:Polar_contacts_small.png "Enlarge")

Polar Contacts in PyMol

Using the actions [A] button for an object or selection you can display Hydrogen bonds and Polar Contacts. [A]->find->polar contacts-><select from menu>

The command behind the menus is the **dist** ance command called with the additional argument mode=2. 

Parameters that control the the identification of H-bonds are defined as 
    
    
    set h_bond_cutoff_center, 3.6
    

with ideal geometry and 
    
    
    set h_bond_cutoff_edge, 3.2
    

with minimally acceptable geometry. 

These settings can be changed *before* running the detection process (dist command mode=2 or via the menus). 

Note that the hydrogen bond geometric criteria used in PyMOL was designed to emulate that used by [DSSP](http://swift.cmbi.kun.nl/gv/dssp/). 

#### Hydrogen bonds between specific atoms
    
    
    dist name, sele1, sele2, mode=2
    

  


#### Hydrogen bonds where find->polar contacts doesn't do what you need

You can show H-bonds between two objects using atom selections so long as hydrogens are present in both molecules. If you don't have hydrogens, you can use [h_add](/index.php/H_add "H add") on the proteins, or provide ligands with valence information and then use h_add. 

Two examples are below. For clarity, they draw dashes between the heavy atoms and hide the hydrogens. 
    
    
    # EXAMPLE 1: Show hydrogen bonds between protein 
    # and docked ligands (which must have hydrogens)
    
    load target.pdb,prot
    load docked_ligs.sdf,lig
    
    # add hydrogens to protein
    
    h_add prot
    
    select don, (elem n,o and (neighbor hydro))
    select acc, (elem o or (elem n and not (neighbor hydro)))
    dist HBA, (lig and acc),(prot and don), 3.2
    dist HBD, (lig and don),(prot and acc), 3.2
    delete don
    delete acc
    hide (hydro)
    
    hide labels,HBA
    hide labels,HBD
    
    
    
    # EXAMPLE 2
    # Show hydrogen bonds between two proteins
    
    load prot1.pdb
    load prot2.pdb
    
    h_add prot1
    h_add prot2
    
    select don, (elem n,o and (neighbor hydro))
    select acc, (elem o or (elem n and not (neighbor hydro)))
    dist HBA, (prot1 and acc),(prot2 and don), 3.2
    dist HBD, (prot1 and don),(prot2 and acc), 3.2
    delete don
    delete acc
    hide (hydro)
    
    hide labels,HBA
    hide labels,HBD
    
    # NOTE: that you could also use this approach between two
    # non-overlapping selections within a single object.
    

The "polar contacts" mentioned above are probably better at finding hydrogen bonds than these scripts. "Polar contacts" check geometry as well as distance. 

#### Details

Generally speaking, PyMOL does not have sufficient information to rigorously determine hydrogen bonds, since typical PDB file are ambiguous with respect to charge states, bonds, bond valences, and tautomers. As it stands, all of those things are guessed heuristically. Rigorously determining the location of lone pair electrons and proton coordinates from raw PDB files is a nontrival problem especially when arbitrary small molecule structures are present. In addition, PyMOL would also need to consider the implied coordinate error due to overall structure resolution and local temperature factors before rigorously asserting than any specific hydrogen bond does or does not exist. 

Furthermore, our hydrogen bond detection machinery was originally developed for purposes of mimicking Kabsch and Sander's DSSP secondary structure assignment algorithm (Biopolymers 22, 2577, 1983) which is based on a rather generous notion of hydrogen bonding (see Kabsch Figure 1). 

Although this approximate capability can be accessed via the distance command using mode=2, the criteria applied by our implementation may be based on heavy-atom coordinates (only) and does not necessarily correspond to anything rigorous or published. So the bottom line is that PyMOL merely offers up putative polar contacts and leaves it to the user to determine whether or not the interactions present are in fact hydrogen bonds, salt bridges, polar interactions, or merely artifacts of incorrect assignments (i.e. two carbonyls hydrogen bonding because they're being treated like hydroxyls). 

With respect to the h_bond_* settings, the angle in question for h_bond_cutoff_* and h_bond_max_angle is ADH, assuming H exists. If H does not exist, then PyMOL will guess a hypothetical coordinate which may not actually be valid (in plane, etc.). Tthe hydrogen must also lie within a cone of space with its origin on A (along B->A) and angular width h_bond_cone. Since h_bond_cone in 180 by default, the present behavior is to simply reject any hydrogen bond where the hydrogen would lie behind the plane defined by the acceptor atom (A) in relation to its bonded atom(s) B (if any). In other words, if B is only one atom (e.g. C=O vs. C-O-C), then by default, HAB cannot be less then 90 degrees. 

The two h_bond_power_* settings are merely fitting parameters which enable PyMOL to reproduce a curve shape reflecting Kabsch Figure 1. The endpoints of the effective cutoff curve is a function of the two h_bond_cutoff_* setting. 

## Calculating dihedral angles

The get_dihedral function requires four single-atom selections to work: 
    
    
    get_dihedral prot1///9/C, prot1///10/N, prot1///10/CA, prot1///10/C
    

## Cavities

See [Surfaces_and_Voids](/index.php/Surfaces_and_Voids "Surfaces and Voids"). Also [Caver](/index.php/Caver "Caver") and [CASTp](/index.php/CASTp "CASTp"). 

## Surface-Related

### Surface Area

To calculate the surface area of a selection, see [Get_Area](/index.php/Get_Area "Get Area"). 

### Polar surface area

For a solvent accessible PSA approximation: 
    
    
    set dot_density, 3
    remove hydro
    remove solvent
    show dots
    set dot_solvent, on
    get_area elem N+O
    get_area elem C+S
    get_area all
    

For molecular PSA approximation 
    
    
    set dot_density, 3
    remove hydro
    remove solvent
    set dot_solvent, off
    get_area elem N+O
    get_area elem C+S
    get_area all
    

Showing dots isn't mandatory, but it's a good idea to confirm that you're getting the value for the atom dot surface you think you're using. Please realize that the resulting numbers are only approximate, reflecting the sum of partial surface areas for all the dots you see. To increase accuracy, set dot_density to 4, but be prepared to wait... 

### Display solvent accessible surface

Using the surface display mode, PyMOL doesn't show the solvent accessible surface, rather it shows the solvent/protein contact surface. The solvent accessible surface area is usually defined as the surface traced out by the center of a water sphere, having a radius of about 1.4 angstroms, rolled over the protein atoms. The contact surface is the surface traced out by the vdw surfaces of the water atoms when in contact with the protein. 

PyMOL can show solvent accessible surfaces using the dot or sphere representations: 

for dots: 
    
    
    show dots
    set dot_mode,1
    set dot_density,3
    

for spheres: 
    
    
    alter all,vdw=vdw+1.4
    show spheres
    

Once the Van der Waals radii for the selection have been altered, the surface representation will also be "probe-inflated" to show a pseudo solvent accessible surface, as detailed above. 

for surfaces: 
    
    
    alter all,vdw=vdw+1.4
    show surface
    

[![](/images/d/d0/Solvent-accessible_surface.jpg)](/index.php/File:Solvent-accessible_surface.jpg)

[](/index.php/File:Solvent-accessible_surface.jpg "Enlarge")

Solvent-Accessible Surface Example

  


  


  


  


  


  


  


  


  
Note that to display both the molecular surface and the solvent-accessible surface, the object must be duplicated, as is done for [Surface#Representation-independent Color Control](/index.php/Surface#Representation-independent_Color_Control "Surface"). This also applies if the spheres representation is to be used to display "real" atoms. 

### Contact Potential

See [Protein_contact_potential](/index.php/Protein_contact_potential "Protein contact potential") and [APBS](/index.php/APBS "APBS"). 

[![Prot contact pot.png](/images/6/64/Prot_contact_pot.png)](/index.php/File:Prot_contact_pot.png)

[](/index.php/File:Prot_contact_pot.png "Enlarge")

### Electric Field Lines

[![](/images/6/66/Elines.png)](/index.php/File:Elines.png)

[](/index.php/File:Elines.png "Enlarge")

PyMOL and APBS used to show electronic field lines.

To produce an image with electric field lines, first run APBS. Then, input the following: 
    
    
    gradient my_grad, pymol-generated
    ramp_new my_grad_ramp, pymol-generated
    color my_grad_ramp, my_grad
    

### Residues with functional groups

[![](/images/c/c4/1igt_cys_lys_asp_glu_colored.png)](/index.php/File:1igt_cys_lys_asp_glu_colored.png)

[](/index.php/File:1igt_cys_lys_asp_glu_colored.png "Enlarge")

Whole residues colored (Cys: yellow, Lys: blue, Asp and Glu: red)

Poor man's solution: Display protein as surface, colorize all Lys (-NH2), Asp and Glu (-COOH) and Cys (-SH): 
    
    
    remove resn hoh    # remove water
    h_add              # add hydrogens
    
    as surface
    color grey90
    
    color slate, resn lys       # lysines in light blue
    color paleyellow, resn cys  # cysteines in light yellow
    color tv_red, (resn asp or(resn glu))  # aspartic and glutamic acid in light red
    

[![](/images/3/33/1igt_functional_groups_colored.png)](/index.php/File:1igt_functional_groups_colored.png)

[](/index.php/File:1igt_functional_groups_colored.png "Enlarge")

Only central atoms of functional groups colored (Cys: S, Lys: NH2, Asp and Glu: CO2)

Not-_so_ -poor-man's solution: In order to have the functional groups better localized, only the central atoms can be colored: 

  * the S atom of cystein,
  * the N and H atoms of the free amine of lysine (may be displayed with three H atoms at all three possible positions)
  * the C and two O atoms of free carboxylic groups in aspartic and glutamic acid



In this way, they are better visible through the surface compared to only one colored atom, both amines and carboxylic groups consist of three colored atoms each. 
    
    
    remove resn hoh    # remove water
    h_add              # add hydrogens
    
    as surface
    color grey90
    
    select sulf_cys, (resn cys and (elem S))      # get the sulfur atom of cystein residues
    color yellow, sulf_cys
    
    select nitro_lys, (resn lys and name NZ)              # get the nitrogens of free amines ("NZ" in PDB file)
    select hydro_lys, (elem H and (neighbor nitro_lys))   # get the neighboring H atoms 
    select amine_lys, (nitro_lys or hydro_lys)
    color tv_blue, amine_lys
    
    
    select oxy_asp, (resn asp and (name OD1 or name OD2))             # get the two oxygens of -COOH  ("OD1", "OD2")
    select carb_asp, (resn asp and (elem C and (neighbor oxy_asp)))   # get the connecting C atom
    select oxy_glu, (resn glu and (name OE1 or name OE2))             # oxygens "OE1" and "OE2" in PDB file
    select carb_glu, (resn glu and (elem c and (neighbor oxy_glu)))
    select carboxy, (carb_asp or oxy_asp or carb_glu or oxy_glu)
    color tv_red, carboxy
    

By displaying the protein as non-transparent surface, only the functional groups (colored atoms) at the surface are visible. The visualization of those groups can be pronounced by displaying the corresponding atoms as spheres, e.g. "as spheres, carboxy + amine_lys + sulf_cys", in this way it might become more clear how accessible they are. 

When displaying the protein as cartoon, the functional groups can be shown as spheres, and the whole residues cys, lys, asp and glu as sticks connected to the backbone, with the atoms of the functional groups again as spheres. However, then also the not accessible residues inside the protein are visible. 

## Backbones

### Displaying the C-Alpha trace of proteins
    
    
    hide
    show ribbon
    set ribbon_sampling,1
    

And if your model only contains CA atoms, you'll also need to issue: 
    
    
    set ribbon_trace,1
    

### Displaying the Amino Acid Backbone

The easiest way to see the backbone of the protein is to do 
    
    
    hide all
    show ribbon
    

If you don't like the ribbon representation, you can also do something like 
    
    
    hide all
    show sticks, name C+O+N+CA
    

You can replace **sticks** in the above by other representations like **spheres** or **lines**. 

### Displaying the Phosphate backbone of nucleic acids

#### Native Nucleic Acid Rendering in PyMol

PyMol now better supports viewing nucleic acid structure. [Nuccyl](/index.php/Nuccyl "Nuccyl") still seems to be the reigning champ for image quality, but see PyMol's native [Cartoon](/index.php/Cartoon "Cartoon") command. For more information on representing nucleic acids, please see the [Nucleic Acids](/index.php/Category:Nucleic_Acids "Category:Nucleic Acids") Category. 

[![Cnam 0.png](/images/0/0c/Cnam_0.png)](/index.php/File:Cnam_0.png)

[](/index.php/File:Cnam_0.png "Enlarge")

Should you ever want to show the phosphate trace of a nucleic acid molecule: 
    
    
    def p_trace(selection="(all)"):
        s = str(selection)
        cmd.hide('lines',"("+s+")")
        cmd.hide('spheres',"("+s+")")
        cmd.hide('sticks',"("+s+")")
        cmd.hide('ribbon',"("+s+")")
        cmd.show('cartoon',"("+s+")")
        cmd.set('cartoon_sampling',1,"("+s+")")
        cmd.set('cartoon_tube_radius',0.5,"("+s+")")
    cmd.extend('p_trace',p_trace)
    

and then: 
    
    
    p_trace (selection)
    

## Align proteins with CA fit

If two proteins have significant homology, you can use the [Align](/index.php/Align "Align") command: 
    
    
    align prot1////ca,prot2
    

which will perform a sequence alignment of prot1 against prot2, and then an optimizing fit using the CA positions. I'm not sure if the help text for align got into 0.82, but the next version will definitely have it. 

Retrieved from "[https://pymolwiki.org/index.php?title=Displaying_Biochemical_Properties&oldid=11686](https://pymolwiki.org/index.php?title=Displaying_Biochemical_Properties&oldid=11686)"


---

## Symmetry

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

There are many ways that symmetry can be important/useful/beautiful to look at in macromolecular structures. A typically application is the reconstruction of a symmetric oligomer from a few subunits. PyMOL has tools that can help with this type of analysis or depiction. 

  


## Expanding crystallographic symmetry

Structures determined by X-ray crystallography are typically deposited as files containing the coordinates for one asymmetric unit (ASU). Knowledge of the symmetry operators that describe how the ASUs are arranged relative to each other allows the arrangement of the crystal lattice to be recreated. PyMOL can read this symmetry information from the input coordinate file and recreate the neigbouring copies of the ASU using symmexp. 

PyMOL's built in symmetry expansion functionality is available as A->generate->symmetry mates for an object or as the [symexp](/index.php/Symexp "Symexp") command. 

The [SuperSym](/index.php/SuperSym "SuperSym") plugin has additional unit cell and symmetry axis tools. 

## Displaying symmetry axes

Often it is necessary to be able to find and draw symmetry axes. There are some contributed scripts that help do this including [Symmetry Axis](/index.php/Symmetry_Axis "Symmetry Axis"), for which you need to know the coordinates and direction of the axis you would like to draw, and [RotationAxis](/index.php/RotationAxis "RotationAxis"), which does an excellent job of displaying symmetry relationships between selections. 

Retrieved from "[https://pymolwiki.org/index.php?title=Symmetry&oldid=13341](https://pymolwiki.org/index.php?title=Symmetry&oldid=13341)"


---

## Nonstandard Amino Acids

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This page talks a little about how PyMOL deals with nonstandard amino acids, and the various representations and options available. 

By default, PyMOL considers [nonstandard amino acids](http://en.wikipedia.org/wiki/Amino_acid#Nonstandard_amino_acids) as HETERO atoms. Therefore, when you draw a default [surface](/index.php/Surface "Surface"), the heteroatoms are not included. Also, nonstandard atoms are left out of the backbone representation in the [cartoon](/index.php/Cartoon "Cartoon") and [ribbon](/index.php/Ribbon "Ribbon") representations. 

### Workarounds

  * If you're looking to represent the backbone via [ribbon](/index.php/Ribbon "Ribbon") or [cartoon](/index.php/Cartoon "Cartoon"), then just use the [mutagenesis](/index.php/Mutagenesis "Mutagenesis") to modify the nonstandard amino acid to something standard, like GLU, or LEU. The [draw](/index.php/Draw "Draw")/[ray](/index.php/Ray "Ray") your structure as needed.


  * If you want the nonstandard amino acid to be included in the [surface representation](/index.php/Surface "Surface"), then [set surface_mode](/index.php/Surface_mode "Surface mode") to 1. For example, consider the images below. We have a hetero atom shown in the image at left, not included in the surface, become part of the surface when we set the [surface_mode](/index.php/Surface_mode "Surface mode") to 1.


  * Surface Mode Examples
  * [![surface_mode set to 0, the default. The galactose \(blue\) is not considered part of the surface.](/images/6/6c/Sm0.png)](/index.php/File:Sm0.png "surface_mode set to 0, the default. The galactose \(blue\) is not considered part of the surface.")

[surface_mode](/index.php/Surface_mode "Surface mode") set to 0, the default. The galactose (blue) is not considered part of the surface. 

  * [![surface_mode set to 1 -- now including heteroatoms. The galactose and all heteroatoms \(blue\) are now considered part of the surface and colored blue.](/images/4/46/Sm1.png)](/index.php/File:Sm1.png "surface_mode set to 1 -- now including heteroatoms. The galactose and all heteroatoms \(blue\) are now considered part of the surface and colored blue.")

[surface_mode](/index.php/Surface_mode "Surface mode") set to 1 -- now including heteroatoms. The galactose and all heteroatoms (blue) are now considered part of the surface and colored blue. 



  * Another work around is to try


    
    
    flag ignore, bymol polymer, set
    rebuild
    

Retrieved from "[https://pymolwiki.org/index.php?title=Nonstandard_Amino_Acids&oldid=5906](https://pymolwiki.org/index.php?title=Nonstandard_Amino_Acids&oldid=5906)"


---

## Category:Representations

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Representations

Representations are how PyMOL draws a selection or object. The choices available are shown below. 

## Pages in category "Representations"

The following 15 pages are in this category, out of 15 total. 

### B

  * [Ball and Stick](/index.php/Ball_and_Stick "Ball and Stick")



### C

  * [Cartoon](/index.php/Cartoon "Cartoon")



### D

  * [Dots](/index.php/Dots "Dots")



### E

  * [Ellipsoids](/index.php/Ellipsoids "Ellipsoids")
  * [Examples of nucleic acid cartoons](/index.php/Examples_of_nucleic_acid_cartoons "Examples of nucleic acid cartoons")



### L

  * [Light](/index.php/Light "Light")
  * [Lines](/index.php/Lines "Lines")



### M

  * [Mesh](/index.php/Mesh "Mesh")



### R

  * [Ribbon](/index.php/Ribbon "Ribbon")



### S

  * [Spheres](/index.php/Spheres "Spheres")
  * [Sticks](/index.php/Sticks "Sticks")
  * [Surface](/index.php/Surface "Surface")
  * [Main Page](/index.php/Main_Page "Main Page")
  * [Surface type](/index.php/Surface_type "Surface type")



### V

  * [Volume](/index.php/Volume "Volume")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Representations&oldid=4866](https://pymolwiki.org/index.php?title=Category:Representations&oldid=4866)"


---

## Category:Selector Quick Reference

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page contains links to all of the various ways one can make selections in PyMOL. 

We need to outline the importance of good selection capabilities here, as well. 

## Pages in category "Selector Quick Reference"

The following 4 pages are in this category, out of 4 total. 

### P

  * [Property Selectors](/index.php/Property_Selectors "Property Selectors")



### S

  * [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")
  * [Selection Macros](/index.php/Selection_Macros "Selection Macros")
  * [Single-word Selectors](/index.php/Single-word_Selectors "Single-word Selectors")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Selector_Quick_Reference&oldid=2819](https://pymolwiki.org/index.php?title=Category:Selector_Quick_Reference&oldid=2819)"


---

## Category:Working with Objects

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Working_with_Objects&oldid=2945](https://pymolwiki.org/index.php?title=Category:Working_with_Objects&oldid=2945)"


---

## Category:Working with Selections

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

To learn about making PyMol selections, see 

  1. (E) [ Selection Expression](/index.php/Selection-expressions "Selection-expressions")
  2. (N) [ Named Atom Selections](/index.php/Named_Atom_Selections "Named Atom Selections")



## Pages in category "Working with Selections"

The following 2 pages are in this category, out of 2 total. 

### N

  * [Named Atom Selections](/index.php/Named_Atom_Selections "Named Atom Selections")



### S

  * [Selection-expressions](/index.php/Selection-expressions "Selection-expressions")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Working_with_Selections&oldid=2944](https://pymolwiki.org/index.php?title=Category:Working_with_Selections&oldid=2944)"


---

