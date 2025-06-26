# Category: Biochemical Properties

## APBS

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/f/f8/Rna_surface_apbs.png)](/index.php/File:Rna_surface_apbs.png)

[](/index.php/File:Rna_surface_apbs.png "Enlarge")

APBS-generated electrostatic surface displayed in PyMOL

[APBS](http://www.poissonboltzmann.org/), the Adaptive Poisson-Boltzmann Solver, is a [freely](http://www.oreilly.com/openbook/freedom/) available macromolecular electrostatics calculation program released under a [BSD license](http://apbs-pdb2pqr.readthedocs.io/en/latest/apbs/license.html). PyMOL can display the results of the calculations as an electrostatic potential molecular surface. 

## Contents

  * 1 APBS Plugins for PyMOL
  * 2 Required Dependencies
  * 3 Troubleshooting
  * 4 Using APBS



## APBS Plugins for PyMOL

  * [APBS Electrostatics Plugin](/index.php/APBS_Electrostatics_Plugin "APBS Electrostatics Plugin"), included in [Incentive PyMOL 2.0](https://pymol.org/)
  * [APBS Tools 2.1](/index.php/Apbsplugin "Apbsplugin"), based on the original version by [Michael Lerner](/index.php/User:Mglerner "User:Mglerner").



Both plugins make it possible to run APBS from within PyMOL, and then display the results as a color-coded electrostatic surface (units  K b T / e c {\displaystyle K_{b}T/e_{c}} ![{\\displaystyle K_{b}T/e_{c}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/f6c49dfd8965c68897680da7541c3332a8d04bff)) in the molecular display window (as with the image to the right). 

## Required Dependencies

The plugins require **apbs** and **pdb2pqr**. 

  * [Incentive PyMOL](https://pymol.org) ships preconfigured with apbs and pdb2pqr
  * [Official release downloads](https://github.com/Electrostatics/apbs-pdb2pqr/releases) for all platforms
  * Precompiled packages are also available in Ubuntu, Debian, Gentoo, [MacPorts](https://www.macports.org/), [Fink](http://www.finkproject.org/), and many other Unix-like distributions. Example for installation on Ubuntu:


    
    
    apt-get install apbs pdb2pqr
    

After all components are installed, open the plugin and browse for apbs and pdb2pqr on the "Program Locations" tab (if they haven't been found automatically). 

## Troubleshooting

  * _(this might be outdated information):_ If the B-factor is  ≥ 100 , {\displaystyle \geq 100,} ![{\\displaystyle \\geq 100,}](https://wikimedia.org/api/rest_v1/media/math/render/svg/acff937e615aab509cd3154946a28c26fdeb3ee2) then APBS doesn't properly read in the PDB file and thus outputs garbage (or dies). To fix this, set all b factors to be less than 100. 
        
        alter all, b=min(b,99.9)
        

The problem stems from how to parse a PDB file. The PDB file originally was written when most people used FORTRAN programs, and so the file format was specified by columns, not by the more modern comma separated value format we tend to prefer today. For the latest on the PDB format see the [new PDB format docs](http://www.wwpdb.org/docs.html).
  * APBS has problems, sometimes, in reading atoms with **alternate conformations**. You can remove the alternate locations with a simple script [removeAlt](/index.php/RemoveAlt "RemoveAlt").
  * For pdb2pqr, **RNA** resdiue names must be RA, RC, RG, and RU.
        
        alter polymer & resn A+C+G+U, resn = "R" + resn
        

  * Incomplete Residues: Some truncated PDB files include a single backbone atom of the next residue, e.g. [2xwu](https://www.rcsb.org/structure/2xwu) chain B residue 954 atom N. **pdb2pqr** reports: _Error encountered: Too few atoms present to reconstruct or cap residue LEU B 954 in structure!_. The easiest solution is to remove that atom:
        
        remove /2xwu//B/954
        




[![](/images/3/37/Apbs_ex.png)](/index.php/File:Apbs_ex.png)

[](/index.php/File:Apbs_ex.png "Enlarge")

PyMOL visualizing two maps at once

## Using APBS

See [APBS Electrostatics Plugin](/index.php/APBS_Electrostatics_Plugin "APBS Electrostatics Plugin") and [APBS Tools2.1](/index.php/Apbsplugin "Apbsplugin"). 

Retrieved from "[https://pymolwiki.org/index.php?title=APBS&oldid=12790](https://pymolwiki.org/index.php?title=APBS&oldid=12790)"


---

## Caver

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/loschmidt/caver-pymol-plugin/archive/master.zip>  
Author(s)  | CAVER Team   
License  | [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)  
<http://www.caver.cz>  
  
  
The Caver3 plugin is the successor of [Caver2](/index.php/Caver2 "Caver2"). CAVER is a software tool for analysis and visualization of tunnels and channels in protein structures. 

## Installation

  1. Install [Java](https://java.com/en/download/)
  2. Install [caver-pymol-plugin-master.zip](https://github.com/loschmidt/caver-pymol-plugin/archive/master.zip) using the PyMOL [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



[![PluginManagerCaverPlugin301 fixed.png](/images/b/b2/PluginManagerCaverPlugin301_fixed.png)](/index.php/File:PluginManagerCaverPlugin301_fixed.png)

## See Also

  * [Mole2](/index.php/Mole2 "Mole2")
  * [Caver2](/index.php/Caver2 "Caver2") (old version)



Retrieved from "[https://pymolwiki.org/index.php?title=Caver3&oldid=13180](https://pymolwiki.org/index.php?title=Caver3&oldid=13180)"


---

## CavitOmiX

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://innophore.com/cavitomix>  
Author(s)  | [Georg Steinkellner](https://www.linkedin.com/in/georg-steinkellner-81b2bbb9), [Christian C. Gruber](https://www.linkedin.com/in/christiancgruber), [Karl Gruber](https://www.linkedin.com/in/karl-gruber-2a486135)*, and the [Innophore Team](https://www.linkedin.com/company/innophore)  
License  |   
<https://innophore.com>  
  
  
[![Cavitomix.png](/images/8/8f/Cavitomix.png)](/index.php/File:Cavitomix.png)

## CavitOmiX 1.0

CavitOmiX plugin for Schrodinger’s PyMOL, a nifty tool that allows you to analyze protein cavities from any input structure. You can now dive deep into your proteins, cavities, and binding sites using crystal structures and state-of-the-art AI models from OpenFold (powered by NVIDIA’s BioNeMo service), DeepMind`s AlphaFold and ESMFold by Meta. Even more exciting: just enter any protein sequence and you will get the structure predicted by OpenFold or ESMFold loaded into your PyMOL within seconds. 

  
[CavitOmiX](https://innophore.com/cavitomix). 

  * Catalophore™ Cavities can be calculated for molecules loaded in PyMOL
  * Predict protein structures within seconds for any protein sequence using [OpenFold](https://www.nvidia.com/en-us/gpu-cloud/bionemo) by NVIDIA [BioNeMo](https://www.nvidia.com/en-us/gpu-cloud/bionemo) (coming soon!) and [ESMFold by Meta](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/)
  * [AlphaFold](https://www.deepmind.com/research/highlighted-research/alphafold) models can be retrieved via [UniProt ID](https://www.uniprot.org/)
  * Analyze the hydrophobicity of your Catalophore™ cavities
  * Protein structures can be retrieved from the [PDB](https://www.rcsb.org/) using the PDB code
  * Mix and match all the above in a single entry, align the structures and get a quick overview
  * Catalophore™ cavity calculation settings can be changed
  * Each [Catalophore™ cavity](https://innophore.com) is an "residue" entry and each cavity point is an "atom", so you can select, remove, copy, represent cavities to your liking!



  


## Overview

  *   * 


  


[![](/images/3/38/Cavitomix-plugin-V1.0_screenshot.png)](/index.php/File:Cavitomix-plugin-V1.0_screenshot.png)

[](/index.php/File:Cavitomix-plugin-V1.0_screenshot.png "Enlarge")

CavitOmiX PyMOL plugin with BioNeMo by NVIDIA

## References & Information

  1. [Innophore GmbH](https://innophore.com)
  2. [NVIDIA BioNeMo](https://www.nvidia.com/en-us/gpu-cloud/bionemo)
  3. [ESMFold by Meta](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/)
  4. [AlphaFold2](https://alphafold.ebi.ac.uk)



Retrieved from "[https://pymolwiki.org/index.php?title=CavitOmiX&oldid=13520](https://pymolwiki.org/index.php?title=CavitOmiX&oldid=13520)"


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

## Distancetoatom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/distancetoatom.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/distancetoatom.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar") and [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![get distance to atoms within cutoff](/images/c/c2/Distancetoatom_0.png)](/index.php/File:Distancetoatom_0.png "get distance to atoms within cutoff")

## Contents

  * 1 About distancetoatom
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
  * 4 SEE ALSO



## About distancetoatom

**distancetoatom** prints all distances between a specified atom, coordinate or group selection center and all atoms within cutoff distance that are part of the selection.  
All coordinates and distances can be saved in a csv-style text file report and/or can be stored in a (custom) atom property, if defined. 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage
    
    
        distancetoatom [ origin [, cutoff [, filename [, selection [, state [, property_name [, coordinates [, decimals [, sort [, quiet ]]]]]]]]]]
    

|   
---|---  
  
### Arguments

Arguments for distancetoatom   
---  
Keyword  | Default  | What it does   
origin  | pk1  | "origin" defines the coordinates for the origin and can be: 1\. a list with coordinates [x,y,z]  
2\. a single atom selection string  
3\. a multi-atom selection string (center will be used)   
cutoff  | 10  | "cutoff" defines the maximum distance, atoms further away are not considered   
filename  | None  | "filename" is the name of a report file that will be created (optional). set to e.g. 'report.txt' to create a report. This file is CSV style and can be imported into EXCEL.  
(omit or set to "", None, 0 or False to disable)   
selection  | all  | only atoms within "selection" (and cutoff distance) will be used for calculation   
state  | 0  | "state" is the state that will be used for calculations use 0 (omit), to automatically use the current state   
property_name  | p.dist  | "property_name" defines a (custom) property in which the calculated distance will be stored. This can be used e.g. for labeling or spectrum coloring etc. (cf. examples)  
NB! older PyMOL versions only support b or q and the distance will overwrite this property if set.   
coordinates  | 0  | "coordinates" toggles whether, besides distance, the atom coordinates will be included in the report.   
decimals  | 3  | "decimals" is the number of decimal places calculated. Note that PDB files do not support a higher resolution than 3.   
sort  | 1  | "sort" defines how the output will be sorted: 1: ascending (default)  
0: no sorting (by selection)  
-1: descending   
quiet  | 1  | toggle verbosity   
  
  


## Examples
    
    
    # make function available to PyMOL
    import distancetoatom
    
    ##############################
    # Example 1: print all distances
    frag LYS
    edit name CA
    distancetoatom
    
    ##############################
    # Example 2: print all distances to file
    fetch 1cll, async=0
    distancetoatom origin=/1cll//A/CA`150/CA, filename=distances.txt, selection=elem O, property_name=b
    # the file is CSV-format and can be imported, e.g. to EXCEL
    
    # format
    hide everything
    select byres (resi 150 expand 5 and not resn hoh) 
    show_as sticks, sele
    util.cbaw sele
    show_as spheres, resi 150
    zoom sele
    
    ##########
    # Label by stored distance (property_name)
    label sele and elem O, "%.2f"%b
    
    
    ##############################
    # Example 3: color an object by distance to a coordinate
    
    fetch 1cll, async=0
    distancetoatom [0,0,0], cutoff=1000, property_name=b
    # the distance is stored in the b-factor in this case
    # newer PyMOL versions support custom properties, e.g. "p.dist"
    spectrum b, rainbow_rev
    

|  [![Example 2](/images/5/5c/Distancetoatom_2.png)](/index.php/File:Distancetoatom_2.png "Example 2")  
  
[![Example 3](/images/3/3a/Distancetoatom_1.png)](/index.php/File:Distancetoatom_1.png "Example 3")  
---|---  
  
  


## SEE ALSO

  * [Distance](/index.php/Distance "Distance")
  * [Get Distance](/index.php/Get_Distance "Get Distance")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")



Retrieved from "[https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939](https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939)"


---

## Dss

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**dss** _d_ efines _s_ econdary _s_ tructure based on backbone geometry and hydrogen bonding patterns. 

With PyMOL, heavy emphasis is placed on cartoon aesthetics, and so both hydrogen bonding patterns and backbone geometry are used in the assignment process. Depending upon the local context, helix and strand assignments are made based on geometry, hydrogen bonding, or both. 

This command will generate results which differ slightly from DSSP and other programs. Most deviations occur in borderline or transition regions. Generally speaking, PyMOL is more strict, thus assigning fewer helix/sheet residues, except for partially distorted helices, which PyMOL tends to tolerate. 

WARNING: This algorithm has not yet been rigorously validated. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES
  * 4 NOTES
    * 4.1 Secondary Structure Determination



### USAGE
    
    
    dss [selection [, state]]
    

### ARGUMENTS

  * selection = string: atom selection {default: (all)}
  * state = integer: state-index if positive number or any of these: 
    * state = 0: consensus over all states {default}
    * state = -1: current state
    * state = -2: consensus over effective object states (implemented?)
    * state = -3: no consensus, residue gets assignment from earliest state that is not 'L'
    * state = -4: consensus over first and last state



### EXAMPLES
    
    
    # determine secondary structures in 
    # all loaded objects in PyMOL
    dss
    

### NOTES

If you dislike one or more of the assignments made by dss, you can use the alter command to make changes (followed by "rebuild"). For example: 
    
    
    # set residues 123-125 as being loops
    alter 123-125/, ss='L'
    
    # set the secondary structure of selection (pk1) to beta sheet
    alter pk1, ss='S'
    
    # set residue 90 to be alpha-helical
    alter 90/, ss='H'
    
    # update the scene in PyMOL to reflect the changes.
    rebuild
    

#### Secondary Structure Determination

As is typical with PyMOL, the secondary structure assignment engine is ad hoc and empirically tuned to produce desirable aesthetics. Though there are some phi/psi's that are clearly helix/sheet and others that are clearly not, there are certain regions of phi/psi space were the assignment is subjective or arbitrary. In my experience, algorithms based on strict definitions tend to operate poorly in such regions, and so PyMOL's algorithm is "fuzzy" in that there is a grey area where residues may be accepted or rejected as helix/sheet depending upon the surrounding context. 

There aren't any hard & fast definitions. But you are welcome to check out the collection of settings beginning with "ss_helix" and "ss_strand", noting that the include and exclude settings are deviations around the target in degrees. If you think PyMOL is incorrectly assigning secondary structure, then you might try varying these. 

Retrieved from "[https://pymolwiki.org/index.php?title=Dss&oldid=9168](https://pymolwiki.org/index.php?title=Dss&oldid=9168)"


---

## Get area

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_area** calculates the surface area in square Angstroms of the selection given. Note that the accessibility is assessed in the context of the object(s) that the selection is part of. So, to get the surface areas of e.g. a component of a complex, you should make a new object containing a copy of just that component and calculate its area. 

The "get_area selection" command will return the effective surface area of the dots that you would see from "[show](/index.php/Show "Show") [dots](/index.php/Dots "Dots"), selection". This is a discrete approximation -- not an exact calculation. 

**Attention:** Atoms with the "[ignore](/index.php/Surface#Calculating_a_partial_surface "Surface")" [flag](/index.php/Flag "Flag") may lead to unexpected results. Clear the "ignore" flag first or exclude those atoms from the calculation. _Fixed in PyMOL 2.5,[ignored atoms are now excluded](https://github.com/schrodinger/pymol-open-source/commit/1659fde83d2a125f86cad13414cead6b8f74abe9). Previously, their entire sphere surface was added to the surface area._

## Contents

  * 1 Usage
    * 1.1 Arguments
    * 1.2 Settings
    * 1.3 Python API
  * 2 Examples
    * 2.1 Example 1 - starting with a complex in a single file
    * 2.2 Example 2 - starting with two components in separate files
    * 2.3 Example 3 - using load_b to get surface area per atom
    * 2.4 Example 4 - using a Python script to compute the SASA for individual residues
  * 3 See Also



# Usage
    
    
    get_area [ selection [, state [, load_b ]]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **state** = int: object state {default: 1}
  * **load_b** = 0/1: Load the surface area per atom into the b-factor {default: 0}



## Settings

The following PyMOL settings control how **get_area** works: 

Setting  | Description   
---|---  
[dot_solvent](/index.php?title=Dot_solvent&action=edit&redlink=1 "Dot solvent \(page does not exist\)") | **0** = calculate **molecular surface** area {default} **1** = calculate **solvent accessible surface** area (SASA)   
[dot_density](/index.php?title=Dot_density&action=edit&redlink=1 "Dot density \(page does not exist\)") | 1-4 {default: 2}. Sampling density. Higher density (more dots) means higher accuracy but slower performance.   
[solvent_radius](/index.php/Solvent_radius "Solvent radius") | The solvent radius {default: 1.4}   
  
The [dot_hydrogens](/index.php?title=Dot_hydrogens&action=edit&redlink=1 "Dot hydrogens \(page does not exist\)") setting is **not** used. 

For example: 
    
    
    PyMOL> load $TUT/1hpv.pdb
    PyMOL> show dots, resn arg
    PyMOL> get_area resn arg
    cmd.get_area: 1147.956 Angstroms^2.
    PyMOL>set dot_solvent, on
    PyMOL>get_area resn arg
     cmd.get_area: 673.084 Angstroms^2.
    PyMOL>set dot_density, 3
    PyMOL>get_area resn arg
     cmd.get_area: 674.157 Angstroms^2.
    PyMOL>set dot_density, 4
    PyMOL>get_area resn arg
     cmd.get_area: 672.056 Angstroms^2.
    PyMOL>get_area all
     cmd.get_area: 13837.804 Angstroms^2.
    

This code has not been recently validated though it was checked a couple years back. We suggest that people perform some kind of independent check on their system before trusting the results. 

## Python API
    
    
    float cmd.get_area(str selection="(all)", int state=1, int load_b=0)
    

# Examples

## Example 1 - starting with a complex in a single file
    
    
    # load complex
    # Haemoglobin in this example illustrates careful use of selection algebra
    load 2HHB.pdb
    
    # create objects for alpha1, beta1 and alpha1,beta1 pair of subunits
    create alpha1, 2HHB and chain A
    create beta1, 2HHB and chain B
    create ab1, 2HHB and chain A+B
    
    # get hydrogens onto everything (NOTE: must have valid valences on e.g. small organic molecules)
    h_add
    
    # make sure all atoms including HETATM within an object occlude one another, but ignore solvent
    flag ignore, none
    flag ignore, solvent
    
    # use solvent-accessible surface with high sampling density
    set dot_solvent, 1
    set dot_density, 3
    
    # measure the components individually storing the results for later
    alpha1_area=cmd.get_area("alpha1")
    beta1_area=cmd.get_area("beta1")
    
    
    # measure the alpha1,beta1 pair
    ab1_area=cmd.get_area("ab1")
    
    # now print results and do some maths to get the buried surface
    print alpha1_area
    print beta1_area
    print ab1_area
    print (alpha1_area + beta1_area) - ab1_area
    

## Example 2 - starting with two components in separate files
    
    
    # load components separately
    load my_ligand.pdb
    load my_target.pdb
    
    # get hydrogens onto everything (NOTE: must have valid valences on the ligand...)
    h_add
    
    # make sure all atoms including HETATM within an object occlude one another but ignore solvent
    flag ignore, none
    flag ignore, solvent
    
    # use solvent-accessible surface with high sampling density
    set dot_solvent, 1
    set dot_density, 3
    
    # measure the components individually
    ligand_area=cmd.get_area("my_ligand")
    target_area=cmd.get_area("my_target")
    
    # create the complex
    create my_complex, my_ligand my_target
    
    # measure the complex
    complex_area=cmd.get_area("my_complex")
    
    # now print results
    print ligand_area
    print target_area
    print complex_area
    print (ligand_area + target_area) - complex_area
    

## Example 3 - using load_b to get surface area per atom
    
    
    # example usage of load_b
    # select some organic small molecule
    select ligand, br. first organic
    # get its area and load it into it's b-factor column
    get_area ligand, load_b=1
    # print out the b-factor/areas per atom
    iterate ligand, print b
    

  


## Example 4 - using a Python script to compute the SASA for individual residues
    
    
    import __main__
    __main__.pymol_argv = ['pymol','-qc']
    import pymol
    from pymol import cmd, stored
    pymol.finish_launching()
    
    cmd.set('dot_solvent', 1)
    cmd.set('dot_density', 3)
    
    cmd.load('file.pdb')  # use the name of your pdb file
    stored.residues = []
    cmd.iterate('name ca', 'stored.residues.append(resi)')
    
    sasa_per_residue = []
    for i in stored.residues:
        sasa_per_residue.append(cmd.get_area('resi %s' % i))
    
    print sum(sasa_per_residue)
    print cmd.get_area('all')  # just to check that the sum of sasa per residue equals the total area
    

  


# See Also

  * For an example of **load_b** in use check out [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues").
  * [Surface](/index.php/Surface "Surface"), most notably [Surface#Calculating_a_partial_surface](/index.php/Surface#Calculating_a_partial_surface "Surface").



Retrieved from "[https://pymolwiki.org/index.php?title=Get_area&oldid=13461](https://pymolwiki.org/index.php?title=Get_area&oldid=13461)"


---

## Modeling and Editing Structures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Modeling in PyMOL
    * 1.1 Saving with transformed coordinates
    * 1.2 Translate or rotate individual objects
    * 1.3 Moving one segment relative to the rest
    * 1.4 Split states to objects
    * 1.5 Altering secondary structures
    * 1.6 Altering van der Waals radii
    * 1.7 Altering atom coordinates
    * 1.8 Deleting bonds
    * 1.9 Converting D- to L- amino acids
    * 1.10 Adding disulfide bonds
    * 1.11 Adding hydrogen bonds
    * 1.12 Protonating ligands
    * 1.13 Superposition of two molecules
    * 1.14 Manual superposition of two molecules
    * 1.15 Adding and using your own fragments



# Modeling in PyMOL

## Saving with transformed coordinates

Here is a simple script that saves the molecule with coordinates from the current orientation. (invoke it with 'run save_transformed.py' and type the new save_transformed.py command thereafter). 
    
    
    # Adds the command save_transformed
    # Usage: save_transformed object, file
    def save_transformed(object,file):
        m = cmd.get_view(0)
        ttt = [m[0], m[1], m[2], 0.0,
               m[3], m[4], m[5], 0.0,
               m[6], m[7], m[8], 0.0,
               0.0,   0.0,  0.0, 1.0]
        cmd.transform_object(object,ttt,transpose=1)
        cmd.save(file,object)
    
    cmd.extend('save_transformed',save_transformed)
    

  


## Translate or rotate individual objects

There is a "[translate](/index.php/Translate "Translate")" function similar to "[rotate](/index.php/Rotate "Rotate")", feel free to use them in the following forms: 
    
    
    translate vector,object-name,state
       (vector needs to be something like [x,y,z])
    
       translate [1,0,0],pept
    
    rotate axis,angle,object-name,state
       (axis can be either the letter x,y,z or a 3D vector [x,y,z])
       rotate x,90,pept
       rotate [1,1,1],10,pept
    

## Moving one segment relative to the rest

This means moving two parts of one object into different directions. The easiest way to do this is to split the objects and then use the rotate command. 
    
    
    load 1FJ1.pdb
    
    # split PDB file
    
    create anti=(chain F) 
    create fab=(chain A,B)
    
    # delete original object
    delete 1FJ1
    
    # color objects
    color green,fab
    color pink,anti
    
    # color interface
    select inter = (byres ((fab within 5 of anti)\
       or (anti within 5 of fab)))
    
    color yellow,inter
    
    # splay apart
    orient
    origin fab
    rotate y,60,fab
    origin anti
    rotate y,-60, anti
    
    # zoom interface region
    zoom inter
    show sph,inter
    disable inter
    

## [Split states](/index.php/Split_States "Split States") to objects

There is also a new command in the 0.95 series: 
    
    
       split_states object-name
    

which will spread a PDB "biological unit" (or any multi-state object -- including SD files) over a series of independent objects. This makes it possible to interact with such objects more naturally than with "all_states = 1". 

  


## Altering secondary structures

Examples: 
    
    
     alter A/10:34/, ss='H'
     alter A/35:40/, ss='L'
     alter A/41:60/, ss='S'
    

## Altering van der Waals radii

Example: 
    
    
    alter (elem Fe),vdw=1.8
    rebuild
    

(The value for Fe is wrecked in PyMOL at the moment, so running the above line might be a good idea). 

  


## Altering atom coordinates

Example: 
    
    
    alter_state 1,(pdb1cse),x=x-10.0
    

The latter section can contain formulae involving at least the xyz coordinates, lots of constants and the (+-*/) operators. 

  


## Deleting bonds

Select the bond using Ctrl-right-click, then either 
    
    
    unbond pk1,pk2
    

or hit Ctrl-D. 

  


## Converting D- to L- amino acids

The inversion function was changed in version 0.95 to take advantage of multiple picked atoms. To invert a center, Ctrl-middle-click to pick the center atom as pk1 and two stationary atoms as pk2 and pk3. Then type Ctrl-E to invert. 

  


## Adding disulfide bonds

You can use the [bond](/index.php/Bond "Bond") command to attach them: 
    
    
    bond 24/sg,26/sg
    bond 56/sg,99/sg
    unpick
    

(unpick will hide the bond baton which gets displayed.) Additionally, the residue names can be changed for bonded cysteines: 
    
    
    alter cys/,name='CYX'
    

or for specific residues 
    
    
    alter 24+26+56+99/,name='CYX'
    

  


## Adding hydrogen bonds

    _See[Displaying biochemical properties](/index.php/Displaying_Biochemical_Properties#Hydrogen_bonds_and_Polar_Contacts "Displaying Biochemical Properties")._

## Protonating ligands

If your ligands come in with valid valencies and formal charges, PyMOL's h_add command can protonate ligands. (NOTE that there is a minor technical hiccup with SD-files which are loaded by default as immutable "discrete" objects.) Suffice it to say that in order to make changes to the chemical structure, an object must be loaded with the "discrete" flag set to zero. Unfortunately, much of the molecular editing stuff remains to be documented. Here's an example sequence, but I'm not sure it will help to much...as indicated in the manual, this is immature functionality with some major gaps. Attach in particular is very limited... 
    
    
    # show valences
    set valence=0.05
    
    # load cysteine fragment
    fragment cys
    
    # remove hydrogens
    remove (hydro)
    
    # edit gamma S
    edit cys////sg
    
    # add hydrogen
    attach H,1,1
    
    # add planer, trivalent nitrogen onto C terminus
    edit cys////C
    attach N,3,3
    
    # edit that nitrogen
    edit (elem N and neighbor cys////C)
    
    # attach a tetrahedral methyl (note random position)
    attach C,4,4
    
    # here's an example of adding a whole residue from the library
    edit cys////N
    editor.attach_amino_acid("pk1","ace")
    
    # now restore missing hydrogens (note that the names are off...)
    h_add
    

  


## Superposition of two molecules

Using pair_fit requires that you specify a set of paired atoms in each structure. Fortunately, you no longer have to specify each pair separately, so long as the ordering is the same in each selection (almost always true). 
    
    
    pair_fit ( trna10 and resid 10:15 and name P ), ( ref4 and resid 10:15 and name P )
    

Another example: 
    
    
    pair_fit prot1///11-26/CA, prot2///34-49/CA
    

would superimpose prot1 on prot2 using C-alphas from residues 11-26 in prot1 and 34-49 in prot2. 

## Manual superposition of two molecules

You can also align to structures using mouse rotation/translation. For this, you need to protect those molecules you don't want to move with (action menu -> movement -> protect) in the selection menu. 

Protect one object, deprotect the other, grab the deprotected object and move with Shift-Mouse. Don't forget to switch to Mouse Editing mode. 

  


## Adding and using your own fragments

Pymol has some build-in fragments (amino acids and simple functional groups). You can add your own fragments, eg. sugars, in this way: 

Create the molecule you want to use as a fragment. Save it as a .pkl file in <pymol_path>/data/chempy/fragments. 

How to use the fragment: 

Pick the atom (ctrl-middle) where you want to add the fragment. This will usually be a hydrogen atom (which will be removed). Then use the command: 
    
    
    editor.attach_fragment('pk1','my_fragment_name',11,0)
    

where my_fragment_name is the name of the pkl-file (w/o .pkl extension) and 11 is the number of the connecting (hydrogen) atom in the fragment. To determine this number, press '[L]abel' -> 'atom identifiers' -> 'index' and choose the hydrogen atom you want. 

If you want a menu item for your fragment, you can probably put it in <pymol_path>/modules/pmg_tk/skins/normal/__init__.py, but I haven't tried this. 

Retrieved from "[https://pymolwiki.org/index.php?title=Modeling_and_Editing_Structures&oldid=12249](https://pymolwiki.org/index.php?title=Modeling_and_Editing_Structures&oldid=12249)"


---

## Mole

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

#### There is a new version of MOLE [check it out](/index.php/MOLE_2.0:_advanced_approach_for_analysis_of_biomacromolecular_channels "MOLE 2.0: advanced approach for analysis of biomacromolecular channels")!

[MOLE](http://mole.chemi.muni.cz) is an universal toolkit for rapid and fully automated location and characterization of channels, tunnels and pores in molecular structures. The core of MOLE algorithm is a Dijsktra path search algorithm, which is applied to a Voronoi mesh. MOLE is a powerful software (overcomming some limitations of [Caver](/index.php/Caver "Caver") tool) for exploring large molecular channels, complex networks of channels and molecular dynamics trajectories (AMBER ascii traj and parm7 are supported) in which analysis of a large number of snapshots is required. 

Program is available as standalone application for Linux(Unix), MacOS and Windows. Plugin for pyMol enables calling the MOLE program via an user friendly interface and allows results visualization. The on-line interface provides another interactive and easy-to-use interface to analyze molecular channels. 

[Tutorial](http://mole.chemi.muni.cz/web/docs/doc_pymolplugin.php) is also available. 

References: 

Petrek M., Kosinova P., Koca J., Otyepka M.: [MOLE: A Voronoi Diagram-Based Explorer of Molecular Channels, Pores, and Tunnels. Structure 2007, 15, 1357-1363](http://dx.doi.org/10.1016/j.str.2007.10.007)

[![KscA.png](/images/8/81/KscA.png)](/index.php/File:KscA.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Mole&oldid=11428](https://pymolwiki.org/index.php?title=Mole&oldid=11428)"


---

## MOLE 2.0: advanced approach for analysis of biomacromolecular channels

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![Molelogo.png](/images/d/d6/Molelogo.png)](/index.php/File:Molelogo.png)

[![](/images/a/ad/1jj2.png)](/index.php/File:1jj2.png)

[](/index.php/File:1jj2.png "Enlarge")

Example of channel system of 1JJ2

[MOLE 2](http://mole.chemi.muni.cz) is a successor of a popular software tool [Mole](http://lcc.ncbr.muni.cz/whitezone/development/mole/web/index.php) for detection and characterization of channels in biomacromolecules (proteins, nucleic acids and glycans). This completely redesigned version allows user rapid and accurate **analysis of channels and transmembrane pores** even **in large structures** (hundreds of thousands of atoms). As a new feature, MOLE 2 **estimates** not only **physicochemical properties** of the identified channels, i.e., **hydropathy, hydrophobicity, polarity, charge, and mutability** , but also physicochemical properties of the cavities. For thorough description of the functionality and instructions on working with MOLE 2 please consult the [paper](http://dx.doi.org/10.1186/1758-2946-5-39) or our [manual](http://webchem.ncbr.muni.cz/Platform/AppsBin/Mole2_manual.pdf). 

## MOLE 2.5 Update

Update of MOLE 2 with novel functionality and bug fixes has been released. The binaries and plugins are available for download from the [| main page](https://webchem.ncbr.muni.cz/Platform/App/Mole). 

  * Easily remove parts of the PDB entry with [PatternQuery](https://webchem.ncbr.muni.cz/Wiki/PatternQuery:UserManual).
  * New PDB standard mmCIF is supported and recommended for all calculations.
  * New weight functions for tunnel/pore calculation.



  


## Availability and Requirements

MOLE 2.5 is available as a [GUI application](https://webchem.ncbr.muni.cz/Platform/App/Mole) with in-built molecular browser enabling user interactive work; [command-line application](https://webchem.ncbr.muni.cz/Platform/App/Mole) and [PyMOL & Chimera plugin](https://webchem.ncbr.muni.cz/Platform/App/Mole). Some functionality is also available via online [web service](http://beta.mole.upol.cz). 

  * GUI application is available for Windows with [.NET 4.0 framework](http://www.microsoft.com/en-us/download/details.aspx?id=17851) or newer installed.
  * Command-line application and PyMOL&Chimera plugin are available for Windows, Linux (Unix) and MacOS. For Linux (Unix) and MacOS systems [Mono framework](http://mono-project.com/Main_Page) 2.10 or newer is required.



  


[![](/images/b/b7/Plugin.png)](/index.php/File:Plugin.png)

[](/index.php/File:Plugin.png "Enlarge")

MOLE 2.5 PyMOL plugin

## References

  1. Sehnal D, Svobodová Vařeková R, Berka K, Pravda L, Navrátilová V, Banáš P, Ionescu C-M, Otyepka M, Koča J: [MOLE 2.0: advanced approach for analysis of biomacromolecular channels](http://dx.doi.org/10.1186/1758-2946-5-39). Journal of Cheminformatics 2013, 5:39.
  2. Pravda,L., Berka,K., Svobodová Vařeková,R., Sehnal,D., Banáš,P., Laskowski,R.A., Koča,J. and Otyepka,M. (2014) [Anatomy of enzyme channels](http://dx.doi.org/10.1186/s12859-014-0379-x). BMC Bioinformatics, 15, 379.
  3. Berka K, Hanák O, Sehnal D, Banáš P, Navrátilová V, Jaiswal D, Ionescu C-M, Svobodová Vařeková R, Koča J, Otyepka M: [MOLEonline 2.0: interactive web-based analysis of biomacromolecular channels](http://dx.doi.org/10.1093/nar/GKS363). Nucleic acids research 2012, 40:W222–7.



Retrieved from "[https://pymolwiki.org/index.php?title=MOLE_2.0:_advanced_approach_for_analysis_of_biomacromolecular_channels&oldid=12616](https://pymolwiki.org/index.php?title=MOLE_2.0:_advanced_approach_for_analysis_of_biomacromolecular_channels&oldid=12616)"


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

## Nuccyl

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

See [Nuccyl Home](http://www.biosci.ki.se/groups/ljo/software/nuccyl.html)

## Overview

**nuccyl** is a [Perl](http://www.perl.org/) program that allows [PyMOL](http://www.pymol.org/), a powerful open-source molecular modeling system, to display atomic models of nucleic acids in a highly simplified representation. By depicting the nucleic acid backbone as a cartoon and the bases of nucleotides as cylinders, similarly to the excellent programs [Ribbons](http://sgce.cbse.uab.edu/ribbons/) and [Drawna](http://www-ibmc.u-strasbg.fr/upr9002/westhof/download.html), nuccyl and PyMOL allow to quickly grasp the overall folding of an RNA or DNA molecule 

[![](/images/6/63/Phe_trna_small.jpg)](/index.php/File:Phe_trna_small.jpg)

[](/index.php/File:Phe_trna_small.jpg "Enlarge")

Example Nuccyl Image

[![](/images/5/5a/45s_rna_small.jpg)](/index.php/File:45s_rna_small.jpg)

[](/index.php/File:45s_rna_small.jpg "Enlarge")

Example Nuccyl II

[![](/images/b/be/50s_small.jpg)](/index.php/File:50s_small.jpg)

[](/index.php/File:50s_small.jpg "Enlarge")

Example Nuccyl III

[![](/images/7/7a/P4p6_small.jpg)](/index.php/File:P4p6_small.jpg)

[](/index.php/File:P4p6_small.jpg "Enlarge")

Example Nuccyl IV

nuccyl produces base cylinder coordinates and command files to display them with PyMOL. The base pair analysis programs [RNAView](http://beta-ndb.rutgers.edu/services/download/index.html#rnaview) and [3DNA](http://rutchem.rutgers.edu/%7exiangjun/3DNA/) can also be optionally used to obtain starting input files for nuccyl. 

[RNAView](http://beta-ndb.rutgers.edu/services/download/index.html#rnaview) can be downloaded and installed as [described](http://ndbserver.rutgers.edu/services/help/rnaview-readme.html); a [publication](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12824344&dopt=Abstract) describing the program is also available. Information and download links for [3DNA](http://rutchem.rutgers.edu/%7exiangjun/3DNA/) may be found at the program's [home page](http://rutchem.rutgers.edu/%7exiangjun/3DNA/); a [paper](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12930962&dopt=Abstract) has also been published. 

* * *

The four example images Copyright © 2000-2005 Luca Jovine. Do **NOT** reuse them without contacting him or reading his release license for the images! 

Retrieved from "[https://pymolwiki.org/index.php?title=Nuccyl&oldid=5397](https://pymolwiki.org/index.php?title=Nuccyl&oldid=5397)"


---

## Protein contact potential

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  * Comparing PyMOL-generated electrostatics and APBS' version.
  * [![PyMOL Charge-smoothed potential \(see below\)](/images/6/64/Prot_contact_pot.png)](/index.php/File:Prot_contact_pot.png "PyMOL Charge-smoothed potential \(see below\)")

PyMOL Charge-smoothed potential (see below) 

  * [![APBS-generated potential \(see APBS\).](/images/6/65/Prot_contact_APBS.png)](/index.php/File:Prot_contact_APBS.png "APBS-generated potential \(see APBS\).")

APBS-generated potential (see [APBS](/index.php/APBS "APBS")). 




Protein contact potential is an automated PyMOL representation where the false red/blue charge-smoothed surface is shown on the protein. 

The rule of thumb with respect to PyMOL's internal "protein contact potential" is that if you care enough to be concerned with how it works, then you should instead be using a true Possion-Boltzman electrostatics solver such as [APBS](/index.php/APBS "APBS"). 

Regardless, what PyMOL does to generate a qualitative electrostatic representation (via action **popup- >generate->vacuum electrostatics->...**) amounts to averaging charges over a small region of space using a quasi-Coulombic-shaped convolution function. Nathan Baker has suggested that we simply refer to this approach as "charge smoothing". 

PyMOL's use of the term "contact" relates to the fact that the potential shown by the default coloring approximates the potential that would be felt by a point-charge one solvent radius *above* the protein surface, if we ignore solvent screening and only consider nearby atoms. My personal opinion is that this sort of treatment makes a lot sense when viewing results from Poisson-Boltzmann calculations as well. (Select `Color by potential on sol. acc. surf.` in the visualization tab of the APBS plugin to get that effect.) 

As an aside: It has never made any sense to me whatsoever that electrostatics visualization programs would show potentials located right on the molecular surface given that (1) we employ point charge models that have *only* been parameterized to approximate potential energies computed between partial charges located at atomic centers and that (2) the molecular surface is located in the region of space where artifacts & noise in a PB calculation would be greatest due to discrete discontinuities between the low-dielectric interior and the high-dielectric exterior (solvent region). 

Retrieved from "[https://pymolwiki.org/index.php?title=Protein_contact_potential&oldid=6908](https://pymolwiki.org/index.php?title=Protein_contact_potential&oldid=6908)"


---

## Pucker

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [pucker.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/pucker.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 DESCRIPTION
  * 2 USAGE
  * 3 EXAMPLES
  * 4 PYMOL API



### DESCRIPTION

**pucker.py** is a PyMol script that returns the sugar pucker information (phase, amplitude, pucker) for a given selection. 

This script uses its own dihedral calculation scheme rather than the get_dihedral command. Thus, it is lightning fast! 

If a selection does not contain any ribose sugars then an error message is returned. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    pucker selection
    #Calculates the sugar information for the given selection
    
    pucker selection, state=1
    #Calculates the sugar information for the given selection in state 1
    
    pucker selection, state=0 
    #Iterates over all states and calculates the sugar information for the given selection
    

### EXAMPLES
    
    
    #fetch 1BNA.pdb
    fetch 1bna
    
    #Select DNA only
    #Otherwise, you will get an error for water not having sugars
    select DNA, not solvent
    
    #Execute pucker command
    pucker DNA
    
    #The output should look like this
     Phase     Amp    Pucker  Residue
    161.22   55.32  C2'-endo  A 1
    139.52   41.67   C1'-exo  A 2
     92.82   38.31  O4'-endo  A 3
    166.35   48.47  C2'-endo  A 4
    128.57   46.30   C1'-exo  A 5
    126.92   49.75   C1'-exo  A 6
    101.30   47.32  O4'-endo  A 7
    115.62   49.23   C1'-exo  A 8
    140.44   46.37   C1'-exo  A 9
    145.79   53.36  C2'-endo  A 10
    147.47   47.04  C2'-endo  A 11
    113.80   51.69   C1'-exo  A 12
     
     Phase     Amp    Pucker  Residue
    153.24   43.15  C2'-endo  B 13
    128.49   45.01   C1'-exo  B 14
     67.74   43.84   C4'-exo  B 15
    149.33   41.01  C2'-endo  B 16
    169.27   42.30  C2'-endo  B 17
    147.03   42.30  C2'-endo  B 18
    116.29   47.52   C1'-exo  B 19
    129.62   49.92   C1'-exo  B 20
    113.61   42.93   C1'-exo  B 21
    156.34   50.98  C2'-endo  B 22
    116.89   44.21   C1'-exo  B 23
     34.70   45.55  C3'-endo  B 24
    

  


### PYMOL API
    
    
    from pymol.cgo import *    # get constants
    from math import *
    from pymol import cmd
    
    def pucker(selection,state=1):
    
    	# Author: Sean Law
    	# Institute: University of Michigan
    	# E-mail: seanlaw@umich.edu
    
    	#Comparison to output from 3DNA is identical
    	#Note that the 3DNA output has the second chain
    	#printed out in reverse order
    	state=int(state)
    	if state == 0:
    		for state in range(1,cmd.count_states()+1):
    			sel=cmd.get_model(selection,state)
    			if state == 1:
    				print " " #Blank line to separate chain output
    				print "%5s  %6s  %6s  %8s  Residue" % ("State","Phase","Amp", "Pucker")
    			get_pucker(sel,iterate=state)
    	else:
    		sel=cmd.get_model(selection,state)
    		get_pucker(sel,iterate=0)
    	return
    
    def get_pucker (sel,iterate=0):
    	iterate=int(iterate)
    	first=1
    	old=" "
    	oldchain=" "
    	residue={}
    	theta={}
    	count=0
    	for atom in sel.atom:
    		new=atom.chain+" "+str(atom.resi)
    		newchain=atom.chain+" "+atom.segi
    		if oldchain != newchain and first:
    			if iterate == 0:
    				print " " #Blank line to separate chain output
    				print "%6s  %6s  %8s  Residue" % ("Phase", "Amp", "Pucker")
    		if new != old and not first:
    			#Check that all 5 atoms exist
    			if len(residue) == 15:
    				#Construct planes
    				get_dihedrals(residue,theta)
    				#Calculate pucker
    				info = pseudo(residue,theta)
    				print info+"  "+old
    			else:
    				print "There is no sugar in this residue"
    			if oldchain != newchain:
    				print " " #Blank line to separate chain output
    				print "%6s  %6s  %8s  Residue" % ("Phase", "Amp", "Pucker")
    			#Clear values
    			residue={}
    			dihedrals={}
    			theta={}
    			#Store new value
    			store_atom(atom,residue)
    		else:
    			store_atom(atom,residue)
    		first=0
    		old=new
    		oldchain=newchain
    	#Final Residue
    	#Calculate dihedrals for final residue
    	if len(residue) == 15:
    		#Construct planes
    		get_dihedrals(residue,theta)
    		#Calculate pucker for final residue
    		info = pseudo(residue,theta)
    		if iterate == 0:
    			print info+"  "+old
    		else:
    			print "%5d  %s" % (iterate,info+"  "+old)
    	else:
    		print "There is no sugar in this residue"
    	return
    
    def sele_exists(sele):
    	return sele in cmd.get_names("selections");
    
    def pseudo(residue,theta):
    	other=2*(sin(math.radians(36.0))+sin(math.radians(72.0)))
    	
    	#phase=atan2((theta4+theta1)-(theta3+theta0),2*theta2*(sin(math.radians(36.0))+sin(math.radians(72.0))))
    	phase=atan2((theta['4']+theta['1'])-(theta['3']+theta['0']),theta['2']*other)
    	amplitude=theta['2']/cos(phase)
    	phase=math.degrees(phase)
    	if phase < 0:
    		phase+=360 # 0 <= Phase < 360
    	#Determine pucker
    	if phase < 36:
    		pucker = "C3'-endo"
    	elif phase < 72:
    		pucker = "C4'-exo"
    	elif phase < 108:
    		pucker = "O4'-endo"
    	elif phase < 144:
    		pucker = "C1'-exo"
    	elif phase < 180:
    		pucker = "C2'-endo"
    	elif phase < 216:
    		pucker = "C3'-exo"
    	elif phase < 252:
    		pucker = "C4'-endo"
    	elif phase < 288:
    		pucker = "O4'-exo"
    	elif phase < 324:
    		pucker = "C1'-endo"
    	elif phase < 360:
    		pucker = "C2'-exo"
    	else:
    		pucker = "Phase is out of range"
    	info = "%6.2f  %6.2f  %8s" % (phase, amplitude, pucker)
    	return info
    	
    
    def store_atom(atom,residue):
    	if atom.name == "O4'" or atom.name == "O4*":
    		residue['O4*X'] = atom.coord[0]
    		residue['O4*Y'] = atom.coord[1]
    		residue['O4*Z'] = atom.coord[2]
    	elif atom.name == "C1'" or atom.name == "C1*":
    		residue['C1*X'] = atom.coord[0]
    		residue['C1*Y'] = atom.coord[1]
    		residue['C1*Z'] = atom.coord[2]
    	elif atom.name == "C2'" or atom.name == "C2*":
    		residue['C2*X'] = atom.coord[0]
    		residue['C2*Y'] = atom.coord[1]
    		residue['C2*Z'] = atom.coord[2]
    	elif atom.name == "C3'" or atom.name == "C3*":
    		residue['C3*X'] = atom.coord[0]
    		residue['C3*Y'] = atom.coord[1]
    		residue['C3*Z'] = atom.coord[2]
    	elif atom.name == "C4'" or atom.name == "C4*":
    		residue['C4*X'] = atom.coord[0]
    		residue['C4*Y'] = atom.coord[1]
    		residue['C4*Z'] = atom.coord[2]
    	return
    
    def get_dihedrals(residue,theta):
    
    	C = []
    	ribose = residue.keys()
    	ribose.sort()
    	
    	shift_up(ribose,6)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['0']=calcdihedral(C)
    	
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['1']=calcdihedral(C)
    
    
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['2']=calcdihedral(C)
    
    
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['3']=calcdihedral(C)
    
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['4']=calcdihedral(C)
    	
    	return
    
    def shift_up(list,value):
    	for i in range (0,value):
    		list.insert(0,list.pop())
    	return
    
    def shift_down(list,value):
    	for i in range (0,value):
    		list.append(list.pop(0))
    	return
    
    def calcdihedral(C):
    	
    	DX12=C[0]-C[3]
    	DY12=C[1]-C[4]
    	DZ12=C[2]-C[5]
    	
    	DX23=C[3]-C[6]
    	DY23=C[4]-C[7]
    	DZ23=C[5]-C[8]
    	
    	DX43=C[9]-C[6];
    	DY43=C[10]-C[7];
    	DZ43=C[11]-C[8];
    
    	#Cross product to get normal
    	
    	PX1=DY12*DZ23-DY23*DZ12;
    	PY1=DZ12*DX23-DZ23*DX12;
    	PZ1=DX12*DY23-DX23*DY12;
    
    	NP1=sqrt(PX1*PX1+PY1*PY1+PZ1*PZ1);
    
    	PX1=PX1/NP1
    	PY1=PY1/NP1
    	PZ1=PZ1/NP1
    
    	PX2=DY43*DZ23-DY23*DZ43;
    	PY2=DZ43*DX23-DZ23*DX43;
    	PZ2=DX43*DY23-DX23*DY43;
    
    	NP2=sqrt(PX2*PX2+PY2*PY2+PZ2*PZ2);
    	
    	PX2=PX2/NP2
    	PY2=PY2/NP2
    	PZ2=PZ2/NP2
    
    	DP12=PX1*PX2+PY1*PY2+PZ1*PZ2
    
    	TS=1.0-DP12*DP12
    
    	if TS < 0:
    		TS=0
    	else:
    		TS=sqrt(TS)
    	
    	ANGLE=math.pi/2.0-atan2(DP12,TS)
    
    	PX3=PY1*PZ2-PY2*PZ1
    	PY3=PZ1*PX2-PZ2*PX1
    	PZ3=PX1*PY2-PX2*PY1
    
    	DP233=PX3*DX23+PY3*DY23+PZ3*DZ23
    
    	if DP233 > 0:
    		ANGLE=-ANGLE
    
    	ANGLE=math.degrees(ANGLE)
    
    	if ANGLE > 180:
    		ANGLE-=360
    	if ANGLE < -180:
    		ANGLE+=360
    
    	return ANGLE
    
    cmd.extend("pucker",pucker)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pucker&oldid=11145](https://pymolwiki.org/index.php?title=Pucker&oldid=11145)"


---

## Pytms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/pytms.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/pytms.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [GPL-2.0](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![PyTMs example: PTMs in red](/images/d/dd/Pytms.gif)](/index.php/File:Pytms.gif "PyTMs example: PTMs in red")

**PyTMs** is a PyMOL plugin that enables to easily introduce a set of common post-translational modifications (PTMs) into protein models.  
Latest version: **1.2 (November 2015)**  


## Contents

  * 1 Installation
  * 2 Using the GUI
  * 3 PyMOL API functions
  * 4 Selections and custom property selectors
  * 5 PyTMs in python scripts
  * 6 Basic residue-based optimization, animation, vdW clashes, surface selections
  * 7 Hints
  * 8 Update notes
  * 9 References
  * 10 SEE ALSO



## Installation

  * For instruction on setting up plugins see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



PyTMS can be [downloaded separately](https://github.com/andwar/Pymol-script-repo/blob/master/plugins/pytms.py), or together with the [pymol-script-repository](https://github.com/Pymol-Scripts/Pymol-script-repo)  
The **pytms.py** file should be placed in the **plugins** folder. 

## Using the GUI

In PyMOL, navigate to: _P_ lugin --> PyTMs  
[![PyTMs GUI](/images/6/6a/Pytms_menu.png)](/index.php/File:Pytms_menu.png "PyTMs GUI")  
The supported PTM can be selected from the left panel. The options will appear on the right. For help, click the help button.  
Note that available objects/selection can be selected using the buttons or be entered into the field. Target atoms are automatically sub-selected and do not need to be specified. 

  


## PyMOL API functions

List of functions   
---  
PTM  | PyTMs API command  | p.PTM property selector *  | Target amino acid(s)  | Modifyable N-terminus? **   
acetylation  | acetylate  | acetylation  | Lysine  | Yes, biologically relvant   
carbamylation  | carbamylate  | carbamylation  | Lysine  | Yes, implemented   
citrullination  | citrullinate  | citrullination  | Arginine  | No   
cysteine oxidations  | oxidize_cys  | oxidation_cys  | Cysteine /(Seleno-)  | No   
malondialdehyde adducts  | mda_modify  | MDA  | Lysine  | Yes, biologically relevant   
methionine oxidation  | oxidize_met  | oxidation_met  | Methionine /(Seleno-)  | No   
methylation  | methylate  | methylation  | Lysine  | Yes, implemented   
nitration  | nitrate  | nitration  | Tyrosine (Tryptophan)  | No   
nitrosylation  | nitrosylate  | nitrosylation  | Cysteine  | No   
proline hydroxylation  | hydroxy_pro  | hydroxylation_pro  | Proline  | No   
phosphorylation  | phosphorylate  | phosphorylation  | Serine, Threonine, Tyrosine  | No   
  
  * (*) requires incentive PyMOL version 1.7, or higher
  * (**) N-terminal modifiactions excluding Proline
  * Note that each function has an **individual help description** , callable by entering, e.g.:


    
    
    help citrullinate
    

  


## Selections and custom property selectors

PTM amino acids have altered names, and altered/added atoms are identifyable by unique names.  
The nomenclature is adapted mostly from [RCSB](http://www.rcsb.org/pdb/home/home.do). 
    
    
    Examples
    select resn CIR # citrullines
    select resn NIY # nitro-tyrosines
    # etc. ...
    

  * Note that more information can be found in the associated help of individual functions
  * Hint: using the [label](/index.php/Label "Label") command is a convenient way to get hold of the names of atoms or residues



  
Incentive PyMOL version (>1.7) support custom property selectors.  
PyTMs supports these by intoducting the property **p.PTM** selector (cf. table above), which will select the altered PTM atoms:  

    
    
    # To select a PTM, use e.g.:
    select p.PTM in nitration #selects nitro groups
    select p.PTM in * #selects any PTM
    

  


## PyTMs in python scripts

In simple API scripts, the API commands can be used directly.  
For use in python scripts PyTMs can also be imported, e.g.: 
    
    
    import pmg_tk.startup.pytms
    pmg_tk.startup.pytms.mda_modify()
    

  


## Basic residue-based optimization, animation, vdW clashes, surface selections

Some functions (MDA & phosphorylation) support basic structure optimization that is used to position PTMs in more favorable locations. Note that this optimization is based solely on sterical vdW strain and currently ignores charge. Due to interation of testing, there may be a significant calculation time associated with this procedure. The optimization can be followed in the console and graphically in the PyMOL window (only in case the GUI is used). There is an option to animate the states from original to final after calculation for reference (cf. example).  
Note that the presence of hydrogens will significantly affect both the calculation during optimization and the display of vdW clashes!  
  
All functions support the display of vdW clashes after modification to allow assessment of potential steric interactions/displacements.  
This is essentially an integrated version of the original [ show_bumps script](/index.php/Show_bumps "Show bumps") by Thomas Holder.  
Note that a dedicated **display vdW strain** function is available to perform this in retrospect, or on unmodified models.   
PyTMs also enables the user to sub-select surface accessible residues. This corresponds to an integrated version of [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") by Jason Vertrees. This option is enables by providing the required cutoff (in A^2).  | **Animation of optimization with steric vdW clashes:**  
[![pytms example: optimization, animation, clashes](/images/0/08/Pytms_optimization_animation.gif)](/index.php/File:Pytms_optimization_animation.gif "pytms example: optimization, animation, clashes")  
---|---  
  
  


  


## Hints

  * console output appears first in the PyMOL console prior to the API interface
  * when using the GUI menu, the modification process can be followed in the display window



  


## Update notes
    
    
    '''
        0.90 pre release testing
        1.00 first release
            * minor changes and typo correction
            * fixed non-incentive users receiving
              error messages related to the 'alter' command
            * new feature: integrated surface selection
            * new feature: integrated display of vdW clashes for all PTMs
            * new function: display of vdW clashes indpendent of modification
        1.1 Minor fixes
            * changed boolean processing of 'delocalized' keyword
        1.2 Additional function and minor improvements
            * fixed crashed related to random selections
            * new function: nitrosylate (Cysteine S-Nitrosylation)
            * updated usage descriptions
            * the user interface window is now resizable
            * changes in nitrate-function:
                * the torsion angle can now be defined for the nitro group and has
                  a default of ~22.352 (slight angle);
                  this feature is currently only supported for Tyrosines!
                * selecting a surface cutoff in conjuction with random placement
                  of CE1 or CE2 for Nitro-tyrosines will now select the most
                  accessible CE atom rather than a random one
            * font size can now be adjusted from the Main menu
    '''
    

## References

The original citation can be found on [PubMED](http://www.ncbi.nlm.nih.gov/pubmed/25431162) A current (non-repository) version of PyTMs can be downloaded here: [pytms.py](https://github.com/andwar/Pymol-script-repo/blob/master/plugins/pytms.py)  
In case of suggestions, questions or feedback please feel free to [ contact me](/index.php/User:Andwar "User:Andwar"). 

  


## SEE ALSO

  * [Git intro](/index.php/Git_intro "Git intro")
  * [Plugin_manager](/index.php/Plugin_manager "Plugin manager")
  * [Plugins](/index.php/Plugins "Plugins")



Retrieved from "[https://pymolwiki.org/index.php?title=Pytms&oldid=12231](https://pymolwiki.org/index.php?title=Pytms&oldid=12231)"


---

## Quickdisplays

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/quickdisplays.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/quickdisplays.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The module **quickdisplays** is a package containing several functions for quick standard displays:  


List of functions   
---  
Function  | What it does  | Usage   
disp_list  | prints a list of functions  | disp_list   
disp_ss  | secondary structure  | disp_ss [ selection [, colors [, only ]]]   
disp_ball_stick  | balls and sticks  | disp_ball_stick [ selection [, hydrogens [, only ]]]   
disp_mesh  | surface mesh  | disp_mesh [ selection [, color_m [, hydrogens [, only [, limits ]]]]]   
disp_surf  | surface  | disp_surf [ selection [, color_s [, transparency [, hydrogens [, solvent [, ramp_above [, only [, limits ]]]]]]]]   
disp_putty  | putty b-factor sausage  | disp_putty [ selection [, limits [, only ]]]   
  
  * Note that each function has an **individual help description** , callable by entering, e.g.:


    
    
    help disp_surf
    

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



  


Examples   
---  
[![Quickdisplays disp surf.png](/images/3/3d/Quickdisplays_disp_surf.png)](/index.php/File:Quickdisplays_disp_surf.png)**disp_surf**  
 _surface_ | [![Quickdisplays disp ss.png](/images/a/a8/Quickdisplays_disp_ss.png)](/index.php/File:Quickdisplays_disp_ss.png)**disp_ss**  
 _secondary structure_ | [![Quickdisplays disp putty.png](/images/5/51/Quickdisplays_disp_putty.png)](/index.php/File:Quickdisplays_disp_putty.png)**disp_putty**  
 _putty b-factor sausage_ |  The _combined displays_ example was produced like this: 
    
    
    import quickdisplays
    fetch 1hpv, async=0
    
    disp_putty all, limits=10
    disp_surf color_s=lightblue, transparency=0.8
    disp_ss chain B
    disp_ball_stick hetatm
    util.cbao hetatm
    set mesh_mode, 1 # make sure hetams are meshed
    set mesh_cutoff, 4.5
    disp_mesh resn 478, color_m=green
      
  
[![Quickdisplays disp mesh.png](/images/0/07/Quickdisplays_disp_mesh.png)](/index.php/File:Quickdisplays_disp_mesh.png)**disp_mesh**  
 _mesh_ | [![Quickdisplays disp ball stick.png](/images/4/4b/Quickdisplays_disp_ball_stick.png)](/index.php/File:Quickdisplays_disp_ball_stick.png)**disp_ball_stick**  
 _balls and sticks_ | [![Quickdisplays combo.png](/images/f/fc/Quickdisplays_combo.png)](/index.php/File:Quickdisplays_combo.png)_combined displays_  
see example   
  
## Some notes on the arguments

  * **selection** can be used to restrict the display to certain objects or selections, the default is 'all'
  * **only** can be _True_ or _False_ and will toggle whether the display will be added (_show_) or replace others (_show_as_)
  * **disp_mesh** and **disp_surf** support the color **putty** , which will color the object/selection by b-factor
  * **limits** defines the b-factor color range: 
    * a list entry will define absolute values **[min;max]**
    * a float value will calculate the corresponding percentiles **±value%** , _default=5_
  * setting **hydrogens** will add hydrogen to the model; if set to _=False_ , hydrogens are removed, if omitted the state of hydogens will be 'as is'
  * **colors** in **disp_ss** is flexible: 
    * set e.g. three colors, e.g. 'red green blue' for sheets, helices and loops, respectively (cf. example)
    * alternatively certain 'util.' functions can be used (cf. example)
    * setting _False_ once or for individual colors will suppress coloring



  


## Applied example

Enter the following lines one-by-one and observe how the display will be affected 
    
    
    # import function to PyMOL
    import quickdisplays
    
    # prep
    fetch 1hpv, async=0
    bg black
    orient
    disp_list # list of functions
    
    help disp_ss # check out help
    disp_ss
    disp_ss only=T
    disp_ss colors=smudge orange grey
    disp_ss colors=chartreuse False False # False suppresses coloring
    disp_ss colors=util.chainbow # You can use util.cbc, util.rainbow, ...
    
    help disp_ball_stick # check out help
    disp_ball_stick hetatm
    disp_ball_stick all, only=T
    util.cbaw
    disp_ball_stick hydrogens=1
    disp_ball_stick hydrogens=-1
    disp_stick_ball # will this work?
    
    help disp_mesh # check out help
    disp_mesh
    disp_mesh color_m=green, only=False
    disp_mesh color_m=green, only=T
    disp_ball_stick hetatm
    set mesh_skip, 2 # omits lines in mesh
    disp_mesh color_m=default, hydrogens=1, only=T
    disp_mesh color_m=putty, hydrogens=0, limits=20
    disp_mesh color_m=putty, limits=30
    disp_mesh color_m=putty, limits=[15,50]
    
    help disp_putty # check out help
    disp_putty chain A, only=False, limits=[15,50] #only default is True
    disp_putty all, limits=0
    disp_putty all, limits=10
    
    help disp_surf # check out help
    disp_surf color_s=white, solvent=1
    disp_surf color_s=white# solvent=0
    disp_surf color_s=lightblue, transparency=0.8
    

  


## SEE ALSO

[Git intro](/index.php/Git_intro "Git intro"), [Displaying Biochemical Properties](/index.php/Displaying_Biochemical_Properties "Displaying Biochemical Properties")

Retrieved from "[https://pymolwiki.org/index.php?title=Quickdisplays&oldid=13947](https://pymolwiki.org/index.php?title=Quickdisplays&oldid=13947)"


---

## Solvent radius

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This defines the solvent radius. Changing the solvent radius in PyMOL affects how PyMOL interprets [surfaces](/index.php/Surface_solvent "Surface solvent"), [surface area](/index.php/Displaying_Biochemical_Properties#Surface-Related "Displaying Biochemical Properties"). 

The default solvent radius is 1.4 Angstroms. 

# Syntax
    
    
    # set the solvent radius
    set solvent_radius, float
    
    # for example, set it to 2.0 Angstroms
    set solvent_radius, 2.0
    

# See Also

[Surface_solvent](/index.php/Surface_solvent "Surface solvent")

Retrieved from "[https://pymolwiki.org/index.php?title=Solvent_radius&oldid=6195](https://pymolwiki.org/index.php?title=Solvent_radius&oldid=6195)"


---

## Surfaces and Voids

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Help with integrating output from programs that generate/analyse surfaces and voids in molecules. 

See also: 

  1. [http://sourceforge.net/search/?group_id=4546&words=cavity&type_of_search=mlists&limit=30](http://sourceforge.net/search/?group_id=4546&words=cavity&type_of_search=mlists&limit=30)
  2. <http://loschmidt.chemi.muni.cz/caver/>
  3. <http://alpha2.bmc.uu.se/usf/voidoo.html>
  4. <http://sts.bioengr.uic.edu/castp/>
  5. <http://www.ccl.net/cca/software/UNIX/pass/overview.shtml>



Retrieved from "[https://pymolwiki.org/index.php?title=Surfaces_and_Voids&oldid=5396](https://pymolwiki.org/index.php?title=Surfaces_and_Voids&oldid=5396)"


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

## Vdw

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

vdw is an atomic property available in PyMOL that defines the atomic radius. It can be changed using the [alter](/index.php/Alter "Alter") command. 

Retrieved from "[https://pymolwiki.org/index.php?title=Vdw&oldid=6275](https://pymolwiki.org/index.php?title=Vdw&oldid=6275)"


---

