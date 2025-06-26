# Category: Electrostatics

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

## APBS Electrostatics Plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![APBS-Electrostatics-Plugin-Options.png](/images/6/63/APBS-Electrostatics-Plugin-Options.png)](/index.php/File:APBS-Electrostatics-Plugin-Options.png)

The APBS Electrostatics Plugin integrates the [APBS](http://www.poissonboltzmann.org) software package into PyMOL. Its primary purpose is electrostatic surface visualization. It superseedes the [APBS Tool2.1](/index.php/Apbsplugin "Apbsplugin") plugin. 

_New in Incentive PyMOL 2.0_

## Contents

  * 1 Usage
  * 2 Procedure
    * 2.1 Prepare Molecule
    * 2.2 Calculate Map with APBS
    * 2.3 Molecular Surface Visualization
  * 3 Troubleshooting
    * 3.1 RNA
    * 3.2 Incomplete Residues
  * 4 APBS Template
  * 5 Examples
    * 5.1 RNA Molecule
    * 5.2 Virus Assembly
  * 6 Python API
  * 7 See Also



## Usage

_A molecule must already be loaded into PyMOL._

For regular protein and DNA molecules, the plugin should be a one-click solution (just click "Run"). 

For RNA, the residue names might need adjustment, see troubleshooting RNA. 

Non-protein/DNA/RNA molecules or modified residues need special preparation, see #Prepare Molecule. 

## Procedure

The default procedure involves three steps: 

  1. preparing the molecule by assigning partial charges and adding hydrogens and other missing atoms
  2. calculating the electrostatic map
  3. visualization of the charged molecular surface



Each step can be run individually by unchecking the checkboxes on the other steps. 

### Prepare Molecule

APBS requires **partial charges** and **atom radii**. Since PDB files don't provide this information, they have to be preprocessed first. 

**pdb2pqr** is the default preparation method, it is limited to protein and nucleic acid, unless the user provides his own forcefield parameters. 

The **prepwizard** method is available if the user has the [Schrodinger Suite](https://schrodinger.com) installed. It can handle a variety of molecules. 

You can also use any **third-party tool** that assigns partial charges and load the result into PyMOL. Supported file formats include PQR, MOL2, and MAE. In that case, uncheck "Prepare Molecule", or use "Method: use vdw" if you only have charges but no radii. Tools to consider are for example [Open Babel](http://openbabel.org/) or the [PDB2PQR web server](http://www.poissonboltzmann.org/). 

### Calculate Map with APBS

  * Focus Selection: If you only care about e.g. a binding site, you can specify a [PyMOL selection](/index.php/Selection_Algebra "Selection Algebra") here to speed up the calculation.


  * Grid spacing (Angstrom): Increase to speed up the calculation and to reduce memory usage. Decrease to get a higher map resolution.



[![](/images/5/57/Ramp-levels-menu.png)](/index.php/File:Ramp-levels-menu.png)

[](/index.php/File:Ramp-levels-menu.png "Enlarge")

changing ramp level after APBS calculation

### Molecular Surface Visualization

The most common visualization is to show the solvent excluded surface. PyMOL uses the [potential at the solvent accessible surface](/index.php/Surface_ramp_above_mode "Surface ramp above mode") for coloring in that case. 

The **range** (color intensity) can also be changed later once the result is in PyMOL. Click "run01 > apbs_ramp01 > A > levels > ..." in the object menu panel. 

## Troubleshooting

### RNA

For **pdb2pqr** , RNA resdiue names must be RA, RC, RG, and RU. 
    
    
    alter polymer & resn A+C+G+U, resn = "R" + resn
    

### Incomplete Residues

Electrostatics should be analyzed on a complete molecule without missing sidechains, and ideally without missing loops. **pdb2pqr** will automatically model up to 10% of missing sidechain atoms. If your structure is incomplete, consider creating a homology model first using a tool of your choice (e.g. [MODELLER](https://salilab.org/modeller/)). 

Some truncated PDB files include a single backbone atom of the next residue, e.g. [2xwu](https://www.rcsb.org/structure/2xwu) chain B residue 954 atom N. **pdb2pqr** reports: 

> _Error encountered: Too few atoms present to reconstruct or cap residue LEU B 954 in structure!_

The easiest solution is to remove that atom: 
    
    
    remove /2xwu//B/954
    

## APBS Template

Parameters like temperature, protein and solvent dielectric, or ion concentrations can be changed directly in the "APBS Template". Please read the [APBS input file documentation](http://apbs-pdb2pqr.readthedocs.io/en/latest/apbs/input/). 

## Examples

### RNA Molecule
    
    
    fetch 1rna, async=0
    alter polymer & resn A+C+G+U, resn = "R" + resn
    

  * Open plugin
  * click "Run"



### Virus Assembly
    
    
    set assembly, 1
    fetch 3j7l, async=0
    split_states 3j7l
    delete 3j7l
    

  * Open plugin
  * Selection: 3j7l_*
  * Calculate Map with APBS > Options >> Grid Spacing: 2.0
  * click "Run" _(takes about 20 minutes!)_



## Python API

_Disclaimer: The plugin code is not part of the official**cmd** API and could be subject to change._

Example with PDB [1ubq](https://www.rcsb.org/structure/1ubq): 
    
    
    from pmg_tk.startup.apbs_gui.creating import pdb2pqr_cli
    from pmg_tk.startup.apbs_gui.electrostatics import map_new_apbs
    
    cmd.fetch("1ubq")
    
    pdb2pqr_cli("prepared01", "1ubq", options=["--ff", "amber"])
    map_new_apbs("apbs_map01", "prepared01")
    
    cmd.ramp_new("apbs_ramp01", "apbs_map01", [-5, 0, 5])
    cmd.set("surface_ramp_above_mode", 1, "prepared01")
    cmd.set("surface_color", "apbs_ramp01", "prepared01")
    cmd.show("surface", "prepared01")
    

## See Also

  * [apbsplugin](/index.php/Apbsplugin "Apbsplugin")
  * [APBS](/index.php/APBS "APBS")
  * [ramp_new](/index.php/Ramp_new "Ramp new") command



Retrieved from "[https://pymolwiki.org/index.php?title=APBS_Electrostatics_Plugin&oldid=13556](https://pymolwiki.org/index.php?title=APBS_Electrostatics_Plugin&oldid=13556)"


---

## H Add

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[h_add](/index.php/H_add "H add") uses a primitive algorithm to add hydrogens onto a molecule. 

# Usage and Algorithm Notes

PyMOL fills missing valences but does no optimization of hydroxyl rotamers. Also, many crystal structures have bogus or arbitrary ASN/GLN/HIS orientations. Getting the proper amide rotamers and imidazole tautomers & protonation states assigned is a nontrivial computational chemistry project involving electrostatic potential calculations and a combinatorial search. There's also the issue of solvent & counter-ions present in systems like aspartyl proteases with adjacent carboxylates . 

## Accessing Through GUI

[A]->hydrogens->add 

## Syntax
    
    
    # normal usage
    h_add (selection)
    
    # API usage
    cmd.h_add( string selection="(all)" )
    

see also [H_Fill](/index.php/H_Fill "H Fill")

Retrieved from "[https://pymolwiki.org/index.php?title=H_Add&oldid=6849](https://pymolwiki.org/index.php?title=H_Add&oldid=6849)"


---

## H Fill

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

H_Fill removes and replaces hydrogens on the atom or bond picked for editing. 

## Contents

  * 1 USAGE
  * 2 NOTES
  * 3 PYMOL API
  * 4 SEE ALSO



# USAGE
    
    
    h_fill
    

# NOTES

This is useful for fixing hydrogens after changing bond valences. 

# PYMOL API
    
    
     
    cmd.h_fill()
    

# SEE ALSO

  * [Edit](/index.php/Edit "Edit")
  * [Cycle_Valence](/index.php/Cycle_Valence "Cycle Valence")
  * [H_add](/index.php/H_add "H add")



Retrieved from "[https://pymolwiki.org/index.php?title=H_Fill&oldid=7424](https://pymolwiki.org/index.php?title=H_Fill&oldid=7424)"


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

