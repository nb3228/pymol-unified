# Category: Plugins

## 3Dstrut

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Plugin   
---|---  
Download  | <https://github.com/raghuyennamalli/3Dstrut>  
Author  | Tony Kamenick, Stuart Ballard   
Other authors  | George N Phillips Jr.   
Other contributors  | Shashank Ravichandran, Yatindrapravanan Narasimhan, Gautham Krishna, [Ragothaman Yennamalli](/index.php/User:Ragothaman "User:Ragothaman")  
  
3Dstrut is an external plugin for the molecular visualization program PyMOL. The plugin generates struts at optimal positions and thus enables 3D printing of any structure which can be visualized using PyMOL. This software is distributed freely for Windows and MacOS platforms. 

## Contents

  * 1 Website
  * 2 Version history
  * 3 Installation
  * 4 Tutorial



## Website

<https://github.com/raghuyennamalli/3Dstrut>

## Version history

Version 2.0 

## Installation

The following steps are applicable to both windows and MacOS systems 

I. Download the ‘.zip’ file and save it in your home directory. Once you unzip the file, you should see a directory named ‘3D_strut_plugin’ containing three files: 
    
    
       __init__.py, struts.py, build_struts.py
    

  
II. Open PyMOL and go to Plugin → Plugin manager. 

III. Select the ‘Install new Plugin’ tab and click on the choose button, a window will popup. 

IV. Navigate to where you extracted the ‘.zip’ file and select the ‘__init__.py’ file. Click open. 

V. Make sure that the path in the ‘Select Plugin directory’ popup box is similar to the one in the figure below and then click ok. 

VI. A new window will pop up with the message ‘Plugin “3D_Strut_plugin” has been installed.” and a new menu named “3-D printer” will be added to the menu bar. If you do not see the new menu immediately, close the plugin manager window quit and re-open PyMOL. 

## Tutorial

  * Open a PDB file of your choice


  * Navigate to 3D Printer → Build Struts (CA)


  * A popup will ask you to enter a value for the maximum length a strut. Click OK to proceed to next step/popup


  * Then, enter the value for the minimum distance between struts


  * If you want to build struts for everything opened in PyMOL, type “all”. Otherwise, if you want only one chain to be built, select that particular chain, make a new object, and hide the rest of the protein. In this case, you will type “visible”


  * You should be able to see the structs for the model in white color


  * You should look for spots where building one more strut would make the printed model structurally stable. For example, in this case this loop will be unstable after printing and will crumble during processing of the model


  * For the additional struts to be built, change the selection mode from “Residues” to “Atoms”


  * Pick the two atoms


  * And press F1. (In some keyboards the “fn” button needs to be pressed while pressing F1)


  * You can see view the same tutorial with screenshots here: <https://github.com/raghuyennamalli/3Dstrut/blob/main/Readme.pdf>



Retrieved from "[https://pymolwiki.org/index.php?title=3Dstrut&oldid=13454](https://pymolwiki.org/index.php?title=3Dstrut&oldid=13454)"


---

## Annocryst

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/annocryst.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/annocryst.py)  
Author(s)  | Anna Gerber   
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Introduction

AnnoCryst for PyMOL extends the functionality of PyMOL to support collaborative annotation of 3D crystallographic models. It is implemented as a plugin to PyMOL that provides an interface to browse and create annotations on structures loaded from Protein Data Bank (PDB) files. The annotations are retrieved from and stored on an annotation server using the W3C's Annotea protocol. 

Website: [AnnoyCryst](http://www.itee.uq.edu.au/eresearch/projects/annocryst-pymol) online. 

## Example Pymol Script

See userguide <http://www.itee.uq.edu.au/eresearch/downloads/annocryst-pymol-user-guide.pdf>

Standard settings of **AnnoCryst Settings**
    
    
    annotationsServerURL: http://maenad.itee.uq.edu.au:8080/Annotea/AnnoteaServlet
    keywordOntologyNamespace: http://www.proteinontology.info/po.owl
    keywordOntologyURL: http://maenad.itee.uq.edu.au/agerber/po.owl
    pdbRepositoryURL: http://maenad.itee.uq.edu.au:8080/harvanapdb/au.edu.uq.itee.eresearch.harvana.gwt.Main/pdb/
    uploadServerURL: http://maenad.itee.uq.edu.au:8080/Annotea/FileUploadServlet
    username: Anonymous
    

Change following settings, and save them. 
    
    
    pdbRepositoryURL: <http://www.rcsb.org/pdb/files/>
    

Then in the **PyMOL command window** , write 
    
    
     remotepdb 3ait
    

In the **Annocryst** , is now written: <http://maenad.itee.uq.edu.au:8080/harvanapdb/au.edu.uq.itee.eresearch.harvana.gwt.Main/pdb/3ait.pdb>

Push Open. 

Retrieved from "[https://pymolwiki.org/index.php?title=Annocryst&oldid=12516](https://pymolwiki.org/index.php?title=Annocryst&oldid=12516)"


---

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

## Autodock plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/autodock_plugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/autodock_plugin.py)  
Author(s)  | [Daniel Seeliger](/index.php/User:Dseelig "User:Dseelig")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Installation
    * 2.1 Linux
    * 2.2 Windows
  * 3 Tutorial Video
  * 4 Important
  * 5 Example 1 - Rename of funny atom names
  * 6 Example 2 - Alternative configurations of a Ligand
  * 7 References and LICENSES
    * 7.1 Vina
    * 7.2 MGLTOOLS



## Description

This plugin should help to set up docking runs with Autodock and view docking results. It has two features: 

  1. Setup of a docking grid for Autodock with PyMOL.
  2. View the docking results.



## Installation

### Linux

This plugin is ready "out-of-box" for Linux users through the project [ Pymol-script-repo](/index.php/Git_intro "Git intro")

### Windows

This plugin is ready "out-of-box" for win users through the project [ Pymol-script-repo](/index.php/Git_intro "Git intro")

You can download it manually 

  1. Download [plugin](http://wwwuser.gwdg.de/~dseelig/adplugin.html)
  2. **PyMOL > Plugin > Install Plugin**
  3. [MGLTools/AutoDockTools ](http://autodock.scripps.edu/downloads/resources/adt/index_html)
  4. [AutoDock4.2](http://autodock.scripps.edu/downloads/autodock-registration/autodock-4-2-download-page/)
  5. [AutoDock Vina](http://vina.scripps.edu/download.html)



## Tutorial Video

Watch Dan Seeliger's [autodock plugin tutorial](http://wwwuser.gwdg.de/~dseelig/plugin_tutorial.html). 

## Important

  1. The Autodock tools, does not like funny atoms names like "C1, N13, O28" and so on. Rename them! See example 1.
  2. Ligands can not be in alternative configuration. Create by: **create ROC_A, 1HXB and resn ROC and alt a** . See example 2.



## Example 1 - Rename of funny atom names

Read about the protein here: <http://www.proteopedia.org/wiki/index.php/3ig7>   
The example of video <http://wwwuser.gwdg.de/~dseelig/plugin_tutorial.html>   
Note, that the module "prepare_ligand4.py" does not like funny names of atoms, so we have to rename them 

Download: [examples/autodock_plugin_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/autodock_plugin_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/autodock_plugin_1.pml>" highlight="python" />

Open the Autodock/Vina plugin. Check/set the executable to folders (they should be automatically found): 
    
    
    Autodock Tools:       /home/you/another/pymol/Pymol-script-repo/modules/ADT/AutoDockTools/Utilities24
    autogrid4 executable: /home/you/another/pymol/Pymol-script-repo/modules/autodock_423/ia64Linux2/autogrid4                     # Platform dependent
    autogrid4 executable: /home/you/another/pymol/Pymol-script-repo/modules/autodock_423/ia64Linux2/autodock4                     # Platform dependent
    vina executable:      /home/you/another/pymol/Pymol-script-repo/modules/autodock_vina/autodock_vina_1_1_2_linux_x86/vina      # Platform dependent
    Working Directory:    /home/you/another/pymol/workingdir
    

Enlarge the window and click **Save Plugin Configuration file**. Plugin will not work, if you do not save. 

Click tab **Grid settings**. In **Calculate Grid Center by Selection** , write **EFP**.   
Try clicking **Wired Box** in top right, and then **Show Box** , **Hide Box** , **Cylindric Box** , **Show Box**.  
Then push **Select binding site**. In PyMOL, write 
    
    
    show sticks, binding_site
    cmd.disable("binding_site")
    hide sticks, binding_site
    

Click tab **Receptor**. Mark **cdk** , and click **Generate Receptor**.  
Click cdk2 to the "right", and see where the files are located.  
Click tab **Ligands**. Mark **EFP** , and click **EFP** in the Ligand list. Here is the files stored. Click **Generate Ligand**  
Click tab **Docking** , and click **Run Vina**. Wait 3-5 minutes, until **Writing output ... done.**  
Click tab **View poses** , Browse to **EFP.docked.pdbqt**. Click **Load** , **Show best 10**. Experiment with "Display options, and click the candidates"  
Congratulations, you did it! 

## Example 2 - Alternative configurations of a Ligand

Read about the protein here: <http://www.proteopedia.org/wiki/index.php/1hxb>   
Ligands can not be in alternative configuration. Create by: create ROC_A, 1HXB and resn ROC and alt a . 

Download: [examples/autodock_plugin_2.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/autodock_plugin_2.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/autodock_plugin_2.pml>" highlight="python" />

## References and LICENSES

### Vina

If you used AutoDock Vina in your work, [please cite:](http://vina.scripps.edu/)  
O. Trott, A. J. Olson, _AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading_ , **Journal of Computational Chemistry** 31 (2010) 455-461 

### MGLTOOLS

LICENSE of MGLTOOLS: 
    
    
    This software is copyrighted by Michel F. Sanner (sanner@scripps.edu) and TSRI.
    The following terms apply to all files associated with the software 
    unless explicitly disclaimed in individual files.
    
           MGLTOOLS SOFTWARE LICENSE AGREEMENT.
    
      1. Grant Of Limited License; Software Use Restrictions. The programs
         received by you will be used only for NON COMMERCIAL purposes.
         This license is issued to you as an individual.
    
         For COMMERCIAL use done with the software please contact Michel F. 
         Sanner for details about commercial usage license agreements.
         
         For any question regarding license agreements, please contact
                 Michel Sanner:
    	     TSRI, Molecular Biology Department, TCP 26,
    	     10550 North Torrey Pines Road, La Jolla, CA 92037
    	     sanner@scripps.edu
    	     tel (858) 784-7742
    	     fax (858) 784-2341
    
      2. COMMERCIAL USAGE is defined as revenues generating activities. These
         include using this software for consulting activities and selling
         applications built on top of, or using this software. Scientific 
         research in an academic environment and teaching are considered 
         NON COMMERCIAL.
    
      3. Copying Restrictions. You will not sell or otherwise distribute commercially 
         these programs or derivatives to any other party, whether with or without 
         consideration.
    
      4. Ownership of Software. You will not obtain, and will not attempt to 
         obtain copyright coverage thereon without the express purpose written 
         consent of The Scripps Research Institute and Dr. Sanner.
    
      5. IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
         FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
         ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
         DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
         POSSIBILITY OF SUCH DAMAGE.
    
      6. THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
         INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
         FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
         IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
         NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
         MODIFICATIONS.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Autodock_plugin&oldid=10429](https://pymolwiki.org/index.php?title=Autodock_plugin&oldid=10429)"


---

## Azahar

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/BIOS-IMASL/Azahar/archive/v0.8.8-beta.zip>  
Author(s)  | Agustina Arroyuelo and [Osvaldo Martin](/index.php/User:OsvaldoMartin "User:OsvaldoMartin")  
License  | [MIT](http://opensource.org/licenses/MIT)  
  
## Contents

  * 1 Description
  * 2 Installation
  * 3 Usage
    * 3.1 Glycan Creation
    * 3.2 Visualization
    * 3.3 Calculations
    * 3.4 Monte Carlo
  * 4 New templates
  * 5 Citation
  * 6 Change log



## Description

[![](/images/c/c6/Azahar_GUI.png)](/index.php/File:Azahar_GUI.png)

[](/index.php/File:Azahar_GUI.png "Enlarge")

Azahar GUI.

Azahar (pronounced /ɑːsɑːˈɑːr/) is a plugin that extends the PyMOL's capabilities to visualize, analyze and model glycans and glycoconjugated molecules. 

## Installation

The plugin have been tested on Linux and Windows and should also works on Mac Os X. To install the plugin just download [this](https://github.com/BIOS-IMASL/Azahar/archive/v0.8.8-beta.zip) zip file and install it using the plugin manager. 

Azahar can be used without any dependencies (i.e. you don't need to install anything else) but if you have installed OpenBabel (and OpenBabel Python bindings) Azahar will use it to perform certain calculations like optimize the geometry of your newly created molecule. OpenBabel is also used during the Monte Carlo with Minimization routine. 

If you still want to install OpenBabel on your computer see instructions [here](http://openbabel.org/wiki/Get_Open_Babel). Users of 64-bit Windows OS, please install OpenBabel 64-bit from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#openbabel)

## Usage

The plugin's GUI displays three interactive tabs: Creation, Visualization and Calculations. Usage is described by tab: 

### Glycan Creation

The Creation tab contains: 

  * Three option menus, two of which are used to specify residues to be connected from a list of predefined templates. And the third one, selects the chemical bond to be used.


  * Three entry boxes, the first one it's called "Repetitions", where the user can enter the number of residues to be connected. The "Position" entry box indicates the index of the residue to which the next bond will be added. On a typical glycan molecule constructed using this plugin, the user would use the "Position" entry box and the bonds menu to specify ramifications. And in the "New Molecule" entry box the user can write the name of the PyMOL object corresponding to the newly created glycan molecule.


  * Finally, The "add" button, when clicked, writes a text file containing the connectivity of the glycan molecule to be created. The user can edit this file or provide it as input. The connectivity instructions are also printed in PyMOL's external GUI. When the "create" button is clicked, the molecule is constructed and displayed in the PyMOL's GUI. The "reset" button resets the entry box's parameters.



### Visualization

On this tab the user will find options to visualize glycans with 3 different _cartoon-like_ representations. The _cartoon_ representation is build using PyMOL's [CGO](/index.php?title=CGO&action=edit&redlink=1 "CGO \(page does not exist\)") and hence the color of the object has to be define it before creating it. Azahar provides a color menu for that purpose, you can choose from a list of colors or use the _auto_ option (default), that reads the atom colors from the object. For example if you want to color a glycan cartoon using the b-factors, you should first choose that option from the PyMOL's menu and then tell Azahar to create a cartoon representation with the _auto_ color option. 

You can _cartoonize_ glycan molecules created by this plugin or loaded from other sources. 

### Calculations

Radius of gyration and Ramachandran Plot are useful tools for the analysis of glycan molecules. Through this plugin's GUI the user can easily select a range of states to compute these calculations, and in case of the RG, the result will be displayed in the PyMOL's external GUI. 

### Monte Carlo

This feature is still experimental 

A Monte Carlo with Minimization (MCM) is a global optimization method that have proven very efficient to sample the conformational states of complex molecules like proteins, It was first published by [Li and Scheraga](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC299132/), and is some circles is also know as basin hopping. Briefly, the conformational search is done by perturbing a torsional angle (psi, psi, omega or chi), follow by and energy minimization and the new (minimized) conformation is accepted or rejected following the Metropolis-criteria. If the conformation is accepted the next iteration start from the new conformation, otherwise the old conformation is perturbed again. 

From the GUI the user can choose the number of iterations to run, whether or not to use a SASA model to take into account the solvent contribution, and whether or not to start from a random conformation. 

After the MCM procedure is completed the user will get: 

  * one .pdb file for each accepted conformation.
  * one multi-state .pdb file (automatically uploaded at the end of the run)
  * one .pdb file for the conformation with the minimum energy.
  * a plain text listing the energy of each saved conformation



  
The internal energy is computed using the MMFF94 force field (as implemented on open-babel) and the contribution of the solvent is approximated using a SASA model (the validity of the SASA model for glycans and the relative weight of the internal and salvation energies needs further test). 

## New templates

You can add your own templates by adding a .pdb file to the Azahar/db_glycans folder, that was created during the plugin installation. The residue number of the template should be 0 and atoms names should have certain names. When in doubt just take a look at one of the already available templates and do not hesitate to contact us. 

## Citation

If you find this plugin useful please cite it as: 

Arroyuelo, A., Vila, J.A. & Martin, O.A. J Comput Aided Mol Des (2016) 30: 619. doi:10.1007/s10822-016-9944-x 

## Change log

  * 2017-06-17 (Version 0.8.7) 
    1. Python3 compatible
    2. improved error messages when no HB are found



  


  * 2015-10-2 (Version 0.8.5) 
    1. Fix templates: reorder atoms, standardize nomenclature and minimize conformations
    2. Add new templates
    3. Add Monte Carlo with Minimization Routine
    4. fix cartoonize for 1-6 and 6-1 bonds
    5. add 2 new cartoon-like representations
    6. fix minor bugs, clean code and GUI


  * 2015-08-31 (Version 0.8.1) 
    1. Several bugs fixed, including problems when creating branches.


  * 2015-08-28 (Version 0.8) 
    1. First beta version and first public version (several bugs expected!).



Retrieved from "[https://pymolwiki.org/index.php?title=Azahar&oldid=13304](https://pymolwiki.org/index.php?title=Azahar&oldid=13304)"


---

## Bnitools

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/bnitools.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/bnitools.py)  
Author(s)  | [Georg Steinkellner](/index.php/User:Steinkeg "User:Steinkeg")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![](/images/d/d4/Bnitools.png)](/index.php/File:Bnitools.png)

[](/index.php/File:Bnitools.png "Enlarge")

Screenshot of the BNI Tools plugin menu

**BNI Tools** is a plug in for PyMOL which adds additional functionalities and presets to the PyMOL GUI. 

## Contents

  * 1 Installation
  * 2 Example images created using BNI Tools
  * 3 BNI Tools dropdown menu (Plugin > BNI Tools > ...)
  * 4 BNI tools integrated in PyMOL sidebar



## Installation

Install it using the 'Install Plugin' menu within PyMOL:  
**PyMOL > Plugin > Plugin Manager**  


## Example images created using BNI Tools

[![Create a box around a selection](/images/6/65/Bni_box.jpg)](/index.php/File:Bni_box.jpg "Create a box around a selection") [![Load and show DelPhi or other maps](/images/8/83/Bni_delphimap.jpg)](/index.php/File:Bni_delphimap.jpg "Load and show DelPhi or other maps") [![Create a box](/images/7/7a/Bni_box2.jpg)](/index.php/File:Bni_box2.jpg "Create a box") [![Exclude parts of surface](/images/0/02/Bni_excludesurf.jpg)](/index.php/File:Bni_excludesurf.jpg "Exclude parts of surface") [![Create Plane](/images/c/c0/Bni_plane.jpg)](/index.php/File:Bni_plane.jpg "Create Plane") [![Create a main chain track](/images/1/10/Bni_track.jpg)](/index.php/File:Bni_track.jpg "Create a main chain track")

## BNI Tools dropdown menu (Plugin > BNI Tools > ...)
    
    
        [+]Load Files->
         [-]modeller files
           * load multiple modeller files and sort by 
           * objective function energy 
           * and add energy to object title
         [-]autodock files
           * load autodock .dlg or .dlg.pdb file into
           * different states with energy in object 
           * title
         [-]amber minimization
           * load amber minimization file with energy in
           * title if the .info file is present in same
           * directory and has the same name as the pdb file.
         [-]delphi phi,dx map
           * load delphi map and corresponding pdb file
           * simultaneously and show surface colored by 
           * PHI or DX map. (show the surface to see the 
           *  effect)
         [-]casox map
           * load casox map (cavity calculation) and 
           * show ligsite accessibility 
           * value maps. (7 closed cavity to 1 open)
         [-]multiple files into states
           * load multiple pdb files (e.g. MD simulation 
           * snapshots) into one state. (object is named
           * by first object loaded)
        [+]Fetch->
         [-]Density View (EDS)
           * load density and pdb file from EDS
           * (Electron Density Server)
           * if available, and show density with density
           * wizard
         [-]RCSB Biol. Assembly
           * load biological assembly from RCSB
           * protein database
         [-]2FoFc map(s)
           * load (multiple) 2FoFc maps
           * from EDS density server
           * if available
         [-]FoFc maps(s)
           * load (multiple) FoFc maps
           * from EDS density server
           * if available
        [+]Edit->
         [-]HIS --> HID,HIE,HIP
           * change histidine residues to HID,HIE,HIP
           * depending on hydrogens on histidine
         [-]HID,HIE,HIP --> HIS
           * change altered histidine residues 
           * back to HIS
         [-]Poly-Alanine Chain
           * create a poly alanine chain (GLY and ALA)
           * for molecular replacement
         [-]MSE --> MET
           * change selenomethionine to methionine
         [-]del alternates
           * delete alternates in selection
         [-]Unbond- >
           * unbond atoms in selection
        [+]Images->
          * create ray traced images depending on
          * x size and resolution (dpi)
        [+]Create->
          * create compiled graphics objects (CGO)
          * these objects can be altered in color or
          * transparency, and they can be dragged
          * and rotated in space by the
          * "action->drag" command and using
          * "shift" and mouse buttons.
         [-]Plane
           * create a plane (with certain cushion)
           * using a selection of three atoms
         [-]Box
           * create a box around selection
           * the whole box can be altered as group
           * or by side planes separately
         [-]Triangle
           * create a triangle using a selection of three
           * atoms
        [+]pseudo center atom
            * create pseudo atom at the center of atoms
            * in selection
            * move atoms with editing->"shift"
            * and middle mouse button
    

## BNI tools integrated in PyMOL sidebar
    
    
        [+]Action on "all"
         [-]delete enabled
           * delete all enabled objects or 
           * selections
         [-]invert enabled/disabled
           * disable currently enabled objects
         [-]combine selections
           * combine all enabled or disabled 
           * selections to selection (sele)
        [+]action
         [-]sequence
           * show sequence in different formats
           * of selection
           * copy and paste to text file for later use
           [.]fasta
             * show sequence of selection
             * in fasta format
           [.]pir
             * show sequence of selection
             * in pir format
           [.]modeller
             * show sequence of selection
             * in modeller pir format
           [.]list
             * create residue or atom lists
             * of selection
        [+]preset
         [-]track main chain
           * create a new object which tracks the
           * main chain atoms and shows main chain
           * and side chain polar contacts
         [-]symmetry surface
           * create a symmetry view of the selected
           * atoms showing the contact surface as well
           * as a selection entry of the atoms in 
           * contact with symmetry mates
           * (only includes atoms of the initial selection).
         [-]hydrophobic residues
           * show hydrophobic residues
           * depending on hydrophobic residue scales
           * by KandD (Kyte & Doolittle 
           *           J Mol Biol 157:105, 1982)
           *    Rose  (Rose et al. 
           *           Science, 229, 834-838,1985)
           *    GES   (Engelman Engelman et al. 
           *           Annu Rev Biophys Biophys Chem,
           *           15, 321-353(1986)
           * (no window selection; just raw categories
           * are colored by: blue-hydrophile
           *                green-neutral
           *                red-hydrophobe )
         [-]surface inspection
           * create selections for surface
           * inspections
      # BNI tools additional settings in PyMOL sidebar
         [+]color
          [-]by ss
            * color helix,sheet and loop separately
            [.]Helix
            [.]Sheet
            [.]Loop
          # this section is replaced by "[-] by_rep" in PyMOL versions >1.6
          [-]surface
            * color surface separately from atoms
            [.]by atom
              * set surface color to standard
            [.]by map
              * if a map or ramp is loaded
              * color surface by ramp/map
            [.]my color
              * color surface by own defined color
          [-]mesh
            * color mesh separately from atoms
            [.]by atom
              * set mesh color to standard
            [.]by map
              * if a map or ramp is loaded
              * color mesh by ramp/map
            [.]my color
              * color mesh by own defined color
          [-]label
            * color labels separately
            [.]by atom
              * set label color by atom
            [.]my color
              * color labels by own defined color
          [-]stick
            * color sticks separately
            [.]standard
              * set stick color to standard
            [.]my color
              * color sticks by own defined color
          [-]my colors
            * use/append own defined colors
            * own colors can be defined by
            * Setting->Colors..->New
            * to keep color settings for
            * other pymol sessions
            * you have to set the colors
            * in .pymolrc or similar
            * pymol setting file
            * like
            * set_color mycolor,[ 1.00, 1.00, 1.00]
         [+]show
          [-]surface flag
            * set surface flag of atoms to show hide
            * or ignore fro surface calculation
          [-]transparency
            * set different transparency types
            * on selections or atoms
    

Retrieved from "[https://pymolwiki.org/index.php?title=Bnitools&oldid=12367](https://pymolwiki.org/index.php?title=Bnitools&oldid=12367)"


---

## Bondpack

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/rasbt/BondPack>  
Author(s)  | Sebastian Raschka   
License  | GNU GENERAL PUBLIC LICENSE   
  
## Contents

  * 1 BondPack
  * 2 Introduction
  * 3 HydroBond
  * 4 BondVis
  * 5 BondLab
  * 6 Github



## BondPack

A collection of PyMOL plugins to visualize atomic bonds.  


  


## Introduction

PyMOL is without any doubt a great tool to visualize protein and ligand molecules.  
However, drawing interactions between atoms can be often quite cumbersome when done manually.  
For the sake of convenience, I developed three plugins for PyMOL that will make our life as protein biologists a little bit easier.  
All three PyMOL plugins can be installed and used separately; they don't depend on each other, but rather complement each other.  
At the end of this article, you will find brief instructions on how to install plugins in PyMOL - a very quick and simple process.  


## HydroBond

HydroBond visualizes all potential polar contacts between protein and ligand molecules within a user-specified distance.  
The underlying function is based on the different atom types, such as hydrogen bond acceptor and donor atoms,  
and thus it is required to have hydrogen atoms present in the structure.   
If your structure file doesn't contain hydrogen atoms already, you can add them directly in PyMOL as shown in the image below.  
  
[![Add hydrogens.png](/images/e/e0/Add_hydrogens.png)](/index.php/File:Add_hydrogens.png)   
  
HydroBond is related to PyMOLs "[A]->find->polar contacts" command, however,  
it doesn't consider geometry criteria and angle thresholds,  
but is rather based on atom types.  
When you select HydroBond from the "Plugin" menu, you will be prompted to enter the name of the protein object,  
the ligand object, and a distance cutoff as shown in the figure below.  
If HydroBond was able to detect hydrogen bond and acceptor atoms within the  
specified distance, potential interactions will be visualized as yellow dashed lines.  
  
[![Hydrobond action.png](/images/8/8b/Hydrobond_action.png)](/index.php/File:Hydrobond_action.png)   
  


## BondVis

The BondVis plugin lets you visualize interactions between any pair of atoms you specified.  
Often I find it helpful for my understanding (and for verification) to visualize the bonds between certain atoms  
that were assigned in docking or any other prediction software.  
Most software will provide you with information about the atoms that were "connected" to perform the analysis.  
If you run BondVis from the "Plugin" menu, it will ask you to select a "bond info file."  
  
[![Bondinfo.png](/images/3/33/Bondinfo.png)](/index.php/File:Bondinfo.png)   
  
This is just a simple text file that contains the atom numbers of connected atoms in pairs.  
Those can be separated by spaces, tabs, or commas. An example file with bond information could look like this:  


`
    
    
    1174		1357
    1175		1358
    1176		1359
    

`

  
When you selected a "bond info" file, BondVis will connect all atom pairs by yellow dashed lines  
and print out the connected atoms in the command field.  
  
[![Bondvis.png](/images/d/db/Bondvis.png)](/index.php/File:Bondvis.png)   
  


## BondLab

If you are not happy with the looks of the lines that were drawn to visualize connections between atoms,  
you will like to use BondLab.  
This plugin offers a simple interface that allows you to change the diameter,  
gap width, and color of the lines.   
  
[![Bondlab.png](/images/b/b3/Bondlab.png)](/index.php/File:Bondlab.png)  
The following video demonstrates the different features of BondLab:  
<http://youtu.be/14UZctxtK3w>   


## Github

If you are interested, you can follow the BondPack Github repository  
<https://github.com/rasbt/BondPack> for updates 

Retrieved from "[https://pymolwiki.org/index.php?title=Bondpack&oldid=11553](https://pymolwiki.org/index.php?title=Bondpack&oldid=11553)"


---

## CASTp

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/castp.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/castp.py)  
Author(s)  | Joe Dundas   
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview

The CASTpyMOL plugin allows for the visualization of surfaces and voids identified from CASTp using pyMOL. For more information see the [CASTp Website](http://sts-fw.bioengr.uic.edu/castp/pymol.php). 

Retrieved from "[https://pymolwiki.org/index.php?title=CASTp&oldid=10430](https://pymolwiki.org/index.php?title=CASTp&oldid=10430)"


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

## Caver2

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[A new version of the Caver plugin is available](/index.php/Caver3 "Caver3")

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/Caver2_1_2.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/Caver2_1_2.py)  
Author(s)  | [CAVER 2.0](/index.php/User:Hci "User:Hci")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[CAVER 2.0](http://loschmidt.chemi.muni.cz/caver/) is a software tool for protein analysis leading to detection of channels. Channels are in fact void pathways leading from buried cavities to the outside solvent of a protein structure. Studying of these pathways is highly important in the area of drug design and molecular enzymology. 

[![Ukazky-plugin.png](/images/c/cb/Ukazky-plugin.png)](/index.php/File:Ukazky-plugin.png)

CAVER 2.0 provides rapid, accurate and fully automated calculation of channels not only for static molecules, but also for dynamic molecular simulations.  
CAVER 2.0 facilitates analysis of any molecular structure including proteins, nucleic acids or inorganic materials. 

[![Ukazky.png](/images/1/1a/Ukazky.png)](/index.php/File:Ukazky.png)

CAVER 2.0 can be used in two possible ways. If you are accustomed with the PyMOL application, you can download the CAVER plugin for PyMOL. Then the standard functionality of PyMOL application for the subsequent exploration of computed channels can be used. The other possibility is to take advantage of our CAVER Viewer application, which has been designed exactly for the problem of detection and mainly for further exploration of channels. It includes the CAVER 2.0 algorithm for computation of channels in static molecule as well as in molecular dynamics. CAVER Viewer can be downloaded and installed as a standalone version or you can use the online version based on the Java Web Start technology. 

References: 

Medek, Petr and Beneš, Petr and Kozlíková, Barbora and Chovancová, Eva and Pavelka, Antonín and Szabó, Tibor and Zamborský, Matúš and Andres, Filip and Klvaňa, Martin and Brezovský, Jan and Sochor, Jiří and Damborský, Jiří. "[CAVER 2.0](http://loschmidt.chemi.muni.cz/caver.html)", 2008. 

Medek, Petr and Beneš, Petr and Sochor, Jiří,"[Computation of tunnels in protein molecules using delaunay triangulation](http://decibel.fi.muni.cz/download/papers/medek08.pdf)". Journal of WSCG. Volume 15(1-3), 2007, Pages: 107-114. 

Retrieved from "[https://pymolwiki.org/index.php?title=Caver2&oldid=12624](https://pymolwiki.org/index.php?title=Caver2&oldid=12624)"


---

## Caver3

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

## Cheshift

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/aloctavodia/cheshift/archive/v3.6.zip>  
Author(s)  | [Osvaldo Martin](/index.php/User:OsvaldoMartin "User:OsvaldoMartin")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
  
  


## Contents

  * 1 Description
    * 1.1 Version
  * 2 Installation
    * 2.1 Linux
    * 2.2 Windows
    * 2.3 Mac OsX
  * 3 Using the Plugin for predicting 13Cα and 13Cβ chemical shifts
  * 4 Using the Plugin for protein structure validation
  * 5 License
  * 6 Change log
  * 7 Citation
  * 8 References



## Description

[![](/images/e/ed/Cheshift.png)](/index.php/File:Cheshift.png)

[](/index.php/File:Cheshift.png "Enlarge")

Result of a CheShift analysis.  
Colors indicate the difference between predicted and observed 13Cα and/or 13Cβchemical shifts values averaged over all uploaded conformers. Green, yellow and red colors represent small, medium and large differences, respectively. White is used if either the prediction failed or the observed value is missing.  
Blue is used to highlight residues for which the agreement between observed and predicted 13Cα and 13Cβchemical shifts can be improved, i.e., if the (χ1 /χ2 ) side-chain

CheShift (pronounced /tʃeʃɪft/) is a software for prediction of 13Cα and 13Cβ chemical shifts and validation of proteins structures. It is based on the idea that the differences between observed and predicted 13Cα and 13Cβ chemical shifts can be used as a sensitive probe with which to detect possible local flaws in protein structures. A [Web Server](http://www.cheshift.com) is also available. 

This plugin provides a way to use PyMOL to validate a protein model using observed chemical shifts. 

### Version

The current version of this plugin is 3.6 

## Installation

### Linux

1) The plugin can be downloaded from here [[Here](https://github.com/aloctavodia/cheshift/archive/v3.6.zip)] 

2) If you are using the incentive version of PyMOL skip the next step and go directly to step 4. 

3) You should install NumPy and SciPy. These Python packages are available from the repositories of the main Linux distributions. Just use your default package manager (or command line) to install it. In Ubuntu/Debian you could do "sudo apt-get install python-numpy python-scipy". 

4) Open PyMOL and then go to plugin --> plugin manager --> Install New Plugin --> Install from local file --> Choose file. and choose the zip file you download in step 1. Now the plugin should be installed. If you have an old version of PyMOL without the plugin manager then, unzip the file downloaded in step 1 and copy the "cheshift" folder to the plugin directory, probably something like "/usr/lib/python2.7/dist-packages/pmg_tk/startup" 

### Windows

This plugin have not been extensively tested on Windows machines, but it passed all the test I have done... 

1) The plugin can be downloaded from here [[Here](https://github.com/aloctavodia/cheshift/archive/v3.6.zip)] 

2) If you are using the incentive version of PyMOL skip the next step and go directly to step 4. 

3) You should install NumPy and SciPy. Probably the easiest way to do this on a Windows machine is to install a Python Distribution like [[Anaconda](https://store.continuum.io/cshop/anaconda/)] or [[Canopy](https://www.enthought.com/products/canopy/)] 

4) Open PyMOL and then go to plugin --> plugin manager --> Install New Plugin --> Install from local file --> Choose file. and choose the zip file you download in step 1. On Windows machine this step could take a couple of minutes, just wait for the confirmation message. Now the plugin should be installed. If you have an old version of PyMOL without the plugin manager then, unzip the file downloaded in step 1 and copy the "cheshift" folder to the plugin directory, probably something like "C:\Python27\Lib\site-packages\pmg_tk\startup\" 

### Mac OsX

This plugin have not been extensively tested on Windows machines, but it passed all the test I have done... 

1) The plugin can be downloaded from here [[Here](https://github.com/aloctavodia/cheshift/archive/v3.6.zip)] 

2) If you are using the incentive version of PyMOL skip the next step and go directly to step 4. 

3) You should install NumPy and SciPy. Probably the easiest way to do this on a Windows machine is to install a Python Distribution like [[Anaconda](https://store.continuum.io/cshop/anaconda/)] or [[Canopy](https://www.enthought.com/products/canopy/)] 

4) Open PyMOL and then go to plugin --> plugin manager --> Install New Plugin --> Install from local file --> Choose file. and choose the zip file you download in step 1. Now the plugin should be installed. If you have an old version of PyMOL without the plugin manager then, unzip the file downloaded in step 1 and copy the "cheshift" folder to the plugin directory, probably something like "/usr/lib/python2.7/dist-packages/pmg_tk/startup" 

## Using the Plugin for predicting 13Cα and 13Cβ chemical shifts

1) Launch PyMOL and select a PDB file   
2) Select "Cheshift" from the plugin menu.   
3) Click the "Run" button.  
4) The results will be saved in your working directory as a .txt file  


## Using the Plugin for protein structure validation

1) Launch PyMOL and select a PDB file   
2) Select "Cheshift" from the plugin menu.   
3) Select a file with the experimental chemical shift values.  
4) Click the "Run" button.  
5) Wait until results are displayed (this could take from seconds to a few minutes, depending on the size of your protein and the speed of your computer  


**Note:**

  * If the PDB file has more than one chain only the first one will be analysed.
  * The PDB file must have no missing residues.
  * Missing observed 13Cα and 13Cβ chemical shifts are tolerated.
  * A file with observed 13Cα and/or 13Cβ chemical shift values is needed. The format of this file should be the one used in the BMRB or the PDB. Alternatively, you can provide a file with the following format the first line should contain the name of the reference compound i.e DSS, TSP or TMS. The following lines should have four columns. The first column should be the residue number, the second the residue name (three-letter code) and the third column the 13Cα experimental chemical shifts and the last column the 13Cβ chemical shifts. Spaces should be used to separates the columns.



DSS  
1 MET 55.63 32.95  
2 TYR 62.81 39.27  
3 ALA 53.78 18.97  
4 GLY 47.24 999.00  
5 LYS 57.55 32.77  
6 ILE 56.38 38.59  


  


## License

_Che_ Shift plugin is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License. A complete copy of the GNU General Public License can be accessed [here](http://www.gnu.org/licenses). 

_Che_ Shift uses data derived from the [Neighbor-Dependent Ramachandran Distributions](http://dunbrack.fccc.edu/ndrd/) obtained by the Dunbrack Lab. This derived data is also released under the GPL license, with the permission of Professor Roland Dunbrack. The full NDRD data is released using a different [license](http://dunbrack.fccc.edu/ndrd/license/index.html). 

## Change log

  * 2014-07-02 (Version 3.6)



Fix bug with relative path and add computation of the conformationally averaged rmsd (ca_rmsd). The ca_rmsd is the root mean square deviation of the theoretical chemical shift and the observed chemical shifts. The ca_rmsd is calculated for both nuclei (13Cα and 13Cβ). If more than one conformation is used, the theoretical chemical shift for each residue is the average along all conformations. 

  * 2014-05-23 (Version 3.5)



Previous versions required an Internet connection. This is the first stand alone version, i.e. all computations are performed on the local machine. This version of the plugin performs the same computations the _Che_ Shift Server does. For details please read [1] 

  * 2013-08-28 (Version 3.0)



This version was able to establish a connection to the version of the _Che_ Shift Server described in [1] 

  * 2012-02-14 (Version 2.0)



This version was able to establish a connection to the version of the _Che_ Shift Server described in [2] 

## Citation

If you find this plugin useful please cite it as: 

Martin O.A. Vila J.A. and Scheraga H.A. (2012). _Che_ Shift-2: Graphic validation of protein structures. Bioinformatics 2012. 28(11), 1538-1539. 

## References

[1] Martin O.A. Arnautova Y.A. Icazatti A.A. Scheraga H.A. and Vila J.A. A Physics-Based Method to Validate and Repair Flaws in Protein Structures. Proc Natl Acad Sci USA 2013, vol. 110 16826-16831. 

[2] Vila J.A. Arnautova Y.A. Martin O.A. and Scheraga, H.A. (2009). Quantum-mechanics-derived 13C chemical shift server (_Che_ Shift) for protein structure validation. PNAS, 106(40), 16972-16977. 

[3] ﻿Vila, J.A. and Scheraga H.A. (2009). Assessing the accuracy of protein structures by quantum mechanical computations of 13Cα chemical shifts. Accounts of chemical research, 42(10), 1545-53. 

Retrieved from "[https://pymolwiki.org/index.php?title=Cheshift&oldid=12506](https://pymolwiki.org/index.php?title=Cheshift&oldid=12506)"


---

## Cluster mols

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![Cluster mols py pymol.png](/images/4/4c/Cluster_mols_py_pymol.png)](/index.php/File:Cluster_mols_py_pymol.png)

cluster_mols is a PyMOL plugin that allows the user to quickly select compounds from a virtual screen to be purchased or synthesized. 

It helps the user by automatically clustering input compounds based on their molecular fingerprints [[1]](http://openbabel.org/wiki/FP2) and loading them into the PyMOL window. cluster_mols also highlights both good and bad polar interactions between the ligands and a user specified receptor. Additionally there are a number of keyboard controls for selecting and extracting compounds, as well as functionality for searching online to see if there are vendors for a selected compound. 

## Contents

  * 1 Description
  * 2 Download
  * 3 Installation
  * 4 Usage
    * 4.1 GUI Options
    * 4.2 Keyboard Controls
    * 4.3 show_contacts
  * 5 Citing ClusterMols
  * 6 Authors



# Description

The basic work flow of cluster_mols.py can be broken up into three parts. 

  1. Computing a similarity matrix from the input compounds
  2. Performing hierarchical clustering on the results from 1)
  3. Cutting the tree at a user-specified height and creating and sorting clusters



The results of 1 and 2 are saved to python pickle files so you do not have to recompute them in subsequent runs. 

In addition, it also highlights both good and bad polar contacts between the ligand and a user specified protein using the 'show_contacts' module described below. 

This script also integrates keyboard controls which allows for WASD movement through the clusters as well as keyboard shortcuts for pulling out compounds. See below for usage. 

# Download

The most up to date version (recommended) of cluster_mols is available through BitBucket at: <https://bitbucket.org/mpb21/cluster_mols_py/overview>

# Installation

This plugin has a number of dependencies that are required. And it is currently only supported on Linux and OSX. 

Python packages (install using easy_install or pip) 

  1. openbabel
  2. numpy
  3. scipy
  4. Tkinter
  5. fastcluster
  6. Pmw-py3 **Important:** Pmw 2.0.1 does not work; install the Pmw-py3 package instead of Pmw to get version 2.1



Command line tools (These must be accessible through your PATH environment variable): 

  1. obabel -- from <http://openbabel.org>



Recent versions of cluster_mols do not require sdsorter, but it is still a very useful tool for dealing with sdf files. 

  1. sdsorter -- <https://sourceforge.net/projects/sdsorter/>



Once you have the required dependencies, install it through PyMOL's Plugin menu. 

PyMOL > Plugin > Install Plugin 

  


# Usage

The GUI is relatively straightforward, if you follow it from top to bottom, and then then left to right through the tabs. 

The program requires that the input be a '.sdf' or '.sdf.gz' file. If your compounds are not in that format, use the 'babel' tool from OpenBabel to convert them. 

## GUI Options

[![Cluster mols screen 1 desc.png](/images/e/ed/Cluster_mols_screen_1_desc.png)](/index.php/File:Cluster_mols_screen_1_desc.png)

[](/index.php/File:Cluster_mols_screen_1_desc.png "Enlarge")

[![Cluster mols screen 2 desc.png](/images/6/6d/Cluster_mols_screen_2_desc.png)](/index.php/File:Cluster_mols_screen_2_desc.png)

[](/index.php/File:Cluster_mols_screen_2_desc.png "Enlarge")

In the 'Compute Similarities' tab, there are options for selecting a new ligand and for specifying how many CPUs you want to run the similarity calculation on. Clicking the 'Compute Similarity' button will start the similarity calculations. If you check the 'Ignore saved results?' box it will ignore any saved intermediate results files. This could be useful if you change the contents of the original input file while keeping the file name the same. 

Depending on how many compounds there are, the similarity calculations may take between 1 and 10 minutes. If you launched PyMOL from the command line, you will be able to see the progress printing out in the console. The similarity results are saved to a file so if you want to re-cluster the same input file, you do not need to wait to recompute the similarities. 

  
The first option on the Cluster Compounds tab defines how the clusters will be sorted. The default is to sort by the 'minimizedAffinity' which is inserted into the output sdf file after minimization with 'smina' (An enhanced version of AutoDock Vina. Available at: <http://www.smina.sf.net>). You can also sort the clusters by any SD tag that exists in the input file, or by the Title (alphabetically) or by the size of the cluster. 

The second option is the height at which the hierarchical clustering tree is cut. The units are arbitrary, but a higher number leads to a small number of large clusters of less similar compounds, and lower cutoffs lead to more small clusters of more similar compounds. Play around with the cutoff until you get a clustering that you like. The third option is a check box for whether to group clusters with only one compound into one ‘singletons’ cluster. The forth option enables the show_contacts tool that is described below. There is also a field to enter a PyMOL selection string to compute the hydrogen bonds to. Finally, there is a button to create the clusters and load them into PyMOL. 

## Keyboard Controls

Once you have finished the similarity calculations and clustering mentioned above, you can navigate the clusters using the keyboard. Familiar to gamers, you can move through clusters using the WASD keys, (W for up, S for down, A for left, D for right). The one important caveat is that due to [limitations](/index.php/Set_Key#KEYS_WHICH_CAN_BE_REDEFINED "Set Key") in PyMOL, the WASD movement needs to be used with the Control (or Alt) key. Meaning Ctrl-W moves up. It seems weird, but you quickly get used to it. 

  
Navigation Controls 

Ctrl-W – Move up a cluster 

Ctrl-S – Move down a cluster 

Ctrl-A – Move to the previous compound in a cluster 

Ctrl-D – Move to the next compound in the cluster 

Ctrl-F -- Check for vendors 

If you acquired your compounds from ZINCPharmer (<http://zincpharmer.csb.pitt.edu/>) and/or your compounds have title that start with a ZINC ID (<http://www.docking.zinc.org>) or a MolPort ID (<http://www.molport.com>), you can hit 'Ctrl-F' to see if there are any vendors available. 

  
Compound selection 

In addition to moving through the clusters, you can also extract compounds that you like for later viewing using the following controls. Pressing F3 will append the current compounds into a new object with the suffix '_selected'. 

F1 – Print title of currently selected molecule 

F2 – Remove most recently added compound 

**F3 – Add currently visible compound to list** (Most commonly used) 

F4, F12 – Print List 

## show_contacts

show_contacts is an expanded version of list_hbonds[[2]](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/) that shows both favorable and unfavorable contacts between ligands and a protein receptor. show_contacts has been integrated into cluster_mols as a function and is executed automatically when clustering. It can also be run by itself, not in the context of cluster_mols. In the standalone case, the usage is as follows: 

show_contacts(selection,selection2,result="contacts",cutoff=3.6, bigcutoff = 4.0): 

The arguments are as follows: 

  1. selection -- pymol selection string for the protein
  2. selection2 -- pymol selection string for the ligands
  3. results -- prefix of the object that the distances should be shown in. (Default "contacts")
  4. cutoff -- Distance cutoff for what is considered an ideal hydrogen bond.
  5. bigcutoff -- Distance cutoff for a non-ideal hydrogen bond.



  
Output: The output of show_contacts are a set of pymol distance objects. They are color-coded and size coded to indicate different interactions between the ligand and protein. They are controlled by the parameter indicated. 

  1. thin-purple lines -- all possible polar contacts (acc-acc, don-don, acc-don) -- bigcutoff
  2. thick-yellow lines -- All ideal hydrogen bonds -- cutoff
  3. thin-yellow lines -- Non ideal hydrogen bonds -- bigcutoff
  4. thick-red lines -- Polar clashes, i.e. Donor-Donor, Acceptor-Acceptor -- cutoff



  


# Citing ClusterMols

If you use ClusterMols in your work, please cite the following. 

Baumgartner, Matthew (2016) IMPROVING RATIONAL DRUG DESIGN BY INCORPORATING NOVEL BIOPHYSICAL INSIGHT. Doctoral Dissertation, University of Pittsburgh. 

# Authors

The main cluster_mols.py script was conceived of by Matthew P Baumgartner (mpb21 [at] pitt.edu) and Dr. David Koes while working in the lab of Dr. Carlos Camacho at the University of Pittsburgh. The cluster_mols.py script was implemented (and later rewritten) by MPB. The show_contacts functionality and the first version of the objectfocus.py keyboard controls was written by DK. 

Please send questions/comments/bug reports to matthew.p.baumgartner [at] gmail.com. 

Retrieved from "[https://pymolwiki.org/index.php?title=Cluster_mols&oldid=13488](https://pymolwiki.org/index.php?title=Cluster_mols&oldid=13488)"


---

## CMPyMOL

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Software   
---|---  
Download  | <https://github.com/emptyewer/CMPyMOL>  
Author(s)  | [Venkatramanan Krishnamani](/index.php/User:Venky "User:Venky")  
License  | The MIT License (MIT)   
  
CMPyMOL is an add-on software to molecular visualization program PyMOL. It combines the protein 3D visualization capabilities of PyMOL and the protein's 2D contact map with an interactive interface for scientific analysis. This software is freely distributed under the MIT license for Linux and Mac OS X platforms. 

## Contents

  * 1 Website
  * 2 Version History
  * 3 Prerequisites
    * 3.1 Mac OS X
    * 3.2 Linux
  * 4 Installation
  * 5 Usage
  * 6 Software
    * 6.1 Inputs
    * 6.2 Overlays
    * 6.3 Plots
    * 6.4 Word of Caution
  * 7 Requests and Disclaimer
  * 8 Main Window
  * 9 Pairwise Aminoacid Heatmap
  * 10 License



## Website

<http://emptyewer.github.io/CMPyMOL/>

## Version History

  * Added support for reading multi-model PDB files. (Multi-Model PDB file for NMR structure or Trajectory from MD simulations.)
  * Supports displaying variance of contact points from a series of contact maps generated from multi-model PDB files.
  * CMPyMOL stores the calculated contact maps, heat maps and contact density information in a local SQLite database for fast and easy subsequent access.
  * Cleaner GUI.
  * Parallelized the code for contact map calculation for multi-frame PDB files.



## Prerequisites

  * Python 2.7
  * Python module dependency 
    * [wxpython](http://www.wxpython.org)
    * [matplotlib](http://matplotlib.org)
    * [python imaging library (PIL)](http://www.pythonware.com/products/pil/)
    * [numpy](http://www.numpy.org)
  * PyMOL. (It is recommended that the user add the PyMOL installation directory to the $PATH environment variable.)
  * Stride secondary structure assignment tool. This program can be downloaded from <http://webclu.bio.wzw.tum.de/stride/> and compiled into a stand-alone executable. It is recommended that the Stride executable or its installation directory is added to the $PATH environment variable. NOTE: If this executable is not detected in the $PATH variable, the secondary structure calculation will be disabled in CMPyMOL.



### Mac OS X

  * Users can install the python libraries using "easy_install" or "pip". It is recommended that the user use Enthought Canopy python distribution and management package downloaded from <https://www.enthought.com/products/canopy/>. This package includes a robust python library management software and a python IDE.
  * PyMOL 1.5.x can be installed using MacPorts <http://www.macports.org>. NOTE: This automatically adds the executable into the $PATH environment variable.



### Linux

  * The python dependencies and PyMOL can be installed using apt-get (aptitude) or a similar package management system.



## Installation

There is no need for installation of the script. Optionally, a standalone executable can be complied using "pyinstaller" or "py2exe" or "py2app" package, depending on the users operating system. 

## Usage
    
    
    python /<path to CMPyMOL directory>/CMPyMOL.py
    

  * This command will automatically invoke the PyMOL executable and the user is led through the rest of the program with a series of pop-up windows.



## Software

Clicking (left) and dragging a selection of contact points on the displayed contact map will highlight the corresponding residues in the PyMOL window (as red and blue colored atoms in spheres representation). In addition, several structural/biochemical properties can be overlaid on top of the contact map. The contact-map data can also be plotted in other representations. The calculated contact-map, heat-map and contact density information is stored in a local SQL database. Any subsequent access of the same PDB with matching parameters will be read from the database for fast access. The code for calculating contact map for trajectory files is parallelized for efficiency. 

### Inputs

  * Single-frame PDB files (local)
  * Multi-frame PDB trajectory files (local)
  * Multi-frame trajectory files should have the following format.


    
    
    MODEL X
    .
    .
    .
    ATOM ...
    ATOM ...
    .
    .
    .
    ENDMDL
    

NOTE: The PDB can include REMARKS, CRYST and other standard PDB information entries. The MODEL line is essential for the software to work properly (ENDMDL is optional). 

### Overlays

  * Secondary structure of the protein is overlaid as translucent strips over the contact map. This button won't be active if secondary structure calculation program stride is not found in the system path ($PATH). (Button: Secondary Structure)
  * Contact points where a Charge-Charge interaction occurs are highlighted. (Button: Charged Interactions)
  * Residues that interact via hydrophobic interaction are highlighted. (Button: Hydrophobic Interactions)
  * Contact regions that have a B-factor that is higher than a certain cutoff are highlighted (Button: B-factor). The b-factor cutoff can be varied using a slider (Slider).
  * Highlights a contact point/region where the pair of selected residues are in contact (selected by checking the checkboxes). If only one aminoacid is selected from the list, interaction site of the selected aminoacid with another one of the same type is highlighted. (List of checkboxes for each aminoacid)



### Plots

  * Pairwise Heat Map - Plots a 20x20 matrix of pairwise aminoacid interaction count.
  * Contacts Histogram - Plots the number of contacts around a given residue. Selecting a particular bar highlights the corresponding residue in the PyMOL window.
  * Variance Contact Map - For Multi-frame PDB files (trajectory), this button toggles the displays the variance contact map starting from the initial frame until the current frame. This view can be used to identifying the dynamic regions in a protein.



### Word of Caution

When using a multi-frame PDB file, the contact-map for the next frame(s) are being pre-calculated in the background (depending on the number of free CPU cores available). Clicking on "Next Frame" in rapid succession may lead to undesired results and/or a crash. 

In the event of a crash, delete the database that is created in the working directory and relaunch the program. 

## Requests and Disclaimer

Users are welcome to send me an email to request the addition of a specific feature or to report a bug. 

## Main Window

[![CMPyMOL.PNG](/images/d/dd/CMPyMOL.PNG)](/index.php/File:CMPyMOL.PNG)

[](/index.php/File:CMPyMOL.PNG "Enlarge")

## Pairwise Aminoacid Heatmap

[![Heatmap-CMPyMOL.PNG](/images/3/39/Heatmap-CMPyMOL.PNG)](/index.php/File:Heatmap-CMPyMOL.PNG)

[](/index.php/File:Heatmap-CMPyMOL.PNG "Enlarge")

## License
    
    
    #  The MIT License (MIT)
    # =======================
    # 
    # The PyMOL Plugin source code in this file is copyrighted, but you are
    # free to use and copy it as long as you don't change or remove any of
    # the copyright notices.
    # 
    # -----------------------------------------------------------------------------------
    # CMPyMOL
    # Copyright (C) 2013 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy of this 
    # software and associated documentation files (the "Software"), to deal in the Software 
    # without restriction, including without limitation the rights to use, copy, modify, 
    # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
    # permit persons to whom the Software is furnished to do so, subject to the following 
    # conditions:
    #
    # The above copyright notice and this permission notice shall be included in all copies 
    # or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
    # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    # OTHER DEALINGS IN THE SOFTWARE.
    

Retrieved from "[https://pymolwiki.org/index.php?title=CMPyMOL&oldid=12359](https://pymolwiki.org/index.php?title=CMPyMOL&oldid=12359)"


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

## Contact map visualizer

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/contact_map_visualizer.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/contact_map_visualizer.py)  
Author(s)  | [Venkatramanan Krishnamani](/index.php/User:Venky "User:Venky")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**Enhanced version of this plugin is now available at[CMPyMOL](/index.php/CMPyMOL "CMPyMOL")**

The **contact map visualizer** plugin can link contact map images to the residues in PyMOL in a interactive way. Contact maps are pixel graphics where each protein residue corresponds to one line and one row of pixels. Thus for a 100 residue protein, such a image has 100x100 pixels. A common tool to generate such images is **g_mdmat** from the [gromacs](http://www.gromacs.org) package. 

## Contents

  * 1 Usage
  * 2 Required Dependencies
  * 3 Installation
  * 4 Generate Contact Map
  * 5 Screenshot
  * 6 Copyright



## Usage
    
    
    contact_map_visualizer [ image_file [, selection ]]
    

## Required Dependencies

  * [pygame](http://pygame.org/)
  * [Tkinter](http://wiki.python.org/moin/TkInter) (optional and usually included with PyMOL)
  * [PIL](http://www.pythonware.com/products/pil/) (optional, for automatically converting XPM images)



Example for installing all dependencies on a Ubuntu like system: 
    
    
    sudo apt-get install python-tk python-imaging python-pygame
    

## Installation

  * Navigate to plugins > install
  * Locate the downloaded **contact_map_visualizer.py** file in the dialogbox and select 'OK'
  * Quit and Restart 'pymol'



Alternative way: Just [run](/index.php/Run "Run") the script, it will provide a command but no menu plugin entry. 

## Generate Contact Map

[![](/images/e/ea/Contact-Map-of-a-Trajectory.png)](/index.php/File:Contact-Map-of-a-Trajectory.png)

[](/index.php/File:Contact-Map-of-a-Trajectory.png "Enlarge")

Mean contact map of a protein trajectory generated from g_mdmat tool in the Gromacs analysis package.

Use the command [g_mdmat](http://manual.gromacs.org/online/g_mdmat.html) from [Gromacs](http://www.gromacs.org) analysis package. A typical contact map looks like the figure on the right. 

To generate contact map of a single PDB. For example contact map for a PDB from RCSB, use the following command 
    
    
    g_mdmat -f <protein.pdb> -s <protein.pdb> -mean contact-map.xpm
    

To generate a mean contact map form a protein trajectory 
    
    
    g_mdmat -f <trajectory.pdb> -s <starting-frame.pdb> -mean contact-map.xpm
    

For Gromacs 2021 To generate contact map of a single PDB. For example contact map for a PDB from RCSB, use the following command 
    
    
    gmx g_mdmat -f <protein.pdb> -s <protein.pdb> -mean contact-map.xpm
    

To generate a mean contact map form a protein trajectory 
    
    
    gmx g_mdmat -f <trajectory.pdb> -s <starting-frame.pdb> -mean contact-map.xpm
    

To convert XPM to PNG format 
    
    
    convert contact-map.xpm contact-map.png
    

## Screenshot

[![Screenshot.png](/images/e/e6/Screenshot.png)](/index.php/File:Screenshot.png)

[](/index.php/File:Screenshot.png "Enlarge")

## Copyright
    
    
    # Copyright Notice
    # ================
    # 
    # The PyMOL Plugin source code in this file is copyrighted, but you are
    # free to use and copy it as long as you don't change or remove any of
    # the copyright notices.
    # 
    # -----------------------------------------------------------------------------------
    # This PyMOL Plugin Contact Maps Visualizer is
    # Copyright (C) 2012 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
    # 
    #                        All Rights Reserved
    # 
    # Permission to use, copy, modify, distribute, and distribute modified
    # versions of this software and its documentation for any purpose and
    # without fee is hereby granted, provided that the above copyright
    # notice appear in all copies and that both the copyright notice and
    # this permission notice appear in supporting documentation, and that
    # the name(s) of the author(s) not be used in advertising or publicity
    # pertaining to distribution of the software without specific, written
    # prior permission.
    # 
    # THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
    # INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
    # NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
    # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
    # USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
    # OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    # PERFORMANCE OF THIS SOFTWARE.
    #
    

Retrieved from "[https://pymolwiki.org/index.php?title=Contact_map_visualizer&oldid=13459](https://pymolwiki.org/index.php?title=Contact_map_visualizer&oldid=13459)"


---

## Dehydron

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/dehydron.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/dehydron.py)  
Author(s)  | [Osvaldo Martin](/index.php/User:OsvaldoMartin "User:OsvaldoMartin")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 Installation
      * 1.1.1 Linux
      * 1.1.2 Windows
      * 1.1.3 Mac OsX
    * 1.2 Usage
  * 2 Methods
    * 2.1 Acknowledgement
    * 2.2 Change log
    * 2.3 References



# Introduction

A dehydron is a protein backbone hydrogen bond incompletely shielded from water attack. A desolvated hydrogen bond is energetically more favourable than one exposed to the solvent and hence dehydrons are _sticky_ , since they promote the removal of surrounding water through protein associations or ligand binding. 

Dehydrons are less conserved than other structural motifs, hence identification of dehydrons could help to increase specificity during the rational drug design process. Certain proteins are enriched in dehydrons such as membrane proteins, toxic proteins and proteins that have a strong tendency to aggregate. Dehydrons have been invoked to explain biological processes above the molecular level such as the dosage imbalance effect in duplicated genes and the high connectivity of the protein interactomes of higher organisms. 

A putative dehydron can be detected by counting the number of _wrappers_ that surround a hydrogen bond. A wrapper is defined as a carbon atom not bonded directly to an oxygen or nitrogen atom, i.e. a non-polar carbon atom. A dehydron is defined by the number of wrappers inside two overlapping spheres centred at the Cα carbon of the donor and acceptor residues. If the number of wrappers around an hydrogen bond is below a certain cut-off value that hydrogen bond is identified as a dehydron. 

  


## Installation

### Linux

This plugin is ready "out-of-box" for Linux users through the project [ Pymol-script-repo](/index.php/Git_intro "Git intro")

### Windows

This plugin is ready "out-of-box" for Windows users through the project [ Pymol-script-repo](/index.php/Git_intro "Git intro")

### Mac OsX

This plugin have not been tested on a Mac OsX machine, but it should work. 

## Usage

The plugin can be accessed using the following command: 
    
    
    wrappy [ selection [, angle_range [, max_distance [, desolv [, min_wrappers ]]]]]
    

[![](/images/c/cf/Dehydrons.png)](/index.php/File:Dehydrons.png)

[](/index.php/File:Dehydrons.png "Enlarge")

**Figure 0** : Dehydrons calculated and displayed in PyMOL.

  
Or using a graphical environment (see figure 0) 

  
There are five parameters the user can change: 

Two of them control the hydrogen bonds detection. 

  * Angle range: deviation in degrees from the optimal hydrogen bond angle (C=0 N-H).
  * Max distance: maximum donor-acceptor distance in angstroms.



Although no thoroughly optimized, the default values for the hydrogen bonds parameters were adjusted to get a close agreement between the hydrogen bond listed by wrappy and the ones listed by the [What-If](http://swift.cmbi.ru.nl/servers/html/index.html) server. 

Another two control the dehydron detection. 

  * Desolvatation sphere radius: this parameter controls the radius of the two spheres centred at the Cα carbon of the donor and acceptor residues. A dehydron is defined by the number of "wrappers" inside this two spheres.  

  * Min wrappers: a hydrogen bond surrounded with less "wrappers" than "min_wrappers" is a dehydron. Setting this parameter to a "high" value, something like 100, will return all main chain hydrogen bonds (according to the angle range and max distance parameters). The default value (19) was taken from a statistical analysis of ~7400 high-quality proteins form the Protein Data Bank (see below).
  * Max wrappers: a hydrogen bond surrounded with more "wrappers" than "max_wrappers" is a buried hydrogen bond. The default value (35) was taken from a statistical analysis of ~7400 high-quality proteins form the Protein Data Bank (see below).



A wrapper is defined as a carbon atom not bonded directly to an oxygen or nitrogen atom, i.e. a non-polar carbon atom. The plug-in count as wrappers any non-polar carbon from any protein chain, organic ligand or other type of molecule, if the atoms belong "selection" (see below). This means that if you have, for example, a dimeric protein you will probably get different results for the dimer and for the isolated monomers. Instead, if you upload two (or more) different files the results will be independent because the plug-in does not count atoms from other objects 

  * Selection: This parameter allows the user to select which part of system is used to calculate dehydrons. This parameter is useful, for example, to calculate dehydronds for different objects independently or to easily calculate dehydrons with and without an organic ligand. The default selection is "all". Compute dehydrons for an specific selection is equivalent to delete all but the selected atoms and then compute dehydrons.



# Methods

[![](/images/c/c0/Wrappers_histogram.png)](/index.php/File:Wrappers_histogram.png)

[](/index.php/File:Wrappers_histogram.png "Enlarge")

**Figure 1** :Histogram of the number of wrappers per hydrogen bond (blue bars), the distribution of wrappers approximate a Gaussian distribution (red line). Parameters for the Gaussian fit (mean and standard deviation)are show in the grey box

An analysis of 7476 high quality X-ray proteins was performed in order to estimate the number of wrappers that should be used as a cut-off to determine whether to call an hydrogen bond a dehydron (i.e. the _min wrappers_ parameter). Although this value have been already estimated in the literature; differences in the algorithms used before and the ones used by wrappy could lead to differences in the exact value of the cut-off and hence this parameter was re-estimated in order to obtain reliable calculation of the dehydrons. To compute the numbers of wrappers the following parameters were used _angle range_ = 40°, _Max distance_ = 3.5 Å and _desolvatation sphere_ = 6.5 Å. A non-redundant set of 7476 proteins were obtained from the Protein Data Bank. Each protein in this set conforms with the following criteria: Resolution < 2.0 Å, R-factors <= 0.25, not containing DNA and/or RNA molecules. Additionally, proteins with a sequence identity of 30% were removed. The frequency of wrappers in this set of proteins approximate a Gaussian distribution (see figure 1) with a mean of ~27 and a standard deviation of ~8, hence an hydrogen bond with less than 19 wrappers is defined as a dehydron (27-8 = 19). In the same fashion and over-wrapped hydrogen bond is and hydrogen bond with more than 35 wrappers (27+8). 

Using the same set of 7476 proteins it was obtained that on average a high quality and globular protein should have 0.62 hydrogen bonds per residue (with a standard deviation of 0.06) and 17 wrappers per residue (with a standard deviation of 2). Hydrogen bonds per residue and wrappers per residue could be used as indicators of the global protein structure quality. Wrappy reports such indicators as z-score, i.e. the number of standard deviations an observation is above or below the expected mean value. 

  


## Acknowledgement

The H-bond detection code is based on the list_mc_hbonds.py script from Robert L. Campbell <http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/>

  


## Change log

  * 2012-01-14 (Version 1.0) 
    1. First public version was released and put under version control. In the project, [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro).


  * 2012-01-28 (Version 1.1) 
    1. Minor changes in the code most of them not visible for the end-user.


  * 2012-02-28 (Version 1.5) 
    1. The code was cleaned (e.g. remove global variables and other ugly stuff)
    2. The code was made available as a PyMOL command
    3. Better support for multiple objects



All features in this version and most of the code was provided by [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3"), thanks Thomas! :-) 

  * 2012-03-14 (Version 1.6) 
    1. Representation is not changed to "cartoon" after each calculation
    2. Total control over the selection from which dehydrons are calculated


  * 2013-03-26 (Version 1.7) 
    1. Wrappers were not correctly counted for structures with hydrogen atoms.



Thanks Shafqat Rasool for reporting the bug. 

  * 2013-07-19 (Version 2.0) 
    1. The plug-in was renamed to wrappy.
    2. Wrappy reports, now, the hydrogen bonds per residue and the wrappers per residue as z-scores. The values of the mean and standard deviation, necessary to compute zscores, where taken from the analysis of ~7400 high quality X-ray proteins from the PDB.


  * 2014-03-22 (Version 2.1) 
    1. Wrappy reports all the hydrogens bonds. Dehydrons are displayed using red dashes, average wrapped hydrogens bonds are yellow, and over-wrapped hydrogens bonds are green.



## References

Citation for Dehydrons:  
De Simone, A., Dodson, G. G., Verma, C. S., Zagari, A., and Fraternali, F. (2005). Prion and water: tight and dynamical hydration sites have a key role in structural stability. Proceedings of the National Academy of Sciences of the United States of America, 102(21), 75357540. 

Fernández, A. and Berry, R. S. (2003). Proteins with h-bond packing defects are highly interactive with lipid bilayers: Implications for amyloidogenesis. Proceedings of the National Academy of Sciences, 100(5), 2391–2396. 

Fernández, A. and Crespo, A. (2008). Protein wrapping: a molecular marker for association, aggregation and drug design. Chemical Society Reviews, 37(11), 2373. 

Fernández, A. and Lynch, M. (2011). Non-adaptive origins of interactome complexity. Nature, 474(7352), 502–505. 

Fernández, A., Rogale, K., Scott, R., and Scheraga, H. A. (2004a). Inhibitor design by wrapping packing defects in HIV-1 proteins. Proceedings of the National Academy of Sciences of the United States of America, 101(32), 1164011645. 

Fernández, A., Scott, R., and Berry, R. S. (2004b). The nonconserved wrapping of conserved protein folds reveals a trend toward increasing connectivity in proteomic networks. Proceedings of the National Academy of Sciences of the United States of America, 101(9), 28232827. 

Kabsch, W. and Sander, C. (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637. PMID: 6667333. 

Liang, H., Plazonic, K. R., Chen, J., Li, W.-H., and Fernández, A. (2008). Protein under-wrapping causes dosage sensitivity and decreases gene duplicability. PLoS Genetics, 4(1), e11. 

Retrieved from "[https://pymolwiki.org/index.php?title=Dehydron&oldid=12748](https://pymolwiki.org/index.php?title=Dehydron&oldid=12748)"


---

## DivScore

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A quote from the page: 

_The DivScore package is intended to be used with AMBER and DivCon for preprocessing biological structures (primarily proteins and small molecules) quick geometry optimization, charge calculation, creating job scripts and streamlining workflow within the PyMOL environment. This package was developed to aid my research and hence works within the confines of assumptions that were necessary for streamlining my workflow. However for those who have access to the AMBER and DivCon, it could potentially be useful since it combines the power of PyMOL with that of AMBER for molecular simulation and DivCon for Quantum Mechanics calculations on large biomolecular systems. An example of the capability of DivScore is that, it can perform a quick minimization of a ligand in the active site of a protein based on the AMBER force field and upload the resulting coordinates into PyMOL for visualization. More examples of things one could do will follow in the examples sections._

Web link: [DivScore](http://shoichetlab.compbio.ucsf.edu/~raha/research/divscore.html)

Retrieved from "[https://pymolwiki.org/index.php?title=DivScore&oldid=8490](https://pymolwiki.org/index.php?title=DivScore&oldid=8490)"


---

## Dockingpie

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Description
  * 2 Requirements
  * 3 DockingPie Download
  * 4 DockingPie Install
  * 5 How to use DockingPie
  * 6 Configuration
  * 7 How to Cite DockingPie
  * 8 References



## Description

**DockingPie** : a Molecular Docking plugin for PyMOL [1]. [Tutorial:Plugins:Threads](/index.php?title=Tutorial:Plugins:Threads&action=edit&redlink=1 "Tutorial:Plugins:Threads \(page does not exist\)") DockingPie implements: **Smina** , **Autodock Vina** , **RxDock** and **ADFR** [2-5]. 

Providing an easy interface to four docking programs, DockingPie is particularly suited as a platform to carry out consensus docking and scoring analyses. 

[![](/images/7/79/Dockingpie.png)](/index.php/File:Dockingpie.png)

[](/index.php/File:Dockingpie.png "Enlarge")

DockingPie

## Requirements

Minimal requirement: a recent version of PyMOL installed on your computer. 

DockingPie is compatible with incentive PyMOL builds distributed by Schrodinger (required PyMOL version >= 2.3.4) and open source builds (required PyMOL version >= 2.3.0). 

DockingPie is distributed freely to the public and it has been tested and runs on Windows, macOS and Linux versions of PyMOL. 

(Some incompatibilities may arise with the usage of PyMOL version 2.5.x if ‘undo’ function is enabled, which in PyMOL 2.5.2 still shows some shortcomings. Therefore, when the plugin is opened, the ‘undo’ function is automatically disabled and it is strongly suggested to keep it disabled when using the plugin.) 

## DockingPie Download

Download: <https://github.com/paiardin/DockingPie/archive/refs/heads/main.zip>

## DockingPie Install

DockingPie is installed, as any other PyMOL [1] plugin, via the PyMOL plugin manager: 

First download the latest version of the plugin ZIP file [here](https://github.com/paiardin/DockingPie/archive/refs/heads/main.zip)

Launch PyMOL and use the Plugin → Plugin Manager command from the main menu of PyMOL. The plugin manager window of PyMOL will open. 

Click on Install New Plugin and press the Choose File… button. Select the DockingPie ZIP file which you have downloaded before. You will be asked to give the path of the directory in which to install the plugin files. Just select the default option if you are unsure about what to do (the location of the plugin files does not make any difference when running the plugin). 

## How to use DockingPie

A detailed explanation on how to use DockingPie, some videos and tutorials can be found at the following links. 

  * User's Guide: [download from here](https://github.com/paiardin/DockingPie/releases/download/versioning/DockingPie_User_Guide.pdf)


  * DockingPie GitHub Wiki: [Home](https://github.com/paiardin/DockingPie/wiki) [Tutorials](https://github.com/paiardin/DockingPie/wiki/Tutorials)


  * Website: [Structural Bioinformatics Group at Sapienza](http://schubert.bio.uniroma1.it)



## Configuration

DockingPie, at the current and first release, integrates four different docking programs: RxDock [2], Vina [3], Smina [4] and ADFR [5]; several chemo-informatics python modules (i.e. AutoDockTools [6], Openbabel [7], sPyRMSD [8]) and other external tools like sdsorter [9]. The CONFIGURATION tab provides an easy way for the installation of the needed tools from within the plugin in two steps, as reported next. 

  * Configure external tools: CONFIGURATION tab → Configure → Start Download → Finish Download


  * Install external tools: If the needed tools are not currently installed on the user’s machine, the *Install* button is enabled and it can be used to install the external components.



## How to Cite DockingPie

If you find DockingPie useful, please cite: 

Serena Rosignoli and Alessandro Paiardini, DockingPie: a consensus docking plugin for PyMOL, Bioinformatics, 2022, btac452, [DOI](https://doi.org/10.1093/bioinformatics/btac452)

## References

[1] DeLano, WL. (2002). “The PyMOL Molecular Graphics System on World Wide Web.” CCP4 Newsletter On Protein Crystallography 

[2] Ruiz-Carmona, S., Alvarez-Garcia, D., Foloppe, N., Garmendia-Doval, A. B., Juhos, S., Schmidtke, P., Barril, X., Hubbard, R. E., & Morley, S. D. (2014). rDock: a fast, versatile and open source program for docking ligands to proteins and nucleic acids. PLoS computational biology, 10(4), e1003571. <https://doi.org/10.1371/journal.pcbi.1003571>

[3] O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461 

[4] Koes, David Ryan et al. “Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise.”Journal of chemical information and modeling vol. 53,8 (2013): 1893-904. doi:10.1021/ci300604z 

[5] Ravindranath PA, Forli S, Goodsell DS, Olson AJ, Sanner MF. AutoDockFR: Advances in Protein-Ligand Docking with Explicitly Specified Binding Site Flexibility. PLoS Comput Biol. 2015;11(12):e1004586. Published 2015 Dec 2. doi:10.1371/journal.pcbi.1004586 

[6] Morris, G. M., Huey, R., Lindstrom, W., Sanner, M. F., Belew, R. K., Goodsell, D. S., & Olson, A. J. (2009). AutoDock4 and AutoDockTools4: Automated docking with selective receptor flexibility. Journal of Computational Chemistry, 30(16), 2785–2791. 

[7] O'Boyle, N.M., Banck, M., James, C.A. et al. Open Babel: An open chemical toolbox. J Cheminform 3, 33 (2011). <https://doi.org/10.1186/1758-2946-3-33>

[8] Meli, R., Biggin, P.C. spyrmsd: symmetry-corrected RMSD calculations in Python. J Cheminform 12, 49 (2020). <https://doi.org/10.1186/s13321-020-00455-2>

[9] <https://sourceforge.net/projects/sdsorter/>

[10] Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from <https://docs.anaconda.com/>

Retrieved from "[https://pymolwiki.org/index.php?title=Dockingpie&oldid=13558](https://pymolwiki.org/index.php?title=Dockingpie&oldid=13558)"


---

## DSSP Stride

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/dssp_stride.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/dssp_stride.py)  
Author(s)  | [Hongbo Zhu](/index.php/User:Hongbo_zhu "User:Hongbo zhu")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[DSSP](http://swift.cmbi.ru.nl/gv/dssp/) and [Stride](http://webclu.bio.wzw.tum.de/stride/) are popular tools for assigning secondary structures for proteins. DSSP & Stride Plugin for PyMOL provides a graphical user interface for coloring proteins according to their secondary structures assigned by DSSP or Stride. Most recent DSSP & Stride Plugin code can be obtained from [this link](https://github.com/hongbo-zhu-cn/Pymol-script-repo/blob/master/plugins/dssp_stride.py). 

## News

  * 2014_01_24_Update: The plugin is able to cope with empty chain name now.


  * 2012_10_22_Update: The most recent version of DSSP (v2.0.4) does not consider residues after the first [`TER` record](http://www.wwpdb.org/documentation/format33/sect9.html#TER) from the same chain. The plugin has been updated to cope with this new feature.



## External links

[![](/images/9/9b/Demo_DSSP_plugin_1pyg.png)](/index.php/File:Demo_DSSP_plugin_1pyg.png)

[](/index.php/File:Demo_DSSP_plugin_1pyg.png "Enlarge")

Demonstration of DSSP plugin (pdb:1pyg).

  * Most recent code on github: [DSSP & Stride Plugin for PyMOL @ GitHub](https://github.com/hongbo-zhu-cn/Pymol-script-repo/blob/master/plugins/dssp_stride.py)


  * Original site (not updated anymore): [DSSP & Stride Plugin for PyMOL](http://www.biotec.tu-dresden.de/~hongboz/dssp_pymol/dssp_pymol.html)



## See Also

  * [dssp](/index.php/Dssp "Dssp") (psico)
  * [dss](/index.php/Dss "Dss")



Retrieved from "[https://pymolwiki.org/index.php?title=DSSP_Stride&oldid=12654](https://pymolwiki.org/index.php?title=DSSP_Stride&oldid=12654)"


---

## EMovie

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/emovie.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/emovie.py)  
Author(s)  | [Israel Structural Proteomics Center](http://www.weizmann.ac.il/ISPC/home.html)  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 eMovie
    * 1.1 Program Features
    * 1.2 Screenshots
      * 1.2.1 eMovie in use
      * 1.2.2 eMovie menu bar
    * 1.3 Availability
    * 1.4 Authors
    * 1.5 Reference



# eMovie

eMovie is a plug-in for PyMOL that makes the creation of molecular movies both easy and intuitive via a breakthrough storyboard interface, similar in nature to what is used in the creation of traditional movies. The eMovie homepage is accessible at [www.weizmann.ac.il/ISPC/eMovie.html](http://www.weizmann.ac.il/ISPC/eMovie.html). 

## Program Features

eMovie is currently the most user-friendly way for users to create movies in PyMOL (even inexperienced users). Users interact with a user-friendly eMovie GUI that does not require typing commands into PyMOL. 

Modular actions such as zooms, rotations, fadings, and morphs (morphs require incentive PyMOL) can be inserted to any frame in the movie and the actions comprising the movie can be reviewed in list-format by viewing the eMovie storyboard. The storyboard also allows for deletion and reinsertion of actions. 

Movies can be saved, loaded, and exported as a series of image files to be later merged into a traditional movie format such as .mov using an external program like GraphicConverter. 

## Screenshots

### eMovie in use

[![EMovie in use.jpg](/images/3/38/EMovie_in_use.jpg)](/index.php/File:EMovie_in_use.jpg)

### eMovie menu bar

[![EMovie menubar.jpg](/images/3/3d/EMovie_menubar.jpg)](/index.php/File:EMovie_menubar.jpg)

## Availability

eMovie is freely available for download at the [eMovie homepage](http://www.weizmann.ac.il/ISPC/eMovie.html). 

## Authors

eMovie was created at the [Israel Structural Proteomics Center](http://www.weizmann.ac.il/ISPC/home.html) (ISPC) at the Weizmann Institute of Science. 

## Reference

Hodis, E., Schreiber, G., Rother, K., Sussman, J.L., eMovie: a storyboard-based tool for making molecular movies, Trends in Biochemical Sciences 32, 199-204 (2007). 

Retrieved from "[https://pymolwiki.org/index.php?title=EMovie&oldid=10436](https://pymolwiki.org/index.php?title=EMovie&oldid=10436)"


---

## EZ-Viz

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## EZ-Viz

EZ-Viz was developed as an assistance tool for the difficult to understand user interface of PyMOL. The standard interface is very complex takes a long time to learn, especially for students unfamiliar with the chemistry and busy researchers with little time to learn the commands necessary to use the PyMOL software to its full extent. EZ-Viz allows for users to navigate through PyMOL using a series of clickable buttons and tabs to complete all functions of the PyMol software without ever having to type a single command.  
  
For more information and downloads, please visit our site [ProMOL.org](http://www.rit.edu/~ez-viz/). 

Retrieved from "[https://pymolwiki.org/index.php?title=EZ-Viz&oldid=4020](https://pymolwiki.org/index.php?title=EZ-Viz&oldid=4020)"


---

## Geo Measures Plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/lkagami/geo_measures_pymol/archive/master.zip>  
Author(s)  | Luciano Porto Kagami, Gustavo Machado das Neves, Luís Fernando Saraiva Macedo Timmers, Rafael Andrade Caceres and Vera Lucia Eifler-Lima   
License  | GNU General Public License v3.0   
  
The Geo Measures Plugin can carry out geometric analysis on protein structures.In addition, it makes other trajectory analyzes such as Probability Density Function (PDF), Root Mean Square Deviation (RMSD), Radius of Gyration (RG), Free Energy Landscape (FEL), Principal Component Analysis (PCA), Ramachandran map, Root Mean Square Fluctuation (RMSF), Define Secondary Structure of Proteins (DSSP), and Modevectors. Altogether there are 14 tools, which can be easily used by the graphical interface. 

[![](/images/6/6f/Geo_measures.png)](/index.php/File:Geo_measures.png) [](/index.php/File:Geo_measures.png "Enlarge")Geo-Measures Graphical User Interface | [![](/images/7/79/Figure_4.png)](/index.php/File:Figure_4.png) [](/index.php/File:Figure_4.png "Enlarge")Free Energy Landscape | [![](/images/a/a6/Ramachandran.png)](/index.php/File:Ramachandran.png) [](/index.php/File:Ramachandran.png "Enlarge")Ramachandran Plot  
---|---|---  
[![](/images/a/a6/PDF.png)](/index.php/File:PDF.png) [](/index.php/File:PDF.png "Enlarge")Probability Density Function | [![](/images/c/c4/Modevector.png)](/index.php/File:Modevector.png) [](/index.php/File:Modevector.png "Enlarge")Modevectors | [![](/images/f/fa/Pincer_Angle.png)](/index.php/File:Pincer_Angle.png) [](/index.php/File:Pincer_Angle.png "Enlarge")Pincer Angle  
  
## Video Tutorial

* * *

<https://youtu.be/YKmIcU5weI0>

* * *

## Installation

* * *

Requirements: 

PyMOL version higher than 2.0. 

It is highly recommended to install PyMOL using the Ananconda3 package [[1]](https://www.anaconda.com/distribution/). 

You need to install python packages. 

To do this use the command: 

pip install pandas matplotlib scipy mdtraj sklearn 

* * *

  
Install [geo_measures_pymol-master.zip](https://github.com/lkagami/geo_measures_pymol/archive/master.zip) using the PyMOL [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")

## References

* * *

Kagami LP, das Neves GM, Timmers LFSM, Caceres RA,Eifler-Lima VL, Geo-Measures: A PyMOL plugin for protein structure ensembles analysis,Computational Biology and Chemistry(2020) doi.org/10.1016/j.compbiolchem.2020.107322 

  


* * *

**Universidade Federal do Rio Grande do Sul - Laboratório de Síntese Orgânica Medicinal - LaSOM**

Retrieved from "[https://pymolwiki.org/index.php?title=Geo_Measures_Plugin&oldid=13343](https://pymolwiki.org/index.php?title=Geo_Measures_Plugin&oldid=13343)"


---

## GROMACS Plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [pymol_plugin_dynamics.py](https://raw.githubusercontent.com/makson96/Dynamics/master/pymol_plugin_dynamics.py)  
Author(s)  | [Laboratory of Biomolecular Systems Simulations](http://biotech.ug.edu.pl/dzialalnosc_naukowa/zespoly_badawcze/pracownia_symulacji_ukladow_biomolekularnych)  
License  | [GPLv3](https://www.gnu.org/licenses/gpl.html)  
  
[Dynamics PyMOL Plugin](https://github.com/makson96/Dynamics) is plugin for PyMOL, which add molecular dynamics simulation feature. It is meant to be easy to use. The plugin uses GROMACS tools as a back-end. Project is developed as an open source and as such create full open source stack together with PyMOL and GROMACS. Software works on Linux, MacOS X and Windows/Cygwin. 

## Contents

  * 1 Website
  * 2 Features
  * 3 Installation
    * 3.1 Ubuntu
    * 3.2 Gentoo
    * 3.3 Other Linux distributions and macOS
    * 3.4 Windows/Cygwin
    * 3.5 Latest Snapshots
  * 4 References
  * 5 License



## Website

<https://github.com/makson96/Dynamics>

## Features

[![Menu1.png](/images/5/55/Menu1.png)](/index.php/File:Menu1.png)

  * Easy to use GUI, to take advantage of complex software GROMACS.
  * Work directly on molecules loaded to PyMOL.
  * Display results of calculations directly in PyMOL.
  * Set restraints, choose water models, force fields and many more.
  * Save your work and finish calculations later or on the other machine.
  * Do your work free of charge and without any restrictions. Feel free to modify any component of the stack.
  * Minimum dependencies. Plugin use the same graphic libraries as PyMOL, so working PyMOL and GROMACS installations are enough to make plugin work.



Note that exact instruction of how to use program is in the [software manual](https://github.com/makson96/Dynamics/blob/master/manual.odt), which is available together with the plugin. 

## Installation

There are at least few ways to install the plugin. In this section we will describe the most common. 

### Ubuntu

The easiest way to install the plugin is to use PPA repositories: 
    
    
    sudo add-apt-repository ppa:makson96/dynamics
    sudo apt-get update
    sudo apt-get install pymol-plugin-dynamics
    

All required libraries will be downloaded as dependencies. Supported versions of Ubuntu are latest stable and latest LTS. 

### Gentoo

The plugin is also present in [Gentoo ebuild system](http://packages.gentoo.org/package/sci-chemistry/pymol-plugins-dynamics?arches=all). To make the installation use (as root): 
    
    
    emerge --ask sci-chemistry/pymol-plugins-dynamics
    

### Other Linux distributions and macOS

The second way for all other platforms, is to download latest release of the plugin directly from its GitHub webpage:  
<https://github.com/makson96/Dynamics/releases>   
Then run PyMOL as a root. On the top Menu choose **Plugin- >Manage Plugins->Install...** Then choose **pymol_plugin_dynamics.py** file. After installation restart PyMOL with normal user privilege. On macOS you will need to use [X11 Hybrid](https://pymolwiki.org/index.php/MAC_Install#X11_Hybrid) version of the PyMOL. 

Required dependencies: 

  * PyMOL
  * GROMACS
  * ProDy (optional)



In order to take advantage of latest features you will need to have [ProDy](http://www.csb.pitt.edu/prody/) library installed. 

### Windows/Cygwin

  * Download and install the latest version of Cygwin including appropriate code development packages.
  * Download, compile, and install the latest version of GROMACS 2016.3 under Cygwin. Follow the standard compilation, installation and testing instructions to build, compile and install GROMACS 2016.3 under Cygwin. For example:


    
    
    % cd gromacs-2016.3
    % mkdir build
    % cmake .. -DBUILD_SHARED_LIBS=off -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/cygdrive/c/<YOUR_GROMACS_ROOT>
    % make >& make.log
    % make install >& make.install.log
    

  * Test GROMACS from Cywin command line.
  * Download and install 64 bit PyMOL 1.8.2.0, which comes with Python 2.7.9.
  * Test PyMOL.
  * Update PATH environment variable on Windows to include the following bin directories corresponding to GROMACS and Cygwin:


    
    
    C:\<YOUR_GROMACS_ROOT>\bin
    C:\<YOUR_CYGWIN_PATH>\bin
    

  * Install the GROMACS plugin following standard PyMOL instructions. Download latest plugin version from link: <https://github.com/makson96/Dynamics/archive/master.zip> and unpack it. Run PyMOL. On the top Menu choose Plugin->Manage Plugins->Install... Then choose pymol_plugin_dynamics.py file. After installation restart PyMOL



### Latest Snapshots

You can download latest master snapshot by command: 
    
    
    git clone <git://github.com/makson96/Dynamics.git>
    

Then run PyMOL as a root. On the top Menu choose **Plugin- >Manage Plugins->Install...** Then choose **pymol_plugin_dynamics.py** file. After installation restart PyMOL with normal user privilege. 

Required dependencies: 

  * PyMOL
  * GROMACS
  * ProDy (optional)
  * PLUMED (optional, compiled into GROMACS)



In order to take advantage of latest features you will need to have [ProDy](http://www.csb.pitt.edu/prody/) library installed and [PLUMED](http://www.plumed.org/) compiled into GROMACS. 

## References

[Molecular Dynamics Simulation by GROMACS Using GUI Plugin for PyMOL Tomasz Makarewicz and Rajmund Kaźmierkiewicz Journal of Chemical Information and Modeling 2013 53 (5), 1229-1234](http://pubs.acs.org/doi/abs/10.1021/ci400071x)  
[Improvements in GROMACS plugin for PyMOL including implicit solvent simulations and displaying results of PCA analysis Tomasz Makarewicz, Rajmund Kaźmierkiewicz Journal of Molecular Modeling May 2016, 22:109](http://link.springer.com/article/10.1007/s00894-016-2982-4)

## License

This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3". Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.  
Contributors:  
\- Tomasz Makarewicz   
\- Ajit B. Datta   
\- Sara Boch Kminikowska   
\- Manish Sud   


Retrieved from "[https://pymolwiki.org/index.php?title=GROMACS_Plugin&oldid=13266](https://pymolwiki.org/index.php?title=GROMACS_Plugin&oldid=13266)"


---

## Gyration tensor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <http://emptyewer.github.io/Gyration_tensor>  
Author(s)  | [Venkatramanan Krishnamani](/index.php/User:Venky "User:Venky")  
License  | \-   
  
## Contents

  * 1 Introduction
    * 1.1 How to get it
  * 2 Required Dependencies on Linux and Mac OS X
    * 2.1 Installing on Linux
    * 2.2 Usage
    * 2.3 Screenshot
    * 2.4 Requests for features are welcome
  * 3 Copyright



## Introduction

This program calculates and visualizes the gration tensors of a protein (based on geometry) for each chain in the PDB file. 

* * *

#### How to get it
    
    
    git clone <https://github.com/emptyewer/Gyration_tensor.git>
    

* * *

## Required Dependencies on Linux and Mac OS X

  * [Tkinter](http://wiki.python.org/moin/TkInter) and
  * [Numpy](http://numpy.scipy.org/)


    
    
    sudo apt-get install python-tk python-numpy
    

* * *

### Installing on Linux

  * Navigate to plugins > install
  * Locate the downloaded gyration_tensor.py file in the dialogbox and select 'OK'
  * Quit and Restart 'pymol'



* * *

### Usage

  * Launch PyMOL
  * Navigate to Plugins > Contact Map Visualizer in PyMOL
  * When the first dialog box pops up, select PDB file and select 'OK'.
  * The plugin calculates and outputs the Magnitude and Direction of the tensors in the debug window and draws the corresponding tensors as a cgo object in the PyMOL opengl window.



* * *

### Screenshot

[![Gyration tensor.png](/images/9/90/Gyration_tensor.png)](/index.php/File:Gyration_tensor.png)

[](/index.php/File:Gyration_tensor.png "Enlarge")

* * *

### Requests for features are welcome

Please see your user account [User:Venky](/index.php/User:Venky "User:Venky"). 

* * *

## Copyright
    
    
    # Copyright Notice
    # ================
    # 
    # The PyMOL Plugin source code in this file is copyrighted, but you are
    # free to use and copy it as long as you don't change or remove any of
    # the copyright notices.
    # 
    # -----------------------------------------------------------------------------------
    # This PyMOL Plugin Gyration Tensor is
    # Copyright (C) 2012 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
    # 
    #                        All Rights Reserved
    # 
    # Permission to use, copy, modify, distribute, and distribute modified
    # versions of this software and its documentation for any purpose and
    # without fee is hereby granted, provided that the above copyright
    # notice appear in all copies and that both the copyright notice and
    # this permission notice appear in supporting documentation, and that
    # the name(s) of the author(s) not be used in advertising or publicity
    # pertaining to distribution of the software without specific, written
    # prior permission.
    # 
    # THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
    # INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
    # NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
    # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
    # USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
    # OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    # PERFORMANCE OF THIS SOFTWARE.
    #
    

* * *

Retrieved from "[https://pymolwiki.org/index.php?title=Gyration_tensor&oldid=11227](https://pymolwiki.org/index.php?title=Gyration_tensor&oldid=11227)"


---

## Helicity check

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**helicity_check** shows the evolution of O—N distances over an amino acid sequence 

See for further info: 

["Models for the 3(10)-helix/coil, pi-helix/coil, and alpha-helix/3(10)-helix/coil transitions in isolated peptides"](http://www.proteinscience.org/cgi/content/abstract/5/8/1687). _Protein Sci_ ROHL and DOIG 5 (8) 1687. 

Uses: 

  * in the pymol console:


    
    
     >run pymol_helicity_check.py
       ----> select some consecutive amino acids
              - this is nicely done with the Display->Sequence tool
     >helicity_check()
    

  * installing helicity_check


    
    
     copy pymol_helicity_check.py in $PYMOL_INSTALL_DIR/modules/pmg_tk/startup
     launch Pymol: you now have a new option in the Plugin menu
    
    

helicity_check uses [gnuplot](http://www.gnuplot.info) to display its results. As a consequence gnuplot needs to be installed. 

This plugin was tested on linux only, it might need some modifications to run on other OSes. (hints: launching gnuplot and path to dumpfile) 

  

    
    
    # pymol_helicity_check.py
    # Copyright (c) 2006-2007 Julien Lefeuvre <lefeuvrejulien@yahoo.fr>
    #
    
    """
    Pymol plugin for checking helicity type
    
    helicity_check() takes as input a selection ('sele' by default)
    of at least 5 amino acids and computes the distances between
    O(i) - N(i+3)
    O(i) - N(i+4)
    O(i) - N(i+5)
    See for further info:
    Protein Sci ROHL and DOIG 5 (8) 1687
    'Models for the 3(10)-helix/coil, pi-helix/coil,
    and alpha-helix/3(10)-helix/coil transitions in isolated peptides.'
    
    uses:
    *in the pymol console:
      >run pymol_helicity_check.py
        ----> select some consecutive amino acids
               - this is nicely done with the Display->Sequence tool
      >helicity_check()
    *installing helicity_check
      copy pymol_helicity_check.py in $PYMOL_INSTALL_DIR/modules/pmg_tk/startup
      launch Pymol: you now have a new option in the Plugin menu
    
    helicity_check uses gnuplot (http://www.gnuplot.info) to display its results
    As a consequence gnuplot needs to be installed.
    
    This plugin was tested on linux only, it my need some modifications to run on
    other OSes (hints: launching gnuplot and path to dumpfile)
    """
    
    __author__ =    "Julien Lefeuvre <lefeuvrejulien@yahoo.fr>"
    __version__ =   "1.0"
    __date__ =      "2007-04-02"
    __copyright__ = "Copyright (c) 2007 %s. All rights reserved." % __author__
    __licence__ =   "BSD"
    
    from pymol import cmd
    from math import sqrt
    import sys
    import os
    import subprocess
    import time
    
    def __init__(self):
        """init function in order to have a nice menu option in Pymol"""
        self.menuBar.addmenuitem('Plugin', 'command', 'Helicity Check',
                 label='Helicity Check', command = lambda: helicity_check())
    
    
    class Residue(object):
    
        def __init__(self):
            self.name=None
            self.index=None
            self.Ocoord=None
            self.Ncoord=None
    
    
    def calc_distON(Ocoord,Ncoord):
        """return the distance between 2 atoms given their coordinates"""
        sum = 0
        for o, n in zip(Ocoord, Ncoord):
            sum += (o - n)**2
        return sqrt(sum)
    
    
    def helicity_check(selection='sele'):
        """calcultate distance O[res i]-N[res i+3]
                               O[res i]-N[res i+4]
                               O[res i]-N[res i+5]
        """
        seq_model = cmd.get_model(selection) #get info from selection
        res_lim = seq_model.get_residues()
    
        if len(res_lim)<5:
            sys.stderr.write("\nPlease select at least 5 residues\n")
            return
    
        atom_list = seq_model.atom
        res_data=[]
    
        for start,end in res_lim:   #extract the data we are interested in
            res=Residue()
            for atom in atom_list[start:end]:
                if atom.name == 'N':
                    res.name = atom.resn
                    res.index = atom.resi
                    res.Ncoord = atom.coord
                elif atom.name == 'O':
                    res.Ocoord = atom.coord
            if res.Ocoord and res.Ncoord and res.name and res.index:
                res_data.append(res)
            else:
                sys.stderr.write("\nPlease select complete protein residues\n")
                return
    
        res_list = [int(res.index) for res in res_data]
    
        if res_list != range(res_list[0], res_list[-1]+1):
            sys.stderr.write("\nPlease select a unbrocken residue sequence\n")
            return
    
        distON3 = []
        distON4 = []
        distON5 = []
        distONs = [distON3, distON4, distON5]
    
        for i,res in enumerate(res_data[:-5]): #distances calculations
            resis = res_data[i+3:i+6]
            for resi, distONi in zip(resis, distONs):
                distONi.append(calc_distON(res.Ocoord, resi.Ncoord))
    
        dump = os.tmpnam()+'.dat'
        dumpfile = file(dump, 'w')
    
        sys.stdout.write('\n#Distances O(i)---N(i+n)\n'
               '#ResNum , d(O(i)-N(i+3)) , d(O(i)-N(i+4)) , d(O(i)-N(i+4))\n')
        for i, d3, d4, d5 in zip(res_list, distON3, distON4, distON5):
            #writing console output
            sys.stdout.write(
                  '  %i ,      %f ,       %f ,       %f \n'%(i, d3, d4, d5))
            #writing data to a dump file for use by gnuplot
            dumpfile.write(
                  '  %i       %f        %f        %f \n'%(i, d3, d4, d5))
        dumpfile.flush()
    
        #launch a gnuplot window to show the distances
        gnuplotcmd = subprocess.Popen(['/usr/bin/gnuplot'], shell=True,
                                   stdin=subprocess.PIPE)
        gnuplotcmd.stdin.write('set autoscale\n')
        gnuplotcmd.stdin.write("plot "
             "'%s' using 1:2 title 'd(O(i)-N(i+3))' with lines, "
             "'%s' using 1:3 title 'd(O(i)-N(i+4))' with lines, "
             "'%s' using 1:4 title 'd(O(i)-N(i+5))' with lines\n'"
                              % (dump, dump, dump))
        time.sleep(3)
        dumpfile.close()
        os.remove(dump)
    

## Download

[Helicity_check-1.0.tar.zip](/images/f/fb/Helicity_check-1.0.tar.zip "Helicity check-1.0.tar.zip")

Retrieved from "[https://pymolwiki.org/index.php?title=Helicity_check&oldid=12200](https://pymolwiki.org/index.php?title=Helicity_check&oldid=12200)"


---

## Isoslider

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [scripts/isoslider.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/isoslider.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Isocontour slider ("density slider") plugin. It opens a dialog with isolevel sliders for each [isomesh](/index.php/Isomesh "Isomesh") and [isosurface](/index.php/Isosurface "Isosurface") object that is currently loaded in PyMOL. 

If the mouse is above the slider, the mouse wheel is bound to increasing/decreasing the contour level. 

## See Also

  * [isolevel](/index.php/Isolevel "Isolevel")
  * Similar plugin: <http://www.ebi.ac.uk/~gareth/pymol/pymol.shtml#density_slider>



Retrieved from "[https://pymolwiki.org/index.php?title=Isoslider&oldid=13913](https://pymolwiki.org/index.php?title=Isoslider&oldid=13913)"


---

## Lisica

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/lisica.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/lisica.py)  
Author(s)  | Janez Konc   
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Important
  * 2 Description
  * 3 LiSiCA Plugin Features
  * 4 Requirements
  * 5 Installation
  * 6 Usage
    * 6.1 Inputs Tab
    * 6.2 Load Project Tab
    * 6.3 Output Tab
    * 6.4 About Tab
  * 7 Reference
  * 8 Licensing



## Important

Currently, the most recent version can be obtained at the [Insilab web page](http://insilab.org/lisica-plugin). 

## Description

LiSiCA is a software for 2D and 3D ligand based virtual screening. It uses a fast maximum clique algorithm to find two- and three- dimensional similarities between reference compound and a database of target compounds in Mol2 format. The similarities are expressed using Tanimoto coefficients and the target compounds are ranked accordingly. LiSiCA is developed and maintained at **National Institute of Chemistry, Slovenia**. 

LiSiCA Plugin is **free for Academic (NON-COMMERCIAL) use**. However, for **COMMERCIAL use the potential users have to write to[Janez Konc](mailto:konc@cmm.ki.si)**. 

  


## LiSiCA Plugin Features

  1. Graphical User Interface which facilitate the choice of parameters for LiSiCA.
  2. Results displayed according to the ranking in similarities as measured by Tanimoto coefficients.
  3. Structural similarity between molecules visualized using PyMOL viewer.



  


## Requirements

The plugin should work on Windows and all variety of Unix like systems. Most of the tests were performed on Windows 7, Ubuntu 14.04 and Linux Mint. To work properly, the plugin requires Python with Tkinter and PyMOL. The plugin was mainly tested with PyMOL 1.7.x. 

  


## Installation

This plugin is ready "out-of-box" in most Linux and Windows operating systems and is available on GitHub in the Pymol-script-repo repository. 

The lisica.py script initializes the installation of the plugin. The script downloads and installs all the required files. On successful installation, a directory named .lisicagui is downloaded and saved in the home directory. This folder contains the executables, log files, icon files, python modules etc. 

Once the files are properly installed, LiSiCA plugin version will be ready to use. 

**Important note** : If installation fails with "pmg_tk" error or similar, one or both of the tkinter and ttk libraries are not installed. To install them open terminal and type as root "yum install tkinter python-pip" and "pip install pyttk" (in CentOS) or "apt-get install tkinter python-pip" and then "pip install pyttk" (in Ubuntu). 

**Important note** : For complete installation of the plugin, a stable INTERNET connection is required. 

## Usage

The plugin window has four tabs. 

#### Inputs Tab

Both the reference and the target compound files need to be in the Tripos mol2 format. The reference file should contain only one (reference) compound; the target file may contain many compounds. If the reference file contains more than one compound, the first molecule is used as a reference. The molecules in the target file having the same name are considered as different conformers of the same molecule. By default only the best-scoring (by Tanimoto coefficient) conformer will be shown in the final output for the 3D screening option. 

In the mol2 input files, the name of the molecule (ZINC ID) must be specified under the line for @<TRIPOS>MOLECULE tag. For example, 
    
    
      @<TRIPOS>MOLECULE
      **ZINC73655097**
      46    48     0     0     0
      SMALL
      USER_CHARGES
      
      @<TRIPOS>ATOM
      ...
    

LiSiCA checks similarities based on the mol2 atom types. Hence **SYBYL** atom types have to be specified in the mol2 files. For example, under the @<TRIPOS>ATOM tag, each atom specification should include the SYBYL atom type, as in: 
    
    
      1 C1         -0.0647    1.4496   -0.0592   **C.3**       1 <0>        -0.167
    

In the above line from a mol2 file, in bold, in column 5, the SYBYL atom type is specified. 

  
In 2D screening, the option **Maximum Allowed Shortest path Difference** corresponds to the maximum allowed difference in shortest-path length between atoms of the two compared product graph vertices. Lesser values correspond to a more rigorous screening. By default this value is unit bond. 

In 3D screening, the option **Maximum Allowed Spatial Distance Difference** corresponds to the maximum allowed difference in distances between atoms of the two compared product graph vertices. Lesser values correspond to a more rigorous screening. By default this value is 1 Å. The **Number of Conformations** option corresponds to the maximum number of outputted files of one molecule in different conformations and is to be used only for 3D screening. 

  
[![InputTab2DLiSiCALinuxv1.0.0.png](/images/a/ad/InputTab2DLiSiCALinuxv1.0.0.png)](/index.php/File:InputTab2DLiSiCALinuxv1.0.0.png) [![InputTab3DLiSiCALinuxv1.0.0.png](/images/5/5e/InputTab3DLiSiCALinuxv1.0.0.png)](/index.php/File:InputTab3DLiSiCALinuxv1.0.0.png)

  
According to the value of **Number of highest ranked molecules to be written to the output (say W)** , LiSiCA will create mol2 files of that many (W) highest scoring target molecules with a comment section at the end of the file where the matching atom pairs are displayed. This value is by default 100. The resulting mol2 files are written into a time-stamp directory in the folder specified in **Save results in:**. By default this folder is the user's home directory. Also, a text file named lisica_results.txt with the target molecules with (in the descending order of) Tanimoto coefficients is written to the time-stamp folder. 

The **Number of CPU cores to be used** allows selection of CPU threads used for LiSiCA. By default, it tries to detect the number of CPUs available. The **Consider Hydrogen** options lets the user to choose if the hydrogen atoms are to be considered for the calculation of the similarity using the maximum clique algorithm. By default, hydrogen atoms are not considered in finding the largest substructure common to the reference and target molecules, so as to obtain faster results. 

#### Load Project Tab

The plugin also has a feature to load saved results. On the **Load Project** tab, the user can choose the directory with the saved results (mol2 files of each target and the reference) and the lisica_results.txt file. When the load button is clicked, the results will be loaded onto the output tab and the PyMOL Viewer window. 

#### Output Tab

In the output tab, there are two listboxes: 

  * One contains ZINC ID and Tanimoto Coefficients of target molecules in the decreasing order of the Tanimoto coefficient values.



Any single target molecule can be selected on this listbox using a mouse click or using up/down arrow keys. The selected target molecule is displayed with the reference molecule on the PyMOL viewer window. 

  * Depending on the target molecule chosen on the first listbox, the corresponding atoms from reference and the target molecules are displayed on the other listbox.



Any single pair of corresponding atoms can be selected on this listbox using a mouse click or using up/down arrow keys. The selected pair is highlighted on the PyMOL viewer window. 

  
For **2D Screening** , the two molecules (the reference and the selected target) are visualized side by side on the PyMOL viewer screen. 

[![2D LiSiCA output- Note the reference and target molecules are aligned side by side](/images/a/a7/2DResultOutputtabPyMOLViewerLinux1.png)](/index.php/File:2DResultOutputtabPyMOLViewerLinux1.png "2D LiSiCA output- Note the reference and target molecules are aligned side by side")

For **3D Screening** , the 3D structures of the two molecules (the reference and the selected target) are superimposed on one another to visualize the similarity on the PyMOL viewer screen. 

[![3D LiSiCA output-Note that the reference and target molecules are superimposed](/images/6/64/3DResultsOutputtabPyMOLViewer1.png)](/index.php/File:3DResultsOutputtabPyMOLViewer1.png "3D LiSiCA output-Note that the reference and target molecules are superimposed")

#### About Tab

The users can get information on new updates if available on the About tab. 

  


## Reference

**If you are using LiSiCA in your work, please cite:**

S. Lesnik, T. Stular, B. Brus, D. Knez, S. Gobec, D. Janezic, J. Konc, LiSiCA: A Software for Ligand-Based Virtual Screening and Its Application for the Discovery of Butyrylcholinesterase Inhibitors, J. Chem. Inf. Model., 2015, 55, 1521–1528. 

[Read the article](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00136)

[PDF Link](http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00136)

  


## Licensing

LiSiCA software is copyrighted by : 

    **Janez Konc**
    National Institute of Chemistry
    Laboratory for Molecular Modeling
    Hajdrihova 19
    SI-1000 Ljubljana
    Slovenia.

The terms stated in the following agreement apply to all files associated with the software unless explicitly disclaimed in individual files. 
    
    
       **LiSiCA SOFTWARE LICENSE AGREEMENT**
    
     1. Grant Of Limited License; Software Use Restrictions. The programs received by you will be used only for NON COMMERCIAL purposes.
        This license is issued to you as an individual.
        
        For COMMERCIAL use of the software, please contact Janez Konc for details about commercial usage license agreements.
        For any question regarding license agreements, please contact:
        
        Janez Konc
        National Institute of Chemistry
        Laboratory for Molecular Modeling
        Hajdrihova 19
        SI-1000 Ljubljana
        Slovenia.
        
     2. COMMERCIAL USAGE is defined as revenues generating activities. These
        include using this software for consulting activities and selling
        applications built on top of, or using this software. Scientific 
        research in an academic environment and teaching are considered 
        NON COMMERCIAL.
    
     3. Copying Restrictions. You will not sell or otherwise distribute commercially 
        these programs or derivatives to any other party, whether with or without 
        consideration.
    
     4. Ownership of Software. You will not obtain, and will not attempt to 
        obtain copyright coverage thereon without the express purpose written 
        consent of Janez Konc.
    
     5. IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
        FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
        ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
        DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
        POSSIBILITY OF SUCH DAMAGE.
    
     6. THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
        INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
        IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
        NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
        MODIFICATIONS.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Lisica&oldid=12429](https://pymolwiki.org/index.php?title=Lisica&oldid=12429)"


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

## MSMS

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/msms.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/msms.py)  
Author(s)  | [Hongbo Zhu](/index.php/User:Hongbo_zhu "User:Hongbo zhu")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[MSMS](http://mgltools.scripps.edu/packages/MSMS/) is an excellent tool for computing protein solvent excluded surface (SES). MSMS Plugin for PyMOL provides a graphical user interface for running MSMS and displaying its results in PyMOL. See the original site (not updated anymore) for more details: [MSMS Plugin for PyMOL](http://www.biotec.tu-dresden.de/~hongboz/msms_pymol/msms_pymol.html)

Most recent MSMS Plugin code can be obtained from [this link](https://github.com/hongbo-zhu-cn/Pymol-script-repo/blob/master/plugins/msms.py). 

## External links

[![](/images/7/70/Example.png)](/index.php/File:Example.png)

[](/index.php/File:Example.png "Enlarge")

Demonstration of MSMS plugin.

  * [Molecular Surfaces](http://en.wikipedia.org/wiki/Molecular_surface)
  * Most recent code on github: [MSMS Plugin for PyMOL @ GitHub](https://github.com/hongbo-zhu-cn/Pymol-script-repo/blob/master/plugins/msms.py)
  * Original site (not updated anymore): [MSMS Plugin for PyMOL](http://www.biotec.tu-dresden.de/~hongboz/msms_pymol/msms_pymol.html)



## See Also

  * [msms_surface](/index.php/Msms_surface "Msms surface") (psico)



Retrieved from "[https://pymolwiki.org/index.php?title=MSMS&oldid=12653](https://pymolwiki.org/index.php?title=MSMS&oldid=12653)"


---

## NsSNP Loader

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Description
  * 2 Comments
  * 3 Download
  * 4 Author



## Description

This plugin allows the user to enter a SNP (rs) or a Swiss-Prot id and load the structural model from the large scale non-synonymous SNP database [LS-SNP](http://alto.compbio.ucsf.edu/LS-SNP/) if one exists. If more than one nsSNP is present for the chosen protein, all nsSNPs are listed and the user is able to choose which nsSNPs they wish to view. Once selected, the modeled structure is loaded, each residue is mutated to the nsSNP to allow easy comparison to the base structure. 

## Comments

As this plugin requires an internet connection it also requires the [Proxy Config](/index.php/Proxy_Config "Proxy Config") plugin. 

## Download

Download by right clicking the link and choosing 'Save Link As...'.  
[Download here](http://vbc.med.monash.edu.au/~fauxn/pymol_plugin/nsSNP_loader.py)  
Then run PyMol and under the 'Plugin' menu option choose 'Install Plugin...'. 

## Author

[Noel Faux](/index.php/User:Fauxn "User:Fauxn")

Retrieved from "[https://pymolwiki.org/index.php?title=NsSNP_Loader&oldid=3457](https://pymolwiki.org/index.php?title=NsSNP_Loader&oldid=3457)"


---

## Optimize

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/optimize.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/optimize.py)  
Author(s)  | [Osvaldo Martin](/index.php/User:OsvaldoMartin "User:OsvaldoMartin")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Installation
  * 3 Usage
    * 3.1 Local optimization
      * 3.1.1 Usage
      * 3.1.2 PyMOL API
    * 3.2 Global optimization
  * 4 Change log
  * 5 Citation
  * 6 Known Issues
  * 7 See Also



# Introduction

_Optimize_ provides a PyMOL graphical interface to some of the molecular mechanics features available in [openbabel](http://openbabel.org), allowing the user to optimize (minimize) the energy of any molecule uploaded on PyMOL. 

# Installation

The plugin can be downloaded through the project [ Pymol-script-repo](/index.php/Git_intro "Git intro"). 

_Optimize_ needs OpenBabel (and OpenBabel Python bindings) to be installed on your computer (see instructions [here](http://openbabel.org/wiki/Get_Open_Babel)). 

**Note to Windows users:** If you are using a 64-bit Windows, please install OpenBabel 64-bit from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#openbabel)

# Usage

The plugin can be accessed using a graphic user interface (see figure 0) or from the PyMOL`s terminal. There are 5 types of optimization routines available now, 2 local and 3 global optimization routines. 

## Local optimization

Local optimization can be done using the `minimize` command. 

### Usage
    
    
    minimize [selection [, forcefield [, method [, nsteps [, conv [, cutoff [, cut_vdw [, cut_elec]]]]]]]]
    

[![](/images/d/de/Optimize.png)](/index.php/File:Optimize.png)

[](/index.php/File:Optimize.png "Enlarge")

**Figure 0** : Optimize plugin GUI

  * **selection** : The name of the object to be minimized. From the command line, the default value is 'all'. When using the GUI, the default value is the first uploaded object.
  * **forcefield** : The forcefield used to compute the Internal Energy. Available options are GAFF, MMFF94s (default), MMFF94, UFF and Ghemical.
  * **method** : The method used to find the local minimum. Available options are "conjugate gradients" (default) and "steepest descent".
  * **nsteps** : Number of iteration steps during the minimization (default = 500).
  * **conv** : Criteria used to judge minimization convergence (default = 0.0001).
  * **cutoff** : Control if cut-off are used or not to compute non-bonded interactions, possible values are True or False (default).
  * **cut_vdw** : If cutoff is True, then this parameter sets the distance (in Angstroms) beyond which two atoms do not interact through Van der Waals forces (default = 6.0).
  * **cut_elec** : If cutoff is True, then this parameter sets the distance (in Angstroms) beyond which two atoms do not interact through electrostatic forces (default = 8.0).



  


### PyMOL API
    
    
    cmd.minimize(string selection="all", string forcefield="MMFF94s", string method="conjugate gradients", 
    int nsteps=500, float conv=0.0001, bool cutoff=False, float cut_vdw=6.0, float cut_elec=8.0)
    

## Global optimization

Global optimization can be done using the conf_search command from the PyMOL`s terminal: 
    
    
    conf_search [selection string [, forcefield string [, method string [, nsteps int [, conformers int [, lowest_confor int]]]]]]
    

  
Where: 

  * selection: The name of teh object that is going to be minimized. The default value is 'all'. Using the GUI the default value is the first uploaded object.
  * forcefield: Choose the forcefield used to compute the Internal Energy, options available are GAFF, MMFF94s, MMFF94, UFF and Ghemical.
  * method: Choose the method used to find the global minimum. The methods available are:



Systematic: Systematically iterate through all possible conformers according to Open Babel’s torsion library. This approach guarantee to find the local minimum (according to the forcefield in use). This approach scales to the power of N, where N is the number of rotatable bonds, hence it is only applicable to molecules with very few rotatable bonds. 

Random: Conformations are generated by randomly choosing from the allowed torsion angles. 

Weighted: This method uses an iterative procedure to find a global minimum. As with the Random method, it randomly choses from the allowed torsion angles but the choice is re-weighted based on the energy of the generated conformer. For molecules with to many rotatable bonds, that are not suitable for for the _systematic_ , this method is generally the best option. 

  * nsteps: Number of iteration steps during the minimization.
  * conformers: Total number of conformers to be analysed. This option is not available when using the _systematic_ method because all possible conformers are analysed.
  * lowest_conf: This options sets how many of the low-energy conformers are retrieved as the result of a conformational search. Conformers are ordered from low to high energy. This option is not available when using the _systematic_ method because this method return only the lowest energy conformer.



In general, it is a good idea to minimize the initial conformation before doing a conformational search. 

# Change log

  * 2013-10-06 (Version 0.1) 
    1. First version. (More features coming soon!)
  * 2013-10-24 (Version 0.2) 
    1. now openbabel add the hydrogen and not PyMOL.
  * 2014-01-27 (Version 0.6) 
    1. now it is possible to perform global optimization.
  * 2014-08-22 (Version 0.8) 
    1. Bug fixed. In version 0.6 the 'selection' and 'steps' parameters for the local minimization were incorrectly taken from the gobal optimization tab



# Citation

If you find optimize useful please consider citing this work Noel M. O'Boyle , Michael Banck , Craig A. James , Chris Morley , Tim Vandermeersch and Geoffrey R. Hutchison. "Open Babel: An open chemical toolbox." Journal of Cheminformatics 2011, 3:33. <http://www.jcheminf.com/content/3/1/33>

# Known Issues

  * Resizing the dialog on Windows with PyMOL 2.x (Qt interface) leads to a crash.



# See Also

  * [minimize_ob](/index.php/Minimize_ob "Minimize ob")



Retrieved from "[https://pymolwiki.org/index.php?title=Optimize&oldid=12907](https://pymolwiki.org/index.php?title=Optimize&oldid=12907)"


---

## DYNMAP

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## DYNMAP

DYN-MAP is a PyMol v0.99 plug-in that generates maps of functional groups present in a protein and visualize them using dynamic parameters. Why DYNMAP could be useful: 

  * Visualize distribution of functional groups in a protein, their solvent accessibility and their spatial fluctuation
  * Compare functional groups (dynmaps) among homologous proteins or mutants
  * Run automatically GROMACS molecular dynamics



## Availability

This plugin can be downloaded on this [website](http://www.giacomobastianelli.com/work.html)

## Author

[Giacomo Bastianelli](/index.php?title=User:Ghepardolesto&action=edit&redlink=1 "User:Ghepardolesto \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=DYNMAP&oldid=4076](https://pymolwiki.org/index.php?title=DYNMAP&oldid=4076)"


---

## PDB Index Search

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Description
  * 2 Availability
  * 3 Author
    * 3.1 See Also



## Description

The [Protein Data Bank (PDB)](http://www.rcsb.org/pdb/) is a repository of large molecule structures (proteins and nucleic acids). The website offers also a variety of tools and resources for studying these structures. The structures are available in PDB and mmCIF format, that can be read by PyMOL. In order to exploit these resources, some plugins have been developed. 

  
This plugin has been developed for PyMol v0.97. It works also with never version of PyMOL. This piece of software has been tested successfully on Linux, Mac OS X and Microsoft Windows. This plugin permits you to: 

  * search for any text into the PDB Index File ("Human Hemoglobin", i.e.)
  * retrieve the results into a spreadsheet. By double clicking on the row, the structure of the entry is downloaded.



## Availability

This plugin is Free Software and can be downloaded on this [website](http://pansanel.adlp.org/index.php?page=chemistry&sub_page=software). 

## Author

[Jerome Pansanel](/index.php/User:Pansanel "User:Pansanel")

### See Also

[fetch](/index.php/Fetch "Fetch")

Retrieved from "[https://pymolwiki.org/index.php?title=PDB_Index_Search&oldid=12119](https://pymolwiki.org/index.php?title=PDB_Index_Search&oldid=12119)"


---

## PDB Loader Service

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Description

Originally bundled with PyMol 0.96, the PDB Loader Service plugin dynamically downloads a requested PDB from [RCSB](http://www.rcsb.org/pdb/) and displays it in PyMol. 

## Comments

Although the plugin is completely functional, I have written newer versions which never touch the filesystem and as a result are much cleaner. See the [Tutorial on writing plugins](/index.php/Plugins_Tutorial "Plugins Tutorial") for one such example. 

## Author

[Charlie](/index.php/User:Cmoad "User:Cmoad") Moad 

Retrieved from "[https://pymolwiki.org/index.php?title=PDB_Loader_Service&oldid=12124](https://pymolwiki.org/index.php?title=PDB_Loader_Service&oldid=12124)"


---

## PDB plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 PDB PyMOL plugin
    * 1.1 Installation
      * 1.1.1 PyMOL 1.6
      * 1.1.2 PyMOL 1.7 or later
    * 1.2 PDB Analysis - All
    * 1.3 PDB Analysis - Molecules
    * 1.4 PDB Analysis - Domains
    * 1.5 PDB Analysis - Validation
    * 1.6 PDB Analysis - Assemblies
    * 1.7 View assemblies for your own mmCIF file
    * 1.8 Acknowledgements



## PDB PyMOL plugin

PDBe’s PyMOL plugin provides an easy way to visualise PDB data and annotations in PyMOL. This includes displaying individual molecules, Pfam, SCOP and CATH domains and viewing geometric validation for PDB entries based on PDBe's API. 

### Installation

The PDB PyMOL plugin works on PyMOL 1.6 or later 

#### PyMOL 1.6

Download the plugin from the following URL: 

<http://www.ebi.ac.uk/pdbe/pdb-component-library/pymol_plugin/PDB_plugin.py>

#### PyMOL 1.7 or later

Open the “plugin manager” Enter the following url in the “Install from PyMOL Wiki or URL 

<http://www.ebi.ac.uk/pdbe/pdb-component-library/pymol_plugin/PDB_plugin.py>

This will add the following items to the Plugins menu: 

PDB Analysis - All 

PDB Analysis - Molecules 

PDB Analysis - Domains 

PDB Analysis - Validation 

PDB Analysis - Assemblies 

  
Screenshot of the menu 

  
[![Screenshot of the PyMOL plugin menu with the PDB plugin installed](/images/e/e1/PDB_plugin.png)](/index.php/File:PDB_plugin.png "Screenshot of the PyMOL plugin menu with the PDB plugin installed")

  
To use each component click the menu item in the Plugins menu and then enter and PDB ID into the dialogue box. 

### PDB Analysis - All

This displays chemically distinct molecules, domains which make up a PDB entry. It also displays the assemblies within a PDB entry. 

Chemically distinct molecules which make up a PDB entry, including proteins, nucleic acids and ligands, are highlighted. Each molecule is given a unique colour and selected as a different object. Polymers (protein, DNA and RNA chains) are shown as cartoon or ribbon representation and ligands are shown as spheres. Pfam, SCOP and CATH domains within a PDB entry are highlighted. Each domain is coloured in a different colour and selected as a different object. This enables each domain to be turned on and off with PyMOL. Assemblies annotated within a PDB entry are highlighted. This requires pymol 1.76 or later. 

All PDB annotations show in PDB entry 3unb 

  
[![All PDB annotations shown on PDB entry 3unb](/images/2/26/Pdbe_3unb_all_annotation.png)](/index.php/File:Pdbe_3unb_all_annotation.png "All PDB annotations shown on PDB entry 3unb")

  


### PDB Analysis - Molecules

This highlights the chemically distinct molecules which make up a PDB entry, including proteins, nucleic acids and ligands. Each molecule is given a unique colour and selected as a different object. Polymers (protein, DNA and RNA chains) are shown as cartoon or ribbon representation and ligands are shown as spheres. 

In the below example the distinct molecules in PDB entry 3l2p are shown. 

  
[![Distinct molecules in PDB entry 3l2p](/images/b/b8/Pdbe_3l2p_entity.png)](/index.php/File:Pdbe_3l2p_entity.png "Distinct molecules in PDB entry 3l2p")

  


### PDB Analysis - Domains

This highlights the Pfam, Rfam, SCOP and CATH domains within a PDB entry. Each domain is coloured in a different colour and selected as a different object. This enables each domain to be turned on and off with PyMOL. The Domains Plugin also highlights chemically distinct molecules and domains are overlaid on these molecules. 

The CATH domains in PDB entry 3b43 are highlighted in the figure below. 

  
[![CATH domains shown on PDB entry 3b43](/images/0/0f/PDBe_3b43_domains_2.png)](/index.php/File:PDBe_3b43_domains_2.png "CATH domains shown on PDB entry 3b43")

### PDB Analysis - Validation

This overlays geometric validation on a PDB entry. Geometry outliers are coloured as per the PDB validation report: Residues with no geometric outliers are shown in green Residues with one geometric outliers are shown in yellow Residues with two geometric outliers are shown in orange Residues with three or more geometric outliers are shown in red. 

Geometric validation is shown on PDB entry 2gc2 

  
[![Geometric validation shown on PDB entry 2gc2](/images/7/78/Pdbe_2gc2_validation.png)](/index.php/File:Pdbe_2gc2_validation.png "Geometric validation shown on PDB entry 2gc2")

### PDB Analysis - Assemblies

This highlights the assemblies within a PDB entry. This requires pymol 1.76 or later. 

Assemblies are shown within PDB entry 5j96, with each assembly as an object in pymol and coloured in a distinct colour. 

  
[![Assemblies shown for PDB entry 5j96](/images/8/85/PDB_5j96_assemblies.png)](/index.php/File:PDB_5j96_assemblies.png "Assemblies shown for PDB entry 5j96")

  


### View assemblies for your own mmCIF file

To view annotated assemblies in your own mmCIF file - for example after PDB annotation - run pymol with the following command: 

pymol -r PDB_plugin.py -- mmCIF_file=YOUR_MMCIF_FILE 

where YOUR_MMCIF_FILE is the filename you want display the assemblies for. 

  


### Acknowledgements

This plugin was developed by PDBe [![PDBe logo](/images/7/75/Pdbe_logo.png)](http://pdbe.org "PDBe logo") who are a founding member of the wwPDB [![wwPDB logo](/images/d/de/Wwpdb-logo.png)](http://www.wwpdb.org/ "wwPDB logo")

Retrieved from "[https://pymolwiki.org/index.php?title=PDB_plugin&oldid=12772](https://pymolwiki.org/index.php?title=PDB_plugin&oldid=12772)"


---

## PICv

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [**PICv**](http://vidyaniranjan.co.in/PICv/PICv.py)  
Author(s)  | [Akshay Uttarkar](https://akshayuttarkar.wixsite.com/akshay), [Vasanth Kumar Desai](https://vpdesai2020.github.io/vasanth.github.io/), Namitha P, John Berrisford, [Sameer Velankar](https://www.ebi.ac.uk/about/people/sameer-velankar), [Vidya Niranjan](http://www.vidyaniranjan.co.in/)*   
License  | GNU Free Documentation License 1.2   
  
  


  


  


## Contents

  * 1 **About PICv**
  * 2 **Code**
  * 3 **Installation and Demo video**
  * 4 **Developed by**
  * 5 **Copyright**



## **About PICv**

[![PICv.png](/images/f/f2/PICv.png)](/index.php/File:PICv.png)

Protein interaction clustering and visualization is an pioneer attempt in understanding protein-protein interaction at a residue level. For any given protein the interaction is purely dependent on its charges and surface-structural modifications. The clustering of proteins based on there preferential amino acid interactions provides a biological insight on both the above mentioned aspects. The clusters such obtained can be used to infer the interaction behavior for a class or family of proteins. Such interpretation can be useful in understanding structural protein chemistry. The interactions also provide information on the crucial amino acids required for interactions to remain stable. This information can be used to design antibody or induce mutations to depreciate its functionality 

[![PICv home screen.PNG](/images/a/a7/PICv_home_screen.PNG)](/index.php/File:PICv_home_screen.PNG)

## **Code**

Download the plugin from the following URL: 

[**PICv Source Code**](http://vidyaniranjan.co.in/PICv/PICv.py)

## **Installation and Demo video**

Visit the below mentioned link for detailed installation followed by demo. 

[**PICv (Video)**](https://youtu.be/3bxpASy-FV0)

## **Developed by**

The plugin was developed by **Center of Excellence Computational Genomics** , R V College of Engineering, Bangalore, India in collaboration with **Protein Data Bank in Europe (PDBe), UK**

**For more details visit**

[**Center of Excellence Computational Genomics**](http://www.vidyaniranjan.co.in/)

## **Copyright**
    
    
    # Copyright Notice
    # ================
    # 
    # The PyMOL Plugin source code in this file is copyrighted, but you are
    # free to use and copy it as long as you don't change or remove any of
    # the copyright notices.
    # 
    # -----------------------------------------------------------------------------------
    # This PyMOL Plugin PICv is
    # Copyright (C) 2021 by Vidya Niranjan <vidya.n@rvce.edu.in>
    # 
    #                        All Rights Reserved
    # 
    # Permission to use, copy, modify, distribute, and distribute modified
    # versions of this software and its documentation for any purpose and
    # without fee is hereby granted, provided that the above copyright
    # notice appear in all copies and that both the copyright notice and
    # this permission notice appear in supporting documentation, and that
    # the name(s) of the author(s) not be used in advertising or publicity
    # pertaining to distribution of the software without specific, written
    # prior permission.
    # 
    # THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
    # INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
    # NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
    # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
    # USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
    # OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    # PERFORMANCE OF THIS SOFTWARE.
    #
    

Retrieved from "[https://pymolwiki.org/index.php?title=PICv&oldid=13433](https://pymolwiki.org/index.php?title=PICv&oldid=13433)"


---

## Plugin Manager

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Since version 1.5.0.5 PyMOL comes with a **Plugin Manager** , which can be used to load plugins such as those in the [PyMOL Script Repo](https://github.com/Pymol-Scripts/Pymol-script-repo#).  


## Contents

  * 1 Features
    * 1.1 Appending additional paths to the Plugin Manager
  * 2 Screenshots
  * 3 See Also



## Features

  * Install/Uninstall plugins
  * Disable plugins and load them on demand to optimize PyMOL's startup time
  * Configure the plugin search path



### Appending additional paths to the Plugin Manager

Should your scripts be located in several alternative locations, it is possible to append additional directories to the **Plugin Manager**.  


    Plugin > Plugin Manager > Settings > Add new directory...

## Screenshots

[![Plugin manager installed.png](/images/3/3d/Plugin_manager_installed.png)](/index.php/File:Plugin_manager_installed.png) [![Plugin manager install new.png](/images/1/1b/Plugin_manager_install_new.png)](/index.php/File:Plugin_manager_install_new.png) [![Plugin manager settings.png](/images/7/7b/Plugin_manager_settings.png)](/index.php/File:Plugin_manager_settings.png)

## See Also

  * [Pymol-script-repo](/index.php/Pymol-script-repo "Pymol-script-repo")



Retrieved from "[https://pymolwiki.org/index.php?title=Plugin_Manager&oldid=12823](https://pymolwiki.org/index.php?title=Plugin_Manager&oldid=12823)"


---

## PluginArchitecture

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes PyMOL's plugin architechture since version 1.5.0.6. 

## Contents

  * 1 API Entry Points
  * 2 __init_plugin__
    * 2.1 Example for a PyQt5 GUI (PyMOL 2.x)
    * 2.2 Example for a legacy Tkinter GUI (PyMOL 1.x)
  * 3 More than one file per plugin
  * 4 Config files
  * 5 PyMOL OS Fellowship Project 2011/2012
  * 6 See Also



## API Entry Points

A plugin is a Python module which uses PyMOL's API. The following entrypoints add functionality to PyMOL: 

  1. [`pymol.cmd.extend(func)`](/index.php/Extend "Extend") registers a function as a PyMOL command
  2. `pymol.plugins.addmenuitemqt(label, callback)` adds a plugin menu item, should be called inside `__init_plugin__`
  3. `MyPlugin.__init_plugin__(app)` is called during plugin initialization (if defined by the plugin)



## __init_plugin__

A plugin can implement an `__init_plugin__(app)` function (_previously`__init__`, deprecated_) which is called during plugin initialization. It is passed an object which is compatible with the legacy `PMGApp` instance and provides access to the Tkinter parent (`app.root`). 

### Example for a PyQt5 GUI (PyMOL 2.x)

_See also the[Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")_
    
    
    def __init_plugin__(app=None):
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('My Qt Plugin', myqtdialog)
    
    def myqtdialog():
        from pymol.Qt import QtWidgets
        QtWidgets.QMessageBox.information(None, 'My Plugin Title', 'Hello World')
    

### Example for a legacy Tkinter GUI (PyMOL 1.x)

It is important that all Tkinter objects are created with the `app.root` parent, otherwise legacy plugins won't work in PyMOL 2.0. 
    
    
    def __init_plugin__(app):
        app.menuBar.addmenuitem('Plugin', 'command',
            label='My Tk Plugin',
            command=lambda: mytkdialog(app.root))
    
    def mytkdialog(parent):
        try:
            import tkMessageBox  # Python 2
        except ImportError:
            import tkinter.messagebox as tkMessageBox  # Python 3
        tkMessageBox.showinfo(parent=parent, title='My Plugin Title', message='Hello World')
    

## More than one file per plugin

Plugins can be either a single Python file, or a directory with an `__init__.py` file. 

**Single file layout:**
    
    
    MyPlugin.py (defines "__init_plugin__(app=None)" function)
    

**Directory layout:**
    
    
    MyPlugin/
    ├── data
    │   ├── datasheet.txt
    │   └── image.png
    ├── submodule1.py
    ├── submodule2.py
    └── __init__.py (defines "__init_plugin__(app=None)" function)
    

They can be zipped (`.zip`) for distribution, the [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager") handles the unzipping during installation. 
    
    
    zip -r MyPlugin-1.0.zip MyPlugin
    

## Config files

PyMOL offers a basic API to store custom settings in `~/.pymolpluginsrc.py`. Only basic types (`str`, `int`, `float`, `list`, `tuple`, `dict`, etc.) are suppored (those who produce executable code with `repr()`). 

**Store:**
    
    
    pymol.plugins.pref_set(key, value)
    
    # needed if pymol.plugins.pref_get("instantsave") == False
    pymol.plugins.pref_save()
    

**Load:**
    
    
    value = pymol.plugins.pref_get(key, default=None)
    

**Example:**
    
    
    apbsbinary = pymol.plugins.pref_get("APBS_BINARY_LOCATION", "/usr/bin/apbs")
    

## PyMOL OS Fellowship Project 2011/2012

[Thomas Holder](/index.php/User:Speleo3 "User:Speleo3") was working on the new plugin system as part of his [2011-2012 Fellowship](http://pymol.org/fellowship). The improvements were incorporated into PyMOL 1.5.0.5. 

  * graphical [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * multiple files per plugin/module (see #More than one file per plugin)
  * user plugin directory
  * rarely used plugins can be disabled for better performance and enabled on demand
  * metadata support (version, author, description, citation, tags, ...)
  * install from online repositories
  * settings API for plugins (see #Config files)



## See Also

  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")
  * [Plugins](/index.php/Plugins "Plugins")



Retrieved from "[https://pymolwiki.org/index.php?title=PluginArchitecture&oldid=12829](https://pymolwiki.org/index.php?title=PluginArchitecture&oldid=12829)"


---

## Plugins

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Plugins are external modules that add functionality to PyMOL. 

Plugins can be single Python files (*.py) or directories with an `__init__.py` file (see [PluginArchitecture](/index.php/PluginArchitecture "PluginArchitecture")). 

Looking for available Plugins? Browse the [plugins category](/index.php/Category:Plugins "Category:Plugins"). 

## Contents

  * 1 Installing Plugins
  * 2 Writing Plugins
  * 3 Legacy Notes
    * 3.1 Version 3.0
    * 3.2 Before Version 2.0
    * 3.3 Before Version 1.5.0.6
    * 3.4 MacPyMOL
  * 4 See Also



## Installing Plugins

To install a plugin, open PyMOL and then go to 

    **Plugin > [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager") > Install new plugin**

You will have the option to install from a file on your computer (*.py, *.zip) or directly from the PyMOLWiki by providing the wiki page url. 

[![Plugin manager install new.png](/images/1/1b/Plugin_manager_install_new.png)](/index.php/File:Plugin_manager_install_new.png)

## Writing Plugins

See the [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial"). 

## Legacy Notes

### Version 3.0

Currently, support for plugins written in Tkinter is considered deprecated with full expectation of removal by PyMOL 4.0. Plugin developers should consider migrating legacy plugins to PyQt. 

### Before Version 2.0

Plugins written with PyQt5 do not run in PyMOL 1.x. 

PyMOL 2.x has a legacy layer for Tkinter to support old plugins, but the preferred toolkit for new plguins is PyQt5. 

### Before Version 1.5.0.6

Before PyMOL's Plugin Manager was added in version 1.5.0.6, plugins could not be installed per user, the plugin search path was a single directory inside the PyMOL installation: **$PYMOL_PATH/modules/pmg_tk/startup**

Plugins were limited to single Python files, directories were not supported yet. 

### MacPyMOL

MacPyMOL was the native macOS version of PyMOL before version 2.0. It supported Plugins only in the [X11 Hybrid](/index.php/MAC_Install#X11_Hybrid "MAC Install") mode. 

## See Also

  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")
  * [PluginArchitecture](/index.php/PluginArchitecture "PluginArchitecture")



Retrieved from "[https://pymolwiki.org/index.php?title=Plugins&oldid=13574](https://pymolwiki.org/index.php?title=Plugins&oldid=13574)"


---

## Plugins Tutorial

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This tutorial is about writing your own plugin for PyMOL 2.x. 

See the [plugins](/index.php/Plugins "Plugins") page for how to install and use exisiting plugins. 

## Contents

  * 1 Writing Plugins: Learn By Example
    * 1.1 Plugin Files
    * 1.2 Registering your Plugin
    * 1.3 Creating a GUI
    * 1.4 Filling the GUI with Widgets
    * 1.5 Make Buttons do something
    * 1.6 Deploy the final plugin
    * 1.7 Full Source
  * 2 Extending Plugins to the Command Line
  * 3 See Also



## Writing Plugins: Learn By Example

This tutorial shows how to write a PyMOL plugin with PyQt. **The full source of the demo plugin is[available on github](https://github.com/Pymol-Scripts/pymol2-demo-plugin)**. 

The demo plugin adds a dialog to render images at a custom resolution. 

### Plugin Files

We will create a plugin which consists of [multiple files inside a directory](/index.php/PluginArchitecture#More_than_one_file_per_plugin "PluginArchitecture"): 
    
    
    pymol2-demo-plugin/
    ├── demowidget.ui
    └── __init__.py (defines "__init_plugin__(app=None)" function)
    

### Registering your Plugin

First you must add your plugin to the _Plugins_ menu. This is done in the `__init_plugin__` function of your plugin. A callback (here: `run_plugin_gui`) is added. 
    
    
    # file pymol2-demo-plugin/__init__.py
    def __init_plugin__(app=None):
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Demo "Render" Plugin', run_plugin_gui)
    

_Legacy note_ : The `__init_plugin__` function takes one argument, a reference to the main Tk app to [support legacy plugins witten with Tkinter](/index.php/PluginArchitecture#init_plugin "PluginArchitecture") (unused with PyQt plugins). 

### Creating a GUI

The callback may do arbitrary stuff. Here we're going to create a dialog and show it to the user. 
    
    
    # global reference to avoid garbage collection of our dialog
    dialog = None
    
    def run_plugin_gui():
        # pymol.Qt provides the PyQt5 interface, but may support PyQt4 and/or PySide as well
        from pymol.Qt import QtWidgets
    
        global dialog
    
        if dialog is None:
            # create a new (empty) Window
            dialog = QtWidgets.QDialog()
    
            # TODO: FILL DIALOG WITH WIDGETS HERE
    
        dialog.show()
    

### Filling the GUI with Widgets

[![](/images/8/81/Designer-demo-plugin.png)](/index.php/File:Designer-demo-plugin.png)

[](/index.php/File:Designer-demo-plugin.png "Enlarge")

Screenshot of Qt Designer

The most convenient and maintainable way to create user interfaces with Qt is by using the [Qt Designer](http://doc.qt.io/qt-5/qtdesigner-manual.html). Follow these steps to get started: 

  1. Create a new "Widget" form (New Form > templates/forms > Widget > Create)
  2. Drag widgets like "Push Button" or "Line Edit" into your form
  3. Name your widgets ("objectName" in the "Property Editor"), this name will become a Python variable in your code
  4. Save as a UI file (`*.ui`) inside the `pymol2-demo-plugin` directory



PyMOL provides a utility function to load a UI file into a parent widget: `pymol.Qt.utils.loadUi(filename, widget)`
    
    
        # filename of our UI file
        uifile = os.path.join(os.path.dirname(__file__), 'demowidget.ui')
    
        # load the UI file into our dialog
        from pymol.Qt.utils import loadUi
        form = loadUi(uifile, dialog)
    

### Make Buttons do something

We need to connect the [signals](http://doc.qt.io/qt-5/signalsandslots.html) of those widgets we created with the Qt Designer, like our buttons' `clicked` signal. Our example has a "Ray" button (_objectName=button_ray_) that we'll connect to the `run()` function. 
    
    
        def run():
            # get form data
            height = form.input_height.value()
    
            # some debugging feedback
            print('User entered height', height)
    
            # TODO: DO SOMETHING WITH FORM DATA
    
        form.button_ray.clicked.connect(run)
    

### Deploy the final plugin

The `pymol2-demo-plugin` directory can be zipped for deployment (see [PluginArchitecture#More than one file per plugin](/index.php/PluginArchitecture#More_than_one_file_per_plugin "PluginArchitecture")). 

### Full Source

The full source of the demo plugin is [available on github](https://github.com/Pymol-Scripts/pymol2-demo-plugin)**.**

## Extending Plugins to the Command Line

While not really applicable to our simple "Render" plugin (it only wraps existing `ray` and `png` commands), most plugins should also provide their functionality as new PyMOL commands using [cmd.extend()](/index.php/Extend "Extend"). 

## See Also

  * [PluginArchitecture](/index.php/PluginArchitecture "PluginArchitecture")
  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins](/index.php/Plugins "Plugins")
  * [cmd.extend()](/index.php/Extend "Extend")



Retrieved from "[https://pymolwiki.org/index.php?title=Plugins_Tutorial&oldid=13173](https://pymolwiki.org/index.php?title=Plugins_Tutorial&oldid=13173)"


---

## Pmlbeta

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://gitlab.com/awacha/pmlbeta/-/raw/binaries/pmlbeta_latest.zip>  
Author(s)  | András Wacha   
License  | [BSD 3-clause](https://choosealicense.com/licenses/bsd-3-clause/)  
<https://pmlbeta.readthedocs.io>  
  
  
Pmlbeta is primarily a plugin for creating and manipulating molecular models of β-peptides, but works for α-peptides as well. In other words, the [Fab](/index.php/Fab "Fab") command on steroids. Its central piece is the betafab2 command, but a graphical user interface is also supplied. 

## Installation

Install the current version from [pmlbeta_latest.zip](https://gitlab.com/awacha/pmlbeta/-/raw/binaries/pmlbeta_latest.zip) using the PyMOL [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")

## On-line documentation

On-line documentation is available at <https://pmlbeta.readthedocs.io>

## See Also

  * [Fab](/index.php/Fab "Fab")



Retrieved from "[https://pymolwiki.org/index.php?title=Pmlbeta&oldid=13402](https://pymolwiki.org/index.php?title=Pmlbeta&oldid=13402)"


---

## PocketPicker

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PocketPicker -- binding site prediction 

Website: [PocketPicker](http://gecco.org.chemie.uni-frankfurt.de/pocketpicker/index.html)

Retrieved from "[https://pymolwiki.org/index.php?title=PocketPicker&oldid=8491](https://pymolwiki.org/index.php?title=PocketPicker&oldid=8491)"


---

## Polarpairs

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Find polar pairs with the [cmd.find_pairs](/index.php/Find_pairs "Find pairs") function and return them as a list of atom pairs. It finds (more or less) the same contacts like [cmd.distance(mode=2)](/index.php/Distance "Distance"). 

## Example
    
    
    pairs = polarpairs('chain A', 'chain B')
    for p in pairs:
        dist = cmd.get_distance('(%s`%s)' % p[0], '(%s`%s)' % p[1])
        print(p, 'Distance: %.2f' % (dist))
    

## The Script
    
    
    '''
    (c) 2011 Thomas Holder, MPI for Developmental Biology
    '''
    
    from pymol import cmd
    
    def polarpairs(sel1, sel2, cutoff=4.0, angle=63.0, name='', state=1, quiet=1):
        '''
    ARGUMENTS
    
        sel1, sel2 = string: atom selections
    
        cutoff = float: distance cutoff
    
        angle = float: h-bond angle cutoff in degrees. If angle="default", take
        "h_bond_max_angle" setting. If angle=0, do not detect h-bonding.
    
        name = string: If given, also create a distance object for visual representation
    
    SEE ALSO
    
        cmd.find_pairs, cmd.distance
        '''
        cutoff = float(cutoff)
        quiet = int(quiet)
        state = int(state)
        if angle == 'default':
            angle = cmd.get('h_bond_max_angle', cmd.get_object_list(sel1)[0])
        angle = float(angle)
        mode = 1 if angle > 0 else 0
        x = cmd.find_pairs('(%s) and donors' % sel1, '(%s) and acceptors' % sel2,
                state, state,
                cutoff=cutoff, mode=mode, angle=angle) + \
            cmd.find_pairs('(%s) and acceptors' % sel1, '(%s) and donors' % sel2,
                state, state,
                cutoff=cutoff, mode=mode, angle=angle)
        x = sorted(set(x))
        if not quiet:
            print('Settings: cutoff=%.1fangstrom angle=%.1fdegree' % (cutoff, angle))
            print('Found %d polar contacts' % (len(x)))
        if len(name) > 0:
            for p in x:
                cmd.distance(name, '(%s`%s)' % p[0], '(%s`%s)' % p[1])
        return x
    
    cmd.extend('polarpairs', polarpairs)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Polarpairs&oldid=12897](https://pymolwiki.org/index.php?title=Polarpairs&oldid=12897)"


---

## ProBiS H2O

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <http://insilab.org/probis-h2o/>  
Author(s)  | Marko Jukic   
License  |   
  
## Contents

  * 1 Update & Download
  * 2 Description
  * 3 ProBiS H2O Plugin Features
  * 4 Requirements & Dependencies
  * 5 Installation
  * 6 First Time Usage
  * 7 Usage
    * 7.1 Define System Tab
    * 7.2 Cluster Analysis Tab
    * 7.3 About/Help Tab
  * 8 Reference
  * 9 Licensing



## Update & Download

The most recent version can be obtained at the [Insilab web page](http://insilab.org/probis-h2o/). 

## Description

In order to help in identification of **conserved water sites** , we developed **ProBiS H2O** workflow that supports the complete process of data collection, extraction, identification of conserved water locations and result visualisation. **ProBiS H2O** workflow adopts the available experimental pdb data deposited in online databases or local user experimental data collections, collects similar macromolecular systems, performs local/specific binding site superimposition with the computational speed and accuracy of ProBiS algorithm, collects the experimental water location data reported in the parent macromolecule systems and transposes gathered data to the examined system, its specific chain, binding site or individual water. Water location data is then clustered to identify discrete spaces with high conservation of water molecules and visualised in the context of studied system. **ProBiS H2O** workflow is a robust, transparent and fast methodology with emphasis on experimental data that can identify trends in water localization intra/inter macromolecule or macromolecule-small molecule interfaces. For further info please refer to the relevant articles by the authors of ProBiS H2O. 

  


  
ProBiS H2O Plugin is **free for Academic (NON-COMMERCIAL) use**. For **COMMERCIAL use the potential users have to write to[Marko Jukic](mailto:marko.jukic@ffa.uni-lj.si)**. 

  


## ProBiS H2O Plugin Features

  1. Graphical User Interface
  2. lean and fast workflow
  3. fast iteration
  4. data logging and export
  5. usage of custom datasets



  


## Requirements & Dependencies

Python 2.7, NumPy (>= 1.6.1), SciPy (>= 0.9), scikit-learn (>= 1.18), PyMOL software, ProBiS (installed during ProBiS H2O setup). 

The plugin works in Unix like systems. Testing environment was Ubuntu 16.04. 

## Installation

ProBiS H2O is installed as a PyMOL plugin: 

  1. Run PyMOL
  2. Install the plugin in PyMOL by following the path: Plugin > Plugin Manager > Install New Plugin
  3. Restart PyMOL!



  
After plugin is installed, it can be run from PyMOL by following Plugin>ProBiS H2O and used through PyMOL Tcl-Tk GUI where plugin window supplements PyMOL default control and display windows. 

**Important note** : Before first time usage use SETUP DB button to set up th working environment. 

  


## First Time Usage

In the lower left corner of the ProBiS H2O main window resides SETUP DB button. Before usage ProBiS H2O plugin requires up-to-date information on RCSB sequence clustering and ProBiS algorithm binary file. When SETUP DB button is pressed ProBiS H2O plugin creates /Probis_H2O/ folder in default working directory of PyMOL (default is user home directory) and downloads recent sequence clustering data from RCSB PDB Database. ProBiS H2O plugin also creates /Probis_H2O/.pro/ folder where ProBiS algorithm binary file is downloaded to achieve a final workspace: 
    
    
      ~/Probis_H2O/.pro/probis
      ~/Probis_H2O/clusters50.txt
      ~/Probis_H2O/clusters70.txt
      ~/Probis_H2O/clusters90.txt
      ~/Probis_H2O/clusters95.txt
      ~/Probis_H2O/bc-30.out
      ~/Probis_H2O/bc-40.out
      ~/Probis_H2O/bc-50.out
      ~/Probis_H2O/bc-70.out
      ~/Probis_H2O/bc-90.out
      ~/Probis_H2O/bc-95.out
      ~/Probis_H2O/bc-100.out
    

Final step is to grant probis binary executable rights: 
    
    
      $sudo chmod +x ~/Probis_H2O/.pro/probis
    

## Usage

Plugin is composed of three respective tabs: 

  


#### Define System Tab

Default workflow proceeds in following steps (all the data manipulation is performed in /Probis_H2O/ working directory that was configured using SETUP DB button): 

[![Define system tab.PNG](/images/8/83/Define_system_tab.PNG)](/index.php/File:Define_system_tab.PNG)

  


  1. User inputs a desired PDB ID as examined system
  2. User presses “Find” button to display identified structures in cluster
  3. User downloads all relevant structures by pressing “Download” button (this step is in parentheses because this step is optional if the user already analysed a particular protein cluster and has the structures already downloaded. If it is a first time experiment, the “Download” button downloads the set from the PDB website)
  4. User identifies binding sites or chains of the examined system (1.)
  5. User selects desired binding site or chain from a list
  6. User proceeds with conserved water identification by pressing GO button



  


#### Cluster Analysis Tab

After the calculation is finished ProBiS H2O switches to the cluster analysis tab and displays the results 

[![Cluster analysis tab.PNG](/images/5/50/Cluster_analysis_tab.PNG)](/index.php/File:Cluster_analysis_tab.PNG)

  1. User displays the studied system in PyMOL display window by clicking “fetch/reset” button. This button can be used for resetting the display for alternative clustering representation
  2. User selects the cluster from calculated clusters window
  3. User displays the relevant cluster with display button



[![Feature 4duh.PNG](/images/a/ad/Feature_4duh.PNG)](/index.php/File:Feature_4duh.PNG)

#### About/Help Tab

The users can get reference information on the About/Help tab. 

  


## Reference

**If you are using ProBiS H2O in your work, please cite:**

Jukič, Marko, Janez Konc, Stanislav Gobec, Dušanka Janežič. Identification of Conserved Water Sites in Protein Structures for Drug Design. J. Chem. Inf. Model., 2017, doi: 10.1021/acs.jcim.7b00443 [[1]](https://pubs.acs.org/doi/10.1021/acs.jcim.7b00443)

## Licensing

ProBiS H2O software is copyrighted by : 

    **Marko Jukic**
    Fakulteta za farmacijo
    Univerza v Ljubljani
    Askerceva 7
    SI-1000 Ljubljana
    Slovenia.

The terms stated in the following agreement apply to all files associated with the software unless explicitly disclaimed in individual files. 
    
    
      # Copyright Notice
      # ================
      #
      # The PyMOL Plugin source code in this file is copyrighted, but you can
      # freely use and copy it as long as you don't change or remove any of
      # the copyright notices.
      #
      # ----------------------------------------------------------------------
      # This PyMOL Plugin is Copyright (C) 2017 by Marko Jukic <marko.jukic@ffa.uni-lj.si>
      #
      #                        All Rights Reserved
      #
      # Permission to use, copy and distribute
      # versions of this software and its documentation for any purpose and
      # without fee is hereby granted, provided that the above copyright
      # notice appear in all copies and that both the copyright notice and
      # this permission notice appear in supporting documentation, and that
      # the name(s) of the author(s) not be used in advertising or publicity
      # pertaining to distribution of the software without specific, written
      # prior permission.
      #
      # THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
      # INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
      # NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
      # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
      # USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
      # OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
      # PERFORMANCE OF THIS SOFTWARE.
      # ----------------------------------------------------------------------
    

Retrieved from "[https://pymolwiki.org/index.php?title=ProBiS_H2O&oldid=12760](https://pymolwiki.org/index.php?title=ProBiS_H2O&oldid=12760)"


---

## ProMOL

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://pymol.org/tmp/promol-plugin-54r461u.zip> (modified to work with PyMOL 2.x. For the unmodified version, see <http://www.promol.org>)  
Author(s)  | Paul Craig, Rochester Institute of Technology   
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
<http://www.promol.org>  
  
  
## ProMOL

ProMOL was developed as a program that searches for catalytic binding sites by using the advanced selection algebra from Warren DeLano's PyMOL software. ProMOL is capable of accurately recognizing catalytic sites nearly every time with few erroneous results that are easily distinguishable.  
  
For more information and downloads, please visit our site [ProMOL.org](http://www.promol.org/). 

Retrieved from "[https://pymolwiki.org/index.php?title=ProMOL&oldid=12742](https://pymolwiki.org/index.php?title=ProMOL&oldid=12742)"


---

## Proxy Config

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_This plugin is no longer available. For easy proxy configuration, see[proxy setting](/index.php/Fetch#Proxy_Setting "Fetch")._

## Description

This plugin which stores the users internet setting, i.e. whether they connect directly to the internet or through a proxy. 

## Download

Download by right clicking the link and choosing 'Save Link As...'.  
[Download here](http://vbc.med.monash.edu.au/~fauxn/pymol_plugin/ProxyConfig.py)  
Then run PyMol and under the 'Plugin' menu option choose 'Install Plugin...'. 

## Author

[Noel Faux](/index.php/User:Fauxn "User:Fauxn")

Retrieved from "[https://pymolwiki.org/index.php?title=Proxy_Config&oldid=11165](https://pymolwiki.org/index.php?title=Proxy_Config&oldid=11165)"


---

## PyANM

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/pyanm.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/pyanm.py)  
Author(s)  | [Yuan Wang](/index.php/User:YuanWang "User:YuanWang")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Dependencies
  * 3 Installation
  * 4 Usage
    * 4.1 Structure Handling
    * 4.2 ANM Parameters
    * 4.3 Build ANM
    * 4.4 Show Springs
    * 4.5 Color by MSFs
    * 4.6 Make Movies
    * 4.7 Draw Arrows
    * 4.8 Export Data
    * 4.9 Control Panel
  * 5 Contact Author



# Introduction

Elastic Network Models (ENM) have been successful in reproducing fluctuations for proteins of native conformations. ENM is a coarse-grained method for modeling protein dynamics, meaning in generally in ENM each residue is represented by one bead in its Alpha Carbon position. These beads are then connected by elastic springs if the distance between two beads fall under a cutoff value (usually within the range of 7 to 15 angstroms). For more details please read [here](http://www.cell.com/biophysj/abstract/S0006-3495%2801%2976033-X). 

PyANM is developed as a cross-platform Pymol Plugin to allow its users to build and visualize Anisotropic Network Models, a member of the ENM family. This plugin allows its users to draw arrows or make movies based on calculated mode motions, to draw all springs used to build the ANM within the protein, to color the protein based on its Mean Square Fluctuations (MSF) calculated from ANM and to export data from ANM for future processing. PyANM also allows its users to change the colors and scales (sizes) for movies or arrows generated by PyANM for better visualization. 

  


# Dependencies

  * PyANM requires Numerical Python (NumPy) to do all the matrix calculations. If you wish to use PyANM but don't have NumPy on your machine, please follow the download and build instructions [here](http://www.numpy.org/). If you will be using PyANM on a Windows machine, installing a Python distribution such as [Anaconda](https://store.continuum.io/cshop/anaconda/) might be easier. After installing NumPy, type the following line inside Pymol and if you don't see any error messages, you now have NumPy ready for Pymol! 
        
        import numpy
        




  


  * Having [SciPy](http://www.scipy.org/index.html#) might improve the performance of PyANM but it is not required.


  * PyANM has only been tested for Pymol 1.1 and higher, using an earlier version of Pymol might cause problems.



# Installation

Download pyanm.py from this page, save it on your system, then go to Pymol-->Plugin-->Manage Plugins-->Install..., select the file you have just downloaded and then restart Pymol. Things might be different depending on your Pymol version and your local system, but this should be fairly similar. 

[![](/images/e/ea/PyANM_interface.png)](/index.php/File:PyANM_interface.png)

[](/index.php/File:PyANM_interface.png "Enlarge")

**Figure 1** : PyANM Interface

After a successful installation, go to Pymol-->Plugin and you will see PyANM, click on it and you should see PyANM's interface just like in figure 1. 

# Usage

Now that you have PyANM installed, let's build one ANM together! 

### Structure Handling

  * Type the name of your protein structure under 'Select Structure', this could either be a valid 4-letter PDB (Protein Data Bank) ID or the name of a object you have already loaded into Pymol.
  * You could also choose to upload a local file from your system using the 'Browse' button.
  * If the check box 'Use only Alpha Carbons' is checked, even if you choose an all-atom structure, PyANM will coarse-grain the structure for you. Otherwise PyANM will build an all-atom ANM if an all-atom structure is used.
  * If the check box 'Include HETATMs' is checked, PyANM will include atoms in PDB that starts with HETATM as well.



### ANM Parameters

  * So far there are two ANMs available for you to choose in PyANM : the cutoff model as described [here](http://www.cell.com/biophysj/abstract/S0006-3495%2801%2976033-X) and the Parameter Free Model described [here](http://www.pnas.org/content/106/30/12347.short), where spring constants are adjusted based on distance between two residues.
  * If you choose to use **Cutoff Model** , you will be able to choose your desired cutoff value. The default value is 12 angstroms, personally I would recommend using any value between 7 angstroms to 15 anstroms.
  * If you choose to use **Parameter Free Model** , you will be able to choose the power t for calculating spring constants for all residue pairs, where spring constant k for residue i and j is calculated as 1/distance(i,j) to the power of t. The default value for power t is 6.



### Build ANM

Once you have selected the structure you wish to use as well as the model and its parameters, you should be able to build your ANM simply by clicking 'Build ANM'.A window will pop up telling you your ANM was built successfully once the calculation has finished. For a protein with less than 200 residues, this should take less than 10 seconds. 

[![](/images/d/de/Pyanm_example_of_springs.png)](/index.php/File:Pyanm_example_of_springs.png)

[](/index.php/File:Pyanm_example_of_springs.png "Enlarge")

**Figure 2** : All springs of HIV protease (PDB ID: 1T3R) when cutoff value is set as 12 angstroms

### Show Springs

After building your ANM, you will be able to visualize all the springs in Pymol by simply clicking 'Show Springs' in PyANM. 

See Figure 2 for a demonstration of all springs for HIV Protease (PDB ID: 1T3R) using cutoff model with cutoff value set as 12 angstroms. 

### Color by MSFs

[![](/images/f/f3/Pyanm_example_of_msf.png)](/index.php/File:Pyanm_example_of_msf.png)

[](/index.php/File:Pyanm_example_of_msf.png "Enlarge")

**Figure 3** : HIV protease (PDB ID: 1T3R) colored by its calculated MSFs

One can also calculate the Mean Square Fluctuations (MSF) with all the modes calculated by ANM. ANM has been successful in getting high correlations between its calculated MSFs and experimental b-factors for most cases. Here PyANM allows you to color the structure based on its calculated MSFs just like Pymol colors a structure with its b-factors. Here the spectrum is from blue to red corresponding to from low values to high values. 

Figure 3 shows the same HIV protease in figure 2 colored by its calculated MSFs. 

### Make Movies

[![](/images/e/ee/Pyanm_example_of_arrows.png)](/index.php/File:Pyanm_example_of_arrows.png)

[](/index.php/File:Pyanm_example_of_arrows.png "Enlarge")

**Figure 3** : Arrows showing mode motions for HIV Protease (PDB ID: 1T3R)

You will be able to make harmonic fluctuation movies for different modes when you click 'Make Movies' after building your ANM. A dialogue will pop up asking for the desired modes you wish to make movies with. Modes are sorted with an ascending frequencies and the first 6 rigid-body modes are already discarded. So modes 1,2,3 will be the first 3 slowest modes. 

### Draw Arrows

Instead of making movies for different modes, you can also draw arrows for each residue for different modes. Each arrow's direction will point to the direction where the residue will move and each arrow's length will indicate the moving scale for each residue. 

Figure 4 shows the arrows for mode 2 of HIV protease produced by PyANM. 

### Export Data

This will allow you to save data calculated from PyANM to your local system for future processing. You can choose to save eigenvalues (frequencies for modes), eigenvectors (modes), MSFs, contact matrix or hessian matrix from **Control Panel**. 

### Control Panel

This is where PyANM allows you to change the scales (sizes) of the movies (arrows), or change color of the arrows (The color of movie objects could simply be changed like you change colors for other objects in Pymol). You could either predefine scales (sizes) or colors you would like before building your ANM or you could change these settings after you have made movies or arrows and these changes will happen real-time. 

You could also select what attributes you would like to export to your local system here inside control panel. These data will be written as text files and Numpy-ready for future processing. 

# Contact Author

Thank you for using PyANM, please feel free to contact me at yuanwang (at) iastate (dot) edu about anything, especially if you are having problems with PyANM. 

Retrieved from "[https://pymolwiki.org/index.php?title=PyANM&oldid=11742](https://pymolwiki.org/index.php?title=PyANM&oldid=11742)"


---

## PyDet

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A PyMOL plugin for visualizing geometric concepts. 

Reference: [PyDet at _Bioinformation_. 2008; 2(8): 346–347.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2478735/)

Download: <http://pitgroup.org/downloads/>

Retrieved from "[https://pymolwiki.org/index.php?title=PyDet&oldid=9173](https://pymolwiki.org/index.php?title=PyDet&oldid=9173)"


---

## PyMod

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Description
  * 2 Requirements
  * 3 PyMod Download
  * 4 PyMod Install
  * 5 How to use PyMod
  * 6 References



## Description

**PyMod 3** is a PyMOL plugin, designed to act as simple and intuitive interface between PyMOL and several bioinformatics tools (i.e., PSI-BLAST, Clustal Omega, HMMER, MUSCLE, CAMPO, PSIPRED, and MODELLER). The current PyMod release, PyMod 3, has been extended with a rich set of functionalities that substantially improve it over its predecessor (PyMod 2.0), particularly in its ability to build homology models through the popular MODELLER package. 

Starting from the amino acid sequence of a target protein, users may take advantage of PyMod 3 to carry out the three steps of the homology modeling process (that is, template searching, target-template sequence alignment and model building) in order to build a 3D atomic model of a target protein (or protein complex). Additionally, PyMod 3 may also be used outside the homology modeling context, in order to extend PyMOL with numerous types of functionalities. Sequence similarity searches, multiple sequence-structure alignments and evolutionary/phylogeny analyses can all be performed in the PyMod 3/PyMOL environment. 

[![Pymod3](/images/9/93/Pymod3.png)](/index.php/File:Pymod3.png "Pymod3")

## Requirements

PyMod 3 runs on Windows, macOS and Linux. It is compatible with PyMOL versions >= 2.3 and has been developed and tested on both incentive PyMOL builds (distributed by Schrodinger) and open source builds. 

## PyMod Download

Download: <https://github.com/pymodproject/pymod/releases/download/v3.0/pymod3.zip>

## PyMod Install

Please refer to the PyMod 3 [User's Guide](https://github.com/pymodproject/pymod/releases/download/v3.0/pymod_users_guide.pdf) to learn how to install PyMod 3 on your system. **Long story short** : PyMod 3 can be installed as any other PyMOL plugin: 

1\. Download the [PyMod 3 plugin file](https://github.com/pymodproject/pymod/releases/download/v3.0/pymod3.zip) (a ZIP file named pymod3.zip) 

2\. Use the PyMOL plugin manager (Plugin -> Plugin Manager from the menu of its main window) to install it 

3\. The external tools of PyMod 3 can be obtained and configured through an easy-to-use installer dialog which can be launched by the PyMod 3 plugin (Help -> Install PyMod Components from the main menu of the PyMod plugin). The way to configure MODELLER in PyMod may vary according to your PyMOL version (in the User's Guide we cover different scenarios). The easiest way is to install Modeller using the Pymod 3 installer dialog which can be launched by the PyMod 3 plugin (a Modeller licence key is required). 

## How to use PyMod

A series of introductory tutorials and accurate, in-depth information about all the PyMod 3 functionalities can be found in its [User's Guide](https://github.com/pymodproject/pymod/releases/download/v3.0/pymod_users_guide.pdf). Documentation is also available for PyMod3 at the [GitHub site](https://github.com/pymodproject/pymod). 

Home page: <http://schubert.bio.uniroma1.it/pymod/index.html>

Download: <https://github.com/pymodproject/pymod/releases/download/v3.0/pymod3.zip>

GitHub page: <https://github.com/pymodproject/pymod>

## References

PyMod 3: a complete suite for structural bioinformatics in PyMOL. Janson G, Paiardini A. Bioinformatics. 2020 Oct 3:btaa849. doi: 10.1093/bioinformatics/btaa849. Online ahead of print. PMID: 33010156 

  
PyMod 2.0: improvements in protein sequence-structure analysis and homology modeling within PyMOL. Janson G, Zhang C, Prado MG, Paiardini A. Bioinformatics. 2017 Feb 1;33(3):444-446. doi: 10.1093/bioinformatics/btw638. PMID: 28158668 

PyMod: sequence similarity searches, multiple sequence-structure alignments, and homology modeling within PyMOL. Bramucci E, Paiardini A, Bossa F, Pascarella S. BMC Bioinformatics. 2012 Mar 28;13 Suppl 4(Suppl 4):S2. doi: 10.1186/1471-2105-13-S4-S2. PMID: 22536966 

Retrieved from "[https://pymolwiki.org/index.php?title=PyMod&oldid=13334](https://pymolwiki.org/index.php?title=PyMod&oldid=13334)"


---

## PyMOLProbity

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/jaredsampson/pymolprobity/raw/master/pymolprobity.tar.gz>  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [MIT](http://opensource.org/licenses/MIT)  
  
## Contents

  * 1 Description
  * 2 Installation
  * 3 Getting Started
    * 3.1 GUI
    * 3.2 Command-Line Interface



## Description

[![](/images/d/d9/PyMOLProbity_GUI.png)](/index.php/File:PyMOLProbity_GUI.png)

[](/index.php/File:PyMOLProbity_GUI.png "Enlarge")

The PyMOLProbity GUI can be used to inspect and adjust clashes and flip orientations of flippable side chain groups.

[PyMOLProbity](https://github.com/jaredsampson/pymolprobity) is a plugin allows the user to produce MolProbity-style visualization of atomic interactions within a structure (e.g. H-bonds, van der Waals interactions and clashes) directly within a PyMOL session. The plugin runs local copies of several executable programs from the [Richardson Lab](http://kinemage.biochem.duke.edu/) at Duke University, authors of the [MolProbity](http://molprobity.biochem.duke.edu/) software, parses the output, and displays the results in the PyMOL viewport. There are both a graphical user interface (GUI) for general point-and-click use, and a command-line interface (CLI) suitable for scripting. 

## Installation

For installation instructions, please see the [README](https://github.com/jaredsampson/pymolprobity/blob/master/README.md) file in the [repository](https://github.com/jaredsampson/pymolprobity). 

If you use [Anaconda](https://anaconda.org), you can install the dependencies with: 
    
    
    conda install -c speleo3 reduce probe flipkin prekin
    

## Getting Started

### GUI

Once PyMOLProbity is installed, it should appear as an option in PyMOL's _Plugin_ menu. Load or fetch a structure, and launch the GUI by selecting _Plugin > PyMOLProbity_. 

  * Use the _Add Hydrogens_ tab to add hydrogens with Reduce. This will also calculate which N/Q/H residue side chains should be flipped, and perform those flips. Note that this should be done even if the model already includes explicit hydrogens.
  * To examine these more closely, select the _Review Flips_ tab. Here, you can zoom to inspect each flippable residue and choose the ones you wish to keep or change. Save any changes using the _**Save Selections**_ button.
  * Finally, use the _Visualize Contacts_ tab to run Probe on the modified coordinates and generate contact dots and clash vectors for all the atoms in your object.



### Command-Line Interface

The plugin makes the following functions available: 

  * **reduce_obj**(obj, flip=1): Run reduce on a loaded PyMOL object (or named selection) with (default) or without making the Asn/Gln/His flips recommended by Reduce.


  * **flipkin_obj**(obj): Run Flipkin to create both NQ and H flipkin kinemage visualization of the Reduce-modified structure.


  * **probe_obj**(obj): Run Probe on either a structure saved from the Flipkin tab of the GUI, or the Reduce-modified structure.



Note that both `flipkin_obj` and `probe_obj` require previously having run `reduce_obj` on the same object. 

Retrieved from "[https://pymolwiki.org/index.php?title=PyMOLProbity&oldid=12911](https://pymolwiki.org/index.php?title=PyMOLProbity&oldid=12911)"


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

## Show contacts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/show_contacts.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/show_contacts.py)  
Author(s)  | [David Ryan Koes](/index.php?title=User:DavidKoes&action=edit&redlink=1 "User:DavidKoes \(page does not exist\)")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments



## Introduction

PyMOL Plugin for displaying polar contacts. Good hydrogen bonds (as determined by PyMOL) are shown in yellow. Electrostatic clashes (donor-donor or acceptor-acceptor) are shown in red. Close (<4.0 A) but not ideal contacts are shown in purple. Cutoffs are configurable. Exports the command `contacts` which takes two selections and an optional name for the generated contacts group. Alternatively, the selections can be chosen using a dialog box accessible from the Plugins menu. 

  


## Usage
    
    
    show_contacts sel1, sel2, result="contacts", cutoff=3.6, bigcutoff=4.0
    

## Required Arguments

  * **sel1** = first selection
  * **sel2** = second selection



  


## Optional Arguments

  * **result** = name of created group containing contacts
  * **cutoff** = cutoff for good contacts (yellow dashes)
  * **bigcutoff** = cutoff for suboptimal contacts (purple dashes)



Retrieved from "[https://pymolwiki.org/index.php?title=Show_contacts&oldid=13265](https://pymolwiki.org/index.php?title=Show_contacts&oldid=13265)"


---

## Spaceball plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <http://divites.de/files/pymol/spaceball.tar.gz>  
Author(s)  | [Raphael Dives](/index.php?title=User:RaphaelDives&action=edit&redlink=1 "User:RaphaelDives \(page does not exist\)")  
License  | [BSD-3-Clause](http://opensource.org/licenses/BSD-3-Clause)  
  
**Note:** Incentive PyMOL 2.1 has built-in SpaceNavigator support for all platforms, see <https://pymol.org/spacenavigator>

## Description

This plugin integrates Spaceball devices on Linux Systems. 

It was developed and tested with the SpaceNavigator(TM) by 3DConnexion. 

But it *should* work with other (Linux based) systems and Spacenav Driver compatible devices as well. 

The used software was: 

    PyMOL version 1.7
    Ubuntu-12.04 / 14.04

  


## Installation

[Download](http://divites.de/files/pymol/spaceball.tar.gz) the plugin and install it using the PyMOL plugin wizard. 

Please consider the following requirements: 

  1. Install the "libspnav0" and "libspnav-dev" packages from the repositories
  2. Obtain and install the current version of the Spacenav driver from <http://www.spacenav.sourceforge.net/>



  


## Important

The name "SpaceMouse" is registered trademark of 3DConnexion. 

This plugin uses the spnav-0.9 package 

    Copyright (c) 2011, Stanley Seibert
    All rights reserved.
    The BSD-License agreement and the Readme for this package are located in the spnav_docs folder.
    Please see the <http://bitbucket.org/seibert/spnav/> for more information.

Retrieved from "[https://pymolwiki.org/index.php?title=Spaceball_plugin&oldid=13401](https://pymolwiki.org/index.php?title=Spaceball_plugin&oldid=13401)"


---

## Stride

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/dssp_stride.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/dssp_stride.py)  
Author(s)  | [Hongbo Zhu](/index.php/User:Hongbo_zhu "User:Hongbo zhu")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[DSSP](http://swift.cmbi.ru.nl/gv/dssp/) and [Stride](http://webclu.bio.wzw.tum.de/stride/) are popular tools for assigning secondary structures for proteins. DSSP & Stride Plugin for PyMOL provides a graphical user interface for coloring proteins according to their secondary structures assigned by DSSP or Stride. Most recent DSSP & Stride Plugin code can be obtained from [this link](https://github.com/hongbo-zhu-cn/Pymol-script-repo/blob/master/plugins/dssp_stride.py). 

## News

  * 2014_01_24_Update: The plugin is able to cope with empty chain name now.


  * 2012_10_22_Update: The most recent version of DSSP (v2.0.4) does not consider residues after the first [`TER` record](http://www.wwpdb.org/documentation/format33/sect9.html#TER) from the same chain. The plugin has been updated to cope with this new feature.



## External links

[![](/images/9/9b/Demo_DSSP_plugin_1pyg.png)](/index.php/File:Demo_DSSP_plugin_1pyg.png)

[](/index.php/File:Demo_DSSP_plugin_1pyg.png "Enlarge")

Demonstration of DSSP plugin (pdb:1pyg).

  * Most recent code on github: [DSSP & Stride Plugin for PyMOL @ GitHub](https://github.com/hongbo-zhu-cn/Pymol-script-repo/blob/master/plugins/dssp_stride.py)


  * Original site (not updated anymore): [DSSP & Stride Plugin for PyMOL](http://www.biotec.tu-dresden.de/~hongboz/dssp_pymol/dssp_pymol.html)



## See Also

  * [dssp](/index.php/Dssp "Dssp") (psico)
  * [dss](/index.php/Dss "Dss")



Retrieved from "[https://pymolwiki.org/index.php?title=DSSP_Stride&oldid=12654](https://pymolwiki.org/index.php?title=DSSP_Stride&oldid=12654)"


---

## SuperSym

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/SuperSymPlugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/SuperSymPlugin.py)  
Author(s)  | [Stuart Ballard](/index.php?title=User:Srballard&action=edit&redlink=1 "User:Srballard \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![](/images/3/34/SuperSymExample.png)](/index.php/File:SuperSymExample.png)

[](/index.php/File:SuperSymExample.png "Enlarge")

Symmetry partners for 1hpv showing 6-1 screw axis

[![](/images/d/da/SuperSymExample2.png)](/index.php/File:SuperSymExample2.png)

[](/index.php/File:SuperSymExample2.png "Enlarge")

Full cell of symmetry partners with symmetry axes displayed

SuperSym is a PyMOL plugin providing a large number of tools for visualization of space groups; unit cells; and symmetry axes, operators, and partners. 

The original source code is available from <https://sourceforge.net/projects/supersym/>

## Contents

  * 1 Dependencies
  * 2 Bugs
  * 3 Acknowledgments
  * 4 Installing SuperSym
  * 5 Feedback
  * 6 The Menu



## Dependencies

  * [cctbx](/index.php/CCTBX "CCTBX")
  * numpy



**PyMOL, cctbx and numpy must all be compiled with the same Python distribution!** See [CCTBX](/index.php/CCTBX "CCTBX"). 

This plugin was developed in 2010 for PyMOL version 1.2r1 and has not been updated since. 

## Bugs

Symmetry axes are not defined for all space groups, and do not display properly for some. 

## Acknowledgments

Primary coding and development was done by Stuart Ballard. All comments, questions, and issues should be directed to him at srballard@wisc.edu. 

Code for unit cell and symmetry axis building is borrowed from scripts created by Robert Campbell and Ralf W. Grosse-Kunstleve, available at <http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/>. Some of this code has been modified for use in SuperSym. 

[FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") is utilized for some of SuperSym's graphics generation, with some modifications. 

## Installing SuperSym

To install SuperSym v1.2, download SuperSymPlugin12.py from <https://sourceforge.net/projects/supersym/>. In PyMOL, go to: 

  * Plugin > Manage Plugins > Install...



A file selector dialog will appear. Select SuperSymPlugin12.py. PyMOL will direct you to restart, and upon doing so SuperSym will be accessible through the Plugin menu. 

To use functions of SuperSym directly, without creating a drop-down menu, use the run command in PyMOL on SuperSymPlugin12.py. 

Note: previous errors resulting from incorrect naming of the plugin file have been resolved in v1.2. 

## Feedback

Please post any comments, complaints, bug fix requests, useful tricks, or cool adaptations of SuperSym here. 

## The Menu

  * **Default Symmetry Partner Set**
    * See **Build Symmetry Partners > Cell [0,0,0] (default)**
  * **Draw Unit Cell**
    * Creates a cgo object with unit cell axes as cylinders. This functions similarly to _show cell_ , but the cell axes are cylinders instead of lines, allowing for printing.
  * **Build Symmetry Partners >**
    * All options in this submenu generate sets of symmetry partners
    * **Cell [0,0,0] (default)**
      * Generates a suite of symmetry partners for a given object for the default unit cell, which is lattice position [0,0,0]
    * **Cell [x,y,z] (custom)**
      * Generates a suite of symmetry partners for a given object for a lattice position which you specify
    * **2x2x2 Block**
      * Generates 8 sets of symmetry partners for a given object, filling lattice positions [0,0,0] through [1,1,1]
    * **3x3x3 Block**
      * Generates 27 sets of symmetry partners for a given object, filling lattice positions [-1,-1,-1] through [1,1,1]. This option may take a long time to execute
    * **By Partner**
      * Generates only those symmetry partners which the user specifies by their defining symmetry operators
  * **Coloring >**
    * **Default Rainbow**
      * Colors all symmetry objects with a specified by their symmetry operations automatically
    * **Select color for each operation**
      * Select symmetry partners to color by their defining symmetry operation and select the color for each
    * **Select one color for custom set of operations**
      * Select a set of symmetry partners defined by symmetry operations and select one color for all of them
  * **Graphics >**
    * **Lines**
      * Convenience function to display symmetry partners as lines
    * **Ribbon**
    * Convenience function to display symmetry partners as ribbons
    * **Cartoon**
      * Convenience function to display symmetry partners as cartoons
    * **Sphere Surface (best for printing)**
      * Uses the findSurfaceResidues function and shows surface residues as spheres. If printing, this option saves at least 60% of materials relative to regular surfaces, with minimal loss in resolution
    * **Surface (high load render)**
      * Displays symmetry partners as surfaces. This option may take a very long time to execute
  * **Symmetry Axes >**
    * **Build Axes**
      * Builds all symmetry axes for the given object. This functionality will be customizable and extended in future versions
  * **Move symmetry partners**
    * Merely displays instructions for using built in hotkeys to move symmetry partners
  * **About**
    * Developer info
  * **Help**
    * Reference to this page



Retrieved from "[https://pymolwiki.org/index.php?title=SuperSym&oldid=13315](https://pymolwiki.org/index.php?title=SuperSym&oldid=13315)"


---

## SuperSymSource

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page contains source files for the [SuperSym](/index.php/SuperSym "SuperSym") plugin. For documentation and use instructions, see [SuperSym](/index.php/SuperSym "SuperSym"). 

## Source Code

See <https://github.com/Pymol-Scripts/Pymol-script-repo/blob/master/plugins/SuperSymPlugin.py>

Retrieved from "[https://pymolwiki.org/index.php?title=SuperSymSource&oldid=13314](https://pymolwiki.org/index.php?title=SuperSymSource&oldid=13314)"


---

## TCONCOORD

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Description

tCONCOORD is a program that predicts conformational flexibility based on geometrical considerations. The focus is on proteins, however, most non-protein residues like ligands or other organic compounds are treated correctly. You can use it for: 

1) Generating ensembles of protein structures 

2) Generating alternate conformations of protein parts (e. g. loops) 

3) Generating ensembles of organic compounds. 

The PyMOL plugin helps to setup tCONCOORD runs and provides visual support for the structure analysis in tCONCOORD. 

## Installation

1) You need a working gromacs installation (<ftp://ftp.gromacs.org/pub/gromacs/gromacs-3.3.3.tar.gz>) 

2) the tCONCOORD binaries (<http://wwwuser.gwdg.de/~dseelig/tconcoord.html>) 

3) and the plugin itself (<http://wwwuser.gwdg.de/~dseelig/cncplugin.html>) 

## Author

[Daniel Seeliger](/index.php/User:Dseelig "User:Dseelig")

Retrieved from "[https://pymolwiki.org/index.php?title=TCONCOORD&oldid=7900](https://pymolwiki.org/index.php?title=TCONCOORD&oldid=7900)"


---

## VASCo

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

Generation of protein-protein specific surface properties. 

VASCo Website: <http://genome.tugraz.at/VASCo/>

# Reference

<http://www.biomedcentral.com/1471-2105/10/32>

Retrieved from "[https://pymolwiki.org/index.php?title=VASCo&oldid=8487](https://pymolwiki.org/index.php?title=VASCo&oldid=8487)"


---

## Voronoia

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A quote from their site: _Voronoia is a program suite to analyse and visualize the atomic packing of protein structures. It is based on the Voronoi Cell method and can be used to estimate the quality of a protein structure, e.g. by comparing the packing density of buried atoms to a reference data set or by highlighting protein regions with large packing defects. Voronoia is also targeted to detect locations of putative internal water or binding sites for ligands. Accordingly, Voronoia is beneficial for a broad range of protein structure approaches. It is applicable as a standalone version coming with a user friendly GUI or, alternatively, as a Pymol Plugin. Finally, Voronoia is also available as an easy to use webtool to process user defined PDB-files or to asses precalculated packing files from DOPP, the regularly updated Dictionary of Packing in Proteins._

Website: [Voronoia](http://141.42.202.21/voronoia/index.php?site=home)

Retrieved from "[https://pymolwiki.org/index.php?title=Voronoia&oldid=8492](https://pymolwiki.org/index.php?title=Voronoia&oldid=8492)"


---

## Xpyder

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[xPyder](http://pubs.acs.org/doi/abs/10.1021/ci300213c) is a PyMOL plugin for analyzing coupled residued and their networks in proteins. 

You can download it [here](http://linux.btbs.unimib.it/xpyder/). 

Retrieved from "[https://pymolwiki.org/index.php?title=Xpyder&oldid=11359](https://pymolwiki.org/index.php?title=Xpyder&oldid=11359)"


---

