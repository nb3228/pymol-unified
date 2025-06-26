# Category: Surfaces and Voids

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

## SURFNET

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A recipe for reading surfaces from [Roman Laskowski](http://www.ebi.ac.uk/~roman/)'s [SURFNET](https://www.ebi.ac.uk/thornton-srv/software/SURFNET/) program (for finding cavities in macromolecules) into PyMol for visualisation. 

1\. Create your surfaces in "CCP4" format in SURFNET. 
    
    
    Asides: 
    A. The "endianness" of SURFNET is set to big endian by default - see the 
    remarks about the SGI flag.  Change this if you're on a little endian machine, 
    e.g. LINUX/i386.
    B.SURFNET can be compiled against ccp4 version 5 and 6 libraries 
    by following the instructions in the SURFNET distribution and modifiying 
    the link lines at the end of ccp4link.scr to replace 
    
    $CLIB/libccp4.a 
    
    with 
    
    $CLIB/libccp4f.a $CLIB/libccp4c.a
    

2\. Use [Gerard Kleywegt](http://xray.bmc.uu.se/gerard)'s mapman from the [USF](http://structure.usc.edu/usf/) [RAVE](http://structure.usc.edu/usf/rave.html) package to convert the CCP4 density map to XPLOR format 

e.g. in a shell on LINUX: 
    
    
    $ lx_mapman
    
    MAPMAN > READ map1 gaps.den
    
    MAPMAN > WRITE map1 gaps.xplor XPLOR 
    

3\. Open the XPLOR map in PyMol 

4\. Generate a mesh or surface object from the map using isomesh or isosurface. 

e.g. on the PyMol command line: 
    
    
    isomesh gaps_mesh, gaps, 100.0
    

Retrieved from "[https://pymolwiki.org/index.php?title=SURFNET&oldid=12721](https://pymolwiki.org/index.php?title=SURFNET&oldid=12721)"


---

