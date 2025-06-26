# Category: PyMol Integration

## CCTBX

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes how to use the [Computational Crystallography Toolbox (cctbx)](https://github.com/cctbx/cctbx_project) with PyMOL. 

cctbx and PyMOL need to be compiled with the same Python distribution, otherwise they won't be compatible. 

### Incentive PyMOL

Create a [conda](https://repo.continuum.io/miniconda/) environment with [schrodinger::pymol](https://anaconda.org/schrodinger/pymol) and [conda-forge::cctbx-base](https://anaconda.org/conda-forge/cctbx-base). 

On Linux and Windows: 
    
    
    conda create -n env1 schrodinger::pymol conda-forge::cctbx-base
    

On macOS, you also need [schrodinger::tk](https://anaconda.org/schrodinger/pymol) and [XQuartz](https://www.xquartz.org/): 
    
    
    conda create -n env1 schrodinger::pymol conda-forge::cctbx-base schrodinger::tk
    

### Open-Source PyMOL

See [CCTBX-fedora32](/index.php/CCTBX-fedora32 "CCTBX-fedora32"). 

Retrieved from "[https://pymolwiki.org/index.php?title=CCTBX&oldid=13326](https://pymolwiki.org/index.php?title=CCTBX&oldid=13326)"


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

## Category:PyMol Integration Arbitrary Graphics Objects

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  * [PL: The Protein Library](http://protlib.uchicago.edu/index.html)



_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:PyMol_Integration_Arbitrary_Graphics_Objects&oldid=3154](https://pymolwiki.org/index.php?title=Category:PyMol_Integration_Arbitrary_Graphics_Objects&oldid=3154)"


---

