# Category: Pymol-script-repo

## AAindex

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/aaindex.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/aaindex.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.aaindex](http://pymol.org/psicoredirect.php?psico.aaindex)  
  
## Contents

  * 1 Introduction
  * 2 Python Example
  * 3 PyMOL Example
  * 4 See Also



## Introduction

[![](/images/4/44/AAindexExample.png)](/index.php/File:AAindexExample.png)

[](/index.php/File:AAindexExample.png "Enlarge")

Hydrophobicity coloring with KYTJ820101

AAindex is a database of numerical indices representing various physicochemical and biochemical properties of amino acids and pairs of amino acids. See <http://www.genome.jp/aaindex/>

This script is a python parser for the AAindex flat files which will be downloaded from <ftp://ftp.genome.jp/pub/db/community/aaindex/>

The script provides two PyMOL commands (but can also be used without PyMOL). 

  * aaindex2b: Loads numerical indices from aaindex1 as b-factors into your structure
  * pmf: Potential of Mean Force (aaindex3)



## Python Example

Consider the script is called aaindex.py, it is placed somewhere in your PYTHONPATH and the aaindex flatfiles are found in the current directory. 
    
    
    import aaindex
    aaindex.init(path='.')
    
    aaindex.grep('volume')
    x = aaindex.get('KRIW790103')
    print(x)
    print(x.get('A'))
    
    aaindex.init(index='2')
    aaindex.grep('blosum')
    x = aaindex.get('HENS920102')
    print(x.get('A', 'K'))
    

## PyMOL Example

Solvent-accessible surface coloring by [Hydropathy index (Kyte-Doolittle, 1982)](https://www.genome.jp/dbget-bin/www_bget?aaindex:KYTJ820101)
    
    
    run aaindex.py
    aaindex2b KYTJ820101
    spectrum b, white yellow forest
    set surface_solvent
    show surface
    

## See Also

  * Protscale from [rTools](http://www.rubor.de/pymol_extensions_de.html) does a similar job in coloring by amino acid properties



Retrieved from "[https://pymolwiki.org/index.php?title=AAindex&oldid=13884](https://pymolwiki.org/index.php?title=AAindex&oldid=13884)"


---

## AKMT Lys pred

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/aKMT_Lys_pred.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/aKMT_Lys_pred.py)  
Author(s)  | [Troels Schwarz-Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Data background
    * 2.1 PDB structures in Sulfolobus Islandicus
    * 2.2 Methylated peptide
  * 3 Example
  * 4 References



## Overview

Chu _et al._ has compiled a large information of peptides which are methylated in _Sulfolobus Islandicus._ The mass-spec data is available at: <http://www.mcponline.org/content/15/9/2908/suppl/DC1> in 

  * Supplemental Table S2a revisied (.xlsx, 6.5 MB) - Supplemental Table S2a revisied : mcp.M115.057778-3.xlsx
  * Supplemental Table S2a revisied (.xlsx, 4.8 MB) - Supplemental Table S2a revisied : mcp.M115.057778-4.xlsx



The authors has tried to find correlation with amino acids in 1D, with 7 residues up and downstream to the lysines. The authors also tried to first predict secondary structure in the chain, and then find correlation. They made a weblogo figure <http://weblogo.threeplusone.com> with their findings. There seems to be a preference for helical structure and Glutamic acid (E), Isoleucin (I), Leucin (L) and Valine (V) seems to be dominant with 7 residues up/downstream to the methylated lysines. 

This script should try to make an analysis this for 3D, with residues up and downstream in angstrom distance to the residue of interest. 

## Data background

### PDB structures in Sulfolobus Islandicus

Following structures was found June 2017 

  * www.rcsb.org -> Search advanced
  * Biology -> Source organism -> Sulfolobus organism

PDB id | Info | Text | Authors | Released | Method | Residue Count | Macromolecule | Gene Name(s) | UniProtKB AC | NCBI - Accession Numbers   
---|---|---|---|---|---|---|---|---|---|---  
4G2D | Entity 1 containing Chain A | Crystal structure of the hyperthermophilic Sulfolobus islandicus PLL SisLac | Hiblot, J., Gotthard, G., Chabriere, E., Elias, M. | 10/24/2012 | X-ray Diffraction | 315 | Aryldialkylphosphatase | M164_0332 | C4KKZ9 | ADX84410   
5EWT | Entity 1 containing Chain A | Crystal structure of ExoIII endonuclease from Sulfolobus islandicus | Yan, Z., Yuan, Z.L., Gu, L.C., Shen, Y.L. | 11/23/2016 | X-ray Diffraction | 247 | Exodeoxyribonuclease III Xth | SiRe_0100 | F0NG49  | ADX84195   
5JWJ | Entity 1 containing Chain A | NMR solution structure of a thermophilic lysine methyl transferase from Sulfolobus islandicus | de Lichtenberg, C., Stiefler-Jensen, D., Schwarz-Linnet, T., She, Q., Teilum, K. |  5/24/2017 | Solution NMR | 172 | Protein-lysine N-methyltransferase |  |  |   
5F4H | Entity 1 containing Chain A, B, C, D, E, F | Archael RuvB-like Holiday junction helicase | Zhai, B., DuPrez, K.T., Doukov, T.I., Shen, Y., Fan, L. | 12/21/2016 | X-ray Diffraction | 505/3030 | Nucleotide binding protein PINc | LS215_1665 | C3MQK6 | ADX85498   
5FA8 | Entity 1 containing Chain A | SAM complex with aKMT from the hyperthermophilic archaeon Sulfolobus islandicu | Chu, Y., Zhu, Y., Chen, Y., Li, W., Zhang, Z., Liu, D., Wang, T., Ma, J., Deng, H., Liu, Z.J., Ouyang, S., Huang, L. | 6/29/2016 | X-ray Diffraction | 161 | Ribosomal protein L11 methyltransf |  |  |   
5FAD | Entity 1 containing Chain A | SAH complex with aKMT from the hyperthermophilic archaeon Sulfolobus islandicus | Chu, Y., Zhu, Y., Chen, Y., Li, W., Zhang, Z., Liu, D., Wang, T., Ma, J., Deng, H., Liu, Z.J., Ouyang, S., Huang, L. |  6/29/2016 | X-ray Diffraction | 161 | Ribosomal protein L11 methyltransf |  |  |   
3O27 | Entity 1 containing Chain A, B | The crystal structure of C68 from the hybrid virus-plasmid pSSVx | Contursi, P., D'Ambrosio, K., Pirone, L., Pedone, E.M., Aucelli, T., She, Q., De Simone, G., Bartolucci, S. | 1/19/2011 | X-ray Diffraction | 136 | Putative uncharacterized protein |  |  |   
3M1M | Entity 1 containing Chain A | Crystal structure of the primase-polymerase from Sulfolobus islandicus | Beck, K., Vannini, A., Cramer, P., Lipps, G. | 6/16/2010 | X-ray Diffraction | 335 | ORF904 |  |  |   
2K9I | Entity 1 containing Chain A, B | NMR structure of plasmid copy control protein ORF56 from sulfolobus islandicus | Weininger, U., Zeeb, M., Neumann, P., Low, C., Stubbs, M.T., Lipps, G., Balbach, J. | 10/20/2009 | Solution NMR | 110 | Uncharacterized protein ORF56 |  |  |   
3FT7 | Entity 1 containing Chain A, B | Crystal structure of an extremely stable dimeric protein from sulfolobus islandicus | Weininger, U., Zeeb, M., Neumann, P., Low, C., Stubbs, M.T., Lipps, G., Balbach, J. | 10/20/2009 | X-ray Diffraction | 110 | Uncharacterized protein ORF56 |  |  |   
1RNI | Entity 1 containing Chain A | Bifunctional DNA primase/polymerase domain of ORF904 from the archaeal plasmid pRN1 | Lipps, G., Weinzierl, A.O., Von Scheven, G., Buchen, C., Cramer, P. | 1/27/2004 | X-ray Diffraction | 216 | ORF904 |  |  |   
1RO0 | Entity 1 containing Chain A | Bifunctional DNA primase/polymerase domain of ORF904 from the archaeal plasmid pRN1- Triple mutant F50M/L107M/L110M SeMet remote | Lipps, G., Weinzierl, A.O., Von Scheven, G., Buchen, C., Cramer, P. | 1/27/2004 | X-ray Diffraction | 216 | ORF904 |  |  |   
1RO2 | Entity 1 containing Chain A | Bifunctional DNA primase/polymerase domain of ORF904 from the archaeal plasmid pRN1- Triple mutant F50M/L107M/L110M manganese soak | Lipps, G., Weinzierl, A.O., Von Scheven, G., Buchen, C., Cramer, P. | 1/27/2004 | X-ray Diffraction | 216 | hypothetical protein ORF904 |  |  |   
  
### Methylated peptide

Data from: 

  * mcp.M115.057778-3.xlsx
  * mcp.M115.057778-4.xlsx

PDB id | NCBI - Accession Number | 'Methylated peptides' 'Sequence' in -3.xlsx | 'Methylated peptides' ' Modification' in -3.xlsx s | 'Methylated peptides' 'Sequence' in -4.xlsx | 'Methylated peptides' ' Modification' in -4.xlsx s   
---|---|---|---|---|---  
4G2D | ADX84410 Sequence  | SIDEIADLFIHDIkEGIQATSNK   
kIADKGSFIGLDR   
ILmEEGVDPGkILIGHLGDTDNTDYIK   
kNGmSEEVIDIIFK   
SIDEIADLFIHDIkEGIQATSNK   
SIDEIADLFIHDIkEGIQATSNK  | K14(Methyl)   
K1(Methyl)   
M3(Oxidation); K11(Trimethyl)   
K1(Methyl); M4(Oxidation)   
K14(Dimethyl)   
K14(Trimethyl)  | mRIPLVGkEPIEAEDmGFTLIHEHLR  | M1(Oxidation); K8(Trimethyl); M16(Oxidation)   
  
## Example
    
    
    reinitialize
    viewport 1100,700
    # Imports
    import aKMT_Lys_pred
    
    # Load
    fetch 4g2d, async=0
    remove hydrogens
    remove solvent
    
    # Display
    preset.publication(selection='4g2d')
    center 4g2d
    
    # Get the residues with lysine
    return_list_resn_resi, return_list_resn_resi_sel = aKMT_Lys_pred.get_resis_from_resn("4g2d", "LYS")
    print(return_list_resn_resi)
    
    # Match resi with peptide sequence
    python
    peptides = [
        ["SIDEIADLFIHDIkEGIQATSNK", "K14"],
        ["kIADKGSFIGLDR", "K1"],
        ["ILmEEGVDPGkILIGHLGDTDNTDYIK", "K11"],
        ["kNGmSEEVIDIIFK", "K1"],
        ["mRIPLVGkEPIEAEDmGFTLIHEHLR", "K8"],
    ]
    python end
    
    # Find peptides
    return_list_resn_match, return_list_resn_match_sel = aKMT_Lys_pred.match_peptides(target_sel="4g2d", peptides=peptides)
    print(return_list_resn_match)
    
    # Make list of lysines not methylated
    lys_not_meth = [x for x in return_list_resn_resi if x not in return_list_resn_match]
    print(lys_not_meth)
    
    # Make stats
    aKMT_Lys_pred.get_resi_stats(target_sel="4g2d", residues=lys_not_meth, group_id="not_meth", verb=True)
    aKMT_Lys_pred.get_resi_stats(target_sel="4g2d", residues=return_list_resn_match, group_id="methylated", verb=True)
    

## References

For Mass-Spec data   
_aKMT Catalyzes Extensive Protein Lysine Methylation in the Hyperthermophilic Archaeon Sulfolobus islandicus but is Dispensable for the Growth of the Organism'  
_ Chu, Yindi; Zhu, Yanping; Chen, Yuling; Li, Wei; Zhang, Zhenfeng; Liu, Di; Wang, Tongkun; Ma, Juncai; Deng, Haiteng; Liu, Zhi-Jie; Ouyang, SongyingHuang, Li   
_Molecular & Cellular Proteomics_ **2016** , Vol 15, Issue 9, 2908-2923, DOI: 10.1074/mcp.M115.057778 

Retrieved from "[https://pymolwiki.org/index.php?title=AKMT_Lys_pred&oldid=13937](https://pymolwiki.org/index.php?title=AKMT_Lys_pred&oldid=13937)"


---

## AngleBetweenHelices

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/anglebetweenhelices.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/anglebetweenhelices.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.orientation](http://pymol.org/psicoredirect.php?psico.orientation)  
  
  * [![Orientation of a helix](/images/a/ad/Helix_orientation1.png)](/index.php/File:Helix_orientation1.png "Orientation of a helix")

Orientation of a helix 

  * [![Orientation of two helices, their angle that separates the two is 145.08 degrees as reported by the script](/images/5/5b/Helix_orientation2.png)](/index.php/File:Helix_orientation2.png "Orientation of two helices, their angle that separates the two is 145.08 degrees as reported by the script")

Orientation of two helices, their angle that separates the two is 145.08 degrees as reported by the script 




Calculate angle between alpha-helices or beta-sheets. There are four different methods implemented to fit a helix, two of them also work for sheets or loops. 

# Commands
    
    
    angle_between_helices selection1, selection2 [, method [, visualize ]]
    
    
    
    helix_orientation selection [, visualize [, sigma_cutoff ]]
    
    
    
    helix_orientation_hbond selection [, visualize [, cutoff ]]
    
    
    
    loop_orientation selection [, visualize ]
    
    
    
    cafit_orientation selection [, visualize ]
    

# Example
    
    
    import anglebetweenhelices
    
    fetch 2x19, async=0
    select hel1, /2x19//B/23-36/
    select hel2, /2x19//B/40-54/
     
    # just calculate/visualize orientation of single alpha-helix
    helix_orientation_hbond hel1
     
    # get angle between two helices
    angle_between_helices hel1, hel2
    angle_between_helices hel1, hel2, method=1
    angle_between_helices hel1, hel2, method=2
     
    # get angle between beta-sheets
    select sheet1, A/47-54/
    select sheet6, A/146-149/
    angle_between_helices sheet1, sheet6, method=loop_orientation
    angle_between_helices sheet1, sheet6, method=cafit_orientation
    

Output: 
    
    
    PyMOL>angle_between_helices hel1, hel2, method=cafit_orientation
    Using method: cafit_orientation
    Angle: 145.08 deg
    

## See Also

  * [angle_between_domains](/index.php/Angle_between_domains "Angle between domains")



Retrieved from "[https://pymolwiki.org/index.php?title=AngleBetweenHelices&oldid=13885](https://pymolwiki.org/index.php?title=AngleBetweenHelices&oldid=13885)"


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

## Annotate v

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/annotate_v.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/annotate_v.py)  
Author(s)  | [Xin Yu](/index.php/User:XinYu "User:XinYu")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Annotation Schemes
  * 3 Dependencies and Limitations
  * 4 How to use
  * 5 Contact



## Description

_Version: 1.0_

A simple-to-use script for annotation of VH or VL sequence of an antibody. It creates a Pymol object for each FR or CDR region. It utilizes the REST API interface of _Abnum_ (<http://www.bioinf.org.uk/abs/abnum/>) from Dr Andrew Martin's group at UCL 

## Annotation Schemes

Currently supports **Kabat, Chothia, Contact, IMGT** 1

1Definitions for Kabat, Chothia, Contact, and IIMGT are same as listed in the table (<http://www.bioinf.org.uk/abs/info.html#kabatnum>), except that for IMGT, H-CDR2 is defined as **H51-H57** in this script, as opposed to of **H51-H56** in the table. This slight revision generates result that matches that from _IMGT website_ (<http://www.imgt.org/>) 

## Dependencies and Limitations

1\. Import **request** module 

2\. Relies on internet connection to _Abnum_

3\. Incomplete VH or VL sequence might not be annotated 

## How to use

1\. Create a selection object for the V region, (such as **VH** and **VL** shown here). Only 1 V-region (VH OR VL) can be annotated each time. Alternatively, simply select the residues into the default **< sele>** group. 

  


[![](/images/8/82/Annotate_v_pic1.png)](/index.php/File:Annotate_v_pic1.png)

[](/index.php/File:Annotate_v_pic1.png "Enlarge")

Create a selection group for input seq, or use the default sele

2\. Copy the entire script and create a `.py` file. Run script. 

3\. The general syntax is: 
    
    
    annotate_v("selection_group_name", "scheme")
    

**selection_group_name:** name of the selection group where you have the input sequence. Must be single-letter amino acid coded. Either VH or VL but not both 

**scheme:** currently supports Kabat, Chothia, Contact, IMGT. Must be lowercase 

For example: 
    
    
    annotate_v("VH", "kabat") #input sequence from VH selection group, annotate using "kabat"
    annotate_v("sele", "imgt") #input sequence from the default sele group, annotate using "imgt". You can also try "contact", "chothia"
    

  
4\. In the output window, the script will print the FR and CDR regions (if found). It also automatically creates a selection group for each of the FR and CDR region. 

[![](/images/6/63/Annotate_v_pic2.png)](/index.php/File:Annotate_v_pic2.png)

[](/index.php/File:Annotate_v_pic2.png "Enlarge")

example output from command lines

[![](/images/7/72/Annotate_v_pic3.png)](/index.php/File:Annotate_v_pic3.png)

[](/index.php/File:Annotate_v_pic3.png "Enlarge")

new FR and CDR selection groups are created

  


## Contact

For bugs & questions, emai @ xinyu18018@gmail.com. Enjoy! 

Retrieved from "[https://pymolwiki.org/index.php?title=Annotate_v&oldid=13886](https://pymolwiki.org/index.php?title=Annotate_v&oldid=13886)"


---

## Apbsplugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes the "APBS Tools2.1" plugin, based on the [original version](https://sourceforge.net/projects/pymolapbsplugin/) by [Michael G. Lerner](/index.php/User:Mglerner "User:Mglerner"). 

For the new plugin in Incentive PyMOL 2.0, see [APBS Electrostatics Plugin](/index.php/APBS_Electrostatics_Plugin "APBS Electrostatics Plugin"). 

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/apbsplugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/apbsplugin.py)  
Author(s)  | [Michael G. Lerner](/index.php/User:Mglerner "User:Mglerner")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The APBS Tool2.1 plugin integrates the [APBS](http://www.poissonboltzmann.org) software package into PyMOL. Its primary purpose is electrostatic surface visualization. 

## Contents

  * 1 Software requirements
  * 2 Examples
    * 2.1 Example 1
  * 3 Saving default locations
    * 3.1 Via environment variables
    * 3.2 Via editing the plugin
  * 4 Troubleshooting
  * 5 References and LICENSES
    * 5.1 APBS
    * 5.2 pdb2pqr
    * 5.3 MALOC



## Software requirements

  * [apbsplugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/apbsplugin.py) (Install via PyMOL's plugin manager)
  * [APBS](https://github.com/Electrostatics/apbs-pdb2pqr/releases) (Ubuntu users: `apt-get install apbs`)
  * [pdb2pqr](https://github.com/Electrostatics/apbs-pdb2pqr/releases) (Ubuntu users: `apt-get install pdb2pqr`)



Please register at <http://www.poissonboltzmann.org/> to ensure continued support for the APBS-PDB2PQR software. 

## Examples

There is a nice tutorial on the APBS homepage: <http://apbs-pdb2pqr.readthedocs.io/en/latest/examples/using-pymol.html>

For further help, there is a [mailing list](https://sourceforge.net/projects/apbs/lists/apbs-users). 

### Example 1

Read about the protein here: <http://www.proteopedia.org/wiki/index.php/3ig7>   


Download: [examples/apbsplugin_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/apbsplugin_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/apbsplugin_1.pml>" highlight="python" />

Open the **APBS Tool2.1...** plugin. The executable should be found for itself: 

Click in **Selection to use** , and write **cdk2**.   
Click **Set grid**   
Click **Run APBS**   
To the left in **Molecular Surface** , click **Show**.   
Tadaaa! Hide it again, now show both _Positive_ and _Negative_ Isosurface 

## Saving default locations

### Via environment variables

Set the environment variables APBS_BINARY_DIR, APBS_WEB_DIR, APBS_PSIZE_DIR, APBS_PDB2PQR_DIR, and/or TEMP, and things should work properly as long as you're using the most recent version of the plugin. 

### Via editing the plugin

Open up the python file associated with the plugin (typically apbsplugin.py), look for the section near the top labeled "Global config variables" and change the values from None to a string. 

  


## Troubleshooting

See: [APBS#Troubleshooting](/index.php/APBS#Troubleshooting "APBS")

## References and LICENSES

### APBS

When using the APBS, you are "obligated" to register your use of the software. This will the authors able to require funding for further development.  
Register here, it takes 1 minutes: <http://www.poissonboltzmann.org/apbs/downloads>

Please acknowledge your use of APBS by citing:   
Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. _Electrostatics of nanosystems: application to microtubules and the ribosome._ **Proc. Natl. Acad. Sci. USA 98** , 10037-10041 2001. [doi:10.1073/pnas.181342398](http://dx.doi.org/10.1073/pnas.181342398)

Please acknowledge your use of the Holst group software by citing:   
M. Holst and F. Saied, _Multigrid solution of the Poisson-Boltzmann equation._ **J. Comput. Chem. 14** , 105-113, 1993.   
M. Holst and F. Saied, _Numerical solution of the nonlinear Poisson-Boltzmann equation: Developing more robust and efficient methods._ **J. Comput. Chem.** 16, 337-364, 1995.   
For PMG (the multigrid solver):   
M. Holst, _Adaptive numerical treatment of elliptic systems on manifolds._ **Advances in Computational Mathematics 15** , 139-191, 2001. [doi:10.1023/A:1014246117321](http://dx.doi.org/10.1023/A:1014246117321)   
For FEtk (the finite element solver):   
R. Bank and M. Holst, _A New Paradigm for Parallel Adaptive Meshing Algorithms._ **SIAM Review 45** , 291-323, 2003. [doi:10.1137/S003614450342061](http://dx.doi.org/10.1137/S003614450342061)

### pdb2pqr

When using the pdb2pqr, you are "strongly encourage" to register your use of the software. This will the authors able to require funding for further development.  
Register here, it takes 1 minutes: <http://www.poissonboltzmann.org/pdb2pqr/d/downloads>

### MALOC

Please acknowledge the use of MALOC and FETK by citing:   
M. Holst, _Adaptive numerical treatment of elliptic systems on manifolds._ **Advances in Computational Mathematics** , 15 (2001), pp. 139-191. 

Copyright and Terms of Use:   
GNU General Public License: <http://www.gnu.org/copyleft/gpl.html>
    
    
    This version of MALOC is distributed under the following guidelines:
    MALOC (Minimal Abstraction Layer for Object-oriented C) 
    Copyright (C) 1994-2010 Michael Holst 
    
    This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. 
    
    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
    
    You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Apbsplugin&oldid=12789](https://pymolwiki.org/index.php?title=Apbsplugin&oldid=12789)"


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

## B2transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/b2transparency.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/b2transparency.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
b2transparency can set transparency settings like [transparency](/index.php/Transparency "Transparency") or [sphere_transparency](/index.php/Sphere_transparency "Sphere transparency") for each atom scaled by b-factor (or any other numeric atomic property). 

[![B2.png](/images/d/d1/B2.png)](/index.php/File:B2.png)

## Usage
    
    
    b2transparency [ selection [, setting [, minimum [, maximum [, var ]]]]]
    

## Example
    
    
    # some dummy molecule with increasing b-factor along the chain
    fab AAAAAAAAAAAAA
    alter all, b=resv
    
    b2transparency all
    
    as sticks
    show surface
    
    set surface_color, gray
    

## See Also

  * [spectrum](/index.php/Spectrum "Spectrum")



Retrieved from "[https://pymolwiki.org/index.php?title=B2transparency&oldid=13887](https://pymolwiki.org/index.php?title=B2transparency&oldid=13887)"


---

## BbPlane

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/bbPlane.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/bbPlane.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Blaine Bell](/index.php?title=User:Bell&action=edit&redlink=1 "User:Bell \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO plane between the backbone atoms of two neighboring residues. This is to show the planarity of the atoms. The image style this is meant to represent can be found many places, like "Introduction to Protein Structure" by Branden and Tooze (2nd ed. pp. 8). 

  * [![Close up of planar atoms](/images/8/8d/BbPlane3.png)](/index.php/File:BbPlane3.png "Close up of planar atoms")

Close up of planar atoms 

  * [![A few more](/images/8/89/BbPlane1.png)](/index.php/File:BbPlane1.png "A few more")

A few more 

  * [![Global view](/images/b/b8/BbPlane2.png)](/index.php/File:BbPlane2.png "Global view")

Global view 




# Examples
    
    
    # download the source and save as bbPlane.py
    run bbPlane.py
    fetch 1cll
    # make planes for residues 4-9
    bbPlane i. 4-10
    

Retrieved from "[https://pymolwiki.org/index.php?title=BbPlane&oldid=13888](https://pymolwiki.org/index.php?title=BbPlane&oldid=13888)"


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

## Ccp4 contact

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ccp4_contact.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ccp4_contact.py)  
Author(s)  | [Gerhard Reitmayr and Dalia Daujotyte](/index.php?title=User:Dalyte&action=edit&redlink=1 "User:Dalyte \(page does not exist\)")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Example 1
  * 4 Getting a CONTACT file
    * 4.1 Install CCP4 - for Linux
    * 4.2 Use of CONTACT - for Linux



## Overview

[![](/images/6/64/HhaExample.png)](/index.php/File:HhaExample.png)

[](/index.php/File:HhaExample.png "Enlarge")

Interface residues (at cutoff <4A) in the 2c7r.pdb were found using NCONT, but similar results can be obtained using this script and CONTACT output. Usage of ccp4_contact script in PyMOL allows easy selection of residues and atoms listed in the output file. Interacting protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

The script selects residues and atoms from the list of the contacts found by CONTACT from CCP4 Program Suite (CONTACT analyses contacts between subsets of atoms in a PDB file). First, we run CONTACT on our pdb file to find interface residues. Then by using the ccp4_contact script in PyMOL we separately select residues and atoms listed in the output file. This generates two selections (atoms and residues) for each interacting chain, allowing quick manipulation of (sometimes) extensive lists in CONTACT log file. 

## Usage

ccp4_contact( contactsfile, selName1 = "source", selName2 = "target" ) 

  


## Example 1

First use CONTACT to find interface residues/atoms in the pdb file. Once you have the log file proceed to PyMOL. Make sure you import the ccp4_contact script first. 
    
    
    fetch 2c7r
    ccp4_contact 2c7r.contact, selName1=prot, selName2=dna
    

[![](/images/8/81/HhaI20example.png)](/index.php/File:HhaI20example.png)

[](/index.php/File:HhaI20example.png "Enlarge")

Quick and easy selection of interacting residues and atoms listed in the CONTACT log file. Protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

Download: [examples/ccp4_contact_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_contact_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_contact_1.pml>" highlight="python" />

## Getting a CONTACT file

### Install CCP4 - for Linux

Goto: <http://www.ccp4.ac.uk/download.php>   
Click: automated Downloads Pages   
Select: Linux, generic linux (x86)   
Select: Customized installation   
Select: Only CCP4 Program Suite, Executables -> Continue   
No additional packages -> Continue   
Download   


Extract for example to: **/home/YOU/Software/CCP** 4   
Then run:   

    
    
    $ /home/YOU/Software/CCP4/install.sh
    

write yes, read agreement, push y to agree license   
For sourcing scripts, say yes.   
See the changes to your environmental virables:   

    
    
    $ less ~/.bashrc
    

### Use of CONTACT - for Linux

See here for the CONTACT program and options: <http://www.ccp4.ac.uk/html/contact.html>   
Locate the pdb, and now run in terminal:   

    
    
    $ contact XYZIN 2c7r.pdb >> 2c7r.contact << eof   (#press enter)
    > MODE ISUB  (#press enter)
    > ATYPE NON-CARBON  (#press enter)
    > eof        (#press enter, and now the program runs, and shell saves to 2c7r.contact)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ccp4_contact&oldid=13889](https://pymolwiki.org/index.php?title=Ccp4_contact&oldid=13889)"


---

## Ccp4 ncont

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ccp4_ncont.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ccp4_ncont.py)  
Author(s)  | [Gerhard Reitmayr and Dalia Daujotyte](/index.php?title=User:Dalyte&action=edit&redlink=1 "User:Dalyte \(page does not exist\)")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 Getting a NCONT file
    * 4.1 Install CCP4 - for Linux
    * 4.2 Use of NCONT - for Linux



## Overview

[![](/images/6/64/HhaExample.png)](/index.php/File:HhaExample.png)

[](/index.php/File:HhaExample.png "Enlarge")

Interface residues (at cutoff <4A) in the 2c7r.pdb were found using NCONT. Usage of ccp4_ncont script in PyMOL allows easy selection of residues and atoms listed in ncont.log file. Interacting protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

The script selects residues and atoms from the list of the contacts found by NCONT from CCP4 Program Suite (NCONT analyses contacts between subsets of atoms in a PDB file). First, we run NCONT on our pdb file to find interface residues. Then by using the ccp4_ncont script in PyMOL we separately select residues and atoms listed in a ncont.log file. This generates two selections (atoms and residues) for each interacting chain, allowing quick manipulation of (sometimes) extensive lists in NCONT log file. 

This script works best for intermolecular contacts (when NCONT target and source selections don't overlap). If crystal contacts (NCONT parameter cell = 1 or 2) are included then additional coding is required to distinguish inter from intramolecular contacts. 

## Usage

ccp4_ncont( contactsfile, selName1 = "source", selName2 = "target" ) 

  


## Examples

First use NCONT to find interface residues/atoms in the pdb file. Once you have ncont.log file proceed to PyMOL. Make sure you import the ccp4_ncont script first. 
    
    
    fetch 2c7r
    ccp4_ncont 2c7r.ncont, selName1=prot, selName2=dna
    

[![](/images/8/81/HhaI20example.png)](/index.php/File:HhaI20example.png)

[](/index.php/File:HhaI20example.png "Enlarge")

Quick and easy selection of interacting residues and atoms listed in the NCONT log file. Protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

Download: [examples/ccp4_ncont_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_ncont_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_ncont_1.pml>" highlight="python" />

## Getting a NCONT file

### Install CCP4 - for Linux

Goto: <http://www.ccp4.ac.uk/download.php>   
Click: automated Downloads Pages   
Select: Linux, generic linux (x86)   
Select: Customized installation   
Select: Only CCP4 Program Suite, Executables -> Continue   
No additional packages -> Continue   
Download   


Extract for example to: **/home/YOU/Software/CCP** 4   
Then run:   

    
    
    $ /home/YOU/Software/CCP4/install.sh
    

write yes, read agreement, push y to agree license   
For sourcing scripts, say yes.   
See the changes to your environmental virables:   

    
    
    $ less ~/.bashrc
    

### Use of NCONT - for Linux

See here for the NCONT program and options:   
<http://www.ccp4.ac.uk/html/ncont.html>   
<http://www.ccp4.ac.uk/html/pdbcur.html#atom_selection>   
Locate the pdb, and now run in terminal:   

    
    
    $ ncont XYZIN 2c7r.pdb >> 2c7r.ncont << eof   (#press enter)
    > source A    (#press enter)
    > target C,D  (#press enter)
    > eof         (#press enter, and now the program runs, and shell saves to 2c7r.ncont)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ccp4_ncont&oldid=13890](https://pymolwiki.org/index.php?title=Ccp4_ncont&oldid=13890)"


---

## Ccp4 pisa

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ccp4_pisa.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ccp4_pisa.py)  
Author(s)  | [Gerhard Reitmayr and Dalia Daujotyte](/index.php?title=User:Dalyte&action=edit&redlink=1 "User:Dalyte \(page does not exist\)")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Overview

The script selects atoms from the list of the contacts found by PISA. First, we run PISA on our pdb file to find the interfaces. Then by using the ccp4_pisa script in PyMOL we separately select atoms for all interface types and individual interfaces. This generates many selections, two for each interface, allowing quick manipulation of (sometimes) extensive lists in PISA log file. 

## Usage

ccp4_pisa( pisafile ) 

  


## Example 1

The script parses the XML output files from the PISA service or command line tool. A short description of how to download the XML output files is available here <http://www.ebi.ac.uk/msd-srv/prot_int/pi_download.html>. 

(For example, the following URL downloads the interfaces in 2c7r.pdb <http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?2c7r>) 

Make sure you import the ccp4_pisa script first. 

Download: [examples/ccp4_pisa_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_pisa_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_pisa_1.pml>" highlight="python" />

Retrieved from "[https://pymolwiki.org/index.php?title=Ccp4_pisa&oldid=13891](https://pymolwiki.org/index.php?title=Ccp4_pisa&oldid=13891)"


---

## Center of mass

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/center_of_mass.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/center_of_mass.py)  
Author(s)  | [Henschel](/index.php?title=User:Henschel&action=edit&redlink=1 "User:Henschel \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Examples
  * 4 See Also



### Description

This script calculates the true center-of-mass (COM) or the center-of-geometry (COG) for a given selection and returns the x, y, z values in the form of a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") (rather than a CGO sphere). The benefit of using a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") is that it can be selected and used in calculations. In addition, this script also iteratively goes through all states of a selection if more than one state exists and appends the corresponding COM/COG values as states into the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). The script itself is quite simple yet robust enough to be applied in different settings. As well, the calculation of the COM/COG is handled independently from the formation of the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") and can be called as an independent function where applicable. 

### Usage
    
    
    com selection [,state=None [,mass=None [,object=None]]]
    

  


### Examples
    
    
    import center_of_mass
    fetch 1c3y, finish=1, multiplex=0
    
    com 1c3y, state=1
    #Create a pseudoatom representing the 1c3y COG and store it as "1c3y_COM"
    #The "1c3y_COM" object will contain 1 state only
    
    com 1c3y, state=1, object=COG
    #Create a pseudoatom representing the 1c3y COG and store it as "COG"
    #The "COG" object will contain 1 state only
    
    com 1c3y, state=1, object=COM, mass=1
    #Create a pseudoatom representing the 1c3y COM and store it as "COM"
    #The "COM" object will contain 1 state only
    
    com 1c3y, object=COM, mass=1
    #Create a single pseudoatom containing the COM for each state found in 1c3y and store it as "COM"
    #The "COM" object will contain MULTIPLE states!
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [COM](/index.php/COM "COM")
  * [get_extent](/index.php/Get_extent "Get extent")
  * [get_position](/index.php/Get_position "Get position")



Retrieved from "[https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892](https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892)"


---

## Centroid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/centroid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/centroid.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
    * 3.1 See Also



## Overview

Centroid is a small script that returns the value of the geometric center (or centroid) of your selection. It also can translate the object of that selection to the origin. 

## Syntax
    
    
    centroid (selection=PyMOLSelection), [center=boolean]
    

## Examples
    
    
    # get the centroid of the polymer
    import centroid
    fetch 4ins, async=0
    centroid polymer
    
    # move some 'ligand' to the origin
    centroid ligand, center=1
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [Center_Of_Mass](/index.php/Center_Of_Mass "Center Of Mass")



Retrieved from "[https://pymolwiki.org/index.php?title=Centroid&oldid=13893](https://pymolwiki.org/index.php?title=Centroid&oldid=13893)"


---

## Cgo arrow

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/cgo_arrow.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_arrow.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO arrow between two picked atoms. 

## Usage
    
    
    cgo_arrow [ atom1 [, atom2 [, radius [, gap [, hlength [, hradius [, color [, name ]]]]]]]]
    

## Example
    
    
    run cgo_arrow.py
    
    fetch 1rx1, async=0
    preset.pretty("*")
    
    cgo_arrow [34.9, 68.4, 19.1], A/164/O3X, gap=1.0
    

[![Cgo arrow example.png](/images/1/17/Cgo_arrow_example.png)](/index.php/File:Cgo_arrow_example.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_arrow&oldid=13894](https://pymolwiki.org/index.php?title=Cgo_arrow&oldid=13894)"


---

## Cgo grid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/cgo_grid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_grid.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![cgo_grid creates flowing mesh objects](/images/2/27/Cgo_grid.gif)](/index.php/File:Cgo_grid.gif "cgo_grid creates flowing mesh objects")

## Contents

  * 1 About cgo_grid
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 SEE ALSO



## About cgo_grid

**cgo_grid** will generate a flowing mesh object using the points provided or the current view. By default is will generate a flowing membrane. The shape is affected substantially by the arguments!  


  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage

**cgo_grid** has many arguments, but not all need to be set necessarily (see arguments or examples). 
    
    
    cgo_grid [ pos1 [, pos2 [, pos3 [, length_x [, length_z [, npoints_x [, npoints_z
    [, nwaves_x [, nwaves_z [, offset_x [, offset_z [, gain_x [, gain_z
    [, thickness [, color [, nstates [, startframe [, endframe
    [, mode [, view [, name [, quiet ]]]]]]]]]]]]]]]]]]]]]]
    

|   
---|---  
  
  


## Arguments
    
    
        pos1 = single atom selection (='pk1') or list of 3 floats {default: [0,0,0]}
    
        pos2 = single atom selection (='pk2') or list of 3 floats {default: [1,0,0]}
    
        pos3 = single atom selection (='pk3') or list of 3 floats {default: [0,0,1]}
    
        --> the plane is defined by pos1 (origin) and vectors to pos2 and pos3, respectively
    
        length_x = <float>: length of membrane {default: 30}
        length_z = <float>: length of membrane {default: ''} # same as length_x
    
        npoints_x = <int>: number of points(lines) along x-direction
                    {default: ''} #will be set to give a ~1 unit grid
        npoints_z = <int>: number of points(lines) along z-direction
                    {default: ''} #will be set to give a ~1 unit grid
                    {minimum: 1 # automatic}
    
        nwaves_x =   <float>: number of complete sin waves along object x-axis
                     {default: 2}
        nwaves_z =  <float>: number of complete sin waves along object z-axis
                    {default: ''} # same as nwaves_x
                    define separately to adjust number of waves in each direction
    
        offset_x = <float> phase delay (in degrees) of sin wave in x-axis
                 can be set to affect shape and starting amplitude {default: 0}
        offset_z = <float> phase delay (in degrees) of sin wave in z-axis
                 can be set to affect shape and starting amplitude
                 {default: ''} # same as  offset_x
        offset_x and offset_z can be used together to phase
        otherwise identical objects
    
        gain_x = <float>: multiplication factor for y-amplitude for x-direction
                 {default: 1}
        gain_z = <float>: multiplication factor for y-amplitude for z-direction
                 {default: ''} #=gain_x
    
        thickness = <float>: line thickness {default: 2}
    
        color = color name <string> (e.g. 'skyblue') OR
                rgb-value list of 3 floats (e.g. [1.0,1.0,1.0]) OR
                {default: ''} // opposite of background
                input illegal values for random coloring
    
        nstates =  <int>: number of states; {default: 60}
                   this setting will define how many states
                   the object will have (per wave) and how fluent and fast the
                   animation will be.
                   Higher values will promote 'fluent' transitions,
                   but decrease flow speed.
                       Note: Frame animation cycles thought the states one at a time
                       and needs to be set accordingly. Can also be used to phase
                       otherwise identical objects.
                   Set to 1 for static object {automatic minimum}
    
        startframe: specify starting frame <int> or set (='') to use current frame
                    set to 'append' to extend movie from the last frame {default: 1}
          endframe: specify end frame <int> or set (='') to use last frame
                    if 'append' is used for startframe,
                    endframe becomes the number of frames to be appended instead
                    {default: 1}
                    Note: if start- and endframe are the same, movie animation will
                    be skipped, the object will be loaded and can be used afterwards
    
        mode: defines positioning {default: 0}:
        0: pos1 is center
        1: pos1 is corner
    
        view {default: 0}:
        '0': off/ uses provided points to create CGO
        '1': overrides atom selections and uses current orienatation for positioning
             - pos1 = origin/center
             - pos2 = origin +1 in camera y
             - pos3 = origin +1 in camera z
    
        name: <string> name of cgo object {default: ''} / automatic
    
        quiet: <boolean> toggles output
    

| Concept sketch:  
[![cgo_grid concept sketch](/images/d/d1/Cgo_grid.png)](/index.php/File:Cgo_grid.png "cgo_grid concept sketch")  
---|---  
  
  


## Examples

The behaviour or shape of the cgo_grid object are substantially influenced by the arguments 
    
    
    delete all
    set_view (\
         0.263772517,   -0.113038681,    0.957937598,\
        -0.040853567,    0.990910411,    0.128179103,\
        -0.963716805,   -0.072944686,    0.256756991,\
         0.000000000,    0.000000000, -131.816467285,\
         0.000000000,    0.000000000,    0.000000000,\
       -50.008331299,  353.641235352,  -20.000000000 )
    
    #membrane
    cgo_grid color=blue
    
    #swimming worm, random color
    cgo_grid \
    pos1=[0,-5,0], pos2=[1,-5,1], pos3=[0,-5,1],\
    length_x=15,\
    npoints_z=1,\
    gain_x=2,\
    gain_z=0,\
    thickness=20,\
    color=3,\
    mode=1,\
    name="worm"
    
    #Moving Ladder
    cgo_grid \
    length_x=15,\
    pos1=[0,10,0], pos2=[0,10,1], pos3=[0,9,0],\
    npoints_x=2, npoints_z=30,\
    name="ladder"
    
    #Roof
    cgo_grid \
    nstates=1,\
    npoints_x=15,\
    npoints_z=15,\
    gain_x=20,\
    gain_z=20,\
    nwaves_x=0.5,\
    thickness=5,\
    color=cyan,\
    name="roof"
    
    #Boxes
    cgo_grid \
    pos1=[0,-10,0], pos2=[1,-10,0], pos3=[0,-10,1],\
    nstates=1,\
    npoints_x=50,\
    npoints_z=50,\
    nwaves_x=0,\
    color=[0.00 , 0.53 , 0.22],\
    thickness=5,\
    name="bottom"
    
    cgo_grid \
    nstates=1,\
    npoints_x=2,\
    npoints_z=2,\
    nwaves_x=0,\
    color=gray60,\
    thickness=10,\
    name="top"
    
    cgo_grid \
    pos1=[-15,-10,15], pos2=[-14,-10,15], pos3=[-15,-9,15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=-2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="front"
    
    cgo_grid \
    pos1=[-15,-10,-15], pos2=[-14,-10,-15], pos3=[-15,-9,-15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="back"
    
    set ray_trace_frames, 0
    set movie_loop,1
    mplay
    
    # play around with the ARGUMENTS! :-)
    

|   
---|---  
  
  


## SEE ALSO

  * [CgoCircle](/index.php/CgoCircle "CgoCircle")



Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938](https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938)"


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

## Color by conservation

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/color_by_conservation.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/color_by_conservation.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | Free   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview

This script reads an alignment object and colors the protein objects in the alignment by the sequence conservation found in the alignment. 

[![Cbcon.png](/images/9/94/Cbcon.png)](/index.php/File:Cbcon.png)

# Example Usage
    
    
    reinitialize
    import color_by_conservation
     
    # get some kinases
    fetch 1opk 3dtc 3p86 2eva 3efw, async=0
     
    # turn on the sequence viewer
    set seq_view
    
    # align them into the "algn" object
    for x in cmd.get_names(): cmd.align(x, "3efw and c. A", object="algn")
     
    # color
    color_by_conservation aln=algn, as_putty=1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_by_conservation&oldid=13895](https://pymolwiki.org/index.php?title=Color_by_conservation&oldid=13895)"


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

## Colorbydisplacement

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/colorbydisplacement.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorbydisplacement.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Bug in code
    * 2.1 Examples
  * 3 Example 1



## Introduction

This script allows you to color two structures by distance displacement between an Open and Closed form of a protein, as calculated by PyMol's internal distance command. The pairwise distance is calculated between all-atoms. The distance displacement values are stored as B-factors of these residues, which are colored by a _rainbow_ color spectrum, with blue specifying minimum and red indicating maximum. 

Do keep in mind, all original B-factors values are overwritten! 

There exist one version.   
ColorByDisplacement**All** is between All atoms in residues and is quite slow => 3-5 mins for a run. Ideal for sticks representation. 

**You have to specify which residues should be used in the alignment procedure, or it will take all residues as standard**

V.2 is implemented the 2011.01.06 - Due to a bug in coloring. 

## Bug in code

A bug in the boolean operator of the spectrum command has been found. This versions work for version 1.3 Educational product.   
For other versions of pymol, try to change (comment/uncomment) the **cmd.spectrum** line. The other spectrum line works for Open-Source PyMOL 1.2r3pre, Incentive product 

### Examples
    
    
    ColorByDisplacementAll O5NT, C5NT, super1=resi 26-355, super2=resi 26-355, doColor=t, doAlign=t
    
    ColorByDisplacementAll O5NT, C5NT, super1=resi 26-355, super2=resi 26-355, doColor=t, doAlign=t, AlignedWhite='no'
    

  * [![ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.](/images/a/a2/ColorByDisplacement-All-1.png)](/index.php/File:ColorByDisplacement-All-1.png "ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.")

ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement. 

  * [![ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.](/images/6/69/ColorByDisplacement-All-2.png)](/index.php/File:ColorByDisplacement-All-2.png "ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.")

ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement. 




Dark blue is low displacement, higher displacements are in orange/yellow/red.   
Residues used for alignment is colored white. Can be turned off in top of algorithm. Residues not in both pdb files is colored black 

## Example 1

Download: [examples/colorbydisplacement_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/colorbydisplacement_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/colorbydisplacement_1.pml>" highlight="python" />

Retrieved from "[https://pymolwiki.org/index.php?title=Colorbydisplacement&oldid=13897](https://pymolwiki.org/index.php?title=Colorbydisplacement&oldid=13897)"


---

## ColorByRMSD

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/colorbyrmsd.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorbyrmsd.py)  
Author(s)  | [Shivender Shandilya](/index.php/User:Shiven "User:Shiven"), [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Arguments
  * 4 Examples



## Introduction

This script allows you to color two structures by Root Mean Square Deviation (RMSD). The distances between aligned C-alpha atom pairs are stored as B-factors of these residues, which are colored by a color spectrum, with blue specifying the minimum pairwise RMSD and red indicating the maximum. Unaligned residues are colored gray. 

## Usage
    
    
    colorbyrmsd mobile, target [, doAlign [, doPretty [, guide [, method ]]]]
    

## Arguments

  * **mobile** = string: atom selection for mobile atoms
  * **target** = string: atom selection for target atoms
  * **doAlign** = 0 or 1: Superpose selections before calculating distances {default: 1}
  * **doPretty** = 0 or 1: Show nice representation and colors {default: 1}
  * **guide** = 0 or 1: Only use C-alpha atoms {default: 1}
  * **method** = align or super: Method to match atoms {default: super}



## Examples
    
    
    # example #1
    colorbyrmsd 1cbs, 1hmt, doAlign=1, doPretty=1
    # example #2
    colorbyrmsd 1eaz, 1fao, doAlign=1, doPretty=1
    

  * [![1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.](/images/5/57/ColorByRMSD_1cbs_1hmt.png)](/index.php/File:ColorByRMSD_1cbs_1hmt.png "1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.")

1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white. 

  * [![1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.](/images/b/b6/ColorByRMSD_1eaz_1fao.png)](/index.php/File:ColorByRMSD_1eaz_1fao.png "1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.")

1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white. 




Retrieved from "[https://pymolwiki.org/index.php?title=ColorByRMSD&oldid=13898](https://pymolwiki.org/index.php?title=ColorByRMSD&oldid=13898)"


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

## Cubes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/cubes.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cubes.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the **cubes** and **tetrahedra** commands which create CGOs with cube or tetrahedron representations for all atoms in selection. 

## Usage
    
    
    cubes [ selection [, name [, state [, scale [, atomcolors ]]]]]
    
    
    
    tetrahedra [ selection [, name [, state [, scale [, atomcolors ]]]]]
    

## Example
    
    
    fetch 1rx1, async=0
    cubes solvent & b < 50
    tetrahedra solvent & b > 50
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cubes&oldid=13899](https://pymolwiki.org/index.php?title=Cubes&oldid=13899)"


---

## Cyspka

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/cyspka.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cyspka.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 Overview
    * 1.2 Algorithm development
    * 1.3 Correction to article
  * 2 Example of use
  * 3 References



## Introduction

This script is an experimental surface cysteine pKa predictor.  
The script is solely based on the work by: 

_Predicting Reactivities of Protein Surface Cysteines as Part of a Strategy for Selective Multiple Labeling_.  
Maik H. Jacob, Dan Amir, Vladimir Ratner, Eugene Gussakowsky, and Elisha Haas.  
**Biochemistry**. Vol 44, p. 13664-13672, [doi:10.1021/bi051205t](http://dx.doi.org/10.1021/bi051205t)

Questions to the article should be send to [Maik Jacob](http://www.jacobs-university.de/node/1816). 

[![Cyspka.png](/images/2/2b/Cyspka.png)](/index.php/File:Cyspka.png)

[](/index.php/File:Cyspka.png "Enlarge")

### Overview

The authors Jacob _et al_. were able to describe a computational algorithm that could predict the reactivity of surface cysteines. The algorithm was based on reaction rates with Ellmans reagent, Riddles _et al_., on 26 single cysteine mutants of adenylate kinase. The authors could predict the reactivity of the cysteines with a pearson correlation coefficient of 0.92. The algorithm was based on predicting the pKa values of cysteines by a calculation of electrostatic interactions to the backbone and sidechains of the protein and a energetic solvation effect from the number of atom neighbours. The algorithm is different from other pKa algorithms, since it calculates a Boltzmann energy distribution for the rotational states of cysteine. The reaction rate with Ellman's reagent was set proportional to the fraction of negatively charged cysteines, Bulaj _et al_. 

The authors Ratner _et al_., used the prediction of cysteine reactivity to selectively label several two cysteine mutants of adenylate kinase. A double mutant was selected with a high and a low reactive cysteine. A dye was first reacted with the most reactive cysteine and purified. Subsequent the second dye was attached to the low reactive cysteine under unfolding conditions. The strategy is interesting since in meet the double challenge of both site-specificity and double labeling of proteins. 

The prediction algorithm was tested against a popular, (see Sanchez _et al_., Sundaramoorthy _et al_.), web-based pKa prediction program [PROPKA](http://propka.ki.ku.dk/). [ A script](/index.php/Propka "Propka") was developed to to send a structure from PyMOL to the server and fetch and display the calculated pKa values in PyMOL. The adenylate kinase protein was virtual mutated at the 26 positions described by Jacob et al. and showed a lower pearson correlation coefficient of 0.7. In collaboration with Dr. Maik Jacob, the computational algorithm was developed as a general cysteine prediction algorithm for PyMOL. The script is here published at the pymolwiki, but the phase of validating the algorithm against other experimental pKa values has not been finished. The validation phase was hindered by the limited amount of available cysteine pKa data. For example revealed a search in the "[PPD a database of protein ionization constants](http://www.ddg-pharmfac.net/ppd/PPD/pKasearch.htm)" only 12 canditates, where several was dimerized and only 2 candidate cysteines were surface exposed. 

### Algorithm development

The algorithm is based on electrostatic calculations, where some parameters have been fine-tuned.   
The distance from the sulphur atom (SG) of the cysteine to the nearest backbone amide groups and residues with a partial charge, is considered in the electrostatic model.   
The model is including a evalution of Boltzman distribution of the rotation of the SG atom around the CA->CB bond. 

Twenty-six mutants of Escherichia coli adenylate kinase (4AKE) were produced, each containing a single cysteine at the protein surface, and the rates of the reaction with [Ellman's reagent](http://en.wikipedia.org/wiki/Ellman's_reagent) were measured. The reaction rate was set proportional to the pKa, to fine-tune the parameters in the electro static model. 

### Correction to article

There is a type error in equation 6. There is missing a minus "-". The equations should read:  
W M C , S C ( i ) = − ( ∑ W M C ( i ) + ∑ W S C ( i ) ) {\displaystyle W_{MC,SC(i)}=-\left(\sum W_{MC(i)}+\sum W_{SC(i)}\right)} ![{\\displaystyle W_{MC,SC\(i\)}=-\\left\(\\sum W_{MC\(i\)}+\\sum W_{SC\(i\)}\\right\)}](https://wikimedia.org/api/rest_v1/media/math/render/svg/7d3970750d5609b67ea693104acdba3f8e6910f5)

## Example of use

[Escherichia coli adenylate kinase.](http://www.proteopedia.org/wiki/index.php/4ake)   

    
    
    reinitialize
    import cyspka
    
    fetch 4AKE, async=0
    create 4AKE-A, /4AKE//A and not resn HOH
    delete 4AKE
    hide everything
    show cartoon, 4AKE-A
    cyspka 4AKE-A, A, 18
    
    ### You can loop over several residues. 
    loopcyspka 4AKE-A, A, residue=18.25.41-42
    
    ### OR for the original 26 residues. Takes a long time, so not to many at the time.
    #loopcyspka 4AKE-A, A, residue=18.25.41-42.55.73.90.113.162.188-189.203.28.58.75.102.138.142.148.154.169.214.3.24.86.109
    

## References

  * _Ellman's Reagent: 5,5'-Dithiobis(2-nitrobenzoic Acid) a Reexamination_. Peter W. Riddles, Robert L. Blakeley, and Burt Zerner. **Analytical Biochemistry**. Vol 94, p. 75-81, 1979
  * _Ionization Reactivity Relationships for Cysteine Thiols in Polypeptides_. Grzegorz Bulaj, Tanja Kortemme, and David P. Goldenberg. **Biochemistry**. Vol 37, p. 8965-8972, 1998
  * _A General Strategy for Site-Specific Double Labelling of Globular Proteins for Kinetic FRET Studies_. V. Ratner, E. Kahana, M. Eichler, and E. Haas. **Bioconjugate Chem.**. Vol 13, p. 1163-1170, 2002
  * _Prediction of reversibly oxidized protein cysteine thiols using protein structure properties_. Ricardo Sanchez, Megan Riddle, Jongwook Woo, and Jamil Momand. **Protein Science**. Vol 17, p. 473-481, 2008
  * _Predicting protein homocysteinylation targets based on dihedral strain energy and pKa of cysteines_. Elayanambi Sundaramoorthy, Souvik Maiti, Samir K. Brahmachari, and Shantanu Sengupta. **Proteins**. December, p. 1475-1483, 2007 [doi:10.1002/prot.21846](http://dx.doi.org/10.1002/prot.21846)



Retrieved from "[https://pymolwiki.org/index.php?title=Cyspka&oldid=13900](https://pymolwiki.org/index.php?title=Cyspka&oldid=13900)"


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

## Displacementmap

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/displacementmap.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/displacementmap.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Bugs
  * 3 Example
  * 4 Output
  * 5 Text output
  * 6 Example 1
  * 7 References



## Overview

DisplacementMap is made for easy investigations of suitable positions for site-directed mutagenesis of amino residues into cysteines and FRET/EPR pair labelling. 

A Open and Closed form of a protein should be loaded. New objects should be created for the selected asymmetric unit. Parts of the protein should be aligned, leaving the flexible part in two different positions. 

The input is the objects, Open (molecule1) and Closed (molecule2).   
Further is the criteria for selecting which atom the distance should be calculated between. Standard is atom='CA' (atom).   
Then one selects the Förster distance R0 (mindist). This is the minimum distance between the residues. This depends of the selection of the FRET pair and protein at hand. But usually in the range 40 - 80 Angstrom is suitable.   
Then one defines the minimum displacement that is accepted. Usually R0/2 (mindelta).   
The script will find the 5 best (listlength=5) positive and negative distance displacement between the two objects. 

It parses the results back to Pymol, that is standard set to show it as sticks (showsticks='yes').   
If one is looking for a particular residue range for the FRET pair, this can be specified with two input. resi1=24.45-47.86 resi2=100-105.107 resi1 is "from" and resi2 is "to". Individual residues are split by a ".", and ranges are defined with "-".  
In the end, it makes a large data-matrix with all the distances. It also produces a gnuplot file, for easy visualisation. Just drag the .plt file for win gnuplot command window and it plots your datamatrix. 

## Bugs

If the criterion is set to low, the memory gets flooded in the data-matrix file, making the file unreadable. No solutions found yet. 

## Example
    
    
    dispmap(molecule1="NIL", molecule2="NIL", mindist=30.0, mindelta=15.0, resi1=str(0), resi2=str(0), atom='CA', listlength=5, showsticks='yes'):
    

Use of functions 
    
    
    import displacementmap
    dispmap Open5NT, Closed5NT, 40.0, 15.0, resi1=206, resi2=1-512.515-550 
    dispmap Open5NT, Closed5NT, 30.0, 15.0, resi2=1-512.515-550, atom=CA, listlength=10
    

## Output

Suggestions are created in pymol, and gnuplot file is created for easy visualisation of pair data-matrix and the general backbone displacement. 

  * [![O5NT-C5NT-CA-dist.png](/images/1/1a/O5NT-C5NT-CA-dist.png)](/index.php/File:O5NT-C5NT-CA-dist.png)

  * [![O5NT-C5NT-CA-dist-menu.png](/images/c/cd/O5NT-C5NT-CA-dist-menu.png)](/index.php/File:O5NT-C5NT-CA-dist-menu.png)

  * [![O5NT-C5NT-CA-dist-gnuplot.png](/images/6/66/O5NT-C5NT-CA-dist-gnuplot.png)](/index.php/File:O5NT-C5NT-CA-dist-gnuplot.png)

  * [![O5NT-C5NT-CA-dist-backbone.png](/images/c/ce/O5NT-C5NT-CA-dist-backbone.png)](/index.php/File:O5NT-C5NT-CA-dist-backbone.png)




## Text output

In the data-matrix.txt file, you find the best suggestions 
    
    
    # Input 1: Open5NT  and Input 2: Closed5NT
    # Find for: CA  with min. residue-residue dist: 30.0 Angstrom
    # Looking for min. displacement dist: 15.0 Angstrom
    # I give nr# suggestions: 5, and do I show sticks in pymol?: yes
    # I look for suggestions in the range: ([0]=>means all residues)
    # for Input 1: ['0'] and for Input 2: ['0']
    # Mutation info is BLOSUM62 log-odds likelihood score and PAM250 is probability in % for evolutionary distance
    ###########################################################################################################
    # Max Negative and positive distances                                       #       Mutation info         #
    ###########################################################################################################
    # Obj.1   Obj.2   Delta   Op-Op Cl-Cl # Obj.1   Obj.2   Delta   Op-Op Cl-Cl # Res.1  Res.2 # Res.1  Res.2 #
    # Res.1   Res.2   -Dist   Dist  Dist  # Res.1   Res.2   +Dist   Dist  Dist  # B62/PAM250%  # B62/PAM250%  #
    ###########################################################################################################
    # PRO241  ASP456  -25.7   59.1   33.4 # PRO274  PRO513   26.8   31.2   58.0 # -3/ 2 -3/ 1  # -3/ 2 -3/ 2  #
    # LYS197  ASP456  -25.6   57.3   31.7 # THR273  PRO513   26.1   31.6   57.7 # -1/ 1 -3/ 1  # -1/ 2 -3/ 2  #
    # PRO513  ASP456  -25.4   32.4    7.0 # PRO274  GLY514   24.8   32.9   57.6 # -3/ 2 -3/ 1  # -3/ 2 -3/ 2  #
    # LEU198  ASP456  -25.3   59.0   33.7 # PRO274  LYS512   24.7   30.3   55.0 # -1/ 1 -3/ 1  # -3/ 2 -1/ 1  #
    # GLN201  ASP456  -25.2   62.8   37.6 # ASN311  PRO513   24.7   35.6   60.3 # -3/ 1 -3/ 1  # -3/ 1 -3/ 2  #
    

The script also automatically make the gnuplot plot file (.plt), with all the defined variables, for easy visualisation of the data-matrix.txt and the backbone displacement. 

## Example 1

Download: [examples/displacementmap_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/displacementmap_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/displacementmap_1.pml>" highlight="python" />

## References

For EPR considerations   
_Conformation of T4 Lysozyme in Solution. Hinge-Bending Motion and the Substrate-Induced Conformational Transition Studied by Site-Directed Spin Labeling_   
Hassane S. Mchaourab, Kyoung Joon Oh, Celia J. Fang, and Wayne L. Hubbell   
_Biochemistry_ **1997** , 36, 307-316 

_Probing Single-Molecule T4 Lysozyme Conformational Dynamics by Intramolecular Fluorescence Energy Transfer_   
Yu Chen, Dehong Hu, Erich R. Vorpagel, and H. Peter Lu   
J. Phys. Chem. B **2003** , 107, 7947-7956 

For FRET pair selection and considerations   
_Fluorescent probes and bioconjugation chemistries for single-molecule fluorescence analysis of biomolecules_   
Achillefs N. Kapanidisa and Shimon Weiss   
_Journal of chemical physics_ VOLUME 117, Number 24 22 **December 2002**

**For inspiration to DisplacementMap. Fig: 6, Difference-distance matrix for the difference in CA-CA distances.**   
_Structure of a Hinge-bending Bacteriophage T4 Lysozyme Mutant, Ile3 - > Pro_   
M. M. Dixon, H. Nicholsont, L. Shewchuk W. A. Baase and B. W. Matthews1   
J. Mol. Biol. (**1992**) 227. 917-933 

Retrieved from "[https://pymolwiki.org/index.php?title=Displacementmap&oldid=13901](https://pymolwiki.org/index.php?title=Displacementmap&oldid=13901)"


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

## Draw Protein Dimensions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/Draw_Protein_Dimensions.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/Draw_Protein_Dimensions.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw the dimensions of a protein based on an Inertia Axis Aligned Bounding Box (IABB). 

  
The idea behind this script is to calculate an approximate minimal bounding box (MBB) to extract the cell dimensions of a protein. To calculate the MBB is not trivial and usually the Axis Aligned Bounding Box (AABB) does not show up the real dimensions of the protein. This script calculates the inertia tensor of the object, extract the eigenvectors and use them to rotate the molecule (using as rotation matrix the transpose of the eigenvector matrix). The result is a molecule oriented with the inertia axis aligned with the cartesian axis. A new Bounding Box is calculated, which is called Inertia Axis Aligned Bounding Box (IABB), with a volume always lower than the AABB volume, and in many cases may correspond with the MBB. 

As always with these type of things, you have to use at your own risk. I did not try all the possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado (at) gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run path-to-the-script/Draw_Protein_Dimensions.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_Protein_Dimensions _selection_**

This will draw the cell dimensions of your selection based on a IABB box. It also generates the IABB box and the inertia axis, you just need to do "show cgo" to display them. 

You could also try: 

**draw_BB _selection_**

This will draw the AABB and IABB boxes with their cell dimensions and show up their volumes, you can compare them. 

  


# Examples
    
    
    # download the source and save as Draw_Protein_Dimensions.py
    run Draw_Protein_Dimensions.py
    fetch 2vak
    # calculate the dimensions of the full ASU
    draw_Protein_Dimensions 2vak
    

  * [![Dimensions of 2vak ASU](/images/7/77/2vak_ASU.png)](/index.php/File:2vak_ASU.png "Dimensions of 2vak ASU")

Dimensions of 2vak ASU 



    
    
    # you can extract only one chain and calculates the dimensions
    draw_Protein_Dimensions obj01
    

  * [![Dimensions of one protomer \(chain A - 2vak\)](/images/e/ed/2vak_A.png)](/index.php/File:2vak_A.png "Dimensions of one protomer \(chain A - 2vak\)")

Dimensions of one protomer (chain A - 2vak) 



    
    
    # You can also draw the Bounding boxes (AABB and IABB) to compare them.
    fetch 2vak
    draw_BB 2vak
    

  


  * [![Axis-aligned bounding box](/images/4/46/2vak_AABB.png)](/index.php/File:2vak_AABB.png "Axis-aligned bounding box")

Axis-aligned bounding box 

  * [![Inertia-axis-aligned bounding box](/images/4/41/2vak_IABB.png)](/index.php/File:2vak_IABB.png "Inertia-axis-aligned bounding box")

Inertia-axis-aligned bounding box 




Retrieved from "[https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902](https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902)"


---

## DrawGridBox

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/drawgridbox.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/drawgridbox.py)  
Author(s)  | [Cunliang Geng](/index.php?title=User:Clgeng&action=edit&redlink=1 "User:Clgeng \(page does not exist\)")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the [drawgridbox](/index.php?title=Drawgridbox&action=edit&redlink=1 "Drawgridbox \(page does not exist\)") command, which draw grid boxes around a selection. 

## Usage
    
    
     drawgridbox [ selection,  [nx, [ny, [nz, [padding, [lw, [r, [g, [b]]]]]]]]]
    

  


## Prameters

  * **selection** = str: atom selection {default: all}
  * **nx** = int: number of grids on axis X {default: 10}
  * **ny** = int: number of grids on axis Y {default: 10}
  * **nz** = int: number of grids on axis Z {default: 10}
  * **padding** = float: default 0.0
  * **lw** = float: line width {default: 2.0}
  * **r** = float: red color component, valid range is [0.0, 1.0] {default 1.0}
  * **g** = float: green color component, valid range is [0.0, 1.0] {default 1.0}
  * **b** = float: blue color component, valid range is [0.0, 1.0] {default 1.0}



  


## Examples
    
    
    load drawgridbox.py
    fetch 1cbh
    show surface
    drawgridbox 1cbh, nx=5, ny=5, nz=5,  lw=1, g=0, b=0
    

[![Drawgridbox.png](/images/5/5a/Drawgridbox.png)](/index.php/File:Drawgridbox.png)

  

    
    
    drawgridbox 1cbh, nx=5, ny=5, nz=5, padding=5, lw=1, r=0, g=0
    

[![Drawgridbox padding.png](/images/3/32/Drawgridbox_padding.png)](/index.php/File:Drawgridbox_padding.png)

Retrieved from "[https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904](https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904)"


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

## Dssr block

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [script/dssr_block.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/script/dssr_block.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![Trna.png](/images/d/dd/Trna.png)](/index.php/File:Trna.png)

[](/index.php/File:Trna.png "Enlarge")

This script adds the dssr_block command, which is a simple wrapper for the **DSSR** program and can create "block" shaped cartoons for nucleic acid bases and base pairs. 

_Requires the**x3dna-dssr** executable, obtainable from <http://x3dna.org>_

Instead of installing **x3dna-dssr** and the plugin locally, you can also use the web server at <http://skmatic.x3dna.org/>

See also the blog post by DSSR author Xiang-Jun Lu: <http://x3dna.org/highlights/dssr-base-blocks-in-pymol-interactively>

## Contents

  * 1 Recommended Setup
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 See Also



## Recommended Setup

Place the **x3dna-dssr** executable into **$PATH**. Examples: 

  * Linux/Mac: **/usr/bin/x3dna-dssr**
  * Windows: **C:\Program Files\PyMOL\PyMOL\x3dna-dssr.exe**



The script can be [run](/index.php/Run "Run"), imported as a Python module, or installed with the [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager"). 

## Usage
    
    
    dssr_block [ selection [, state [, block_file [, block_depth
           [, block_color [, name [, exe ]]]]]]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **state** = int: object state (0 for all states) {default: -1, current state}
  * **block_file** = face|edge|wc|equal|minor|gray: Corresponds to the **\--block-file** option (see [DSSR manual](http://x3dna.org/)). Values can be combined, e.g. "wc-minor". {default: face}
  * **block_depth** = float: thickness of rectangular blocks {default: 0.5}
  * **block_color** = str: Corresponds to the **\--block-color** option (new in DSSR v1.5.2) {default: }
  * **name** = str: name of new CGO object {default: dssr_block##}
  * **exe** = str: path to "x3dna-dssr" executable {default: x3dna-dssr}



## Examples

Combining DSSR block representation with regular PyMOL cartoons: 
    
    
    fetch 1ehz, async=0
    as cartoon
    set cartoon_ladder_radius, 0.1
    set cartoon_ladder_color, gray
    set cartoon_nucleic_acid_mode, 1
    dssr_block
    

Joined base-pair blocks (block_file=wc): 
    
    
    fetch 1ehz, async=0
    dssr_block block_file=wc
    

Multi-state Example: 
    
    
    fetch 2n2d, async=0
    dssr_block 2n2d, 0
    set all_states
    

Custom coloring: 
    
    
    fetch 1msy, async=0
    dssr_block block_color=N red | minor 0.9 | major yellow
    

## See Also

  * [3DNA](/index.php/3DNA "3DNA") (old, obsoleted by x3dna-dssr)
  * [Overview of nucleic acid cartoons](/index.php/Overview_of_nucleic_acid_cartoons "Overview of nucleic acid cartoons")
  * <http://skmatic.x3dna.org/>



Retrieved from "[https://pymolwiki.org/index.php?title=Dssr_block&oldid=13940](https://pymolwiki.org/index.php?title=Dssr_block&oldid=13940)"


---

## Dynamic mesh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/dynamic_mesh.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/dynamic_mesh.py)  
Author(s)  | [Takanori Nakane](/index.php?title=User:TakanoriNakane&action=edit&redlink=1 "User:TakanoriNakane \(page does not exist\)")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**dynamic_mesh** displays [isomesh](/index.php/Isomesh "Isomesh") around the center of the view. When the view is moved, the isomesh will be updated automatically. You can also change contour level by PageDown/PageUp keys. This script is intended to implement interface similar to Coot for examing electron density maps. 

Note: PyMOL's [Density Wizard](/index.php?title=Density_Wizard&action=edit&redlink=1 "Density Wizard \(page does not exist\)") (Menu > Wizard > Density) provides similar functionality. It is implemented using [wizard](/index.php/Wizard "Wizard") framework, while this uses [CallBack](/index.php?title=CallBack&action=edit&redlink=1 "CallBack \(page does not exist\)") object. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    dynamic_mesh map_name [, level [, radius [, name [ sym_source ]]]]
    

## Arguments

  * map_name = string: name of volumetric object(map) to display
  * level = float: contour level of isomesh {default: 1.0}
  * radius = float: radius of isomesh around the center of the view {default: 8}
  * name = string: name of mesh object {default: dynamic_mesh}
  * sym_source = string: name of object from which symmetry information is derived {default: map_name}



## Example
    
    
    fetch 1HWK, async=1
    fetch 1HWK, 1hwk_map, type=2fofc, async=1
    run dynamic_mesh.py
    dynamic_mesh 1hwk_map, sym_source=1hwk
    show sticks, resn 117
    show ribbon
    zoom chain A and resn 117
    

Note: On PyMOL <= 1.4, you have to download the electron density map from the Uppsala Electron Density Server manually. 

## See Also

  * [isomesh](/index.php/Isomesh "Isomesh")
  * [get_position](/index.php/Get_position "Get position")
  * [Density Wizard](/index.php?title=Density_Wizard&action=edit&redlink=1 "Density Wizard \(page does not exist\)")
  * [map_auto_expand_sym](/index.php/Map_auto_expand_sym "Map auto expand sym")
  * Coot <http://lmb.bioch.ox.ac.uk/coot/>



Retrieved from "[https://pymolwiki.org/index.php?title=Dynamic_mesh&oldid=13903](https://pymolwiki.org/index.php?title=Dynamic_mesh&oldid=13903)"


---

## DynoPlot

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/dynoplot.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/dynoplot.py)  
Author(s)  | [Dan Kulp](/index.php/User:Tmwsiy "User:Tmwsiy")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 IMAGES
    * 1.2 SETUP
    * 1.3 NOTES / STATUS
    * 1.4 USAGE
    * 1.5 EXAMPLES



## Introduction

This script was setup to do generic plotting, that is given a set of data and axis labels it would create a plot. Initially, I had it setup to draw the plot directly in the PyMol window (allowing for both 2D and 3D style plots), but because I couldn't figure out how to billboard CGO objects (Warren told me at the time that it couldn't be done) I took a different approach. The plot now exists in it's own window and can only do 2D plots. It is however interactive. I only have here a Rama.(phi,psi) plot, but the code can be easily extended to other types of data. For instance, I had this working for an energy vs distance data that I had generated by another script. 

This script will create a Phi vs Psi(Ramachandran) plot of the selection given. The plot will display data points which can be dragged around Phi,Psi space with the corresponding residue's Phi,Psi angles changing in the structure (PyMol window). 

### IMAGES

  * [![Initial Ramachandran plot of 1ENV](/images/e/e0/RamaPlotInitComposite.png)](/index.php/File:RamaPlotInitComposite.png "Initial Ramachandran plot of 1ENV")

Initial Ramachandran plot of 1ENV 

  * [![Modified pose and plot of 1ENV](/images/c/c7/RamaPlotBentComposite.png)](/index.php/File:RamaPlotBentComposite.png "Modified pose and plot of 1ENV")

Modified pose and plot of 1ENV 




### SETUP

Install from the plugins menu with _Plugin > Manage Plugins > Install ..._ or just [run](/index.php/Run "Run") the script. 

### NOTES / STATUS

  * Tested on Linux, PyMol version 1.4
  * Left, Right mouse buttons do different things; Right = identify data point, Left = drag data point around
  * Post comments/questions or send them to: dwkulp@mail.med.upenn.edu



### USAGE
    
    
    rama [ sel [, name [, symbols [, filename ]]]]
    

### EXAMPLES
    
    
    fetch 1ENV, async=0 # (download it or use the PDB loader plugin)
    select sel01, resi 129-136
    rama sel01
    rock   # the object needs to be moving in order for the angles to be updated.
    

Don't create callback object, use symbols by secondary structure and dump canvas as postscript file: 
    
    
    fetch 2x19, async=0
    color yellow, chain A
    color forest, chain B
    rama polymer, none, ss, /tmp/canvasdump.ps
    rama ss H,    none, aa, /tmp/canvasdump_helix.ps
    rama ss S,    none, aa, /tmp/canvasdump_sheet.ps
    

Retrieved from "[https://pymolwiki.org/index.php?title=DynoPlot&oldid=10435](https://pymolwiki.org/index.php?title=DynoPlot&oldid=10435)"


---

## Elbow angle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/elbow_angle.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/elbow_angle.py)  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  


## Contents

  * 1 Introduction
  * 2 Syntax
  * 3 Examples
    * 3.1 Basic usage
    * 3.2 Included Example
  * 4 Implementation



## Introduction

This script allows you to calculate the elbow angle of an antibody Fab fragment object and optionally draw a graphical representation of the vectors used to calculate the elbow angle. 

  


  


## Syntax
    
    
    elbow_angle object [, light=L [, heavy=H [, limit_l=107 [, limit_h=113 [, draw=0 ]]]]]
    

The first argument should be a PyMOL object. You can calculate the elbow angles of multiple Fabs from the same object (one at a time) by specifying the chains manually for each Fab. 

The `light` and `heavy` arguments are for PDB chain IDs, and `limit_l` and `limit_h` are the last residue of the light and heavy chain variable domains, respectively. (For Kabat-numbered structures, these limits will be 107 and 113, respectively. For structures with different numbering schemes, the limits can be estimated by visual inspection of the PDB file.) Setting `draw=1` will draw the "dumbbells" seen in the images below. 

## Examples

### Basic usage
    
    
    # load an antibody Fab fragment from the PDB
    fetch 3ghe, async=0
    # calculate the elbow angle and draw the vectors
    elbow_angle 3ghe, draw=1
    

This will produce something like the first image below. 

Note that if you don't specify the light/heavy chain IDs or the variable domain limit residue numbers, the default values (L/H and 107/113, respectively) will be used. If your antibody is not Kabat or Chothia numbered, or has different chain names, this will result in an incorrect value or, in the case of wrong chain IDs, may cause the script to fail entirely due to an empty selection. Have a look at Dr. Andrew Martin's [Abnum](http://www.bioinf.org.uk/abs/abnum/) for more information on antibody numbering. 

  * [![Fab fragment 3ghe shown with draw=1.](/images/1/11/Elbow_angle_3ghe.png)](/index.php/File:Elbow_angle_3ghe.png "Fab fragment 3ghe shown with draw=1.")

Fab fragment 3ghe shown with draw=1. 

  * [![5 PDB examples from Stanfield, et al., JMB 2006, shown in the same orientations as in Figure 1 of that paper.](/images/9/95/Stanfield.png)](/index.php/File:Stanfield.png "5 PDB examples from Stanfield, et al., JMB 2006, shown in the same orientations as in Figure 1 of that paper.")

5 PDB examples from Stanfield, et al., JMB 2006, shown in the same orientations as in Figure 1 of that paper. 




The black "dumbbells" pass through the centers of mass of the variable and constant domains of each Fab (as determined using [com.py](/index.php/Com "Com")). The green and red dumbbells denote the residues used to split the variable and constant domains, with a green ball for the light chain, and a red ball for the heavy chain. 

### Included Example

There is also an example .pml file in the Pymol-script-repo/examples directory which can be run by the following: 
    
    
    import ex
    ex elbow_angle
    

This will run the following .pml file: 

Type  | [PyMOL Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [examples/elbow_angle.pml](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/examples/elbow_angle.pml)  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/elbow_angle.pml>" highlight="python" />

and will produce the second image above. 

## Implementation

The elbow angle is the angle between the pseudo-twofold axes between the light and heavy chain variable and constant domains, respectively. The rotation matrix to superpose VL onto VH and CL onto CH are calculated and the vectors corresponding to the rotation axes are determined (using Christoph Gohlke's [transformations.py](/index.php/Transformations "Transformations")). The elbow angle is the obtuse angle obtained by taking the arccos of the dot product of the two vectors. For consistency, the elbow angle is reported as less than 180° when the cross product of the Variable and Constant domain rotation axis vectors ( **V** ⨯ **C** ) is a vector pointing the same direction as that from the limit_h residue alpha carbon atom to the limit_l alpha carbon atom. 

Retrieved from "[https://pymolwiki.org/index.php?title=Elbow_angle&oldid=13905](https://pymolwiki.org/index.php?title=Elbow_angle&oldid=13905)"


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

## Ex

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ex.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ex.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Introduction

A little script, to execute "@" example scripts from **Pymol-script-repo/examples** , no matter what directory you are in. 

### Example of use
    
    
    import ex
    
    ex propka_1
    ex rotkit_1
    ex colorbydisplacement_1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ex&oldid=13906](https://pymolwiki.org/index.php?title=Ex&oldid=13906)"


---

## Extra fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

extra_fit aligns multiple objects to a reference object. 

_New in PyMOL 1.7.2_

It can use any of PyMOL's pairwise alignment methods ([align](/index.php/Align "Align"), [super](/index.php/Super "Super"), [cealign](/index.php/Cealign "Cealign"), [fit](/index.php/Fit "Fit")...). More precisely it can use any function which takes arguments **mobile** and **target** , so it will for example also work with [tmalign](/index.php/Tmalign "Tmalign"). Additional keyword arguments are passed to the used method, so you can for example adjust outlier cutoff or create an alignment object. 

There are several similar commands/scripts for this job, like "A > align > all to this" from the PyMOL panel, the "[alignto](/index.php?title=Alignto&action=edit&redlink=1 "Alignto \(page does not exist\)")" command, [align_all.py and super_all.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/) scripts from Robert Campbell. 

## Usage
    
    
    extra_fit [ selection [, reference [, method ]]]
    

## Example

This will align 4 structures on CA atoms using the [super](/index.php/Super "Super") method. It will also create an alignment object. 
    
    
    fetch 1e4y 1ake 4ake 3hpq, async=0
    remove not chain A
    
    extra_fit name CA, 1ake, super, object=aln_super
    

Same, but with tmalign method (see [TMalign](/index.php/TMalign "TMalign")) 
    
    
    import tmalign
    extra_fit name CA, 1ake, tmalign, object=aln_tmalign
    

## See Also

  * [align_all.py and super_all.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/)



Retrieved from "[https://pymolwiki.org/index.php?title=Extra_fit&oldid=12909](https://pymolwiki.org/index.php?title=Extra_fit&oldid=12909)"


---

## Findseq

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/findseq.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findseq.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview & Motivation

Anyone ever give you a protein and then say, find the sequence "FLVEW"? Well, this script will find any string or regular expression in a given object and return that selection for you. Here's an example, 
    
    
    reinitialize
    import findseq
    # fetch two sugar-binding PDB
    fetch 1tvn, async=0
    # Now, find FATEW in 1tvn, similarly
    findseq FATEW, 1tvn
    # lower-case works, too
    findseq fatew, 1tvn
    # how about a regular expression?
    findseq F.*W, 1tvn
     
    # Find the regular expression:
    #  ..H[TA]LVWH
    # in the few proteins loaded.
    # I then showed them as sticks and colored them to highlight matched AAs
    for x in cmd.get_names(): findseq.findseq("..H[TA]LVWH", x, "sele_"+x, firstOnly=1)
    

[![](/images/c/c3/SeqFinder.png)](/index.php/File:SeqFinder.png)

[](/index.php/File:SeqFinder.png "Enlarge")

Red residues were those matching the regular expression '..H[TA]LVWH'.

# Usage

I built this to be rather flexible. You call it as: 
    
    
    findseq needle, haystack[, selName[, het[, firstOnly ]]]
    

where the options are: 

    

    **needle** the sequence of amino acids to find. Should be a string of one letter amino acid abbreviations. Can also be a string-style regular expression (eg. FW.*QQ).
    **haystack** the PyMOL object or selection in which to search
    **selName** the name of the returned selection. If you leave this blank, it'll be _foundSeqXYZ_ where XYZ is some random integer (eg. foundSeq1435); if you supply _sele_ then the usual PyMOL _(sele)_ is used; and, finally, if it's anything else, then that will be used verbatim. Defaults to _foundSeqXYZ_ so as not to overwrite any selections you might have in _sele_.
    **het** 0/1 -- if 0 then heteroatoms are not considered; if 1 then they are; defaults to 0
    **firstOnly** 0/1 -- if 0 then all matches are selected and returned; if 1 then only the first is returned

# See Also

select_pepseq@[Psico](/index.php/Psico "Psico")

Retrieved from "[https://pymolwiki.org/index.php?title=Findseq&oldid=13907](https://pymolwiki.org/index.php?title=Findseq&oldid=13907)"


---

## FindSurfaceCharge

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/findSurfaceCharge.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findSurfaceCharge.py)  
Author(s)  | [Teddy Warner](/index.php/User:TeddyWarner "User:TeddyWarner")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Drawing upon the [findSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") script, the findSurfaceCharge script will identify and output a list of all charged residues on the surface of a selectionand calculates the ionization state of a protein at a given pH. The charge can be calculated for either a folded or denatured protein. This function is intended to be used to give buffer conditions for mass spectrometry. 

# Usage
    
    
       findSurfaceCharge [selection, [pH, [folded ,[cutoff]]]]
    

# Arguments

  * **selection** = str: The object or selection for which to find exposed residues {default: all}
  * **pH** = float: The pH to calculate the surface charge at {default: 7.0}
  * **folded** = bool: Whether the program should calculate the charge of a folded (True) or denatured (False) protein.
  * **cutoff** = float: The cutoff in square Angstroms that defines exposed or not. Those atoms with > cutoff Ang^2 exposed will be considered _exposed_ {default: 2.5 Ang^2}



# Examples

  * [![Result of 4FIX.pdb at pH 7.0.](/images/2/2a/SurfaceCharge.PNG)](/index.php/File:SurfaceCharge.PNG "Result of 4FIX.pdb at pH 7.0.")

Result of 4FIX.pdb at pH 7.0. 



    
    
    run findSurfaceResiduesListCharged.py
    fetch 4FIX
    
    findSurfaceResiduesListCharged
    
    # see how pH changes the protein surface charge:
    findSurfaceCharge("4fix",7.0)
        Exposed charged residues: 
            ERREDRKEE...
        The expected surface charge of 4fix at pH 7.0 is: +3.24
    
    findSurfaceCharge("4fix",7.0,False)
        Charged residues: 
            HHHHHHRHERREDRKEE...
        The expected charge of denatured 4fix at pH 7.0 is: +0.13
    
    findSurfaceCharge("4fix",10.0)
        Charged residues: ...
        The expected surface charge of 4fix at pH 10.0 is: -3.86
    

# See Also

  * [findSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [get_area](/index.php/Get_Area "Get Area")



Retrieved from "[https://pymolwiki.org/index.php?title=FindSurfaceCharge&oldid=13951](https://pymolwiki.org/index.php?title=FindSurfaceCharge&oldid=13951)"


---

## FindSurfaceResidues

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/findSurfaceResidues.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findSurfaceResidues.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The findSurfaceResidues script will select (and color if requested) surface residues and atoms on an object or selection. See the options below. 

Each time, the script will create two new selections called, **exposed_res_XYZ** and **exposed_atm_XYZ** where _XYZ_ is some random number. This is done so that no other selections/objects are overwritten. 

# Usage
    
    
    findSurfaceResidues [ selection=all [, cutoff=2.5 [, doShow=0 ]]]
    

# Arguments

  * **selection** = str: The object or selection for which to find exposed residues {default: all}
  * **cutoff** = float: The cutoff in square Angstroms that defines exposed or not. Those atoms with > cutoff Ang^2 exposed will be considered _exposed_ {default: 2.5 Ang^2}
  * **doShow** = 0/1: Change the visualization to highlight the exposed residues vs interior {default: 0}



# Examples

  * [![Result of $TUT/1hpv.pdb at 2.5 Ang cutoff.](/images/f/fe/FindExRes.png)](/index.php/File:FindExRes.png "Result of $TUT/1hpv.pdb at 2.5 Ang cutoff.")

Result of $TUT/1hpv.pdb at 2.5 Ang cutoff. 

  * [![Example coloring of surface residues](/images/8/83/Surface_residues_ex.png)](/index.php/File:Surface_residues_ex.png "Example coloring of surface residues")

Example coloring of surface residues 



    
    
    run findSurfaceResidues.py
    fetch 1hpv, async=0
    
    findSurfaceResidues
    
    # now show the exposed
    findSurfaceResidues doShow=1
    
    # watch how the visualization changes:
    findSurfaceResidues doShow=1, cutoff=0.5
    findSurfaceResidues doShow=1, cutoff=1.0
    findSurfaceResidues doShow=1, cutoff=1.5
    findSurfaceResidues doShow=1, cutoff=2.0
    findSurfaceResidues doShow=1, cutoff=2.5
    findSurfaceResidues doShow=1, cutoff=3.0
    

# See Also

  * [get_sasa_relative](/index.php/Get_sasa_relative "Get sasa relative")
  * [get_area](/index.php/Get_Area "Get Area")



Retrieved from "[https://pymolwiki.org/index.php?title=FindSurfaceResidues&oldid=13953](https://pymolwiki.org/index.php?title=FindSurfaceResidues&oldid=13953)"


---

## Flatten obj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/flatten_obj.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/flatten_obj.py)  
Author(s)  | [Spencer Bliven](/index.php/User:Sbliven "User:Sbliven")  
License  | Public Domain   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Arguments
  * 4 Notes
  * 5 Examples
  * 6 See Also



# Description

The `flatten_obj` python script combines multiple objects or states into a single object, renaming chains where required. 

This is particularly useful for dealing with biological assemblies, which are loaded as multi-state objects when fetched using `[fetch](/index.php/Fetch "Fetch") PDBID, type=pdb1`. It can also be used as a quick way to combine multiple objects without causing collisions between chain identifiers. 

The command re-letters chains to avoid collisions. Older versions of PyMOL restrict the chain id to a single character, so the script will fail for assemblies with >62 chains. With more recent versions, this problem is solved with multi-character chain IDs. Several options are available for how re-lettering should occur. 

# Usage
    
    
       flatten_obj name, selection[, state[, rename[, quiet[, chain_map]]]]
    

# Arguments

  * name = a unique name for the flattened object {default: flat}


  * selection = the set of objects to include in the flattening. The selection will be expanded to include all atoms of objects. {default: all}


  * state = the source state to select. Use 0 or -1 to flatten all states {default: 0}


  * rename = The scheme to use for renaming chains: {default: 0} 
    * (0) preserve chains IDs where possible, rename other chains alphabetically
    * (1) rename all chains alphabetically
    * (2) rename chains using the original chain letter, object name, and state


  * quiet = If set to 0, print some additional information about progress and chain renaming {default: 1}


  * chain_map = An attribute name for the 'stored' scratch object. If specified, `stored.<chain_map>` will be populated with a dictionary mapping the new chain names to a tuple giving the originated object, state, and chainID. {default: ""}



# Notes

Like the select command, if name is omitted then the default object name ("flat") is used as the name argument. 

Chain renaming is tricky. PDB files originally limited chains to single letter identifiers containing [A-Za-z0-9]. When this was found to be limiting, multi-letter chains (ideally < 4 chars) were allowed. This is supported as of PyMOL 1.7. Earlier versions do not accept rename=2, and will raise an exception when flattening a structure with more than 62 chains. 

# Examples
    
    
       flatten_obj flat, nmrObj
       flatten_obj ( obj1 or obj2 )
    

# See Also

  * [split_states](/index.php/Split_states "Split states")



Retrieved from "[https://pymolwiki.org/index.php?title=Flatten_obj&oldid=13954](https://pymolwiki.org/index.php?title=Flatten_obj&oldid=13954)"


---

## Format bonds

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/format_bonds.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/format_bonds.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The script **format_bonds** will automatically format bonds in amino acids. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Notes
  * 4 SEE ALSO



## Usage
    
    
    format_bonds [ selection [, bonds ]]
    

## Examples

  * [![format_bonds bonds=1](/images/1/10/PHE_valence_0.png)](/index.php/File:PHE_valence_0.png "format_bonds bonds=1")

format_bonds bonds=1 

  * [![format_bonds bonds=2](/images/9/9f/PHE_valence_1_mode_1.png)](/index.php/File:PHE_valence_1_mode_1.png "format_bonds bonds=2")

format_bonds bonds=2 

  * [![format_bonds](/images/c/ce/PHE_delocalized.png)](/index.php/File:PHE_delocalized.png "format_bonds")

format_bonds 



    
    
    import format_bonds
    
    frag PHE
    format_bonds
    
    format_bonds bonds=2
    

  


## Notes

  * Remember to correctly configure plugin import (see: [Git intro](/index.php/Git_intro "Git intro"))
  * **format_bonds** will introduce delocalized bonds by default or when _bonds_ is larger than 2.  

  * Setting _bonds=1_ will simply disable valence display (globally!)
  * The _selection_ argument is 'all' by default and can be used to restrict editing to selected residues.  

  * Note that **format_bonds** will also format acidic residues, the C-terminus, arginine and nitro groups.
  * Tip: press the TAB key after entering format_bonds to get argument suggestions



  


# SEE ALSO

[Git intro](/index.php/Git_intro "Git intro"), [Valence](/index.php/Valence "Valence"), [Modeling and Editing Structures](/index.php/Modeling_and_Editing_Structures "Modeling and Editing Structures")

Retrieved from "[https://pymolwiki.org/index.php?title=Format_bonds&oldid=13941](https://pymolwiki.org/index.php?title=Format_bonds&oldid=13941)"


---

## Forster distance calculator

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/forster_distance_calculator.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/forster_distance_calculator.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Spectre input
  * 3 Getting spectres
  * 4 Input values
  * 5 How the script works
  * 6 How to run the script



## Introduction

This script can be handsome, if one is working with Förster resonance energy transfer in proteins. The script calculates the [Förster distance](http://en.wikipedia.org/wiki/F%C3%B6rster_resonance_energy_transfer): R0, from two dyes excitation/emission spectres.   
This script is very handsome, if one want's to pick two dyes from different companies, and the Förster distance is not provided in a table. 

This script does no calculation of proteins in pymol, but is made as a "hack/shortcut" for a python script.   
We use the python part of pymol to do the calculations, so a student would not need to install python at home, but simply pymol.   


## Spectre input

**There should be provided the path to four datafiles:**   
Excitation and Emission spectre for the Donor dye: D_Exi="path.txt" D_Emi="path.txt   
Excitation and Emission spectre for the Acceptor dye: A_Exi="path.txt" A_Emi="path.txt   


Each of the files shall be a two column file. Separated by space. Numbers are "." dot noted and not by "," comma. One can for example prepare the files by search-and-replace in a small text editor.   
The first column is the wavelength in nanometers "nm" or centimetres "cm". Second column is the arbitrary units of excitation/emission. **For example:**
    
    
    Absorption
    Wavelength "Alexa Fluor 488 H2O (AB)" 
    
    300.00 0.1949100000 
    301.00 0.1991200000 
    302.00 0.2045100000 
    303.00 0.2078800000 
    ....
    

## Getting spectres

For example provides the company [ATTO-TEC](http://www.atto-tec.com/index.php?id=128&L=1&language=en) spectre for their dyes in excel format [Spectre(excel-file)](https://www.atto-tec.com/index.php?id=113&no_cache=1&L=1). These can easily be copied from excel to a flat two column file. 

Some of the most cited dyes are the Alexa Fluor dyes from the company [invitrogen](http://www.invitrogen.com/site/us/en/home/References/Molecular-Probes-The-Handbook/Fluorophores-and-Their-Amine-Reactive-Derivatives/Alexa-Fluor-Dyes-Spanning-the-Visible-and-Infrared-Spectrum.html). Invitrogen does not provide spectres for their dyes in dataformat, but as flat picture files. 

Luckily, a group on University of Arizona have traced several spectre of dyes from the literature with a graphics program. They have made these spectre easily public at <http://www.spectra.arizona.edu/>. I highly recommend this homepage. With these Spectra, the script can calculate the Forster Distance for different dyes from different companies.   


Download one spectrum at the time by "deselecting" in the right side of the graph window. Then get the datafile with the button in the lower right corner. 

## Input values

The calculation needs input for: 

  * **Qd** : Fluorescence quantum yield of the donor in the absence of the acceptor. 

    This is normally provided from the manufacture, and this information is normally also available in the [database](http://www.spectra.arizona.edu/).
    It is recognized in the literature, that the value for **Qd** , changes on behalf on which position its located on the protein.
    Changes in 5-10 % is normal, but much larger values can be expected if the dye make hydrophobic interactions with the protein. This is dye-nature dependant.
    This value can only be determined precisely with experiments for the protein at hand.
    But anyway, for a "starting" selection for a FRET pair, the manufacture **Qd** should be "ok".
  * **Kappa2 = 2.0/3.0** : Dipole orientation factor. 

    In the literature, I often see this equal 2/3=0,667 This corresponds for free diffusing dye.
    For a dye attached to a protein, and restricted in one direction, the value is equal Kappa2=0.476
    But since we take the 6th root of the value, the difference ends up only being 5% in relative difference for the R0 value.
    Kappa2=2/3 can be a valid fast approximation, for the purpose of ordering dyes. But be careful and check the literature for discussions.
  * **n = 1.33** : 

    water=1.33, protein=1.4, n(2M GuHCl)=1.375
    Again, we take the sixth root, so the differences are "not so important".



## How the script works

  * The script integrate the Donor emission spectre, to divide with it. 

    This is done by simple numerical integration by the [Rectangle method](http://en.wikipedia.org/wiki/Rectangle_method)
  * All datapoints for the Acceptor excitation spectre are scaled, so that **e**(molar extinction coefficient) fits the wavelength where it is determined.
  * For the multiplication of the spectre, the script test that both data points exist before multiplication. 

    This can be troublesome, if there is 1-2 nm between each data point. Or if they are "un-even". Like Donor has 100.2 100.4 100.4 and Acceptor has 100.1 100.3 100.5
    But most spectre have much better resolution. Mostly 0.2 nm, and "even spaced"
    So this should work quite good as well.



The output is a text file, with all the datapoints. 

  1. Column: Donor: The input Emission wavelength.
  2. Column: Donor: The input Emission, arbitrary units of excitation.
  3. Column: Donor: The input Emission, arbitrary units of excitation, divided by the integrated area.
  4. Column: Acceptor: The input excitation wavelength.
  5. Column: Acceptor: The input excitation **e**(molar extinction coefficient).
  6. Column: Acceptor: The input excitation **e**(molar extinction coefficient), scaled correctly.
  7. Column: The calculated overlap for this datapoint.
  8. Column: The summed values for calculated overlap, until this datapoint.



  
It also automatically generates a [gnuplot](http://www.gnuplot.info/) .plt file, for super fast making graphs of the spectres. In linux, gnuplot is normally part of the installation. Just write: gnuplot FILENAME.plt If you are in windows, [download gnuplot](http://sourceforge.net/projects/gnuplot/files/) and open the .plt file. 

Gnuplot makes three graphs in .eps format. They can be converted to pdf with the program: epstopdf or eps2pdf. They are for example part of LaTeX: C:\Program Files (x86)\MiKTeX 2.9\miktex\bin or you can [download it here](http://tinyurl.com/eps2pdf). The format of .eps is choses, since gnuplot then allows for math symbols in the graphs. 

  * [![All graphs plotted together, but not scaled correctly.](/images/7/77/1-ALEXA488Emi-ALEXA633Exi-overlap-all-spectre.png)](/index.php/File:1-ALEXA488Emi-ALEXA633Exi-overlap-all-spectre.png "All graphs plotted together, but not scaled correctly.")

All graphs plotted together, but not scaled correctly. 

  * [![The Donor emission \(y1\) and Acceptor excitation \(y2\) scaled correctly.](/images/7/7d/2-ALEXA488Emi-ALEXA633Exi-overlap-normalized-spectre.png)](/index.php/File:2-ALEXA488Emi-ALEXA633Exi-overlap-normalized-spectre.png "The Donor emission \(y1\) and Acceptor excitation \(y2\) scaled correctly.")

The Donor emission (y1) and Acceptor excitation (y2) scaled correctly. 

  * [![The Rectangle integral datapoint \(y1\) and the summed integration \(y2\).](/images/0/09/3-ALEXA488Emi-ALEXA633Exi-overlap-integral.png)](/index.php/File:3-ALEXA488Emi-ALEXA633Exi-overlap-integral.png "The Rectangle integral datapoint \(y1\) and the summed integration \(y2\).")

The Rectangle integral datapoint (y1) and the summed integration (y2). 




## How to run the script

Make a pymol .pml file like this. Write in the required information, and execute/run the script with pymol. Then open the .plt file afterwards with gnuplot. 
    
    
    ## Change to your directory
    cd /homes/YOU/Documents/Atto-dyes/Spectre/ALEXA488-ALEXA633
    import forster_distance_calculator
    forster D_Exi=ALEXA488Exi.txt, D_Emi=ALEXA488Emi.txt, A_Exi=ALEXA633Exi.txt, A_Emi=ALEXA633Emi.txt, A_e_Max_Y=159000, A_e_Max_X=621, Qd=0.92
    

Retrieved from "[https://pymolwiki.org/index.php?title=Forster_distance_calculator&oldid=13908](https://pymolwiki.org/index.php?title=Forster_distance_calculator&oldid=13908)"


---

## Frame slider

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/frame_slider.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/frame_slider.py)  
Author(s)  | Matthew Baumgartner   
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  
The Frame slider plugin provides a simple slider bar that allows you to quickly skip through frames in PyMOL. It also has a text field that you can type the desired frame in and it will automatically go to that frame. 

[![Frame slider.png](/images/a/ad/Frame_slider.png)](/index.php/File:Frame_slider.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Frame_slider&oldid=13909](https://pymolwiki.org/index.php?title=Frame_slider&oldid=13909)"


---

## Get colors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/get_colors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/get_colors.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  * **get_colors** contains two functions that can be useful when working with colors 
    * _get_colors_ : returns all available PyMOL colors
    * _get_random_color_ : returns a random available PyMOL color


  * Note that basic colors can be accessed manually without this script from the PyMOL menu under **Setting** \--> **Colors...**
  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



[![1LSD colored randomly by residue](/images/1/13/1LSD_random_colors.png)](/index.php/File:1LSD_random_colors.png "1LSD colored randomly by residue")

  


## Usage
    
    
    get_colors [ selection [, quiet ]]
    get_random_color [ selection [, quiet ]]
    

## Examples
    
    
    # basic example
    get_colors # basic colors
    get colors all # larger range with intermediates
    
    
    
    #get disco funky
    import get_colors
    from get_colors import get_colors
    from get_colors import get_random_color
    
    cmd.delete('all')
    cmd.fetch('1LSD', async=0) # :P
    cmd.hide('everything')
    cmd.show_as('sticks','not hetatm')
    cmd.orient()
    
    python # start a python block
    from pymol import stored
    stored.atom_list=[]
    cmd.iterate('name CA','stored.atom_list.append([model, resi])')
    resi_list=["model %s and resi %s"%(value[0],value[1]) for value in stored.atom_list]
    for resi in resi_list: cmd.color(get_random_color(),resi)
    python end # end python block
    cmd.ray()
    

|   
---|---  
  
  


## SEE ALSO

  * [Color](/index.php/Color "Color")
  * [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_colors&oldid=13942](https://pymolwiki.org/index.php?title=Get_colors&oldid=13942)"


---

## Get raw distances

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/get_raw_distances.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/get_raw_distances.py)  
Author(s)  | [Takanori Nakane](/index.php?title=User:TakanoriNakane&action=edit&redlink=1 "User:TakanoriNakane \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.querying](http://pymol.org/psicoredirect.php?psico.querying)  
  
get_raw_distances can dump distance objects, created with [distance](/index.php/Distance "Distance"). 

This script also provides the command **select_distances** , which selects atoms from distance objects. 

_Warning: the atoms are hashed by coordinates; this could cause issues if coordinates are altered after distance objects have been created (see also[dynamic_measures](/index.php?title=Dynamic_measures&action=edit&redlink=1 "Dynamic measures \(page does not exist\)") setting)._

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Example Atom Selection
  * 5 See Also



## Usage
    
    
    get_raw_distances [ names [, state [, selection ]]]
    
    
    
    select_distances [ names [, name [, state [, selection [, cutoff ]]]]]
    

## Arguments

  * **names** = string: names of distance objects (no wildcards!) {default: all measurement objects}
  * **state** = integer: object state {default: 1}
  * **selection** = string: atom selection {default: all}



## Example
    
    
    fetch 2xwu, async=0
    
    # interface polar contacts
    distance iface_hbonds, chain A, chain B, mode=2
    
    # dump (model,index) information
    get_raw_distances iface_hbonds
    

## Example Atom Selection
    
    
    # select atoms
    select_distances iface_hbonds, hbsele
    
    # nice representation
    set cartoon_side_chain_helper
    set sphere_scale, 0.5
    as cartoon, 2xwu
    show sticks, byres hbsele
    show spheres, hbsele
    

## See Also

  * [distance](/index.php/Distance "Distance")
  * [dynamic_measures](/index.php?title=Dynamic_measures&action=edit&redlink=1 "Dynamic measures \(page does not exist\)")
  * [find_pairs](/index.php/Find_pairs "Find pairs")
  * [get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_raw_distances&oldid=13910](https://pymolwiki.org/index.php?title=Get_raw_distances&oldid=13910)"


---

## Grepset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/grepset.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/grepset.py)  
Author(s)  | [Ezequiel Panepucci](/index.php?title=User:Zac&action=edit&redlink=1 "User:Zac \(page does not exist\)")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.helping](http://pymol.org/psicoredirect.php?psico.helping)  
  
Use this little script to explore PyMOL's myriad settings. 

Usefull for newbies and those with not so good memory skills... 

To use: 

  1. put the script in a file called **grepset.py**
  2. from within PyMOL execute **run grepset.py**. If you have the Pymol-script library configured, then **import grepset**
  3. try it out, see examples below.



Example 1: **grepset light**
    
    
    cartoon_highlight_color        default
    dot_lighting                   on
    light                          [ -0.40000, -0.40000, -1.00000 ]
    mesh_lighting                  off
    two_sided_lighting             off
    5 settings matched
    

Example 2: **grepset trans**
    
    
    cartoon_transparency           0.00000
    ray_transparency_contrast      1.00000
    ray_transparency_shadows       1
    ray_transparency_spec_cut      0.90000
    ray_transparency_specular      0.40000
    sphere_transparency            0.00000
    stick_transparency             0.00000
    transparency                   0.00000
    transparency_mode              2
    transparency_picking_mode      2
    10 settings matched
    

Example 3: **grepset ^trans**
    
    
    transparency                   0.00000
    transparency_mode              2
    transparency_picking_mode      2
    3 settings matched
    

Retrieved from "[https://pymolwiki.org/index.php?title=Grepset&oldid=13911](https://pymolwiki.org/index.php?title=Grepset&oldid=13911)"


---

## Hbplus

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [scripts/hbplus.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/hbplus.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
hbplus is a PyMOL wrapper for the [HBPLUS](http://www.biochem.ucl.ac.uk/bsm/hbplus/home.html) program. It creates distance objects for the hydrogen bonds found by HBPLUS. 

## Usage
    
    
    hbplus [ selection [, exe [, prefix [, state ]]]]
    

## See Also

  * [distance](/index.php/Distance "Distance")



Retrieved from "[https://pymolwiki.org/index.php?title=Hbplus&oldid=13912](https://pymolwiki.org/index.php?title=Hbplus&oldid=13912)"


---

## Inertia tensor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/inertia_tensor.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/inertia_tensor.py)  
Author(s)  | [Mateusz Maciejewski](http://www.mattmaciejewski.com)  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![It.png](/images/d/dd/It.png)](/index.php/File:It.png)

This script will draw the eigenvectors of the inertia tensor of the selected part of the molecule. 

## Usage
    
    
    tensor selection, name, model
    

## Examples

To draw the inertia tensor of the first model in PDB file 4B2R do: 
    
    
    fetch 4b2r, async=0
    hide lines
    show cartoon
    run inertia_tensor.py
    tensor 4b2r & n. n+ca+c & i. 569-623, tensor_dom11, 1
    

## Reference

A figure generated using this script can be seen in the following reference: 

  * _Estimation of interdomain flexibility of N-terminus of factor H using residual dipolar couplings_. Maciejewski M, Tjandra N, Barlow PN. **Biochemistry**. 50(38) 2011, p. 8138-49, fig. 1 [doi:10.1021/bi200575b](http://dx.doi.org/10.1021/bi200575b)



Retrieved from "[https://pymolwiki.org/index.php?title=Inertia_tensor&oldid=13943](https://pymolwiki.org/index.php?title=Inertia_tensor&oldid=13943)"


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

## Load img stack

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/load_img_stack.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/load_img_stack.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the load_img_stack command, which loads a set of images as a map object. The images can either be individual files, or contained in a multi-page TIFF file. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    load_img_stack pattern [, name [, grid [, channel [, normalize [, extent ]]]]]
    

## Arguments

  * **pattern** = str: image filename or pattern
  * **name** = str: map object name to create
  * **grid** = float: grid spacing in Angstrom {default: 1.0}
  * **channel** = int: color channel for RGB images {default: 0}
  * **normalize** = 0 or 1: normalize data {default: 1}
  * **extent** = 3-float: (a,b,c) edge lengths in Angstrom, overwrites "grid" arguments if given {default: }



## Example
    
    
    load_img_stack img*.png, name=map, extent=(10.0, 10.0, 5.0)
    volume vol, map, rainbow
    

## See Also

  * [tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")
  * [volume](/index.php/Volume "Volume")



Retrieved from "[https://pymolwiki.org/index.php?title=Load_img_stack&oldid=13944](https://pymolwiki.org/index.php?title=Load_img_stack&oldid=13944)"


---

## Modevectors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/modevectors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/modevectors.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
### DESCRIPTION

[![](/images/3/34/Modevectors.png)](/index.php/File:Modevectors.png)

[](/index.php/File:Modevectors.png "Enlarge")

Example showing modevectors in action. (See Examples below).

**modevectors.py** is a PyMol script that was originally written to visualize results obtained from Normal Mode Analysis (NMA) by drawing arrows or vectors from a starting structure to a final structure. However, this script can be used to visualize the direction of motion between two specified states (e.g. targeted MD, comparing open and closed structures, etc). The strength of this script is that the arrows are highly customizable with numerous options to choose from (see script below). It is important to note that the two states MUST be identical except for the X, Y, and Z coordinates. That is, you cannot draw vectors between homologous structures. The default settings sets the background to white and displays the arrows along with the first object frame (in cartoon representation). 

  
Note: Ray tracing these CGO arrows appears to be a memory intensive process so please save a session before ray tracing and use "pymol -qc" where necessary. If anybody can come up with a method to improve this please e-mail me (see address in script). 

Update: A new way of drawing the arrows has been implemented for PyMOL v.1.1 and has been updated in this code (which automatically detects the version that you are using). The new method uses the built in cone primitive rather than concentric circles decreasing in size. Thus, the new method is much less memory intensive and is a significant improvement over the old method. Additionally, you may want to turn off [ray shadow](/index.php/Ray_shadow "Ray shadow") depending on the [scene](/index.php/Scene "Scene") that you are viewing as shadows can be confusing to others when seen in publications. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    modevectors first_obj_frame, last_obj_frame [,first_state=1 [,last_state=1 [,outname=modevectors \
      [,head=1.0 [,tail=0.3 \[,head_length=1.5 [,headrgb=(1.0,1.0,1.0) [,tailrgb=(1.0,1.0,1.0) [,cutoff=4.0 
      [,skip=0 [,cut=0.5 [,atom=CA [,stat=show [,factor=1.0 [,notail=0]]]]]]]]]]]]]]
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "modevectors" (which can be changed through the options). 

  


### EXAMPLES
    
    
    modevectors 1E3M, 1W7A
    modevectors 1E3M, 1W7A, outname="arrows"
    modevectors 1E3M, 1W7A, atom="P"
    

Copy/paste the following code to see an example of modevectors. This uses a multistate protein and the arrows are connected between the first and last states. 
    
    
    import modevectors
    # fetch the PDBs from pdb.org
    fetch 1c3y, finish=1, multiplex=0, async=0
    # separate the first and last states of the NMR ensemble to individual objects
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
    # run the modevectors code
    modevectors 1c3y_0001, 1c3y_0023
    # just setup a nice representation
    as cartoon, 1c3y_0001 or 1c3y_0023
    show cgo, modevectors
    color marine, 1c3y_0001
    color purpleblue, 1c3y_0023
    

The following set of examples will illustrate the power of each optional argument. Each example should be compared to the default figure in the table below. 

[![](/images/3/35/Mv_default.png)](/index.php/File:Mv_default.png) [](/index.php/File:Mv_default.png "Enlarge")Fig.1 - Default Settings | [![](/images/b/b0/Mv_fig2.png)](/index.php/File:Mv_fig2.png) [](/index.php/File:Mv_fig2.png "Enlarge")Fig.2 - Arrow direction drawn from CA to CA | [![](/images/4/4f/Mv_fig3.png)](/index.php/File:Mv_fig3.png) [](/index.php/File:Mv_fig3.png "Enlarge")Fig.3 - Arrow head radius  
---|---|---  
[![](/images/4/4f/Mv_fig4.png)](/index.php/File:Mv_fig4.png) [](/index.php/File:Mv_fig4.png "Enlarge")Fig.4 - Arrow tail radius | [![](/images/f/f1/Mv_fig5.png)](/index.php/File:Mv_fig5.png) [](/index.php/File:Mv_fig5.png "Enlarge")Fig.5 - Arrow head length | [![](/images/e/e4/Mv_fig6.png)](/index.php/File:Mv_fig6.png) [](/index.php/File:Mv_fig6.png "Enlarge")Fig.6 - Arrow head RGB  
[![](/images/f/fb/Mv_fig7.png)](/index.php/File:Mv_fig7.png) [](/index.php/File:Mv_fig7.png "Enlarge")Fig.7 - Arrow tail RGB | [![](/images/b/bc/Mv_fig8.png)](/index.php/File:Mv_fig8.png) [](/index.php/File:Mv_fig8.png "Enlarge")Fig.8 - Arrow length cutoff | [![](/images/7/7f/Mv_fig9.png)](/index.php/File:Mv_fig9.png) [](/index.php/File:Mv_fig9.png "Enlarge")Fig.9 - Skip every other arrow  
[![](/images/d/d7/Mv_fig10.png)](/index.php/File:Mv_fig10.png) [](/index.php/File:Mv_fig10.png "Enlarge")Fig.10 - Subtract value from vector length | [![](/images/4/42/Mv_fig11.png)](/index.php/File:Mv_fig11.png) [](/index.php/File:Mv_fig11.png "Enlarge")Fig.11 - Draw arrows from CB atoms | [![](/images/b/b1/Mv_fig12.png)](/index.php/File:Mv_fig12.png) [](/index.php/File:Mv_fig12.png "Enlarge")Fig.12 - Shorten by 50%  
| [![](/images/a/a3/Mv_fig13.png)](/index.php/File:Mv_fig13.png) [](/index.php/File:Mv_fig13.png "Enlarge")Fig.13 - Multiple options being used (see final example below)  
      
    
    reinitialize
    import modevectors
    
    fetch 1c3y, async=0
    
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
     
    #This is the default setting (Fig.1)
    modevectors 1c3y_0001, 1c3y_0023
     
    #This is shows that the direction of the arrows is drawn from the 1c3y_0001 towards 1c3y_0023 (Fig.2)
    show cartoon, 1c3y_0023
    color red, 1c3y_0023
     
    #This controls the base radius of the cone/arrow head in angstroms (Fig.3)
    modevectors 1c3y_0001, 1c3y_0023, head=2.5
     
    #This controls the radius of the cylinders/arrow tails in angstroms (Fig.4)
    modevectors 1c3y_0001, 1c3y_0023, tail=0.75
     
    #This controls the length of the cone/arrow head in angstroms (Fig.5)
    #Note that this option does NOT increase the vector length but simply changes the tail length
    modevectors 1c3y_0001, 1c3y_0023, head_length=3.0
     
    #This controls the colour of the cone/arrow head using RGB values (Fig.6)
    modevectors 1c3y_0001, 1c3y_0023, headrgb=(1.0,0.0,0.0)
     
    #This controls the colour of the cylinder/arrow tails using RGB values (Fig.7)
    modevectors 1c3y_0001, 1c3y_0023, tailrgb=(1.0,0.0,0.0)
     
    #This controls the which vectors to show based on a specific cutoff value in angstroms.  Vector lengths that are less 
    #than the cutoff value will not be displayed (Fig.8)
    modevectors 1c3y_0001, 1c3y_0023, cutoff=30.0
     
    #This controls how many vectors to skip (integer value) and is useful when there are too many vectors showing.  
    #Skip=1 will show every other arrow (Fig.9)
    modevectors 1c3y_0001, 1c3y_0023, skip=1
     
    #This controls how much to cut from each vector (in angstroms).  Note that arrows will point in the opposite directions 
    #when too much is cutoff (resulting in a negative vector length) (Fig.10) and should be used with caution!
    modevectors 1c3y_0001, 1c3y_0023, cut=15.0
     
    #This controls which atom to draw a vector from (Fig.11).  Note that this is case-sensitive and is really only useful 
    #when atom=CA or when atom=P (for DNA backbones)
    modevectors 1c3y_0001, 1c3y_0023, atom=CB
     
    #This controls how much to multiple the length of each vector by (percentage increase/decrease) (Fig.12)
    #This example halves the length of each vector (50%)
    modevectors 1c3y_0001, 1c3y_0023, factor=0.5
     
    #This hides the statistics which count atoms skipped, atoms counted (number of arrows showing), and number of atoms 
    #that did not meet the cutoff and are not shown
    modevectors 1c3y_0001, 1c3y_0023, stat=hide
     
    #This example shows multiple options being used (Fig.13)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, tail=1.0, head_length=2.0, headrgb=(1.0,0.0,0.0), tailrgb=(0.0,0.0,1.0),
    cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8
     
    #Finally, this example hides all arrow tails and only uses arrow heads via the notail option(No Figure)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8, notail=1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Modevectors&oldid=13914](https://pymolwiki.org/index.php?title=Modevectors&oldid=13914)"


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

## Movie color fade

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/movie_color_fade.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/movie_color_fade.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![movie_color_fade example](/images/0/05/Movie_color_fade_example_1.gif)](/index.php/File:Movie_color_fade_example_1.gif "movie_color_fade example")

**movie_color_fade** is like [movie_fade](/index.php/Movie_fade "Movie fade"), but will fade colors in a movie. Simply specify the arguments (see below). 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Examples
  * 3 SEE ALSO



## Usage
    
    
     movie_color_fade [ startframe [, startcolor [, endframe [, endcolor [, selection ]]]]]
     help movie_color_fade
    

|   
---|---  
  
### Arguments

  * startframe, endframe = beginning and end movie frame for fading 
    * if omitted, the current or last frame will be respectively set automatically
  * startcolor, endcolor = coloring at start and end
  * selection: target selection



## Examples

This example yields the movie to the top right 
    
    
    # make function available to PyMOL
    import movie_color_fade
    
    # object prep
    bg black
    delete all
    fetch 1hpv, async=0
    as cartoon
    orient
    color yellow, chain A
    color skyblue, chain B
    # movie prep
    mset 1x120
    
    # color fading
    movie_color_fade 1, yellow, 60, skyblue, chain A
    movie_color_fade 60, skyblue, 120, yellow, chain A
    movie_color_fade 1, skyblue, 60, yellow, chain B
    movie_color_fade 60, yellow, 120, skyblue, chain B
    

|   
---|---  
  
  


## SEE ALSO

  * [Movie_fade](/index.php/Movie_fade "Movie fade")
  * [mappend](/index.php/Mappend "Mappend")
  * [Color](/index.php/Color "Color")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_color_fade&oldid=13945](https://pymolwiki.org/index.php?title=Movie_color_fade&oldid=13945)"


---

## Movie fade

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/movie_fade.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/movie_fade.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will help fade in and out settings in a movie. Just specify the setting, it's initial value at an initial frame and it's ending value and frame. 

## Usage
    
    
    movie_fade setting, startFrame, startVal, endFrame, endVal [, selection ]
    

## Examples

To fade in sticks from fully transparent to fully opaque across 60 frames do: 
    
    
    mset 1x60
    movie_fade stick_transparency, 1, 1., 60, 0.
    

More complex example which involves camera motion: 
    
    
    fetch 1rx1, async=0
    as cartoon
    show surface
    mset 1x80
    movie.roll
    movie_fade transparency,  1, 0., 40, 1.
    movie_fade transparency, 41, 1., 80, 0.
    

**Example file:**

[![](/images/9/98/Movie_fade_example.gif)](/index.php/File:Movie_fade_example.gif)

[](/index.php/File:Movie_fade_example.gif "Enlarge")

Result
    
    
    import movie_fade
    
    fetch 1hpv, async=0
    orient
    
    #format
    bg black
    show_as cartoon
    color marine, ss s
    color red, ss h
    color white, ss l+""
    show sticks, resn 478
    util.cbao resn 478
    set surface_color, white, 1hpv
    show surface
    
    #movie
    mset 1x120
    movie_fade transparency,1, 0, 60, 1, 1hpv
    movie_fade transparency,61, 1, 120, 0, 1hpv
    

## See Also

  * [mdo](/index.php/Mdo "Mdo")
  * [mappend](/index.php/Mappend "Mappend")
  * [set](/index.php/Set "Set")
  * [movie_color_fade](/index.php/Movie_color_fade "Movie color fade")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_fade&oldid=13946](https://pymolwiki.org/index.php?title=Movie_fade&oldid=13946)"


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

## Nmr cnstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/nmr_cnstr.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/nmr_cnstr.py)  
Author(s)  | [Evangelos Papadopoulos](/index.php?title=User:Evangelos&action=edit&redlink=1 "User:Evangelos \(page does not exist\)")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  
This script will display the NMR constrains used for a structure calculation atop a structure. Usage: Save this as "NMRcnstr.py" load your protein in PyMOL, and run the script. type upl('fname') or cns('fname') where fname is the filename with the NMR constrains you want to display. It is still a very preliminary version. 

If you generated the structure by CYANA type: 

cyana> read final.pdb 

to input the structure in cyana then: 

cyana>pseudo=1 

before exporting the structure again by: 

cyana> write final.pdb 

this way the structure will contain the appropriate pseudoatoms nomeclature. 

Welcome to contact me if you need some help to set it up. 

[![Cns.png](/images/b/ba/Cns.png)](/index.php/File:Cns.png) [![Upl.png](/images/8/8f/Upl.png)](/index.php/File:Upl.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Nmr_cnstr&oldid=13915](https://pymolwiki.org/index.php?title=Nmr_cnstr&oldid=13915)"


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

## Perp maker

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/perp_maker.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/perp_maker.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The **perp_maker.py** script creates perpendicular planes. 

Nothing to do with cops. Given a simple PyMol scene, attempts to create a CGO background triangle perpendicular to the vector created - which is parallel to the line segment drawn through the camera point and current center of mass - as obtained by "get_position," or "get_view." 

## Usage

  * Load your scene
  * Orient the scene as you wish
  * [Run](/index.php/Running_Scripts "Running Scripts") the script



Could it be any simpler? 

You can rotate and move the plane using the editing features (A > drag, Shift+Drag with right, middle or left mouse button). 

Retrieved from "[https://pymolwiki.org/index.php?title=Perp_maker&oldid=13916](https://pymolwiki.org/index.php?title=Perp_maker&oldid=13916)"


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

## Plot noe

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/plot_noe.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/plot_noe.py)  
Author(s)  | [Justin Lorieau](http://www.lorieau.com/) and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
plot_noe reads a XPLOR NOE distance restraints file and shows them as distance objects in PyMOL. 

This is a modified version of [Justin's original code](http://lorieau.com/software/biophysics-software/40-plot-xplor-noes-in-pymol.html). You may visit his site for more description and examples. 

## Usage
    
    
    plot_noe filename [, selection [, line_color [, line_width ]]]
    

## See Also

  * [distance](/index.php/Distance "Distance")



Retrieved from "[https://pymolwiki.org/index.php?title=Plot_noe&oldid=13917](https://pymolwiki.org/index.php?title=Plot_noe&oldid=13917)"


---

## PoseView

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/poseview.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/poseview.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**poseview** is a PyMOL wrapper for the [PoseView](http://www.biosolveit.de/poseview/) program. PoseView generates 2D structure-diagrams of protein-ligand complexes. 

## Setup

You need the poseview executable and a valid license file. Example for a launcher script which should be saved to **/usr/bin/poseview** so that PyMOL can find it: 
    
    
    #!/bin/sh
    POSEVIEW_PATH="/opt/poseview-1.1.2-Linux-x64"
    export BIOSOLVE_LICENSE_FILE="$POSEVIEW_PATH/poseview.lic"
    exec "$POSEVIEW_PATH/poseview" "$@"
    

## Usage
    
    
    poseview [ ligand [, protein [, width [, height [, filename [, exe [, state ]]]]]]]
    

## Example
    
    
    fetch 1a00, async=0
    poseview chain C and resn HEM
    

[![Poseviewscreenshot.png](/images/c/c7/Poseviewscreenshot.png)](/index.php/File:Poseviewscreenshot.png)

Retrieved from "[https://pymolwiki.org/index.php?title=PoseView&oldid=13918](https://pymolwiki.org/index.php?title=PoseView&oldid=13918)"


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

## Propka

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/propka.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/propka.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Acknowledgement
  * 2 Introduction
  * 3 Dependency of python module: mechanize
    * 3.1 Examples
  * 4 Use of script
    * 4.1 Input paramaters
  * 5 Examples
    * 5.1 Example 1 - Mutagenesis analysis - part 1
    * 5.2 Example 2 - Mutagenesis analysis - part 2
    * 5.3 Example 3 - Scan a range of proteins
  * 6 ScriptVersion
  * 7 Changelog
    * 7.1 Known bugs



## Acknowledgement

propka.py contact and relies on the result from the [propka](http://propka.ki.ku.dk) server 

The PROPKA method is developed by the [Jensen Research Group](http://propka.ki.ku.dk/~luca/wiki/index.php/JansPage) , Department of Chemistry, University of Copenhagen.   


## Introduction

This script can fetch the pka values for a protein from the [propka](http://propka.ki.ku.dk/) server. The "propka" function downloads the results and processes them.   
It also automatically writes a pymol command file and let pymol execute it. This command file make pka atoms, rename them, label them and color them according to the pka value. 

If you put the mechanize folder and the propka.py script somewhere in your pymol search path, then getting the pka values is made super easy. By the way, did you know, that you don't have to prepare the .pdb file by adding/removing hydrogens? The [propka](http://propka.ki.ku.dk/) server uses its own internal hydrogen placement algorithm. 
    
    
    import propka
    fetch 4ins, async=0
    propka
    

If there is no web connection, it is possible to process a result file from a previous run or from a downloaded propka webpage result. This can be a handsome feature in a teaching/seminar situation, since it speeds up the pymol result or that an available web connection can be doubtful. Just point to the .pka file: Remember the dot "." which means "current directory". 
    
    
    import propka
    load 4ins.pdb
    propka pkafile=./Results_propka/4ins"LOGTIME".pka, resi=18.25-30, resn=cys
    

The last possibility, is just to ask for the pka values of a recognized PDB id. This is done with the "getpropka" function. 
    
    
    import propka
    getpropka source=ID, PDBID=4ake, logtime=_, showresult=yes
    

## Dependency of python module: mechanize

The script needs mechanize to run. The module is included in the project, [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro). 

  * On windows, it is not easy to make additional modules available for pymol. So put in into your working folder.
  * The easy manual way:


  1. Go to: <http://wwwsearch.sourceforge.net/mechanize/download.html>
  2. Download mechanize-0.2.5.zip. <http://pypi.python.org/packages/source/m/mechanize/mechanize-0.2.5.zip>
  3. Extract to .\mechanize-0.2.5 then move the in-side folder "mechanize" to your folder with propka.py. The rest of .\mechanize-0.2.5 you don't need.


  * You can also see other places where you could put the "mechanize" folder. Write this in pymol to see the paths where pymol is searching for "mechanize"
  * import sys; print(sys.path)



### Examples

Read about the proteins here:  
<http://www.proteopedia.org/wiki/index.php/4ins>   
<http://www.proteopedia.org/wiki/index.php/1hp1>   

    
    
    import propka
    
    fetch 4ins, async=0
    propka        OR
    propka 4ins   OR
    propka 4ins, resi=19.20, resn=ASP.TYR, logtime=_, verbose=yes
    
    import propka
    fetch 1hp1, async=0
    propka molecule=1hp1, chain=A, resi=305-308.513, resn=CYS, logtime=_
    
    import propka
    getpropka source=ID, PDBID=4ins, logtime=_, server_wait=3.0, verbose=yes, showresult=yes
    

  * [![propka used on 1HP1.](/images/9/94/Propka1HP1.png)](/index.php/File:Propka1HP1.png "propka used on 1HP1.")

propka used on 1HP1. 

  * [![propka used on 1HP1, zoom on ATP ligand.](/images/2/29/Propka1HP1-ATP.png)](/index.php/File:Propka1HP1-ATP.png "propka used on 1HP1, zoom on ATP ligand.")

propka used on 1HP1, zoom on ATP ligand. 

  * [![propka used on 1HP1, zoom on ZN ligand metal center.](/images/a/af/Propka1HP1-ZN.png)](/index.php/File:Propka1HP1-ZN.png "propka used on 1HP1, zoom on ZN ligand metal center.")

propka used on 1HP1, zoom on ZN ligand metal center. 



  * [![propka used on 4INS.](/images/7/75/Propka4INS.png)](/index.php/File:Propka4INS.png "propka used on 4INS.")

propka used on 4INS. 

  * [![The easy menu does it easy to click on/off ligands and bonds](/images/f/f7/Pymolpropkamenu.png)](/index.php/File:Pymolpropkamenu.png "The easy menu does it easy to click on/off ligands and bonds")

The easy menu does it easy to click on/off ligands and bonds 

  * [![An appending logfile ./Results_propka/_Results.log saves the input commands and the result for the residues specified with resi= and resn=](/images/a/a6/Resultslog.png)](/index.php/File:Resultslog.png "An appending logfile ./Results_propka/_Results.log saves the input commands and the result for the residues specified with resi= and resn=")

An appending logfile ./Results_propka/_Results.log saves the input commands and the result for the residues specified with resi= and resn= 




pka atoms are created and renamed for their pka value. That makes it easy to "click" the atom in pymol and instantly see the pka value. 

The atoms b value are also altered to the pka value, and the atoms are then spectrum colored from pka=0-14. 

The pka value of 99.9 represent a di-sulphide bond, and is colored gold and the sphere size is set a little bigger. 

If one wants to see the specified result, the logfile ./Results_propka/_Results.log saves the link to the propka server. Here one can see in an interactive Jmol appp, the interactions to the pka residues. 

## Use of script
    
    
    cd /home/tlinnet/test
    
    import propka
    
    ### The fastest method is just to write propka. Then the last pymol molecule is assumed and send to server. verbose=yes makes the script gossip mode.
    fetch 4ins, async=0
    propka
    
    ### Larger protein
    fetch 1hp1, async=0
    propka logtime=_, resi=5-10.20-30, resn=CYS.ATP.TRP, verbose=yes
    
    ### Fetch 4ins from web. async make sure, we dont execute script before molecule is loaded. The resi and resn prints the interesting results right to command line.
    fetch 4ins, async=0
    propka chain=*, resi=5-10.20-30, resn=ASP.CYS, logtime=_
    
    ### If there is no web connection, one can process a local .pka file. Either from a previous run or from a downloaded propka webpage result.
    ### Then run and point to .pka file with: pkafile=./Results_propka/pkafile.pka Remember the dot "." in the start, to make it start in the current directory.
    load 4ins.pdb
    propka pkafile=./Results_propka/4ins_.pka, resi=18.25-30, resn=cys,
    
    ### Some more examples. This molecule has 550 residues, so takes a longer time. We select to run the last molecule, by writing: molecule=1hp1
    fetch 4ins, async=0
    fetch 1hp1, async=0
    propka molecule=1hp1, chain=A, resi=300-308.513, resn=CYS.ATP.TRP, logtime=_, verbose=no, showresult=no
    propka molecule=1hp1, pkafile=./Results_propka/1hp1_.pka, verbose=yes
    

### Input paramaters
    
    
    ############################################Input parameters: propka############################################
    ############# The order of input and changable things:
    propka(molecule="NIL",chain="*",resi="0",resn="NIL",method="upload",logtime=time.strftime("%m%d",time.localtime()),server_wait=3.0,version="v3.1",verbose="no",showresult="no",pkafile="NIL")
    # method : method=upload is default. This sends .pdb file and request result from propka server.
    ## method=file will only process a manual .pka file, and write a pymol command file. No use of mechanize.
    ## If one points to an local .pka file, then method is auto-changed to method=file. This is handsome in off-line environment, ex. teaching or seminar.
    # pkafile: Write the path to .pka file. Ex: pkafile=./Results_propka/4ins_.pka
    # molecule : name of the molecule. Ending of file is assumed to be .pdb
    # chain : which chains are saved to file, before molecule file is send to server. Separate with "." Ex: chain=A.b
    # resi : Select by residue number, which residues should be printed to screen and saved to the log file: /Results_propka/_Results.log.
    ## Separate with "." or make ranges with "-". Ex: resi=35.40-50
    # resn : Select by residue name, which residues should be printed to screen and saved to the log file: /Results_propka/_Results.log.
    ## Separate with "." Ex: resn=cys.tyr
    # logtime : Each execution give a set of files with the job id=logtime. If logtime is not provided, the current time is used.
    ## Normal it usefull to set it empty. Ex: logtime=_
    # verbose : Verbose is switch, to turn on messages for the mechanize section. This is handsome to see how mechanize works, and for error searching.
    # showresult : Switch, to turn on all results in pymol command window. Ex: showresult=yes
    # server_wait=10.0 is default. This defines how long time between asking the server for a result. Set no lower than 3 seconds.
    # version=v3.1 is default. This is what version of propka which would be used.
    ## Possible: 'v3.1','v3.0','v2.0'. If a newer version is available than the current v3.1, a error message is raised to make user update the script.
    ############################################Input parameters: getpropka############################################
    ############# The order of input and changable things:
    getpropka(PDB="NIL",chain="*",resi="0",resn="NIL",source="upload",PDBID="",logtime=time.strftime("%Y%m%d%H%M%S",time.localtime()),server_wait=3.0,version="v3.1",verbose="no",showresult="no")
    # PDB: points the path to a .pdb file. This is auto-set from propka function.
    # source : source=upload is default and is set at the propka webpage.
    # source=ID, PDBID=4ake , one can print to the command line, the pka value for any official pdb ID. No files are displayed in pymol.
    # PDBID: is used as the 4 number/letter pdb code, when invoking source=ID.
    

## Examples

### Example 1 - Mutagenesis analysis - part 1

This script was developed with the intention of making analysis of possible mutants easier. For example, the reactivity of Cysteines in FRET maleimide labelling is determined by the fraction of the Cysteine residue which is negatively charged (C-). This fraction is related to its pKa value and the pH of the buffer: f(C-)=1/(10(pK-pH)+1). So, one would be interested in having the lowest possible pKa value as possible. Ideally lower than the pH of the buffer. To analyse where to make the best mutant in your protein, you could do the following for several residues. We do the mutagenesis in the command line, since we then could loop over the residues in the protein.  
So, in a loop with defined residues, this could look like the following code. Note, now we are quite happy for the result log file, since it collects the pka for the mutants. 

Download: [examples/propka_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_1.pml>" highlight="python" />

### Example 2 - Mutagenesis analysis - part 2

To only loop over surface residues, you might want to find these with the script [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues"). 

Download: [examples/propka_2.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_2.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_2.pml>" highlight="python" />

A little warning though. You need to be carefull about the rotamer you're choosing. It can happen that the first rotamer ends up being in physically non-reasonable contact distance to other residues, so atoms become overlayed. Also, the mutagenesis wizard can have the funny habit of sometimes not adding hydrogens to terminal -C or -N after mutating a residue. 

### Example 3 - Scan a range of proteins

Could be done with this script 

Download: [examples/propka_3.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_3.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_3.pml>" highlight="python" />

## ScriptVersion

Current_Version=20111202 

## Changelog

  * 20111202 
    1. This code has been put under version control. In the project, [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro).


  * 20110823 
    1. Fixed some issues with selection algebra of Ligands.
    2. Now colors Ligands automatically to purple scheme.


  * 20110822 
    1. Made the naming scheme consistent, so one can work with multiple proteins, and the grouping still works.
    2. Bonds to N-terminal and C-terminal did not show up. Fixed.
    3. If one just write "propka", the "last" molecule in the pymol object list is now assumed, instead of the first. This makes mutagenesis analysis easier.
    4. The pka difference from assumed standard values are now also displayed. Standard values are set to: pkadictio = {'ASP':3.9, 'GLU':4.3, 'ARG':12.0, 'LYS':10.5, 'HIS':6.0, 'CYS':8.3, 'TYR':10.1}
    5. The menu size is made bigger, so it can fit the long names for the bonding partners. There's also a slider that you can drag using the mouse 

    \--it's the tiny tab between the command line and the movie controls. Drag that left or right.


  * 20110821 
    1. The resn function was made wrong, only taking the last item of list. Fixed.
    2. resn= can now also be print to screen/log for ligands like ATP. Ex: resn=CYS.ATP.TRP
    3. Ligands are now also written to the "stripped-file". "./Results_propka/PDB.stripped"
    4. Bonds are now also made for Ligands.


  * 20110820 
    1. If one points to a result .pka file, then "method" is automatically set to "method=file".
    2. Added the ability for pka values for ligands
    3. Bonds are now generated for the pka atoms.
    4. The color scheme is changed from "rainbow" to "red_white_blue". This is easier to interpret.
    * CC: COULOMBIC INTERACTION. Color is red.
    * SH: SIDECHAIN HYDROGEN BOND. Color is brightorange.
    * BH: BACKBONE HYDROGEN BOND. Color is lightorange.


  * 20110817 
    1. If just invoking with "propka", it will select the first molecule. And now it is possible also to write. "propka all".
    2. Removed the "raise UserWarning" if the script is oudated. Only a warning message is printed.


  * 20110816 
    1. Made the execution of the pymol command script silent, by only using cmd.API. This will raise a Warning, which can be ignored.
    2. Built-in a Script version control, to inform the user if the propka script is updated on this page
    3. The alternate attribute for the labeling atoms are reset. It was found pymol objected altering names for atoms which had alternate positions/values.
    4. Reorganized the input order, which means that: molecule=all is default.



### Known bugs

  * Bonds can be multiplied from amino acids to Ligands like ZN or ATP. Assume the shortest bond to be correct. 
    1. ZN: This is caused, since ZN has no ID number and when there are several in the same chain. This can be reproduced for 1HP1, bond: A41ZNCC. Here it shows to bonds, where only the shortes is correct. No fix possible.
    2. ATP: This is caused, since propka sometimes eats the ending of the atom name. In the .pdb file is written O3', which propka represents like O3. Therefore a wildcard "*" is generel inserted, which can cause multiple bonds to the ATP molecule. This can be reproduced for 1HP1, bond: A504ATPSH. Here multiple bonds are made to all the O3*,O2* atoms. No fix possible.
  * The propka server sometimes use another naming scheme for Ligands. These Ligands will not be pka labelled or shown in pymol. The results will still downloaded to "/.Results_propka/file.pka" 
    1. This can be reproduced with 1BJ6. Here the propka webpage predicts the pka values for the ligands: AN7,GN1,CN3,CO2,GO2, but the .pdb file names these ligands: DA,DA,DC,DG.
  * Alternative configurations of a Ligand is at the moment a problem, and will not be shown. For example 1HXB and the ligand ROC. 

    The script extraxt AROC and BROC from pymol, which does match with ROC in the .pka file. Try save the protein with only one configuration of the Ligand.
    * >In pymol: 

    fetch 1hxb, async=0
    create 1hxbA, 1hxb and not alt B
    save 1hxbA.pdb, 1hxbA
    * >Quit pymol
    * >Manually replace with text editor in .pdb file "AROC" with " ROC" Remember the space " "!
    * >Then in pymol 

    load 1hxbA.pdb
    import propka
    propka verbose=yes



Retrieved from "[https://pymolwiki.org/index.php?title=Propka&oldid=13919](https://pymolwiki.org/index.php?title=Propka&oldid=13919)"


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

## Pymol2glmol

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/pymol2glmol.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/pymol2glmol.py)  
Author(s)  | Takanori Nakane   
License  | [LGPL3](http://opensource.org/licenses/LGPL-3.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Example of use
  * 3 System requirements
    * 3.1 Fix explained for Mint 12, gnome shell 3, chrome 16



## Introduction

A script to export a scene in pymol to GLmol. GLmol is a molecular viewer for Web browsers written in WebGL/Javascript. 

With **pymol2glmol** , you can publish your pymol scene to a Web page. Visitors can rotate, zoom the molecule on the page. 

Compared to export of polygon coordinates (VRML or Object3D), the published web page contain only atomic coordinates so that the file size is much smaller and visitors can even change representation. 

[Examples and script can be downloaded from my web page.](http://webglmol.sourceforge.jp/pymol_exporter/index.html)

This script uses **cmd.get_session** to extract which representations is enabled on each part of the molecule. I think this technique is useful for many purposes, for example, writing exporters, copying representations between aligned structures, etc. 

Comments and suggestions are welcome. 

Best regards,   
Takanori Nakane 

[![Pymol2glmol.png](/images/d/d8/Pymol2glmol.png)](/index.php/File:Pymol2glmol.png)

[](/index.php/File:Pymol2glmol.png "Enlarge")

## Example of use
    
    
    reinitialize
    import pymol2glmol
    
    fetch 1TQN, async=0
    preset.pretty_solv("1TQN")
    select heme, organic
    
    pymol2glmol 1TQN
    import webbrowser
    webbrowser.open("1TQN.html")
    

## System requirements

GLmol runs on newer versions of Firefox, Chrome, Safari or Opera.   
[See here how to activate in windows.](http://www.windows7hacker.com/index.php/2010/02/how-to-enable-webgl-in-firefox-chrome-and-safari/) Internet Explorer is not supported because IE doesn't implement WebGL. 

GLmol also runs on Sony Ericsson's Android devices which support WebGL.   
Support for Firefox Mobile is currently underway.   
Reportedly, GLmol also runs on WebGL enabled safari in iOS. 

If you see only black screen and you are using 

  * Internet Explorer: sorry. IE doesn't support WebGL.
  * Firefox (version 4 or later): [try force enable WebGL.](https://wiki.mozilla.org/Blocklisting/Blocked_Graphics_Drivers#How_to_force-enable_blocked_graphics_features)
  * Chrome: [try force enable WebGL.](http://www.google.com/support/forum/p/Chrome/thread?tid=4b9244822aa2f2e0&hl=en)
  * Safari: [enable WebGL.](https://discussions.apple.com/thread/3300585?start=0&tstart=0)



### Fix explained for Mint 12, gnome shell 3, chrome 16

Install menu editor, and run it 
    
    
    sudo apt-get install alacarte
    sudo alacarte
    

Locate internet menu. Edit the "Google Chrome" entry from->to 
    
    
    /opt/google/chrome/google-chrome %U
    /opt/google/chrome/google-chrome --enable-webgl %U
    

Then find out what handles your default opening of **text/html** , and edit that .desktop file.  
If you use the repository version of chrome, then search for **chromium-browser.desktop** instead. 
    
    
    less ~/.local/share/applications/mimeapps.list | grep 'text/html'
    sudo find / -name 'google-chrome.desktop' | xargs sudo gedit
    

Then locate the line **Exec=** and insert in line **\--enable-webgl** : 
    
    
    Exec=/opt/google/chrome/google-chrome --enable-webgl %U
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pymol2glmol&oldid=13920](https://pymolwiki.org/index.php?title=Pymol2glmol&oldid=13920)"


---

## Pymolscriptrepo

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Pymolscriptrepo&oldid=10420](https://pymolwiki.org/index.php?title=Pymolscriptrepo&oldid=10420)"


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

## Removealt

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/removealt.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/removealt.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | Free   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[removeAlt](/index.php/RemoveAlt "RemoveAlt") removes all atoms from **obj** that have alternate locations but aren't altloc **keep**. 

## Usage
    
    
    removealt [ obj [, keep ]]
    

## Example
    
    
    import removealt
    
    fetch 1hxb, async=0
    
    select ligA, alt a
    select ligB, alt b
    
    count_atoms ligA    # 49 atoms
    count_atoms ligB    # 49 atoms
    
    removealt 1hxb, b
    
    count_atoms ligA    # 0 atoms
    count_atoms ligB    # 49 atoms
    

## See Also

[Handling Alternate Locations](/index.php/Property_Selectors#Alternate_Locations "Property Selectors") and [alter](/index.php/Alter "Alter")

Retrieved from "[https://pymolwiki.org/index.php?title=Removealt&oldid=13921](https://pymolwiki.org/index.php?title=Removealt&oldid=13921)"


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

## Renumber

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/renumber.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/renumber.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
renumber sets new residue numbers (resi) for a polymer based on connectivity. 

## Example

This examples takes a pdb structure with insertion codes and sets a new, linear numbering based on [Q8N2U3_HUMAN](http://www.uniprot.org/uniprot/Q8N2U3_HUMAN). 
    
    
    fetch 1h4w, async=0
    
    # move everything which is not polymer to another chain
    alter not polymer, chain="B"
    
    # renumber polymer, first 27 residues of Q8N2U3_HUMAN missing.
    renumber chain A, 28
    

This example fixes numbering after concatenating two chains with [fuse](/index.php/Fuse "Fuse"). Note that the cartoon representation and the [sequence viewer](/index.php/Seq_view "Seq view") need [sorting](/index.php/Sort "Sort") to display correctly. 
    
    
    fab ACDEFG, chain1
    fab HIKLMN, chain2
    
    disable chain1
    as cartoon
    
    fuse last (chain1 and name C), first (chain2 and name N)
    renumber chain2
    sort chain2
    

## Troubleshooting

For large molecules it might be necessary to increase the [python recursion limit](http://docs.python.org/library/sys.html#sys.getrecursionlimit): 
    
    
    import sys
    sys.setrecursionlimit(10**5)
    

## See Also

  * [alter](/index.php/Alter "Alter")
  * [rename](/index.php/Rename "Rename")
  * [pdb_retain_ids](/index.php/Pdb_retain_ids "Pdb retain ids")



Retrieved from "[https://pymolwiki.org/index.php?title=Renumber&oldid=13922](https://pymolwiki.org/index.php?title=Renumber&oldid=13922)"


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

## RmsdByResidue

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [Zhenting Gao](/index.php/User:Zhentg "User:Zhentg")  
License  |   
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Examples
  * 6 The Code



## Introduction

RMSD between two structures of the same protein. The concept is similar as RMSF between two structures. 

  * programming 
    * Platform 
      * PyMOL 
        * rms_cur
    * Feature 
      * Output 
        * RMSD of all atoms of each residues pairs
        * Least RMSD of all atoms of each residues pairs 
          * symmetry of Phe, Tyr, His, Asp, Glu, Gln, Asn, Arg, Leu and Val needs to be considered 
            * switch the atom name and then calculate the RMSD again
            * Selected least RMSD of a residue pair for report 
              * RMSD of backbone atoms of each residues pairs
              * RMSD of C alpha atoms of each residues pairs
          * With defined residues pairs 
            * Residue pair can be limited to within binding site
    * Workflow 
      * Read reference and target pdb files 
        * two structures should be superposed before using this function
  * Note 
    * Python
    * PyMOL 
      * Clean attributes 
        * otherwise rms_cur will fail
      * How to get residue name? 
        * residue name, residue index and etc. can only be read from an atom
  * Reference 
    * <https://sourceforge.net/p/pymol/mailman/message/32710889/>



## Usage

  * Open PyMOL
  * Load PDB files
  * run this Python script inside PyMOL
  * call the function 
    * rmsdByRes pdb1, resSelection, pdb2



## Required Arguments

Text 

## Optional Arguments

Text 

## Examples

Text 

## The Code
    
    
    #load library
    from pymol import cmd, stored
    import sys
    import os
    
    #function to judge if a file exists
    def is_non_zero_file(fpath):
        return os.path.isfile(fpath) and os.path.getsize(fpath) > 0
    
    #function to switch atom names within a residue
    def switchName(residueSelection, atomName1, atomName2):
      """
      switch the name of two atoms
      """
      cmd.alter(residueSelection+" and name "+atomName1, 'name="gzt"')
      cmd.alter(residueSelection+" and name "+atomName2, 'name="'+atomName1+'"')
      cmd.alter(residueSelection+" and name gzt", 'name="'+atomName2+'"')
    
    
    #function to change atom names of some residues with symetric sidechain
    def flipAtomName(targetResidueSelection):
      """
      switch the atom names of specific residues
      """
    # Create flipped residue
      cmd.create("flippedRes",targetResidueSelection+" and not alt B")
      targetResidueCa=cmd.get_model("flippedRes and name CA")
      for g in targetResidueCa.atom:
    #   print g.resn
        if g.resn=='ARG':
         switchName("flippedRes", "NH1", "NH2")
        elif g.resn=='HIS':
         switchName("flippedRes", "ND1", "CD2")
         switchName("flippedRes", "CE1", "NE2")
        elif g.resn=='ASP':
         switchName("flippedRes", "OD1", "OD2")
        elif g.resn=='PHE':
         switchName("flippedRes", "CD1", "CD2")
         switchName("flippedRes", "CE1", "CE2")
        elif g.resn=='GLN':
         switchName("flippedRes", "OE1", "NE2")
        elif g.resn=='GLU':
         switchName("flippedRes", "OE1", "OE2")
        elif g.resn=='LEU':
         switchName("flippedRes", "CD1", "CD2")
        elif g.resn=='ASN':
         switchName("flippedRes", "OD1", "ND2")
        elif g.resn=='TYR':
         switchName("flippedRes", "CD1", "CD2")
         switchName("flippedRes", "CE1", "CE2")
        elif g.resn=='VAL':
          switchName("flippedRes", "CG1", "CG2")
      cmd.sort()
    # cmd.label("flippedRes","name")
      return "flippedRes"
    
      
    
    #main function
    def rmsdByRes(referenceProteinChain,sel, targetProteinChain):
      """
      Update
        Zhenting Gao on 7/28/2016
    
      USAGE
    
        rmsf referenceProteinChain, targetProteinChain, selection [,byres=0], [reference_state=1]
    
        Calculate the RMSD for each residue pairs from two chains of the same protein from two crystal structures.
    
      Workflow
        Read reference and target pdb files
        Align two structures
            sel target, proA and chain A
                #define target protein chain
            sel refrence, proB and chain A
                #define reference protein chain
            align target, reference
                #automatical alignment
        Clean attributes
            otherwise rms_cur will fail
    
    
      """
    # Create temporary objects, exclude alternative conformation B
      cmd.create("ref_gzt", referenceProteinChain+" and polymer and not alt B")
      cmd.alter("ref_gzt", "chain='A'")
      cmd.alter("ref_gzt", "segi=''")
      cmd.create("target_gzt", targetProteinChain+" and polymer and not alt B")
      cmd.alter("target_gzt", "chain='A'")
      cmd.alter("target_gzt", "segi=''")
    #  cmd.align("target_gzt","ref_gzt",object="align")
    # parameters
      outputText=""
      res2Check=['HIS','ASP','ARG','PHE','GLN','GLU','LEU','ASN','TYR','VAL']
           
    # select alpha carbon of selected residues in reference structure
      calpha=cmd.get_model(sel+" and name CA and not alt B")
      
      for g in calpha.atom:
    #  print g.resi+g.resn
       if cmd.count_atoms("ref_gzt and polymer and resi "+g.resi)==cmd.count_atoms("target_gzt and polymer and resi "+g.resi):
        rmsdRes=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi,"target_gzt and polymer and resi "+g.resi)
        rmsdResCa=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi+" and name ca","target_gzt and polymer and resi "+g.resi+" and name ca")
        rmsdResBackbone=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi+" and name ca+n+c+o","target_gzt and polymer and resi "+g.resi+" and name ca+n+c+o")
    #  calculate minimum rmsd
        rmsdResMin=rmsdRes
        if g.resn in res2Check:
          flippedRes=flipAtomName("target_gzt and polymer and resi "+g.resi)
          rmsdFlippedRes=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi,flippedRes)
          if rmsdFlippedRes<rmsdRes:
           rmsdResMin=rmsdFlippedRes
    
    #    print cmd.count_atoms("ref_gzt and polymer and resi "+g.resi),cmd.count_atoms("target_gzt and polymer and resi "+g.resi)
        outputText+="%s,%s,%s,%.3f,%.3f,%.3f,%.3f\n" % (targetProteinChain,g.resn,g.resi,rmsdRes,rmsdResCa,rmsdResBackbone,rmsdResMin)
      
      print outputText
    # Destroy temporary objects
      cmd.delete("ref_gzt target_gzt align res_gzt "+flippedRes)
      
    # Save data into csv
      outputFile='rmsdByRes_'+sel+'.csv'
      f=open(outputFile,'a')
      if not is_non_zero_file(outputFile):
       f.write("targe,residueName,residueId,allAtomRMSD,rmsdResCa,rmsdResBackbone,allAtomRMSDMin\n")
      f.write(outputText)
      f.close()
      print "Results saved in "+outputFile
      
    cmd.extend("rmsdByRes",rmsdByRes)
    

Retrieved from "[https://pymolwiki.org/index.php?title=RmsdByResidue&oldid=12811](https://pymolwiki.org/index.php?title=RmsdByResidue&oldid=12811)"


---

## RotationAxis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/draw_rotation_axis.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/draw_rotation_axis.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO cylinder representing a rotation axis for a given transformation. It's very useful for drawing the axes of rotational symmetry in an oligomeric assembly. 

The idea is to align two molecules/domains/chains/selections (using cmd.super) and extract the trasformation (TTT) matrix (T). The direction of the rotation axis, and a point are obtained from T and are used to create a cgo object representing the axis. The script generates two measures: one in the graphical screen (the rotation axis value, and the norm of the translation vector along the rotation axis) and some basic information in the command-line (the transformation matrix, the rotation angle, distance between centroids, and some pml lines that you can use to reproduce the axis...) 

As always with these type of things, you have to use it at your own risk. I did not try all possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado@gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run _path-to-the-script_ /draw_rotation_axis.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_axis('selection1', 'selection2')**

Please, pay attention to the apostrophes around the selections, you MUST use them. Also works with chains: 

**draw_axis('chain A', 'chain B')**

or objects: 

**draw_axis('obj01', 'obj02')**

Also, you can play a bit with the length, width and color of the axis. 

**draw_axis('selection1', 'selection2', scale_factor, width, r1, g1, b1, r2, g2, b2)**

**scale_factor** = to control the length of the axis, the default is 20 

**width** = to control the width of the axis. Default is 0.6 

**r1, g1, b1** = first RGB color code. Default is 1, 1, 1 

**r2, g2, b2** = second RGB color code to create the gradient. Default is 1, 0, 0. 

To create a single color axis, just made r1,g1,b1=r2,g2,b2 

# Examples
    
    
    # download the source and save as draw_rotation_axis.py
    run draw_rotation_axis.py
    fetch 2vak
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/d/df/2vak4.png)](/index.php/File:2vak4.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/3/37/2vak2.png)](/index.php/File:2vak2.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Another example, calculate the rotation axis of an homodimer
    run draw_rotation_axis.py
    fetch 3rfg
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/5/55/3rfg1.png)](/index.php/File:3rfg1.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/d/d6/3rfg4.png)](/index.php/File:3rfg4.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Clearly, each of the domains are made up with motifs with internal symmetry
    draw_axis('resi 200-236 and chain A', 'resi 238-274 and chain A', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # Also, you can create first the selections and use them to calculate the axis
    sele selection1, resi 200-236 and chain A
    sele selection2, resi 238-274 and chain A
    draw_axis('selection1', 'selection2', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # will produce the same result.
    

  * [![Internal rotation axis](/images/6/60/3rfg2.png)](/index.php/File:3rfg2.png "Internal rotation axis")

Internal rotation axis 

  * [![Some basic information is printed in the screen](/images/0/03/3rfg3.png)](/index.php/File:3rfg3.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 




Retrieved from "[https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923](https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923)"


---

## Rotkit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/rotkit.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/rotkit.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 Functions available in PyMOL
  * 2 Example of use
    * 2.1 Example 1 - Make a rotation of domain
    * 2.2 Example 2 - Simulate dye freedom
    * 2.3 Example 3 - Create distance distribution histogram
    * 2.4 Example 4 - A tutorial file



## Introduction

This script-kit is a collection of small script to be able to precisely to put a molecule (like a dye) where you want in relation to a protein. 

You can also create rotational states of a domain or simulate a dye freedom. 

It simply makes the [ PyMOL TTT matrixes](/index.php/Object_Matrix "Object Matrix"), in a easy and user friendly way. The calls to the functions available in PyMOL, takes care of all the conversion of input and such. 

If you are interested in this, you might also want to check out the PyMOL [Chempy](/index.php/Chempy "Chempy") module that is included in PyMOL. It provides handy vector and matrix functions. 

### Functions available in PyMOL

  * rotateline(Pos1,Pos2,degangle,molecule): 

    "Pos1->Pos2" define a line whereabout "molecule" will be rotated "degangle" degrees
    **rotateline Pos1=P513C_CA, Pos2=P513C_CB, degangle=5, molecule=Atto590**
    **rotateline Pos1=dyeatom87, Pos2=dyeatom85, degangle=10, molecule=Atto590**
  * mutate(molecule,chain,resi,target="CYS",mutframe="1"): 

    Mutate a /molecule//chain/resi into a target, and selecting most probable frame 1
    **mutate 1HP1, chain=A, resi=515, target=CYS, mutframe=1**
  * toline(Pos1,Pos2,atom,molecule,dist=1): 

    Translate molecule atom, 1 angstrom away in the same direction Pos1->Pos2 specify
    **toline Pos1=P513C_CA, Pos2=P513C_CB, atom=dyeatom87, molecule=Atto590, dist=3**



**Available through rotkit.functionname**

  * printMat(matrix): 

    prints the TTT matrix in a readable format. (4X4)
  * getxyz(Sel): 

    output is a list [x,y,z] in float. The input can be a list, a string(list) or a selection.
  * vector(Sel1,Sel2): 

    Finds the vector between points. Gets the xyz list from getxyz, so input can be anything.
  * vectorstr(vector): 

    turn a vector in list format into string. No real function actually.
  * transmat(vector,dist=1): 

    Makes a TTT translation matrix for according to the input vector. The vector is multiplied with dist.
  * unitvector(vector): 

    Make a vector a unitvector.
  * radangle(angle): 

    Convert degree to radians. Not that all input are assumed to be in degrees, and are converted automatically.
  * rotmat(angle,vectornorm,pointcoord): 

    This function is the most important. That makes the TTT matrix that rotates a molecule around a normalized vector, which goes through a coordinate point.
  * crossprod(Vector1, Vector2): 

    Makes a crossproduct between two vectors
  * crosspoint(Pos1, crossprod): 

    Returns the endpoint for the Position plus the crossproduct vector. Suitable if one would like to rotate around a crossvector.


  * [![Place your dye how you want. You always get same ending position.](/images/1/16/Rotkitex1.png)](/index.php/File:Rotkitex1.png "Place your dye how you want. You always get same ending position.")

Place your dye how you want. You always get same ending position. 

  * [![Mutate a residue through command line, and quickly put your dye there.](/images/d/de/Rotkitex2.png)](/index.php/File:Rotkitex2.png "Mutate a residue through command line, and quickly put your dye there.")

Mutate a residue through command line, and quickly put your dye there. 




## Example of use

### Example 1 - Make a rotation of domain

Download: [examples/rotkit_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_1.pml>" highlight="python" />

### Example 2 - Simulate dye freedom

Download: [examples/rotkit_2.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_2.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_2.pml>" highlight="python" />

### Example 3 - Create distance distribution histogram

Download: [examples/rotkit_3.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_3.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_3.pml>" highlight="python" />

### Example 4 - A tutorial file

To understand how the functions works, read through the tutorial. Hash/Unhash "##" each step at the time to see the effect. To be able to follow the tutorial, you need the dye molecule, which is loaded from the Pymol-script-repository. 

Download: [examples/rotkit_4.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_4.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_4.pml>" highlight="python" />

Retrieved from "[https://pymolwiki.org/index.php?title=Rotkit&oldid=13924](https://pymolwiki.org/index.php?title=Rotkit&oldid=13924)"


---

## Save Mopac

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/save_mopac.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/save_mopac.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**save_mopac** attempts to save a PDB file in the MOPAC file format. 

## Usage
    
    
    save_mopac filename [, selection [, zero [, state ]]]
    

## See Also

  * [thread on the pymol-users mailing list](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10678.html)
  * [openmopac.net](http://openmopac.net/)



Retrieved from "[https://pymolwiki.org/index.php?title=Save_Mopac&oldid=13948](https://pymolwiki.org/index.php?title=Save_Mopac&oldid=13948)"


---

## Save pdb with anisou

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/save_pdb_with_anisou.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/save_pdb_with_anisou.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Save in PDB format including ANISOU records. 

**Warning:** PyMOL does not rotate anisotropic temperature factors if you rotate an object (for example with [align](/index.php/Align "Align") or [rotate](/index.php/Rotate "Rotate")). See also related bugs on Sourceforge. 

## Usage
    
    
    save_pdb_with_anisou filename [, selection [, state ]]
    

## See Also

  * [save](/index.php/Save "Save")
  * [ellipsoids](/index.php/Ellipsoids "Ellipsoids")
  * On Sourceforge Tracker: 
    * [#3371398](https://sourceforge.net/tracker/?func=detail&aid=3371398&group_id=4546&atid=354546)
    * [#2112541](https://sourceforge.net/tracker/?func=detail&aid=2112541&group_id=4546&atid=104546)



Retrieved from "[https://pymolwiki.org/index.php?title=Save_pdb_with_anisou&oldid=13925](https://pymolwiki.org/index.php?title=Save_pdb_with_anisou&oldid=13925)"


---

## Save settings

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/save_settings.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/save_settings.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
save_settings is a prototype implementation to save all settings with non-default values to a file. 

By default settings will be stored to **~/.pymolrc-settings.py** which is automatically recognised and loaded by PyMOL on startup. You have to manually delete the file if you no more want those settings to be loaded. 

This answers a [feature request on sourceforge](https://sourceforge.net/tracker/?func=detail&aid=1009951&group_id=4546&atid=354546). 

## Usage
    
    
    save_settings [ filename ]
    

## See Also

  * [pymolrc](/index.php/Pymolrc "Pymolrc")



Retrieved from "[https://pymolwiki.org/index.php?title=Save_settings&oldid=13926](https://pymolwiki.org/index.php?title=Save_settings&oldid=13926)"


---

## Select sites

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/select_sites.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/select_sites.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3") & [Troels Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
select_sites make named selections from SITE, LINK, SSBOND, HET records.  
A super fast way to annotate the protein, according to the authors input i pdb file. 

See:   
<http://www.wwpdb.org/procedure.html#toc_10>

For instance, trypsin (<http://www.rcsb.org/pdb/files/1SGT.pdb?headerOnly=YES>): 
    
    
    SITE     1 CAT  3 HIS A  57  ASP A 102  SER A 195
    SITE     1 AC1  6 ASP A 165  ALA A 177A GLU A 180  GLU A 230
    SITE     2 AC1  6 HOH A 259  HOH A 261
    

This requires that the authors have performed annotation in the pdb file. 

## Usage
    
    
    import select_sites
    select_sites [ selection [, filename [, prefix [, nice [, quiet ]]]]]
    

nice = 0 or 1: make colored sticks representation for sites {default :1} 

## Example
    
    
    fetch 1sgt, async=0
    import select_sites
    select_sites
    

Or fetch and select_sites in one go 
    
    
    import select_sites
    sites 1sgt
    

## See Also

  * [uniprot_features](/index.php/Uniprot_features "Uniprot features")



Retrieved from "[https://pymolwiki.org/index.php?title=Select_sites&oldid=13927](https://pymolwiki.org/index.php?title=Select_sites&oldid=13927)"


---

## Show bumps

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/show_bumps.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/show_bumps.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The red discs in the image below show vdW overlaps or steric clashing. This script will create those bumps for your chosen molecule. 

[![Sb.png](/images/a/a3/Sb.png)](/index.php/File:Sb.png)

This is a script to visualize VDW clashes, just like the [mutagenesis](/index.php/Mutagenesis "Mutagenesis") wizard does. 

Posted on [pymol-users mailing list](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg09190.html). 

## Usage
    
    
    show_bumps [ selection [, name ]]
    

## Related Settings

  * sculpt_vdw_vis_max (default: 0.3)
  * sculpt_vdw_vis_mid (default: 0.1)
  * sculpt_vdw_vis_min (default: -0.1)



Retrieved from "[https://pymolwiki.org/index.php?title=Show_bumps&oldid=13928](https://pymolwiki.org/index.php?title=Show_bumps&oldid=13928)"


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

## Spectrumbar

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/spectrumbar.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrumbar.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
    * 1.1 USAGE
    * 1.2 EXAMPLES
  * 2 Script



## Description

Often times one may color their structure based upon some predefined color spectrum. However, this alone is not helpful when dealing with publication images. The purpose of this script is to allow the user to draw their own custom spectrum bar to accompany their structural images. 

Further work will be done to improve the functionality of this script. Please feel free to contact the author for function requests. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "spectrumbar" (which can be changed through the options) which is a cylinder colored according to the user-specified colors. 

Note: It is strongly recommended to turn off the specular reflections before ray tracing 

### EXAMPLES

[![](/images/b/ba/Sb_red.png)](/index.php/File:Sb_red.png) [](/index.php/File:Sb_red.png "Enlarge")Fig.1 - spectrumbar red | [![](/images/b/ba/Sb_red.png)](/index.php/File:Sb_red.png) [](/index.php/File:Sb_red.png "Enlarge")Fig.2 - spectrumbar 1,0,0 | [![](/images/f/f5/Sb_red_orange_yellow.png)](/index.php/File:Sb_red_orange_yellow.png) [](/index.php/File:Sb_red_orange_yellow.png "Enlarge")Fig.3 - spectrumbar red, orange, yellow  
---|---|---  
[![](/images/0/00/Sb_default.png)](/index.php/File:Sb_default.png) [](/index.php/File:Sb_default.png "Enlarge")Fig.4 - spectrumbar red,orange,yellow,green,blue,purple | [![](/images/0/00/Sb_default.png)](/index.php/File:Sb_default.png) [](/index.php/File:Sb_default.png "Enlarge")Fig.5 - spectrumbar 1,0,0, orange, yellow, 0,1,0, 0,0,1, purple | [![](/images/5/58/Sb_radius.png)](/index.php/File:Sb_radius.png) [](/index.php/File:Sb_radius.png "Enlarge")Fig.6 - spectrum red,radius=2  
[![](/images/6/63/Sb_long_red_orange_long_yellow.png)](/index.php/File:Sb_long_red_orange_long_yellow.png) [](/index.php/File:Sb_long_red_orange_long_yellow.png "Enlarge")Fig 7 - spectrumbar red, red, red, orange, yellow, yellow | [![](/images/b/ba/Sb_rounded.png)](/index.php/File:Sb_rounded.png) [](/index.php/File:Sb_rounded.png "Enlarge")Fig 8 - spectrumbar red, radius=2, ends=rounded  
  
## Script

load the script using the [run](/index.php/Run "Run") command 
    
    
    from pymol.cgo import *
    from math import *
    from pymol import cmd
    from re import *
    
    def spectrumbar (*args, **kwargs):
    
        """
        Author Sean M. Law
        University of Michigan
        seanlaw_(at)_umich_dot_edu
    
        USAGE
    
        While in PyMOL
    
        run spectrumbar.py
    
        spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)
    
        Parameter     Preset         Type     Description
        RGB_Colors    [1.0,1.0,1.0]  N/A      RGB colors can be specified as a
                                              triplet RGB value or as PyMOL
                                              internal color name (i.e. red)
        radius        1.0            float    Radius of cylindrical spectrum bar
        name          spectrumbar    string   CGO object name for spectrum bar
        head          (0.0,0.0,0.0)  float    Starting coordinate for spectrum bar
        tail          (10.0,0.0,0.0) float    Ending coordinate for spectrum bar
        length        10.0           float    Length of spectrum bar
        ends          square         string   For rounded ends use ends=rounded
    
        Examples:
    
        spectrumbar red, green, blue
        spectrumbar 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0
    
        The above two examples produce the same spectrumbar!
    
        spectrumbar radius=5.0
        spectrumbar length=20.0
    
        """
    
        rgb=[1.0, 1.0, 1.0]
        name="spectrumbar"
        radius=1.0
        ends="square"
        x1=0
        y1=0
        z1=0
        x2=10
        y2=0
        z2=0
        num=re.compile('[0-9]')
        abc=re.compile('[a-z]')
    
        for key in kwargs:
            if (key == "radius"):
                radius = float(kwargs["radius"])
            elif (key == "name"):
                name=kwargs["name"]
            elif (key == "head"):
                head=kwargs["head"]
                head=head.strip('" []()')
                x1,y1,z1=map(float,head.split(','))
            elif (key == "tail"):
                tail=kwargs["tail"]
                tail=tail.strip('" []()')
                x2,y2,z2=map(float,tail.split(','))
            elif (key == "length"):
                if (abc.search(kwargs["length"])):
                    print "Error: The length must be a value"
                    return
                else:
                    x2=float(kwargs["length"]);
            elif (key == "ends"):
                ends=kwargs["ends"]
            elif (key != "_self"):
                print "Ignoring unknown option \""+key+"\""
            else:
                continue
    
        args=list(args)
        if (len(args)>=1):
            rgb=[]
        while (len(args)>=1):
            if (num.search(args[0]) and abc.search(args[0])):
                if (str(cmd.get_color_tuple(args[0])) != "None"):
                    rgb.extend(cmd.get_color_tuple(args.pop(0)))
                else:
                    return
            elif (num.search(args[0])):
                rgb.extend([float(args.pop(0))])
            elif (abc.search(args[0])):
                if (str(cmd.get_color_tuple(args[0])) != "None"):
                    rgb.extend(cmd.get_color_tuple(args.pop(0)))
                else:
                    return
            else:
                print "Error: Unrecognized color format \""+args[0]+"\""
                return
        
        if (len(rgb) % 3):
            print "Error: Missing RGB value"
            print "Please double check RGB values"
            return
    
        dx=x2-x1
        dy=y2-y1
        dz=z2-z1
        if (len(rgb) == 3):
            rgb.extend([rgb[0]])
            rgb.extend([rgb[1]])
            rgb.extend([rgb[2]])
        t=1.0/(len(rgb)/3.0-1)
        c=len(rgb)/3-1
        s=0
        bar=[]
        
        while (s < c):
            if (len(rgb) >0):
                r=rgb.pop(0)
                g=rgb.pop(0)
                b=rgb.pop(0)
            if (s == 0 and ends == "rounded"):
                bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz, radius])
            bar.extend([CYLINDER])
            bar.extend([x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz])
            bar.extend([x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz])
            bar.extend([radius, float(r), float(g), float(b)])
            if (len(rgb) >= 3):
                bar.extend([float(rgb[0]), float(rgb[1]), float(rgb[2])])
                r=rgb[0]
                g=rgb[1]
                b=rgb[2]
            else:
                bar.extend([float(r), float(g), float(b)])
            if (s == c-1 and ends == "rounded"):
                bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz, radius])
            s=s+1
    
        cmd.delete(name)
        cmd.load_cgo(bar,name)
        
    
        return
    cmd.extend("spectrumbar",spectrumbar)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spectrumbar&oldid=13949](https://pymolwiki.org/index.php?title=Spectrumbar&oldid=13949)"


---

## Stereo ray

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/stereo_ray.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/stereo_ray.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This scripts provides the stereo_ray command, which saves two ray traced images in PNG format. The left image is rotated +3 degrees and the right image is rotated -3 degrees. 

_Note that the end result is almost identical to using "walleye"[stereo](/index.php/Stereo "Stereo") mode, which might be more convenient and is readily available in PyMOL._

## Contents

  * 1 Usage
  * 2 Example
  * 3 Manually
  * 4 See Also



## Usage
    
    
    stereo_ray filename [, width [, height ]]
    

## Example
    
    
    import stereo_ray 
    stereo_ray output, 1000, 600
    

This will create two images, "output_l.png" and "output_r.png". Just paste the two images next to each other (in some image editing program) and you're done. 

[![](/images/a/ad/Stereo_ex_l.png)](/index.php/File:Stereo_ex_l.png)

[](/index.php/File:Stereo_ex_l.png "Enlarge")

Left Image

[![](/images/4/49/Stereo_ex_r.png)](/index.php/File:Stereo_ex_r.png)

[](/index.php/File:Stereo_ex_r.png "Enlarge")

Right Image

[![](/images/7/7b/Stereo_complete.png)](/index.php/File:Stereo_complete.png)

[](/index.php/File:Stereo_complete.png "Enlarge")

Combined Images

  


## Manually

To obtain the same result without the script, use the [ray](/index.php/Ray "Ray") and the [png](/index.php/Png "Png") command like this: 
    
    
    ray angle=+3
    png left-image.png
    ray angle=-3
    png right-image.png
    

Obtain almost the same result using the [stereo](/index.php/Stereo "Stereo") command: 
    
    
    stereo walleye
    png combined.png, ray=1
    

## See Also

  * [stereo](/index.php/Stereo "Stereo")
  * [stereo_mode](/index.php/Stereo_Mode "Stereo Mode")



Retrieved from "[https://pymolwiki.org/index.php?title=Stereo_ray&oldid=13931](https://pymolwiki.org/index.php?title=Stereo_ray&oldid=13931)"


---

## TMalign

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/tmalign.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/tmalign.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.fitting](http://pymol.org/psicoredirect.php?psico.fitting)  
  
**tmalign** is a python module (or script) that provides wrappers to TMalign, TMscore and MMalign. The executables can be downloaded from <http://zhanglab.ccmb.med.umich.edu/TM-align/> and should be saved to any directory in [PATH](http://en.wikipedia.org/wiki/PATH_\(variable\)). 

The module also provides the command **alignwithanymethod** which is useful to quickly test different alignment methods (with their respective default values). 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Examples
  * 4 See Also



## Installation

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c speleo3 tmalign
    

## Usage
    
    
    tmalign mobile, target [, args [, exe [, ter [, transform [, object ]]]]]
    

**tmscore** and **mmalign** usage is equivalent. 
    
    
    alignwithanymethod mobile, target [, methods ]
    

## Examples
    
    
    fetch 2xwu 2x19 3gjx, async=0
    
    # TMscore example
    tmscore 2x19 and chain B, 2xwu and chain B
    
    # TMalign example with alignment object
    tmalign 3gjx and chain A, 2xwu and chain B, object=aln
    
    # full path to executable
    mmalign 3gjx, 2x19, exe=/usr/local/src/TM/MMalign
    

Example for **alignwithanymethod** (tmalign and [cealign](/index.php/Cealign "Cealign") will nicely align the beta sheet, but [align](/index.php/Align "Align") and [super](/index.php/Super "Super") will fail to find a nice superposition): 
    
    
    fetch 1C0M  1BCO, async=0
    remove 1C0M and not chain B
    alignwithanymethod 1C0M, 1BCO
    

  * [![1C0M_align01](/images/5/59/1C0M_align01.png)](/index.php/File:1C0M_align01.png "1C0M_align01")

1C0M_align01 

  * [![1C0M_super01](/images/4/44/1C0M_super01.png)](/index.php/File:1C0M_super01.png "1C0M_super01")

1C0M_super01 

  * [![1C0M_cealign01](/images/0/08/1C0M_cealign01.png)](/index.php/File:1C0M_cealign01.png "1C0M_cealign01")

1C0M_cealign01 

  * [![1C0M_tmalign01](/images/2/26/1C0M_tmalign01.png)](/index.php/File:1C0M_tmalign01.png "1C0M_tmalign01")

1C0M_tmalign01 




## See Also

  * [align](/index.php/Align "Align")
  * [super](/index.php/Super "Super")
  * [cealign](/index.php/Cealign "Cealign")



Retrieved from "[https://pymolwiki.org/index.php?title=TMalign&oldid=13932](https://pymolwiki.org/index.php?title=TMalign&oldid=13932)"


---

## ToGroup

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/togroup.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/togroup.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Examples
  * 3 The Code
  * 4 See Also



# Overview

toGroup will convert a multistate object into a group of single-state objects. _Be warned, by default it deletes your original object (since it's extracting a copy)._

PyMOL does a great job at handling multistate objects and grouping them together. One thing that I found myself doing over and over again was 

  * loading a multistate object (say a PDBQT file with 100 ligand poses)
  * splitting that object into all 100 states, with some given prefix
  * then grouping them into their own group
  * and then finally removing the original.



This became tedious, so I automated that with this script. 

# Examples
    
    
    # A multistate object (20 NMR states)
    fetch 1nmr
    
    # Create the group called,  "nmrEnsemble"
    # from '1nmr' and name all the new states state1,
    # state2, state3, etc.
    toGroup nmrEnsemble, 1nmr, prefix=state
    

# The Code
    
    
    import pymol
    from pymol import cmd
    
    def toGroup(groupName,sel,prefix="",delOrig=True):
        """
        DESCRIPTION
        toGroup will take a multistate object and extract it
        to a group with N objects all in state #1.  It essentially
        performs the following:
    
        split_states myObj, prefix=somePrefix
        group newGroup, somePrefix*
        delete myObj
    
        PARAMETERS:
        
        groupName (string)
            The name of the group to create
    
        sel (string)
            The name of the selection/object from which
            to make the group
    
        prefix (string)
            The prefix of the names of each of split states.
            For example, if your prefix is ''obj'' and is in
            states 1 through 100 then the states will be labeled
            obj1, obj2, obj3, ..., obj100.
    
        delOrig (string/boolean)
            If true then delete the original selection, otherwise not.
    
        RETURN
    
        Nothing, it makes a new group.
    
        """
        if prefix=="":
            prefix="grouped"
    
        cmd.split_states(sel, prefix=prefix)
        cmd.group(groupName,prefix+"*")
        
        if delOrig:
            cmd.delete(sel)
    
    cmd.extend("toGroup", toGroup)
    

# See Also

[group](/index.php/Group "Group"), [saveGroup](/index.php/SaveGroup "SaveGroup"), [select](/index.php/Select "Select"),, [split_states](/index.php/Split_states "Split states"), [delete](/index.php/Delete "Delete"), [extend](/index.php/Extend "Extend"). 

Retrieved from "[https://pymolwiki.org/index.php?title=ToGroup&oldid=13950](https://pymolwiki.org/index.php?title=ToGroup&oldid=13950)"


---

## Transformations

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/transformations.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/transformations.py)  
Author(s)  | Christoph Gohlke   
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Introduction

This script provides a library of functions for manipulation of matrices and vectors, useful in comparing alignments of various objects. It was written by Christoph Gohlke and is available along with many other scripts from his [website](http://www.lfd.uci.edu/~gohlke/). 

## Usage

The script must be imported in any other script that uses it. The easiest way to do this is to use the [Pymol-script-repo](/index.php/Git "Git"). 
    
    
    import transformations
    

## See also

The following scripts use transformations.py: 

  * [elbow_angle](/index.php/Elbow_angle "Elbow angle")



Retrieved from "[https://pymolwiki.org/index.php?title=Transformations&oldid=13933](https://pymolwiki.org/index.php?title=Transformations&oldid=13933)"


---

## Uniprot features

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/uniprot_features.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/uniprot_features.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
uniprot_features fetches the sequence annotation (features) from [uniprot.org](http://www.uniprot.org) and create named selections for each feature. 

Requires residue numbering (resi) to match uniprot sequence! 

## Usage
    
    
    uniprot_features uniprot_id [, selection [, withss [, prefix ]]]
    

## Example

Show the active site of [human pepsin A](http://www.uniprot.org/uniprot/PEPA_HUMAN): 
    
    
    fetch 1flh, async=0
    
    # fix residue numbering to match uniprot sequence
    alter all, resi=resv+62
    
    # create selections
    uniprot_features PEPA4_HUMAN
    
    # show sticks for active site
    as cartoon
    set cartoon_side_chain_helper
    show sticks, feature_active_site
    zoom feature_active_site
    

### See also

[select_sites](/index.php/Select_sites "Select sites")

Retrieved from "[https://pymolwiki.org/index.php?title=Uniprot_features&oldid=13934](https://pymolwiki.org/index.php?title=Uniprot_features&oldid=13934)"


---

## Wfmesh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/wfmesh.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/wfmesh.py)  
Author(s)  | [Dan Kulp](/index.php/User:Tmwsiy "User:Tmwsiy")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 DESCRIPTION
  * 2 IMAGES
  * 3 SETUP
  * 4 NOTES / STATUS
  * 5 EXAMPLES
  * 6 REFERENCES
  * 7 ADDITIONAL RESOURCES



### DESCRIPTION

This script will create an object for any Wavefront(.OBJ) mesh file. This is a way to extend the number of objects you can use. Also, you have more control over the coloring, transformations, etc than the CGOs. Although there are a number of these obj files on the web, you can also easily created them with open source tools ([OpenFX](http://www.openfx.org), [Crossroads3D](http://synapses.mcg.edu/tools/xroads/xroads.stm)). It takes literally, 2 min to get an object created and then loaded into pymol. Simply open OpenFX Designer, click File->Insert->Model, then choose any of the models (or create your own of course!), then export it as .3ds file. Then open the .3ds file from Crossroads3D and export as Wavefront OBJ. 

  * createWFMesh - create a mesh object from Wavefront (*.obj) formated file



### IMAGES

[![](/images/3/38/Starwars_pymol.png)](/index.php/File:Starwars_pymol.png)

[](/index.php/File:Starwars_pymol.png "Enlarge")

Star Wars Anyone?

[![](/images/8/8e/Torus_pymol.png)](/index.php/File:Torus_pymol.png)

[](/index.php/File:Torus_pymol.png "Enlarge")

A Torus, as an example shape you could do. Notice polygon normals are being used...need smoothing!

### SETUP

Simply "import wfmesh.py" 

### NOTES / STATUS

  * Tested on Pymolv0.97, Windows platform, should work on linux as well.
  * Coloring is fixed for grey and sections of mesh are stored, but not used.
  * Simple opengl calls; not optimized (display lists, etc) or anything.
  * Vertex Normal code is broken, so normals are per polygon right now.
  * Post problems in the discussion page, on 'my talk' page or just email me : dwkulp@mail.med.upenn.edu



### EXAMPLES
    
    
    import wfmesh
    cd /home/tlinnet/Software/pymol/Pymol-script-repo/files_for_examples
    wfmesh.createWFObj('torus.obj','Torus',flip=1)    # Flip = 1, if OBJ created by openFX, crossroads3D combination
    wfmesh.createWFObj("torus.obj","Torus2",translate=[5,5,0], flip=1)
    wfmesh.createWFObj("ship.obj","Ship")
    

### REFERENCES

[OpenFX](http://www.openfx.org) [Crossroads3D](http://synapses.mcg.edu/tools/xroads/xroads.stm)

### ADDITIONAL RESOURCES

Torus.obj [Torus.zip](http://pymolwiki.org/images/9/98/Torus.zip)

Retrieved from "[https://pymolwiki.org/index.php?title=Wfmesh&oldid=13935](https://pymolwiki.org/index.php?title=Wfmesh&oldid=13935)"


---

