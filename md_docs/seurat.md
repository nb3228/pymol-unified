# seurat

seurat

### Table of Contents

  * SEURAT Integration

    * Development Version (future releases)

    * SEURAT 5.0

    * Older SEURAT Releases




# SEURAT Integration

PyMOL directly integrates with the SEURAT collaborative data-mining tool from [Synaptic Science LLC](http://synapticscience.com "http://synapticscience.com"). For information on evaluating SEURAT, please complete [this form](http://synapticscience.com/contact.html "http://synapticscience.com/contact.html"). 

The documentation below is organized by SEURAT version since each subsequent release of SEURAT will likely include changes and improvements as we work with Synaptic Science and our mutual customers to implement an effective interface between the two packages. 

## Development Version (future releases)

Incorporating customer feedback, this version contains a new “Complex” Result Type for consolidating co-crystal structure complexes into a single column of a Seurat spreadsheet. Use of this Result Type requires consolidation of all atom selection information into a single string so that PyMOL can appropriately break the structure into its constituent parts. 

The format of the Expt Result Disc field for the new **Complex** Result Type is: 
    
    
    filename | display-name | display-sele | align-sele | ligand-sele | target-sele | other-sele

An updated excerpt from the pri_mapping.properties file is shown below: 
    
    
    # syn_pdb_rules allows customers to control which
    # "assay names" are treated as special columns for controlling
    # the shipping of information to and from 3D molecular modelling
    # tools like PyMol, JMol and MarvinSpace
    #
    syn_pdb_rules={"XTAL-ABL (Complex)", "short", "|"},\
    	{"XTAL-CDK2 (Complex)", "short", "|"},\
    	{"XTAL-CHK1 (Complex)", "short", "|"},\
    	{"XTAL-LCK (Complex)", "short", "|"},\
    	{"XTAL-PI3K (Complex)", "short", "|"},\
    	{"XTAL-PKA (Complex)", "short", "|"},\
    	{"XTAL-SRC (Complex)", "short", "|"},\
    	{"XTAL-SYK (Complex)", "short", "|"},
    
    syn_pdb_viewer=com.delsci.integration.Pymol3DViewerImpl

Additional regarding this result type will be made available when the next version ships. 

## SEURAT 5.0

This version introduces automatic superposition capabilities for PDB files exposed through SEURAT. Selected portions of the structure, such as ligands only, can be visualized and compared independently. 

  1. PDB files containing the structures to expose should be placed in the “sandbox” folder of the Seurat Server installation. It it not necessary to name these files using corporate identifers.

  2. The **pri_mapping.properties** file should be modified to contain a syn_pdb_rules entry for the Assay Type(s) you choose to associate with registered structures.

  3. Assay results should be then loaded to associate the various structural visualizations with specific corporate identifiers. 




For purposes of registration within SEURAT, crystal structures are simply another kind of assay result available for query and retrieval inside the Seurat Client interface. The key fields to specify when registering such data are as described here: 

  * The Corporate ID should match the structure of the compound as registered in your database

  * The Assay Name associated with the result should typically be something like XTAL-CMET (for cMet crystal structuces)

  * The Result Type should be one of Target, Ligand, or Other. 

  * A multi-part text string (detailed below) should be loaded into the Expt Result Desc field.




A representative excerpt from the pri_mapping.properties file is shown below. Note that you need to add an entry for each combination of Assay Name and Result Type you wish to expose. 
    
    
    # syn_pdb_rules allows customers to control which
    # "assay names" are treated as special columns for controlling
    # the shipping of information to and from 3D molecular modelling
    # tools like PyMol, JMol and MarvinSpace
    #
    syn_pdb_rules={"XTAL-CDK2 (Target)", "short", "|"},\
    	{"XTAL-CDK2 (Ligand)", "short", "|"},\
    	{"XTAL-CDK2 (Other)", "short", "|"},\
    	{"XTAL-CHK1 (Target)", "short", "|"},\
    	{"XTAL-CHK1 (Ligand)", "short", "|"},\
    	{"XTAL-CHK1 (Other)", "short", "|"},\
    	{"XTAL-LCK (Target)", "short", "|"},\
    	{"XTAL-LCK (Ligand)", "short", "|"},\
    	{"XTAL-LCK (Other)", "short", "|"}
    
    syn_pdb_viewer=com.delsci.integration.Pymol3DViewerImpl

The format of the Expt Result Desc string for the **Target** , **Ligand** , and **Other** Result Types is as follows. Note that a vertical bar character `|` is used to separate the fields of the string. 
    
    
    filename|display-name|display-selection|align-selection

For example, the following two Expt Result Desc strings were registered for a Result Type of **Target** for the 2C5N PDB file, which happens to contains two copies of the protein-ligand complex. 
    
    
    2c5n.pdb|2C5Na|polymer&A//|polymer&A//
    2c5n.pdb|2C5Nb|polymer&C//|polymer&C//

Translated into plain English, these strings convey that the structures to expose derive from the 2c5n.pdb file and should be presented in the Seurat spreadsheet with the labels “2C5Na” and “2C5Nb”. For 2C5Na, the atoms to be displayed map to protein chain A, whereas for 2C5Nb, they map to protein chain C. In this case, the same atom selections used for display are repeated and used for the alignments. 

Next, the following strings were registered for a Result Type of **Ligand** : 
    
    
    2c5n.pdb|2C5Na|C/CK8`1297/|polymer&C//
    2c5n.pdb|2C5Nb|C/CK8`1298/|polymer&C//

As above, the structures derive from 2c5n.pdb and should be displayed in the spreadsheet as 2C5Na and 2C5Nb. This time however, the atom selections to display correspond to residue CK8 with identifier 1297 in chain C (for 2C5Na), and identifier 1298 in chain C (for 2C5Nb). However, the structure alignments are performed using the same protein atom selections used in the prior example above. 

Finally, for a result type of **Other** , the following strings were registered: 
    
    
    2c5n.pdb|2C5Na|Z//&!(polymer or resn CK8)|polymer&A//
    2c5n.pdb|2C5Nb|X//&!(polymer or resn CK8)|polymer&C//

where the strings indicate that only the atoms in chain X (for 2C5Na) or Z (for 2C5Nb) not matching either protein atoms or residues CK8 should be displayed. The “Other” selection would typically map to solvent atoms, inorganic salts, and anything else in the complex. As before, protein chain A (or C) are used the perform the alignments. 

An excerpt from an Excel CSV file used to register this data via the Seurat client interface can be [downloaded here](http://pymol.org/synsci/seurat50example.csv "http://pymol.org/synsci/seurat50example.csv") for your perusal. Note that reasonable default values are required for some of the other fields (Lot Number, Assay Protocol, and Expt Result Units) when uploading data via this method. 

If you have a lot of PDB files to register, then we recommend that you seek assistance from Schrödinger in developing one or more PyMOL Python scripts which may be able to automatically generate the result description strings from your PDB files. 

## Older SEURAT Releases

  1. PDB files and/or PyMOL session (.PSE) files are placed in the “sandbox” folder inside the Seurat Server with a filename prefix matching the corporate identifier.

  2. An asterisk will appear next to the corporate indentifier the SEURAT spreadsheet indicating that a structure is available.

  3. Users can right-click on the corporate identifier to open these files for visualization inside of PyMOL.




seurat.txt · Last modified: 2013/08/19 21:00 (external edit)
