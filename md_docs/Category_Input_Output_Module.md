# Category: Input Output Module

## Delete

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**delete** removes objects or selections matching an expression _name_ , which can include wildcards. 

### USAGE
    
    
    delete name  
    delete all   # deletes all objects
    

name = name of object or selection 

### PYMOL API
    
    
    cmd.delete(string name = object-or-selection-name )
    

Note that special care needs to be taken to ensure that the object or selection name does not contain any quotes when passed as an argument. 

### SEE ALSO

[Remove](/index.php/Remove "Remove")

Retrieved from "[https://pymolwiki.org/index.php?title=Delete&oldid=9227](https://pymolwiki.org/index.php?title=Delete&oldid=9227)"


---

## Load

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**load** reads several file formats. The file extension is used to determine the format. PDB files must end in ".pdb", MOL files must end in ".mol", Macromodel files must end in ".mmod", XPLOR maps must end in ".xplor", CCP4 maps must end in ".ccp4", Raster3D input (Molscript output) must end in ".r3d", PyMOL session files must end in ".pse", and pickled ChemPy models with a ".pkl" can also be directly read. 

If an object is specified, then the file is loaded into that object. Otherwise, an object is created with the same name as the file prefix. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 ARGUMENTS
  * 4 NOTES
  * 5 User Comments/Examples
  * 6 SEE ALSO



### USAGE
    
    
    load filename [,object [,state [,format [,finish [,discrete [,multiplex ]]]]]]
    

### PYMOL API
    
    
    cmd.load( filename [,object [,state [,format [,finish [,discrete [,multiplex ]]]]]] )
    

### ARGUMENTS

  * **filename : string** Path or URL to the file to load.
  * **object : string** Name of Pymol object to store the structure in. Defaults to the filename prefix.
  * **state : integer** State number to store the object in, or 0 to append (default:0)
  * **format : string** Format for the file (see notes). Defaults to the filename extension.
  * **finish : integer**
  * **discrete : integer** For multi-model structures, a value of 0 indicates that models have the same set of atoms (e.g. trajectory files or NMR structures), allowing memory savings, while a value of 1 forces the creation of independent atom sets for each model {default: -1 (file type dependent)} (see [discrete objects](/index.php/Discrete_objects "Discrete objects"))
  * **quiet : integer** (default 1)
  * **multiplex : integer** Load a multi-model file as separate objects instead of states (see also [split_states](/index.php/Split_states "Split states"))
  * **zoom : integer** {default: -1 = use [auto_zoom](/index.php/Auto_zoom "Auto zoom") setting}
  * **partial : integer** For session files (.pse). partial=0: discard current session. partial=1: merge with current session (will not load global settings, selections, movie, camera). partial=2: like 1, but also load camera view {default: 0}
  * **mimic : integer** For .mae files, match style from file as close as possible, uses atom-level settings (like [cartoon_color](/index.php/Cartoon_color "Cartoon color")) {default: 1}
  * **object_props : string** = Space delimited list of property names (or * wildcard) to load from .sdf or .mae files {default: use [load_object_props_default](/index.php?title=Load_object_props_default&action=edit&redlink=1 "Load object props default \(page does not exist\)") setting} _Incentive PyMOL 1.6+_
  * **atom_props : string** = Space delimited list of property names (or * wildcard) to load from .mae files {default: use [load_atom_props_default](/index.php?title=Load_atom_props_default&action=edit&redlink=1 "Load atom props default \(page does not exist\)") setting} _Incentive PyMOL 1.6+_



### NOTES

  * You can override the file extension by giving a format string:


    
    
    'pdb' : PDB,  'mmod' : Macromodel, 'xyz' : Tinker, 'cc1' : ChemDraw3D  
    'mol' : MDL MOL-file, 'sdf' : MDL SD-file
    'xplor' : X-PLOR/CNS map, 'ccp4' : CCP4 map,
    'callback' : PyMOL Callback object (PyOpenGL)
    'cgo' : compressed graphics object (list of floats)
    'trj' : AMBER trajectory (use load_traj command for more control)
    'top' : AMBER topology file 'rst' : AMBER restart file
    'cex' : Metaphorics CEX format
    'pse' : PyMOL Session file
    'pqr' : PQR (a modified PDB file with charges and radii)
    'mol2' : MOL2
    

  * A new feature has been added to load. You can specify an URL to a PDB and PyMOL will download it. This is a very handy feature for loading experimental/theoretical data from servers across the web. Just copy the link they give you and do,


    
    
    load http://www.theurl.com/path/toYourData
    

or you can open a remote file just from the command line 
    
    
    # load a PDB file; I placed one example file on the PyMOL Wiki
    pymol http://www.pymolwiki.org/1rsy.pdb
    
    # PyMOL can also handle the gzipped files on the PDB.  :-)
    pymol ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/00/pdb100d.ent.gz
    

### User Comments/Examples

  * Load xyz.pdb using the PyMOL API:


    
    
    cmd.load("xyz.pdb")
    

  * Loading multiple PDBs into one object with many states:


    
    
    load trj0001.pdb, mov
    load trj0002.pdb, mov
    load trj0003.pdb, mov
    etc.
    

or, if you have too many states to do that by hand, 
    
    
    for idx in range(1,1001):cmd.load("trj%04d.pdb"%idx,"mov")
    

or, you can use "glob" from Python, 
    
    
    from glob import glob
    lst = glob("trj*.pdb")
    lst.sort()
    for fil in lst: cmd.load(fil,"mov")
    

  * Load a NAMD multi-PDB file. These are just concatenated PDB files.


    
    
    load NAMDtrajFile.pdb, multiplex=0
    

  * Hint: You can save some time & a lot of memory by loading the file and removing the atoms in a single-line compound statement (with a semicolon



after the load statement). 
    
    
    load 1E3M.pdb; remove not A-C+F//
    

  * Decorating the load command to include technical info about the loaded object


    
    
    def load_with_props(fileName, objName):
      # store whatever info you want, like the filename
      obj_info[objName] = fileName
      # ... do more recording of properties you choose
      # ask PyMOL to now load the file
      cmd.load(fileName,objName)
    

Then on can query obj_info based on the object name: 
    
    
    for obj in cmd.get_names():
      print obj_info[obj]
    

### SEE ALSO

[Fetch](/index.php/Fetch "Fetch") [Save](/index.php/Save "Save") [Load Traj](/index.php/Load_Traj "Load Traj")

Retrieved from "[https://pymolwiki.org/index.php?title=Load&oldid=12717](https://pymolwiki.org/index.php?title=Load&oldid=12717)"


---

## Load CGO

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**cmd.load_cgo()** loads Compiled Graphics Objects into the PyMOL session. CGOs are defined by lists of floating point numbers that define either a series of triangles/vertices and colors, or one of several predefined sequences that describe mathematical shapes. 

### PYMOL API
    
    
    cmd.load_cgo(object, name)
    

### ARGUMENTS

  * **object : list of floats** Coordinate data for the CGO.
  * **name : string** Name of the object (for the internal GUI)



### SEE ALSO

[CGO Shapes](/index.php/CGO_Shapes "CGO Shapes")

Retrieved from "[https://pymolwiki.org/index.php?title=Load_CGO&oldid=13474](https://pymolwiki.org/index.php?title=Load_CGO&oldid=13474)"


---

## Multifilesave

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The multifilesave command saves each object and/or state in a selection to a separate file. 

_New in PyMOL 2.1_

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 See Also



## Usage
    
    
    multifilesave filename [, selection [, state [, format ]]]
    

The **filename** argument must contain one or more placeholders: 
    
    
    {name}  : object name
    {state} : state number
    {title} : state title
    {num}   : file number
    {}      : object name (first) or state (second)
    

## Arguments

  * **filename** = str: file path with placeholders
  * **selection** = str: atom selection {default: all}
  * **state** = int: object state (-1 for current, 0 for all states) {default: -1}
  * **format** = str: file format, e.g. `pdb` {default: guess from file extension}



## Examples
    
    
    multifilesave /tmp/{name}.pdb
    multifilesave /tmp/{name}-{state}.cif, state=0
    multifilesave /tmp/{}-{}.cif, state=0
    multifilesave /tmp/{}-{title}.sdf, state=0
    

## See Also

  * [save](/index.php/Save "Save")
  * [loadall](/index.php/Loadall "Loadall")



Retrieved from "[https://pymolwiki.org/index.php?title=Multifilesave&oldid=12791](https://pymolwiki.org/index.php?title=Multifilesave&oldid=12791)"


---

## Quit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

quit terminates the program. 

_Changed in PyMOL 1.6: Added return code_

## Usage
    
    
    quit [code]
    

## Arguments

  * **code** = int: Return code {default: 0}



## Python API
    
    
    cmd.quit(code: int)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Quit&oldid=13342](https://pymolwiki.org/index.php?title=Quit&oldid=13342)"


---

## Save

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**save** writes selected atoms to a file. 

## Contents

  * 1 Details
    * 1.1 USAGE
    * 1.2 EXAMPLES
    * 1.3 PYMOL API
    * 1.4 NOTES
    * 1.5 SEE ALSO



# Details

  * The file format is autodetected if the extension is .pdb, .pqr, .mol, .sdf, .pkl, .pkla, .mmd, .out, .dat, .mmod, .pmo, .pov, .png, .pse, .psw, .aln, .fasta, .obj, .mtl, .wrl, .idtf, .dae, or .mol2.
  * If the file extension is ".pse" (PyMOL Session), the complete PyMOL state is always saved to the file (the selection and state parameters are thus ignored).
  * CLUSTALW formatted alignments can be written by PyMOL as well. Once you perform an alignment like the following,



    

    
    
    
    align proteinA, proteinB, object=A_on_B
    

    the alignment can be written using:
    
    
    
    save A_aligned_with_B.aln, A_on_B
    

### USAGE
    
    
    save file [,(selection) [,state [,format]] ]
    

### EXAMPLES
    
    
    # save only the alpha carbons
    save onlyCAs.pdb, n. CA
    
    # save my MD trajectory file to disk
    save myTraj.pdb, myMDTrajectory, state=0
    
    # save a PyMOL session
    save thisSession.pse
    

### PYMOL API
    
    
    cmd.save(filename[, selection[, state[, format]]])
    

### NOTES

  * When saving a session file, then "state" has no effect.
  * Default is state = -1, which saves only the current state.
  * When state = 0, all states in the file are written. If you have more than one state, this produces a multi-state PDB file.



### SEE ALSO

[Load](/index.php/Load "Load"), [Get Model](/index.php/Get_Model "Get Model")

Retrieved from "[https://pymolwiki.org/index.php?title=Save&oldid=12445](https://pymolwiki.org/index.php?title=Save&oldid=12445)"


---

