# Category: Commands

## Abort

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Abort is a dummy function which returns [Python's None element](/index.php?title=Python:None&action=edit&redlink=1 "Python:None \(page does not exist\)"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Abort&oldid=7430](https://pymolwiki.org/index.php?title=Abort&oldid=7430)"


---

## Alias

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**alias** allows you to bind a commonly used command to a single PyMOL keyword. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    alias name, command-sequence
    

### PYMOL API
    
    
    cmd.alias(string name,string command)
    

### EXAMPLES
    
    
    alias go,load $TUT/1hpv.pdb; zoom 200/; show sticks, 200/ around 8 go
    

### NOTES

For security reasons, new PyMOL commands created using "extend" are not saved or restored in sessions. 

### SEE ALSO

[extend](/index.php/Cmd_extend "Cmd extend"), [api](/index.php?title=Cmd_api&action=edit&redlink=1 "Cmd api \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Alias&oldid=7431](https://pymolwiki.org/index.php?title=Alias&oldid=7431)"


---

## Align

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/6/6e/After_alignment.png)](/index.php/File:After_alignment.png)

[](/index.php/File:After_alignment.png "Enlarge")

Two proteins after structure alignment

**align** performs a sequence alignment followed by a structural superposition, and then carries out zero or more cycles of refinement in order to reject structural outliers found during the fit. align does a good job on proteins with decent sequence similarity (identity >30%). For comparing proteins with lower sequence identity, the [super](/index.php/Super "Super") and [cealign](/index.php/Cealign "Cealign") commands perform better. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Alignment Objects
  * 4 RMSD
  * 5 Examples
  * 6 PyMOL API
  * 7 Notes
  * 8 See Also



## Usage
    
    
    align mobile, target [, cutoff [, cycles
        [, gap [, extend [, max_gap [, object
        [, matrix [, mobile_state [, target_state
        [, quiet [, max_skip [, transform [, reset ]]]]]]]]]]]]]
    

## Arguments

  * **mobile** = string: atom selection of mobile object
  * **target** = string: atom selection of target object
  * **cutoff** = float: outlier rejection cutoff in RMS {default: 2.0}
  * **cycles** = int: maximum number of outlier rejection cycles {default: 5}
  * **gap, extend, max_gap** : sequence alignment parameters
  * **object** = string: name of alignment object to create {default: (no alignment object)}
  * **matrix** = string: file name of substitution matrix for sequence alignment {default: BLOSUM62}
  * **mobile_state** = int: object state of mobile selection {default: 0 = all states}
  * **target_state** = int: object state of target selection {default: 0 = all states}
  * **quiet** = 0/1: suppress output {default: 0 in command mode, 1 in API}
  * **max_skip** = ?
  * **transform** = 0/1: do superposition {default: 1}
  * **reset** = ?



## Alignment Objects

An alignment object can be created with the **object=**_somename_ argument. An alignment object provides: 

  * aligned [sequence viewer](/index.php/Seq_view "Seq view")
  * graphical representation of aligned atom pairs as lines in the 3D viewer
  * can be [saved](/index.php/Save "Save") to a clustalw sequence alignment file



## RMSD

The RMSD of the aligned atoms (after outlier rejection!) is reported in the text output. The **all-atom RMSD** can be obtained by setting **cycles=0** and thus not doing any outlier rejection. The RMSD can also be captured with a python script, see the API paragraph below. Note that the output prints "RMS" but it is in fact "RMSD" and the units are Angstroms. 

## Examples
    
    
    fetch 1oky 1t46, async=0
    
    # 1) default with outlier rejection
    align 1oky, 1t46
    
    # 2) with alignment object, save to clustalw file
    align 1oky, 1t46, object=alnobj
    save alignment.aln, alnobj
    
    # 3) all-atom RMSD (no outlier rejection) and without superposition
    align 1oky, 1t46, cycles=0, transform=0
    

## PyMOL API
    
    
    cmd.align( string mobile, string target, float cutoff=2.0,
               int cycles=5, float gap=-10.0, float extend=-0.5,
               int max_gap=50, string object=None, string matrix='BLOSUM62',
               int mobile_state=0, int target_state=0, int quiet=1,
               int max_skip=0, int transform=1, int reset=0 )
    

This returns a list with 7 items: 

  1. RMSD after refinement
  2. Number of aligned atoms after refinement
  3. Number of refinement cycles
  4. RMSD before refinement
  5. Number of aligned atoms before refinement
  6. Raw alignment score
  7. Number of residues aligned



## Notes

  * The molecules you want to align need to be in **two different objects**. Else, PyMOL will answer with: _ExecutiveAlign: invalid selections for alignment._ You can skirt this problem by making a temporary object and aligning your original to the copy.
  * By defaults, **all states** (like in NMR structures or trajectories) are considered, this might yield a bad or suboptimal alignment for a single state. Use the **mobile_state** and **target_state** argument to be explicit in such cases.



## See Also

  * [super](/index.php/Super "Super"), [cealign](/index.php/Cealign "Cealign"), [fit](/index.php/Fit "Fit"), [pair_fit](/index.php/Pair_fit "Pair fit")
  * [rms](/index.php/Rms "Rms"), [rms_cur](/index.php/Rms_cur "Rms cur"),
  * [intra_fit](/index.php/Intra_fit "Intra fit"), [intra_rms](/index.php/Intra_rms "Intra rms"), [intra_rms_cur](/index.php/Intra_rms_cur "Intra rms cur")
  * [extra_fit](/index.php/Extra_fit "Extra fit")
  * [align_all.py and super_all.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/)
  * [tmalign](/index.php/Tmalign "Tmalign")
  * [Color_by_conservation](/index.php/Color_by_conservation "Color by conservation")
  * [Get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")
  * [mcsalign](/index.php/Mcsalign "Mcsalign") (psico)



Retrieved from "[https://pymolwiki.org/index.php?title=Align&oldid=12703](https://pymolwiki.org/index.php?title=Align&oldid=12703)"


---

## Cealign

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/d/d8/Cealign_ex1.png)](/index.php/File:Cealign_ex1.png)

[](/index.php/File:Cealign_ex1.png "Enlarge")

cealign superposition of 1c0mB and 1bco

cealign aligns two proteins using the CE algorithm. It is very robust for proteins with little to no sequence similarity (twilight zone). For proteins with decent structural similarity, the [super](/index.php/Super "Super") command is preferred and with decent sequence similarity, the [align](/index.php/Align "Align") command is preferred, because these commands are much faster than cealign. 

_This command is new in PyMOL 1.3, see the[cealign plugin](/index.php/Cealign_plugin "Cealign plugin") for manual installation._

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    cealign target, mobile [, target_state [, mobile_state
        [, quiet [, guide [, d0 [, d1 [, window [, gap_max
        [, transform [, object ]]]]]]]]]]
    

## Arguments

**Note** : The **mobile** and **target** arguments are swapped with respect to the [align](/index.php/Align "Align") and [super](/index.php/Super "Super") commands. 

  * **target** = string: atom selection of target object (CE uses only CA atoms)
  * **mobile** = string: atom selection of mobile object (CE uses only CA atoms)
  * **target_state** = int: object state of target selection {default: 1}
  * **mobile_state** = int: object state of mobile selection {default: 1}
  * **quiet** = 0/1: suppress output {default: 0 in command mode, 1 in API}
  * **guide** = 0/1: only use "guide" atoms (CA, C4') {default: 1}
  * **d0, d1, window, gap_max** : CE algorithm parameters
  * **transform** = 0/1: do superposition {default: 1}
  * **object** = string: name of alignment object to create {default: (no alignment object)}



## Example
    
    
    fetch 1c0mB 1bco, async=0
    as ribbon
    cealign 1bco, 1c0mB, object=aln
    

## See Also

  * [super](/index.php/Super "Super")
  * [align](/index.php/Align "Align")
  * [cealign plugin](/index.php/Cealign_plugin "Cealign plugin")



Retrieved from "[https://pymolwiki.org/index.php?title=Cealign&oldid=12745](https://pymolwiki.org/index.php?title=Cealign&oldid=12745)"


---

## Alter

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **iterate** command executes a Python expression for all atoms in a selection. The local namespace exposes all atomic identifiers and properties as read-only Python variables. The global namespace by default is the **pymol** module namespace, and thus exposes objects like **cmd** and **[stored](/index.php/Pymol.stored "Pymol.stored")**. The order of iteration is that of the internal atom ordering. 

The **alter** command is equivalent to **iterate** , but provides read-write access to the exposed variables. This can for example be used to rename a chain or to assign new b-factors. If changing identifiers which would affect atom ordering, calling [sort](/index.php/Sort "Sort") is necessary to reorder based on the new identifiers. 

The **iterate_state** command is similar to **iterate** , but iterates over the coordinates in the given state and selection, and exposes **x, y, z** in addition to the atomic identifiers and properties. 

The **alter_state** command is equivalent to **iterate_state** , but allows modification of **x, y, z** (and atom-state level settings in Incentive PyMOL 1.7+). 

## Contents

  * 1 Exposed Variables
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
    * 3.1 Example: Print b-factors
    * 3.2 Example: Get b-factor list
    * 3.3 Example: Get coordinates
    * 3.4 Example: Modify coordinates
    * 3.5 Example: Rename chain
    * 3.6 Example: Sequence numbering offset
    * 3.7 Example: Net Charge
    * 3.8 Example: Color transfer
  * 4 PYMOL API
    * 4.1 Arguments
    * 4.2 "space" argument
  * 5 See Also



## Exposed Variables

All commands in the **iterate** -family expose the following variables: 

  * **model** (str): the object name (appearing in the selection panel on the right) (cannot be altered)
  * **name** (str): the atom name
  * **resn** (str): the residue name
  * **oneletter** (str): **resn** translated to one letter code, like it's displayed in the sequence viewer (cannot be altered, alter **resn** instead) (PyMOL 1.8.7+)
  * **resi** (str): the residue identifier (residue number) as a string, including optional insertion code
  * **resv** (int): the residue identifier (residue number) as an integer, excluding insertion code
  * **chain** (str): the chain name
  * **alt** (str): alternate location identifier
  * **elem** (str): the chemical element
  * **q** (float): occupancy
  * **b** (float): the B Factor
  * **segi** (str): segment identifier (columns 73-76 in PDB file)
  * **type** (str: ATOM,HETATM): the atom type (PDB convention for canonical polymer residues)
  * **formal_charge** (int): the formal charge of the atom
  * **partial_charge** (float): the partial charge of the atom
  * **numeric_type**
  * **text_type** (str): automatic mol2 atom typing (Incentive PyMOL 1.4+)
  * **stereo** (str): automatic stereochemical R/S label (Incentive PyMOL 1.4+)
  * **ID** (int): PDB atom id (not guaranteed to be unique)
  * **rank** (int): atom index from original file import (not guaranteed to be unique)
  * **index** (int): internal atom index (unique per object, sensitive to [sorting](/index.php/Sort "Sort") and [removing](/index.php/Remove "Remove") of atoms, cannot be altered)
  * **vdw** (float): Van der Waals radius
  * **ss** (str): secondary structure
  * **color** (int): color index
  * **reps** (int): numeric mask of shown representations (PyMOL 1.7.4+)
  * **protons** (int): (PyMOL 1.7.4+)
  * **p** (object): property object to access user-defined properties like **p.abc** (Incentive PyMOL 1.7+)
  * **s** (object): settings object to access settings, e.g. **s.sphere_scale** (Incentive PyMOL 1.7+, Open-Source PyMOL 1.8.1+)
  * **label** (str): see [label](/index.php/Label "Label")
  * **geom** (int): 1=single, 2=linear, 3=planar, 4=tetrahedral (see [set_geometry](/index.php?title=Set_geometry&action=edit&redlink=1 "Set geometry \(page does not exist\)"))
  * **valence** (int): expected number of bonds (?)
  * **flags** (int): bitmask, see [flag](/index.php/Flag "Flag")
  * **cartoon** (int): see [cartoon](/index.php/Cartoon "Cartoon")



The **iterate_state** and **alter_state** commands in addition expose: 

  * **x** (float): x-coordinate
  * **y** (float): y-coordinate
  * **z** (float): z-coordinate
  * **state** (int): (PyMOL 1.7.0+)



## Usage
    
    
    iterate (selection), expression
    
    
    
    iterate_state state, (selection), expression
    
    
    
    alter (selection), expression
    
    
    
    alter_state state, (selection), expression
    

### Arguments

  * **state** = int: object state, -1 for current state, 0 for all states
  * **selection** = str: atom selection
  * **expression** = str: expression in Python language



## Examples

### Example: Print b-factors

The following prints the atom names (**name** variable) and b-factors (**b** variable) for residue #1. 
    
    
    iterate (resi 1), print(name + " %.2f" % b)
    

### Example: Get b-factor list

The following example fills a list, **[stored](/index.php/Pymol.stored "Pymol.stored").bfactors** with the f-factors (**b** variable) from residue #1. 
    
    
    stored.bfactors = []
    iterate (resi 1), stored.bfactors.append(b)
    print(stored.bfactors)
    

### Example: Get coordinates
    
    
    stored.coords = []
    iterate_state 1, (all), stored.coords.append([x,y,z])
    

### Example: Modify coordinates

This example shifts the selection by 10 Angstrom along the model x-axis. 
    
    
    alter_state 1, (all), x += 10.0
    rebuild
    

### Example: Rename chain

This renames chain **A** to **C**. Note the difference between the selection syntax (first argument) and the Python syntax with the quoted string (second argument). The [sort](/index.php/Sort "Sort") call is necessary if let's say there is also a chain B, and you expect chain B to appear before chain C (formerly A) in the sequence viewer. 
    
    
    alter chain A, chain="C"
    sort
    

### Example: Sequence numbering offset

Assume the residue numbering in a PDB file starts at 100, then the following updates the residue numbers to start at 1. 
    
    
    alter (chain A), resv -= 99
    sort
    

### Example: Net Charge

The following example calculates the net charge of an object. 
    
    
    stored.net_charge = 0
    iterate (all), stored.net_charge += partial_charge
    print('Net charge: ' + str(stored.net_charge))
    

### Example: Color transfer

Copy (transfer) the color from one object to another 
    
    
    stored.colors = {}
    iterate obj1, stored.colors[chain,resi,name] = color
    alter obj2, color = stored.colors.get((chain,resi,name), color)
    recolor
    

## PYMOL API
    
    
    cmd.iterate(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.iterate_state(int state, str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter_state(int state, str selection, str expression, int quiet=1, dict space=None)
    

### Arguments

  * **state** = int: state-index if positive number or any of these: 
    * **state** = 0: all states
    * **state** = -1: current state
  * **selection** = str: atom selection
  * **expression** = str: expression in valid [python syntax](http://en.wikipedia.org/wiki/Python_syntax_and_semantics)
  * **space** = dict: namespace dictionary {default: pymol namespace}
  * **atomic** = 0/1: provide atomic properties as variables if 1, or only x/y/z if 0 (in older PyMOL versions, atomic=0 gives some speed improvement) {default: 1}



### "space" argument

The **space** argument can be used to pass local objects into the expression namespace, instead of evaluating the expression in the global **pymol** module namespace. 

The b-factor list example from above but without the global [pymol.stored](/index.php/Pymol.stored "Pymol.stored") variable: 
    
    
    myspace = {'bfactors': []}
    cmd.iterate('(all)', 'bfactors.append(b)', space=myspace)
    

User defined functions can also be included in the namespace: 
    
    
    def myfunc(resi,resn,name):
        print('%s`%s/%s' % (resn ,resi, name))
    
    myspace = {'myfunc': myfunc}
    cmd.iterate('(all)', 'myfunc(resi,resn,name)', space=myspace)
    

## See Also

  * [sort](/index.php/Sort "Sort")
  * [get_model](/index.php/Get_Model "Get Model")



Retrieved from "[https://pymolwiki.org/index.php?title=Iterate&oldid=12709](https://pymolwiki.org/index.php?title=Iterate&oldid=12709)"


---

## Alter State

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **iterate** command executes a Python expression for all atoms in a selection. The local namespace exposes all atomic identifiers and properties as read-only Python variables. The global namespace by default is the **pymol** module namespace, and thus exposes objects like **cmd** and **[stored](/index.php/Pymol.stored "Pymol.stored")**. The order of iteration is that of the internal atom ordering. 

The **alter** command is equivalent to **iterate** , but provides read-write access to the exposed variables. This can for example be used to rename a chain or to assign new b-factors. If changing identifiers which would affect atom ordering, calling [sort](/index.php/Sort "Sort") is necessary to reorder based on the new identifiers. 

The **iterate_state** command is similar to **iterate** , but iterates over the coordinates in the given state and selection, and exposes **x, y, z** in addition to the atomic identifiers and properties. 

The **alter_state** command is equivalent to **iterate_state** , but allows modification of **x, y, z** (and atom-state level settings in Incentive PyMOL 1.7+). 

## Contents

  * 1 Exposed Variables
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
    * 3.1 Example: Print b-factors
    * 3.2 Example: Get b-factor list
    * 3.3 Example: Get coordinates
    * 3.4 Example: Modify coordinates
    * 3.5 Example: Rename chain
    * 3.6 Example: Sequence numbering offset
    * 3.7 Example: Net Charge
    * 3.8 Example: Color transfer
  * 4 PYMOL API
    * 4.1 Arguments
    * 4.2 "space" argument
  * 5 See Also



## Exposed Variables

All commands in the **iterate** -family expose the following variables: 

  * **model** (str): the object name (appearing in the selection panel on the right) (cannot be altered)
  * **name** (str): the atom name
  * **resn** (str): the residue name
  * **oneletter** (str): **resn** translated to one letter code, like it's displayed in the sequence viewer (cannot be altered, alter **resn** instead) (PyMOL 1.8.7+)
  * **resi** (str): the residue identifier (residue number) as a string, including optional insertion code
  * **resv** (int): the residue identifier (residue number) as an integer, excluding insertion code
  * **chain** (str): the chain name
  * **alt** (str): alternate location identifier
  * **elem** (str): the chemical element
  * **q** (float): occupancy
  * **b** (float): the B Factor
  * **segi** (str): segment identifier (columns 73-76 in PDB file)
  * **type** (str: ATOM,HETATM): the atom type (PDB convention for canonical polymer residues)
  * **formal_charge** (int): the formal charge of the atom
  * **partial_charge** (float): the partial charge of the atom
  * **numeric_type**
  * **text_type** (str): automatic mol2 atom typing (Incentive PyMOL 1.4+)
  * **stereo** (str): automatic stereochemical R/S label (Incentive PyMOL 1.4+)
  * **ID** (int): PDB atom id (not guaranteed to be unique)
  * **rank** (int): atom index from original file import (not guaranteed to be unique)
  * **index** (int): internal atom index (unique per object, sensitive to [sorting](/index.php/Sort "Sort") and [removing](/index.php/Remove "Remove") of atoms, cannot be altered)
  * **vdw** (float): Van der Waals radius
  * **ss** (str): secondary structure
  * **color** (int): color index
  * **reps** (int): numeric mask of shown representations (PyMOL 1.7.4+)
  * **protons** (int): (PyMOL 1.7.4+)
  * **p** (object): property object to access user-defined properties like **p.abc** (Incentive PyMOL 1.7+)
  * **s** (object): settings object to access settings, e.g. **s.sphere_scale** (Incentive PyMOL 1.7+, Open-Source PyMOL 1.8.1+)
  * **label** (str): see [label](/index.php/Label "Label")
  * **geom** (int): 1=single, 2=linear, 3=planar, 4=tetrahedral (see [set_geometry](/index.php?title=Set_geometry&action=edit&redlink=1 "Set geometry \(page does not exist\)"))
  * **valence** (int): expected number of bonds (?)
  * **flags** (int): bitmask, see [flag](/index.php/Flag "Flag")
  * **cartoon** (int): see [cartoon](/index.php/Cartoon "Cartoon")



The **iterate_state** and **alter_state** commands in addition expose: 

  * **x** (float): x-coordinate
  * **y** (float): y-coordinate
  * **z** (float): z-coordinate
  * **state** (int): (PyMOL 1.7.0+)



## Usage
    
    
    iterate (selection), expression
    
    
    
    iterate_state state, (selection), expression
    
    
    
    alter (selection), expression
    
    
    
    alter_state state, (selection), expression
    

### Arguments

  * **state** = int: object state, -1 for current state, 0 for all states
  * **selection** = str: atom selection
  * **expression** = str: expression in Python language



## Examples

### Example: Print b-factors

The following prints the atom names (**name** variable) and b-factors (**b** variable) for residue #1. 
    
    
    iterate (resi 1), print(name + " %.2f" % b)
    

### Example: Get b-factor list

The following example fills a list, **[stored](/index.php/Pymol.stored "Pymol.stored").bfactors** with the f-factors (**b** variable) from residue #1. 
    
    
    stored.bfactors = []
    iterate (resi 1), stored.bfactors.append(b)
    print(stored.bfactors)
    

### Example: Get coordinates
    
    
    stored.coords = []
    iterate_state 1, (all), stored.coords.append([x,y,z])
    

### Example: Modify coordinates

This example shifts the selection by 10 Angstrom along the model x-axis. 
    
    
    alter_state 1, (all), x += 10.0
    rebuild
    

### Example: Rename chain

This renames chain **A** to **C**. Note the difference between the selection syntax (first argument) and the Python syntax with the quoted string (second argument). The [sort](/index.php/Sort "Sort") call is necessary if let's say there is also a chain B, and you expect chain B to appear before chain C (formerly A) in the sequence viewer. 
    
    
    alter chain A, chain="C"
    sort
    

### Example: Sequence numbering offset

Assume the residue numbering in a PDB file starts at 100, then the following updates the residue numbers to start at 1. 
    
    
    alter (chain A), resv -= 99
    sort
    

### Example: Net Charge

The following example calculates the net charge of an object. 
    
    
    stored.net_charge = 0
    iterate (all), stored.net_charge += partial_charge
    print('Net charge: ' + str(stored.net_charge))
    

### Example: Color transfer

Copy (transfer) the color from one object to another 
    
    
    stored.colors = {}
    iterate obj1, stored.colors[chain,resi,name] = color
    alter obj2, color = stored.colors.get((chain,resi,name), color)
    recolor
    

## PYMOL API
    
    
    cmd.iterate(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.iterate_state(int state, str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter_state(int state, str selection, str expression, int quiet=1, dict space=None)
    

### Arguments

  * **state** = int: state-index if positive number or any of these: 
    * **state** = 0: all states
    * **state** = -1: current state
  * **selection** = str: atom selection
  * **expression** = str: expression in valid [python syntax](http://en.wikipedia.org/wiki/Python_syntax_and_semantics)
  * **space** = dict: namespace dictionary {default: pymol namespace}
  * **atomic** = 0/1: provide atomic properties as variables if 1, or only x/y/z if 0 (in older PyMOL versions, atomic=0 gives some speed improvement) {default: 1}



### "space" argument

The **space** argument can be used to pass local objects into the expression namespace, instead of evaluating the expression in the global **pymol** module namespace. 

The b-factor list example from above but without the global [pymol.stored](/index.php/Pymol.stored "Pymol.stored") variable: 
    
    
    myspace = {'bfactors': []}
    cmd.iterate('(all)', 'bfactors.append(b)', space=myspace)
    

User defined functions can also be included in the namespace: 
    
    
    def myfunc(resi,resn,name):
        print('%s`%s/%s' % (resn ,resi, name))
    
    myspace = {'myfunc': myfunc}
    cmd.iterate('(all)', 'myfunc(resi,resn,name)', space=myspace)
    

## See Also

  * [sort](/index.php/Sort "Sort")
  * [get_model](/index.php/Get_Model "Get Model")



Retrieved from "[https://pymolwiki.org/index.php?title=Iterate&oldid=12709](https://pymolwiki.org/index.php?title=Iterate&oldid=12709)"


---

## Attach

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**attach** adds a single atom onto the picked atom. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 ARGUMENTS
  * 4 EXAMPLE
  * 5 NOTES



## USAGE
    
    
    attach element, geometry, valence
    

## PYMOL API
    
    
    cmd.attach( element, geometry, valence )
    

## ARGUMENTS

  * element = string: element name (symbol) of the added atom. Case-sensitive.


  * geometry = int: geometry (steric number) of added atom.


  * valence = int: valence of the added atom.



## EXAMPLE
    
    
    attach O, 1, 1 # adds an Oxygen atom to picked atom.
    

## NOTES

Immature functionality. See code for details. 

Retrieved from "[https://pymolwiki.org/index.php?title=Attach&oldid=13467](https://pymolwiki.org/index.php?title=Attach&oldid=13467)"


---

## Backward

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**backward** moves the movie back one frame. 

### USAGE
    
    
    backward
    

### PYMOL API
    
    
    cmd.backward()
    

### SEE ALSO
    
    
    [Mset](/index.php/Mset "Mset"), [Forward](/index.php/Forward "Forward"), [Rewind](/index.php/Rewind "Rewind"), [Category:Movies](/index.php/Category:Movies "Category:Movies")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Backward&oldid=7436](https://pymolwiki.org/index.php?title=Backward&oldid=7436)"


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

## Bond

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**bond** creates a new bond between two [selections](/index.php/Selection "Selection"), each of which should contain one atom. You can easily create a new bond by selecting two atoms, each with the CTRL-MIDDLE-MOUSE-BUTTON and typing "bond" on the command line. 

In order to visualize the newly created bond, [lines](/index.php/Lines "Lines") or [sticks](/index.php/Sticks "Sticks") have to be [shown](/index.php/Show "Show") for both atoms. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    bond [atom1, atom2 [,order]]
    

atom1 and atom2 must either be selections (created by the [select](/index.php/Select "Select") command, for example) or chosen using the PkAt mouse action (Ctrl+Middle mouse button). order is an optional integer (default is 1). 

Example usage: 
    
    
    bond resi 1 and name NE2, resi 10 and name ND1
    

### PYMOL API
    
    
    cmd.bond(string atom1, string atom2)
    

where atom1 and atom2 are selections. 

### NOTES

The atoms must both be within the same object. 

The default behavior is to create a bond between the (lb) and (rb) selections. 

### SEE ALSO

[Unbond](/index.php/Unbond "Unbond"), [Fuse](/index.php/Fuse "Fuse"), [Attach](/index.php/Attach "Attach"), [Replace](/index.php/Replace "Replace"), [Remove_picked](/index.php/Remove_picked "Remove picked")

Retrieved from "[https://pymolwiki.org/index.php?title=Bond&oldid=13579](https://pymolwiki.org/index.php?title=Bond&oldid=13579)"


---

## Button

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**button** can be used to redefine what the mouse buttons do. 

### USAGE
    
    
    button <button>,<modifier>,<action>
    

### PYMOL API
    
    
    cmd.button( string button, string modifier, string action )
    

### NOTES

button: 

  * L
  * M
  * R
  * S



modifers: 

  * None
  * Shft
  * Ctrl
  * CtSh



actions: 

  * Rota
  * Move
  * MovZ
  * Clip
  * RotZ
  * ClpN
  * ClpF
  * lb
  * mb
  * rb
  * +lb
  * +lbX
  * -lbX
  * +mb
  * +rb
  * PkAt
  * PkBd
  * RotF
  * TorF
  * MovF
  * Orig
  * Cent



Retrieved from "[https://pymolwiki.org/index.php?title=Button&oldid=13578](https://pymolwiki.org/index.php?title=Button&oldid=13578)"


---

## Cache

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Cache manages storage of precomputed results, such as molecular surfaces. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Notes
  * 4 API



# Usage
    
    
    cache action [, scenes [, state ]]
    

where, the arguments are: 

**action**

    

    string: enable, disable, read_only, clear, or optimize

**scenes**

    

    string: a space-separated list of scene names (default: _)_

**state**

    

    integer: state index (default: -1)

# Examples
    
    
    cache enable
    cache optimize
    cache optimize, F1 F2 F5
    

# Notes

"cache optimize" will iterate through the list of scenes provided (or all defined scenes), compute any missing surfaces, and store them in the cache for later reuse. 

# API
    
    
    cmd.cache(string action, string scenes, int state, int quiet)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cache&oldid=9091](https://pymolwiki.org/index.php?title=Cache&oldid=9091)"


---

## Capture

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Capture** obtains the entire openGL window. 

### USAGE
    
    
    capture
    

### PYMOL API
    
    
    cmd.capture
    

## Notes

It's helpful to [save] the image once you've captured it, 
    
    
    capture
    save img.png
    

Retrieved from "[https://pymolwiki.org/index.php?title=Capture&oldid=7441](https://pymolwiki.org/index.php?title=Capture&oldid=7441)"


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

## Cd

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## DESCRIPTION

Changes PyMOL's current working directory. This useful command allows the user to set the current working directory without having to leave PyMOL. 

See examples below for using spaces in pathnames on Windows. 

## USAGE
    
    
    cd DIR_NAME
    

where **DIR_NAME** is some directory. 

## EXAMPLES
    
    
    # go to my home directory
    cd ~
    
    # go to /tmp
    cd /tmp
    
    # examples of spaces on Windows
    cd \"program files"
    
    cd "\program files"
    
    cd \program?files
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cd&oldid=8498](https://pymolwiki.org/index.php?title=Cd&oldid=8498)"


---

## Center

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**center** translates the window, the clipping slab, and the origin to a point centered within the atom selection. 

## Contents

  * 1 PYMOL API
  * 2 NOTES
  * 3 User Example
    * 3.1 SEE ALSO



### PYMOL API
    
    
    cmd.center( string selection, int state = 0, int origin = 1 )
    

### NOTES

  * state = 0 (default) use all coordinate states
  * state = -1 use only coordinates for the current state
  * state > 0 use coordinates for a specific state


  * origin = 1 (default) move the origin
  * origin = 0 leave the origin unchanged



## User Example

  * Center around any given point


    
    
    # define the point as x=1.0, y=2.0, z=3.0; replace [1.0, 2.0, 3.0] with your coordinate.
    origin position=[1.0,2.0,3.0]
    # center on it
    center origin
    

### SEE ALSO

[Origin](/index.php/Origin "Origin"), [Orient](/index.php/Orient "Orient"), [Zoom](/index.php/Zoom "Zoom")

Retrieved from "[https://pymolwiki.org/index.php?title=Center&oldid=7442](https://pymolwiki.org/index.php?title=Center&oldid=7442)"


---

## Centerofmass

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The centerofmass command calculates the center of mass for an atom selection. It considers atom mass and occupancy. 

_New in PyMOL 1.7.2. See the[center_of_mass](/index.php/Center_of_mass "Center of mass") script for older PyMOL versions._

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Python API
  * 5 See Also



## Usage
    
    
    centerofmass [ selection [, state ]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **state** = int: object state, -1 for current state, 0 for all states {default: -1}



## Example
    
    
    PyMOL> fetch 1ubq, async=0
    PyMOL> centerofmass
     Center of Mass: [  30.004,  28.522,  14.701]
    

## Python API

`cmd.centerofmass()` returns a list of 3 floats. 

## See Also

  * [get_extent](/index.php/Get_extent "Get extent")
  * [center of mass](/index.php/Center_of_mass "Center of mass") (script)
  * [COM](/index.php/COM "COM") (script)



Retrieved from "[https://pymolwiki.org/index.php?title=Centerofmass&oldid=12458](https://pymolwiki.org/index.php?title=Centerofmass&oldid=12458)"


---

## Clip

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**clip** alters the near and far clipping planes 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



### USAGE
    
    
    clip {near|far|move|slab|atoms}, distance [,selection [,state ]]
    

### EXAMPLES
    
    
    clip near, -5           # moves near plane away from you by 5 A
    clip far, 10            # moves far plane towards you by 10 A
    clip move, -5           # moves the slab away from you by 5 A
    clip slab, 20           # sets slab thickness to 20 A
    clip slab, 10, resi 11  # clip 10 A slab about residue 11
    clip atoms, 5, pept     # clip atoms in "pept" with a 5 A buffer
                            # about their current camera positions
    

## PYMOL API
    
    
    cmd.clip( string mode, float distance, string selection = None)
    

## SEE ALSO

  * [Zoom](/index.php/Zoom "Zoom")
  * [Reset](/index.php/Reset "Reset")



Retrieved from "[https://pymolwiki.org/index.php?title=Clip&oldid=7443](https://pymolwiki.org/index.php?title=Clip&oldid=7443)"


---

## Cls

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**cls** clears the output buffer. 

  * [![Messy output buffer](/images/f/f9/No_cls.png)](/index.php/File:No_cls.png "Messy output buffer")

Messy output buffer 

  * [![Cleared output buffer](/images/9/9d/Cls.png)](/index.php/File:Cls.png "Cleared output buffer")

Cleared output buffer 




### USAGE
    
    
    cls
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cls&oldid=7444](https://pymolwiki.org/index.php?title=Cls&oldid=7444)"


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

## Commands

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Commands](/index.php/Category:Commands "Category:Commands") are functions in pymol that are used for viewing, manipulating, and storing molecules. 

You can get a list of some common commands by typing "help commands" into the pymol command line. Typing this in to version 1.1r1 returns the following: 

` `

`
    
    
       INPUT/OUTPUT  load      save      delete    quit
       VIEW          turn      move      clip      rock
                     show      hide      enable    disable
                     reset     refresh   rebuild   
                     zoom      origin    orient   
                     view      get_view  set_view
       MOVIES        mplay     mstop     mset      mdo
                     mpng      mmatrix   frame
                     rewind    middle    ending
                     forward   backward
       IMAGING       png       mpng
       RAY TRACING   ray       
       MAPS          isomesh   isodot
       DISPLAY       cls       viewport  splash    
       SELECTIONS    select    mask   
       SETTINGS      set       button
       ATOMS         alter     alter_state 
       EDITING       create    replace   remove    h_fill   remove_picked
                     edit      bond      unbond    h_add    fuse       
                     undo      redo      protect   cycle_valence  attach
       FITTING       fit       rms       rms_cur   pair_fit  
                     intra_fit intra_rms intra_rms_cur   
       COLORS        color     set_color
       HELP          help      commands
       DISTANCES     dist      
       STEREO        stereo
       SYMMETRY      symexp
       SCRIPTS       @         run
       LANGUAGE      alias     extend
    

```

``

Retrieved from "[https://pymolwiki.org/index.php?title=Commands&oldid=6655](https://pymolwiki.org/index.php?title=Commands&oldid=6655)"


---

## Config Mouse

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Config_mouse** allows one to change the context the mouse; say for example if it has two button instead of three (so the user doesn't have to set the 2-Button Viewing Mode). 

## USAGE

If you have a 2-Button mouse and don't want to have to set it every time, add the following line to your **pymolrc** file. 
    
    
    config_mouse two_button
    
    # put the mouse in the cool-new movie-enabling mode
    config_mouse three_button_motions
    

  


## Three-Button Motions Mode

The [Three_Button_Motions](/index.php/Three_Button_Motions "Three Button Motions") mouse configuration mode has been added to PyMOL. To use it, first prepare the mouse mode: 
    
    
    # put the mouse in the cool-new movie-enabling mode
    config_mouse three_button_motions
    

and then click on the mouse-panel (lower right hand side of PyMOL window) until you see a boxed "M" appear in the openGL menu window above. (It should appear after the "C" for Color). 

## See Also

[Three_Button_Motions](/index.php/Three_Button_Motions "Three Button Motions")

Retrieved from "[https://pymolwiki.org/index.php?title=Config_Mouse&oldid=7446](https://pymolwiki.org/index.php?title=Config_Mouse&oldid=7446)"


---

## Copy

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**copy** creates a new object that is an identical copy of an existing object 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 SEE ALSO
  * 4 User Comments/Examples



### USAGE
    
    
    copy target, source
    

### PYMOL API
    
    
    cmd.copy(string target,string source)
    

### SEE ALSO

  * [create](/index.php/Create "Create")



### User Comments/Examples
    
    
    ## will copy the object "trna" into six new objects with a number suffic
    s = range(6)
    for x in s:
    	cmd.copy("trna%s" %x, "trna")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Copy&oldid=11459](https://pymolwiki.org/index.php?title=Copy&oldid=11459)"


---

## Count Atoms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**count_atoms** returns a count of atoms in a selection. 

### USAGE
    
    
    count_atoms (selection)
    

### PYMOL API
    
    
    cmd.count_atoms(string selection)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Count_Atoms&oldid=7448](https://pymolwiki.org/index.php?title=Count_Atoms&oldid=7448)"


---

## Count Frames

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**count_frames** is an API-only function which returns the number of frames defined for the PyMOL movie. 

### PYMOL API
    
    
    cmd.count_frames()
    

## SEE ALSO

[Frame](/index.php/Frame "Frame"), [Count_States](/index.php/Count_States "Count States")

Retrieved from "[https://pymolwiki.org/index.php?title=Count_Frames&oldid=7449](https://pymolwiki.org/index.php?title=Count_Frames&oldid=7449)"


---

## Count States

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**count_states** is an API-only function which returns the number of states in the selection. 

### PYMOL API
    
    
    cmd.count_states(string selection="(all)")
    

### SEE ALSO

[Frame](/index.php/Frame "Frame")

Retrieved from "[https://pymolwiki.org/index.php?title=Count_States&oldid=7450](https://pymolwiki.org/index.php?title=Count_States&oldid=7450)"


---

## Create

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**create** creates a new molecule object from a selection. It can also be used to create states in an existing object. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 User Comments/Examples
  * 5 EXAMPLES
  * 6 SEE ALSO



### USAGE
    
    
    create name, (selection) [, source_state [, target_state ]]
    

  * name = object to create (or modify)
  * selection = atoms to include in the new object
  * source_state (default: 0 - copy all states)
  * target_state = int: -1 appends after last state {default: 0}



### PYMOL API
    
    
    cmd.create(string name, string selection, int source_state, int target_state,int discrete)
    

### NOTES

If the source and target states are zero (default), all states will be copied. Otherwise, only the indicated states will be copied. 

Target state -1 will append new states. 

### User Comments/Examples

### EXAMPLES
    
    
    # merging two multi-state objects with the create command
    # obj1 has 25 states and obj2 has 40 states 
    
    load multi_state_obj1, obj1
    load multi_state_obj2, obj2
    
    # merge the two multi-state objects back to back into obj3
    # 0 indicates every states from obj1 should be included
    # 1 indicates that the first state of obj3 is the first state from obj1 etc.
    
    create obj3, obj1, 0, 1
    
    # as obj1 has 25 states, the first state of obj2 must be state 26 to create a multi-state object of obj1 and obj2 back-to-back 
    # 26 indicates that the 26th state of obj3 is state nr. 1 of obj2 etc. 
    
    create obj3, obj2, 0, 26
    
    # Obj3 now has 65 states (25 from obj1 and 40 from obj2)
    # save obj3 as a multi-state pdb file
    # use state=0 to get all states into one pdb file
    
    save obj3.pdb, obj3, state=0
    

### SEE ALSO

[Fuse](/index.php/Fuse "Fuse"), [load](/index.php/Load "Load"), [copy](/index.php/Copy "Copy"), [Extract](/index.php/Extract "Extract")

Retrieved from "[https://pymolwiki.org/index.php?title=Create&oldid=12908](https://pymolwiki.org/index.php?title=Create&oldid=12908)"


---

## Cycle Valence

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**cycle_valence** cycles the valence on the currently selected bond. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    cycle_valence [ h_fill ]
    

### PYMOL API
    
    
    cmd.cycle_valence(int h_fill)
    

### EXAMPLES
    
    
    cycle_valence
    cycle_valence 0
    

### NOTES

If the h_fill flag is true, hydrogens will be added or removed to satisfy valence requirements. 

This function is usually connected to the DELETE key or "CTRL-W". (Try CTRL-W first.) 

There are two distinct ways to cycle a bond valence in PyMOL. First, if you start by clicking the Builder button you don't need to worry about Editing Mode, PyMOL will take care of that for you. After clicking Builder, click Cycle. In green text you should see, "Pick bonds to set as Cycle bond..." now just use the LEFT mouse button to click on a bond (not an atom). If you repeatedly click the bond, PyMOL will cycle the valence. The second method requires PyMOL to be in Editing Mode. To ensure you're in Editing Mode, please either click Mouse > 3 Button Editing. Now, using the mouse please CTRL-RIGHT-click on a bond (not an atom). The white cuff appears with a number (the angle of the substituents). Last, type CTRL-w to cycle the bond. 

### SEE ALSO

[remove_picked](/index.php/Remove_picked "Remove picked"), [attach](/index.php/Attach "Attach"), [replace](/index.php/Replace "Replace"), [fuse](/index.php/Fuse "Fuse"), [h_fill](/index.php?title=H_fill&action=edit&redlink=1 "H fill \(page does not exist\)") [Edit_Keys](/index.php/Edit_Keys "Edit Keys")

Retrieved from "[https://pymolwiki.org/index.php?title=Cycle_Valence&oldid=9183](https://pymolwiki.org/index.php?title=Cycle_Valence&oldid=9183)"


---

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

## Delete states

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**delete_states** removes states from a multi-state object like a trajectory. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    delete_states name, 1 2 3  # delete states 1, 2, and 3
    delete_states name, 1-3 9-13  # delete states 1 through 3 and 9 through 13
    

name = name of object or name expression (wildcard supported) 

### PYMOL API
    
    
    delete_states(name: str, states: str) -> None
    
    
    
    cmd.delete_states(string name = object-name, string states = states string)
    

### EXAMPLES
    
    
       delete_state 1nmr, 1-5     # delete states 1 to 5 from 1nmr
       delete_state *, 1-3 10-40  # deletes states 1 to 3 and 10 to 40 from all applicable objects
    

Note that special care needs to be taken to ensure that the object or selection name does not contain any quotes when passed as an argument. 

Note This function currently only applies to non-discrete multistate molecular objects. 

### SEE ALSO

[Remove](/index.php/Remove "Remove") [Delete](/index.php/Delete "Delete")

Retrieved from "[https://pymolwiki.org/index.php?title=Delete_states&oldid=13577](https://pymolwiki.org/index.php?title=Delete_states&oldid=13577)"


---

## Deprotect

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**deprotect** reverses the effect of the [Protect](/index.php/Protect "Protect") command. 

### USAGE
    
    
    deprotect (selection)
    

### PYMOL API
    
    
     cmd.deprotect(string selection="(all)")
    

### SEE ALSO
    
    
    [protect](/index.php/Protect "Protect"), [Mask](/index.php/Mask "Mask"), [Unmask](/index.php/Unmask "Unmask")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Deprotect&oldid=12875](https://pymolwiki.org/index.php?title=Deprotect&oldid=12875)"


---

## Deselect

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**deselect** disables any and all visible selections 

### USAGE
    
    
    deselect
    

### PYMOL API
    
    
    cmd.deselect()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Deselect&oldid=7455](https://pymolwiki.org/index.php?title=Deselect&oldid=7455)"


---

## Disable

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **disable** command toggles off the display of all currently visible representations of an object. It is the equivalent of deselecting the object in the list in the top panel of the [Internal GUI](/index.php/Internal_gui "Internal gui"). 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    disable [name]
    

name
    the name of an object or a named selection (default="all")

### PYMOL API
    
    
    cmd.disable( string name='all' )
    

### EXAMPLES
    
    
    disable my_object
    disable (my_object1 or my_object2)
    disable my_object*
    disable
    

### SEE ALSO

[Show](/index.php/Show "Show"), [Hide](/index.php/Hide "Hide"), [Enable](/index.php/Enable "Enable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates")

Retrieved from "[https://pymolwiki.org/index.php?title=Disable&oldid=11561](https://pymolwiki.org/index.php?title=Disable&oldid=11561)"


---

## Distance

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Distance
  * 2 Usage
  * 3 PYMOL API
    * 3.1 EXAMPLES
  * 4 See Also



# Distance

**distance** creates a new distance object between two selections. It will display all distances within the cutoff. Distance is also used to make hydrogen bonds. Calling distance without arguments will show distances between selections (pk1) and (pk1), which can be set in editing mode or using the PkAt mouse action (usually CTRL-middle-click).  
  
Note: For interactive use, the **measurement** [wizard](/index.php/Wizard "Wizard") (from the PyMOL menu) makes measuring distances easier than using the distance command.  
If you want to measure distance and avoid creating a distance object, use [Get Distance](/index.php/Get_Distance "Get Distance") or [Distancetoatom](/index.php/Distancetoatom "Distancetoatom").  
|  [![Example of distance. See example #1.](/images/e/ea/Dist_ex1.png)](/index.php/File:Dist_ex1.png "Example of distance. See example #1.") |  [![Example of distance. See example #2.](/images/4/4e/Dist_ex2.png)](/index.php/File:Dist_ex2.png "Example of distance. See example #2.") |   
---|---|---|---  
  
# Usage
    
    
    distance [ name [, selection1 [, selection2 [, cutoff [, mode ]]]]]
    

**name**

    

    string: name of the distance object to create

**selection1**

    

    string: first atom selection

**selection2**

    

    string: second atom selection

**cutoff**

    

    float: longest distance to show

**mode**

    

    0: all interatomic distances
    1: only bond distances
    2: only show polar contact distances
    3: like mode=0, but use [distance_exclusion](/index.php?title=Distance_exclusion&action=edit&redlink=1 "Distance exclusion \(page does not exist\)") setting
    4: distance between centroids (_new in 1.8.2_)
    5: pi-pi and pi-cation interactions
    6: pi-pi interactions
    7: pi-cation interactions
    8: like mode=3, but cutoff is the ratio between distance and sum of VDW radii

    

    modes 5 to 8 are available only in incentive PyMOL.

# PYMOL API
    
    
    cmd.distance( string name, string selection1, string selection2,
                  float cutoff, int mode )
       # returns the average distance between all atoms/frames
    

### EXAMPLES

  * Get and show the distance from residue 10's alpha carbon to residue 40's alpha carbon in 1ESR:


    
    
    # make the first residue 0.
    zero_residues 1esr, 0
    distance i. 10 and n . CA, i. 40 and n. CA
    

  * Get and show the distance from residue 10's alpha carbon to residue 35-42's alpha carbon in 1ESR:


    
    
    # make the first residue 0.
    zero_residues 1esr, 0
    distance i. 10 and n . CA, i. 35-42 and n. CA
    

  * This neat example shows how to create distance measurements from an atom in a molecule to all other atoms in the molecule (since PyMol supports wildcards).


    
    
     cmd.distance("(/mol1///1/C)","(/mol1///2/C*)")
    

or written without the PyMolScript code, 
    
    
    dist /mol1///1/C, /mol1///2/C*
    

  * Create multiple distance objects


    
    
    for at1 in cmd.index("resi 10"): \
       for at2 in cmd.index("resi 11"): \
           cmd.distance(None, "%s`%d"%at1, "%s`%d"%at2)
    
    
    
    distance (selection1), (selection2)
    
    # example
    dist i. 158 and n. CA, i. 160 and n. CA 
    
    distance mydist, 14/CA, 29/CA
    distance hbonds, all, all, 3.2, mode=2
    

# See Also

  * [Get Distance](/index.php/Get_Distance "Get Distance") # Avoid creating an object
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom") # Automated script for distances to a point
  * [Measure_Distance](/index.php/Measure_Distance "Measure Distance") # basic script
  * [Lines](/index.php/Lines "Lines")



Retrieved from "[https://pymolwiki.org/index.php?title=Distance&oldid=13867](https://pymolwiki.org/index.php?title=Distance&oldid=13867)"


---

## Draw

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/4/4a/Draw_ex.png)](/index.php/File:Draw_ex.png)

[](/index.php/File:Draw_ex.png "Enlarge")

Showing Draw in Action

**Draw** creates an oversized and antialiased OpenGL image using the current window. It's like [Ray](/index.php/Ray "Ray") but not ray traced. Also, as now with [Ray](/index.php/Ray "Ray") the oversized images are scaled and shown in the viewer window. As Draw doesn't ray trace the shadows of the scene, it is **far** faster than ray. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Match ray tracing appearance
  * 4 Notes
  * 5 See Also



## Usage
    
    
    draw [ width [, height [, antialias ]]]
    

## Examples
    
    
    draw 1600
    

will create an 1600-pixel wide image with an aspect ratio equal to that of the current screen. 
    
    
    draw 2000, 1500, 0
    

will create a 2000 by 1500 pixel image with antialiasing disabled 
    
    
    draw 600, 400, 2
    

will create a 600 by 500 pixel image with maximum (16X) antialiasing 

## Match ray tracing appearance

Since PyMOL 1.6, all "line"-type representations can be rendered as cylinders if [shaders](/index.php?title=Use_shaders&action=edit&redlink=1 "Use shaders \(page does not exist\)") are available and all ***_as_cylinders** settings are set. Example: 
    
    
    set line_as_cylinders
    set nonbonded_as_cylinders
    set ribbon_as_cylinders
    set mesh_as_cylinders
    set dash_as_cylinders
    set render_as_cylinders
    draw 3000, 2000
    

## Notes

This is a quick alternative to ray tracing with antialiasing enabled by default. 

Note that image size and antialiasing may be limited in some cases due to OpenGL hardware limits, such as screen size. For example, high-end ATI and NVidia cards max out at 4096 x 4096. 

## See Also

  * [ray](/index.php/Ray "Ray")
  * [png](/index.php/Png "Png")



Retrieved from "[https://pymolwiki.org/index.php?title=Draw&oldid=12430](https://pymolwiki.org/index.php?title=Draw&oldid=12430)"


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

## Dump

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The dump command writes the geometry of an [isosurface](/index.php/Isosurface "Isosurface"), [isomesh](/index.php/Isomesh "Isomesh"), or [isodot](/index.php/Isodot "Isodot") object to a simple text file. Each line contains one [vertex](https://en.wikipedia.org/wiki/Vertex_\(geometry\)). 

For [surface](/index.php/Isosurface "Isosurface") objects, XYZ coordinates and the [normal](https://en.wikipedia.org/wiki/Normal_\(geometry\)) are exported. Three lines make one triangle (like `GL_TRIANGLES`). 

For [mesh](/index.php/Isomesh "Isomesh") objects, XYZ coordinates are exported (no normals). The vertices form line strips (like `GL_LINE_STRIP`), a blank line starts a new strip. 

For [dot](/index.php/Isodot "Isodot") objects, XYZ coordinates are exported. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    dump filename, object
    

## Arguments

  * filename = str: file that will be written
  * object = str: object name



## Example
    
    
    fetch 1ubq, mymap, type=2fofc, async=0
    
    isosurface mysurface, mymap
    dump surfacegeometry.txt, mysurface
    
    isomesh mymesh, mymap
    dump meshgeometry.txt, mymesh
    
    isodot mydot, mymap
    dump dotgeometry.txt, mydot
    

## See Also

  * [COLLADA](/index.php/COLLADA "COLLADA") export



Retrieved from "[https://pymolwiki.org/index.php?title=Dump&oldid=13192](https://pymolwiki.org/index.php?title=Dump&oldid=13192)"


---

## Edit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**edit** picks an atom or bond for editing. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    edit (selection) [ ,(selection) ]
    

### PYMOL API
    
    
    cmd.edit( string selection  [ ,string selection ] )
    

### NOTES

If only one selection is provided, an atom is picked. If two selections are provided, the bond between them is picked (if one exists). 

### SEE ALSO

[unpick](/index.php/Unpick "Unpick"), [remove_picked](/index.php/Remove_picked "Remove picked"), [cycle_valence](/index.php/Cycle_valence "Cycle valence"), [torsion ](/index.php/Torsion "Torsion")

Retrieved from "[https://pymolwiki.org/index.php?title=Edit&oldid=8379](https://pymolwiki.org/index.php?title=Edit&oldid=8379)"


---

## Edit Keys

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

### EDITING KEYS
    
    
    These are defaults, which can be redefined.  Note that while
    entering text on the command line, some of these control keys take on
    text editing functions instead (CTRL - A, E, and K, and DELETE), so
    you should clear the command line before trying to edit atoms.
    

#### ATOM REPLACEMENT
    
    
      CTRL-C    Replace picked atom with carbon   (C)
      CTRL-N    Replace picked atom with nitrogen (N)
      CTRL-O    Replace picked atom with oxygen   (O)
      CTRL-S    Replace picked atom with sulpher  (S)
      CTRL-G    Replace picked atom with hydrogen (H)
      CTRL-F    Replace picked atom with fluorene (F)
      CTRL-L    Replace picked atom with chlorine (Cl)
      CTRL-B    Replace picked atom with bromine  (Br)
      CTRL-I    Replace picked atom with iodine   (I)
    

#### ATOM MODIFICATION
    
    
      CTRL-J    Set charge on picked atom to -1
      CTRL-K    Set charge on picked atom to +1
      CTRL-D    Remove atom or bond (DELETE works too).
      CTRL-Y    Add a hydrogen to the current atom
      CTRL-R    Adjust hydrogens on atom/bond to match valence.
      CTRL-E    Inverts the picked stereo center, but you must first
                indicate the constant portions with the (lb) and (rb)
                selections.
    
      CTRL-T    Connect atoms in the (lb) and (rb) selections.
      CTRL-W    Cycle the bond valence on the picked bond.
    

UNDO and REDO of conformational changes (not atom changes!) 
    
    
      CTRL-Z    undo the previous conformational change.
                (you can not currently undo atom modifications).
      CTRL-A    redo the previous conformational change.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Edit_Keys&oldid=6700](https://pymolwiki.org/index.php?title=Edit_Keys&oldid=6700)"


---

## Enable

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **enable** command toggles on the display of all currently visible representations of an object. It is the equivalent of selecting the object in the list at the top of the [Internal GUI](/index.php/Internal_gui "Internal gui"). 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    enable [name [, parents]]
    

name
    the name of an object or a named selection (default="all")

parents
    when set to 1, enables the named object and all parent objects of the named object (recursively), so it won't be hidden by having a parent object (or [group](/index.php/Group "Group")) that is [disabled](/index.php/Disable "Disable") (default=0)

### PYMOL API
    
    
    cmd.enable( string name='all', int parents=0 )
    

### EXAMPLES
    
    
    enable my_object
    enable (my_object1 or my_object2)
    enable my_object*
    enable
    enable my_object, parents=1
    

### SEE ALSO

[Show](/index.php/Show "Show"), [Hide](/index.php/Hide "Hide"), [Disable](/index.php/Disable "Disable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates")

Retrieved from "[https://pymolwiki.org/index.php?title=Enable&oldid=11562](https://pymolwiki.org/index.php?title=Enable&oldid=11562)"


---

## Ending

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**ending** goes to the end of the movie. 

### USAGE
    
    
    ending
    

### PYMOL API
    
    
    cmd.ending()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ending&oldid=7468](https://pymolwiki.org/index.php?title=Ending&oldid=7468)"


---

## Fab

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**fab** builds peptide entities from sequence. The sequence must be specified in one-letter code. Several fragments will be created if the sequence contains spaces (chain breaks). To specify chain-id and residue number, an advanced syntax pattern like _chain/resi/_ can be put before each sequence fragment. 

Similar functionality is also provided by the graphical [Builder](/index.php/Builder "Builder") in "Protein" mode. 

_New in PyMOL version 1.2_

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLE
  * 4 SEE ALSO



## USAGE
    
    
    fab input [, name [, mode [, resi [, chain [, segi [, state [, dir [, hydro [, ss [, async ]]]]]]]]]]
    

## ARGUMENTS

  * input = string: Amino acid sequence in one-letter code


  * name = string: name of object to create or modify {default: obj??}


  * mode = string: Only supported mode is "peptide" {default: peptide}


  * resi = integer: Residue number to start numbering from {default: 1}


  * chain, segi = string: Chain id and segment id to assign


  * state = integer: {default: -1}


  * dir = 0/1: 0=append to N-terminus, 1=append to C-terminus {default: 1}


  * hydro = 0/1: With or without hydrogens (BROKEN! Use [auto_remove_hydrogens](/index.php/Settings "Settings")}


  * ss = int: Secondary structure 1=alpha helix, 2=antiparallel beta, 3=parallel beta, 4=flat {default: 0}



## EXAMPLE

The two examples produce the same result: 
    
    
    # use resi and chain attributes
    fab KVRISAEL, myprot1, resi=10, chain=B, ss=2
    
    # use advanced syntax
    fab B/10/ KVRISAEL, myprot2, ss=2
    

## SEE ALSO

  * [Builder](/index.php/Builder "Builder")
  * [Peptide Sequence](/index.php/Peptide_Sequence "Peptide Sequence")
  * [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure")
  * [[1]](https://github.com/zigeuner/robert_campbell_pymol_scripts/tree/master/work_pymol/build_seq.py) by Robert Campbell



Retrieved from "[https://pymolwiki.org/index.php?title=Fab&oldid=13844](https://pymolwiki.org/index.php?title=Fab&oldid=13844)"


---

## Feedback

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  
feedback allows you to change the amount of information output by pymol. 

## Contents

  * 1 USAGE
  * 2 NOTES
  * 3 PYMOL API
  * 4 EXAMPLES



## USAGE
    
    
    feedback action,module,mask
    

  * action is one of ['set','enable','disable']
  * module is a space-separated list of strings or simply "all"
  * mask is a space-separated list of strings or simply "everything"



## NOTES

feedback alone will print a list of the available module choices 

## PYMOL API
    
    
    cmd.feedback(string action,string module,string mask)
    

## EXAMPLES
    
    
    feedback enable, all , debugging
    feedback disable, selector, warnings actions
    feedback enable, main, blather
    

Retrieved from "[https://pymolwiki.org/index.php?title=Feedback&oldid=7470](https://pymolwiki.org/index.php?title=Feedback&oldid=7470)"


---

## Fetch

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Fetch** retrieves a protein structure from the PDB and loads it into PyMOL. The PDB file is saved in [fetch_path](/index.php/Fetch_Path "Fetch Path"), which defaults to the current working directory for PyMOL. 

To download a so-called [biological assembly or biological unit](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies), use the [assembly](/index.php/Assembly "Assembly") setting or use _type=pdb1_ , _type=pdb2_ and so on. 

## Contents

  * 1 ChangeLog
  * 2 Usage
  * 3 Arguments
  * 4 Proxy Setting
  * 5 Examples
    * 5.1 Example
    * 5.2 Example 2
    * 5.3 Example 3 - pdb1
    * 5.4 Example 4 - multistate
  * 6 See Also



## ChangeLog

_Changed in PyMOL 2.3.0: Default async=0_

_New in PyMOL 1.8.6: Support type=mmtf and[fetch_type_default](/index.php/Fetch_type_default "Fetch type default") setting_

_Changed in PyMOL 1.8.0: Default type=cif_

_New in PyMOL 1.8.0: Support type=cc to download a chemical component by 3-letter code_

_New in PyMOL 1.7.4: Support type=cif_

_New in PyMOL 1.7.2: Support type=emd to download maps from[EMDataBank](http://www.emdatabank.org/)_

_New in PyMOL 1.7.0: Support type=cid and type=sid to download from[PubChem](https://pubchem.ncbi.nlm.nih.gov/)_

_New in PyMOL 1.6.0: Support 5 letter codes to download a single chain (4-letter pdb + 1-letter chain)_

_New in PyMOL 1.3: Support type=2fofc and type=fofc to download electron density maps_

## Usage
    
    
    fetch codes [, name [, state [, finish [, discrete [, multiplex [, zoom [, type [, async ]]]]]]]]
    

## Arguments

  * **codes** = str: one or more accession codes, separated by spaces
  * **name** = str: new object name {default: accession code}
  * **type** = cif|mmtf|pdb|pdb1|2fofc|fofc|emd|cid|sid|cc: file type and/or accession code type {default: negotiated by code or use [fetch_type_default](/index.php/Fetch_type_default "Fetch type default")}
  * **async** = 0/1: Download structures in the background and do not block the PyMOL command line. **For scripting, you typically need async=0.** {default: 0 since PyMOL 2.3, before that: !quiet, which means 1 for the PyMOL command language, and 0 for the Python API}



_Other arguments: See[load](/index.php/Load "Load") command._

## Proxy Setting

If your network requires a proxy server, you can specify it by 'http_proxy' and 'ftp_proxy' environmental variables. At least in Mac OS X, these values are setup automatically. Otherwise, add to your [pymolrc](/index.php/Pymolrc "Pymolrc") script: 
    
    
    import os
    os.environ["http_proxy"] = "http://username:password@proxy.server.name:port/"
    os.environ["ftp_proxy"] = os.environ["http_proxy"]
    

## Examples

### Example
    
    
    # fetch them singly
    fetch 1kao
    fetch 1ctq
    
    # fetch them at once
    fetch 1kao 1ctq
    
    # fetch them at once, load them into PyMOL all at once (synchronously)
    fetch 1kao 1ctq, async=0
    
    # multiple commands in one line is accepted
    fetch 1kao, async=0; as cartoon
    
    # Example loading a protein and its electron density map
    fetch 1cll
    fetch 1cll, type=2fofc
    # focus on residues 30-40
    map_trim *, i. 30-40, 4
    zoom i. 30-40
    

### Example 2
    
    
    # fetch PDB files and process each of them
    # using async=0, PyMOL will wait for fetch to finish before executing the next command
    fetch 1kao, async=0
    remove not (alt ''+A)
    alter all, alt=''
    save 1koa_clean.pdb,1koa
    delete 1koa
    fetch 1ctq, async=0
    remove not (alt ''+A)
    alter all, alt=''
    save 1ctq_clean.pdb,1ctq
    

### Example 3 - pdb1
    
    
    # fetch the biological assembly of 1avd
    # the assembly is composed of asymmetric units (ASUs) stored in different MODELs
    # split the biological assembly using split_state
    fetch 1avd, type=pdb1
    split_state 1avd
    

### Example 4 - multistate

This will fetch the biological assembly (type=pdb1) from the pdb, split the 60 states into separate objects (multiplex=1), and tell PyMOL to wait for all this to be completed before moving to the next command in a .pml script (async=0). 
    
    
    # fetch the biological assembly of 2buk
    fetch 2buk, type=pdb1, multiplex=1, async=0
    

## See Also

  * [fetch_path](/index.php/Fetch_Path "Fetch Path")
  * [fetch_type_default](/index.php/Fetch_type_default "Fetch type default")
  * [fetch_host](/index.php/Fetch_Host "Fetch Host")
  * [FetchLocal](/index.php/FetchLocal "FetchLocal")



Retrieved from "[https://pymolwiki.org/index.php?title=Fetch&oldid=12878](https://pymolwiki.org/index.php?title=Fetch&oldid=12878)"


---

## Fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Fit superimposes the model in the first selection on to the model in the second selection. Only _matching atoms_ in both selections will be used for the fit. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    fit mobile, target [, mobile_state [, target_state [, quiet [, matchmaker [, cutoff [, cycles [, object ]]]]]]]
    

### ARGUMENTS

  * **mobile** = string: atom selection
  * **target** = string: atom selection
  * **mobile_state** = integer: object state {default=0, all states)
  * **target_state** = integer: object state {default=0, all states)
  * **matchmaker** = integer: how to match atom pairs {default: 0} 
    * -1: assume that atoms are stored in the identical order
    * 0/1: match based on all atom identifiers (segi,chain,resn,resi,name,alt)
    * 2: match based on ID
    * 3: match based on rank
    * 4: match based on index (same as -1 ?)
  * **cutoff** = float: outlier rejection cutoff (only if cycles>0) {default: 2.0}
  * **cycles** = integer: number of cycles in outlier rejection refinement {default: 0}
  * **object** = string: name of alignment object to create {default: None}



### EXAMPLES
    
    
    fit ( mutant and name ca ), ( wildtype and name ca )
    

If atom identifiers (like segi, chain, ...) in mobile and target do not match, you need to alter them (or use [Pair_Fit](/index.php/Pair_Fit "Pair Fit") instead): 
    
    
    fetch 1a00, async=0
    extract hbaA, chain A
    extract hbaC, chain C
    # hbaA and hbaC are the same protein, but have different chain identifiers
    alter hbaC, chain='A'
    alter hbaC, segi='A'
    # now both have identical atom identifiers
    fit hbaC, hbaA
    

### SEE ALSO

[Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

Retrieved from "[https://pymolwiki.org/index.php?title=Fit&oldid=13194](https://pymolwiki.org/index.php?title=Fit&oldid=13194)"


---

## Flag

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Flag is a command to set or clear _flags_ on atom sets. A flag is just some atom-specific property. A flag is either on or off for a residue. Possible flags are: 

Flag name  | Value  | Description  | Notes   
---|---|---|---  
**focus** | 0  | Atoms of Interest (i.e. a ligand in an active site)  | Reserved for molecular modeling. Affects [sculpting](/index.php/Molecular_Sculpting "Molecular Sculpting").   
**free** | 1  | Free Atoms (free to move subject to a force-field)   
**restrain** | 2  | Restrained Atoms (typically harmonically constrained)   
**fix** | 3  | Fixed Atoms (no movement allowed)   
**exclude** | 4  | Atoms which should not be part of any simulation   
**study** | 5  |   
| 6  | Protein (`polymer.protein` selector)  | See [auto_classify_atoms](/index.php?title=Auto_classify_atoms&action=edit&redlink=1 "Auto classify atoms \(page does not exist\)") and [auto_show_classified](/index.php/Auto_show_classified "Auto show classified")  
| 7  | Nucleic acid (`polymer.nucleic` selector)   
| 8-15  | _Free for end users to manipulate_ |   
| 16-23  | Reserved for external GUIs and linked applications   
**exfoliate** | 24  | **Deprecated**. Remove surface from atoms when surfacing. Redundant with excluding those atoms from the selection in `show surface, sele` | Affects [surface](/index.php/Surface "Surface") (with [surface_mode](/index.php/Surface_mode "Surface mode")=0), [dots](/index.php/Dots "Dots") (with [trim_dots](/index.php?title=Trim_dots&action=edit&redlink=1 "Trim dots \(page does not exist\)")=on), [get_area](/index.php/Get_Area "Get Area")  
**ignore** | 25  | Ignore atoms altogether when surfacing   
**no_smooth** | 26  | Do not smooth atom position  | Affects [cartoon](/index.php/Cartoon "Cartoon")  
| 27  | Polymer  | See [auto_classify_atoms](/index.php?title=Auto_classify_atoms&action=edit&redlink=1 "Auto classify atoms \(page does not exist\)") and [auto_show_classified](/index.php/Auto_show_classified "Auto show classified") See [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra") "Chemical classes"   
| 28  | Solvent   
| 29  | Organic   
| 30  | Inorganic   
| 31  | Guide atom (e.g. CA in proteins)   
  
## Usage
    
    
    flag flag, selection [, action [, quiet ]]
    

If the [auto_indicate_flags](/index.php?title=Auto_indicate_flags&action=edit&redlink=1 "Auto indicate flags \(page does not exist\)") setting is true, then PyMOL will automatically create a selection called "indicate" which contains all atoms with that flag after applying the command. 

## Arguments

  * **flag** = int or str: Flag name or value
  * **selection** = str: atom selection
  * **action** = set|clear|reset: _Note that "reset" will set the flag on the given selection, and clear it on all other atoms_ {default: reset}



## Examples

[![](/images/1/10/Flags-ignore-exfoliate.png)](/index.php/File:Flags-ignore-exfoliate.png)

[](/index.php/File:Flags-ignore-exfoliate.png "Enlarge")

CYS residue with flag "ignore" (left) and flag "exfoliate" (right)
    
    
    fab AC
    
    # Image on the left
    flag ignore, resn CYS
    show surface
    
    # Image on the right
    flag ignore, all, clear
    flag exfoliate, resn CYS
    rebuild surface
    
    
    
    # in sculpting, ensure the newMethyl group just added doesn't move around
    flag fix, newMethyl
    
    # Introspect the flags bitmask
    iterate all, print(hex(flags))
    
    # Select atoms with "fix" flag
    select fixedatoms, flag 3
    

Retrieved from "[https://pymolwiki.org/index.php?title=Flag&oldid=13318](https://pymolwiki.org/index.php?title=Flag&oldid=13318)"


---

## Fnab

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**fnab** builds nucleic acid entities from sequence. The sequence must be specified in one-letter code. Only one fragment can be created at a type (no spaces allowed unlike the 'fab' command). 

Similar functionality is also provided by the graphical [Builder](/index.php/Builder "Builder") in "Nucleic Acid" mode. 

_New in PyMOL version 2.3_

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLE
  * 4 SEE ALSO



## USAGE
    
    
    fnab input [, name [, type [, form [, dbl_helix ]]]]
    

## ARGUMENTS

  * input = string: Nucleic Acid sequence in one-letter code


  * name = string: name of object to create or modify {default: obj??}


  * mode = string: "DNA" or "RNA" {default: "DNA"}


  * form = string: "A" or "B" {default: "B"}


  * dbl_helix = bool (0/1): flag for using double helix in DNA



## EXAMPLE
    
    
    # Create default (B-DNA) chain
    fnab ATGCGATAC
    
    # Create B-DNA chain with specific name 
    fnab ATGCGATAC, name=myDNA, mode=DNA, form=B, dbl_helix=1
    
    # Create RNA chain
    fnab AAUUUUCCG, mode=RNA
    

  * [![Example 1](/images/f/ff/Fnab_example.png)](/index.php/File:Fnab_example.png "Example 1")

Example 1 




## SEE ALSO

  * [Builder](/index.php/Builder "Builder")
  * [Fab](/index.php/Fab "Fab")



Retrieved from "[https://pymolwiki.org/index.php?title=Fnab&oldid=12854](https://pymolwiki.org/index.php?title=Fnab&oldid=12854)"


---

## Fork

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**spawn** launches a Python script in a new thread which will run concurrently with the PyMOL interpreter. It can be run in its own namespace (like a Python module, default), a local name space, or in the global namespace. 

**fork** is an alias for **spawn**. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    spawn python-script [, ( local | global | module | main | private )]
    

### PYMOL API

Not directly available. Instead, use cmd.do("spawn ..."). 

### NOTES

The default mode for spawn is "module". 

Due to an idiosyncracy in Pickle, you can not pickle objects directly created at the main level in a script run as "module", (because the pickled object becomes dependent on that module). Workaround: delegate construction to an imported module. 

The best way to spawn processes at startup is to use the -l option (see "help launching"). 

### SEE ALSO

  * [run](/index.php/Run "Run")



Retrieved from "[https://pymolwiki.org/index.php?title=Fork&oldid=9442](https://pymolwiki.org/index.php?title=Fork&oldid=9442)"


---

## Forward

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**forward** moves the movie one frame forward. 

### USAGE
    
    
    forward
    

### PYMOL API
    
    
    cmd.forward()
    

### SEE ALSO

[Mset](/index.php/Mset "Mset"), [Backward](/index.php/Backward "Backward"), [Rewind](/index.php/Rewind "Rewind"), [Category:Movies](/index.php/Category:Movies "Category:Movies")

Retrieved from "[https://pymolwiki.org/index.php?title=Forward&oldid=7475](https://pymolwiki.org/index.php?title=Forward&oldid=7475)"


---

## Fragment

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

fragment retrieves a 3D structure from the fragment library, which is currently pretty meager (just amino acids). Use the three-letter abbreviations. 

  * [![Result of "fragment trp".](/images/3/3c/Fragment_ex.png)](/index.php/File:Fragment_ex.png "Result of "fragment trp".")

Result of "fragment trp". 




## USAGE
    
    
    # generates ''name''.
    fragment name
    
    # for example, generate an alanine:
    fragment ala
    

## See also

[Modeling_and_Editing_Structures#Adding_and_using_your_own_fragments](/index.php/Modeling_and_Editing_Structures#Adding_and_using_your_own_fragments "Modeling and Editing Structures")

Retrieved from "[https://pymolwiki.org/index.php?title=Fragment&oldid=10342](https://pymolwiki.org/index.php?title=Fragment&oldid=10342)"


---

## Frame

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**frame** sets the viewer to the indicated movie frame. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    frame frame-number
    

### PYMOL API
    
    
    cmd.frame( int frame_number )
    

### NOTES

Frame numbers are 1-based 

### SEE ALSO

[Cmd count_states](/index.php/Cmd_count_states "Cmd count states")

Retrieved from "[https://pymolwiki.org/index.php?title=Frame&oldid=7477](https://pymolwiki.org/index.php?title=Frame&oldid=7477)"


---

## Full Screen

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**full_screen** enables or disables PyMOL's full screen mode. This does not work well on all platforms. 

### USAGE
    
    
    full_screen on
    full_screen off
    

### PYMOL API
    
    
    cmd.full_screen('on')
    cmd.full_screen('off')
    

Retrieved from "[https://pymolwiki.org/index.php?title=Full_Screen&oldid=7478](https://pymolwiki.org/index.php?title=Full_Screen&oldid=7478)"


---

## Full screen

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This determines if PyMOL takes up the entire desktop screen or not. 

You can also use the "F" button in the bottom-right corner of the [internal GUI](/index.php/Internal_Gui "Internal Gui") to toggle full screen. 

# Syntax
    
    
    # make PyMOL run full screen
    set full_screen, on
    
    # go back to normal mode
    set full_screen, off
    

# Notes

Didn't do anything in my Mac version of PyMOL. 

Retrieved from "[https://pymolwiki.org/index.php?title=Full_screen&oldid=9411](https://pymolwiki.org/index.php?title=Full_screen&oldid=9411)"


---

## Fuse

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**fuse** joins two objects into one by forming a bond. A copy of the object containing the first atom is moved so as to form an approximately reasonable bond with the second, and is then merged with the first object. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Notes
  * 5 See Also



## Usage
    
    
    fuse [ selection1 [, selection2 [, mode [, recolor [, move ]]]]]
    

## Arguments

  * selection1 = str: single atom selection (will be copied to object 2) {default: (pk1)}
  * selection2 = str: single atom selection {default: (pk2)}
  * mode = int: {default: 0} 
    * mode=0:
    * mode=1: adopt segi/chain/resn`resi
    * mode=2: adopt segi/chain
    * mode=3: don't move and don't create a bond, just combine into single object
  * recolor = bool: recolor C atoms to match target {default: 1}
  * move = bool: {default: 1}



## Example

Phosphorylate a serine residue: 
    
    
    fetch 1ubq, async=0
    
    fragment phosphite
    remove phosphite & elem H
    
    fuse phosphite & elem P, /1ubq///65/OG, mode=1
    delete phosphite
    
    show sticks, resi 65
    orient resi 65
    unpick
    

## Notes

Each selection must include a single atom in each object. The atoms can both be hydrogens, in which case they are eliminated, or they can both be non-hydrogens, in which case a bond is formed between the two atoms. 

## See Also

[Bond](/index.php/Bond "Bond"), [Unbond](/index.php/Unbond "Unbond"), [Attach](/index.php/Attach "Attach"), [Replace](/index.php/Replace "Replace"), [Remove_picked](/index.php/Remove_picked "Remove picked")

Retrieved from "[https://pymolwiki.org/index.php?title=Fuse&oldid=12294](https://pymolwiki.org/index.php?title=Fuse&oldid=12294)"


---

## Get

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

get returns the value of a setting. 

This command is very useful for determining the any setting(s) when writing a script. For example, with this command you can find out if the background is opaque, where the light source is, etc. 

## Contents

  * 1 USAGE
  * 2 Examples
  * 3 PYMOL API
  * 4 Notes
  * 5 SEE ALSO



# USAGE
    
    
    get name [, selection [, state ]]
    

# Examples
    
    
    get opaque_background
    
    get line_width
    

# PYMOL API
    
    
    print cmd.get(string name, string object, int state, int quiet)
    

# Notes

  * The API command will not print out and should be stored or used for comparison
  * "get" currently only works with global, per-object, and per-state settings. There is currently no way to retrieve per-atom settings.



# SEE ALSO

[Set](/index.php/Set "Set")

Retrieved from "[https://pymolwiki.org/index.php?title=Get&oldid=7467](https://pymolwiki.org/index.php?title=Get&oldid=7467)"


---

## Get Angle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_angle** returns the angle between three atoms. By default, the coordinates used are from the current state, however an alternate state identifier can be provided. 

### USAGE

get_angle atom1, atom2, atom3, [,state ] 

### EXAMPLES
    
    
    get_angle 4/n,4/c,4/ca
    get_angle 4/n,4/c,4/ca,state=4
    

### PYMOL API
    
    
    cmd.get_angle(atom1="pk1",atom2="pk2",atom3="pk3",state=0)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Angle&oldid=7466](https://pymolwiki.org/index.php?title=Get_Angle&oldid=7466)"


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

## Get chains

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This command will return a python list of chains for the given selection. 

Proteins structures often have more than one chain of residues present in a structure file. This command will fetch the names of each chain present. Because of the lack of standards, sometimes chain A will also be blank "". This isn't a problem for PyMOL, but be warned sometimes you may get back, [""]. 

# Syntax
    
    
    # using get_chains for the object or selection called objSel
    get_chains objSel
    

# Examples
    
    
    # examples
    fetch 1cll
    print "1CLL has the following chains:", cmd.get_chains("1cll")
    
    # list all chains in all proteins loaded in PyMOL:
    for x in cmd.get_names():
      for ch in cmd.get_chains(x):
        print x, " has chain ", ch
    

# See Also

[Get_Names](/index.php/Get_Names "Get Names")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_chains&oldid=7503](https://pymolwiki.org/index.php?title=Get_chains&oldid=7503)"


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

## Get Dihedral

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_dihedral** returns the dihedral angle between four atoms. By default, the coordinates used are from the current state, however an alternate state identifier can be provided. 

By convention, positive dihedral angles are right-handed (looking down the atom2-atom3 axis). 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



### USAGE
    
    
    get_dihedral atom1, atom2, atom3, atom4 [,state ]
    

### EXAMPLES
    
    
    get_dihedral 4/n,4/c,4/ca,4/cb
    get_dihedral 4/n,4/c,4/ca,4/cb,state=4
    

### PYMOL API
    
    
    cmd.get_dihedral(atom1,atom2,atom3,atom4,state=0)
    

### SEE ALSO

  * [Set_Dihedral](/index.php/Set_Dihedral "Set Dihedral")
  * [DynoPlot](/index.php/DynoPlot "DynoPlot")
  * [Displaying_Biochemical_Properties#Calculating_dihedral_angles](/index.php/Displaying_Biochemical_Properties#Calculating_dihedral_angles "Displaying Biochemical Properties")
  * [Rotamer_Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_Dihedral&oldid=12175](https://pymolwiki.org/index.php?title=Get_Dihedral&oldid=12175)"


---

## Get Distance

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_distance** returns the distance between two atoms. By default, the coordinates used are from the current state, however an alternate state identifier can be provided. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



### USAGE
    
    
    get_distance atom1, atom2, [,state ]
    

### EXAMPLES
    
    
    get_distance 4/n,4/c
    get_distance 4/n,4/c,state=4
    

### PYMOL API
    
    
    cmd.get_distance(atom1="pk1",atom2="pk2",state=0)
    

### SEE ALSO

  * [Distance](/index.php/Distance "Distance") # create a distance object
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom") # Automated script for distances to a point
  * [Measure_Distance](/index.php/Measure_Distance "Measure Distance") # basic script



Retrieved from "[https://pymolwiki.org/index.php?title=Get_Distance&oldid=11647](https://pymolwiki.org/index.php?title=Get_Distance&oldid=11647)"


---

## Get Extent

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_extent** returns the minimum and maximum XYZ coordinates of a selection as an array: 
    
    
    [ [ min-X , min-Y , min-Z ],[ max-X, max-Y , max-Z ]]

Typing this command returns the coordinates for all atoms. To return the coordinates of the default selection type: 
    
    
    get_extent sele

### PYMOL API
    
    
    cmd.get_extent(string selection="(all)", state=0 )
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Extent&oldid=7463](https://pymolwiki.org/index.php?title=Get_Extent&oldid=7463)"


---

## Get Model

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_model** returns a model object. 

## Contents

  * 1 PYMOL API
  * 2 USAGE
  * 3 NOTES
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.get_model(string "selection", integer "state" )
    

### USAGE
    
    
    cmd.get_model("chain A")
    

### NOTES

It can be useful to loop through all the atoms of a selection (rather than using the iterate command) 
    
    
    atoms = cmd.get_model("chain A")
    for at in atoms.atom:
        print("ATOM DEFINITION: "+at.model+" "\
                                 +at.chain+" "\
                                 +at.resn+" "\
                                 +str(at.resi)+" "\
                                 +at.name+" "\
                                 +str(at.index)+" "\
                                 +"%.2f " % (at.b)\
                                 +str(at.coord[0])+" "\
                                 +str(at.coord[1])+" "\
                                 +str(at.coord[2]))
    

### SEE ALSO

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Model&oldid=13576](https://pymolwiki.org/index.php?title=Get_Model&oldid=13576)"


---

## Get Names

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_names** returns a list of object and/or selection names. 

## Contents

  * 1 PYMOL API
  * 2 ARGUMENTS
  * 3 NOTES
  * 4 EXAMPLES
  * 5 SEE ALSO



### PYMOL API
    
    
    cmd.get_names(type,enabled_only,selection)
    

### ARGUMENTS

  * **type : string** determines the type of objects to be returned. It can take any of the following values: 
    * **objects** Structure objects
    * **selections** All selection
    * **all** Objects and Selections
    * **public_objects** (default)
    * **public_selections**
    * **public_nongroup_objects**
    * **public_group_objects**
    * **nongroup_objects**
    * **group_objects**
  * **enabled_only : boolean** If 1, only enabled objects are returned. Default 0 (all objects)
  * **selection**



### NOTES

The default behavior is to return only object names. 

### EXAMPLES

Multiple alignments 
    
    
    # structure align all proteins in PyMOL to the protein named "PROT".  Effectively a 
    # poor multiple method for multiple structure alignment.
    for x in cmd.get_names("all"): cealign( "PROT", x)
    

Determine whether or not an object (objName) is enabled: 
    
    
    if objName in cmd.get_names(enabled_only=1):
        print objName, "is enabled"
    

### SEE ALSO

[get_type](/index.php/Get_type "Get type"), [get_names_of_type](/index.php/Get_names_of_type "Get names of type"), [count_atoms](/index.php?title=Count_atoms&action=edit&redlink=1 "Count atoms \(page does not exist\)"), [count_states](/index.php/Count_states "Count states")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Names&oldid=12292](https://pymolwiki.org/index.php?title=Get_Names&oldid=12292)"


---

## Get names of type

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_names_of_type** returns a list of object and/or selection names. 

## Contents

  * 1 PYMOL API
  * 2 NOTES
  * 3 EXAMPLES
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.get_names_of_type(string type)
    

### NOTES

The object types are strings such as 

  * object:molecule
  * object:map
  * object:mesh
  * object:slice
  * object:surface
  * object:measurement
  * object:cgo
  * object:group
  * object:volume



### EXAMPLES

Truncate names of all molecules 
    
    
    # Get names for all molecules.
    for x in cmd.get_names_of_type("object:molecule"): cmd.set_name(x,x[:5])
    

### SEE ALSO

[get_names](/index.php/Get_names "Get names"), [get_type](/index.php/Get_type "Get type"), [count_atoms](/index.php?title=Count_atoms&action=edit&redlink=1 "Count atoms \(page does not exist\)"), [count_states](/index.php/Count_states "Count states")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_names_of_type&oldid=13185](https://pymolwiki.org/index.php?title=Get_names_of_type&oldid=13185)"


---

## Get object matrix

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Example
    * 3.1 See Also



# Overview

Pymol stores a transformation matrix for each object relative to it's initial position when loaded in. **get_object_matrix** will return a list of floats with that matrix for a named object. 

The matrix is 4X4, with the upper left 3x3 forming a rotation matrix, the fourth column and row representing pre-rotation and post-rotation translation vectors respectively, and the 16th element always being 1.0. 

According the pymol source code, this is an "unsupported command". 

# Syntax
    
    
    cmd.get_object_matrix(object, state=1)
    

# Example
    
    
    cmd.load("prot1.pdb", "prot1")
    cmd.load("prot2.pdb", "prot2")
    cmd.super("prot1", "prot2") #align prot1 to prot 2
    transformation = cmd.get_object_matrix("prot1") #translation and rotation to align the two proteins
    

## See Also

[Object_Matrix](/index.php/Object_Matrix "Object Matrix") [Transform_selection](/index.php/Transform_selection "Transform selection") [Transform_odb](/index.php/Transform_odb "Transform odb") [Matrix_copy](/index.php/Matrix_copy "Matrix copy")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_object_matrix&oldid=8449](https://pymolwiki.org/index.php?title=Get_object_matrix&oldid=8449)"


---

## Get Pdbstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_pdbstr** is an API-only function which returns a pdb corresponding to the atoms in the selection provided and that are present in the indicated state 

### PYMOL API ONLY
    
    
    cmd.get_pdbstr( string selection="all", int state=0 )
    

#### NOTES

  1. **state** is a 1-based state index for the object.
  2. if state is zero, then current state is used.



### EXAMPLES

The following example lists the residues in the selection called, **near**. 
    
    
    import cmd # -- if required
    select near, sel01 around 6
    cmd.get_pdbstr("near")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Pdbstr&oldid=7460](https://pymolwiki.org/index.php?title=Get_Pdbstr&oldid=7460)"


---

## Get Position

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Get_Position gets the 3D coordinates of the center of the viewer window.. 

## Syntax
    
    
    # print the coordinates for center of mass
    zoom
    print cmd.get_position()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Position&oldid=8374](https://pymolwiki.org/index.php?title=Get_Position&oldid=8374)"


---

## Get session

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Returns a dictionary representation the currently loaded PyMOL session.The session file (.pse) is a compressed version of its output. 

By using this API, user scripts can access many interesting properties which are otherwise inaccessible. Examples include [Pymol2glmol](/index.php/Pymol2glmol "Pymol2glmol") and [get_raw_distances](/index.php/Get_raw_distances "Get raw distances") scripts. 

WARNING: This API is undocumented API, intended for internal use. Use it only when it is necessary. 

## Contents

  * 1 Usage
  * 2 Return value
  * 3 Examples
    * 3.1 Object array
    * 3.2 Internals of a map object



## Usage
    
    
    dict cmd.get_session(names='', partial=0)
    

## Return value

The returned dictionary has following key-value pairs. Some entries are missing if using _partial=1_. 

  * **main** : An array encoding window size. For example, [640, 480]
  * **color_ext** :
  * **settings** : An array of PyMOL's global settings. They are dumped by SettingAsPyList in layer1/Setting.c. Further details will not be discussed because scripts can access these values from [cmd.get](/index.php/Get "Get") API.
  * **colors** : If you have defined color names by [Set Color](/index.php/Set_Color "Set Color"), they are described here. Default color names (red, blue, etc...) will not appear. This dict can also be obtained by [cmd.get_color_indices()](/index.php/Get_Color_Indices "Get Color Indices"). To get RGB definition of a color, get its pallet ID by cmd.get_color_index and convert it to RGB tuple by cmd.get_color_tuple.
  * **view** : Same as [cmd.get_view()](/index.php/Get_View "Get View")
  * **version** : Version number of PyMOL
  * **view_dict** :
  * **names** : This is the most complex but interesting part. Produced by ExecutiveGetNamedEntries, this array describes internal C-objects of the current session. Each element is an array of seven elements.



  
(documentation in progress... Please feel free to expand) 

## Examples

### Object array
    
    
    x = cmd.get_session('myobj', 1)['names'][0]
    x = [
        str name ('myobj'),
        int 0=object/1=selection,
        int enabled,
        list representations,
        int object_type,
        list object_data,
        str group_name,
    ]
    

### Internals of a map object

The map data as a numpy array: 
    
    
    import numpy
    from pymol import cmd
    
    cmd.map_new('mymap', 'gaussian', 0.5)
    
    mymap = cmd.get_session('mymap', 1)['names'][0]
    field_data = mymap[5][2][0][14][2]
    
    values = numpy.reshape(field_data[6], field_data[4])
    

The dumped map datastructure: 
    
    
    mymap[5][2][0] = [
        int Active,
        list Symmetry (like cmd.get_symmetry),
        list Origin,
        list Range (ExtentMax - ExtentMin),
        list Dim (?),
        list Grid (grid spacing),
        list Corner,
        list ExtentMin,
        list ExtentMax,
        int MapSource (?, 4),
        list Div (?, Dim - 1),
        list Min (?),
        list Max (?, ncells),
        list FDim (?, (ncells + 1) + [3]),
        [
            list field->dimensions (FDim[:3]),
            int field->save_points (?, 1),
            [# field->data
                int type (?, 0),
                int n_dim,
                int base_size (?, 4),
                int size (?, len * 4),
                list dim,
                list stride,
                list data,
            ],
            [# field->points
                ...
            ]   
        ],
        [None],
    ]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_session&oldid=10855](https://pymolwiki.org/index.php?title=Get_session&oldid=10855)"


---

## Get Symmetry

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_symmetry** can be used to obtain the crystal and spacegroup parameters for a molecule or map object. 

## Contents

  * 1 USAGE
  * 2 DESCRIPTION
  * 3 PYMOL API
  * 4 See Also



### USAGE
    
    
    get_symmetry object-name-or-selection
    

### DESCRIPTION

Returns a tuple containing the following 7 values: 

  * The unit cell lengths (a,b,c)
  * The unit cell angles (alpha, beta, gamma)
  * The space group name (e.g. "P 21 21 21")



### PYMOL API
    
    
    cmd.get_symmetry(string selection)
    

### See Also

  * [set_symmetry](/index.php/Set_Symmetry "Set Symmetry")
  * [symmetry_copy](/index.php?title=Symmetry_copy&action=edit&redlink=1 "Symmetry copy \(page does not exist\)")
  * [Supercell](/index.php/Supercell "Supercell")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_Symmetry&oldid=13389](https://pymolwiki.org/index.php?title=Get_Symmetry&oldid=13389)"


---

## Get Title

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_title** retrieves a text string to the state of a particular object which will be displayed when the state is active. This is useful for printing the names of objects (given a state). 

## USAGE
    
    
    get_title object,state
    

## PYMOL API
    
    
     cmd.get_title(string object, int state)
    

## Example

Print out all the object names of the ensemble of states loaded in: 
    
    
    for x in range(numStates):
      print cmd.get_title("objName", x)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Title&oldid=13528](https://pymolwiki.org/index.php?title=Get_Title&oldid=13528)"


---

## Get Type

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_type** returns a string describing the named object or selection or the string "nonexistent" if the name in unknown. 

### PYMOL API
    
    
    cmd.get_type(string object-name)
    

### NOTES

Possible return values are 

  1. "object:molecule"
  2. "object:map"
  3. "object:mesh"
  4. "object:distance"
  5. "selection"



### SEE ALSO

[Cmd get_names](/index.php/Cmd_get_names "Cmd get names")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Type&oldid=7501](https://pymolwiki.org/index.php?title=Get_Type&oldid=7501)"


---

## Get type

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

PyMOL objects have a base type. get_type will retrieve that base type name. 

# Usage
    
    
    # create a selection
    
    cmd.select("foo", "all")
    
    # print the type
    
    print cmd.get_type("foo")
    

# See Also

[get_names](/index.php/Get_names "Get names")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_type&oldid=9747](https://pymolwiki.org/index.php?title=Get_type&oldid=9747)"


---

## Get Version

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_version** returns the version string of the current PyMOL program. 

## PYMOL API
    
    
    # get the version
    cmd.get_version()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Version&oldid=7500](https://pymolwiki.org/index.php?title=Get_Version&oldid=7500)"


---

## Get View

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_view** returns and optionally prints out the current view information in a format which can be embedded into a command script and can be used in subsequent calls to **set_view**. 

If a log file is currently open, get_view will not write the view matrix to the screen unless the "output" parameter is 2. 

This command is very useful for saving the orientation of a scene for later. Authors of molecular movies may find this command very powerful. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 API USAGE
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    get_view
    

### PYMOL API
    
    
    cmd.get_view(output=1, quiet=1)
       
    my_view = cmd.get_view()
    

output control: 

  * 0 = output matrix to screen
  * 1 = don't output matrix to screen
  * 2 = force output to screen even if log file is open
  * 3 = return formatted string instead of a list



### API USAGE
    
    
    cmd.get_view(0) # zero option suppresses output (LEGACY approach)
    cmd.get_view(quiet=1) # suppresses output using PyMOL's normal "quiet" parameter.
    

### NOTES

Contents of the view matrix 

  * 0 - 8 = column-major 3x3 matrix which rotates model axes to camera axes
  * 9 - 11 = origin of rotation relative to the camera in camera space
  * 12 - 14 = origin of rotation in model space
  * 15 = front plane distance from the camera
  * 16 = rear plane distance from the camera
  * 17 = orthoscopic flag (not implemented in older versions)



### SEE ALSO

[Set View](/index.php/Set_View "Set View"), [View](/index.php/View "View"), [Model_Space_and_Camera_Space](/index.php/Model_Space_and_Camera_Space "Model Space and Camera Space")

Retrieved from "[https://pymolwiki.org/index.php?title=Get_View&oldid=13321](https://pymolwiki.org/index.php?title=Get_View&oldid=13321)"


---

## Get sasa relative

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

get_sasa_relative calculates the [relative per-residue solvent accessible surface area](https://en.wikipedia.org/wiki/Relative_accessible_surface_area) and optionally labels and colors residues. The value is relative to full exposure of the residue, calculated by removing all other residues except its two next neighbors, if present. 

The command loads a value beteween 0.0 (fully buried) and 1.0 (fully exposed) into the b-factor property, available in [iterate](/index.php/Iterate "Iterate"), [alter](/index.php/Alter "Alter") and [label](/index.php/Label "Label") as **b**. 

_New in Incentive PyMOL 1.8 and Open-Source PyMOL 2.1_

_Changed in PyMOL 2.6: Added**subsele** argument_

## Usage
    
    
    get_sasa_relative [ selection [, state [, vis [, var [, quiet [, outfile [, subsele ]]]]]]]
    

## Example
    
    
    fetch 1ubq
    get_sasa_relative polymer
    

Side-chain exposure, including C-alpha atom, using **subsele** argument (new in PyMOL 2.6): 
    
    
    get_sasa_relative polymer, subsele=sidechain guide
    

## See Also

  * [get_area](/index.php/Get_area "Get area")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_sasa_relative&oldid=13465](https://pymolwiki.org/index.php?title=Get_sasa_relative&oldid=13465)"


---

## Gradient

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

gradient creates a gradient ([field lines](http://en.wikipedia.org/wiki/Field_line)) object from a [map](/index.php?title=Map&action=edit&redlink=1 "Map \(page does not exist\)") object. 

## Usage
    
    
    gradient name, map [, minimum [, maximum [, selection
        [, buffer [, state [, carve [, source_state ]]]]]]]
    

## Arguments

  * **map** = string: the name of the map object to use.
  * **minimum, maximum** = float: minimum and maximum levels (default: full map range)
  * **selection** = string: an atom selection about which to display the mesh with an additional "buffer" (if provided).



## See Also

  * [slice](/index.php/Slice "Slice"), [ramp_new](/index.php/Ramp_new "Ramp new"), [isomesh](/index.php/Isomesh "Isomesh")
  * [APBS](/index.php/APBS "APBS")



Retrieved from "[https://pymolwiki.org/index.php?title=Gradient&oldid=10887](https://pymolwiki.org/index.php?title=Gradient&oldid=10887)"


---

## Group

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The Group command creates or updates a "group" object. The grouped objects are collected underneath a **+** sign in the object tree (see images) in the Pymol Internal Gui. 

Group is tremendously helpful with multi-state or multi-structure sessions. Wildcards work great, for example: 
    
    
    # put all of objState into the group "ensemble".
    group ensemble, objState*
    

  * [![Three EF-Hand proteins loaded into PyMOL](/images/d/d2/Group_off.png)](/index.php/File:Group_off.png "Three EF-Hand proteins loaded into PyMOL")

Three EF-Hand proteins loaded into PyMOL 

  * [![Applied the group command to the proteins via: "group efHand, *"](/images/0/00/Group_on1.png)](/index.php/File:Group_on1.png "Applied the group command to the proteins via: "group efHand, *"")

Applied the group command to the proteins via: "group efHand, *" 

  * [![The plus was clicked and expanded to show the hierarchy of objects.](/images/7/7f/Group_on2.png)](/index.php/File:Group_on2.png "The plus was clicked and expanded to show the hierarchy of objects.")

The plus was clicked and expanded to show the hierarchy of objects. 




## Contents

  * 1 Usage
  * 2 Examples
    * 2.1 Creating, opening and closing
    * 2.2 More advanced usage of groups and naming
  * 3 Notes
  * 4 See Also



## Usage
    
    
    group name, members, action
    

Actions: 

  * **add** \- Add member to group
  * **remove** \- Remove members from group
  * **open** \- Open the group in the panel so objects can be dragged in
  * **close** \- Close the group in the panel so nothing can be dragged in
  * **toggle** \- Switch between open or close based on current state
  * **auto** \- (Deprecated)
  * **ungroup** \- (Deprecated) use the ungroup command instead
  * **empty** \- Move members to top level but do not delete groups
  * **purge** \- Delete members but do not delete groups
  * **excise** \- Delete groups but do not delete members
  * **raise** \- (Incentive 3.1+ only) Move the specified group to the top level, relevant for groups within groups



## Examples

### Creating, opening and closing
    
    
    group efHand, 1cll 1ggz 1sra
    
    # allow addition and removal from the group
    # If a group is open, objects can be added to or removed from 
    # it by right-click+drag from the control panel
    group efHand, open
    # disallow addition/removal from the group
    group efHand, close
    

### More advanced usage of groups and naming
    
    
    # names with dots are treated special
    
    set group_auto_mode, 2
    
    # load the example protein
    
    load $TUT/1hpv.pdb, 1hpv.other
    
    # create the new entry called ".protein" in group 1hpv
    
    extract 1hpv.protein, 1hpv.other and polymer
    
    # create ".ligand in the 1hpv group
    
    extract 1hpv.ligand, 1hpv.other and organic
    
    # supports wildcards
    
    show sticks, *.ligand
    
    hide lines, *.protein
    
    show surface, *.protein within 6 of *.ligand
    
    show lines, byres *.protein within 4 of *.ligand
    
    set two_sided_lighting
    
    set transparency, 0.5
    
    set surface_color, white
    
    # Also, to lexicographically sort the names in the control panel:
    
    order *, yes
    

## Notes

Group objects can usually be used as arguments to commands. It can be processed as a group or as a selection, in which case all the atoms from all objects in the group will be used. 

  


## See Also

[ungroup](/index.php?title=Ungroup&action=edit&redlink=1 "Ungroup \(page does not exist\)") [order](/index.php/Order "Order") [select](/index.php/Select "Select")

Retrieved from "[https://pymolwiki.org/index.php?title=Group&oldid=13873](https://pymolwiki.org/index.php?title=Group&oldid=13873)"


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

## Hide

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**hide** conceals atom and bond representations for a certain selection or other graphical objects like distances. 

  * [![Some normal scene, notice the waters shown as spheres](/images/6/6b/Show1.png)](/index.php/File:Show1.png "Some normal scene, notice the waters shown as spheres")

Some normal scene, notice the waters shown as spheres 

  * [![hide spheres issues and the spheres are now hidden.](/images/e/ea/Show2.png)](/index.php/File:Show2.png "hide spheres issues and the spheres are now hidden.")

_hide spheres_ issues and the spheres are now hidden. 




The available representations are: 

  * lines
  * spheres
  * mesh
  * ribbon
  * cartoon
  * sticks
  * dots
  * surface
  * labels
  * nonbonded
  * nb_spheres



## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    hide representation [,object]
    hide representation [,(selection)]
    hide (selection)
    

### PYMOL API
    
    
    cmd.hide( string representation="", string selection="")
    

### EXAMPLES
    
    
    # hides all lines
    hide lines,all
    
    # hides all ribbons
    hide ribbon
    
    # hides all distances
    hide dashes
    
    # hides sticks in protA and all residues that aren't in the range of 40-65
    hide sticks, protA and not i. 40-65
    
    # hide hydrogen atoms
    hide (hydro)  # or hide (h.)
    

### SEE ALSO

[Show](/index.php/Show "Show"), [Enable](/index.php/Enable "Enable"), [Disable](/index.php/Disable "Disable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates")

Retrieved from "[https://pymolwiki.org/index.php?title=Hide&oldid=11564](https://pymolwiki.org/index.php?title=Hide&oldid=11564)"


---

## Id Atom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**id_atom** returns the original source id of a single atom, or raises and exception if the atom does not exist or if the selection corresponds to multiple atoms. 

### PYMOL API
    
    
    list = cmd.id_atom(string selection)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Id_Atom&oldid=7507](https://pymolwiki.org/index.php?title=Id_Atom&oldid=7507)"


---

## Identify

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**identify** returns a list of atom IDs corresponding to the ID code of atoms in the selection. 

## PYMOL API
    
    
    list = cmd.identify(string selection="(all)",int mode=0)
    

### NOTES

  * mode 0: only return a list of identifiers (default)
  * mode 1: return a list of tuples of the object name and the identifier



## See Also

[Selection](/index.php/Selection "Selection"), [Property_Selectors](/index.php/Property_Selectors "Property Selectors")

Retrieved from "[https://pymolwiki.org/index.php?title=Identify&oldid=12161](https://pymolwiki.org/index.php?title=Identify&oldid=12161)"


---

## Index

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**index** returns a list of tuples corresponding to the object name and index of the atoms in the selection. 

### PYMOL API
    
    
    list = cmd.index(string selection="(all)")
    

### NOTE

Atom indices are fragile and will change as atoms are added or deleted. Whenever possible, use integral atom identifiers instead of indices. 

Retrieved from "[https://pymolwiki.org/index.php?title=Index&oldid=7509](https://pymolwiki.org/index.php?title=Index&oldid=7509)"


---

## Indicate

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**indicate** shows a visual representation of an atom selection. 

The image on the left has nothing indicated. The image on the right has the helices "indicated," which PyMOL represents as small red dots. 

  * [![Ind0.png](/images/e/ed/Ind0.png)](/index.php/File:Ind0.png)

  * [![Ind1.png](/images/d/db/Ind1.png)](/index.php/File:Ind1.png)




### USAGE
    
    
    indicate (selection)
    

### PYMOL API
    
    
    cmd.count(string selection)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Indicate&oldid=7511](https://pymolwiki.org/index.php?title=Indicate&oldid=7511)"


---

## Intra fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**intra_fit** fits all states of an object to an atom selection in the specified state. It returns the rms values to python as an array. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
    * 3.1 Simple Selection
    * 3.2 Fitting NMR Ensembles
  * 4 PYTHON EXAMPLE
  * 5 USER EXAMPLES
  * 6 USER COMMENTS
  * 7 SEE ALSO



### USAGE
    
    
    intra_fit (selection),state
    

### PYMOL API
    
    
    cmd.intra_fit( string selection, int state )
    

### EXAMPLES

#### Simple Selection
    
    
    intra_fit ( name ca )
    

#### Fitting NMR Ensembles

Warren provided a great example on the PyMOL list. If the NMR ensemble has all the structures loaded as multiple states (which is the default behavoir (see [split_states](/index.php/Split_states "Split states")) for multimodel PDBs) then we can simply call: 
    
    
    print cmd.intra_fit(selection, state)
    

which will fit all of the other states to the indicated states and return the RMS values as a list. 

A simple example is: 
    
    
    fetch 1i8e, async=0
    print cmd.intra_fit("1i8e////CA", 1)
    
    # results:
    [-1.0, 1.1934459209442139, 1.2950557470321655, 0.71329855918884277,
    0.76704370975494385, 0.78973227739334106, 0.99323123693466187,
    1.0165935754776001, 0.6535714864730835, 0.95591926574707031,
    1.1299723386764526, 0.28637325763702393, 0.69836461544036865,
    0.40816938877105713, 1.1637680530548096]
    

Note that a RMS value of -1.0 is returned for the target state. 

### PYTHON EXAMPLE
    
    
    from pymol import cmd
    rms = cmd.intra_fit("(name ca)",1)
    

### USER EXAMPLES

### USER COMMENTS

See [Rms](/index.php/Rms "Rms") for selection caveats for this group of commands. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit"), [Extra_fit](/index.php/Extra_fit "Extra fit"), [Align](/index.php/Align "Align")

Retrieved from "[https://pymolwiki.org/index.php?title=Intra_fit&oldid=12573](https://pymolwiki.org/index.php?title=Intra_fit&oldid=12573)"


---

## Intra rms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**intra_rms** calculates rms fit values for all states of an object over an atom selection relative to the indicated state. Coordinates are left unchanged. The rms values are returned as a python array. 

## Contents

  * 1 PYMOL API
  * 2 PYTHON EXAMPLE
  * 3 USER COMMENTS
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.intra_rms( string selection, int state)
    

### PYTHON EXAMPLE
    
    
    from pymol import cmd
    rms = cmd.intra_rms("(name ca)",1)
    

### USER COMMENTS

See [Rms](/index.php/Rms "Rms") for selection caveats for this group of commands. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

Retrieved from "[https://pymolwiki.org/index.php?title=Intra_rms&oldid=11476](https://pymolwiki.org/index.php?title=Intra_rms&oldid=11476)"


---

## Intra rms cur

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**intra_rms_cur** calculates rms values for all states of an object over an atom selection relative to the indicated state without performing any fitting. The rms values are returned as a python array. 

## Contents

  * 1 PYMOL API
  * 2 PYTHON EXAMPLE
  * 3 USER EXAMPLES/COMMENTS
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.intra_rms_cur( string selection, int state)
    

### PYTHON EXAMPLE
    
    
    from pymol import cmd
    rms = cmd.intra_rms_cur("(name ca)",1)
    

### USER EXAMPLES/COMMENTS

See [Rms](/index.php/Rms "Rms") for selection caveats for this group of commands. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

Retrieved from "[https://pymolwiki.org/index.php?title=Intra_rms_cur&oldid=11477](https://pymolwiki.org/index.php?title=Intra_rms_cur&oldid=11477)"


---

## Invert

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**invert** inverts the stereo-chemistry of atom (pk1), holding attached atoms (pk2) and (pk3) immobile. 

### USAGE
    
    
    invert 
    

### PYMOL API
    
    
    cmd.invert( )
    

### NOTE

The invert function is usually bound to CTRL-E in Editing Mode. 

Retrieved from "[https://pymolwiki.org/index.php?title=Invert&oldid=7538](https://pymolwiki.org/index.php?title=Invert&oldid=7538)"


---

## Isodot

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**isodot** creates a dot isosurface object from a map object. 

### USAGE
    
    
    isodot name = map, level [,(selection) [,buffer [, state ] ] ] 
    map = the name of the map object to use.
    level = the contour level.
    selection = an atom selection about which to display the mesh with an additional "buffer" (if provided).
    

### NOTES

If the dot isosurface object already exists, then the new dots will be appended onto the object as a new state. 

### SEE ALSO

[load](/index.php/Load "Load"), [isomesh](/index.php/Isomesh "Isomesh"), [Dynamic_mesh](/index.php/Dynamic_mesh "Dynamic mesh")

Retrieved from "[https://pymolwiki.org/index.php?title=Isodot&oldid=10859](https://pymolwiki.org/index.php?title=Isodot&oldid=10859)"


---

## Isolevel

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**isolevel** changes the contour level of a isodot, isosurface, or isomesh object. 

### USAGE
    
    
    isolevel name, level, state
    

Retrieved from "[https://pymolwiki.org/index.php?title=Isolevel&oldid=8564](https://pymolwiki.org/index.php?title=Isolevel&oldid=8564)"


---

## Isomesh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**isomesh** creates a mesh isosurface object from a map object. 

## Contents

  * 1 Usage
  * 2 Example
  * 3 Details
    * 3.1 Selection Argument
    * 3.2 State Arguments
    * 3.3 MAP AROUND THE CENTER
    * 3.4 MAP OUTSIDE THE CALCULATED AREA
    * 3.5 MAP LEVELS
  * 4 See Also



## Usage
    
    
    isomesh name, map, level [,(selection) [,buffer [,state [,carve ]]]] 
    

  * name = str: the name for the new mesh isosurface object.
  * map = str: the name of the map object to use for computing the mesh.
  * level = float: the contour level (in sigma units) {default: 1.0}
  * selection = str: an atom selection about which to display the mesh with an additional "buffer" (if provided).
  * buffer = float: (see selection)
  * state = int: the state into which the object should be loaded (default=1) (set state=0 to append new mesh as a new state) {default: 1}
  * carve = float: a radius about each atom in the selection for which to include density. If "carve" is not provided, then the whole brick is displayed.



## Example
    
    
    fetch 6sps
    fetch 6sps, type=2fofc
    
    # mesh for entire map object
    isomesh mesh_all, 6sps_2fofc
    
    # mesh within bounding box of ligand, enlarged by 2 Angstrom
    isomesh mesh_ligand, 6sps_2fofc, selection=(resn LR5), buffer=2
    
    # mesh only within 2 Angstrom radius of any ligand atom
    isomesh mesh_ligand_carved, 6sps_2fofc, selection=(resn LR5), carve=2
    
    set_view (\
         0.001600198,   -0.993020296,    0.117938228,\
        -0.999629617,    0.001603606,    0.027057989,\
        -0.027055763,   -0.117936097,   -0.992654920,\
         0.000000000,    0.000000000,  -55.829845428,\
        12.989342690,    2.425159931,   19.217729568,\
        51.427371979,   60.232315063,  -19.999994278 )
    

mesh_all  | mesh_ligand  | mesh_ligand_carved   
---|---|---  
[![6sps-mesh-all.png](/images/1/19/6sps-mesh-all.png)](/index.php/File:6sps-mesh-all.png) | [![6sps-mesh-ligand.png](/images/e/e4/6sps-mesh-ligand.png)](/index.php/File:6sps-mesh-ligand.png) | [![6sps-mesh-ligand-carved.png](/images/d/df/6sps-mesh-ligand-carved.png)](/index.php/File:6sps-mesh-ligand-carved.png)  
  
## Details

#### Selection Argument

The arguments `selection`, `buffer` and `carve` can limit the mesh display to a selected area, and/or extend the area by symmetry operators if the selection is located outside the map bounding box itself. 

[![Isomesh-buffer-carve.png](/images/0/04/Isomesh-buffer-carve.png)](/index.php/File:Isomesh-buffer-carve.png)

#### State Arguments

If the mesh object already exists, then the new mesh will be appended onto the object as a new state (unless you indicate a state). 

  * state > 0: specific state
  * state = 0: all states
  * state = -1: current state


  * source_state > 0: specific state
  * source_state = 0: include all states starting with 0
  * source_state = -1: current state
  * source_state = -2: last state in map



#### MAP AROUND THE CENTER

You can create mesh around the center of the view by specifying "center" as the selection argument. 
    
    
    isomesh normal, fake_map, 1.0, center
    

#### MAP OUTSIDE THE CALCULATED AREA

When [map_auto_expand_sym](/index.php/Map_auto_expand_sym "Map auto expand sym") is ON, you can create mesh beyond the precalculated volume. In this case, symmetry information (lattice constants, space group) of the model specified in the selection argument if available, or (new in 1.7) from the map object. 

#### MAP LEVELS

  * Generally speaking there is some ambiguity with visualization tools as to how map data is to treated: Some map file formats are normalized by convention (in the file data itself) and others do not. Some visualization tools automatically normalize maps upon reading, others do not. PyMOL's default behavior is dependent upon map file type: CCP4 and O/BRIX/DSN6 maps are automatically normalized upon reading (disable via **normalize_*** settings), other maps types are not. PyMOL's normalization is a straight statistical average of all map points -- this may or may not be what you want. If migrating to PyMOL from another tool, then it is definitely worth comparing how the maps are being represented by creating an equivalent figure in both, making sure that they match, and if they do not, then figuring out why not. _From the PyMOL list. Author: Warren DeLano._



## See Also

  * [isodot](/index.php/Isodot "Isodot")
  * [load](/index.php/Load "Load")
  * [fetch](/index.php/Fetch "Fetch")
  * [dynamic_mesh](/index.php/Dynamic_mesh "Dynamic mesh")



Retrieved from "[https://pymolwiki.org/index.php?title=Isomesh&oldid=13332](https://pymolwiki.org/index.php?title=Isomesh&oldid=13332)"


---

## Isosurface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**isosurface** creates a new surface object from a map object. 

## Contents

  * 1 Usage
  * 2 Algorithm
  * 3 Examples
  * 4 Notes
  * 5 See Also



## Usage
    
    
     isosurface name, map, [,level [,(selection) [,buffer [,state [,carve [,source_state [,side [,mode ]]]]]]]]
    

  * name = the name for the new mesh isosurface object.
  * map = the name of the map object to use for computing the mesh.
  * level = the contour level. (default=1.0)
  * selection = an atom selection about which to display the mesh with an additional "buffer" (if provided). (default="")
  * state = the state into which the object should be loaded (default=1) (set state=-2 to append new surface as a new state)
  * carve = a radius about each atom in the selection for which to include density. If "carve= not provided, then the whole brick is displayed. (default=None)
  * source_state = the state of the map from which the object should be loaded. (default=0)
  * Front or back face. Triangle-winding/normal direction. (default=1)
  * mode = surface geometry (0: dots; 1: lines; 2: triangle triangle-normals; 3: triangle gradient-normals) (default=3)



## Algorithm

PyMOL offers three different algorithms for isosurface generation. Each of these can be activated by the `isosurface_algorithm` setting 

  * 0: Marching Cubes via VTKm (default) (requires VTKm)
  * 1: Marching Cubes basic (fallback if VTKm not installed)
  * 2: Marching tetrahedra (legacy)



## Examples
    
    
    fetch 1oky, type=2fofc, async=0
    isosurface 1okySurf, 1oky_2fofc
    

With carving at 2 Angstrom around the molecular model: 
    
    
    fetch 1oky, async=0
    fetch 1oky, type=2fofc, async=0
    isosurface 1okySurf, 1oky_2fofc, 1.0, (1oky), carve=2.0
    

## Notes

If there exists a non-map object with the same name, then the new surface will overwrite that object. Surface objects can be appended onto existing surface objects using the aforementioned state argument. 

## See Also

  * [volume](/index.php/Volume "Volume")
  * [isodot](/index.php/Isodot "Isodot")
  * [isomesh](/index.php/Isomesh "Isomesh")
  * [load](/index.php/Load "Load")
  * [Dynamic_mesh](/index.php/Dynamic_mesh "Dynamic mesh")



Retrieved from "[https://pymolwiki.org/index.php?title=Isosurface&oldid=13432](https://pymolwiki.org/index.php?title=Isosurface&oldid=13432)"


---

## Iterate

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **iterate** command executes a Python expression for all atoms in a selection. The local namespace exposes all atomic identifiers and properties as read-only Python variables. The global namespace by default is the **pymol** module namespace, and thus exposes objects like **cmd** and **[stored](/index.php/Pymol.stored "Pymol.stored")**. The order of iteration is that of the internal atom ordering. 

The **alter** command is equivalent to **iterate** , but provides read-write access to the exposed variables. This can for example be used to rename a chain or to assign new b-factors. If changing identifiers which would affect atom ordering, calling [sort](/index.php/Sort "Sort") is necessary to reorder based on the new identifiers. 

The **iterate_state** command is similar to **iterate** , but iterates over the coordinates in the given state and selection, and exposes **x, y, z** in addition to the atomic identifiers and properties. 

The **alter_state** command is equivalent to **iterate_state** , but allows modification of **x, y, z** (and atom-state level settings in Incentive PyMOL 1.7+). 

## Contents

  * 1 Exposed Variables
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
    * 3.1 Example: Print b-factors
    * 3.2 Example: Get b-factor list
    * 3.3 Example: Get coordinates
    * 3.4 Example: Modify coordinates
    * 3.5 Example: Rename chain
    * 3.6 Example: Sequence numbering offset
    * 3.7 Example: Net Charge
    * 3.8 Example: Color transfer
  * 4 PYMOL API
    * 4.1 Arguments
    * 4.2 "space" argument
  * 5 See Also



## Exposed Variables

All commands in the **iterate** -family expose the following variables: 

  * **model** (str): the object name (appearing in the selection panel on the right) (cannot be altered)
  * **name** (str): the atom name
  * **resn** (str): the residue name
  * **oneletter** (str): **resn** translated to one letter code, like it's displayed in the sequence viewer (cannot be altered, alter **resn** instead) (PyMOL 1.8.7+)
  * **resi** (str): the residue identifier (residue number) as a string, including optional insertion code
  * **resv** (int): the residue identifier (residue number) as an integer, excluding insertion code
  * **chain** (str): the chain name
  * **alt** (str): alternate location identifier
  * **elem** (str): the chemical element
  * **q** (float): occupancy
  * **b** (float): the B Factor
  * **segi** (str): segment identifier (columns 73-76 in PDB file)
  * **type** (str: ATOM,HETATM): the atom type (PDB convention for canonical polymer residues)
  * **formal_charge** (int): the formal charge of the atom
  * **partial_charge** (float): the partial charge of the atom
  * **numeric_type**
  * **text_type** (str): automatic mol2 atom typing (Incentive PyMOL 1.4+)
  * **stereo** (str): automatic stereochemical R/S label (Incentive PyMOL 1.4+)
  * **ID** (int): PDB atom id (not guaranteed to be unique)
  * **rank** (int): atom index from original file import (not guaranteed to be unique)
  * **index** (int): internal atom index (unique per object, sensitive to [sorting](/index.php/Sort "Sort") and [removing](/index.php/Remove "Remove") of atoms, cannot be altered)
  * **vdw** (float): Van der Waals radius
  * **ss** (str): secondary structure
  * **color** (int): color index
  * **reps** (int): numeric mask of shown representations (PyMOL 1.7.4+)
  * **protons** (int): (PyMOL 1.7.4+)
  * **p** (object): property object to access user-defined properties like **p.abc** (Incentive PyMOL 1.7+)
  * **s** (object): settings object to access settings, e.g. **s.sphere_scale** (Incentive PyMOL 1.7+, Open-Source PyMOL 1.8.1+)
  * **label** (str): see [label](/index.php/Label "Label")
  * **geom** (int): 1=single, 2=linear, 3=planar, 4=tetrahedral (see [set_geometry](/index.php?title=Set_geometry&action=edit&redlink=1 "Set geometry \(page does not exist\)"))
  * **valence** (int): expected number of bonds (?)
  * **flags** (int): bitmask, see [flag](/index.php/Flag "Flag")
  * **cartoon** (int): see [cartoon](/index.php/Cartoon "Cartoon")



The **iterate_state** and **alter_state** commands in addition expose: 

  * **x** (float): x-coordinate
  * **y** (float): y-coordinate
  * **z** (float): z-coordinate
  * **state** (int): (PyMOL 1.7.0+)



## Usage
    
    
    iterate (selection), expression
    
    
    
    iterate_state state, (selection), expression
    
    
    
    alter (selection), expression
    
    
    
    alter_state state, (selection), expression
    

### Arguments

  * **state** = int: object state, -1 for current state, 0 for all states
  * **selection** = str: atom selection
  * **expression** = str: expression in Python language



## Examples

### Example: Print b-factors

The following prints the atom names (**name** variable) and b-factors (**b** variable) for residue #1. 
    
    
    iterate (resi 1), print(name + " %.2f" % b)
    

### Example: Get b-factor list

The following example fills a list, **[stored](/index.php/Pymol.stored "Pymol.stored").bfactors** with the f-factors (**b** variable) from residue #1. 
    
    
    stored.bfactors = []
    iterate (resi 1), stored.bfactors.append(b)
    print(stored.bfactors)
    

### Example: Get coordinates
    
    
    stored.coords = []
    iterate_state 1, (all), stored.coords.append([x,y,z])
    

### Example: Modify coordinates

This example shifts the selection by 10 Angstrom along the model x-axis. 
    
    
    alter_state 1, (all), x += 10.0
    rebuild
    

### Example: Rename chain

This renames chain **A** to **C**. Note the difference between the selection syntax (first argument) and the Python syntax with the quoted string (second argument). The [sort](/index.php/Sort "Sort") call is necessary if let's say there is also a chain B, and you expect chain B to appear before chain C (formerly A) in the sequence viewer. 
    
    
    alter chain A, chain="C"
    sort
    

### Example: Sequence numbering offset

Assume the residue numbering in a PDB file starts at 100, then the following updates the residue numbers to start at 1. 
    
    
    alter (chain A), resv -= 99
    sort
    

### Example: Net Charge

The following example calculates the net charge of an object. 
    
    
    stored.net_charge = 0
    iterate (all), stored.net_charge += partial_charge
    print('Net charge: ' + str(stored.net_charge))
    

### Example: Color transfer

Copy (transfer) the color from one object to another 
    
    
    stored.colors = {}
    iterate obj1, stored.colors[chain,resi,name] = color
    alter obj2, color = stored.colors.get((chain,resi,name), color)
    recolor
    

## PYMOL API
    
    
    cmd.iterate(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.iterate_state(int state, str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter_state(int state, str selection, str expression, int quiet=1, dict space=None)
    

### Arguments

  * **state** = int: state-index if positive number or any of these: 
    * **state** = 0: all states
    * **state** = -1: current state
  * **selection** = str: atom selection
  * **expression** = str: expression in valid [python syntax](http://en.wikipedia.org/wiki/Python_syntax_and_semantics)
  * **space** = dict: namespace dictionary {default: pymol namespace}
  * **atomic** = 0/1: provide atomic properties as variables if 1, or only x/y/z if 0 (in older PyMOL versions, atomic=0 gives some speed improvement) {default: 1}



### "space" argument

The **space** argument can be used to pass local objects into the expression namespace, instead of evaluating the expression in the global **pymol** module namespace. 

The b-factor list example from above but without the global [pymol.stored](/index.php/Pymol.stored "Pymol.stored") variable: 
    
    
    myspace = {'bfactors': []}
    cmd.iterate('(all)', 'bfactors.append(b)', space=myspace)
    

User defined functions can also be included in the namespace: 
    
    
    def myfunc(resi,resn,name):
        print('%s`%s/%s' % (resn ,resi, name))
    
    myspace = {'myfunc': myfunc}
    cmd.iterate('(all)', 'myfunc(resi,resn,name)', space=myspace)
    

## See Also

  * [sort](/index.php/Sort "Sort")
  * [get_model](/index.php/Get_Model "Get Model")



Retrieved from "[https://pymolwiki.org/index.php?title=Iterate&oldid=12709](https://pymolwiki.org/index.php?title=Iterate&oldid=12709)"


---

## Iterate State

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **iterate** command executes a Python expression for all atoms in a selection. The local namespace exposes all atomic identifiers and properties as read-only Python variables. The global namespace by default is the **pymol** module namespace, and thus exposes objects like **cmd** and **[stored](/index.php/Pymol.stored "Pymol.stored")**. The order of iteration is that of the internal atom ordering. 

The **alter** command is equivalent to **iterate** , but provides read-write access to the exposed variables. This can for example be used to rename a chain or to assign new b-factors. If changing identifiers which would affect atom ordering, calling [sort](/index.php/Sort "Sort") is necessary to reorder based on the new identifiers. 

The **iterate_state** command is similar to **iterate** , but iterates over the coordinates in the given state and selection, and exposes **x, y, z** in addition to the atomic identifiers and properties. 

The **alter_state** command is equivalent to **iterate_state** , but allows modification of **x, y, z** (and atom-state level settings in Incentive PyMOL 1.7+). 

## Contents

  * 1 Exposed Variables
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
    * 3.1 Example: Print b-factors
    * 3.2 Example: Get b-factor list
    * 3.3 Example: Get coordinates
    * 3.4 Example: Modify coordinates
    * 3.5 Example: Rename chain
    * 3.6 Example: Sequence numbering offset
    * 3.7 Example: Net Charge
    * 3.8 Example: Color transfer
  * 4 PYMOL API
    * 4.1 Arguments
    * 4.2 "space" argument
  * 5 See Also



## Exposed Variables

All commands in the **iterate** -family expose the following variables: 

  * **model** (str): the object name (appearing in the selection panel on the right) (cannot be altered)
  * **name** (str): the atom name
  * **resn** (str): the residue name
  * **oneletter** (str): **resn** translated to one letter code, like it's displayed in the sequence viewer (cannot be altered, alter **resn** instead) (PyMOL 1.8.7+)
  * **resi** (str): the residue identifier (residue number) as a string, including optional insertion code
  * **resv** (int): the residue identifier (residue number) as an integer, excluding insertion code
  * **chain** (str): the chain name
  * **alt** (str): alternate location identifier
  * **elem** (str): the chemical element
  * **q** (float): occupancy
  * **b** (float): the B Factor
  * **segi** (str): segment identifier (columns 73-76 in PDB file)
  * **type** (str: ATOM,HETATM): the atom type (PDB convention for canonical polymer residues)
  * **formal_charge** (int): the formal charge of the atom
  * **partial_charge** (float): the partial charge of the atom
  * **numeric_type**
  * **text_type** (str): automatic mol2 atom typing (Incentive PyMOL 1.4+)
  * **stereo** (str): automatic stereochemical R/S label (Incentive PyMOL 1.4+)
  * **ID** (int): PDB atom id (not guaranteed to be unique)
  * **rank** (int): atom index from original file import (not guaranteed to be unique)
  * **index** (int): internal atom index (unique per object, sensitive to [sorting](/index.php/Sort "Sort") and [removing](/index.php/Remove "Remove") of atoms, cannot be altered)
  * **vdw** (float): Van der Waals radius
  * **ss** (str): secondary structure
  * **color** (int): color index
  * **reps** (int): numeric mask of shown representations (PyMOL 1.7.4+)
  * **protons** (int): (PyMOL 1.7.4+)
  * **p** (object): property object to access user-defined properties like **p.abc** (Incentive PyMOL 1.7+)
  * **s** (object): settings object to access settings, e.g. **s.sphere_scale** (Incentive PyMOL 1.7+, Open-Source PyMOL 1.8.1+)
  * **label** (str): see [label](/index.php/Label "Label")
  * **geom** (int): 1=single, 2=linear, 3=planar, 4=tetrahedral (see [set_geometry](/index.php?title=Set_geometry&action=edit&redlink=1 "Set geometry \(page does not exist\)"))
  * **valence** (int): expected number of bonds (?)
  * **flags** (int): bitmask, see [flag](/index.php/Flag "Flag")
  * **cartoon** (int): see [cartoon](/index.php/Cartoon "Cartoon")



The **iterate_state** and **alter_state** commands in addition expose: 

  * **x** (float): x-coordinate
  * **y** (float): y-coordinate
  * **z** (float): z-coordinate
  * **state** (int): (PyMOL 1.7.0+)



## Usage
    
    
    iterate (selection), expression
    
    
    
    iterate_state state, (selection), expression
    
    
    
    alter (selection), expression
    
    
    
    alter_state state, (selection), expression
    

### Arguments

  * **state** = int: object state, -1 for current state, 0 for all states
  * **selection** = str: atom selection
  * **expression** = str: expression in Python language



## Examples

### Example: Print b-factors

The following prints the atom names (**name** variable) and b-factors (**b** variable) for residue #1. 
    
    
    iterate (resi 1), print(name + " %.2f" % b)
    

### Example: Get b-factor list

The following example fills a list, **[stored](/index.php/Pymol.stored "Pymol.stored").bfactors** with the f-factors (**b** variable) from residue #1. 
    
    
    stored.bfactors = []
    iterate (resi 1), stored.bfactors.append(b)
    print(stored.bfactors)
    

### Example: Get coordinates
    
    
    stored.coords = []
    iterate_state 1, (all), stored.coords.append([x,y,z])
    

### Example: Modify coordinates

This example shifts the selection by 10 Angstrom along the model x-axis. 
    
    
    alter_state 1, (all), x += 10.0
    rebuild
    

### Example: Rename chain

This renames chain **A** to **C**. Note the difference between the selection syntax (first argument) and the Python syntax with the quoted string (second argument). The [sort](/index.php/Sort "Sort") call is necessary if let's say there is also a chain B, and you expect chain B to appear before chain C (formerly A) in the sequence viewer. 
    
    
    alter chain A, chain="C"
    sort
    

### Example: Sequence numbering offset

Assume the residue numbering in a PDB file starts at 100, then the following updates the residue numbers to start at 1. 
    
    
    alter (chain A), resv -= 99
    sort
    

### Example: Net Charge

The following example calculates the net charge of an object. 
    
    
    stored.net_charge = 0
    iterate (all), stored.net_charge += partial_charge
    print('Net charge: ' + str(stored.net_charge))
    

### Example: Color transfer

Copy (transfer) the color from one object to another 
    
    
    stored.colors = {}
    iterate obj1, stored.colors[chain,resi,name] = color
    alter obj2, color = stored.colors.get((chain,resi,name), color)
    recolor
    

## PYMOL API
    
    
    cmd.iterate(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.iterate_state(int state, str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter(str selection, str expression, int quiet=1, dict space=None)
    
    cmd.alter_state(int state, str selection, str expression, int quiet=1, dict space=None)
    

### Arguments

  * **state** = int: state-index if positive number or any of these: 
    * **state** = 0: all states
    * **state** = -1: current state
  * **selection** = str: atom selection
  * **expression** = str: expression in valid [python syntax](http://en.wikipedia.org/wiki/Python_syntax_and_semantics)
  * **space** = dict: namespace dictionary {default: pymol namespace}
  * **atomic** = 0/1: provide atomic properties as variables if 1, or only x/y/z if 0 (in older PyMOL versions, atomic=0 gives some speed improvement) {default: 1}



### "space" argument

The **space** argument can be used to pass local objects into the expression namespace, instead of evaluating the expression in the global **pymol** module namespace. 

The b-factor list example from above but without the global [pymol.stored](/index.php/Pymol.stored "Pymol.stored") variable: 
    
    
    myspace = {'bfactors': []}
    cmd.iterate('(all)', 'bfactors.append(b)', space=myspace)
    

User defined functions can also be included in the namespace: 
    
    
    def myfunc(resi,resn,name):
        print('%s`%s/%s' % (resn ,resi, name))
    
    myspace = {'myfunc': myfunc}
    cmd.iterate('(all)', 'myfunc(resi,resn,name)', space=myspace)
    

## See Also

  * [sort](/index.php/Sort "Sort")
  * [get_model](/index.php/Get_Model "Get Model")



Retrieved from "[https://pymolwiki.org/index.php?title=Iterate&oldid=12709](https://pymolwiki.org/index.php?title=Iterate&oldid=12709)"


---

## Label

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![Pl.png](/images/b/b9/Pl.png)](/index.php/File:Pl.png)

The Label command controls how PyMOL draws text labels for PyMOL objects. 

## Contents

  * 1 Details
    * 1.1 Built-in Object Properties
  * 2 Syntax
  * 3 Settings
  * 4 Examples
  * 5 User Comments
    * 5.1 Labels Using ID Numbers
    * 5.2 Labels Using One Letter Abbreviations
    * 5.3 Labels and defer_builds_mode
  * 6 See Also



# Details

Labeling is important so there are many options for your fine tuning needs. You can change the [label size](/index.php/Label_size "Label size"), [label color](/index.php/Label_color "Label color"), positioning, [font](/index.php/Label_font_id "Label font id"), the [label outline color](/index.php/Label_outline_color "Label outline color") that masks the font and much, much more. 

You can have PyMOL label atoms by properties or arbitrary strings as you want; you can even use Unicode fonts for special symbols like,  α , β , ± , \mathrm {\AA}  {\displaystyle \alpha ,\beta ,\pm ,{\textrm {\mathrm {\AA} }}} ![{\\displaystyle \\alpha ,\\beta ,\\pm ,{\\textrm {\\mathrm {\\AA} }}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/07cd06e9fdbee0a625881a9a0ee07dadcde15a66), etc. 

The following gallery shows some examples of how extensible the Label command is. 

  * [![Simple label](/images/7/78/Label_pre.png)](/index.php/File:Label_pre.png "Simple label")

Simple label 

  * [![Example showing usage of Unicode fonts for special characters, see label_font_id.](/images/e/e0/New_fonts.jpeg)](/index.php/File:New_fonts.jpeg "Example showing usage of Unicode fonts for special characters, see label_font_id.")

Example showing usage of Unicode fonts for special characters, see [label_font_id](/index.php/Label_font_id "Label font id"). 

  * [![Another example with Unicode fonts](/images/2/2f/Font_ex.png)](/index.php/File:Font_ex.png "Another example with Unicode fonts")

Another example with Unicode fonts 

  * [![Example label](/images/c/c2/Label_ex.png)](/index.php/File:Label_ex.png "Example label")

Example label 

  * [![Label shadows turned off](/images/0/05/Ls0.png)](/index.php/File:Ls0.png "Label shadows turned off")

Label shadows turned off 

  * [![Label shadows turned on](/images/8/89/Ls2.png)](/index.php/File:Ls2.png "Label shadows turned on")

Label shadows turned on 




## Built-in Object Properties

Aside from arbitrary string labels, like "This is the catalytic residue" for an atom/residue you can also use the following built-in molecular properties: 

  * **name** , the atom name
  * **resn** , the residue name
  * **resi** , the residue number/identifier
  * **chain** , the chain name
  * **q** , charge
  * **b** , the occupancy/b-factor
  * **segi** , the segment identifier
  * **type** _(ATOM,HETATM)_ , the type of atom
  * **formal_charge** , the formal charge
  * **partial_charge** , the partial charge
  * **numeric_type** , the numeric type
  * **text_type** , the text type



You can use one of these properties as: 
    
    
    # simple example: label residue 22's atoms with their names
    label i. 22, name
    
    # Label residue #44's alpha carbon with it's residue name, residue number and B-factor.
    label n. CA and i. 44, "(%s, %s, %s)" % (resn, resi, b)
    

See the syntax and examples below for more info. 

# Syntax

To use the label command follow this syntax: 
    
    
    # labeling syntax
    label [ selection[, expression]]
    

where **selection** is some object/selection you want to label and **expression** is some string (or set of strings) which PyMOL is to use to label the given selection. 

We have plenty of examples. See the examples below. 

# Settings

Here are all the label settings and their general effect. For each label setting, see the respective web page for more details. 

**[label_angle_digits](/index.php/Label_angle_digits "Label angle digits")**

    

    sets the number of decimals in angle label.

**[label_distance_digits](/index.php/Label_distance_digits "Label distance digits")**

    

    sets the number of decimals in distance label.

**[label_shadow_mode](/index.php/Label_shadow_mode "Label shadow mode")**

    

    sets whether or not PyMOL will ray trace shadows for your label text. Eg: 
    
    
    set label_shadow_mode, 2
    

**[label_color](/index.php/Label_color "Label color")**

    

    sets the color of the label text. Note that you can have labels of different colors for different objects or selections. Some examples:
    
    
    # per-object:
    set label_color, color-name, object-name  #eg, set label-color, magenta, /protein
    
    # per-atom:
    set label_color, color-name, selection    #eg, set label-color, marine, /protein/A/A/23/CA
    
    # another example
    fragment arg
    label all, name
    set label_color, yellow, arg
    set label_color, red, elem c
    

**[label_font_id](/index.php/Label_font_id "Label font id")**

    

    sets the font to render your label. There are 12 different fonts from 5—16. Numbers 15 and 16 are special for unicode. Eg: 
    
    
    set label_font_id, 12
    

. See the [label_font_id](/index.php/Label_font_id "Label font id") page for explicit examples on how to use unicode characters in PyMOL labels.

**[label_size](/index.php/Label_size "Label size")**

    

    sets the size of the text. You can use positive numbers 2, 3, 4, etc for point sizes, or negative numbers for Angstrom-based sizes. Default is 14 points. Labels in Angstrom-size scale with the distance from the front plane, labels in point-size don't. Eg: 
    
    
    set label_size, -2  #results in a size of 2 Angstroms
    

**[label_digits](/index.php?title=Label_digits&action=edit&redlink=1 "Label digits \(page does not exist\)")**

    

    sets the number of decimals in label. It affects all digits only if label_distance_digits or label_dihedral_digits or label_angle_digits are set to -1.

**[label_outline_color](/index.php/Label_outline_color "Label outline color")**

    

    each label is outlined (so you can do white-on-white labels, for example). This options sets the color of the label outline. Eg. 
    
    
    set label_outline_color, orange
    

**[label_dihedral_digits](/index.php/Label_dihedral_digits "Label dihedral digits")**

    

    sets the number of decimals in dihedral label.

**[label_position](/index.php/Label_position "Label position")**

    

    sets any offset from the original X,Y,Z coordinates for the label. If you like to use the mouse, you can enter [edit_mode](/index.php?title=Edit_mode&action=edit&redlink=1 "Edit mode \(page does not exist\)") and **ctrl-left_click** to drag labels around; **ctrl-shift-left_click** will let you move the labels in the z-direction. **"Save labels"-workaround** If you want to save the position of your labels, the best way might be to create a new object and move the atoms in this object. Since the labels are positioned from the atom positions this is an indirect way of moving the labels and being able to save them.

# Examples
    
    
    #1.
    # make a very simple label on the 14th alpha carbon.
    label n. CA and i. 14, "This is carbon 14."
    
    #2.
    # make a fake scene label; use this to label entire scenes, not just atoms/bonds.
    pseudoatom foo
    label foo, "Once upon a time..."
    
    #3.
    # make a huge label
    set label_size, -5
    pseudoatom foo
    label foo, "This is large text"
    
    #4. Partial Charge
    label (chain A),chain
    label (n;ca),"%s-%s" % (resn,resi)
    label (resi 200),"%1.3f" % partial_charge
    
    
    #5. The gallery image above Label_ex.png was created with this code
    #   and finally, some labels were moved around in '''edit_mode'''.
    label (resi 200),"%1.3f" % b
    set label_font_id, 10
    set label_size, 10
    
    #6. This example shows how to label a selection with the 
    #   XYZ coordinates of the atoms 
    from pymol import stored
    stored.pos = []
    # select the carbon atoms in my hetero atoms to label
    select nn, het and e. C
    # get the XYZ coordinates and put them into stored.pos
    # insert at the front because pop() will read the array in reverse
    iterate_state 1, (nn), stored.pos.insert(0,(x,y,z))
    # label all N atoms.  You need the pop() function or else
    # PyMOL will complain b/c you didn't provide enough coords.
    label nn, ("%5.5s, %5.5s, %5.5s") %  stored.pos.pop()
    

# User Comments

## Labels Using ID Numbers

The following commnent, 
    
    
    label SELECTION, " %s" % ID
    

labels the SELECTION with atom ID numbers. 

You can make more complicated selections/lables such as 
    
    
    label SELECTION, " %s:%s %s" % (resi, resn, name)
    

which will give you something like "GLU:139 CG" 

## Labels Using One Letter Abbreviations

  * First, Add this to your $HOME/.pymolrc file:


    
    
    # start $HOME/.pymolrc modification
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}
    # end modification
    

  * Second, instead of:


    
    
    label n. ca, resn
    

use: 
    
    
    label n. ca, one_letter[resn]
    

or: ( to get something like D85) 
    
    
    label n. ca, "%s%s" %(one_letter[resn],resi)
    

## Labels and defer_builds_mode

If You have a weak video card, You might want to set 
    
    
    set defer_builds_mode, 5
    

It helps a lot but breaks labels rendering. You can use 
    
    
    set defer_builds_mode, 4
    

instead. 

# See Also

[Category:Labeling](/index.php/Category:Labeling "Category:Labeling")

All the settings posted above. 

Retrieved from "[https://pymolwiki.org/index.php?title=Label&oldid=12514](https://pymolwiki.org/index.php?title=Label&oldid=12514)"


---

## Talk:Label

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Fixes
  * 2 New Page Content
    * 2.1 New Page Overview
  * 3 Old Page
    * 3.1 DESCRIPTION
    * 3.2 USAGE
    * 3.3 SETTINGS
      * 3.3.1 FONT
        * 3.3.1.1 UTF8 Fonts
        * 3.3.1.2 Unicode Fonts
      * 3.3.2 SIZE
      * 3.3.3 COLOR
      * 3.3.4 EXPRESSION
      * 3.3.5 POSITION
    * 3.4 EXAMPLES
      * 3.4.1 Partial Charge
      * 3.4.2 Example 2
      * 3.4.3 More Advanced
    * 3.5 Users Comments
      * 3.5.1 Labels Using ID Numbers
      * 3.5.2 Labels Using One Letter Abbreviations



## Fixes

  * Updates needed



# New Page Content

## New Page Overview

This is the content for the new labels page. 

  


# Old Page

[![PyMol Labels](/images/7/78/Label_pre.png)](/index.php/File:Label_pre.png "PyMol Labels")

### DESCRIPTION

**label** allows one to configure the appearance of text labels for PyMOL objects. It labels one or more atoms properties over a selection using the python evaluator with a separate name space for each atom. The symbols defined in the name space are: 

  * **name** , the atom name
  * **resn** , the residue name
  * **resi** , the residue number/identifier
  * **chain** , the chain name
  * **q** ,
  * **b** , the occupancy/b-factor
  * **segi** , the segment identifier
  * **type** _(ATOM,HETATM)_ , the type of atom
  * **formal_charge** , the formal charge
  * **partial_charge** , the partial charge
  * **numeric_type** , the numeric type
  * **text_type** , the text type



All strings in the expression must be explicitly quoted. This operation typically takes several seconds per thousand atoms altered. To clear labels, simply omit the expression or set it to _._

[Label](/index.php/Label "Label") is great for labeling atoms, residues and objects. For a scene label, see [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). 

### USAGE
    
    
    label (selection),expression
    

  


### SETTINGS

#### FONT

There are 10 different scalable fonts. 
    
    
    set label_font_id, number
    

where number is 5 through 14. 

##### UTF8 Fonts

[![](/images/e/e0/New_fonts.jpeg)](/index.php/File:New_fonts.jpeg)

[](/index.php/File:New_fonts.jpeg "Enlarge")

New fonts in PyMol. Notice the alpha and beta characters.

Newer versions support UTF8 fonts; use **label_font_id** from above to 15 or 16. The good news about the UTF8 fonts is that they support the alpha and beta characters. (See image.) 

Here's some example code for the image at right: 
    
    
    # roman
    set label_font_id, 15
    set label_shadow_mode, 3
    label 5/CA, "\316\261-Helix"
    label 10/CA, "\316\262-Sheet"
    
    # italic
    set label_font_id, 16
    
    # make bigger
    set label_size, 50
    

##### Unicode Fonts

[![](/images/2/2f/Font_ex.png)](/index.php/File:Font_ex.png)

[](/index.php/File:Font_ex.png "Enlarge")

Notice the Angstrom and superscript 2 characters. You can add other characters as well.

PyMOL gives you the flexibility to use encoded unicode fonts. This allows us to insert various symbols, like the symbol used for Angstrom. Here are the steps to insert a character from the unicode character set. 

  * Find the code for your character at [Unicode Charts](http://www.unicode.org/charts). The Angstrom character,  \mathrm {\AA}  {\displaystyle {\textrm {\mathrm {\AA} }}} ![{\\displaystyle {\\textrm {\\mathrm {\\AA} }}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/fe5b45baabc495a7ef582330f075a45a1d6bb5fc) is u"\u00c5" and  ± {\displaystyle \pm } ![{\\displaystyle \\pm }](https://wikimedia.org/api/rest_v1/media/math/render/svg/869e366caf596564de4de06cb0ba124056d4064b) is u"\u00b1".
  * Label the selection. For simple strings, just type the string in double quote, -- "like this" -- and append to the end of that .encode('utf-8') -- "like this".encode('utf-8'). A working example is shown here,


    
    
    # label residue 30 with "4.1 Ang^2 +/- 0.65 Ang^2; see the image at right
    label i. 30, "4.1" + u"\u00c5\u00b2  \u00b1 0.65 \u00c5\u00b2 ".encode('utf-8')
    

#### SIZE

The font size can be adjusted 
    
    
    set label_size, number
    

where number is the point size (or -number for Angstroms) 

#### COLOR

Set a label's color by 
    
    
    set label_color, color
    

where color is a valid PyMol color. 

If the coloring of the labels is not _exactly_ the same as you'd expect (say black turns out grey, or red turns out pink), then try the following settings: 
    
    
    unset depth_cue
    unset ray_label_specular
    

#### EXPRESSION

To set what the label reads (see above) 
    
    
    label selection, expression
    

For example 
    
    
     label all, name
     label resi 10, b
    

#### POSITION

To position labels 
    
    
    edit_mode
    

then ctrl-middle-click-and-drag to position the label in space. (On Windows systems this appears to be shift-left-click-and-drag, presumably because those mice lack a true middle button.) 

ctrl-shift-left-click-and-drag alters a label's z-plane. (Windows only? This may use the middle button, rather than shift-left, under *NIX / 3-button mice systems.) 

### EXAMPLES

#### Partial Charge
    
    
    label (chain A),chain
    label (n;ca),"%s-%s" % (resn,resi)
    label (resi 200),"%1.3f" % partial_charge
    

#### Example 2

The following image was created with 
    
    
    label (resi 200),"%1.3f" % b
    set label_font_id, 10
    set label_size, 10
    

and finally, some labels were moved around in **edit_mode**. 

[![](/images/c/c2/Label_ex.png)](/index.php/File:Label_ex.png)

[](/index.php/File:Label_ex.png "Enlarge")

Labels.

  


#### More Advanced

This example shows how to label a selection with the XYZ coordinates of the atoms 
    
    
    from pymol import stored
    stored.pos = []
    
    # select the carbon atoms in my hetero atoms to label
    select nn, het and e. C
    
    # get the XYZ coordinates and put htem into stored.pos
    iterate_state 1, (nn), stored.pos.append((x,y,z))
    
    # label all N atoms.  You need the pop() function or else
    # PyMOL will complain b/c you didn't provide enough coords.
    label nn, ("%5.5s, %5.5s, %5.5s") %  stored.pos.pop()
    

### Users Comments

#### Labels Using ID Numbers

The following commnent, 
    
    
    label SELECTION, " %s" % ID 
    

labels the SELECTION with atom ID numbers. 

You can make more complicated selections/lables such as 
    
    
    label SELECTION, " %s:%s %s" % (resi, resn, name)
    

which will give you something like "GLU:139 CG" 

#### Labels Using One Letter Abbreviations

  * First, Add this to your $HOME/.pymolrc file:


    
    
    # start $HOME/.pymolrc modification
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}
    # end modification
    

  * . Second, instead of:


    
    
    label n. ca, resn
    

use: 
    
    
    label n. ca, one_letter[resn]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Talk:Label&oldid=12387](https://pymolwiki.org/index.php?title=Talk:Label&oldid=12387)"


---

## Lines

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Lines** is name of the basic representation for atoms and bonds in PyMOL. **Lines** is a very simple representation, where each atom bond is displayed as a single colored line, and each atom is displayed as the intersection of any two or more non-terminal bonds. 

## Contents

  * 1 Usage
    * 1.1 Examples
      * 1.1.1 Example: Displaying dashed lines between two atoms
    * 1.2 See Also



# Usage
    
    
    # show everything as lines
    show lines
    
    # only show residues 50-80 as lines
    show lines, i.50-80
    

## Examples

### Example: Displaying dashed lines between two atoms

The following commands will create a dashed line between two atoms. 
    
    
    # first, create two named selections
    select a, ///A/501/02
    select b, ///B/229/N
    # calculate & show the distance from selection a to selection b.
    distance d, a, b
    # hide just the distance labels; the 
    # dashed bars should still be shown
    hide labels, d
    

Technically, the object _d_ is a labelled distance, only the label is hidden. When ray-tracing the image, the dashes come out a bit fat. You can slim them with 
    
    
    set dash_gap, 0.5
    set dash_radius, 0.1
    

before the 'ray' command. 

[![](/images/d/d4/Lines_ex.png)](/index.php/File:Lines_ex.png)

[](/index.php/File:Lines_ex.png "Enlarge")

Lines Representation Example

## See Also

Please read about other representations in the **[Representation Category](/index.php/Category:Representations "Category:Representations")**.   
[Measure_Distance](/index.php/Measure_Distance "Measure Distance")   
[Distance](/index.php/Distance "Distance")   


Retrieved from "[https://pymolwiki.org/index.php?title=Lines&oldid=7527](https://pymolwiki.org/index.php?title=Lines&oldid=7527)"


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

## Load Model

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**load_model** reads a ChemPy model into an object. If a trajectory (many models with the same name and different values for state) is being loaded, It can be used with the option discrete=1 to allow changes in the b-factors and Van der Waals radius between snapshots. 

### PYMOL API
    
    
    cmd.load_model(model, object [,state [,finish [,discrete ]]])
    

### See Also

[Get_Model](/index.php/Get_Model "Get Model")

Retrieved from "[https://pymolwiki.org/index.php?title=Load_Model&oldid=7529](https://pymolwiki.org/index.php?title=Load_Model&oldid=7529)"


---

## Load traj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**load_traj** loads a trajectory as "states" into an already loaded molecular object. 

Since version 1.0, PyMOL uses the [Molfile Plugin](http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/) backend, which supports a variety of trajectory file formats. Older versions only supported the ascii AMBER format (".trj" file extension). 

Loading a large trajectory may take up a lot of RAM, unless the [defer_builds_mode](/index.php/Defer_builds_mode "Defer builds mode") is set to 3. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 Notes
  * 5 See Also



## Usage
    
    
    load_traj filename [,object [,state [,format [,interval [,average ]
                       [,start [,stop [,max [,selection [,image [,shift 
                       [, plugin ]
                       ]]]]]]]]]
    

## Arguments

  * **filename** = str: trajectory file path
  * **object** = str: name of the molecular object where the trajectory should be appended as states {default: guess from filename}
  * **state** = int: object state where to start appending states. To discard the currently loaded coordinates, use _state=1_. To append new states, use _state=0_ {default: 0}
  * **format** = str: specify file type instead of guessing from file extension (only affects AMBER .trj format, use "plugin" argument for Molfile Plugin types) {default: }



## Examples
    
    
    # topology from PDB file, trajectory from DCD file
    load      sampletrajectory.pdb
    load_traj sampletrajectory.dcd
    
    # gromacs trajectory, using "mytraj" as object name
    load      sampletrajectory.gro, mytraj
    load_traj sampletrajectory.xtc, mytraj
    
    # desmond trajectory
    load      sample-out.cms, mytraj
    load_traj sample_trj/clickme.dtr, mytraj
    
    # playing through states, memory optimized (but eventually slower)
    set defer_builds_mode, 3
    mplay
    

## Notes

  * PyMOL does not know how to wrap the truncated octahedron used by Amber You will need to use the [cpptraj](http://ambermd.org/tutorials/analysis/tutorial0/index.htm) program first to do this.
  * The average option is not a running average. To perform this type of average, use the [smooth](/index.php/Smooth "Smooth") command after loading the trajectory file.
  * For quickly viewing Trajectories as a movie, use the [mset](/index.php/Mset "Mset") command to map each state to a movie frame.



useful notes from the email list:   
<http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg04272.html>   
<http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10266.html>

in one line convert dcd and psf to pdb : 

catdcd -o all.pdb -otype pdb -s autopsf.psf -stype psf out.dcd 

## See Also

[Load](/index.php/Load "Load"), [defer_builds_mode](/index.php/Defer_builds_mode "Defer builds mode")

Retrieved from "[https://pymolwiki.org/index.php?title=Load_traj&oldid=12888](https://pymolwiki.org/index.php?title=Load_traj&oldid=12888)"


---

## Ls

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**ls** lists contents of the current working directory. PyMOL comes with a few commands that allows the user to traverse a directory structure and list files. This is a powerful tool if you have many directories full of files --- you don't have to quit PyMOL, you can just **cd** to the directory and **ls** to find the contents. 

### USAGE
    
    
    ls [pattern]
    dir [pattern]
    

### EXAMPLES
    
    
    # list all files in the current directory
    ls
    
    # list all files ending in, pml, in the current directory. 
    ls *.pml
    

### SEE ALSO

[Cd](/index.php/Cd "Cd"), [Pwd](/index.php/Pwd "Pwd"), [System](/index.php/System "System")

Retrieved from "[https://pymolwiki.org/index.php?title=Ls&oldid=8252](https://pymolwiki.org/index.php?title=Ls&oldid=8252)"


---

## Map double

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**map_double** resamples a map at twice the current resolution. The amount of memory required to store the map will increase eight-fold. 

  * [![Std. map mesh spacing](/images/a/a8/Map_normal.png)](/index.php/File:Map_normal.png "Std. map mesh spacing")

Std. map mesh spacing 

  * [![Map doubled](/images/c/c2/Map_double.png)](/index.php/File:Map_double.png "Map doubled")

Map doubled 

  * [![Map double, doubled](/images/8/8f/Map_double2.png)](/index.php/File:Map_double2.png "Map double, doubled")

Map double, doubled 

  * [![Map double, double, doubled](/images/6/63/Map_double3.png)](/index.php/File:Map_double3.png "Map double, double, doubled")

Map double, double, doubled 




# Usage
    
    
    map_double map_name [, state ]
    

# Example
    
    
    fetch 1rx1
    fetch 1rx1, type=2fofc
    map_double 1rx1_2fofc
    isomesh mesh, 1rx1_2fofc
    

## See Also

[Map_Halve](/index.php?title=Map_Halve&action=edit&redlink=1 "Map Halve \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Map_double&oldid=13322](https://pymolwiki.org/index.php?title=Map_double&oldid=13322)"


---

## Map new

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

map_new creates a map object using one of the built-in map generation routines. 

This command can be used to create low-resolution surfaces of protein structures. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 Related Settings
  * 5 See Also



## Usage
    
    
    map_new name [, type [, grid [, selection [, buffer [, box [, state ]]]]]]
    

## Arguments

  * name = string: name of the map object to create or modify
  * type = vdw, gaussian, gaussian_max, coulomb, coulomb_neutral, coulomb_local {default: gaussian}
  * grid = float: grid spacing {default: [gaussian_resolution](/index.php?title=Gaussian_resolution&action=edit&redlink=1 "Gaussian resolution \(page does not exist\)")/3.0}
  * selection = string: atoms about which to generate the map {default: (all)}
  * buffer = float: cutoff {default: [gaussian_resolution](/index.php?title=Gaussian_resolution&action=edit&redlink=1 "Gaussian resolution \(page does not exist\)")}
  * state = integer: object state {default: 0} 
    * state > 0: use the indicated state
    * state = 0: use all states independently with independent extents
    * state = -1: use current global state
    * state = -2: use effective object state(s)
    * state = -3: use all states in one map
    * state = -4: use all states independent states by with a unified extent



## Examples

See examples for [huge surfaces](/index.php/Huge_surfaces "Huge surfaces") and [isomesh](/index.php/Isomesh "Isomesh"). 

## Related Settings

  * [gaussian_resolution](/index.php?title=Gaussian_resolution&action=edit&redlink=1 "Gaussian resolution \(page does not exist\)")
  * [gaussian_b_adjust](/index.php?title=Gaussian_b_adjust&action=edit&redlink=1 "Gaussian b adjust \(page does not exist\)")
  * [gaussian_b_floor](/index.php?title=Gaussian_b_floor&action=edit&redlink=1 "Gaussian b floor \(page does not exist\)")



## See Also

  * [map_set](/index.php/Map_set "Map set")
  * [map_trim](/index.php/Map_Trim "Map Trim")



Retrieved from "[https://pymolwiki.org/index.php?title=Map_new&oldid=9509](https://pymolwiki.org/index.php?title=Map_new&oldid=9509)"


---

## Map set

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/2/24/Map_set_ex.png)](/index.php/File:Map_set_ex.png)

[](/index.php/File:Map_set_ex.png "Enlarge")

Map_set created this consensus volume representation of aligned ligands bound to kinases. See example below.

Map_set provides a number of common operations on and between maps. For example, with Map_set you can add two maps together or create a consensus map from a series of maps or even take a difference map. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Detailed Example
  * 4 Notes
  * 5 See Also



# Usage
    
    
    map_set name, operator, operands, target_state, source_state
    

operator may be, 

  * average
  * copy
  * difference
  * maximum
  * minimum
  * sum
  * unique



# Examples
    
    
    # add 3 maps
    map_set my_sum, sum, map1 map2 map3
    
    # calculate the average map
    map_set my_avg, average, map1 map2 map3
    

# Detailed Example

This example shows you how to create a consensus map of the bound ligand in a conserved pocket. 
    
    
    # fetch some similar proteins from the PDB
    fetch 1oky 1h1w 1okz 1uu3 1uu7 1uu8 1uu9 1uvr, async=0
    
    # align them all to 1oky; their ligands
    # should all now be aligned
    alignto 1oky
    
    # highlight the ligands
    as sticks, org
    
    # select one of the atoms in the organic small mol
    sele /1uu3//A/LY4`1374/NAT
    
    # select entire molecules very near the chosen atom
    select bm. all within 1 of (sele)
    
    # remove the proteins; just look at small molecules
    remove not (sele)
    
    # create maps for all the ligands
    python
    for x in cmd.get_names():
      cmd.map_new( "map_" + x, "gaussian", 0.5, x)
    python end
    
    # calculate the average map
    map_set avgMap, average, map*
    
    # show as transparent surface
    set transparency, 0.5
    
    isosurface avgSurface, avgMap, 1.0
    
    orient vis
    

# Notes

source_state = 0 means all states 

target_state = -1 means current state 

This is an experimental function. 

  


# See Also

[Map_new](/index.php/Map_new "Map new")

Retrieved from "[https://pymolwiki.org/index.php?title=Map_set&oldid=8906](https://pymolwiki.org/index.php?title=Map_set&oldid=8906)"


---

## Map Set Border

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**map_set_border** is a function (reqd by PDA) which allows you to set the level on the edge points of a map 

### USAGE
    
    
    map_set_border <name>,<level>
    

### NOTES

unsupported. 

### SEE ALSO

[Load](/index.php/Load "Load")

Retrieved from "[https://pymolwiki.org/index.php?title=Map_Set_Border&oldid=8115](https://pymolwiki.org/index.php?title=Map_Set_Border&oldid=8115)"


---

## Map Trim

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Map Trim** introduced in PyMol 0.99beta11, gives the user the ability to cut out a section of a map and read it back in as a completely new map. This was not previously possible prior to 0.99beta11. (Windows build at [PyMol Beta](http://delsci.com/beta) or compile from source). 

## Usage
    
    
    map_trim map-name, selection-name, buffer
    

Combined with the convenient new object-name wildcards (!!!), you could for example trim all your maps to 3 Angstroms around ligands with one commmand as follows 
    
    
    map_trim *, organic, 3
    

With map size now reduced, the map_double command can come in handy to increase mesh density on your figure. 
    
    
    map_double *
    

Retrieved from "[https://pymolwiki.org/index.php?title=Map_Trim&oldid=8111](https://pymolwiki.org/index.php?title=Map_Trim&oldid=8111)"


---

## Mappend

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  


## Contents

  * 1 USAGE
  * 2 PYMOL API
    * 2.1 EXAMPLE
    * 2.2 NOTES
  * 3 SEE ALSO



### USAGE
    
    
    mappend frame : command
    

### PYMOL API

#### EXAMPLE

#### NOTES

### SEE ALSO

[Cmd mset](/index.php/Cmd_mset "Cmd mset"), [Cmd mplay](/index.php/Cmd_mplay "Cmd mplay"), [Cmd mstop](/index.php/Cmd_mstop "Cmd mstop")

Retrieved from "[https://pymolwiki.org/index.php?title=Mappend&oldid=7524](https://pymolwiki.org/index.php?title=Mappend&oldid=7524)"


---

## Mask

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mask** makes it impossible to select the indicated atoms using the mouse. This is useful when you are working with one molecule in front of another and wish to avoid accidentally selecting atoms in the background. 

### USAGE
    
    
    mask (selection)
    

### PYMOL API
    
    
    cmd.mask( string selection="(all)" )
    

### SEE ALSO

[unmask](/index.php/Unmask "Unmask"), [protect](/index.php/Protect "Protect"), [deprotect](/index.php/Deprotect "Deprotect")

Retrieved from "[https://pymolwiki.org/index.php?title=Mask&oldid=12092](https://pymolwiki.org/index.php?title=Mask&oldid=12092)"


---

## Matrix Copy

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Matrix_copy](/index.php/Matrix_copy "Matrix copy") copies the [object matrix](/index.php/Object_Matrix "Object Matrix") from one object to another. 

This command is often used after a protein structure alignment to bring other related objects into the same frame of reference. 

## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Example
  * 3 Example
  * 4 See Also



# Usage
    
    
    matrix_copy source_name, target_name
    

## Arguments

  * source_name = string: name of object to take matrix from
  * target_name = string: name(s) of object(s) to copy matrix to
  * source_mode = integer: 0: raw coordinates, 1: object TTT matrix, 2: state matrix, 3: camera matrix transformation {default: -1: [matrix_mode](/index.php/Matrix_mode "Matrix mode") setting}
  * target_mode = integer: (see source_mode)
  * source_state = integer: object state {default: 1}
  * target_state = integer: object state {default: 1}
  * target_undo = 1/0: ??? {default: 1}



# Example

Here's a practical example. We grab two proteins and their density maps. We then align one to the other and then use matrix_copy to move over the density map: 
    
    
    # fetch two proteins and their maps
    
    fetch 1rx1 3dau, async=0
    fetch 1rx1 3dau, type=2fofc, async=0
    
    # align them proteins
    
    align 1rx1, 3dau
    
    # copy 1rx1's matrix to 1rx1_2fofc
    # it's density map
    
    matrix_copy 1rx1, 1rx1_2fofc
    
    # show the result
    
    enable *
    show extent, *2fofc
    

  


# Example
    
    
    # load two molecules
    load mol1.pdb
    load mol2.pdb
    
    # load their maps
    load map1.ccp4
    load map2.ccp4
    
    # align the two molecules
    align mol2////CA, mol1////CA
    
    # copy mol2's map to mol2
    matrix_copy mol2, map2
    
    # show the isomesh
    isomesh mesh1, map1
    isomesh mesh2, map2
    

# See Also

[ Object Matrix](/index.php/Object_Matrix "Object Matrix") [Matrix_reset](/index.php/Matrix_reset "Matrix reset"), [align](/index.php/Align "Align"), [fit](/index.php/Fit "Fit"), and [pair_fit](/index.php/Pair_fit "Pair fit"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Matrix_Copy&oldid=9441](https://pymolwiki.org/index.php?title=Matrix_Copy&oldid=9441)"


---

## Matrix reset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The matrix_reset command resets the transformation for an object. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    matrix_reset name [, state [, mode ]]
    

## Arguments

  * **name** = str: object name
  * **state** = int: object state {default: 1}
  * **mode** = int: {defualt: -1 = [matrix_mode](/index.php/Matrix_mode "Matrix mode") or 0} 
    * 0: transformation was applied to coordinates
    * 1: reset [TTT matrix](/index.php/Object_Matrix "Object Matrix") (movie transformation)
    * 2: reset state matrix



## Example
    
    
    fetch 1akeA 4akeA, async=0
    
    # transform 4akeA by superposing it to 1akeA
    align 4akeA, 1akeA
    
    # undo the transformation
    matrix_reset 4akeA
    

## See Also

  * [reset](/index.php/Reset "Reset")
  * [matrix_copy](/index.php/Matrix_copy "Matrix copy")
  * [align](/index.php/Align "Align")



Retrieved from "[https://pymolwiki.org/index.php?title=Matrix_reset&oldid=12548](https://pymolwiki.org/index.php?title=Matrix_reset&oldid=12548)"


---

## Mclear

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mclear** clears the movie frame image cache. 

### USAGE
    
    
    mclear
    

### PYMOL API
    
    
    cmd.mclear()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mclear&oldid=7521](https://pymolwiki.org/index.php?title=Mclear&oldid=7521)"


---

## Mdo

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mdo** sets up a command to be executed upon entry into the specified frame of the movie. These commands are usually created by a PyMOL utility program (such as util.mrock). Command can actually contain several commands separated by semicolons ';' 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLE
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    mdo frame : command
    

### PYMOL API
    
    
     
    cmd.mdo( int frame, string command )
    

### EXAMPLE
    
    
    // Creates a single frame movie involving a rotation about X and Y
    load test.pdb
    mset 1
    mdo 1: turn x,5; turn y,5;
    mplay
    
    
    
    //Show waters within 4 Angstroms around the first residue from a 15 frame simulation trajectory
    load structure.pdb
    load_traj structure.dcd, structure, start=1, stop=15
    mset 1 -15
    for a in range(1,16): cmd.mdo(a,"hide sphere; select waters, (structure & i. 1 around 4) & resn HOH, state="+str(a)+"; show sphere, waters")
    

### NOTES

The **mset** command must first be used to define the movie before "mdo" statements will have any effect. Redefinition of the movie clears any existing mdo statements. 

### SEE ALSO

[Mset](/index.php/Mset "Mset"), [Mplay](/index.php/Mplay "Mplay"), [Mstop](/index.php/Mstop "Mstop")

Retrieved from "[https://pymolwiki.org/index.php?title=Mdo&oldid=10524](https://pymolwiki.org/index.php?title=Mdo&oldid=10524)"


---

## Mdump

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mdump** dumps the current set of movie commands. 

### USAGE
    
    
    mdump
    

### PYMOL API
    
    
    cmd.mdump()
    

### SEE ALSO

[mplay](/index.php/Mplay "Mplay"), [mset](/index.php/Mset "Mset"), [mdo](/index.php/Mdo "Mdo"), [mclear](/index.php/Mclear "Mclear"), [mmatrix](/index.php/Mmatrix "Mmatrix")

Retrieved from "[https://pymolwiki.org/index.php?title=Mdump&oldid=7519](https://pymolwiki.org/index.php?title=Mdump&oldid=7519)"


---

## Mem

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mem** dumps current memory state to standard output. This is a debugging feature, not an official part of the API. 

Retrieved from "[https://pymolwiki.org/index.php?title=Mem&oldid=7518](https://pymolwiki.org/index.php?title=Mem&oldid=7518)"


---

## Meter Reset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**meter_reset** resets the frames per secound counter 

### USAGE
    
    
    meter_reset
    

Retrieved from "[https://pymolwiki.org/index.php?title=Meter_Reset&oldid=7517](https://pymolwiki.org/index.php?title=Meter_Reset&oldid=7517)"


---

## Middle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**middle** goes to the middle of the movie. 

### USAGE
    
    
    middle
    

### PYMOL API
    
    
    cmd.middle()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Middle&oldid=7516](https://pymolwiki.org/index.php?title=Middle&oldid=7516)"


---

## Mmatrix

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mmatrix** sets up a matrix to be used for the first frame of the movie. 

### USAGE
    
    
    mmatrix {clear|store|recall}
    

### PYMOL API
    
    
    cmd.mmatrix( string action )
    

### EXAMPLES
    
    
    mmatrix store
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mmatrix&oldid=7515](https://pymolwiki.org/index.php?title=Mmatrix&oldid=7515)"


---

## Morph

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The morph command is an incentive PyMOL feature (not available in open-source version). Given two or more conformations, either as two objects or one multi-state object, **morph** can create an interpolated trajectory from the first to the second conformation. 

_This command is new in PyMOL 1.6, for older versions see[rigimol.morph](/index.php/Rigimol.morph "Rigimol.morph"), which requires some scripting._

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
    * 3.1 Default Parameters
    * 3.2 Morphing with Ligand
  * 4 Refinement
    * 4.1 Exclude atoms from refinement
    * 4.2 Select the sculpting terms
  * 5 GUI
  * 6 See Also



## Usage
    
    
    morph name, sele1 [, sele2 [, state1 [, state2 [, refinement [, steps [, method [, match ]]]]]]]
    

## Arguments

  * name = string: name of object to create
  * sele1 = string: atom selection of first conformation
  * sele2 = string: atom selection of second conformation {default: <sele1>}
  * state1 = int: sele1 state {default: 1}. If state1=0 and sele1 has N states, create N morphings between all consecutive states and back from state N to 1 (so the morph will have N*steps states). If state2=0, create N-1 morphings and stop at last state.
  * state2 = int: sele2 state {default: 2 if sele1=sele2, else 1}
  * refinement = int: number of sculpting refinement cycles to clean distorted intermediates {default: 3}
  * steps = int: number of states for sele2 object {default: 30}
  * method = string: rigimol or linear {default: rigimol}
  * match = string: method how to match **sele1** and **sele2** {default: align} 
    * in: match atoms by ["in"](/index.php/Selection_Algebra "Selection Algebra") selection operator
    * like: match atoms by ["like"](/index.php/Selection_Algebra "Selection Algebra") selection operator
    * align: match atoms with [align](/index.php/Align "Align") function (with cycles=0)
    * super: match atoms with [super](/index.php/Super "Super") function (with cycles=0)
    * _name of alignment object_ : use given alignment



## Example

### Default Parameters
    
    
    fetch 1akeA 4akeA, async=0
    align 1akeA, 4akeA
    morph mout, 1akeA, 4akeA
    

### Morphing with Ligand

Since the default **match=align** method ignores ligands in most situations, we'll use the "in" match method. In this case, both input structures need matching identifiers. The default post-RigiMOL sculpting refinement tends to let all molecules that are not connected to the protein backbone bounce around; that's why **refinement=0** might give better results in this case. 
    
    
    # get two hemoglobin beta chains, the "ligand" will be the heme hetero group
    fetch 1hbb, async=0
    create conf1, chain B
    create conf2, chain D
    delete 1hbb
    
    # important: identifiers must be the same
    alter conf2, chain="B"
    alter all, segi=""
    
    # optional: styling
    as cartoon
    show sticks, not polymer
    show nb_spheres
    
    # superpose and morph
    align conf1, conf2
    morph mout, conf1, conf2, match=in, refinement=0
    

## Refinement

The default procedure performs two steps: (1) RigiMOL, and (2) refinement of non-backbone atoms by sculpting. The second part can be skipped with "refinement=0". The refinement maintains the local sidechain geometry -- based on start and end conformation as references -- and avoids clashes by VDW repulsion. 

### Exclude atoms from refinement

Consider the previous "Morphing with Ligand" example, where refinement was skipped to avoid that the ligand gets pushed out of place. By "fixing" a set of hand-picked atoms we can still perform reasonable refinement: 
    
    
    # morph without refinement (see "Morphing with Ligand" example)
    morph mout, conf1, conf2, match=in, refinement=0
    
    # exclude the HEM core from moving
    flag fix, /mout///HEM/NA+NB+NC+ND+FE extend 1, set
    
    # do the refinement
    import epymol.rigimol
    epymol.rigimol.refine(5, "mout")
    

### Select the sculpting terms

By default, all sculpting terms are used. Just like the [sculpt_field_mask](/index.php/Sculpt_field_mask "Sculpt field mask") setting, the **refine** function takes a **sculpt_field_mask** argument which is a bitmask to select which sculpting terms should be used. For available terms, see: <https://github.com/schrodinger/pymol-open-source/blob/master/layer2/Sculpt.h>

Example which uses only bond length and angle terms: 
    
    
    cSculptBond = 0x001
    cSculptAngl = 0x002
    bondOrAngle = cSculptBond | cSculptAngl
    
    # do the refinement
    import epymol.rigimol
    epymol.rigimol.refine(5, "mout", sculpt_field_mask=bondOrAngle)
    

## GUI

The morph feature is available from the object menu panel: **A > generate > morph**

## See Also

  * [rigimol.morph](/index.php/Rigimol.morph "Rigimol.morph")
  * [flag](/index.php/Flag "Flag") fix



Retrieved from "[https://pymolwiki.org/index.php?title=Morph&oldid=12837](https://pymolwiki.org/index.php?title=Morph&oldid=12837)"


---

## Move

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**move** translates the camera along one of the three primary view axes. 

## Contents

  * 1 USAGE
    * 1.1 EXAMPLES
  * 2 PYMOL API
  * 3 SEE ALSO



### USAGE
    
    
    move axis,distance
    

#### EXAMPLES
    
    
    move x,3
    move y,-1
    

### PYMOL API
    
    
    cmd.move( string axis, float distance )
    

### SEE ALSO

[turn](/index.php/Turn "Turn"), [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [zoom](/index.php/Zoom "Zoom"), [center](/index.php/Center "Center"), [clip](/index.php/Clip "Clip"), [Translate](/index.php/Translate "Translate")

Retrieved from "[https://pymolwiki.org/index.php?title=Move&oldid=11634](https://pymolwiki.org/index.php?title=Move&oldid=11634)"


---

## Mplay

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mplay** starts the movie. 

### USAGE
    
    
    mplay
    

### PYMOL API
    
    
    cmd.mplay()
    

### SEE ALSO

[Cmd mstop](/index.php/Cmd_mstop "Cmd mstop"), [Cmd mset](/index.php/Cmd_mset "Cmd mset"), [Cmd mdo](/index.php/Cmd_mdo "Cmd mdo"), [Cmd mclear](/index.php/Cmd_mclear "Cmd mclear"), [Cmd mmatrix](/index.php/Cmd_mmatrix "Cmd mmatrix")

Retrieved from "[https://pymolwiki.org/index.php?title=Mplay&oldid=7513](https://pymolwiki.org/index.php?title=Mplay&oldid=7513)"


---

## Mpng

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mpng** writes a series of numbered movie frames to png files with the specified prefix. If the **[ray_trace_frames](/index.php/Ray_trace_frames "Ray trace frames")** variable is non-zero, these frames will be ray-traced. This operation can take several hours for a long movie. 

Be sure to disable **[cache_frames](/index.php?title=Cache_frames&action=edit&redlink=1 "Cache frames \(page does not exist\)")** when issuing this operation on a long movie (>100 frames) to avoid running out of memory. 

## Usage
    
    
    mpng prefix [, first [, last [, preserve [, modal [, mode [, quiet
        [, width [, height ]]]]]]]]
    

Options "first" and "last" can be used to specify an inclusive interval over which to render frames. Thus, you can write a smart Python program that will automatically distribute rendering over a cluster of workstations. If these options are left at zero, then the entire movie will be rendered. 

## Arguments

  * **prefix** = string: filename prefix for saved images -- output files will be numbered and end in ".png"
  * **first** = integer: starting frame {default: 0 (first frame)}
  * **last** = integer: last frame {default: 0 (last frame)}
  * **preserve** = 0/1: Only write non-existing files {default: 0}
  * **modal** = integer: will frames be rendered with a modal draw loop
  * **mode** = int: 2=ray, 1=draw, 0=normal {default: -1, check [ray_trace_frames](/index.php/Ray_trace_frames "Ray trace frames") or [draw_frames](/index.php?title=Draw_frames&action=edit&redlink=1 "Draw frames \(page does not exist\)")}
  * **width** = int: width in pixels {default: current [viewport](/index.php/Viewport "Viewport")}
  * **height** = int: height in pixels {default: current [viewport](/index.php/Viewport "Viewport")}



## Python API
    
    
    cmd.mpng(str prefix, int first=0, int last=0, int preserve=0,
        int modal=0, int mode=-1, int quiet=1,
        int width=0, int height=0)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mpng&oldid=12877](https://pymolwiki.org/index.php?title=Mpng&oldid=12877)"


---

## Mset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mset** sets up a relationship between molecular states and movie frames. This makes it possible to control which states are shown in which frame. 

The related command **madd** appends to the end of a movie. It is identical to the **mset** command, except for the default value of the **frame** argument (0 instead of 1). 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Specification Syntax
  * 4 Examples
  * 5 PyMOL API
  * 6 See Also



## Usage
    
    
    mset specification [, frame ]
    madd specification [, frame ]
    

## Arguments

  * **specification** = str: state sequence (see below for syntax), or empty string to delete movie {default: }
  * **frame** = int: start frame, or 0 to append to the end of an existing movie {default: 1 (0 for madd)}



## Specification Syntax

The state sequence specification consists of numbers, and the operators "x" and "-". Spaces between operators are optional. 

Operator semantic: 

  * _state_ **x** _count_ : Repeat state _state_ for _count_ frames
  * _state1_ **-** _state2_ : Iterate from _state1_ to _state2_ , yields _abs(state1 - state2) + 1_ frames



Formal syntax description: 

  * _specification_ ::= _state_ [**x** _count_] {**-** _specification_} [_specification_]
  * _state_ ::= **number**
  * _count_ ::= **number**



## Examples
    
    
    mset 1         # simplest case, one state -> one frame
    mset 1 x10     # ten frames, all corresponding to state 1
    mset 1 1 1 1 1 1 1 1 1 1 # identical to previous
    mset 1-20      # map a 20 state trajectory to 20 frames
    mset 1 2 3 4 5-20 # identical to previous
    
    mset 1 x30 1 -15 15 x30 15 -1
    # more realistic example:
    # the first thirty frames are state 1
    # the next 15 frames pass through states 1-15
    # the next 30 frames are of state 15
    # the next 15 frames iterate back to state 1
    
    
    
    mset 1 x200 -78 -2 -78 -2 -78 x200 79 -156 157 x200 -234 235 x400 
    # mset 1 x200 makes the first state last for 200 frames
    # -78 -2 takes us FROM state 1 to 78, then back to frame 2.  I've repeated this for dramatic effect.
    # Then we pause at 78 for 200 frames, then go from 79-156 and pause at 157 for 200 frames, etc.
    
    
    
    cmd.mset("1 -%d" % cmd.count_states())
    # this will create a one-to-one mapping of states to movie frames. useful for making movies from trajectory files.
    

## PyMOL API
    
    
    cmd.mset(string specification, int frame=1)
    cmd.madd(string specification, int frame=0)
    

## See Also

[mdo](/index.php/Mdo "Mdo"), [mplay](/index.php/Mplay "Mplay"), [mclear](/index.php/Mclear "Mclear")

Retrieved from "[https://pymolwiki.org/index.php?title=Mset&oldid=12179](https://pymolwiki.org/index.php?title=Mset&oldid=12179)"


---

## Mstop

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**mstop** stops the movie. 

### USAGE
    
    
    mstop
    

### PYMOL API
    
    
    cmd.mstop()
    

### SEE ALSO

[Cmd mplay](/index.php/Cmd_mplay "Cmd mplay"), [Cmd mset](/index.php/Cmd_mset "Cmd mset"), [Cmd mdo](/index.php/Cmd_mdo "Cmd mdo"), [Cmd mclear](/index.php/Cmd_mclear "Cmd mclear"), [Cmd mmatrix](/index.php/Cmd_mmatrix "Cmd mmatrix")

Retrieved from "[https://pymolwiki.org/index.php?title=Mstop&oldid=7499](https://pymolwiki.org/index.php?title=Mstop&oldid=7499)"


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

## Multisave

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The multisave command can save a multi-entry PDB file. 

Every object in the given selection will have a HEADER record, and is terminated with END. Loading such a multi-entry PDB file into PyMOL will load each entry as a separate object. 

This behavior is different to the [save](/index.php/Save "Save") command, where a multi-object selection is written "flat" to a PDB file, without HEADER records. 

_Older versions could also export a proprietary binary "pmo" file format with this command. However, "pmo" format never supported multiple entries, support was partially broken in in the 1.7 and 1.8 versions, and was finally removed in 1.8.4._

## Contents

  * 1 Changes in PyMOL 1.8.4
  * 2 Usage
  * 3 Arguments
  * 4 Example
  * 5 See Also



## Changes in PyMOL 1.8.4

  * CRYST1 records are written if symmetry is defined
  * Interpret "pattern" argument as an atom selection. Previous versions expected an object name pattern or a list of objects.
  * dropped "pmo" format support



## Usage
    
    
    multisave filename [, pattern [, state [, append ]]]
    

## Arguments

  * **filename** = str: file path to be written
  * **pattern** = str: atom selection (before 1.8.4: object name pattern)
  * **state** = int: object state (-1=current, 0=all) {default: -1}
  * **append** = 0/1: append to existing file {default: 0}



## Example
    
    
    fragment ala
    fragment gly
    
    # PyMOL >= 1.8.4 will write CRYST1 records
    set_symmetry ala, 50.840, 42.770,  28.950, 90, 90, 90, P 21 21 21
    set_symmetry gly, 99.930, 99.930, 276.840, 90, 90, 90, P 31 2 1
    
    multisave multi.pdb
    

Exported file "multi.pdb": 
    
    
    HEADER    ala
    CRYST1   50.840   42.770   28.950  90.00  90.00  90.00 P 21 21 21    0
    ATOM      1  N   ALA     2      -0.677  -1.230  -0.491  1.00  0.00           N
    ATOM      2  CA  ALA     2      -0.001   0.064  -0.491  1.00  0.00           C
    ATOM      3  C   ALA     2       1.499  -0.110  -0.491  1.00  0.00           C
    ATOM      4  O   ALA     2       2.030  -1.227  -0.502  1.00  0.00           O
    ATOM      5  CB  ALA     2      -0.509   0.856   0.727  1.00  0.00           C
    ATOM      6  H   ALA     2      -0.131  -2.162  -0.491  1.00  0.00           H
    ATOM      7  HA  ALA     2      -0.269   0.603  -1.418  1.00  0.00           H
    ATOM      8 1HB  ALA     2      -1.605   1.006   0.691  1.00  0.00           H
    ATOM      9 2HB  ALA     2      -0.285   0.342   1.681  1.00  0.00           H
    ATOM     10 3HB  ALA     2      -0.053   1.861   0.784  1.00  0.00           H
    END
    HEADER    gly
    CRYST1   99.930   99.930  276.840  90.00  90.00  90.00 P 31 2 1      0
    ATOM      1  N   GLY    22      -1.195   0.201  -0.206  1.00  0.00           N
    ATOM      2  CA  GLY    22       0.230   0.318  -0.502  1.00  0.00           C
    ATOM      3  C   GLY    22       1.059  -0.390   0.542  1.00  0.00           C
    ATOM      4  O   GLY    22       0.545  -0.975   1.499  1.00  0.00           O
    ATOM      5  H   GLY    22      -1.558  -0.333   0.660  1.00  0.00           H
    ATOM      6 3HA  GLY    22       0.482   1.337  -0.514  0.00  0.00           H
    ATOM      7  HA  GLY    22       0.434  -0.159  -1.479  1.00  0.00           H
    END
    

## See Also

  * [save](/index.php/Save "Save")
  * [load](/index.php/Load "Load")



Retrieved from "[https://pymolwiki.org/index.php?title=Multisave&oldid=12535](https://pymolwiki.org/index.php?title=Multisave&oldid=12535)"


---

## Mview

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The mview command can store and delete movie keyframes. 

Keyframes store a view (camera or object position) and optionally the object [state](/index.php/State "State") and/or a [scene](/index.php/Scene "Scene"). Between keyframes, PyMOL will interpolate views and states, allowing for smooth animations. 

Before using mview, the movie timeline has to be set up with [mset](/index.php/Mset "Mset"). 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
    * 3.1 Ligand zoom
    * 3.2 360° rotation
    * 3.3 360° rotation of a single object
    * 3.4 Object-level state-sweep
    * 3.5 Ligand binding
    * 3.6 Scene based movie
  * 4 See Also



## Usage
    
    
    mview [ action [, first [, last [, power [, bias
        [, simple [, linear [, object [, wrap [, hand
        [, window [, cycles [, scene [, cut [, quiet
        [, auto [, state [, freeze ]]]]]]]]]]]]]]]]]]
    

## Arguments

  * **action** = str: one of store, clear, reset, purge, interpolate, uninterpolate, reinterpolate, toggle, toggle_interp, smooth {default: store}
  * **first** = int: frame number or 0 for current frame {default: 0}
  * **power** = float: slow down animation at keyframe (0.0) or not (1.0) {default: 0.0}
  * **object** = str: name of object for object keyframes, or empty for global (camera) keyframes {default: }
  * **scene** = str: name of scene to store scene with key frame {default: }
  * **cut** = float 0.0-1.0: scene switch moment (0.0: beginning of transition, 1.0: end of transition) {default: 0.5}
  * **auto** = -1/0/1: if freeze=0, then auto reinterpolate after store, clear, or toggle {default: -1 = use [movie_auto_interpolate](/index.php?title=Movie_auto_interpolate&action=edit&redlink=1 "Movie auto interpolate \(page does not exist\)")}
  * **state** = int: if > 0, then store object state {default: 0}
  * **freeze** = 0/1: never auto reinterpolate {default: 0}



## Examples

### Ligand zoom
    
    
    fetch 1rx1, async=0
    as cartoon
    as sticks, organic
    mset 1x70
    orient
    mview store, 1
    mview store, 70
    orient organic
    mview store, 30
    mview store, 40
    mplay
    

### 360° rotation
    
    
    fragment ala
    mset 1x90
    mview store, 1
    mview store, 90
    turn y, 120
    mview store, 30, power=1.0
    turn y, 120
    mview store, 60, power=1.0
    mplay
    

### 360° rotation of a single object
    
    
    set movie_auto_store, 0
    fragment ala
    fragment his
    translate [10, 0, 0], his
    zoom
    mset 1x90
    mview store, 1, object=his
    mview store, 90, object=his
    rotate y, 120, object=his  # keyword >>object<< is absolutely necessary! Otherwise movie will not work.
    mview store, 30, power=1.0, object=his
    rotate y, 120, object=his
    mview store, 60, power=1.0, object=his
    mplay
    

### Object-level state-sweep
    
    
    load <http://pymol.org/tmp/morph.pdb.gz>
    dss
    as cartoon
    mset 1x80
    mview store, 1, object=m, state=1
    mview store, 30, object=m, state=30
    mview store, 50, object=m, state=30
    mview store, 80, object=m, state=1
    mplay
    

### Ligand binding
    
    
    set movie_auto_store, 0
    fetch 1rx1, async=0
    extract ligand, organic
    as cartoon, 1rx1
    as sticks, ligand
    set_view (\
      0.527486444, -0.761115909, -0.377440333,\
      0.736519873, 0.631122172, -0.243357301,\
      0.423434794, -0.149625391, 0.893482506,\
      0.000059791, -0.000049331, -140.287048340,\
      34.670463562, 51.407436371, 17.568315506,\
      111.284034729, 169.290832520, -19.999998093 )
    mset 1x60
    mview store, 60, object=ligand
    translate [10, 0, 0], object=ligand
    mview store, 1, object=ligand
    mplay
    

### Scene based movie
    
    
    fragment ala
    as sticks
    color blue
    scene bluesticks, store
    as spheres
    color red
    turn y, 180
    scene redspheres, store
    mset 1x60
    mview store, 1, scene=bluesticks
    mview store, 30, scene=redspheres
    mplay
    

## See Also

  * [mset](/index.php/Mset "Mset")
  * [frame](/index.php/Frame "Frame")
  * [mdo](/index.php/Mdo "Mdo")
  * [mpng](/index.php/Mpng "Mpng")
  * [movie_loop](/index.php/Movie_loop "Movie loop") setting
  * [movie_auto_interpolate](/index.php?title=Movie_auto_interpolate&action=edit&redlink=1 "Movie auto interpolate \(page does not exist\)") setting
  * [movie_auto_store](/index.php?title=Movie_auto_store&action=edit&redlink=1 "Movie auto store \(page does not exist\)") setting



Retrieved from "[https://pymolwiki.org/index.php?title=Mview&oldid=12641](https://pymolwiki.org/index.php?title=Mview&oldid=12641)"


---

## Order

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 DESCRIPTION
  * 2 USAGE
    * 2.1 EXAMPLES
  * 3 PYMOL API
  * 4 NOTES
  * 5 SEE ALSO



### DESCRIPTION

**order** allows you to change ordering of names in the control panel 

### USAGE
    
    
    order names-list, sort, location
    

#### EXAMPLES
    
    
    # sets the order of these three objects
    order 1dn2 1fgh 1rnd
    
    # sorts all names
    order *,yes
    
    # sorts all names beginning with 1dn2_
    order 1dn2_*, yes
    
    # puts 1frg at the top of the list
    order 1frg, location=top
    

### PYMOL API
    
    
    cmd.order(string names-list, string sort, string location)
    

### NOTES

  1. names-list: a space separated list of names
  2. sort: yes or no
  3. location: top, current, or bottom



### SEE ALSO

[Set_Name](/index.php/Set_Name "Set Name")

Retrieved from "[https://pymolwiki.org/index.php?title=Order&oldid=6769](https://pymolwiki.org/index.php?title=Order&oldid=6769)"


---

## Orient

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**orient** aligns the principal components of the atoms in the selection with the XYZ axes. The function is similar to the orient command in X-PLOR. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 EXAMPLES
  * 5 SEE ALSO



### USAGE
    
    
    orient object-or-selection [, state]
    orient (selection)
    

### PYMOL API
    
    
    cmd.orient( string object-or-selection [, state = 0] )
    

### NOTES
    
    
      state = 0 (default) use all coordinate states
      state = -1 use only coordinates for the current state
      state > 0  use coordinates for a specific state
    

### EXAMPLES

For models with NCS symmetry, orient will align the model with the symmetry axis centered along the viewport's z axis. For example, 
    
    
      fetch 1hiw, async=0
      as cartoon
      remove (!chain A,B,C)
      orient
      util.cbc
    

will produce the first image below. However, if there is a larger symmetry, e.g. two trimers, this will not work. In the above example, leaving out `remove (!chain A,B,C)` from the script results in the second image below. 

  * [![One trimer from 1hiw after "orient" command.](/images/8/8d/1hiw_orient.png)](/index.php/File:1hiw_orient.png "One trimer from 1hiw after "orient" command.")

One trimer from 1hiw after "orient" command. 

  * [![Both trimers from 1hiw after "orient" command.](/images/8/86/1hiw_orient2.png)](/index.php/File:1hiw_orient2.png "Both trimers from 1hiw after "orient" command.")

Both trimers from 1hiw after "orient" command. 




### SEE ALSO

[Zoom](/index.php/Zoom "Zoom"), [Origin](/index.php/Origin "Origin"), [Reset](/index.php/Reset "Reset")

Retrieved from "[https://pymolwiki.org/index.php?title=Orient&oldid=8464](https://pymolwiki.org/index.php?title=Orient&oldid=8464)"


---

## Origin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**origin** sets the center of rotation about a selection. If an object name is specified, it can be used to set the center of rotation for the object's [TTT matrix](/index.php/Object_Matrix "Object Matrix"). 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO
  * 5 Example



### USAGE
    
    
    origin selection [, object [,position, [, state]]]
    origin (selection)
    origin position=[1.0,2.0,3.0]
    

### PYMOL API
    
    
    cmd.origin( string object-or-selection )
    

### NOTES

  * state = 0 (default) use all coordinate states
  * state = -1 use only coordinates for the current state
  * state > 0 use coordinates for a specific state



### SEE ALSO

[Zoom](/index.php/Zoom "Zoom"), [Orient](/index.php/Orient "Orient"), [Reset](/index.php/Reset "Reset")

  


### Example

This example puts the camera 'inside' a protein near the ligand and turns the camera 360 degrees about that location. 
    
    
    load $TUT/1hpv.pdb
    
    set field_of_view, 45
    
    zoom organic
    orient organic
    show stick, polymer
    show surface
    set surface_color, white
    set transparency, 0.5
    clip slab, 30
    set line_width, 5
    
    cmd.move("z",-cmd.get_view()[11])
    
    python
    def spin():
       from pymol import cmd
       for x in range(360):
          cmd.move("z",cmd.get_view()[11])
          cmd.turn("y",1)
          cmd.move("z",-cmd.get_view()[11])
          cmd.refresh()
    
    threading.Thread(target=spin).start()
    python end
    

Retrieved from "[https://pymolwiki.org/index.php?title=Origin&oldid=9541](https://pymolwiki.org/index.php?title=Origin&oldid=9541)"


---

## Pair fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Pair_Fit](/index.php/Pair_Fit "Pair Fit") fits a set of atom pairs between two models. Each atom in each pair must be specified individually, which can be tedious to enter manually. Script files are recommended when using this command. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 NOTES
  * 4 USER EXAMPLES/COMMENTS
  * 5 SEE ALSO



### USAGE
    
    
    pair_fit (selection), (selection), [ (selection), (selection) [ ...] ]
    

### EXAMPLES
    
    
    # superimpose protA residues 10-25 and 33-46 to protB residues 22-37 and 41-54:
    pair_fit protA///10-25+33-46/CA, protB///22-37+41-54/CA
    # superimpose ligA atoms C1, C2, and C4 to ligB atoms C8, C4, and C10, respectively:
    pair_fit ligA////C1, ligB////C8, ligA////C2, ligB////C4, ligA////C3, ligB////C10
    

### NOTES

So long as the atoms are stored in PyMOL with the same order internally, you can provide just two selections. Otherwise, you may need to specify each pair of atoms separately, two by two, as additional arguments to pair_fit. 

Script files are usually recommended when using this command. 

### USER EXAMPLES/COMMENTS

An description of selection caveats for these commands may be found at [Rms](/index.php/Rms "Rms"). 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur")

Retrieved from "[https://pymolwiki.org/index.php?title=Pair_fit&oldid=11478](https://pymolwiki.org/index.php?title=Pair_fit&oldid=11478)"


---

## Png

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**png** writes a png format image file of the current image to disk. 

## Contents

  * 1 Usage
  * 2 Example
  * 3 PyMOL API
  * 4 Comments
    * 4.1 Transparent Backgrounds
    * 4.2 DPI Setting



## Usage
    
    
    png filename[, width[, height[, dpi[, ray[, quiet]]]]]
    

  * **filename** = string: file path to be written
  * **width** = integer or string: width in pixels (integer or string without units), inches (in), or centimeters (cm). If unit suffix is given, `dpi` argument is required as well. If only one of `width` or `height` is given, the aspect ratio of the viewport is preserved. {default: 0 (current)}
  * **height** = integer or string: height (see width) {default: 0 (current)}
  * **dpi** = float: dots-per-inch {default -1.0 (unspecified)}
  * **ray** = 0 or 1: should ray be run first {default: 0 (no)}
  * **quiet** = 0 or 1: if 1, logged output is suppressed. {default: 0}



## Example
    
    
    png ~/Desktop/test.png, width=10cm, dpi=300, ray=1
    

## PyMOL API
    
    
    cmd.png(string filename, int width=0, int height=0, float dpi=-1, int ray=0, int quiet=0)
    

## Comments

#### Transparent Backgrounds

Use the `[ray_opaque_background](/index.php/Ray_opaque_background "Ray opaque background")` setting to output images with transparent backgrounds. 
    
    
    set ray_opaque_background, 0
    

This can be useful for presentations, images that are placed on top of a background of nonuniform color (e.g. gradients), and images that overlap text or other images. 

#### DPI Setting

Use the DPI option to have PyMol set the DPI of your image. Executing the command 
    
    
    png /tmp/ex.png, width=1200, height=1200, dpi=300, ray=1
    

will ouput a four-inch square image at 300dpi. Leaving off the **dpi** parameter would yield a 1200x1200 image at your system's default pixel density (e.g. 72 or 96 dpi). This saves the intermediate step of having to use GIMP/PhotoShop/etc to rescale your photos for publication. 

Retrieved from "[https://pymolwiki.org/index.php?title=Png&oldid=11924](https://pymolwiki.org/index.php?title=Png&oldid=11924)"


---

## Protect

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**protect** protects a set of atoms from tranformations performed using the editing features. This is most useful when you are modifying an internal portion of a chain or cycle and do not wish to affect the rest of the molecule. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 Example
  * 4 SEE ALSO



## USAGE
    
    
    protect (selection)
    

## PYMOL API
    
    
    cmd.protect(string selection)
    

## Example

This example makes a very cool little, and fake, molecular movie. Copy/paste this into PyMOL: 
    
    
    load $PYMOL_PATH/test/dat/pept.pdb, obj
    
    # create the fake trajectory (of states)
    for a in range(2,31): cmd.create("obj","obj",1,a)
    # remove the bond
    unbond 5/C, 6/N
    # This protects everything but all atoms witing 4
    # Ang. of residue 5's carbon and residue 6's nitrogen.
    protect not ((5/C or 6/N) extend 4)
    
    # do some quick sculpting
    sculpt_activate obj, 30
    sculpt_iterate obj, 30, 100
    smooth obj, 30, 3
    
    # then program a bond-break/re-form movie
    mset 1 x30 1 -30 30 x30 30 -1 
    mdo 45: unbond 5/C, 6/N, quiet=1
    mdo 100: bond 5/C, 6/N, quiet=1
    
    frame 100
    
    as sticks
    orient 5-6/N+CA+C
    mplay
    

## SEE ALSO

[Deprotect](/index.php/Deprotect "Deprotect"), [Mask](/index.php/Mask "Mask"), [Unmask](/index.php/Unmask "Unmask")

Retrieved from "[https://pymolwiki.org/index.php?title=Protect&oldid=12275](https://pymolwiki.org/index.php?title=Protect&oldid=12275)"


---

## Pseudoatom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

pseudoatom creates a molecular object with a pseudoatom or adds a pseudoatom to a molecular object if the specified object already exists. Default position is in the middle of the viewing window. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES FOR USE
  * 4 References



## USAGE
    
    
    pseudoatom object [, selection [, name [, resn [, resi [, chain
            [, segi [, elem [, vdw [, hetatm [, b [, q [, color [, label
            [, pos [, state [, mode [, quiet ]]]]]]]]]]]]]]]]]
    

## ARGUMENTS

  * object = string: name of object to create or modify


  * selection = string: optional atom selection. If given, calculate position (pos) as center of this selection (like [center of mass](/index.php/Center_Of_Mass "Center Of Mass") without atom masses).


  * name, resn, resi, chain, segi, elem, vdw, hetatm, b, q = string: [atom properties](/index.php/Property_Selectors "Property Selectors") of pseutoatom


  * color = string: [color](/index.php/Color "Color") of pseudoatom


  * label = string: place a [text label](/index.php/Label "Label") and hide all other representations


  * pos = 3-element tuple of floats: location in space (only if no selection given)


  * state = integer: [state](/index.php/State "State") to modify, 0 for current state, -1 for all states {default: 0}


  * mode = string: determines the vdw property if vdw is not given. Either as RMS distance (mode=rms, like [radius of gyration](/index.php/Radius_of_gyration "Radius of gyration") without atom masses) or as maximum distance (mode=extent) from pseudoatom to selection {default: rms}


    
    
    # create the pseudoatom
    pseudoatom tmpPoint
     ObjMol: created tmpPoint/PSDO/P/PSD`1/PS1
    # show it as a sphere.
    show spheres, tmpPoint
    
    # create another, with more options.
    pseudoatom tmpPoint2, resi=40, chain=ZZ, b=40, color=tv_blue, pos=[-10, 0, 10]
     ObjMol: created tmpPoint2/PSDO/ZZ/PSD`40/PS1
    

## EXAMPLES FOR USE

pseudoatom can be used for a wide variety of tasks where on must place an atom or a label in 3D space, e.g. as a placeholder for distance measurement or distance specifications. 
    
    
    # A pseudoatom as a placeholder for selections according to distance:
    load $TUT/1hpv.pdb
    pseudoatom tmp, pos=[10.0, 17.0, -3.0]
    show sticks, tmp expand 6
    delete tmp
    
    # A pseudoatom as placeholder for distance measurement: 
    # position it at the center of an aromatic ring.  Then 
    # calc the distance from another atom to the pseudoatom.
    load $TUT/1hpv.pdb
    pseudoatom pi_cent,b/53/cg+cz
    dist pi_cent////ps1, b/met`46/ce
    

You can use a pseudoatom to make a label for a scene title. Move the protein to the bottom of the window before the pseudoatom is created. Or move the label after creating it (Shift + Middle mouse button in editing mode). 
    
    
    fetch 1rq5
    pseudoatom forLabel
    label forLabel, "This Protein is a Carbohydrate Active Enzyme"
    set label_color, black
    # png ray=1
    

  * [![You can use a pseudoatom for a movable scene label.](/images/3/32/Patom.png)](/index.php/File:Patom.png "You can use a pseudoatom for a movable scene label.")

You can use a pseudoatom for a movable scene label. 

  * [![Example of using a pseudoatom for distance measurement.](/images/3/33/Pseu1.png)](/index.php/File:Pseu1.png "Example of using a pseudoatom for distance measurement.")

Example of using a pseudoatom for distance measurement. 




  


  


  


## References

PyMOL Mailing List 

Retrieved from "[https://pymolwiki.org/index.php?title=Pseudoatom&oldid=9090](https://pymolwiki.org/index.php?title=Pseudoatom&oldid=9090)"


---

## Push Undo

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**push_undo** stores the currently conformations of objects in the selection onto their individual kill rings. 

### USAGE
    
    
    push_undo (all)
    

### SEE ALSO

[Undo](/index.php/Undo "Undo"), [Redo](/index.php/Redo "Redo")

Retrieved from "[https://pymolwiki.org/index.php?title=Push_Undo&oldid=7491](https://pymolwiki.org/index.php?title=Push_Undo&oldid=7491)"


---

## Pwd

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Print current working directory. 

### USAGE
    
    
    pwd
    

### SEE ALSO

[ls](/index.php/Ls "Ls"), [cd](/index.php/Cd "Cd"), [system](/index.php/System "System")

Retrieved from "[https://pymolwiki.org/index.php?title=Pwd&oldid=7490](https://pymolwiki.org/index.php?title=Pwd&oldid=7490)"


---

## Python

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Issuing the **Python** command will put you into a stateful pseudo-interactive Python session. Or, more simply it's stateful in that you can invoke the Python session write some code, end the session, then restart the session and your data will be saved (see Example 1). It's pseudo-interactive in that you don't get feedback until you type "python end," upon which your code is run the output appears. 

This is a helpful command for testing different scripting or state-editing strategies for movie making. 

## USAGE
    
    
    # start the session
    python
    
    # ...
    # your Python code goes here
    # ...
    
    # end the session
    python end
    

  


## EXAMPLES

  * Start the session. Set x to 10. End the session. Restart the session and see if the value of x is recalled.


    
    
    python
    x = 10
    print(x)
    python end
    
    
    
    python
    print(x)
    python end
    

Output: 
    
    
    10
    

## Python Version

Python scripts and commands used within PyMOL can only be written using the current version of Python that is supported by your version of PyMOL. To determine which version of Python you can use, type the following command into PyMOL: 
    
    
    print(sys.version)
    

Note that this version of Python is not necessarily related to the version that you may have installed on your system. 

This command can also be used to ensure that code you are distributing can be supported by the user's system. 

Retrieved from "[https://pymolwiki.org/index.php?title=Python&oldid=13242](https://pymolwiki.org/index.php?title=Python&oldid=13242)"


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

## Ramp update

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

ramp_update updates range and/or color of a color ramp. 

_New in PyMOL 1.8_

## Usage
    
    
    ramp_update name [, range [, color ]]
    

## Example
    
    
    ramp_new    e_pot_color, e_pot_map, [-10,0,10], [red,white,blue]
    ramp_update e_pot_color, range=[-15,0,15]
    ramp_update e_pot_color, color=[green,white,orange]
    

## See Also

  * [ramp_new](/index.php/Ramp_new "Ramp new")



Retrieved from "[https://pymolwiki.org/index.php?title=Ramp_update&oldid=13485](https://pymolwiki.org/index.php?title=Ramp_update&oldid=13485)"


---

## Ray

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**ray** creates a ray-traced image of the current frame. 

[![](/images/b/bc/Gslike.png)](/index.php/File:Gslike.png)

[](/index.php/File:Gslike.png "Enlarge")

Varying settings to play with rendering options

## Contents

  * 1 Details
    * 1.1 Usage
    * 1.2 PyMol API
    * 1.3 Settings
      * 1.3.1 Modes
      * 1.3.2 Perspective
        * 1.3.2.1 Perspective Example Images
          * 1.3.2.1.1 Notes
      * 1.3.3 Renderer
    * 1.4 Performance
      * 1.4.1 Memory
    * 1.5 Examples
      * 1.5.1 Simple
      * 1.5.2 Specify Image Size
      * 1.5.3 Specify Renderer
      * 1.5.4 High Quality B&W Rendering
      * 1.5.5 High Quality Color
      * 1.5.6 Ray Tracing Stereo Images
    * 1.6 See also
    * 1.7 User comments



# Details

This command is used to make high-resolution photos fit for publication and formal movies. Please note, the **ray** command can take some time (up to several minutes, depending on image complexity and size). 

For those who are making movies with PyMOL, **Ray** is one of the most commonly used last steps before stitching the frames together to compile the movie. **Ray** has many powerful features such as setting the size of the image -- and it still works even if the [Viewport](/index.php/Viewport "Viewport") or screen is smaller than the requested output file dimensions. 

[![](/images/b/b9/No_ray_trace.png)](/index.php/File:No_ray_trace.png)

[](/index.php/File:No_ray_trace.png "Enlarge")

Image, not ray traced.

[![](/images/3/30/Ray_traced.png)](/index.php/File:Ray_traced.png)

[](/index.php/File:Ray_traced.png "Enlarge")

Image, ray traced.

## Usage
    
    
    ray [width,height [,renderer [,angle [,shift ]]]
    

**angle** and **shift** can be used to generate matched stereo pairs 

**width** and **height** can be set to any non-negative integer. If both are set to zero than the current window size is used and is equivalent to just using **ray** with no arguments. If one is set to zero (or missing) while the other is a positive integer, then the argument set to zero (or missing) will be scaled to preserve the current aspect ratio. 

## PyMol API
    
    
    cmd.ray(int width,int height,int renderer=-1,float shift=0)
    

## Settings

### Modes

Setting the **[Ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")** variable in PyMOL changes the way PyMOL's internal renderer represents proteins in the final output. New modes were recently added to give the user more options of molecular representation. New modes are: normal rendering, but with a black outline (nice for presentations); black and white only; quantized color with black outline (also, very nice for presentations; more of a _cartoony_ appearance). 

**Note:** Mode 3, the quantized color one, sort of **burns** the background if you're using this setting. This will make a pure white background somewhat "offwhite"; thus, a poster would look poor because you could see the border for the image. If you'll be using this mode, try the [ray_opaque_background](/index.php/Ray_opaque_background "Ray opaque background") setting. 
    
    
    # normal color
    set ray_trace_mode, 0
    
    # normal color + black outline
    set ray_trace_mode,  1
    
    # black outline only
    set ray_trace_mode,  2
    
    # quantized color + black outline
    set ray_trace_mode,  3
    
    set ray_trace_mode, 1 # (or 2 or 3; best with "bg_color white;set antialias,2")
    # These two new modes -- 2 and 3 -- are cool cartoon looking modes.
    
    # change the color of the outline to a named color, or a hex-code
    set ray_trace_color, magenta
    set ray_trace_color, 0x0033ff
    

Here are the example images for the new modes 

  * [![set ray_trace_mode,1](/images/a/a8/Ray_mode_1_ex.png)](/index.php/File:Ray_mode_1_ex.png "set ray_trace_mode,1")

set ray_trace_mode,1 

  * [![set ray_trace_mode,2](/images/9/95/Ray_mode_2_ex.png)](/index.php/File:Ray_mode_2_ex.png "set ray_trace_mode,2")

set ray_trace_mode,2 

  * [![set ray_trace_mode,3](/images/5/51/Ray_mode_3_ex.png)](/index.php/File:Ray_mode_3_ex.png "set ray_trace_mode,3")

set ray_trace_mode,3 




### Perspective

#### Perspective Example Images

  * [![Example with Perspective Turned Off](/images/3/31/No_persp.png)](/index.php/File:No_persp.png "Example with Perspective Turned Off")

Example with Perspective Turned Off 

  * [![Example with Perspective Turned On](/images/c/cc/Persp1.png)](/index.php/File:Persp1.png "Example with Perspective Turned On")

Example with Perspective Turned On 

  * [![Example with Perspective Turned On and Field of View Set High \(70\).](/images/8/88/Persp2.png)](/index.php/File:Persp2.png "Example with Perspective Turned On and Field of View Set High \(70\).")

Example with Perspective Turned On and Field of View Set High (70). 




##### Notes

PyMol 0.97 and prior used **orthoscopic** rendering -- that is, no perspective. Upon the arrival of 0.98 and later, we get perspective based rendering at a cost of a 4x decrease in render speed. If you want perspective 
    
    
    set orthoscopic, off
    

Otherwise 
    
    
    set orthoscopic, on
    

To magnify the effect of perspective on the scene, 
    
    
    set field_of_view, X
    

where 50<X<70\. Default is 20. 50-70 gives a very strong perspective effect. Nb. the field of view is in Y, not X as one would expect. 

  


### Renderer

**renderer = -1** is default (use value in ray_default_renderer) 

**renderer = 0** uses PyMOL's internal renderer 

**renderer = 1** uses PovRay's renderer. This is Unix-only and you must have "povray" in your path. It utilizes two temporary files: "tmp_pymol.pov" and "tmp_pymol.png". Also works on Mac via Povray37UnofficialMacCmd but it needs to be in your path as "povray". 

## Performance

  * The ray performance depends on distance between camera and molecule.



If the distance is big rendering takes much time. If the distance is too small distant parts of molecule dissolve. 

  * [![Too close to molecule](/images/7/70/Close_ray.png)](/index.php/File:Close_ray.png "Too close to molecule")

Too close to molecule 

  * [![Normal distance](/images/1/1c/Middle_ray.png)](/index.php/File:Middle_ray.png "Normal distance")

Normal distance 



  * Tip: If you have a rather complicated scene that is zoomed into only a part of the molecule, you can speed up the ray tracing by hiding everything else outside of a certain range of the zoomed-on point. For example, if I have a large molecule and I'm looking only at the 30-atom ligand bound to it, then I can do something like the following:


    
    
    # setup your complex scene
    ...
    
    # zoom on the hetero atom (ligand and not water) within 5 Angstroms
    select hh, het and not resn HOH
    zoom hh, 5
    
    # turn on depth cueing
    set depth_cue, 1
    
    # now, select stuff to hide; we select everything that is 
    # farther than 8 Ang from our main selection
    select th, (all) and not ( (all) within 8 of hh) )
    
    hide everything, th
    
    # any additional commands you want
    ...
    
    ray
    

As an example of the efficacy of this method, I ray traced a rather complex scene with all the atoms visible here's the output of ray: 
    
    
    PyMOL>ray
     Ray: render time: 24.50 sec. = 146.9 frames/hour (941.88 sec. accum.).
    

and here is the result when I soft-clipped everything else using the above method: 
    
    
    PyMOL>ray
     Ray: render time: 47.93 sec. = 75.1 frames/hour (989.80 sec. accum.).
    

The two images in the following gallery show the results of the rendering. 

  * [![normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI.](/images/1/1a/Ray_method_off.png)](/index.php/File:Ray_method_off.png "normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI.")

normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI. 

  * [![manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left.](/images/7/7d/Ray_method_on.png)](/index.php/File:Ray_method_on.png "manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left.")

manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left. 




### Memory

If memory is an issue for you in PyMOL, try executing your rendering from a script rather than a PyMOL session file. An unfortunate unavoidable consequence of the fact that we use Python's portable, platform-independent "pickle" machinery for PyMOL session files. Packing or unpacking a Session or Scene file thus requires that there be two simultanous copies of the information to reside in RAM simultaneously: one native and a second in Python itself. 

So when memory is a limiting factor, scripts are recommended over sessions. 

## Examples

### Simple
    
    
    # ray trace the current scene using the default size of the viewport
    ray
    

### Specify Image Size
    
    
    # ray trace the current scene, but scaled to 1024x768 pixels
    ray 1024,768
    

### Specify Renderer
    
    
    # ray trace with an external renderer.
    ray renderer=0
    

### High Quality B&W Rendering

[![](/images/6/6c/1l9l.png)](/index.php/File:1l9l.png)

[](/index.php/File:1l9l.png "Enlarge")

Black and White (ray_trace_mode,2); click to see full image
    
    
    # Black and White Script
    load /tmp/3fib.pdb;
    show cartoon;
    set ray_trace_mode, 2;  # black and white cartoon
    bg_color white;
    set antialias, 2;
    ray 600,600
    png /tmp/1l9l.png
    

### High Quality Color

[![](/images/7/7e/1l9l_2.png)](/index.php/File:1l9l_2.png)

[](/index.php/File:1l9l_2.png "Enlarge")

Color mode (ray_trace_mode,3); click to see full image
    
    
    # Color Script
    load /tmp/thy_model/1l9l.pdb;
    hide lines;
    show cartoon;
    set ray_trace_mode, 3; # color
    bg_color white;
    set antialias, 2;
    remove resn HOH
    remove resn HET
    ray 600,600
    png /tmp/1l9l.png
    

### Ray Tracing Stereo Images

    _See[Stereo_Ray](/index.php/Stereo_Ray "Stereo Ray")_

## See also

  1. "help faster" for optimization tips with the builtin renderer. "help povray" for how to use PovRay instead of PyMOL's built-in ray-tracing engine. For high-quality photos, please also see the [Antialias](/index.php/Antialias "Antialias") command. [Ray shadows](/index.php/Ray_shadows "Ray shadows") for controlling shadows.
  2. See also [Ray Tracing](/index.php/Ray_Tracing "Ray Tracing").
  3. [Desaturation Tutorial](http://www.gimp.org/tutorials/Color2BW) \-- A good resource for making nice B&W images from color images (desaturation).
  4. [Ray Trace Gain](/index.php/Ray_Trace_Gain "Ray Trace Gain")



## User comments

How do I ray trace a publication-ready (~300dpi) image using PyMol?
    This answer is in the [Advanced Issues](/index.php/Category:Advanced_Issues "Category:Advanced Issues") (Image Manipulation Section).

Retrieved from "[https://pymolwiki.org/index.php?title=Ray&oldid=11725](https://pymolwiki.org/index.php?title=Ray&oldid=11725)"


---

## Read Molstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**read_molstr** reads an MDL MOL format file as a string 

### PYMOL API ONLY
    
    
    cmd.read_molstr( string molstr, string name, int state=0, \
      int finish=1, int discrete=1 )
    

### NOTES

  * **state** is a 1-based state index for the object, or 0 to append.


  * **finish** is a flag (0 or 1) which can be set to zero to improve [performance](/index.php/Category:Performance "Category:Performance") when loading large numbers of objects, but you must call **finish_object** when you are done.


  * **discrete** is a flag (0 or 1) which tells PyMOL that there will be no overlapping atoms in the file being loaded. **discrete** objects save memory but can not be edited.



Retrieved from "[https://pymolwiki.org/index.php?title=Read_Molstr&oldid=7602](https://pymolwiki.org/index.php?title=Read_Molstr&oldid=7602)"


---

## Read Pdbstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**read_pdbstr** in an API-only function which reads a pdb file from a Python string. This feature can be used to load or update structures into PyMOL without involving any temporary files. 

### PYMOL API ONLY
    
    
    cmd.read_pdbstr( string pdb-content, string object name 
       [ ,int state [ ,int finish [ ,int discrete ] ] ] )
    

### NOTES

**state** is a 1-based state index for the object. 

**finish** is a flag (0 or 1) which can be set to zero to improve [performance](/index.php/Category:Performance "Category:Performance") when loading large numbers of objects, but you must call "finish_object" when you are done. 

**discrete** is a flag (0 or 1) which tells PyMOL that there will be no overlapping atoms in the PDB files being loaded. **discrete** objects save memory but can not be edited. 

Retrieved from "[https://pymolwiki.org/index.php?title=Read_Pdbstr&oldid=7601](https://pymolwiki.org/index.php?title=Read_Pdbstr&oldid=7601)"


---

## Rebuild

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rebuild** forces PyMOL to recreate geometric objects in case any of them have gone out of sync. 

### USAGE
    
    
    rebuild [selection [, representation ]]
    

### PYMOL API
    
    
    cmd.rebuild(string selection = 'all', string representation = 'everything')
    

### SEE ALSO

[refresh](/index.php/Refresh "Refresh")

Retrieved from "[https://pymolwiki.org/index.php?title=Rebuild&oldid=7600](https://pymolwiki.org/index.php?title=Rebuild&oldid=7600)"


---

## Recolor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

See [Rebuild](/index.php/Rebuild "Rebuild"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Recolor&oldid=6782](https://pymolwiki.org/index.php?title=Recolor&oldid=6782)"


---

## Redo

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_(Incentive PyMOL only feature)_

See [undo](/index.php/Undo "Undo") for guidance on PyMOL's undo feature. 

## Contents

  * 1 Usage
  * 2 PyMOL API
  * 3 Example
  * 4 SEE ALSO



## Usage
    
    
    redo [, steps]
    

  * **steps** = integer: number of steps to undo



## PyMOL API
    
    
    cmd.redo(int steps=1)
    

## Example
    
    
     redo
     redo 5
    

## SEE ALSO

[Undo](/index.php/Undo "Undo")

Retrieved from "[https://pymolwiki.org/index.php?title=Redo&oldid=13417](https://pymolwiki.org/index.php?title=Redo&oldid=13417)"


---

## Refresh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**refresh** causes the scene to be refresh as soon as it is safe to do so. 

### USAGE

refresh 

### PYMOL API
    
    
    cmd.refresh()
    

### SEE ALSO

[rebuild ](/index.php/Rebuild "Rebuild")

Retrieved from "[https://pymolwiki.org/index.php?title=Refresh&oldid=7590](https://pymolwiki.org/index.php?title=Refresh&oldid=7590)"


---

## Reinitialize

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**reinitialize** reinitializes PyMOL. Issuing the command 'reinitialize' during a PyMOL session clears all objects and resets all parameters to defaults; this is equivalent to re-starting the program without having to actually do so. 

  * NOTE: any unsaved work will be lost!
  * NOTE: Multi-CPU machines will sometimes lose SMP abilities in PyMol after running this command.
  * NOTE: settings that usually load with a [.pymolrc](/index.php/Pymolrc "Pymolrc") file will also not work



Alternatively, this command can reset all settings to default values. 

### USAGE
    
    
    reinitialize [what [, object]]
    

### ARGUMENTS

  * what = string: everything|settings {default: everything}


  * object = string: object name for per object settings



### SEE ALSO

[Delete all](/index.php/Delete "Delete")

Retrieved from "[https://pymolwiki.org/index.php?title=Reinitialize&oldid=10600](https://pymolwiki.org/index.php?title=Reinitialize&oldid=10600)"


---

## Remove

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**remove** eleminates a selection of atoms from models. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE

remove (selection) 

### PYMOL API
    
    
    cmd.remove( string selection )
    

### EXAMPLES

remove ( resi 124 ) 

### SEE ALSO

[Delete](/index.php/Delete "Delete")

Retrieved from "[https://pymolwiki.org/index.php?title=Remove&oldid=7592](https://pymolwiki.org/index.php?title=Remove&oldid=7592)"


---

## Remove Picked

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**remove_picked** removes the atom or bond currently picked for editing. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    remove_picked [hydrogens]
    

### PYMOL API
    
    
    cmd.remove_picked(integer hydrogens=1)
    

### NOTES

  * This function is usually connected to the DELETE key and "CTRL-D".
  * By default, attached hydrogens will also be deleted unless hydrogen-flag is zero.



### SEE ALSO

[Attach](/index.php/Attach "Attach"), [Replace](/index.php/Replace "Replace")

Retrieved from "[https://pymolwiki.org/index.php?title=Remove_Picked&oldid=7593](https://pymolwiki.org/index.php?title=Remove_Picked&oldid=7593)"


---

## Rename

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rename** creates new atom names which are unique within residues. 

## Contents

  * 1 USAGE
    * 1.1 CURRENT
    * 1.2 PROPOSED
  * 2 PYMOL API
    * 2.1 CURRENT
    * 2.2 PROPOSED
  * 3 NOTES
  * 4 SEE ALSO



### USAGE

#### CURRENT
    
    
    rename object-name [ ,force ]
    force = 0 or 1 (default: 0)
    

#### PROPOSED
    
    
    rename object-or-selection,force
    

### PYMOL API

#### CURRENT
    
    
    cmd.rename( string object-name, int force )
    

#### PROPOSED
    
    
    cmd.rename( string object-or-selection, int force )
    

### NOTES

To regerate only some atom names in a molecule, first clear them with an "alter (sele),name=_" commmand, then use "rename"_

### SEE ALSO

[alter](/index.php/Alter "Alter")

Or for rename an object [set_name](/index.php/Set_name "Set name")

Retrieved from "[https://pymolwiki.org/index.php?title=Rename&oldid=10340](https://pymolwiki.org/index.php?title=Rename&oldid=10340)"


---

## Replace

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**replace** replaces the picked atom with a new atom. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    replace element, geometry, valence [,h_fill [,name ]]
    

### PYMOL API
    
    
    cmd.replace(string element, int geometry, int valence,
                int h_fill = 1, string name = "" )
    

### NOTES

Immature functionality. See code for details. 

### SEE ALSO

[Remove](/index.php/Remove "Remove"), [Attach](/index.php/Attach "Attach"), [Fuse](/index.php/Fuse "Fuse"), [Bond](/index.php/Bond "Bond"), [Unbond](/index.php/Unbond "Unbond")

Retrieved from "[https://pymolwiki.org/index.php?title=Replace&oldid=7595](https://pymolwiki.org/index.php?title=Replace&oldid=7595)"


---

## Reset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**reset** restores the rotation matrix to identity, sets the origin to the center of mass (approx.) and zooms the window and clipping planes to cover all objects. 

### USAGE
    
    
    reset
    

### PYMOL API
    
    
    cmd.reset()
    

## SEE ALSO

[Clip](/index.php/Clip "Clip") [Zoom](/index.php/Zoom "Zoom")

Retrieved from "[https://pymolwiki.org/index.php?title=Reset&oldid=8323](https://pymolwiki.org/index.php?title=Reset&oldid=8323)"


---

## Rewind

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rewind** goes to the beginning of the movie. 

### USAGE
    
    
    rewind
    

### PYMOL API
    
    
    cmd.rewind()
    

# See Also
    
    
    [Forward](/index.php/Forward "Forward"), [Backward](/index.php/Backward "Backward"), [Category:Movies](/index.php/Category:Movies "Category:Movies")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Rewind&oldid=7597](https://pymolwiki.org/index.php?title=Rewind&oldid=7597)"


---

## Rigimol.morph

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The [morph](/index.php/Morph "Morph") command is an incentive PyMOL feature (not available in open-source version). Given a two-state object, **morph** can create an interpolated trajectory from the first to the second conformation. 

_Notice: There is a new[morph](/index.php/Morph "Morph") command in PyMOL 1.6_

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Troubleshooting
  * 5 PyMOL Command
  * 6 See Also



## Usage
    
    
    morph source, target [, first [, last [, refinement [, async [, steps ]]]]]
    

## Arguments

  * source = string: name of two-state input object
  * target = string: name of output object to create
  * first = int: start state of source {default: 1}
  * last = int: end state of source {default: 2}
  * refinement = int: number of sculpting refinement cycles to clean distorted intermediates {default: 10}
  * steps = int: number of states for target object {default: 30}



_Warning: arguments first, last and steps new in PyMOL 1.5 (up to 1.4 they are always at default values)_

## Example
    
    
    # get open and closed conformation of adenylate kinase
    fetch 1ake 4ake, async=0
    
    # make two state object
    align 1ake and chain A, 4ake and chain A, cycles=0, object=aln
    create rin, 1ake and aln, 1, 1
    create rin, 4ake and aln, 1, 2
    
    # morph
    from epymol import rigimol
    rigimol.morph("rin", "rout")
    

## Troubleshooting

Rigimol likes to fail if the atom identifiers (like chain) in the two input states do not match properly. A more robust way to create the two state object is to [update](/index.php/Update "Update") state two coordinates without matching identifiers: 
    
    
    # make two state object
    align 1ake and chain A, 4ake and chain A, cycles=0, object=aln
    create rin, 1ake and aln, 1, 1
    create rin, rin, 1, 2
    update rin, 4ake and aln, 2, 1, matchmaker=0
    

## PyMOL Command

if you prefer PyMOL syntax over python syntax, add this to your [pymolrc](/index.php/Pymolrc "Pymolrc") file: 
    
    
    import epymol.rigimol
    cmd.extend('morph', epymol.rigimol.morph)
    

## See Also

  * [morph](/index.php/Morph "Morph"), the PyMOL 1.6+ morph command
  * [morpheasy](/index.php/Morpheasy "Morpheasy"), a robust wrapper for the morph command which takes care of matching atom identifiers (requires the **psico** module)
  * [slerpy](/index.php/Slerpy "Slerpy")
  * [eMovie](/index.php/EMovie "EMovie")



Retrieved from "[https://pymolwiki.org/index.php?title=Rigimol.morph&oldid=11544](https://pymolwiki.org/index.php?title=Rigimol.morph&oldid=11544)"


---

## Rms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Rms computes a RMS fit between two atom selections, but does not tranform the models after performing the fit. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 USER COMMENTS
  * 4 SEE ALSO



### USAGE
    
    
    rms (selection), (target-selection)
    

### EXAMPLES
    
    
    fit ( mutant and name ca ), ( wildtype and name ca )
    

### USER COMMENTS

To determine the RMS without any fitting, see [Rms_Cur](/index.php/Rms_Cur "Rms Cur")

[Fit](/index.php/Fit "Fit"), Rms, [Rms_Cur](/index.php/Rms_Cur "Rms Cur") are finicky and only work when all atom identifiers match: segi, chain, resn, name, alt. If they don't then you'll need to use the alter command to change the identifiers to that they do -- typically that means clearing out the SEGI field, renaming chains, and sometimes renumbering. 

I tried made two selections A, and D as 
    
    
    PyMOL>sel A, 1gh2 and n. CA and i. 65-99
    Selector: selection "A" defined with 35 atoms.
    PyMOL>sel D, 1kao and n. CA and i. 64-98
    Selector: selection "D" defined with 35 atoms
    

which as you can see both yield 35 atoms. Now, 
    
    
    rms_cur A, D
    

won't work, due to the aforementioned reason. To fix this, one needs to do, 
    
    
    alter all,segi=""
    alter all,chain=""
    alter D, resi=str(int(resi)+1)  # I don't actually use this line
    

and now 
    
    
    rms_cur A, D
    

should work. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

[Warren DeLano's comment on rms_* and commands.](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg00001.html)

Retrieved from "[https://pymolwiki.org/index.php?title=Rms&oldid=11473](https://pymolwiki.org/index.php?title=Rms&oldid=11473)"


---

## Rms cur

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

rms_cur computes the RMS difference between two atom selections without performing any fitting. 

By default, only matching atoms in both selections will be used for the fit (same chain, residue number, atoms names etc.). Alternate location mess up the match! 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
    * 3.1 Example 1: Identical Identifiers
    * 3.2 Example 2: Homologues
  * 4 See Also



## Usage
    
    
    rms_cur mobile, target [, mobile_state [, target_state [, quiet [, matchmaker [, cutoff [, cycles [, object ]]]]]]]
    

## Arguments

See [fit](/index.php/Fit "Fit"). 

## Examples

### Example 1: Identical Identifiers

Alternate location mess up the match, so remove them first. 
    
    
    fetch 1p36 1kw7, async=0
    remove not alt ""+"A"
    rms_cur 1p36, 1kw7
    

### Example 2: Homologues

Use [align](/index.php/Align "Align") or [super](/index.php/Super "Super") to create an alignment object (without fitting) and then use the alignment object in the atom selection and turn off identifier matching with **matchmaker=-1**. 
    
    
    fetch 1oky 1t46, async=0
    
    # create alignment object
    align 1oky, 1t46, cycles=0, transform=0, object=aln
    
    # RMSD of alignment object
    rms_cur 1oky & aln, 1t46 & aln, matchmaker=-1
    

## See Also

  * [align](/index.php/Align "Align")
  * [fit](/index.php/Fit "Fit")
  * [rms](/index.php/Rms "Rms")
  * [intra_fit](/index.php/Intra_fit "Intra fit")
  * [intra_rms](/index.php/Intra_rms "Intra rms")
  * [intra_rms_cur](/index.php/Intra_rms_cur "Intra rms cur")
  * [pair_fit](/index.php/Pair_fit "Pair fit")



Retrieved from "[https://pymolwiki.org/index.php?title=Rms_cur&oldid=12632](https://pymolwiki.org/index.php?title=Rms_cur&oldid=12632)"


---

## Rock

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rock** toggles Y axis rocking. 

### USAGE
    
    
    rock
    

### PYMOL API
    
    
    cmd.rock()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Rock&oldid=7587](https://pymolwiki.org/index.php?title=Rock&oldid=7587)"


---

## Rotate

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rotate** can be used to rotate the atomic coordinates of a molecular object. Behavior differs depending on whether or not the **object** parameter is specified. 

If object is **None** , then rotate rotates the atomic coordinates according to the axes and angle for the selection and state provided. 

If object is set to an object name, then selection and state are ignored and instead of translating the atomic coordinates, the [object matrix](/index.php/Object_Matrix "Object Matrix") is modified. This option is intended for use in animations. 

The "camera" option controls whether the camera or the model's axes are used to interpret the translation vector. To move the object relative to the camera set **camera=1** (which is default), or to move the molecule relative to the model's geometry, set **camera=0**. 

If this doesn't do what you want, consider [Turn](/index.php/Turn "Turn"). 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 PYMOL API
  * 4 EXAMPLES
    * 4.1 Electrostatic Map Caveat
  * 5 SEE ALSO



### USAGE
    
    
    rotate axis, angle [,selection [,state [,camera [,object [,origin]]]]]
    

### ARGUMENTS

**axis**

    

    x, y, z, or float vector: axis about which to rotate

**angle**

    

    float: degrees of rotation

**selection**

    

    string: atoms whose coordinates should be modified {default: all}

**state**

    

    > 0: only the indicated state is modified
    = 0: all states are modified
    = -1: only the current state is modified {default}

**camera**

    

    = 0 or 1: is the axis specific in camera coordinates? {default: 1}

**object**

    

    string: object name (only if changing object matrix) {default: None}

**origin**

    

    float vector: origin of rotation {default: None}

### PYMOL API
    
    
    cmd.rotate(list-or-string axis, angle=0, string selection = "all", int state = 0,
               int camera = 1, string object = None)
    

### EXAMPLES
    
    
    rotate x, 45, pept
    rotate [1,1,1], 10, chain A
    

##### Electrostatic Map Caveat

If you have an electrostatic map and it's not rotating with the molecule as you expect it to, see the [Turn](/index.php/Cmd_turn "Cmd turn") command. [Turn](/index.php/Cmd_turn "Cmd turn") moves the camera and thus the protein and map will be changed. 

### SEE ALSO

[object matrix](/index.php/Object_Matrix "Object Matrix") [Translate](/index.php/Translate "Translate") [Turn](/index.php/Turn "Turn") [Model_Space_and_Camera_Space](/index.php/Model_Space_and_Camera_Space "Model Space and Camera Space")

Retrieved from "[https://pymolwiki.org/index.php?title=Rotate&oldid=8973](https://pymolwiki.org/index.php?title=Rotate&oldid=8973)"


---

## Run

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**run** executes an external Python script in a local name space, the global namespace, or in its own namespace (as a module). 

### USAGE
    
    
    run python-script [, (local | global | module | main | private ) ]
    

### PYMOL API

Not directly available. Instead, use : 
    
    
    cmd.do("run ...").
    

### NOTES

The default mode for run is **global**. 

Due to an idiosyncrasy in Pickle, you can not pickle objects directly created at the main level in a script run as "module", (because the pickled object becomes dependent on that module). Workaround: delegate construction to an imported module. 

Retrieved from "[https://pymolwiki.org/index.php?title=Run&oldid=8991](https://pymolwiki.org/index.php?title=Run&oldid=8991)"


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

## Select

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**select** creates a named selection from an atom selection. Selections are one of the most powerful aspects of PyMOL and learning to use selections well is paramount to quickly achieving your goals in PyMOL. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 Notes
  * 5 See Also



## Usage
    
    
    select name [, selection [, enable [, quiet [, merge [, state ]]]]]
    

Convenience shortcut to create selection with name "sele": 
    
    
    select (selection)
    

## Arguments

  * **name** = str: a unique name for the selection {default if skipped: sele}
  * **selection** = str: a selection-expression
  * **enable** = 0/1: show selection indicators {default: 1}
  * **quiet** = 0/1 {default: 0 for command, 1 for API}
  * **merge** = 0/1: update existing named selection {default: 0}
  * **state** = int: object state, affects spacial operators like **within** {default: 0 (all states)}



## Examples
    
    
    select near , (ll expand 8)
    select near , (ll expand 8)
    select bb, (name CA+N+C+O)
    
    
    
    cmd.select("%s_%s"%(prefix,stretch), "none")
    

## Notes

Type **help selections** for more information about selections. 

## See Also

  * [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")
  * [Property Selectors](/index.php/Property_Selectors "Property Selectors")



Retrieved from "[https://pymolwiki.org/index.php?title=Select&oldid=12779](https://pymolwiki.org/index.php?title=Select&oldid=12779)"


---

## Set

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

set is one of the most utilized commands. PyMOL representations, states, options, etc. are changed with **set**. Briefly, set changes one of the PyMOL state variables. Currently there are over _600_ [PyMOL settings](/index.php/Category:Settings "Category:Settings")! 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 NOTES
  * 5 SEE ALSO



# USAGE
    
    
    # set '''name''' to '''value'''
    set name, [,value [,object-or-selection [,state ]]]
    
    # alternative way to do the above.
    set name = value  # (DEPRECATED)
    

# PYMOL API
    
    
    cmd.set ( string name, 
        string value=1,
        string selection='',
        int state=0,
        int updates=1,
        quiet=1)
    

# EXAMPLES
    
    
    set surface_color, red
    
    set ray_trace_mode, 3
    
    set ribbon_width, 4
    
    # set the label size to 2Ang.
    set label_size, -2
    

# NOTES

The default behavior (with a blank selection) changes the global settings database. If the selection is 'all', then the settings database in all individual objects will be changed. Likewise, for a given object, if state is zero, then the object database will be modified. Otherwise, the settings database for the indicated state within the object will be modified. 

If a selection is provided, then all objects in the selection will be affected. 

# SEE ALSO

[Get](/index.php/Get "Get")

Retrieved from "[https://pymolwiki.org/index.php?title=Set&oldid=7582](https://pymolwiki.org/index.php?title=Set&oldid=7582)"


---

## Set bond

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Set_bond sets properties on _bonds_. This is usually some atom-connecting property like [Line_width](/index.php/Line_width "Line width"), [line_color](/index.php/Line_color "Line color"), [Stick_radius](/index.php/Stick_radius "Stick radius"). 

[![](/images/b/ba/Lw2.png)](/index.php/File:Lw2.png)

[](/index.php/File:Lw2.png "Enlarge")

Line_width was set to 10 for residues 1-10

# Syntax
    
    
    # set settingName to value for object or selection objSel
    set_bond settingName, value, objSel
    

# Examples
    
    
    # Example
    load $TUT/1hpv.pdb
    remove het
    as lines
    color blue
    set_bond line_with, 5, i. 1-10
    

# See Also

  * [Set](/index.php/Set "Set")
  * [Line_width](/index.php/Line_width "Line width")
  * [Stick_radius](/index.php/Stick_radius "Stick radius")



Retrieved from "[https://pymolwiki.org/index.php?title=Set_bond&oldid=7570](https://pymolwiki.org/index.php?title=Set_bond&oldid=7570)"


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

## Set Colour

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

## Set Dihedral

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_dihedral** sets a given dihedral angle given the four atoms and one angle. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    set_dihedral atom1, atom2, atom3, atom4, angle [,state=1] [,quiet=1]
    

### PYMOL API
    
    
    set_dihedral(string atom1,string atom2,string atom3,string atom4,float angle,state=1,quiet=1):
    

### EXAMPLES
    
    
    set_dihedral resi 40 and name N, resi 40 and name CA, resi 40 and  name CB, resi 40 and name CG, -180
    

### SEE ALSO

  * [Get_Dihedral](/index.php/Get_Dihedral "Get Dihedral")
  * [DynoPlot](/index.php/DynoPlot "DynoPlot")
  * [Displaying_Biochemical_Properties#Calculating_dihedral_angles](/index.php/Displaying_Biochemical_Properties#Calculating_dihedral_angles "Displaying Biochemical Properties")
  * [Rotamer_Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle")



Retrieved from "[https://pymolwiki.org/index.php?title=Set_Dihedral&oldid=8184](https://pymolwiki.org/index.php?title=Set_Dihedral&oldid=8184)"


---

## Set Geometry

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_geometry** changes PyMOL's assumptions about the proper valence and geometry of the picked atom. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    set_geometry geometry, valence
    

### PYMOL API
    
    
    cmd.set_geometry(int geometry,int valence )
    

### NOTES

Immature functionality. See code for details. 

### SEE ALSO

[Remove](/index.php/Remove "Remove"), [Attach](/index.php/Attach "Attach"), [Fuse](/index.php/Fuse "Fuse"), [Bond](/index.php/Bond "Bond"), [Unbond](/index.php/Unbond "Unbond")

Retrieved from "[https://pymolwiki.org/index.php?title=Set_Geometry&oldid=7579](https://pymolwiki.org/index.php?title=Set_Geometry&oldid=7579)"


---

## Set Key

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_key** binds a Python function or PyMOL command to a key press. Such key binding customization is typically done with your [pymolrc](/index.php/Pymolrc "Pymolrc") startup script. 

_Changes with PyMOL version:_

  * 1.8: decorator support
  * 1.7: second argument can also be a string in PyMOL command syntax



## Contents

  * 1 PyMOL API
  * 2 PyMOL Command (since PyMOL 1.7)
  * 3 Examples
  * 4 KEYS WHICH CAN BE REDEFINED
  * 5 See Also



## PyMOL API
    
    
    cmd.set_key(string key, function fn, tuple arg=(), dict kw={})
    

_Also supported since PyMOL 1.7:_
    
    
    cmd.set_key(string key, string command)
    

_Decorator support since PyMOL 1.8:_
    
    
    cmd.set_key(string key)(function fn)
    

## PyMOL Command (since PyMOL 1.7)
    
    
    set_key key, command
    

## Examples
    
    
    from pymol import cmd
    
    # define a custom function which colors a selection blue
    def make_it_blue(selection): cmd.color("blue", selection)
    
    # color "object1" blue when the F1 key is pressed
    cmd.set_key( 'F1' , make_it_blue, [ "object1" ] )
    
    
    
    # zoom on everything
    cmd.set_key( 'CTRL-C' , cmd.zoom )
    
    
    
    # turn camera by 90 degrees
    cmd.set_key( 'ALT-A' , cmd.turn, ('x',90) )
    
    
    
    # zoom on first ligand when pressing F1 (PyMOL 1.7+)
    set_key F1, zoom byres (first organic), animate=1
    
    
    
    # zoom on first ligand when pressing F1 (all PyMOL versions)
    cmd.set_key('F1', cmd.zoom, ['byres (first organic)'], {'animate': 1})
    
    
    
    # decorator syntax (PyMOL 1.8+)
    @cmd.set_key('F1')
    def zoom_ligand():
        cmd.zoom('byres (first organic)', animate=1)
    

## KEYS WHICH CAN BE REDEFINED
    
    
    F1 to F12
    left, right, pgup, pgdn, home, insert
    CTRL-A to CTRL-Z 
    ALT-0 to ALT-9, ALT-A to ALT-Z
    

## See Also

  * [alias](/index.php/Alias "Alias")
  * [cmd.extend](/index.php/Extend "Extend")
  * [Button](/index.php/Button "Button")
  * [Check Key](/index.php/Check_Key "Check Key")



Retrieved from "[https://pymolwiki.org/index.php?title=Set_Key&oldid=12315](https://pymolwiki.org/index.php?title=Set_Key&oldid=12315)"


---

## Set name

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_name** can be used to change the name of an object or selection. 

Not only can you simply rename an object or selection, but this command is also a powerful tool for those who deal with multiple structures in one file --- say a collection of NMR models. The user can execute the [Split_States](/index.php/Split_States "Split States") command and then rename the molecule of choice in the state of choice. For example, if one loads an NMR structure (with, say, 20 states) and aligns it to another structure, the balance of the alignment will (most likely) be off due to the weighting of the 19 other structures you probably don't see. To overcome this problem, one simply executes [Split_States](/index.php/Split_States "Split States") and then renames one of the states and then aligns that newly renamed object. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 User Comments/Examples
    * 3.1 SEE ALSO



### USAGE
    
    
    set_name old_name, new_name
    

### PYMOL API
    
    
    cmd.set_name(string old_name, string new_name)
    

  


## User Comments/Examples
    
    
    cmd.set_name("example", "nicename")
    

As the order of arguments is different from, for example, the [Select](/index.php/Select "Select") command (where the desired name for the selection is the first argument), it may be helpful to think of this as similar to the 'mv' shell command, where the existing object is listed first and the destination second. 

  
To batch-rename several objects with an appended number at once, use: 
    
    
    for a in range(1,999):cmd.set_name("example"+str(a), "nicename"+str(a))
    

### SEE ALSO

[get_names](/index.php/Get_names "Get names"), [get_legal_name](/index.php/Get_legal_name "Get legal name"), [get_unused_name](/index.php/Get_unused_name "Get unused name")

Retrieved from "[https://pymolwiki.org/index.php?title=Set_name&oldid=11439](https://pymolwiki.org/index.php?title=Set_name&oldid=11439)"


---

## Set Symmetry

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_symmetry** can be used to define or redefine the crystal and spacegroup parameters for a molecule or map object. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 Example
  * 4 NOTES



### USAGE
    
    
    set_symmetry selection, a, b, c, alpha, beta, gamma, spacegroup
    

### PYMOL API
    
    
    cmd.set_symmetry(string selection, float a, float b, float c,
         float alpha,float beta, float gamma, string spacegroup)
    

### Example
    
    
    # PyMOL command line
    set_symmetry 1a2p, 60, 60, 80, 90, 90, 120, P6122
    
    # API
    cmd.set_symmetry("1a2p", 60, 60, 80, 90, 90, 120, spacegroup="P6122")
    

### NOTES

The new symmetry will be defined for every object referenced by the selection. 

Retrieved from "[https://pymolwiki.org/index.php?title=Set_Symmetry&oldid=12899](https://pymolwiki.org/index.php?title=Set_Symmetry&oldid=12899)"


---

## Set Title

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_title** attaches a text string to the state of a particular object which can be displayed when the state is active. This is useful for display the energies of a set of conformers. 

### USAGE
    
    
    set_title object,state,text
    

### PYMOL API
    
    
    cmd.set_title(string object,int state,string text)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Set_Title&oldid=7575](https://pymolwiki.org/index.php?title=Set_Title&oldid=7575)"


---

## Set title

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Sets the title of an object. 

# Usage
    
    
    set_title objName, stateNo, newTitle
    
    # for example,
    set_title ligs3d, 1, This is a new title.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Set_title&oldid=7569](https://pymolwiki.org/index.php?title=Set_title&oldid=7569)"


---

## Set View

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_view** sets viewing information for the current scene, including the rotation matrix, position, origin of rotation, clipping planes, and the orthoscopic flag. 

This command is extremely useful for making movies. One may set up the scene to be rendered, then save the exact orientation, with respect to the camera, of the scene using, the [Get_View](/index.php/Get_View "Get View") command. The output from the [Get_View](/index.php/Get_View "Get View") command may then be used by the **set_view** command to restore the orientation of the scene. 

## Contents

  * 1 USAGE
  * 2 EXAMPLE
  * 3 PYMOL API
    * 3.1 Example
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    set_view (...)  where ... is 18 floating point numbers
    

### EXAMPLE

This works in a `.pml` script: 
    
    
       set_view (\
           0.999876618,   -0.000452542,   -0.015699286,\
           0.000446742,    0.999999821,   -0.000372844,\
           0.015699454,    0.000365782,    0.999876678,\
           0.000000000,    0.000000000, -150.258514404,\
           11.842411041,   20.648729324,    8.775371552,\
           118.464958191,  182.052062988,    0.000000000 )
    

### PYMOL API
    
    
    cmd.set_view(string-or-sequence view)
    

#### Example

This works in a `.py` script: 
    
    
    cmd.set_view((
        0.999876618,   -0.000452542,   -0.015699286,
        0.000446742,    0.999999821,   -0.000372844,
        0.015699454,    0.000365782,    0.999876678,
        0.000000000,    0.000000000, -150.258514404,
        11.842411041,   20.648729324,    8.775371552,
        118.464958191,  182.052062988,    0.000000000))
    

The result of [get_view](/index.php/Get_view "Get view") is valid input for **set_view** : 
    
    
    myview = cmd.get_view()
    
    cmd.set_view(myview)
    

### NOTES

Contents of the view matrix 

  * 0 - 8: column-major 3x3 matrix which rotates model axes to camera axes
  * 9 - 11: origin of rotation relative to the camera in camera space
  * 12 - 14: origin of rotation in model space
  * 15: front plane distance from the camera
  * 16: rear plane distance from the camera
  * 17: orthoscopic flag (not implemented in older versions)



### SEE ALSO

[Get View](/index.php/Get_View "Get View"), [View](/index.php/View "View")

Retrieved from "[https://pymolwiki.org/index.php?title=Set_View&oldid=13269](https://pymolwiki.org/index.php?title=Set_View&oldid=13269)"


---

## Show

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Show** displays atom and bond representations for certain selections. 

## Contents

  * 1 Details
    * 1.1 USAGE
    * 1.2 PYMOL API
    * 1.3 EXAMPLES
      * 1.3.1 Example
      * 1.3.2 Example
      * 1.3.3 Example
      * 1.3.4 Example
    * 1.4 NOTES
    * 1.5 SEE ALSO



# Details

  * [![Example of a shown surface.](/images/e/e7/Ray2.png)](/index.php/File:Ray2.png "Example of a shown surface.")

Example of a shown surface. 

  * [![Ball and sticks shown.](/images/7/76/Ray_trace_gain2.png)](/index.php/File:Ray_trace_gain2.png "Ball and sticks shown.")

Ball and sticks shown. 

  * [![A cartoon inside a mesh shown.](/images/6/6b/Mesh_w05.png)](/index.php/File:Mesh_w05.png "A cartoon inside a mesh shown.")

A cartoon inside a mesh shown. 




The **Show** command, is one of the most often used commands in PyMOL. For example, you can _show_ certain atoms as _[Lines](/index.php/Lines "Lines")_ or _[Sticks](/index.php/Sticks "Sticks")_ or _[Cartoons](/index.php/Cartoon "Cartoon")_ or any of the following representations: 

  * [lines](/index.php/Lines "Lines")
  * [spheres](/index.php/Spheres "Spheres")
  * [mesh](/index.php/Mesh "Mesh")
  * [ribbon](/index.php/Ribbon "Ribbon")
  * [cartoon](/index.php/Cartoon "Cartoon")
  * [sticks](/index.php/Sticks "Sticks")
  * [dots](/index.php/Dots "Dots")
  * [surface](/index.php/Surface "Surface")
  * [label](/index.php/Label "Label")
  * [extent](/index.php?title=Extent&action=edit&redlink=1 "Extent \(page does not exist\)")
  * [nonbonded](/index.php?title=Nonbonded&action=edit&redlink=1 "Nonbonded \(page does not exist\)")
  * [nb_spheres](/index.php/Nb_spheres "Nb spheres")
  * [slice](/index.php/Slice "Slice")
  * [cell](/index.php/Cell "Cell")



## USAGE
    
    
    show
    show reprentation [,object]
    show reprentation [,(selection)]
    show (selection)
    

## PYMOL API
    
    
    cmd.show( string representation="", string selection="" )
    

## EXAMPLES

#### Example
    
    
    # show the backbone using lines.
    show lines,(name ca or name c or name n)
    

#### Example
    
    
    # show the ribbon representation for all objects
    show ribbon
    

#### Example
    
    
    # show all hetero atoms as spheres
    show spheres, het
    

#### Example
    
    
    # show only polar hydrogens
    hide everything, ele h
    show lines, ele h and neighbor (ele n+o)
    # hide nonpolar hydrogens
    hide (h. and (e. c extend 1))
    

Note: 

The above code hides all representations of nonpolar hydrogens, including surface representations, resulting in broken surface representations. It might be better to remove the nonpolar hydrogens instead: 
    
    
    # show only polar hydrogens
    hide everything, ele h
    show lines, ele h and neighbor (ele n+o)
    # remove nonpolar hydrogens
    remove (h. and (e. c extend 1))
    

## NOTES

**selection** can be an object name. 

**show** alone will turn on lines and nonbonded for all bonds. 

**show cell** will show the triclinic unit cell, provided it's described in the PDB. For cell packing, you need a script. See [The Script Library](/index.php/Category:Script_Library "Category:Script Library") or [Robert Campbell's Site](http://adelie.biochem.queensu.ca/~rlc/work/pymol/)

## SEE ALSO

[Hide](/index.php/Hide "Hide"), [Enable](/index.php/Enable "Enable"), [Disable](/index.php/Disable "Disable"), [Suspend_updates](/index.php/Suspend_updates "Suspend updates"), [Lines](/index.php/Lines "Lines"), [Cartoon](/index.php/Cartoon "Cartoon"), [Spheres](/index.php/Spheres "Spheres"), [Mesh](/index.php/Mesh "Mesh"), [Ribbon](/index.php/Ribbon "Ribbon"), [Sticks](/index.php/Sticks "Sticks"), [Dots](/index.php/Dots "Dots"), [Surface](/index.php/Surface "Surface"), [Label](/index.php/Label "Label"), [Extent](/index.php?title=Extent&action=edit&redlink=1 "Extent \(page does not exist\)"), [Slice](/index.php/Slice "Slice"), [Cell](/index.php/Cell "Cell"), [Select](/index.php/Select "Select"), [Show_as](/index.php/Show_as "Show as")

Retrieved from "[https://pymolwiki.org/index.php?title=Show&oldid=11744](https://pymolwiki.org/index.php?title=Show&oldid=11744)"


---

## Show as

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

show_as turns on and off atom and bond [representations](/index.php/Category:Representations "Category:Representations"). 

## Contents

  * 1 Details
    * 1.1 USAGE
    * 1.2 ARGUMENTS
    * 1.3 EXAMPLES
    * 1.4 PYMOL API
    * 1.5 NOTES
    * 1.6 SEE ALSO
    * 1.7 References



# Details

The available representations are the usual: 

(History: Python 2.6 came out with the keyword [as](/index.php/As "As") which is also a PyMOL keyword. So, we had to change the PyMOL-named keyword to show_as.) 

  * [lines](/index.php/Lines "Lines")
  * [spheres](/index.php/Spheres "Spheres")
  * [mesh](/index.php/Mesh "Mesh")
  * [ribbon](/index.php/Ribbon "Ribbon")
  * [cartoon](/index.php/Cartoon "Cartoon")
  * [sticks](/index.php/Sticks "Sticks")
  * [dots](/index.php/Dots "Dots")
  * [surface](/index.php/Surface "Surface")
  * [labels](/index.php?title=Labels&action=edit&redlink=1 "Labels \(page does not exist\)")
  * [extent](/index.php?title=Extent&action=edit&redlink=1 "Extent \(page does not exist\)")
  * [nonbonded](/index.php?title=Nonbonded&action=edit&redlink=1 "Nonbonded \(page does not exist\)")
  * [nb_spheres](/index.php/Nb_spheres "Nb spheres")
  * [slice](/index.php/Slice "Slice")



## USAGE
    
    
    show_as representation [, selection ]
    

## ARGUMENTS

  * **representation** = lines, spheres, mesh, ribbon, cartoon, sticks, dots, surface, labels, extent, nonbonded, nb_spheres, slice, extent, slice, dashes, angles, dihedrals, cgo, cell, callback, everything
  * **selection** = string {default: all}



## EXAMPLES
    
    
     
    # show the backbone as lines
    show_as lines, name ca or name c or name n
    
    # show everything as a ribbon
    show_as ribbon
    

## PYMOL API
    
    
     
    cmd.show_as(string representation, string selection)
    

## NOTES

  * **selection** can be an object name
  * **as** alone will turn on lines and nonbonded and hide everything else.



## SEE ALSO

[show](/index.php/Show "Show"), [hide](/index.php/Hide "Hide"), [enable](/index.php/Enable "Enable"), [disable](/index.php/Disable "Disable"), [All PyMOL Representations](/index.php/Category:Representations "Category:Representations")

## References

  * PyMOL Source code



Retrieved from "[https://pymolwiki.org/index.php?title=Show_as&oldid=7571](https://pymolwiki.org/index.php?title=Show_as&oldid=7571)"


---

## Slice

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/5/57/Slice_example_2_screenshot.png)](/index.php/File:Slice_example_2_screenshot.png)

[](/index.php/File:Slice_example_2_screenshot.png "Enlarge")

Screenshot of slice example #2

slice creates a slice object from a [map](/index.php?title=Category:Maps&action=edit&redlink=1 "Category:Maps \(page does not exist\)") object. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES
  * 4 PYMOL API
  * 5 NOTES
  * 6 SEE ALSO
  * 7 References



## USAGE
    
    
    slice_new name, map, [state, [source_state]]
    

_or_
    
    
    slice name, map, [state, [source_state]]
    

## ARGUMENTS

  * **name** = the name for the new slice object (string)
  * **map** = the name of the map object to use for computing the slice (string)
  * **state** = the state into which the object should be loaded (default=1; set state=0 to append new mesh as a new state)
  * **source_state** = the state of the map from which the object should be loaded (default=0)



## EXAMPLES
    
    
     
    # Create a map slice plane perpendicular to current view
    slice a_new_slice, a_map
    
    
    # A more complicated example that shows how to create multiple slices
    # (in this case, 3 slices perpendicular to each other), each colored
    # with a different color ramp and different contour levels:
    
    # Reset the view, to align view on XYZ axes
    reset
    # (Optional: Adjust view direction to your liking)
    # Create a map slice *perpendicular* to the current view.
    # The slice seems to be in the center of the APBS (or other) volmap. Map "tracking" is off by default.
    slice slice_A, apbs_map
    # Rotate camera 90 degrees about the vertical axis
    turn y, 90
    # Second, perpendicular, slice
    slice slice_B, apbs_map
    # Rotate again, this time about the horizontal
    turn x, 90
    # Third slice
    slice slice_C, apbs_map
    
    # Define new color ramps: ramp_name, map_object, list of low/mid/hi values, 3 RGB triplets
    ramp_new ramp1010RWB, apbs_map, [-10,0,10], [ [1,0,0], [1,1,1], [0,0,1] ]
    ramp_new ramp11RYG, apbs_map, [-1,0,1], [ [1,0,0], [1,1,0], [0,1,0] ]
    ramp_new ramp55MltGO, apbs_map, [-5,0,5], [ [0,1,1], [0.5,1,0.5], [1,0.5,0.2] ]
    
    # Color the map slices
    color ramp1010RWB, slice_A
    color ramp11RYG, slice_B
    color ramp55MltGO, slice_C
    
    # Adjust the fineness of the slice color gradations:
    cmd.set('slice_grid',0.1) # normally at 0.3; much finer than 0.05 gets a bit slow
    
    # Map slices can be moved, relative to other fixed objects (e.g., your protein/DNA/RNA), 
    # by turning tracking on (Action menu), and the using the Shift-MouseWheel to move
    # the slice forward in and backward in Z. Adjust fineness of this Z-motion with:
    #   cmd.set('mouse_wheel_scale',0.05) # normally at 0.5
    
    # The result is shown in the image above.
    

## PYMOL API
    
    
     
    cmd.slice_new(string slice_name, string map_name, integer state=0, integer source_state=0)
    

## NOTES

Mis-identified as "slice_map" in documentation. slice or slice_new is the correct command. Documentation also mentions two optional parameters which seem to be no longer supported (opacity = opacity of the new slice, and resolution = the number of pixels per sampling). 

  


## SEE ALSO

[isomesh](/index.php/Isomesh "Isomesh"), [isodot](/index.php/Isodot "Isodot"), [load](/index.php/Load "Load")

## References

  * PyMOL Source code



Retrieved from "[https://pymolwiki.org/index.php?title=Slice&oldid=8954](https://pymolwiki.org/index.php?title=Slice&oldid=8954)"


---

## Slice new

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/5/57/Slice_example_2_screenshot.png)](/index.php/File:Slice_example_2_screenshot.png)

[](/index.php/File:Slice_example_2_screenshot.png "Enlarge")

Screenshot of slice example #2

slice creates a slice object from a [map](/index.php?title=Category:Maps&action=edit&redlink=1 "Category:Maps \(page does not exist\)") object. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES
  * 4 PYMOL API
  * 5 NOTES
  * 6 SEE ALSO
  * 7 References



## USAGE
    
    
    slice_new name, map, [state, [source_state]]
    

_or_
    
    
    slice name, map, [state, [source_state]]
    

## ARGUMENTS

  * **name** = the name for the new slice object (string)
  * **map** = the name of the map object to use for computing the slice (string)
  * **state** = the state into which the object should be loaded (default=1; set state=0 to append new mesh as a new state)
  * **source_state** = the state of the map from which the object should be loaded (default=0)



## EXAMPLES
    
    
     
    # Create a map slice plane perpendicular to current view
    slice a_new_slice, a_map
    
    
    # A more complicated example that shows how to create multiple slices
    # (in this case, 3 slices perpendicular to each other), each colored
    # with a different color ramp and different contour levels:
    
    # Reset the view, to align view on XYZ axes
    reset
    # (Optional: Adjust view direction to your liking)
    # Create a map slice *perpendicular* to the current view.
    # The slice seems to be in the center of the APBS (or other) volmap. Map "tracking" is off by default.
    slice slice_A, apbs_map
    # Rotate camera 90 degrees about the vertical axis
    turn y, 90
    # Second, perpendicular, slice
    slice slice_B, apbs_map
    # Rotate again, this time about the horizontal
    turn x, 90
    # Third slice
    slice slice_C, apbs_map
    
    # Define new color ramps: ramp_name, map_object, list of low/mid/hi values, 3 RGB triplets
    ramp_new ramp1010RWB, apbs_map, [-10,0,10], [ [1,0,0], [1,1,1], [0,0,1] ]
    ramp_new ramp11RYG, apbs_map, [-1,0,1], [ [1,0,0], [1,1,0], [0,1,0] ]
    ramp_new ramp55MltGO, apbs_map, [-5,0,5], [ [0,1,1], [0.5,1,0.5], [1,0.5,0.2] ]
    
    # Color the map slices
    color ramp1010RWB, slice_A
    color ramp11RYG, slice_B
    color ramp55MltGO, slice_C
    
    # Adjust the fineness of the slice color gradations:
    cmd.set('slice_grid',0.1) # normally at 0.3; much finer than 0.05 gets a bit slow
    
    # Map slices can be moved, relative to other fixed objects (e.g., your protein/DNA/RNA), 
    # by turning tracking on (Action menu), and the using the Shift-MouseWheel to move
    # the slice forward in and backward in Z. Adjust fineness of this Z-motion with:
    #   cmd.set('mouse_wheel_scale',0.05) # normally at 0.5
    
    # The result is shown in the image above.
    

## PYMOL API
    
    
     
    cmd.slice_new(string slice_name, string map_name, integer state=0, integer source_state=0)
    

## NOTES

Mis-identified as "slice_map" in documentation. slice or slice_new is the correct command. Documentation also mentions two optional parameters which seem to be no longer supported (opacity = opacity of the new slice, and resolution = the number of pixels per sampling). 

  


## SEE ALSO

[isomesh](/index.php/Isomesh "Isomesh"), [isodot](/index.php/Isodot "Isodot"), [load](/index.php/Load "Load")

## References

  * PyMOL Source code



Retrieved from "[https://pymolwiki.org/index.php?title=Slice&oldid=8954](https://pymolwiki.org/index.php?title=Slice&oldid=8954)"


---

## Smooth

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Smooth performs a window average of coordinate states. 

# Usage
    
    
    smooth [ selection [, passes [, window [, first [, last [, ends]]]]]]
    
    # for example, smooth and object with 30 passes and
    # a window size of 100.
    smooth myObj, 30, 100
    

  * ends = 0 or 1: controls whether or not the end states are also smoothed using a weighted asymmetric window



To see smooth in context see this [example](http://www.pymolwiki.org/index.php/Protect#Example)

# NOTES

  * This type of averaging is often used to suppress high-frequency vibrations in a molecular dynamics trajectory.
  * This function is not memory efficient. For reasons of flexibility, it uses two additional copies of every atomic coordinate for the calculation. If you are memory-constrained in visualizing MD trajectories, then you may want to use an external tool such as ptraj to perform smoothing before loading coordinates into PyMOL.



# See Also

[Load_Traj](/index.php/Load_Traj "Load Traj"), [protect](/index.php/Protect "Protect"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Smooth&oldid=12096](https://pymolwiki.org/index.php?title=Smooth&oldid=12096)"


---

## Sort

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**sort** reorders atoms in the structure. It usually only necessary to run this routine after an **alter** command which has modified the names of atom properties. Without an argument, sort will resort all atoms in all objects. 

For more information please see [alter](/index.php/Alter "Alter") command. 

# USAGE
    
    
    sort [object]
    

# PYMOL API
    
    
    cmd.sort(string object)
    

# SEE ALSO

[alter](/index.php/Alter "Alter")

Retrieved from "[https://pymolwiki.org/index.php?title=Sort&oldid=11419](https://pymolwiki.org/index.php?title=Sort&oldid=11419)"


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

## Splash

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Splash, shows the PyMOL start splash screen. 

Retrieved from "[https://pymolwiki.org/index.php?title=Splash&oldid=7564](https://pymolwiki.org/index.php?title=Splash&oldid=7564)"


---

## Split states

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Split_States** splits and orients multiple models and multimers from the biological unit file into a set of single-state molecular objects. 

## Contents

  * 1 Syntax
  * 2 Using
  * 3 Example
  * 4 See Also



## Syntax
    
    
    split_states object [, first [, last [, prefix ]]]
    

This splits the **object** from **first** to **last** out to the array of objects prefixed by **prefix**. The **prefix** option is very handy if all your states--or a subset of the states--have the same name. 

## Using

To use **split_states** simply Load your molecule 
    
    
    # example usage
    load fileName.pdb1, name
    split_states name
    delete name
    
    # split all the states to objects starting with conf
    fetch 1nmr
    split_states 1nmr, prefix=conf
    

## Example

**1VLS** : A dimer. 
    
    
    load 1vls.pdb1, 1vls
    split_states 1vls
    dele 1vls
    

  * [![1VLS as a monomer. This is the state of 1VLS when I load the molecule \(and select cartoon representation\).](/images/e/ef/1vls1.png)](/index.php/File:1vls1.png "1VLS as a monomer. This is the state of 1VLS when I load the molecule \(and select cartoon representation\).")

1VLS as a monomer. This is the state of 1VLS when I load the molecule (and select cartoon representation). 

  * [![1VLS as a dimer using the split_states command. Notice PyMOL automatically loads and orients the new molecules.](/images/f/fa/1vls1_dimer.png)](/index.php/File:1vls1_dimer.png "1VLS as a dimer using the split_states command. Notice PyMOL automatically loads and orients the new molecules.")

1VLS as a dimer using the split_states command. Notice PyMOL automatically loads and orients the new molecules. 




# See Also

  * [PDB Tutorial Biological Units](http://www.rcsb.org/pdb/static.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html)



Retrieved from "[https://pymolwiki.org/index.php?title=Split_states&oldid=10857](https://pymolwiki.org/index.php?title=Split_states&oldid=10857)"


---

## Super

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**super** aligns two selections. It does a **sequence-independent** (unlike [align](/index.php/Align "Align")) structure-based dynamic programming alignment followed by a series of refinement cycles intended to improve the fit by eliminating pairing with high relative variability (just like [align](/index.php/Align "Align")). **super** is more robust than **align** for proteins with low sequence similarity. 

## Contents

  * 1 Usage
  * 2 Caveats
  * 3 User Scripts
    * 3.1 Write rmsd to file
    * 3.2 Write rmsd to file and looping
  * 4 See Also



## Usage

See [align](/index.php/Align "Align") command. 

## Caveats

  * **Alternative Conformations:** If super ever tells you no matched atoms, then instead of 
        
        super p1, p2
        

try 
        
        super p1 & alt A+'', p2 & alt B+''
        




## User Scripts

### Write rmsd to file

**pymol_rmsd_test.pml**
    
    
    reinitialize
    
    fetch 1F9J, async=0
    fetch 1YX5, async=0
    
    extract 1F9J_A, 1F9J and chain A
    extract 1YX5_B, 1YX5 and chain B
    
    test=cmd.super("1F9J_A","1YX5_B")
    
    python
    writefile=open("rmsd_file.txt","a")
    writefile.write(' '.join('%s' % x for x in test))
    writefile.write('\n')
    writefile.close()
    python end
    

In terminal 
    
    
    pymol -c pymol_rmsd_test.pml ; tail -n 1 rmsd_file.txt
    

### Write rmsd to file and looping

**pymol_pymol_loop.sh**
    
    
    #!/bin/csh -f
    set x = 1
    
    while ( $x <= 20 )
    	set prot="energy_${x}.pdb"
    	pymol -c pymol_super.pml $prot
    	@ x = $x + 1
    end
    

  
**pymol_super.pml**
    
    
    reinitialize
    import sys
    
    python 
    prot1="XXXX"
    prot2="YYYY_trimmed"
    
    prot3=sys.argv[3]
    #prot3="energy_54.pdb"
    prot3name=prot3.split(".pdb")[0]
    print prot3, prot3name
    python end
    
    cd /XXXXX
    
    cmd.load("%s.pdb"%prot1)
    cmd.load("%s.pdb"%prot2)
    cmd.load("./ensemblesize2_numstruct/%s"%prot3)
    #show_as cartoon, all
    
    align1=cmd.super("%s"%prot3name,"%s"%prot1)
    print align1
    
    python
    writefile=open("pymol_rmsd_file.txt","a")
    writefile.write('%s %s '%(prot3name, prot1))
    writefile.write(' '.join('%s' % x for x in align1))
    writefile.write(' ')
    python end
    
    align2=cmd.super("%s"%prot3name,"%s"%prot2)
    print align2
    
    python
    writefile=open("pymol_rmsd_file.txt","a")
    writefile.write('%s %s '%(prot3name, prot2))
    writefile.write(' '.join('%s' % x for x in align2))
    writefile.write('\n')
    writefile.close()
    python end
    

In terminal 
    
    
    chmod +x pymol_loop.sh
    ./pymol_loop.sh
    

## See Also

  * [align](/index.php/Align "Align")
  * [Cealign](/index.php/Cealign "Cealign")
  * [Get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")



Retrieved from "[https://pymolwiki.org/index.php?title=Super&oldid=12180](https://pymolwiki.org/index.php?title=Super&oldid=12180)"


---

## Symexp

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Symexp is used to reconstruct neighboring asymmetric units from the crystallographic experiment that produced the given structure. This is assuming the use of a [PDB](http://www.rcsb.org/pdb/home/home.do) file or equivalent that contains enough information ([CRYST1 record](http://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1)) to reproduce the lattice. 

Symexp creates all symmetry related objects for the specified object that occurs within a cutoff about an atom selection. The new objects are labeled using the prefix provided along with their crystallographic symmetry operation and translation. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLE
  * 4 Naming Scheme
  * 5 See Also



## USAGE
    
    
    # Expand the ''object'' around its ''selection'' by cutoff Angstroms and
    # prefix the new objects withs ''prefix''.
    symexp prefix, object, selection, cutoff [, segi]
    

For one protein: 
    
    
    symexp name_for_new_objects,asymmetric_name,(asymmetric_name),distance
    

## ARGUMENTS

  * **prefix** = string: name prefix for new objects
  * **object** = string: name of the object that you wish to reproduce neighboring crystal partners for; the source of the symmetry operators
  * **selection** = string: atom selection to measure cutoff distance from
  * **cutoff** = float: create all symmetry mates that are within "cutoff" distance from selection (but no more than +/-1 unit cells)
  * **segi** = 0/1: if segi=1 then assign to each symmetry mate a unique 4-character segment identifier {default: 0}



## EXAMPLE

load any .pdb file into PyMOL (here we use 1GVF).  

    
    
    fetch 1GVF
    

[![1GVF assym.png](/images/b/b3/1GVF_assym.png)](/index.php/File:1GVF_assym.png)

At the PyMOL command prompt type the following:  

    
    
    symexp sym,1GVF,(1GVF),1
    

produces three new objects. We now have four objects corresponding to two biologic units (the functional protein in a cell). 

[![1GVF 1A.png](/images/8/83/1GVF_1A.png)](/index.php/File:1GVF_1A.png)

  

    
    
    symexp sym,1GVF,(1GVF),5
    

If we color all of the sym* cyan we will produce the following:  


[![1GVF 5A.jpeg](/images/e/ed/1GVF_5A.jpeg)](/index.php/File:1GVF_5A.jpeg)

As you can see, we can begin to understand the crystal environment of our asymmetric unit. Increasing _distance_ will reveal more of the crystal lattice, but will place in increasing demand on your computer's rendering ability. 

PyMOL is known to exit dramatically (crash) if you provide a scene that is too large or complex. This is a result of the low-level _malloc_ function failing. See [Category:Performance](/index.php/Category:Performance "Category:Performance") for workarounds. 

## Naming Scheme

The created objects have the following naming scheme: 
    
    
    <prefix>AAXXYYZZ
    

  * AA: zero-based [symmetry operator](https://sourceforge.net/p/pymol/code/HEAD/tree/trunk/pymol/modules/pymol/xray.py#l95) index, e.g. 00 always corresponds to "x,y,z". A space group like "P 21 21 21" which has 4 symmetry operators will count up to 03
  * XX: -1, 0, or 1, this is the unit cell offset from the center of the selection along the x-axis. PyMOL never creates more than 3 unit cells along each axis
  * YY: analogous to XX, but along the y-axis
  * ZZ: analogous to XX, but along the z-axis



## See Also

  * [PDB Symmetry Info](http://pdbbeta.rcsb.org/robohelp_f/data_download/biological_unit/pdb_and_mmcif_files_.htm)
  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [Supercell](/index.php/Supercell "Supercell")
  * From within PyMOL, _help symexp_ and _symexp ?_.



Retrieved from "[https://pymolwiki.org/index.php?title=Symexp&oldid=12746](https://pymolwiki.org/index.php?title=Symexp&oldid=12746)"


---

## Sync

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

sync is an API-only function which waits until all current commands have been executed before returning. A timeout can be used to insure that this command eventually returns. 

## PYMOL API
    
    
    cmd.sync(timeout: float = 1.0, poll: float = 0.05)
    

## SEE ALSO

[frame](/index.php/Frame "Frame")

Retrieved from "[https://pymolwiki.org/index.php?title=Sync&oldid=13551](https://pymolwiki.org/index.php?title=Sync&oldid=13551)"


---

## System

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **system** command, executes a command in a subshell under Unix or Windows. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 NOTES
  * 4 SEE ALSO



### USAGE
    
    
    # execute 'command'
    system command
    

### PYMOL API
    
    
    cmd.system(string command,int async=0)
    

### NOTES

async can only be specified from the Python level (not the command language) 

  * if async is 0 (default), then the result code from "system" is returned in r
  * if async is 1, then the command is run in a separate thread whose object is returned



### SEE ALSO

[ls](/index.php/Ls "Ls"), [cd](/index.php/Cd "Cd"), [pwd](/index.php/Pwd "Pwd")

Retrieved from "[https://pymolwiki.org/index.php?title=System&oldid=8245](https://pymolwiki.org/index.php?title=System&oldid=8245)"


---

## Torsion

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**torsion** rotates the torsion on the bond currently picked for editing. The rotated fragment will correspond to the first atom specified when picking the bond (or the nearest atom, if picked using the mouse). 

### USAGE
    
    
    torsion angle
    

### PYMOL API
    
    
    cmd.torsion( float angle )
    

## SEE ALSO

[edit](/index.php/Edit "Edit"), [unpick](/index.php/Unpick "Unpick"), [remove_picked](/index.php/Remove_picked "Remove picked"), [cycle_valence](/index.php/Cycle_valence "Cycle valence")

Retrieved from "[https://pymolwiki.org/index.php?title=Torsion&oldid=7558](https://pymolwiki.org/index.php?title=Torsion&oldid=7558)"


---

## Translate

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**translate** can be used to translate the atomic coordinates of a molecular object (move an object). Behavior differs depending on whether or not the "object" parameter is specified. 

If object is None, then translate translates atomic coordinates according to the vector provided for the selection and in the state provided. 

If object is set to an object name, then selection and state are ignored and instead of translating the atomic coordinates, the [object matrix](/index.php/Object_Matrix "Object Matrix") is modified. This option is intended for use in animations. 

The "camera" option controls whether the camera or the model's axes are used to interpret the translation vector. To move the object relative to the camera set **camera=1** (which is default), or to move the molecule relative to the model's geometry, set, **camera=0**. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 PYMOL API
  * 4 EXAMPLES
    * 4.1 Example 1
    * 4.2 Example 2
  * 5 NOTES
  * 6 SEE ALSO



### USAGE
    
    
    translate vector [,selection [,state [,camera [,object ]]]]
    

### ARGUMENTS

**vector**

    

    float vector: [x, y, z] translation vector

**selection**

    

    string: atoms whose coordinates should be modified {default: all}

**state**

    

    state > 0: only the indicated state is modified
    state = 0: all states are modified
    state = -1: only current state is modified {default}

**camera**

    

    0 or 1: is the vector in camera coordinates? {default: 1 (yes)}

**object**

    

    string: object name (only if rotating object matrix) {default: None}

### PYMOL API
    
    
    cmd.translate(list vector, string selection = "all", int state = 0,
                  int camera = 1, string object = None)
    

### EXAMPLES

#### Example 1

A simple example; just move the alpha-carbons one Ang. along the X-axis. 
    
    
    # move the alpha carbons one Ang. on the X-axis
    translate [1,0,0], name ca
    

#### Example 2

A slightly more complicated example. Load two proteins, align them, and then move one away from the other to see their similarities separately. 
    
    
    # load 1cll and 1ggz
    fetch 1cll 1ggz, async=0
    
    # align the two proteins
    cealign 1cll, 1ggz
    
    # orient the proteins
    orient
    
    # move 1cll above 1ggz so we can 
    # see both proteins separately
    translate [0, 40, 0], 1cll
    

### NOTES

  1. if state = 0, then only visible state(s) are affected.
  2. if state = -1, then all states are affected.



### SEE ALSO

[Rotate](/index.php/Rotate "Rotate"), [Move](/index.php/Move "Move"), [Model_Space_and_Camera_Space](/index.php/Model_Space_and_Camera_Space "Model Space and Camera Space")

Retrieved from "[https://pymolwiki.org/index.php?title=Translate&oldid=11527](https://pymolwiki.org/index.php?title=Translate&oldid=11527)"


---

## Turn

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**turn** rotates the camera about one of the three primary axes, centered at the origin. Also consider [rotate](/index.php/Rotate "Rotate"). 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 SEE ALSO



### USAGE
    
    
    turn axis, angle
    

### EXAMPLES
    
    
    turn x,90
    turn y,45
    

### PYMOL API
    
    
    cmd.turn( string axis, float angle )
    

### SEE ALSO

[move](/index.php/Move "Move"), [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [zoom](/index.php/Zoom "Zoom"), [center](/index.php/Center "Center"), [clip](/index.php/Clip "Clip")

Retrieved from "[https://pymolwiki.org/index.php?title=Turn&oldid=7556](https://pymolwiki.org/index.php?title=Turn&oldid=7556)"


---

## Unbond

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**unbond** removes all bonds between two selections. 

  * [![Atoms bound normally, but not the representation we want.](/images/a/a4/Unbond1.png)](/index.php/File:Unbond1.png "Atoms bound normally, but not the representation we want.")

Atoms bound normally, but not the representation we want. 

  * [![Atoms unboud with the unbound command.](/images/9/9c/Unbond2.png)](/index.php/File:Unbond2.png "Atoms unboud with the unbound command.")

Atoms unboud with the unbound command. 




## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 Example
  * 4 SEE ALSO



# USAGE
    
    
    unbond atom1,atom2
    

# PYMOL API
    
    
    cmd.unbond(selection atom1="(pk1)",selection atom2="(pk2)")
    

# Example
    
    
    # remove all bonds in residue 999 to residue 999
    # this command was used in the examples above in PDB ID 1ACO.
    unbond i. 999, i. 999
    

# SEE ALSO

[Bond](/index.php/Bond "Bond"), [Fuse](/index.php/Fuse "Fuse"), [Remove_picked](/index.php/Remove_picked "Remove picked"), [Attach](/index.php/Attach "Attach"), [Detach](/index.php?title=Detach&action=edit&redlink=1 "Detach \(page does not exist\)"), [Replace](/index.php/Replace "Replace")

Retrieved from "[https://pymolwiki.org/index.php?title=Unbond&oldid=7555](https://pymolwiki.org/index.php?title=Unbond&oldid=7555)"


---

## Undo

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_(Incentive PyMOL only feature)_

Most actions as of **Incentive PyMOL 2.5** are undoable. 

Actions are undoable via GUI menu (Edit -> Undo). The undo feature may be manually disabled by unchecking `Edit -> Undo Enabled`, and can be re-enabled by rechecking the same menu item. 

## Contents

  * 1 Usage
  * 2 PyMOL API
  * 3 Example
  * 4 Comments
  * 5 SEE ALSO



## Usage
    
    
    undo [, steps]
    

  * **steps** = integer: number of steps to undo



## PyMOL API
    
    
    cmd.undo(int steps=1)
    

`undo_enable` is an API-only feature that allows you enable undo. 
    
    
    cmd.undo_enable()
    

`undo_disable` is an API-only feature that allows you disable undo. 
    
    
    cmd.undo_disable()
    

## Example
    
    
     undo
     undo 5
    

## Comments

Currently, the undo stack supports the previous 25 actions or up to 1GB memory used. These will later be configurable in an upcoming release. 

Currently, undo support for builder is not available but will be supported in a later patch release. 

While movie is playing, the undo manager is paused and will save the current session's state when the user ends the movie. 

During wizards, the undo manager is paused, and the current session's state is saved when the user clicks 'Done' for active wizards. 

For larger sessions, it may be suitable to manually turn off undo to save on memory usage (commands such as set will incur a high memory cost). Add `cmd.undo_disable()` in your pymolrc if you need to disable undo automatically upon startup. 

For **Incentive PyMOL** , all previous undo features are deprecated and will now be superseded by the one introduced in **Incentive PyMOL 2.5**. 

## SEE ALSO

[Redo](/index.php/Redo "Redo")

Retrieved from "[https://pymolwiki.org/index.php?title=Undo&oldid=13430](https://pymolwiki.org/index.php?title=Undo&oldid=13430)"


---

## Unmask

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**unmask** reverses the effect of "mask" on the indicated atoms. 

### PYMOL API
    
    
    cmd.unmask( string selection="(all)" )
    

### USAGE
    
    
    unmask (selection)
    

### SEE ALSO

[Cmd mask](/index.php/Cmd_mask "Cmd mask"), [Cmd protect](/index.php/Cmd_protect "Cmd protect"), [Cmd deprotect](/index.php/Cmd_deprotect "Cmd deprotect"), [Cmd mouse](/index.php?title=Cmd_mouse&action=edit&redlink=1 "Cmd mouse \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Unmask&oldid=7553](https://pymolwiki.org/index.php?title=Unmask&oldid=7553)"


---

## Unpick

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**unpick** deletes the special "pk" atom selections (pk1, pk2, etc.) used in atom picking and molecular editing. 

### USAGE
    
    
    unpick
    

### PYMOL API
    
    
    cmd.unpick()
    

### SEE ALSO

[Edit](/index.php/Edit "Edit")

Retrieved from "[https://pymolwiki.org/index.php?title=Unpick&oldid=7552](https://pymolwiki.org/index.php?title=Unpick&oldid=7552)"


---

## Unset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**unset** behaves in two ways. 

If selection is not provided: 

  * Since PyMOL 2.5: Changes the named global setting to the default value.
  * Before PyMOL 2.5: Changes the named global setting to a zero or off value.



If a selection is provided, then "unset" undefines object-specific or state-specific settings so that the global setting will be in effect. 

## Usage
    
    
    unset name [,selection [,state ]]
    

## Python API
    
    
    cmd.unset ( string name, string selection = '',
             int state=0, int updates=1, int log=0 )
    

## See Also

  * [unset_deep](/index.php/Unset_deep "Unset deep")
  * [set](/index.php/Set "Set")



Retrieved from "[https://pymolwiki.org/index.php?title=Unset&oldid=13346](https://pymolwiki.org/index.php?title=Unset&oldid=13346)"


---

## Unset deep

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The unset_deep command clears settings on the object, object-state, atom and bond levels. 

_New in PyMOL 1.8.4_

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Note
  * 5 See Also



## Usage
    
    
    unset_deep [ settings [, object ]]
    

## Arguments

  * **settings** = str: space separated list of setting names or empty string for all settings {default: }
  * **object** = str: name of one object or * for all objects {default: *}



## Example
    
    
    fetch 1rx1, async=0
    as cartoon
    color green
    
    # object-level
    set cartoon_color, blue, 1rx1
    
    # atom-level
    set cartoon_color, red, resi 20-30
    
    # clear on all levels
    unset_deep cartoon_color
    

## Note

Does currently not unset atom-state level settings. 

## See Also

  * [unset](/index.php/Unset "Unset")
  * [set](/index.php/Set "Set")



Retrieved from "[https://pymolwiki.org/index.php?title=Unset_deep&oldid=12446](https://pymolwiki.org/index.php?title=Unset_deep&oldid=12446)"


---

## Update

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**update** transfers coordinates from one selection to another. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 ARGUMENTS
  * 4 EXAMPLES
  * 5 NOTES
  * 6 SEE ALSO



### USAGE
    
    
    update (target-selection),(source-selection)
    

### PYMOL API
    
    
    cmd.update( target ,source [,target_state=0 [,source_state=0 [,matchmaker=1 [,quiet=1 [,_self=cmd ]]]]]] )
    

### ARGUMENTS

  * **matchmaker** = integer: how to match atom pairs {default: 1} 
    * 0: assume that atoms are stored in the identical order
    * 1: match based on all atom identifiers (segi,chain,resn,resi,name,alt)
    * 2: match based on ID
    * 3: match based on rank
    * 4: match based on index



_Note that the meaning of matchmaker=0 differs inupdate and [fit](/index.php/Fit "Fit") (also [rms](/index.php/Rms "Rms"), [rms_cur](/index.php/Rms_Cur "Rms Cur"), ...)!_

### EXAMPLES
    
    
    update target,(variant)
    

### NOTES

Currently, this applies across all pairs of states. Fine control will be added later. 

### SEE ALSO

[Cmd load](/index.php/Cmd_load "Cmd load")

Retrieved from "[https://pymolwiki.org/index.php?title=Update&oldid=12722](https://pymolwiki.org/index.php?title=Update&oldid=12722)"


---

## View

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**view** makes it possible to save and restore viewpoints on a given scene within a single session. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 FUNCTION KEY PRESETS
  * 4 EXAMPLES
  * 5 SEE ALSO



### USAGE
    
    
    view key[,action]
    view *
    

key can be any string action should be 'store' or 'recall' (default: 'recall') 

### PYMOL API
    
    
    cmd.view(string key,string action)
    

### FUNCTION KEY PRESETS

Views F1 through F12 are automatically bound to function keys provided that "set_key" has not been used to redefine the behaviour of the respective key, and that a "scene" has not been defined for that key. 

### EXAMPLES
    
    
    view 0,store
    view 0
    

### SEE ALSO

[Scene](/index.php/Scene "Scene"), [Set View](/index.php/Set_View "Set View"), [Get View](/index.php/Get_View "Get View")

Retrieved from "[https://pymolwiki.org/index.php?title=View&oldid=11633](https://pymolwiki.org/index.php?title=View&oldid=11633)"


---

## Viewport

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**viewport** changes the size of the viewing port--the visible openGL window (and thus the size of all png files subsequently output). 

No API command can reliably retrieve the viewport dimensions under all circumstances. However, it is possible to obtain the dimensions using a third party image viewer like Gimp or OS X Preview: 

  1. Save the current view as a png file ("png imagename.png").
  2. Determine the image dimensions using a viewer program.



These dimensions can be applied directly using the viewport command or the API. 

### USAGE
    
    
    viewport width, height
    

### PYMOL API
    
    
    cmd.viewport(int width, int height)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Viewport&oldid=9467](https://pymolwiki.org/index.php?title=Viewport&oldid=9467)"


---

## Volume

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/b/b1/1okyVol.png)](/index.php/File:1okyVol.png)

[](/index.php/File:1okyVol.png "Enlarge")

Volume visualization of electron density for PDB 1oky

[![](/images/5/5d/1okyVolPanel.png)](/index.php/File:1okyVolPanel.png)

[](/index.php/File:1okyVolPanel.png "Enlarge")

Volume panel for the 1oky volume example. It has the iso-levels on the x-axis and the opacity (1.0 - transparency) on the y-axis.

  
Volume creates a new volume object from a map object. The data (3D scalar fields) are shown as a true 3D object using coloring and transparencies defined by the user to illustrate the data values. This technique supports single and multiple isosurfaces. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Screencasts
  * 5 Changes with PyMOL Version
  * 6 Ray Tracing
  * 7 Known Limitations
  * 8 See Also



## Usage
    
    
    volume name, map [, ramp [, selection [, buffer [, state [, carve ]]]]]
    

## Arguments

  * name = the name for the new volume object.
  * map = the name of the map object to use for computing the volume.
  * ramp = str: named color ramp {default: }
  * selection = an atom selection about which to display the mesh with an additional "buffer" (if provided).
  * carve = a radius about each atom in the selection for which to include density. If "carve" is not provided, then the whole brick is displayed.



## Example
    
    
    fetch 1oky, type=2fofc, async=0
    volume 1okyVol, 1oky_2fofc
    

## Screencasts

  * [Silent demo movie](http://www.youtube.com/watch?v=tuAo_8-_HIc) showing the basics of loading and using a volume in PyMOL. There are more capabilities, but this is the basic functionality.



## Changes with PyMOL Version

  * 1.4.0: first version with volume support
  * 1.7.2: 
    * pre-integrated volume rendering (volume_mode=1) as Incentive-PyMOL-only feature.
    * scripting support with custom color ramp ([volume_color](/index.php/Volume_color "Volume color"), [volume_ramp_new](/index.php?title=Volume_ramp_new&action=edit&redlink=1 "Volume ramp new \(page does not exist\)"))
    * improved volume panel, panel can be opened from the object menu ("C > panel")
    * lots of bugs fixed



## Ray Tracing

**There is no actual ray tracing support**. The volume rendering is implemented exclusively with OpenGL shaders. The recommended way to render a high resolution image is to use the [draw](/index.php/Draw "Draw") command. Example: 
    
    
    # render high resolution image on screen
    draw 4000, 3000, antialias=2
    png highres.png
    

Ray trace specific features like shadows or outlines ([ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")) are not available with this approach. If such features are needed, the [ray_volume](/index.php?title=Ray_volume&action=edit&redlink=1 "Ray volume \(page does not exist\)") setting activates a hybrid solution, which blends the OpenGL rendered volume with the ray traced non-volume objects. Image sizes other than the current window size are not possible. 
    
    
    # compose on-screen volume with ray traced image
    set ray_volume
    ray
    png composed.png
    

Neither of these two solutions work with [headless (batch) mode](/index.php/Launching_PyMOL#Running_PyMOL_in_batch_mode "Launching PyMOL"). 

## Known Limitations

  * No real ray-tracing support yet
  * Multiple volume objects don't blend properly



## See Also

  * <http://pymol.org/volume>
  * <http://pymol.org/d/media:volumevisualization>
  * [volume_color](/index.php/Volume_color "Volume color")
  * [volume_ramp_new](/index.php?title=Volume_ramp_new&action=edit&redlink=1 "Volume ramp new \(page does not exist\)")
  * [map_new](/index.php/Map_new "Map new")
  * [isomesh](/index.php/Isomesh "Isomesh")
  * [isosurface](/index.php/Isosurface "Isosurface")



Retrieved from "[https://pymolwiki.org/index.php?title=Volume&oldid=12536](https://pymolwiki.org/index.php?title=Volume&oldid=12536)"


---

## Volume color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

volume_color set or get the [volume](/index.php/Volume "Volume") colors. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 See Also



## Usage
    
    
    volume_color name [, ramp ]
    

## Arguments

  * name = str: volume object name
  * ramp = str, list or empty: named ramp, space delimited string or list with (x, color, alpha, ...) or (x, r, g, b, alpha, ...) values. If empty, get the current volume colors.



## Examples

Setting volume colors with one line: 
    
    
    fetch 1a00, map, type=2fofc, async=0
    volume vol, map
    volume_color vol, .8 cyan 0. 1. blue .3 2. yellow .3 
    

Using a named color ramp: 
    
    
    volume_ramp_new cyanblueyellow, \
        .8 cyan 0. \
        1. blue .3 \
        2. yellow .3
    volume_color vol, cyanblueyellow
    

Getting the current volume ramp: 
    
    
    PyMOL>volume_color vol
    ### cut below here and paste into script ###
    cmd.volume_ramp_new('ramp399', [\
         0.80, 0.00, 1.00, 1.00, 0.00, \
         1.00, 0.00, 0.00, 1.00, 0.30, \
         2.00, 1.00, 1.00, 0.00, 0.30, \
       ])
    ### cut above here and paste into script ###
    

## See Also

  * [volume](/index.php/Volume "Volume")



Retrieved from "[https://pymolwiki.org/index.php?title=Volume_color&oldid=12729](https://pymolwiki.org/index.php?title=Volume_color&oldid=12729)"


---

## Window

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**window** controls the visibility of PyMOL's output window 

Some actions (maximize, focus, defocus) are not consistent across operating systems. 

_Non-functional in PyMOL 2.0, will be available again in 2.1_

## Contents

  * 1 Usage
  * 2 Actions
  * 3 Example
  * 4 See Also



## Usage
    
    
    window [ action [, ... ]]
    

## Actions

Hide the window (Warning: might make the window inaccessible. Intended for programmed automation): 
    
    
    window hide
    

Show the window (reverse of **hide**): 
    
    
    window show
    

Place the window at **x, y** screen coordinates: 
    
    
    window position, x, y
    

Resize the window: 
    
    
    window size, width, height
    

Place and resize in a single operation: 
    
    
    window box, x, y, width, height
    

If any window corner is not on the visible screen, move the window and if necessary, resize (shrink) to screen dimensions: 
    
    
    window fit
    

Maximize the window: 
    
    
    window maximize
    

Give the OpenGL window focus: 
    
    
    window focus
    

## Example

Place in upper left corner and resize to 1000x500 
    
    
    window box, 0, 0, 1000, 500
    

## See Also

  * [full_screen](/index.php/Full_screen "Full screen")
  * [viewport](/index.php/Viewport "Viewport")
  * -W, -H, -X, -Y [Command Line Options](/index.php/Command_Line_Options "Command Line Options")



Retrieved from "[https://pymolwiki.org/index.php?title=Window&oldid=12726](https://pymolwiki.org/index.php?title=Window&oldid=12726)"


---

## Wizard

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**wizard** launches one of the built-in wizards. There are special Python scripts which work with PyMOL in order to obtain direct user interaction and easily peform complicated tasks. 

### USAGE
    
    
    wizard name
    

### PYMOL API
    
    
    cmd.wizard(string name)
    

### EXAMPLE
    
    
    wizard distance  # launches the distance measurement wizard 
    
    # set a message
    cmd.wizard("message", "Hello, I'm a message.")
    
    # dimiss the message
    cmd.wizard()
    

Retrieved from "[https://pymolwiki.org/index.php?title=Wizard&oldid=7546](https://pymolwiki.org/index.php?title=Wizard&oldid=7546)"


---

## Zoom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**zoom** scales and translates the window and the origin to cover the atom selection. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 PYMOL API
  * 4 NOTES
  * 5 SEE ALSO



### USAGE
    
    
    zoom [ selection [,buffer [, state [, complete ]]]]
    

### EXAMPLES
    
    
    # auto zoom depending on what's loaded in PyMOL
    zoom
    
    #
    zoom complete=1
    
    # zoom on just chain A
    zoom (chain A)
    
    # zoom on residue 142
    zoom 142/
    
    # zoom consistenly 20 Ang from each object at the center
    center prot1
    zoom center, 20
    
    # prot1 and prot2 will have the same exact zoom factor
    center prot2
    zoom center, 20
    

### PYMOL API
    
    
     
    cmd.zoom( string selection, float buffer=0.0,
              int state=0, int complete=0 )
    

### NOTES
    
    
    state = 0 (default) use all coordinate states
    state = -1 use only coordinates for the current state
    state > 0  use coordinates for a specific state
    
    
    
    complete = 0 or 1:
    

Normally the zoom command tries to guess an optimal zoom level for visualization, balancing closeness against occasional clipping of atoms out of the field of view. You can change this behavior by setting the complete option to 1, which will guarantee that the atom positions for the entire selection will fit in the field of an orthoscopic view. To absolutely prevent clipping, you may also need to add a buffer (typically 2 A) to account for the perspective transformation and for graphical representations which extend beyond the atom coordinates. 

### SEE ALSO

  * [Origin](/index.php/Origin "Origin")
  * [Orient](/index.php/Orient "Orient")
  * [Center](/index.php/Center "Center")



Retrieved from "[https://pymolwiki.org/index.php?title=Zoom&oldid=7545](https://pymolwiki.org/index.php?title=Zoom&oldid=7545)"


---

## Category:Atoms Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Atoms Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Atoms_Module&oldid=6887](https://pymolwiki.org/index.php?title=Category:Atoms_Module&oldid=6887)"


---

## Category:Color Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Color Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Color_Module&oldid=6889](https://pymolwiki.org/index.php?title=Category:Color_Module&oldid=6889)"


---

## Category:Colour Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Color Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Color_Module&oldid=6889](https://pymolwiki.org/index.php?title=Category:Color_Module&oldid=6889)"


---

## Category:Display Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Display Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Display_Module&oldid=6884](https://pymolwiki.org/index.php?title=Category:Display_Module&oldid=6884)"


---

## Category:Distances Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Distances Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Distances_Module&oldid=6891](https://pymolwiki.org/index.php?title=Category:Distances_Module&oldid=6891)"


---

## Category:Editing Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Editing Module for PyMol's interface. 

## Pages in category "Editing Module"

The following 16 pages are in this category, out of 16 total. 

### A

  * [Attach](/index.php/Attach "Attach")



### B

  * [Bond](/index.php/Bond "Bond")



### C

  * [Create](/index.php/Create "Create")



### E

  * [Edit](/index.php/Edit "Edit")



### F

  * [Fab](/index.php/Fab "Fab")
  * [Fnab](/index.php/Fnab "Fnab")
  * [Fuse](/index.php/Fuse "Fuse")



### H

  * [H Add](/index.php/H_Add "H Add")
  * [H Fill](/index.php/H_Fill "H Fill")



### P

  * [Protect](/index.php/Protect "Protect")



### R

  * [Redo](/index.php/Redo "Redo")
  * [Remove](/index.php/Remove "Remove")
  * [Remove Picked](/index.php/Remove_Picked "Remove Picked")
  * [Replace](/index.php/Replace "Replace")



### U

  * [Unbond](/index.php/Unbond "Unbond")
  * [Undo](/index.php/Undo "Undo")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Editing_Module&oldid=6839](https://pymolwiki.org/index.php?title=Category:Editing_Module&oldid=6839)"


---

## Category:Fitting Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Fitting Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Fitting_Module&oldid=6888](https://pymolwiki.org/index.php?title=Category:Fitting_Module&oldid=6888)"


---

## Category:Help Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Help Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Help_Module&oldid=6890](https://pymolwiki.org/index.php?title=Category:Help_Module&oldid=6890)"


---

## Category:Imaging Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Imaging Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Imaging_Module&oldid=6881](https://pymolwiki.org/index.php?title=Category:Imaging_Module&oldid=6881)"


---

## Category:Input Output Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Input/Output Module for PyMol's interface. 

## Pages in category "Input Output Module"

The following 6 pages are in this category, out of 6 total. 

### D

  * [Delete](/index.php/Delete "Delete")



### L

  * [Load](/index.php/Load "Load")
  * [Load CGO](/index.php/Load_CGO "Load CGO")



### M

  * [Multifilesave](/index.php/Multifilesave "Multifilesave")



### Q

  * [Quit](/index.php/Quit "Quit")



### S

  * [Save](/index.php/Save "Save")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Input_Output_Module&oldid=6855](https://pymolwiki.org/index.php?title=Category:Input_Output_Module&oldid=6855)"


---

## Category:Language Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Language Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Language_Module&oldid=6895](https://pymolwiki.org/index.php?title=Category:Language_Module&oldid=6895)"


---

## Category:Maps Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Maps Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Maps_Module&oldid=6883](https://pymolwiki.org/index.php?title=Category:Maps_Module&oldid=6883)"


---

## Category:Movies Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Movies Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Movies_Module&oldid=6880](https://pymolwiki.org/index.php?title=Category:Movies_Module&oldid=6880)"


---

## Category:Ray Tracing Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Ray Tracing Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Ray_Tracing_Module&oldid=6882](https://pymolwiki.org/index.php?title=Category:Ray_Tracing_Module&oldid=6882)"


---

## Category:Scripts Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Scripts Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Scripts_Module&oldid=6894](https://pymolwiki.org/index.php?title=Category:Scripts_Module&oldid=6894)"


---

## Category:Selections Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Selections Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Selections_Module&oldid=6885](https://pymolwiki.org/index.php?title=Category:Selections_Module&oldid=6885)"


---

## Category:Settings Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Settings Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Settings_Module&oldid=6886](https://pymolwiki.org/index.php?title=Category:Settings_Module&oldid=6886)"


---

## Category:Stereo Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Stereo Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Stereo_Module&oldid=6892](https://pymolwiki.org/index.php?title=Category:Stereo_Module&oldid=6892)"


---

## Category:Symmetry Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the Symmetry Module for PyMol's interface. 

_This category currently contains no pages or media._

Retrieved from "[https://pymolwiki.org/index.php?title=Category:Symmetry_Module&oldid=6893](https://pymolwiki.org/index.php?title=Category:Symmetry_Module&oldid=6893)"


---

## Category:View Module

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a list of all commands in the View Module for PyMol's interface. 

## Pages in category "View Module"

The following 18 pages are in this category, out of 18 total. 

### C

  * [Clip](/index.php/Clip "Clip")



### D

  * [Disable](/index.php/Disable "Disable")



### E

  * [Enable](/index.php/Enable "Enable")



### G

  * [Get View](/index.php/Get_View "Get View")



### H

  * [Hide](/index.php/Hide "Hide")



### M

  * [Move](/index.php/Move "Move")



### O

  * [Orient](/index.php/Orient "Orient")
  * [Origin](/index.php/Origin "Origin")



### R

  * [Rebuild](/index.php/Rebuild "Rebuild")
  * [Refresh](/index.php/Refresh "Refresh")
  * [Reset](/index.php/Reset "Reset")
  * [Rock](/index.php/Rock "Rock")



### S

  * [Set View](/index.php/Set_View "Set View")
  * [Show](/index.php/Show "Show")
  * [Show as](/index.php/Show_as "Show as")



### T

  * [Turn](/index.php/Turn "Turn")



### V

  * [View](/index.php/View "View")



### Z

  * [Zoom](/index.php/Zoom "Zoom")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:View_Module&oldid=6861](https://pymolwiki.org/index.php?title=Category:View_Module&oldid=6861)"


---

