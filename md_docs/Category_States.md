# Category: States

## All states

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Example
  * 4 User Notes
  * 5 See Also



## Overview

When set "on", this setting causes PyMOL to display all states or in NMR jargon: all the models in the ensemble. The 'default' behavior (OFF) can be overridden by placing the "set all_states, on" statement into your '.pymolrc' file, located in your login directory (under all flavors of unix). 

## Syntax
    
    
    set all_states, on      
    set all_states, off
    

## Example
    
    
    # fetch a PDB and show it in multiple states; this one
    # line does the work of the next 6 lines.
    fetch 1nmr
    
    # this is older code, use the above code which is newer
    import urllib2
    pdbCode = '1BRV'
    pdbUrl = 'http://www.rcsb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId='+pdbCode
    pdbFile = urllib2.urlopen(pdbUrl)
    pdbContent = pdbFile.read()
    cmd.read_pdbstr(pdbContent, pdbCode)
    
    set all_states, on
    

This shows the effect of turning on/off the **all_states** setting used with the script above. 

  * [![set all_states, on](/images/7/76/All_states_on.png)](/index.php/File:All_states_on.png "set all_states, on")

set all_states, on 

  * [![set all_states, off](/images/6/6f/All_states_off.png)](/index.php/File:All_states_off.png "set all_states, off")

set all_states, off 




## User Notes

  * There was an error with importing ensembles of ensembles files (complex MOL2 files) before revision 3541. If you experience this problem, update to 3541 or later.



## See Also

[Split_States](/index.php/Split_States "Split States")

Retrieved from "[https://pymolwiki.org/index.php?title=All_states&oldid=5749](https://pymolwiki.org/index.php?title=All_states&oldid=5749)"


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

## Defer builds mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

The defer_builds_mode setting for improved [performance](/index.php/Category:Performance "Category:Performance") with long trajectories. This now makes it possible to work with files containing thousands of states, and to render impossibly long movies piecewise. This setting, as shown below, stops PyMOL caching the geometry of trajectory data in RAM. 

# Usage
    
    
    # improve PyMOL performance for many-state objects and long movies.
    set defer_builds_mode, 3
    

# See Also

  * [Load_Traj](/index.php/Load_Traj "Load Traj")
  * [async_builds](/index.php/Async_builds "Async builds")



Retrieved from "[https://pymolwiki.org/index.php?title=Defer_builds_mode&oldid=5801](https://pymolwiki.org/index.php?title=Defer_builds_mode&oldid=5801)"


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

## Editing atoms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Altering atom types/elements
  * 2 Altering atom coordinates
  * 3 Translate or rotate individual objects
  * 4 Get coordinates from Python



## Altering atom types/elements

Example: 

To transform a carbon into an oxygen, pick an atom by a _select_ command or CTRL+middle. 
    
    
    alter pk1,elem='O'
    alter pk1,name='O2'
    

Note that the _name_ field should contain something more descriptive than just the element symbol. For reasons unknown (?), util.cbag() and similar methods will not recognize the new element. 

## Altering atom coordinates

Example: 
    
    
    alter_state 1,(pdb1cse),x=x-10.0
    

The latter section can contain formulae involving at least the xyz coordinates, lots of constants and the (+-*/) operators. 

  


## Translate or rotate individual objects

There is a "translate" function similar to "rotate", the docs for these don't exist yet, because the implementation isn't finished. However, feel free to use them in the following forms: 
    
    
    translate vector,object-name,state
       vector needs to be something like [x,y,z]
    
       translate [1,0,0],pept
    
    
    
    rotate axis,angle,object-name,state
       axis can be either the letter x,y,z or a 3D vector [x,y,z]
    
       rotate x,90,pept
       rotate [1,1,1],10,pept
    

## Get coordinates from Python

  1. The actual C-langauge arrays aren't exposed, but there are at least three different ways you can modify coordinates from within Python: You can get a python object which contains the molecular information, modify the coordinates in that object, load the modified molecule into PyMOL, update the modified coordinates to the original model, and then delete the modified object. (link to example)
  2. Another approach is the "alter_state" function, which can perform the same transformation in a single PyMOL command statement:
         
         alter_state 1,pept,(x,y)=(-y,x)
         

Likewise sub-selections can be transformed as well:
         
         alter_state 1,(pept and name ca),(x,y,z)=(x+5,y,z)
         

  3. A third approach is to use alter_state with the global "stored" object: Example in a Python script)



Approaches 2 gives the best [performance](/index.php/Category:Performance "Category:Performance"), approach 3 gives more flexibility, and approach 1 gives you a reusable and fully modifiable Python object. 

Retrieved from "[https://pymolwiki.org/index.php?title=Editing_atoms&oldid=10793](https://pymolwiki.org/index.php?title=Editing_atoms&oldid=10793)"


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

## Main Page

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

hosted by [![SBGridlogo2.jpg](/images/e/ec/SBGridlogo2.jpg)](https://sbgrid.org/)  
---  
Welcome to the PyMOL Wiki!  The community-run support site for the [PyMOL](http://pymol.org) molecular viewer.   
---  
To request a new account, email SBGrid at: accounts (@) sbgrid dot org   
Quick Links  **[Tutorials](/index.php/Category:Tutorials "Category:Tutorials")** | **[Table of Contents](/index.php/TOPTOC "TOPTOC")** | **[Commands](/index.php/Category:Commands "Category:Commands")**  
---|---|---  
**[Script Library](/index.php/Category:Script_Library "Category:Script Library")** | **[Plugins](/index.php/Category:Plugins "Category:Plugins")** | **[FAQ](/index.php/Category:FAQ "Category:FAQ")**  
**[Gallery](/index.php/Gallery "Gallery")** | **[Covers](/index.php/Covers "Covers")** | **[PyMOL Cheat Sheet](/index.php/CheatSheet "CheatSheet")** (_[PDF](/images/7/77/PymolRef.pdf "PymolRef.pdf")_)  | **[Getting Help](/index.php/PyMOL_mailing_list "PyMOL mailing list")**  
News & Updates  | New Setup  | [PyMOL-open-source-windows-setup v3.1](https://github.com/kullik01/pymol-open-source-windows-setup/releases/tag/v3.1.0) has been released on January 20, 2025. More information under [Windows Install](/index.php/Windows_Install "Windows Install").   
---|---  
New Plugin  | [PySSA](/index.php/PySSA "PySSA") aims to combine PyMOL and [ColabFold](https://github.com/sokrypton/ColabFold) to enable the prediction and analysis of 3D protein structures for the scientific end-user. [v1.0 has been released](https://github.com/urban233/PySSA/releases/tag/v1.0.1) on July 10, 2024.   
Official Release  | [PyMOL v3.0 has been released](https://pymol.org) on March 12, 2024.   
New Plugin  | [CavitOmiX](/index.php/CavitOmiX "CavitOmiX") calculate [Catalophore™ cavities](https://innophore.com), predict protein structures with [OpenFold by NVIDIA-BioNeMo](https://www.nvidia.com/en-us/gpu-cloud/bionemo), [ESMFold](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/) and retrieve [Alphafold](https://www.deepmind.com/research/highlighted-research/alphafold) models   
Official Release  | [PyMOL v2.5 has been released](https://pymol.org) on May 10, 2021.   
Python 3  | New [Python 3 compatibility guide](/index.php/2to3 "2to3") for scripts and plugins   
POSF  | [New PyMOL fellowship announced for 2022-2023](https://pymol.org/fellowship)  
Tutorial  | [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial") updated for PyQt5   
New Plugin  | [PICv](/index.php/PICv "PICv") is a new plugin for clustering protein-protein interactions and visualization with available data from PDBe   
Selection keywords  | New [polymer.protein and polymer.nucleic](/index.php/Selection_Algebra "Selection Algebra") selection keywords. Thanks everyone who participated in the [poll](https://goo.gl/forms/r0Ck03VTytZQxN4A2)!   
Plugin Update  | [MOLE 2.5](/index.php/MOLE_2.0:_advanced_approach_for_analysis_of_biomacromolecular_channels "MOLE 2.0: advanced approach for analysis of biomacromolecular channels") is an updated version of channel analysis software in PyMOL   
New Script  | [dssr_block](/index.php/Dssr_block "Dssr block") is a wrapper for DSSR (3dna) and creates block-shaped nucleic acid cartoons   
Older News  | See [Older News](/index.php/Older_News "Older News").   
Did you know...  | 

### [CASTp](/index.php/CASTp "CASTp")

| Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/castp.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/castp.py)  
Author(s)  | Joe Dundas   
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview

The CASTpyMOL plugin allows for the visualization of surfaces and voids identified from CASTp using pyMOL. For more information see the [CASTp Website](http://sts-fw.bioengr.uic.edu/castp/pymol.php).   
  
|  [![](/images/7/74/0712channels.jpg)](/index.php/File:0712channels.jpg) [](/index.php/File:0712channels.jpg "Enlarge")A Random PyMOL-generated Cover. See [Covers](/index.php/Covers "Covers").   
  
  
Retrieved from "[https://pymolwiki.org/index.php?title=Main_Page&oldid=13837](https://pymolwiki.org/index.php?title=Main_Page&oldid=13837)"


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

## Practical Pymol for Beginners

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Although PyMOL has a powerful and flexible interface, it is complex, and can appear daunting to new users. This guide is intended to introduce the PyMOL interface and basic tasks without leaving the mouse behind. 

## Contents

  * 1 The PyMOL Interface
    * 1.1 About the command line
  * 2 Getting started: explore a protein
    * 2.1 Related Commands
  * 3 What else can this thing do?
    * 3.1 Saving an image
    * 3.2 Selecting parts of an object
    * 3.3 Whoops: Getting unstuck
    * 3.4 Related Commands
    * 3.5 Saving work
      * 3.5.1 Sessions
      * 3.5.2 Molecules
      * 3.5.3 Images
      * 3.5.4 Movies
    * 3.6 Scripting
      * 3.6.1 Introduction and Very Simple Scripting
      * 3.6.2 The Python MiniShell
      * 3.6.3 Learning More...
  * 4 Coming Soon
    * 4.1 Logging



## The PyMOL Interface

When PyMOL is opened, two windows appear. The smaller window (called the "External GUI" in PyMOL documentation) contains the menu bar (**File** , **Edit** , **Help** , **Display** , etc), shortcut buttons for common commands, and the command line. 

[![](/images/2/20/Viewer_guide.png)](/index.php/File:Viewer_guide.png)

[](/index.php/File:Viewer_guide.png "Enlarge")

The PyMol Viewer Window

The second window is the PyMOL Viewer, which is where all the magic happens. In the Viewer, 3D models are displayed, and the user interacts (eg rotates) and manipulates the model. 

The objects that PyMOL renders in 3D are loaded from coordinate files that describe (in great detail) locations of individual atoms in the molecule. PyMOL can display more than one object at a time, and provides an Object Control Panel to adjust viewing modes, colors, labels, hiding, and just about anything else relating to objects. After each object name is a set of command buttons which control the object. Here are the buttons and some of their options: 

  * **A** \- _Actions_ : Rename, duplicate, remove, apply presets (like "ball-and-stick" or "publication"), perform computations
  * **S** \- _Show_ : Change the way things appear, eg change to stick or cartoon view.
  * **H** \- _Hide_ : Things that are shown using **S** accumulate, and don't automatically replace the last view. **H** is the opposite of **S** and hides unwanted representations.
  * **L** \- _Label_ : Label atoms, residues, etc.
  * **C** \- _Color_ : Change the color of atoms and groups.



The lower-right corner of the Viewer contains a guide to using the mouse, as well as a powerful selection tool. There is also another command line at the bottom of the Viewer (**PyMOL >**). 

### About the command line

The PyMOL command line is a great tool that lets the experienced user change all sorts of options that simply don't appear in the point-and-click graphical interface. It can also be a lot faster. Combined with scripting, it is a powerful option for automating tasks and making intricate sets of changes. But, it's complex, and page upon page of PyMOL documentation cover these commands, so we're going to ignore them as much as possible. 

Although this guide may include some text commands and links to more advanced documentation, they're purely optional and meant to be informative. 

To run any text command, type it in at a **PyMOL >** command line and hit _[Enter]_. 

## Getting started: explore a protein

PyMOL is great for casual visualization of biological molecules. In this example, a PDB file describing a protein is loaded and its style and color are tweaked. 

[![](/images/7/75/Kchannel-rainbow.png)](/index.php/File:Kchannel-rainbow.png)

[](/index.php/File:Kchannel-rainbow.png "Enlarge")

The end result will look something like this

[![](/images/5/54/Mouse-3view.png)](/index.php/File:Mouse-3view.png)

[](/index.php/File:Mouse-3view.png "Enlarge")

Default buttons for viewing with a 3-button mouse

  1. Obtain a PDB coordinates file for your favorite protein. (The [RCSB Protein Data Bank](http://www.pdb.org/) is a public structure repository containing over 40,000 protein structures in PDB format available for download, not a bad place to look.) For this example, we're using the potassium channel from _Streptomyces Lividans_ ([1BL8](http://www.pdb.org/pdb/files/1bl8.pdb)). 
     * or just type, 
           
           fetch 1bl8
           

and you can skip the next step (see [Fetch](/index.php/Fetch "Fetch") command), as PyMOL will open the file for you.
  2. Open the PDB file using **File** => **Open...** from the menu bar. The protein's structure will appear, probably rendered as simple bonding lines.
  3. The right side of the Viewer shows the loaded PDB as an object, as well as its command buttons. Each button contains a submenu with more options. Click **S** , then **cartoon** to show the protein's secondary structure in popular cartoon form. 
     * Notice that the lines view is still visible on top of the cartoon view. To hide the lines, click **H** then **lines**.
  4. To change the color of each protein chain (as defined in the coordinate file), click **C** then select **chainbows** from the **by chain** menu. "Chainbows" colors residues in each protein chain as a rainbow that begins with blue and ends with red. 
     * Another common coloring method assigns a single color to each chain. Click **C** then select **by chain** from the **by chain** menu.
  5. Click and drag the protein to change the view. A list of mouse buttons is below the object control panel. 
     * Rota: Rotate
     * Move: Translate object along an axis
     * MoveZ: aka Zoom
     * Sele: Select
     * Slab:
     * Cent:
     * PkAt:



### Related Commands

[Load](/index.php/Load "Load"), [Fetch](/index.php/Fetch "Fetch"), [Color](/index.php/Color "Color"), [Show](/index.php/Show "Show"), [Show_as](/index.php/Show_as "Show as"), [Cartoon](/index.php/Cartoon "Cartoon"), [Lines](/index.php/Lines "Lines"), [Rotate](/index.php/Rotate "Rotate"), [Select](/index.php/Select "Select"), [Center](/index.php/Center "Center")

## What else can this thing do?

So, now what? Good question. PyMOL is a powerful program, and everyone uses it for something different. The remainder of this guide is devoted to common tasks that come in handy. 

### Saving an image

You've found the perfect view, and you'd like to [Save](/index.php/Save "Save") it? 

  1. Change the size of the viewer and zoom in/out to make the image the right size. Images are saved exactly as they appear in the viewer, and the size of the viewer determines the size of the image. 
     * For example, images for printing and presenting should be larger than images for a posting on a website.
  2. Because PyMOL was designed to run on older computers, the maximum quality is not enabled by default. To change this, select **Display** , **Quality** , **Maximum Quality**. Notice that round things are rounder, curves are curvier, and color shading is smoother.
  3. Once you've found the appropriate view, save an image using **File** , **Save Image...** An picture of the current view will be saved in PNG format.



**Tip:** Using the _[ray](/index.php/Ray "Ray")_ command before saving an image will create a higher quality version with shadows, etc. This can take time, depending on the size of the image and speed of the computer, but the images created after ray tracing are usually spectacular. However, the ray tracing disappears if the view is changed in any way. 

### Selecting parts of an object

Sometimes it might be useful to select part of an object and modify it directly; for example, selecting active-site residues in a protein to highlight them in another color. 

In the lower-right corner of the Viewer, next to the animation controls, is an **S** button (not to be confused with the **S** how button in the object control panel) that activates the selection tool. The selection tool can be changed to select atoms or residues by clicking _Selecting Residues(or whatever)_ until the right mode appears. 

Once selecting is activated, a list of parts to select appears at the top of the Viewer. Select things clicking or dragging across a range. 

Selections can be controlled individually in the object control panel, just like any other object. To save a selection, select **rename** from the **A** menu. 

### Whoops: Getting unstuck

PyMOL is a program meant to be explored and played with, and sometimes funny things happen in the process. A few common problems: 

  * _The model disappeared:_ Sometimes while rotating and moving a model, it can get lost. Right-click the background of the viewer, and select **reset** from the _Main Pop-Up_. The model should return to view; if it doesn't, make sure the object is being drawn using the **S** menu.
  * _The model has funny colors, labels, etc and they won't go away:_ The **H** menu in the object control panel will remove unwanted details; however, sometimes it's difficult to know exactly what to remove. Select **H** , then **everything** to hide all details and start fresh.
  * _Things are really messed up:_ Use **File** , **Reinitalize** to reset PyMOL to its initial state, but all work will be lost.



### Related Commands

[Save](/index.php/Save "Save"), [Viewport](/index.php/Viewport "Viewport"), [Zoom](/index.php/Zoom "Zoom"), [Save](/index.php/Save "Save"), [Ray](/index.php/Ray "Ray"), [Select](/index.php/Select "Select")

### Saving work

PyMOL supports saving your work in various formats. You can save, images, molecules, sessions, movies, etc. 

#### Sessions

A PyMOL sessions retains the state of your PyMOL instance. You can save the session to disk and reload it later. You can setup a complicated scene, with transitions and more, and simply save it as a PyMOL Session (.pse) file. Use, **File** =>**Save Session** or **Save Session As...**. 

Loading sessions is just as easy: **File** =>**Load** and choose your session file. 

#### Molecules

You can save a molecule by selecting **File** =>**Save Molecule**. A dialog box will appear from which you can select your molecule to save. You can also save an object or selection using the [Save](/index.php/Save "Save") command. It's very easy: 
    
    
    save  fileName, objSel
    

where fileName is something like "1abc.pdb", and objSel can be any object or selection. For example, to save just the waters to disk do: 
    
    
    save wat.pdb, resn HOH
    

#### Images

You can save images that you've rendered (with [Ray](/index.php/Ray "Ray")) or drawn (with [Draw](/index.php/Draw "Draw")) again using the [Save](/index.php/Save "Save") command or by **File** =>**Save Image**. You can save in [Png](/index.php/Png "Png"), VRML-2 and the POVRay formats. 

You can save images to disk, through the command line, by using the [Png](/index.php/Png "Png") command. 

#### Movies

PyMOL allows you to save movies you've created, too. You can automatically save the MPEG or save a series of PNG files--for stitching together later. This is a new command, and I don't know too much about it. Use **File** =>**Save Movie**. 

### Scripting

#### Introduction and Very Simple Scripting

Scripting in PyMOL ranges from very simple to very intricate. You can make a simple loop (say rotating a molecule one-degree at a time) or execute full featured scripts. Once you've got the basics of PyMOL down, learning scripting will greatly enhance your productivity and capabilities. 

Because of the way PyMOL is built on top of the Python interpreter, any command that PyMOL doesn't recognize it passes on to Python to execute. This is a **very** handy feature--you essentially have a live Python interpreter you can interact with, which makes your life much easier. Take the following for example: 
    
    
    f = 10.
    for x in range(0,100,10):
      cmd.set("spec_direct_power", float(float(x) / f))
      cmd.png("spec_dir_power" + str(x) + ".png", ray=1)
    

This simple script of 4 lines will create 10 images each one with the [Spec_direct_power](/index.php/Spec_direct_power "Spec direct power") changed (see the [Spec_direct_power](/index.php/Spec_direct_power "Spec direct power") for the output of this script; the animated GIF). Did you notice that **f** on line 1 and **for** and **x** on line 2 are not PyMOL commands or symbols? You can simply write Python code that interacts with PyMOL. Brilliant! 

#### The Python MiniShell

Taking this one level higher, you can write little code snippets, like 10-20+ lines of code that perform some specific task and wrap these in the `python` and `python end` commands. (If your script ever makes a mistake, you can abort the endeavor with `end` instead of `python end`. The power here is that none of your command are executed until you type `python end`. Confused? Here's the above example in using the wrapper: 
    
    
    python
    f = 10.
    for x in range(0,100,10):
      cmd.set("spec_direct_power", float(float(x) / f))
      cmd.png("spec_dir_power" + str(x) + ".png", ray=1)
    python end
    

The `[python](/index.php/Python "Python")` command gives you complete access to the Python shell and `python end` brings you back into PyMOL's shell. Also, note that PyMOL saves information across instantiations of the `python` command. For example, 
    
    
    # enter mini python shell
    python
    ff = 10.
    python end
    
    # now we're back in the normal PyMOL shell, where PyMOL knows about the value
    print(ff)
    
    # in the mini shell, Python still knows about ff.
    python 
    print(ff)
    python end
    

#### Learning More...

To learn more about scripting check out: 

  * [Biochemistry_student_intro](/index.php/Biochemistry_student_intro "Biochemistry student intro") Basic use of GUI and script
  * [Simple_Scripting](/index.php/Simple_Scripting "Simple Scripting") introduction
  * [Advanced_Scripting](/index.php/Advanced_Scripting "Advanced Scripting") pages
  * [Popular Script Library](/index.php/Category:Script_Library "Category:Script Library").



## Coming Soon

### Logging

Retrieved from "[https://pymolwiki.org/index.php?title=Practical_Pymol_for_Beginners&oldid=12744](https://pymolwiki.org/index.php?title=Practical_Pymol_for_Beginners&oldid=12744)"


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

## Ribbon color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 Advanced details



## Overview

Ribbon_color allows one to explicitly state the color to be applied to a ribbon object. 

## Syntax
    
    
    # set it to a color
    set ribbon_color, color
    
    # examples:
    # default auto-cycles the color for each new object
    set ribbon_color, marine
    
    # apply green to the ribbon only for object obj01
    set ribbon_color, green /obj01
    
    # apply orange to the ribbon only for chain B of object obj01; NB: not fully enabled until version 1.0
    set ribbon_color, orange /obj01//B
    

## Examples

  * [![ribbon_color green](/images/1/17/Ribbon_color_green.png)](/index.php/File:Ribbon_color_green.png "ribbon_color green")

ribbon_color green 

  * [![ribbon_color marine](/images/0/0d/Ribbon_color_marine.png)](/index.php/File:Ribbon_color_marine.png "ribbon_color marine")

ribbon_color marine 

  * [![two chains in same object with different colors](/images/f/f1/Ribbon_color_chains.png)](/index.php/File:Ribbon_color_chains.png "two chains in same object with different colors")

two chains in same object with different colors 




  


## Advanced details

Some notes on general [Set](/index.php/Set "Set") syntax (From PyMOL "help set"): 

DESCRIPTION 
    
    
      "set" changes one of the PyMOL state variables,
    
    

USAGE 
    
    
      set name, [,value [,object-or-selection [,state ]]]
    
      set name = value  # (DEPRECATED)
    
    

PYMOL API 
    
    
      cmd.set ( string name, string value=1,
                string selection=_, int state=0,_
                 int updates=1, quiet=1)
    
    

NOTES 
    
    
      The default behavior (with a blank selection) changes the global
      settings database.  If the selection is 'all', then the settings
      database in all individual objects will be changed.  Likewise, for
      a given object, if state is zero, then the object database will be
      modified.  Otherwise, the settings database for the indicated state
      within the object will be modified.
    
      If a selection is provided, then all objects in the selection will
      be affected.
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ribbon_color&oldid=6230](https://pymolwiki.org/index.php?title=Ribbon_color&oldid=6230)"


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

## Save sep

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Description

Saves all objects as separate PDB files. Useful if you want to do things like combine separate *.pse files. 

## Code
    
    
    from pymol import cmd
    import glob
    import re
    
    def save_sep(prefix=''):
      """
      save_sep <prefix>
    
      saves multiple objects into multiple files using an optional prefix name.
    
      e.g. save_sep prefix
      """
      obj_list = cmd.get_names("all")
    
      if obj_list:
        for i in range(len(obj_list)):
          obj_name = "%s%s.pdb" % (prefix, obj_list[i])
          cmd.save(obj_name, obj_list[i])
          print "Saving %s" %  obj_name
      else:
        print "No objects found"
        
    cmd.extend('save_sep',save_sep)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Save_sep&oldid=6512](https://pymolwiki.org/index.php?title=Save_sep&oldid=6512)"


---

## Scene

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**scene** makes it possible to save and restore multiple scenes scene within a single session. A scene consists of the view, all object activity information, all atom-wise visibility, color, representations, and the global frame index. 

## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Using Scene
    * 2.1 Storing scenes
    * 2.2 Scenes as Movies
    * 2.3 Auto-play through Scenes
  * 3 Examples
  * 4 PyMOL API
  * 5 Notes
  * 6 See Also
  * 7 DEVELOPMENT TO DO



## Usage
    
    
    scene [ key [, action [, message [, view [, color [, active [, rep
        [, frame [, animate [, new_key ]]]]]]]]]]
    

### Arguments

  * **key** = string, new, auto, or *: use new for an automatically numbered new scene, use auto for the current scene (if one exists), and use * for all scenes (clear and recall actions only).
  * **action** = store, recall, insert_after, insert_before, next, previous, update, rename, clear or append: (default = recall). If rename, then a new_key argument must be explicitly defined.
  * **message** = string: a text message to display with the scene.
  * **view** = 1 or 0: controls whether the view is stored {default: 1}
  * **color** = 1 or 0: controls whether colors are stored {default: 1}
  * **active** = 1 or 0: controls whether activity (objects enabled/disabled) is stored {default: 1}
  * **rep** = 1 or 0: controls whether the representations are stored {default: 1}
  * **frame** = 1 or 0: controls whether the frame is stored {default: 1}
  * **animate** = float: animation duration in seconds {default: _scene_animation_duration_}
  * **new_key** = string: the new name for the scene



## Using Scene

The Scene command has quite a few actions/options that can be enabled by using the mouse and the keyboard through the usual Scene command or hot-keys. Also, you can shift the scenes around using the new [Scene_buttons](/index.php/Scene_buttons "Scene buttons") and just dragging the scene names. 

### Storing scenes
    
    
    # store this scene in the next spot, giving it the default name.
    scene auto, store
    

has the hot-key equivalent of **CTRL-PageDown** (FN+CTRL+DownArrow on the Mac). Try turning on [Scene_Buttons](/index.php?title=Scene_Buttons&action=edit&redlink=1 "Scene Buttons \(page does not exist\)") and then doing CTRL-PageDown; see the scene buttons popping up? 

### Scenes as Movies

If you desire to make a movie that only has camera changes or representation changes, then scenes are your friend. Simply setup each view and then when ready you can do Scene->Store from the PyMOL menus (or _scene auto, store_ on the command line or the third method Ctrl+PgDn (Fn+Ctrl+DownArrow on the Mac)). Do this for each view you setup. Once done, you can scroll through your scenes by pushing PgUp/PgDn. PyMOL automatically interpolates when you use the PgUp/PgDn buttons, so you get the desired smooth transitions. Mix this with [AxPyMOL](http://www.pymol.org/ax/) and you have movies in PowerPoint with very little work. 

### Auto-play through Scenes

With this simple trick you can auto-play through scenes. This is similar to "Movie > Program > Scene Loop" but uses only a single frame. 
    
    
    cmd.mset('1x1')
    cmd.set('scene_loop')
    cmd.set('movie_fps', 1.0 / 5.0)
    cmd.mdo(1, 'scene auto, next')
    cmd.mplay()
    

## Examples

Simple Examples. 
    
    
    scene F1, store
    scene F2, store, This view shows you the critical hydrogen bond.
     
    scene F1
    scene F2
    
    scene *
    

This example shows how to use scenes in a movie! 
    
    
    # SUMMARY
    #
    
    # This script demonstrates one way of creating a movie from scenes.
    # It assumes that we have three scenes, each running for 10 seconds
    # (300 frames apiece) including 2-second transitions.
    
    # 1) Load or create content for three scenes (this could just as easily
    #    come from a session file).
    
    load $TUT/1hpv.pdb
    util.cbc
    turn x,180
    orient
    as cartoon
    scene 001, store
    
    show sticks, organic
    orient organic
    scene 002, store
    
    hide cartoon
    show lines, byres organic expand 5
    turn x,45
    turn y,45
    scene 003, store
    
    # 2) Specify a 30-second movie -- state 1, 900 frames at 30 frames per second.
    
    mset 1 x900
    
    # 3) Program scene matrices as movie views at appopriate frames
    #    and also add y-axis rocking between scenes.
    
    scene 001, animate=0
    mview store, 1
    mview store, 240
    
    turn y,-30
    mview store, 70
    turn y,60
    mview store, 170
    
    scene 002, animate=0
    mview store, 300
    mview store, 540
    
    turn y,-30
    mview store, 370
    turn y,60
    mview store, 470
    
    scene 003, animate=0
    mview store, 600
    mview store, 840
    
    turn y,-30
    mview store, 670
    turn y,60
    mview store, 770
    
    # 4) Now interpolate the movie camera.
    
    mview interpolate
    mview smooth
    mview smooth
    
    # 5) Activate scene content at the appropriate movie frames.
     
    mdo 1: scene 001, view=0, quiet=1
    mdo 240: scene 002, view=0, quiet=1
    mdo 540: scene 003, view=0, quiet=1
    mdo 840: scene 001, view=0, quiet=1
    
    # 6) Force frame 1 content to load.
    
    rewind
    
    # 6) And play the movie.
    
    mplay
    

## PyMOL API
    
    
    cmd.scene(str key='auto', str action='recall', str-or-list message=None, bool view=1, bool color=1,
        bool active=1, bool rep=1, bool frame=1, float animate=-1, str new_key=None)
    

## Notes

  * To scroll through your frames, as in a presentation, just use the PG-UP and PG-DN keys. Very handy.
  * Scenes F1 through F12 are automatically bound to function keys provided that "set_key" hasn't been used to redefine the behaviour of the respective key.
  * If you have a script that modifies the representation of the molecules and stores them, quickly, then the stored frames may not be up to date. I suggest calling "refresh" between the commands.



## See Also

[View](/index.php/View "View"), [Set_View](/index.php/Set_View "Set View"), [Get_View](/index.php/Get_View "Get View"), [Movie_from_scenes](/index.php/Movie_from_scenes "Movie from scenes")

## DEVELOPMENT TO DO

Add support for save/restore of a certain global and object-and-state specific settings, such as: state, surface_color, ribbon_color, stick_color, transparency, sphere_transparency, etc. This would probably best be done by defining a class of "scene" settings which are treated in this manner. The current workaround is to create separate objects which are enabled/disabled differentially. 

Retrieved from "[https://pymolwiki.org/index.php?title=Scene&oldid=12266](https://pymolwiki.org/index.php?title=Scene&oldid=12266)"


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

## Set Title

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**set_title** attaches a text string to the state of a particular object which can be displayed when the state is active. This is useful for display the energies of a set of conformers. 

### USAGE
    
    
    set_title object,state,text
    

### PYMOL API
    
    
    cmd.set_title(string object,int state,string text)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Set_Title&oldid=7575](https://pymolwiki.org/index.php?title=Set_Title&oldid=7575)"


---

## Simple Scripting

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_This page discusses writing your own simple scripts for PyMOL. If you're looking to download scripts try the[Script Library](/index.php/Category:Script_Library "Category:Script Library"). If you need help with running a script try [Running_Scripts](/index.php/Running_Scripts "Running Scripts") _

One of the more powerful features of PyMOL is that it supports Python scripting. This gives you the power to use most of the [Python libraries](http://docs.python.org/api/api.html) to write programs and then send the results back into PyMOL. Some useful extensions to PyMOL can be found in our [Script Library](/index.php/Category:Script_Library "Category:Script Library"). 

## Contents

  * 1 General Scripts
    * 1.1 Getting PyMOL Data into your Script
      * 1.1.1 Getting Data From your Script into PyMOL
    * 1.2 Example
    * 1.3 Basic Script Body
  * 2 Python Version Supported by PyMOL



# General Scripts

General PyMOL scripting is done in Python. It's really quite simple, just write your function (following a couple simple rules) and then let PyMOL know about it by using the **[cmd.extend](/index.php/Extend "Extend")** command. Here's the simple recipe for writing your own simple scripts for PyMOL: 

**To write them** : 

  1. Write the function, let's call it **doSimpleThing** , in a Python file, let's call the file **pyProgram.py**.
  2. Add the following command to the end of the **pyProgram.py** file 
         
         cmd.extend("doSimpleThing",doSimpleThing)
         




**To use them** : 

  1. simply import the script into PyMOL: 
         
         run /home/userName/path/toscript/pyProgram.py
         

  2. Then, just type the name of the command: _doSimpleThing_ and pass any needed arguments.



That's it. Your script can, through Python, import any modules you need and also edit modify objects in PyMOL. 

## Getting PyMOL Data into your Script

To get PyMOL data into your script you will need to somehow get access to the PyMOL objects and pull out the data. For example, if you want the atomic coordinates of a selection of alpha carbon atoms your Python function may do something like this (see also [iterate_state](/index.php/Iterate_state "Iterate state")): 
    
    
    # Import PyMOL's stored module.  This will allow us with a 
    # way to pull out the PyMOL data and modify it in our script.
    # See below.
    from pymol import stored
    
    def functionName( userSelection ):
        # this array will be used to hold the coordinates.  It
        # has access to PyMOL objects and, we have access to it.
        stored.alphaCarbons = []
    
        # let's just get the alpha carbons, so make the
        # selection just for them
        userSelection = userSelection + " and n. CA"
    
        # iterate over state 1, or the userSelection -- this just means
        # for each item in the selection do what the next parameter says.
        # And, that is to append the (x,y,z) coordinates to the stored.alphaCarbon
        # array.
        cmd.iterate_state(1, selector.process(userSelection), "stored.alphaCarbons.append([x,y,z])")
    
        # stored.alphaCarbons now has the data you want.
    
        ... do something to your coordinates ...
    

### Getting Data From your Script into PyMOL

Usually this step is easier. To get your data into PyMOL, it's usually through modifying some object, rotating a molecule, for example. To do that, you can use the [alter](/index.php/Alter "Alter") or [alter_state](/index.php?title=Alter_state&action=edit&redlink=1 "Alter state \(page does not exist\)") commands. Let's say for example, that we have translated the molecular coordinates from the last example by some vector (we moved the alpha carbons). Now, we want to make the change and see it in PyMOL. To write the coordinates back we do: 
    
    
    # we need to know which PyMOL object to modify.  There could be many molecules and objects
    # in the session, and we don't want to ruin them.  The following line, gets the object
    # name from PyMOL
    objName = cmd.identify(sel2,1)[0][0]
    
    # Now, we alter each (x,y,z) array for the object, by popping out the values
    # in stored.alphaCarbons.  PyMOL should now reflect the changed coordinates.
    cmd.alter_state(1,objName,"(x,y,z)=stored.alphaCarbons.pop(0)")
    

## Example

Here's a script I wrote for [cealign](/index.php/Cealign "Cealign"). It takes two selections **of equal length** and computes the optimal overlap, and aligns them. See [Kabsch](/index.php/Kabsch "Kabsch") for the original code. Because this tutorial is for scripting and not optimal superposition, the original comments have been removed. 
    
    
    def optAlign( sel1, sel2 ):
            """
            @param sel1: First PyMol selection with N-atoms
            @param sel2: Second PyMol selection with N-atoms
            """
    
            # make the lists for holding coordinates
            # partial lists
            stored.sel1 = []
            stored.sel2 = []
            # full lists
            stored.mol1 = []
            stored.mol2 = []
    
            # -- CUT HERE
            sel1 = sel1 + " and N. CA"
            sel2 = sel2 + " and N. CA"
            # -- CUT HERE
    
            # This gets the coordinates from the PyMOL objects
            cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
            cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
    
            # ...begin math that does stuff to the coordinates...
            mol1 = cmd.identify(sel1,1)[0][0]
            mol2 = cmd.identify(sel2,1)[0][0]
            cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
            cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
            assert( len(stored.sel1) == len(stored.sel2))
            L = len(stored.sel1)
            assert( L > 0 )
            COM1 = numpy.sum(stored.sel1,axis=0) / float(L)
            COM2 = numpy.sum(stored.sel2,axis=0) / float(L)
            stored.sel1 = stored.sel1 - COM1
            stored.sel2 = stored.sel2 - COM2
            E0 = numpy.sum( numpy.sum(stored.sel1 * stored.sel1,axis=0),axis=0) + numpy.sum( numpy.sum(stored.sel2 * stored.sel2,axis=0)
    ,axis=0)
            reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
            if reflect == -1.0:
                    S[-1] = -S[-1]
                    V[:,-1] = -V[:,-1]
            RMSD = E0 - (2.0 * sum(S))
            RMSD = numpy.sqrt(abs(RMSD / L))
            U = numpy.dot(V, Wt)
            # ...end math that does stuff to the coordinates...
    
            # update the _array_ of coordinates; not PyMOL the coords in the PyMOL object
            stored.sel2 = numpy.dot((stored.mol2 - COM2), U) + COM1
            stored.sel2 = stored.sel2.tolist()
    
            # This updates PyMOL.  It is removing the elements in 
            # stored.sel2 and putting them into the (x,y,z) coordinates
            # of mol2.
            cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
    
            print "RMSD=%f" % RMSD
    
            cmd.orient(sel1 + " and " + sel2)
    
    # The extend command makes this runnable as a command, from PyMOL.
    cmd.extend("optAlign", optAlign)
    

## Basic Script Body

Want an easy block of working code to start your function from? Just copy/paste the following into your Python editor and get going! 
    
    
    #
    # -- basicCodeBlock.py
    #
    from pymol import cmd, stored
    
    def yourFunction( arg1, arg2 ):
        '''
    DESCRIPTION
    
        Brief description what this function does goes here
        '''
        #
        # Your code goes here
        #
        print "Hello, PyMOLers"
        print "You passed in %s and %s" % (arg1, arg2)
        print "I will return them to you in a list.  Here you go."
        return (arg1, arg2)
    
    cmd.extend( "yourFunction", yourFunction );
    

# Python Version Supported by PyMOL

Scripts used within PyMOL can only be written using the current version of Python that is supported by your version of PyMOL. To determine which version of Python you can use, type the following command into PyMOL: 
    
    
    print sys.version
    

Note that this version of Python is not necessarily related to the version that you may have installed on your system. 

This command can also be used to ensure that code you are distributing can be supported by the user's system. 

Retrieved from "[https://pymolwiki.org/index.php?title=Simple_Scripting&oldid=9103](https://pymolwiki.org/index.php?title=Simple_Scripting&oldid=9103)"


---

## Single-word Selectors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/f/fb/Solvent-and-chain-A.png)](/index.php/File:Solvent-and-chain-A.png)

[](/index.php/File:Solvent-and-chain-A.png "Enlarge")

show spheres, solvent and chain A

[![](/images/0/0e/Solvent-or-chain-A.png)](/index.php/File:Solvent-or-chain-A.png)

[](/index.php/File:Solvent-or-chain-A.png "Enlarge")

show spheres, solvent or chain A

PyMOL's selection language allows to select atoms based on identifiers and properties. Many commands (like [color](/index.php/Color "Color"), [show](/index.php/Show "Show"), etc.) take an atom selection argument to only operate on a subset of all atoms in the scene. Example: 
    
    
    PyMOL>show spheres, solvent and chain A
    

Selections can be made more precise or inclusive by combining them with logical operators, including the boolean **and** , **or** , and **not**. The boolean **and** selects only those items that have both (or all) of the named properties, and the boolean **or** selects items that have either (or any) of them. 

## Contents

  * 1 Selection Operator/Modifier Table
  * 2 Comparison of distance operators
  * 3 Language Properties
  * 4 Examples
  * 5 See Also



## Selection Operator/Modifier Table

Selection operators and modifiers are listed below. The dummy variables _s1_ and _s2_ stand for selection-expressions such as "chain a" or "hydro." 

Operator  | Aliases  | Description   
---|---|---  
Generic   
all | *  | All atoms currently loaded into PyMOL   
none |  | Empty selection   
enabled |  | Atoms from enabled objects   
Named selections   
sele |  | Named selection or object "sele", but only if it doesn't conflict with the name of another operator   
%sele |  | Named selection or object "sele" Recommended, avoids ambiguity  
?sele |  | Named selection or object "sele", or empty selection if "sele" doesn't exist   
Logical   
not S1 | !  | Inverts selection   
S1 and S2 | & | Atoms included in both S1 and S2   
S1 or S2 | |  | Atoms included in either S1 or S2   
S1 S2 |  | implicit **or**  
S1 and (S2 or S3) |  | Parentheses for evaluation order control   
first S1 |  | First atom in S1 (single atom selection)   
last S1 |  | Last atom in S1 (single atom selection)   
Identifiers (see also [Selection Macros](/index.php/Selection_Macros "Selection Macros"))  
model 1ubq | m.  | Atoms from object "1ubq"   
chain C | c.  | Chain identifier "C"   
segi S | s.  | Segment identifier "S" (**label_asym_id** from mmCIF)   
resn ALA | r.  | Residue name "ALA"   
resi 100-200 | i.  | Residue identifier between 100 and 200   
name CA | n.  | Atom name "CA"   
alt A |  | Alternate location "A"   
index 123 | idx.  | Internal per-object atom index (changes with [sorting](/index.php/Sort "Sort"))   
id 123 |  | ID column from PDB file   
rank 123 |  | Per-object atom index at load time (see also [retain_order](/index.php/Retain_order "Retain order"))   
pepseq ACDEF | ps.  | Protein residue sequence with one-letter code "ACDEF" (see also [FindSeq](/index.php/FindSeq "FindSeq"))   
label "Hello World" |  | Atoms with label "Hello World" _(new in PyMOL 1.9)_  
Identifier matching   
S1 in S2 |  | Atoms in S1 whose identifiers _name, resi, resn, chain_ and _segi_ **all** match atoms in S2   
S1 like S2 |  | Atoms in S1 whose identifiers _name_ and _resi_ match atoms in S2   
Entity expansion   
Important: All "by"-operators have a **weak priority** , so (byres S1 or S2) is actually identical to (byres (S1 or S2)) and **not** to ((byres S1) or S2)   
byobject S1 |  | Expands S1 to complete objects   
bysegi S1 | bs.  | Expands S1 to complete segments   
bychain S1 | bc.  | Expands S1 to complete chains   
byres S1 | br.  | Expands S1 to complete residues   
bycalpha S1 | bca.  | CA atoms of residues with at least one atom in S1   
bymolecule S1 | bm.  | Expands S1 to complete molecules (connected with bonds)   
byfragment S1 | bf.  |   
byring S1 |  | All rings of size ≤ 7 which have at least one atom in S1 _(new in PyMOL 1.8.2)_  
bycell S1 |  | Expands selection to unit cell   
Bond expansion   
bound_to S1 | bto.  | Atoms directly bonded to S1, may include S1   
neighbor S1 | nbr.  | Atoms directly bonded to S1, excludes S1   
S1 extend 3 | xt.  | Expands S1 by 3 bonds connected to atoms in S1   
Proximity (see also comparison of distance operators)  
S1 within 12.3 of S2 | w.  | Atoms in S1 that are within 12.3 Angstroms of any atom in S2   
S1 around 12.3 | a.  | Atoms with centers within 12.3 Angstroms of the center of any atom in S1   
S1 expand 12.3 | x.  | Expands S1 by atoms within 12.3 Angstroms of the center of any atom in S1   
S1 gap 1.2 |  | Atoms whose VDW radii are separated from the VDW radii of S1 by a minimum of 1.2 Angstroms.   
S1 near_to 12.3 of S2 | nto.  | Same as _within_ , but excludes S2 from the selection (and thus is identical to `S1 and S2 around 12.3`)   
S1 beyond 12.3 of S2 | be.  | Atoms in S1 that are at least 12.3 Anstroms away from S2   
Properties   
partial_charge < 1.2 | pc.  |   
formal_charge = 1 | fc.  |   
b < 100.0 |  | B-factor less than 100.0   
q < 1.0 |  | Occupancy less than 1.0   
ss H+S |  | Atoms with secondary structure H (helix) or S (sheet)   
elem C | e.  | Atoms of element C (carbon)   
p.foo = 12 |  |   
p.foo < 12.3 |  |   
p.foo in 12+34 |  |   
stereo R |  | Chiral R/S stereo center with label R _(only[Incentive PyMOL 1.4-1.8](https://pymol.org/d/media:stereochemistry))_  
Flags   
bonded |  | Atoms which have at least one bond   
protected |  | see [protect](/index.php/Protect "Protect")  
fixed | fxd.  | see [flag](/index.php/Flag "Flag")  
restrained | rst.  | see [flag](/index.php/Flag "Flag")  
masked | msk.  | see [mask](/index.php/Mask "Mask")  
flag 25 | f.  | Atoms with flag 25, see [flag](/index.php/Flag "Flag")  
Chemical classes   
organic | org.  | Non-polymer organic compounds (e.g. ligands, buffers)   
inorganic | ino.  | Non-polymer inorganic atoms/ions   
solvent | sol.  | Water molecules   
polymer | pol.  | Protein or Nucleic Acid   
polymer.protein |  | Protein _(New in PyMOL 2.1)_  
polymer.nucleic |  | Nucleic Acid _(New in PyMOL 2.1)_  
guide |  | Protein CA and nucleic acid C4*/C4'   
hetatm |  | Atoms loaded from PDB HETATM records   
hydrogens | h.  | Hydrogen atoms   
backbone | bb.  | Polymer backbone atoms _(new in PyMOL 1.6.1)_  
sidechain | sc.  | Polymer non-backbone atoms _(new in PyMOL 1.6.1)_  
metals |  | Metal atoms _(new in PyMOL 1.6.1)_  
donors | don.  | Hydrogen bond donor atoms   
acceptors | acc.  | Hydrogen bond acceptor atoms   
Style   
visible | v.  | Atoms in enabled objects with at least one visible representation   
rep cartoon |  | Atoms with cartoon representation   
color blue |  | Atoms with atom-color blue (by color index)   
cartoon_color blue |  | Atoms with atom-level cartoon_color setting (by color index)   
ribbon_color blue |  | Atoms with atom-level ribbon_color setting (by color index)   
Non molecular   
center |  | Pseudo-atom at the center of the scene   
origin |  | Pseudo-atom at the origin of rotation   
Coordinates   
state 123 |  | Atoms with coordinates in state 123   
present | pr.  | Atoms with coordinates in the current state   
x < 12.3 |  | Atoms with model-space x coordinate less than 12.3   
y < 12.3 |  | Atoms with model-space y coordinate less than 12.3   
z > 12.3 |  | Atoms with model-space z coordinate greater than 12.3   
Atom typing   
text_type TT | tt.  | _Auto-assigned in[Incentive PyMOL 1.4-1.8](https://pymol.org/d/media:atomtyping))_  
numeric_type 123 | nt.  |   
  
## Comparison of distance operators

There are serveral very similar operators that select by pairwise atom distances. The following table lists the details how they differ. 

**Syntax 1** : _s1_ operator X of _s2_  
**Syntax 2** : _s1_ and (_s2_ operator X) 

operator | distance is ... | measured from | includes s2 | syntax | notes   
---|---|---|---|---|---  
near_to | ≤ X | center | never | 1 | equivalent to "around"   
within | ≤ X | center | if matches s1 | 1 |   
beyond | > X | center | never | 1 |   
gap | > X | center+vdw | never | 2 |   
around | ≤ X | center | never | 2 | equivalent to "near_to"   
expand | ≤ X | center | always | 2 |   
  
## Language Properties

  * names and keywords are case-insensitive unless [ignore_case](/index.php/Ignore_case "Ignore case") is set
  * names and keywords can be abbreviated to non-ambiguous prefixes



**Best practice recommendation:** Only write case-sensitive, non-abbreviated selection expressions. That way your scripts will be robust against run-time configuration and future changes to the language (like addition of new keywords). 

## Examples

Logical selections can be combined. For example, you might select atoms that are part of chain a, but not residue number 125: 
    
    
    # selects atoms that are part of chain A, but not residue number 125.
    select chain A and (not resi 125)
    
    # The following two selections are equivalent, 
    select (name CB or name CG1 or name CG2) and chain A
    
    # select c-beta's, c-gamma-1's and c-gamma-2's 
    # that are in chain A.
    select name CB+CG1+CG2 and chain A
    
    # select all residues within 5 Ang. or any organic small molecules
    select br. all within 5 of organic
    
    # select helices
    select ss 'H'
    
    # select anything shown as a line
    select rep lines
    
    # select all residues with a b-factor less than 20, within 3 angstroms of any water
    select br. b<20 & (all within 3 of resn HOH)
    
    # select anything colored blue
    select color blue
    
    # select the 1st arginine
    select first resn ARG
    
    # select 1foo's segment G's chain X's residue 444's alpha carbon
    select 1foo/G/X/444/CA
    # same thing
    select 1foo and segi G and c. X and i. 444 and n. CA
    
    # select the entire object that residue 23's beta caron is in:
    select bo. i. 23 and n. CA
    
    # select the molecule that chain C is in
    select bm. c. C
    

Like the results of groups of arithmetic operations, the results of groups of logical operations depend on which operation is performed first. They have an order of precedence. To ensure that the operations are performed in the order you have in mind, use parentheses: 
    
    
    byres ((chain A or (chain B and (not resi 125))) around 5)
    

PyMOL will expand its logical selection out from the innermost parentheses. 

## See Also

  * [select](/index.php/Select "Select")
  * [Selection Macros](/index.php/Selection_Macros "Selection Macros")
  * [Property_Selectors](/index.php/Property_Selectors "Property Selectors")
  * [Selection Language Comparison](/index.php/Selection_Language_Comparison "Selection Language Comparison") with other modelling applications
  * [Identify](/index.php/Identify "Identify")



Retrieved from "[https://pymolwiki.org/index.php?title=Selection_Algebra&oldid=12730](https://pymolwiki.org/index.php?title=Selection_Algebra&oldid=12730)"


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

## States

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Introduction to States

Working with States in PyMOL is a very common task. A State is one particular conformation of an object. For example, one could load an NMR ensemble and set [all_states](/index.php/All_states "All states") or use the command [Split_States](/index.php/Split_States "Split States") to see all entries in the NMR file. Another example could be the set of states from a molecular dynamic (MD) simulation. Each conformer in the MD ensemble could be viewed as a state in PyMOL. 

  * [![All states in an ensemble are shown.](/images/7/76/All_states_on.png)](/index.php/File:All_states_on.png "All states in an ensemble are shown.")

All states in an ensemble are shown. 

  * [![The states are hidden.](/images/6/6f/All_states_off.png)](/index.php/File:All_states_off.png "The states are hidden.")

The states are hidden. 




If you are making a movie of a static coordinate set (such as a single crystal structure) then you have only one state. All objects in PyMOL can potentially consist of multiple states. For movies, see [Frames](/index.php/Frame "Frame"). 

# Using States in PyMOL

  * Upon loading an object you can separate each state into its own object with [Split_States](/index.php/Split_States "Split States").
  * States can be [colored](/index.php/Color "Color") individually.
  * One can show all states with the [all_states](/index.php/All_states "All states") setting.
  * One can iterate over states using [Iterate_State](/index.php/Iterate_State "Iterate State") or change properties with states using the [Alter_State](/index.php/Alter_State "Alter State") function.



# See Also

[PyMOL States Related Pages](/index.php/Category:States "Category:States")

Retrieved from "[https://pymolwiki.org/index.php?title=States&oldid=6528](https://pymolwiki.org/index.php?title=States&oldid=6528)"


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

