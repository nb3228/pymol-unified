# Category: Editing Module

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

