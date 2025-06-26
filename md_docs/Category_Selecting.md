# Category: Selecting

## Ignore case

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The ignore_case setting (default: **on** , _except in PyMOL 1.8.0.0 - 1.8.0.4_) controls whether PyMOL does case sensitive matching of atomic identifiers and selection operators in the selection language. Most notably, it affects whether chain identifiers are matched case sensitive, which becomes relevant when using upper and lower case chain identifiers in a structure with more than 26 chains. 

The default value was changed to **off** in PyMOL 1.8.0.0. However, due to undesired side effects, the **on** default was restored in 1.8.0.5. 

The next PyMOL version (expect 1.8.2) will introduce a new [ignore_case_chain](/index.php/Ignore_case_chain "Ignore case chain") setting to address the issue of mixed case chain identifiers. 

See also: <https://sourceforge.net/p/pymol/mailman/message/34815599/>

## Example

Load 1a00 which has chains A, B, C, D 
    
    
    PyMOL>fetch 1a00, async=0
    

1) case insensitive selection language (default for PyMOL <= 1.7.6 and >= 1.8.0.5) 
    
    
    PyMOL>set ignore_case
     Setting: ignore_case set to on.
    PyMOL>count_atoms chain A
     count_atoms: 1164 atoms
    PyMOL>count_atoms chain a
     count_atoms: 1164 atoms
    

2) case sensitive selection language (default for PyMOL 1.8.0.0 - 1.8.0.4) 
    
    
    PyMOL>set ignore_case, off
     Setting: ignore_case set to off.
    PyMOL>count_atoms chain A
     count_atoms: 1164 atoms
    PyMOL>count_atoms chain a
     count_atoms: 0 atoms
    

## Best Practice for Scripting

When writing scripts or plugins, all selection expressions should be strict about case (e.g. "name CA" and not "name ca") to not depend on the users setting of ignore_case. 

Retrieved from "[https://pymolwiki.org/index.php?title=Ignore_case&oldid=12289](https://pymolwiki.org/index.php?title=Ignore_case&oldid=12289)"


---

## Ignore case chain

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The ignore_case_chain setting (default: **off**) controls whether the **chain** and **segi** [selection operators](/index.php/Selection_Algebra "Selection Algebra") are case sensitive. The default is to be case sensitive, so "chain A" is different from "chain a". 

_New In PyMOL 1.8.2 (see[ignore_case](/index.php/Ignore_case "Ignore case") for previous versions)_

## Example

The PDB structure **3tcx** has 28 chains, including "A" and "a". 
    
    
    PyMOL>fetch 3tcx, async=0
    PyMOL>count_atoms chain A
     count_atoms: 641 atoms
    PyMOL>set ignore_case_chain, on
    PyMOL>count_atoms chain A
     count_atoms: 1282 atoms
    

## See Also

  * [ignore_case](/index.php/Ignore_case "Ignore case")
  * [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")



Retrieved from "[https://pymolwiki.org/index.php?title=Ignore_case_chain&oldid=12695](https://pymolwiki.org/index.php?title=Ignore_case_chain&oldid=12695)"


---

## Named Atom Selections

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Atom selections can be named for repeated use by using the select command: 
    
    
    SYNTAX
          
       select  selection-name, selection-expression   
                                               # The  selection-name and  
                                               # the  selection-expression 
                                               # are both arguments to select 
                                               # so they are separated by a comma. 
                                               
    EXAMPLE
    
       PyMOL> select bb, name c+o+n+ca   # Create an atom selection named "bb"
                                         # including all atoms named 
                                         # "C","O","N", or "CA";
       PyMOL> color red, bb              # color the selection red,
       PyMOL> hide lines, bb             # hide the line representation,
       PyMOL> show sticks, bb            # show it using the stick representation,
       PyMOL> zoom bb                    # and zoom in on it.
    

In this case, the selection-expression is the property selector name, which takes the list of identifiers ca+c+n+o to complete the specification. Property selectors and their identifiers are discussed below. 

Named atom selections appear in the PyMOL names list in the control panel. They are distinguished from objects by a surrounding set of parentheses. The control panel options available under the diamond menu differ between objects and atom-selections, because objects and named selections play slightly different roles in PyMOL. Named selections are pointers to subsets of data that are found under an object name. After an object is deleted, the data are no longer available, unless you reload the object. Any named selections that refer to atoms in that object will no longer work. But when named selections are deleted, the data are still available under the object name. Disabling objects eliminates them from the viewer, but disabling named-selections just turns off the pink dots that highlight them in the viewer. 

Atom-selections, named or not, can span multiple objects: 
    
    
    EXAMPLE
    
       PyMOL> load $PYMOL_PATH/test/dat/fc.pdb                                      
       PyMOL> load $PYMOL_PATH/test/dat/pept.pdb
    
       PyMOL> select alpha_c, name ca      # The named selection "alpha_c" 
                                             # is created -- it includes atoms
                                             # in both "fc" and "pept" objects.
       PyMOL> color red, name ca             # "CA" atoms in both objects
                                             # are colored red.                                      
    

Named selections will continue working after you have made changes to a molecular structure: 
    
    
    EXAMPLE
    
       PyMOL> load $PYMOL_PATH/test/dat/pept.pdb                                      
       PyMOL> select bb, name c+o+n+ca     # The named selection "bb" 
                                           # is created.        
                                           
       PyMOL> count_atoms bb               # PyMOL counts 52 atoms in "bb."   
       
       PyMOL> remove resi 5                # All atoms in residue 5 are removed 
                                           # from the object "pept."  
                                           
       PyMOL> count_atoms bb               # Now PyMOL counts
                                           # the remaining 48 atoms in "bb."                                                                                                                                            
    

Named selections are static. Only atoms that exist at the time the selection is defined are included in the selection, even if atoms which are loaded subsequently fall within the selection criterion: 
    
    
    EXAMPLE
    
       PyMOL> load $PYMOL_PATH/test/dat/pept.pdb      
       
       PyMOL> select static_demo, pept    # The named selection "static_demo" 
                                          # is created to reference all atoms.  
                                          
       PyMOL> count_atoms static_demo     # PyMOL counts 107 atoms 
                                          # in "static_demo."   
                                          
       PyMOL> h_add                       # PyMol adds hydrogens in
                                          # the appropriate places
                                          
       PyMOL> count_atoms static_demo     # PyMOL still counts 107 atoms
                                          # in "static_demo,"
       PyMOL> count_atoms                 # even though it counts 200 atoms
                                          # in "pept." 
    

Named selections can also be used in subsequent atom selections: 
    
    
    EXAMPLE
    
       PyMOL> select bb, name c+o+n+ca   # An atom selection named "bb"
                                         # is made, consisting of all 
                                         # atoms named "C","O","N", or "CA."
                                               
       PyMOL> select c_beta_bb, bb or name cb   
                                         # An atom selection named "c_beta_bb"
                                         # is made, consisting of 
                                         # all atoms named "C", "O", "N", "CA" or "CB."  
    

Note that the word "or" is used to select all atoms in the two groups, "bb" and "cb." The word "and" would have selected no atoms because it is interpreted in its boolean logical sense, not its natural language sense. See the subsection on "Selection Algebra" below. 

Retrieved from "[https://pymolwiki.org/index.php?title=Named_Atom_Selections&oldid=4327](https://pymolwiki.org/index.php?title=Named_Atom_Selections&oldid=4327)"


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

## Selection Algebra

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

## Selection-expressions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Selection-expressions 

Selection-expressions stand for lists of atoms in arguments that are subject to PyMOL commands. You can name the selections to facilitate their re-use, or you can specify them anonymously (without names). Object and selection names may include the upper or lower case characters A/a to Z/z, numerals 0 to 9, and the underscore character (_). Characters to avoid include: 

! @ # $ % ^ &* ( ) ' " [ ] { } \ | ~ ` <> . ? / 

Selection-expressions describe the class of atoms you are referencing. Most of them require identifiers to complete the specification. For example, the selector resi references biopolymer residues by sequence number, and the identifier gives the number. The selector name references atoms according to the names described in the PDB, and the identifier gives the name (ca for alpha carbons, cb for beta carbons, etc). A handful of selection-expressions don't require identifiers, but most do. 

PyMOL uses several logical operators to increase the generality or specificity of selection-expressions. Logical combinations of selectors can get complex, so PyMOL accepts short forms and macros that express them with a minimum of keystrokes. This section describes named-selections, and then gives the syntax for making selections in a progression from simple one-word selectors to complex combinations of selectors, using macros and short forms. 

Retrieved from "[https://pymolwiki.org/index.php?title=Selection-expressions&oldid=8445](https://pymolwiki.org/index.php?title=Selection-expressions&oldid=8445)"


---

## Selection Language Comparison

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Many molecular modelling applications have an atom selection language. This page summarizes and compares languages from different tools. 

## Selection Operators

| [PyMOL](/index.php/Selection_Algebra "Selection Algebra") | [Phenix](http://www.phenix-online.org/documentation/reference/atom_selections.html) | [BALL](http://ball-docs.bioinf.uni-sb.de/BALLView/V1.4.0/molecularControl.html#regular_expressions) | [VMD](http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.3/ug/node89.html) | [CNS](http://cns-online.org/v1.3/tutorial/language/atom_sele_basic/text.html) | [RasMol](http://rasmol.org/doc/rasmol.html#chexprs) | [OpenStructure](http://www.openstructure.org/docs/dev/mol/base/query/) | [ProDy](http://prody.csb.pitt.edu/manual/reference/atomic/select.html) | [Maestro](http://gohom.win/ManualHom/Schrodinger/Schrodinger_2015-2_docs/maestro/help_Maestro/misc/asl.html)  
---|---|---|---|---|---|---|---|---|---  
Segment  | segi A  | segid A  |  | segname A  | segid A  |  |  | segment A  |   
Chain  | chain A  | chain A  | chain(A)  | chain A  |  | *A  | cname=A  | chain A  | chain.name A   
Residue ID  | resi 3-8  | resseq 3:8  | residueID(3)  | resid 3 to 8  | residue 3:8  | 3-8  | rnum=3:8  | resnum 3 to 8  | res. 3 - 8   
Residue Name  | resn HIS  | resname HIS  | residue(HIS)  | resname HIS  | resname HIS  | HIS  | rname=HIS  | resname HIS  | res.pt HIS   
Atom Name  | name ND1  | name ND1  | name(ND1)  | name ND1  | name ND1  | *.ND1  | aname=ND1  | name ND1  | atom.pt ND1   
Logical And  | and  | and  | AND  | and  | and  | and  | and  | and  | and   
Logical Or  | or  | or  | OR  | or  | or  | or  | or  | or  | or   
Logical Not  | not  | not  | !  | not  | not  | not  |  | not  | not   
Implicit operator  | or  |  |  | and  |  |  |  |  |   
Atom Symbol  | elem H  | element H  | element(H)  |  | chemical H  | elemno == 1  | ele=H  | element H  | atom.ele H   
Hydrogens  | hydro  |  | element(H)  | hydrogen  | hydrogen  | hydrogen  |  | hydrogen  | atom.ato 1   
Secondary Structure  | ss H  |  | secondaryStruct(...)  | structure H  |  | helix or sheet  | rtype=H  | secondary H  | res.sec h   
Solvent  | solvent  | water  | solvent()  | water  |  | solvent  | water  | water  | water   
Backbone  | backbone (since 1.7)  | backbone  | backbone()  | backbone  |  | backbone  |  | backbone  | backbone   
Sidechain  | sidechain (since 1.7)  | sidechain  |  | sidechain  |  | sidechain  |  | sidechain  | sidechain   
Within Radius  | within 6.0 of (chain A)  | within(6.0, chain A)  |  | within 6.0 of (chain A)  | (chain A) around 6.0  | within(6.0, *A)  | 6.0 <> [cname=A]  |  | within 6.0 (chain. A)   
Protein+Nucleic  | polymer  |  |  |  |  |  |  |  |   
Protein  | polymer.protein  |  |  | protein  |  | protein  | peptide  | protein  | protein   
Peptide Sequence  | pepseq ACDEF  |  |  | protein sequence ACDEF  |  |  |  |  |   
DNA  |  | dna  |  |  |  |  |  |  |   
RNA  |  | rna  |  |  |  |  |  |  |   
Nucleic Acids  | polymer.nucleic  | nucleotide  |  | nucleic  |  | nucleic  |  | nucleic  | nucleic_acids   
B-factor  | b < 50  |  |  |  |  | temperature < 50  | abfac<50  | beta < 50  |   
  
## Selection Macros

| chain A and resi 32 and name CA  | Comment   
---|---|---  
[PyMOL](/index.php/Selection_Macros "Selection Macros") | A/32/CA  |   
Modeller  | mdl.atoms['CA:32:A']  | mdl is model instance   
RasMol  | *32A.CA  | full macro example: SER32A.CA;B/2 (Serine, ..., Alt B, Model 2))   
[Chimera](http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/atom_spec.html) | :32.A@CA  | not tested   
  
## See Also

  * [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")
  * [Selection Macros](/index.php/Selection_Macros "Selection Macros")



Retrieved from "[https://pymolwiki.org/index.php?title=Selection_Language_Comparison&oldid=12905](https://pymolwiki.org/index.php?title=Selection_Language_Comparison&oldid=12905)"


---

## Selection Macros

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/a/ac/Right-click-macro.png)](/index.php/File:Right-click-macro.png)

The atom-right-click context menu title shows a selection macro

Selection Macros allow to represent a long [atom selection](/index.php/Selection_Algebra "Selection Algebra") phrase such as 
    
    
    PyMOL> select pept and segi lig and chain B and resi 142 and name CA
    

in a more compact form: 
    
    
    PyMOL> select /pept/lig/B/142/CA
    

## Contents

  * 1 Syntax and Semantics
  * 2 Details
    * 2.1 Macros Beginning With a Slash
    * 2.2 Macros Not Beginning With a Slash
    * 2.3 Omitting Fields in a Macro
  * 3 See Also



## Syntax and Semantics

An atom selection macro uses slashes to define fields corresponding to identifiers. The macro is used to select atoms using the boolean "and," that is, the selected atoms must have all the matching identifiers: 
    
    
    /object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
    

  * must contain **at least one slash**
  * no spaces allowed
  * **empty** fields are interpreted as **wildcards**
  * starting with a slash: 
    * **yes** : fields from the **right** can be omitted
    * **no** : fields from the **left** can be omitted



## Details

Macros come in two flavors: those that begin with a slash and those that don't. The presence or absence of a slash at the beginning of the macro determines how it is interpreted. If the macro begins with a slash, PyMOL expects to find the fields starting from the top of the hierarchy: the first field to the right of the slash is interpreted as an object-name; the second field as an identifier to segi; the third as an identifier to chain, and so on. It may take any of the following forms: 

### Macros Beginning With a Slash
    
    
      
       /object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
       /object-name/segi-identifier/chain-identifier/resi-identifier
       /object-name/segi-identifier/chain-identifier
       /object-name/segi-identifier
       /object-name
    
    
    EXAMPLES
       PyMOL> zoom /pept
       PyMOL> show spheres, /pept/lig/
       PyMOL> show cartoon, /pept/lig/A
       PyMOL> color pink, /pept/lig/A/10
       PyMOL> color yellow, /pept/lig/A/10/CA
    

### Macros Not Beginning With a Slash

If the macro does not begin with a slash, it is interpreted differently. In this case, PyMOL expects to find the fields ending with the bottom of the hierarchy. Macros that don't start with a slash may take the following forms: 
    
    
                                                 resi-identifier/name-identifier
                                chain-identifier/resi-identifier/name-identifier
                segi-identifier/chain-identifier/resi-identifier/name-identifier
    object-name/segi-identifier/chain-identifier/resi-identifier/name-identifier
    
    
    EXAMPLES
       PyMOL> zoom 10/CB
       PyMOL> show spheres, A/10-12/CB
       PyMOL> show cartoon, lig/B/6+8/C+O
       PyMOL> color pink, pept/enz/C/3/N
    

### Omitting Fields in a Macro

You can also omit fields between slashes. Omitted fields will be interpreted as wildcards, as in the following forms: 
    
    
       
       resi-identifier/
       resi-identifier/name-identifier
       chain-identifier//
       object-name//chain-identifier                
    
       
    EXAMPLES
       PyMOL> zoom 142/                  # Residue 142 fills the viewer. 
       
       PyMOL> show spheres, 156/CA       # The alpha carbon of residue 156
                                         # is shown as a sphere     
                                         
       PyMOL> show cartoon, A//          # Chain "A" is shown as a cartoon.  
       
       PyMOL> color pink, /pept//B       # Chain "B" in object "pept"
                                         # is colored pink.
    

## See Also

  * [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")



Retrieved from "[https://pymolwiki.org/index.php?title=Selection_Macros&oldid=12699](https://pymolwiki.org/index.php?title=Selection_Macros&oldid=12699)"


---

## Property Selectors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Most of this page has been merged into [Selection Algebra](/index.php/Selection_Algebra "Selection Algebra")

## Contents

  * 1 Property Selectors
    * 1.1 General
    * 1.2 Property Selector Table
    * 1.3 Numeric Selector Table
  * 2 User Notes
    * 2.1 Select by PDB Assigned Atom Number
    * 2.2 Selecting via a Negative Identifier
    * 2.3 Alternate Locations
      * 2.3.1 Select All Alternate Locations
      * 2.3.2 Select All Alternate Locations Except the A-Locations
  * 3 See Also



## Property Selectors

PyMOL reads data files written in PDB, MOL/SDF, Macromodel, ChemPy Model, and Tinker XYZ formats. Some of the data fields in these formats allow PyMOL to assign properties to atoms. You can group and select atoms according to these properties using property selectors and identifiers: the selectors correspond to the fields in the data files, and the identifiers correspond to the target words to match, or the target numbers to compare. 

### General

The items in a list of identifiers are separated by plus signs (+) only. Do not add spaces within a list of identifiers. The selector **resi** takes (+)-separated lists of identifiers, as in 
    
    
    PyMOL> select nterm, resi 1+2+3
    

or, alternatively, it may take a range given with a dash 
    
    
    PyMOL> select nterm, resi 1-3
    

The identifier for a blank field in an input file is an empty pair of quotes 
    
    
    PyMOL> select unstruct, ss ""   # A named selection is created
                                    # to contain all atoms that are not assigned 
                                    # a secondary structure.
    

### Property Selector Table

Most property selectors select matches to their identifiers 

Matching Property Selector |  Short Form Selector |  Identifier and Example  
---|---|---  
symbol|  e.| ___chemical-symbol-list_   
list of 1- or 2-letter chemical symbols from the periodic table   

    
    
    PyMOL> select polar, symbol O+N  
  
name|  n.| ___atom-name-list_   
list of up to 4-letter codes for atoms in proteins or nucleic acids   

    
    
    PyMOL> select carbons, name CA+CB+CG+CD  
  
resn|  r.| ___residue-name-list_   
list of 3-letter codes for amino acids   

    
    
    PyMOL> select aas, resn ASP+GLU+ASN+GLN

or list of up to 2-letter codes for nucleic acids   

    
    
    PyMOL> select bases, resn A+G  
  
resi|  i.| ___residue-identifier-list_   
list of up to 4-digit residue numbers   

    
    
    PyMOL> select mults10, resi 1+10+100+1000

_residue-identifier-range_   

    
    
    PyMOL> select nterm, resi 1-10  
  
alt|  alt|  ___alternate-conformation-identifier-list_   
list of single letters 
    
    
    PyMOL> select altconf, alt A+""  
  
chain|  c.| ___chain-identifier-list_   
list of single letters or sometimes numbers   

    
    
    PyMOL> select firstch, chain A  
  
segi|  s.| ___segment-identifier-list_   
list of up to 4 letter identifiers 
    
    
    PyMOL> select ligand, segi lig   
  
flag|  f.| ___flag-number_   
a single integer from 0 to 31   

    
    
    PyMOL> select f1, flag 0  
  
numeric_type |  nt.|  _type-number_   
a single integer   

    
    
    PyMOL> select type1, nt. 5  
  
text_type|  tt.| ___type-string_   
a list of up to 4 letter codes 
    
    
    PyMOL> select subset, text_type HA+HC  
  
id|  id|  ___external-index-number_   
a single integer   

    
    
    PyMOL> select idno, id 23  
  
index|  idx.| ___internal-index-number_   
a single integer   

    
    
    PyMOL> select intid, index 11  
  
ss|  ss|  ___secondary-structure-type_   
list of single letters   

    
    
    PyMOL> select allstrs, ss H+S+L+""  
  
  


### Numeric Selector Table

Other _property selectors_ select by comparison to numeric identifiers 

Numeric Selector |  Short Form | Argument and Example  
---|---|---  
b|  b|  ___comparison-operator b-factor-value_   
a real number 
    
    
    PyMOL> select fuzzy, b > 10  
  
q|  q|  ___comparison-operator occupancy-value_   
a real number 
    
    
    PyMOL> select lowcharges, q <0.50  
  
formal_charge |  fc.|  _comparison-operator formal charge-value_   
an integer 
    
    
    PyMOL> select doubles, fc. = -1  
  
partial_charge |  pc.|  _comparison-operator partial charge-value_   
a real number 
    
    
    PyMOL> select hicharges, pc. > 1  
  
  


## User Notes

### Select by PDB Assigned Atom Number

To select atoms by their PDB ATOM number, use the **id** selector: 
    
    
     select Nt, id 1-30
    

makes a new selection called Nt and puts in it the first 30 atoms (assuming your PDB starts numbering at 1). 

### Selecting via a Negative Identifier

If you have a residue id < 0, eg. resi -55, then you have to escape the minus sign or renumber the atoms: 
    
    
    select i. \-55
    

### Alternate Locations

#### Select All Alternate Locations

  * To select **all** residues with alternate locations (alt, alt loc), simply do:


    
    
    select aa, not alt ""
    

Try it on 1CBN. 

#### Select All Alternate Locations Except the A-Locations

  * To select and remove all atoms with alternate locations (ALT, alt loc) that aren't the **A** record try the following. This removes all alt locations that are not the original **A** location. For example given the following pdb record,


    
    
    ATOM    949  N   ARG A 124       6.039   8.484   9.313  1.00  9.82           N
    ATOM    950  CA  ARG A 124       5.839   7.071   8.980  1.00 10.05           C
    ATOM    951  C   ARG A 124       6.009   6.809   7.493  1.00  9.18           C
    ATOM    952  O   ARG A 124       5.803   7.715   6.648  1.00  9.06           O
    ATOM    953  CB  ARG A 124       4.437   6.660   9.425  1.00 10.61           C
    ATOM    954  CG AARG A 124       4.381   6.539  10.955  0.50 12.54           C
    ATOM    955  CG BARG A 124       4.170   6.743  10.930  0.50 14.58           C
    ATOM    956  CD AARG A 124       2.996   6.418  11.561  0.50 11.86           C
    ATOM    957  CD BARG A 124       2.700   7.091  11.242  0.50 17.41           C
    ATOM    958  NE AARG A 124       2.110   7.506  11.144  0.50 13.01           N
    ATOM    959  NE BARG A 124       2.229   6.746  12.588  0.50 19.94           N
    ATOM    960  CZ AARG A 124       0.912   7.729  11.654  0.50 12.21           C
    ATOM    961  CZ BARG A 124       2.123   5.510  13.067  0.50 21.42           C
    ATOM    962  NH1AARG A 124       0.433   6.958  12.627  0.50 12.93           N
    ATOM    963  NH1BARG A 124       2.470   4.470  12.319  0.50 24.35           N
    ATOM    964  NH2AARG A 124       0.185   8.735  11.188  0.50 11.89           N
    ATOM    965  NH2BARG A 124       1.658   5.310  14.313  0.50 18.28           N
    

the command below will remove the (not A so) B-locations (lines 7, 9, 11, 13, 15, 17). 
    
    
    # select & remove all non A altlocs
    remove not (alt ''+A)
    # reset the PDB information
    alter all, alt=''
    

See also, [removeAlt](/index.php/RemoveAlt "RemoveAlt") for a script to do this for you. 

## See Also

[Single-word Selectors](/index.php/Single-word_Selectors "Single-word Selectors") [Selection Macros](/index.php/Selection_Macros "Selection Macros") [Selection_Algebra](/index.php/Selection_Algebra "Selection Algebra") [Identify](/index.php/Identify "Identify")

Retrieved from "[https://pymolwiki.org/index.php?title=Property_Selectors&oldid=12700](https://pymolwiki.org/index.php?title=Property_Selectors&oldid=12700)"


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

