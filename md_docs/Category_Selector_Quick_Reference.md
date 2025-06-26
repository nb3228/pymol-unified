# Category: Selector Quick Reference

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

