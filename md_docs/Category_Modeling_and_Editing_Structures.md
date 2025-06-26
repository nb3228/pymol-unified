# Category: Modeling and Editing Structures

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

## Homology Modeling

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Align proteins with CA fit

If the proteins have significant homology, then you can use the align command: 
    
    
    align prot1////ca,prot2
    

which will perform a sequence alignment of prot1 against prot2, and then an optimizing fit using the CA positions. 

Retrieved from "[https://pymolwiki.org/index.php?title=Homology_Modeling&oldid=3223](https://pymolwiki.org/index.php?title=Homology_Modeling&oldid=3223)"


---

## Molecular Sculpting

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## What the heck is molecular sculpting?

Warren DeLano's reply to the question "molecular sculpting" in PyMOL mailing list.[[1]](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg00359.html): 

> Molecular sculpting works like a real-time energy minimizer, except that it isn't minimizing the energy. Instead, its just trying to return local atomic geometries (bonds, angles, chirality, planarity) to the configuration the molecules possess when they were first loaded into PyMOL. 
> 
> To actually use this feature: 
> 
>   1. Load a PDB file.
>   2. Configure the mouse for editing (_Mouse_ menu) or click in the mouse/key matrix box.
>   3. Select _Sculpting_ from the _Wizard_ menu. (or: Select _auto-sculpting_ from the _Sculpting_ menu.)
>   4. _Ctrl-middle-click_ on any atom in your protein to activate sculpting. The green part will be free to move. The cyan part will be a fixed cushion to provide context. The grey part will be excluded.
>   5. Now perform any conformational editing operation in the green region such as: 
>      * _ctrl-left-click-and-drag_ on an atom
>      * _ctrl-right-click_ on a bond, then _ctrl-left-click-and-drag_ about that bond.
> 

> 
> You can adjust the radius and cushion using the blue pop-up menus. 
> 
> Right now I'm not sure the sculpting feature is more than entertainment, but my expectation is that it will become part of PyMOL's crystallographic model building system in the future. 

## See Also

  * [sculpt_field_mask](/index.php/Sculpt_field_mask "Sculpt field mask") setting
  * [Optimize](/index.php/Optimize "Optimize") plugin



Retrieved from "[https://pymolwiki.org/index.php?title=Molecular_Sculpting&oldid=12803](https://pymolwiki.org/index.php?title=Molecular_Sculpting&oldid=12803)"


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

