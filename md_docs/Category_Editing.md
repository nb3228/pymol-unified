# Category: Editing

## Builder

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

The Builder is a PyMOL GUI menu that allows you to easily build up structures by hand from various elements like atoms, fragments, rings, amino acids, etc. You can also assign charge, secondary structure, etc. 

To access the Builder simply select the "Builder" option from the PyMOL GUI (see images). 

  * [![Builder button in the GUI](/images/4/44/Builder0.png)](/index.php/File:Builder0.png "Builder button in the GUI")

Builder button in the GUI 

  * [![Builder activated.](/images/6/6d/Builder1.png)](/index.php/File:Builder1.png "Builder activated.")

Builder activated. 




# Nucleic Acid Builder

As of **Incentive PyMOL 2.3** , building can also be done with nucleic acids. 

[![DNA Builder1.png](/images/e/e9/DNA_Builder1.png)](/index.php/File:DNA_Builder1.png)

Click the 'Builder' button in the top-right panel. 

[![DNA Builder2.png](/images/c/ce/DNA_Builder2.png)](/index.php/File:DNA_Builder2.png)

Once selecting "Nucleic Acids", you will presented options for building either DNA or RNA. For DNA, you can specify the custom DNA's form (A-DNA or B-DNA) and whether the DNA should be single- or double-stranded. 

[![DNA Builder3.png](/images/9/96/DNA_Builder3.png)](/index.php/File:DNA_Builder3.png)

When ready, select the base type for the nucleic acid by clicking on the button with its specified one-letter code (e.g. A, T, U, C, or G) and select "Create As New Object". 

[![DNA Builder4.png](/images/0/08/DNA_Builder4.png)](/index.php/File:DNA_Builder4.png)

To extend the DNA, select a terminal 5' phosphate or 3' oxygen. 

You can see the change in the nucleic acid's sequence in real-time (Display-> Sequence). For double-stranded DNA, positive residue indices are base-paired to its negative counterpart (5 base-paired with -5, 12 base-paired with -12, etc...). 

Once finished, click the "Done" button. 

[![DNA Builder5.png](/images/7/7e/DNA_Builder5.png)](/index.php/File:DNA_Builder5.png)

Here is an example finished product from the DNA builder. 

_New in Incentive PyMOL 2.5_

As of Incentive PyMOL 2.5, builder can attach nucleic acids to existing structures. Double stranded structures are automatically detected based on distance cutoff and base pairing. Additionally, a phosphate group will be automatically added to the 5' end if missing from the structure. 

The "Residue: Remove" button will remove the residue currently selected based on PyMOL's Selection Algebra. This means the entire residue of the atom selected will be removed. This works for both proteins and nucleic acids. [![Remove-Residue-Button.png](/images/a/a5/Remove-Residue-Button.png)](/index.php/File:Remove-Residue-Button.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Builder&oldid=13427](https://pymolwiki.org/index.php?title=Builder&oldid=13427)"


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

