# Category: Wizards

## Annotation wizard

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Here's an example of using the annotation wizard (see top-left corner): 

[![Ann.png](/images/3/30/Ann.png)](/index.php/File:Ann.png)
    
    
    import pymol                                                                                                                                                                                
                                                                                                                                                                                                
    # fetch a protein and setup the view
    cmd.fetch("1oky", async=0)                                                                                                                                                                  
    cmd.extract("ligand", "resn STU")                                                                                                                                                           
    cmd.show_as("sticks", "resn STU")                                                                                                                                                           
    cmd.show("surface", "ligand")                                                                                                                                                               
    cmd.flag("ignore", "not rep surface")                                                                                                                                                       
    cmd.set("surface_color", "grey")                                                                                                                                                            
    cmd.set("transparency", 0.3)                                                                                                                                                                
    cmd.distance("(ligand)", "(poly)", quiet=1, mode=2, label=1)                                                                                                                                
    
    # turn on the annotation wizard + prompt                                                                                                                                                                                            
    cmd.wizard("annotation")                                                                                                                                                                    
    cmd.set("wizard_prompt_mode", 1)                                                                                                                                                            
    
    pymol.session.annotation = {}                                                                                                                                                               
                                                                                                                                                                                                
    state_dict = {1: ['\\999Molecules:','  \\459staurosporine (ligand)','  \\459protein (1oky)','\\999Objects:','  \\459Hydrogen bond neighbors',]}                                             
                                                                                                                                                                                                
    pymol.session.annotation["1oky"] = state_dict
    

The following source shows you how to use the annotation wizard for multiple objects. 
    
    
    import pymol                                                                                                                          
    
    # fetch a protein and setup the view
    cmd.fetch("1oky", async=0)                                                                                                           
    cmd.extract("ligand", "resn STU")                                                                                                    
    cmd.show_as("sticks", "resn STU")                                                                                                    
    cmd.show("surface", "ligand")                                                                                                        
    cmd.flag("ignore", "not rep surface")                                                                                                
    cmd.set("surface_color", "grey")                                                                                                     
    cmd.set("transparency", 0.3)                                                                                                         
    cmd.distance("(ligand)", "(poly)", quiet=1, mode=2, label=1)                                                                         
     
    # turn on the annotation wizard + prompt                                                                                             
    cmd.wizard("annotation")                                                                                                             
    
    cmd.set("wizard_prompt_mode", 1)                                                                                                     
    pymol.session.annotation = {}                                                                                                        
    
    # using the annotation wizard for multiple objects:
    
    prot_dict = { 1: ['\\999Protein:',' \\459protein (1oky)',]}
    lig_dict  = { 1: ['\\999Ligand:',' \\459staurosporine (ligand)',]}
    dist_dict = { 1: ['\\999Objects: ',' \\459Hydrogen bond neighbors',]}
    
    pymol.session.annotation["1oky"] = prot_dict
    
    pymol.session.annotation["ligand"] = lig_dict
    
    pymol.session.annotation["dist01"] = dist_dict
    

Retrieved from "[https://pymolwiki.org/index.php?title=Annotation_wizard&oldid=10492](https://pymolwiki.org/index.php?title=Annotation_wizard&oldid=10492)"


---

## Mutagenesis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Mutagenesis

PyMol has a **Mutagenesis Wizard** to make mutagenesis very easy for the end user. 

In [rotkit](/index.php/Rotkit "Rotkit"), a function has been made to call a mutagenesis. 

As of PyMOL version [2.2](https://pymol.org/dokuwiki/?id=media:new22), users may now perform base mutations in nucleotide chains. 

## Walk-through

To mutate a residue follow these easy steps: 

  1. Load a PDB file  
  


[![Mutag1.png](/images/7/78/Mutag1.png)](/index.php/File:Mutag1.png)

  2. Under the **Wizard** menu select **Mutagenesis**  
  


[![Mutag2.png](/images/1/19/Mutag2.png)](/index.php/File:Mutag2.png)

  3. In the PyMol viewer window select a residue  
  


[![Mutag3.png](/images/1/19/Mutag3.png)](/index.php/File:Mutag3.png)

  4. Select **No Mutation** and select resultant residue  
  


[![Mutag4.png](/images/b/bb/Mutag4.png)](/index.php/File:Mutag4.png)

  5. Selecting the rotamer you think better fits your structure.   
  
Several side chain orientations (rotamers) are possible. You can use the back-and-forth movie controls (lower right corner) to display (in white) each of the rotamers available for this residue in PyMOL, whose current and total numbers are shown in the (green) Frame info. The rotamers are ordered according to their frequencies of occurrence in proteins, shown as a red percentage at the mutation object, which exists while mutagenesis is being performed.   
  

  6. Select Apply 
  7. Select Done 



## Explanation of colour codes

From [[1]](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg06196.html): 

The visible disks & colors in the Mutagenesis Wizard indicate pairwise overlap of atomic van der Waals radii. The intent is to provide a qualitative feedback regarding contacts and bumps. 

Short green lines or small green disks are shown when atoms are almost in contact or slightly overlapping. Large red disks indicate signficant van der Waals overlap. Everything else lies between those extremes. 

Retrieved from "[https://pymolwiki.org/index.php?title=Mutagenesis&oldid=12784](https://pymolwiki.org/index.php?title=Mutagenesis&oldid=12784)"


---

