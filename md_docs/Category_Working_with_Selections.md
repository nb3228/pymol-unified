# Category: Working with Selections

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

## Selection-expressions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Selection-expressions 

Selection-expressions stand for lists of atoms in arguments that are subject to PyMOL commands. You can name the selections to facilitate their re-use, or you can specify them anonymously (without names). Object and selection names may include the upper or lower case characters A/a to Z/z, numerals 0 to 9, and the underscore character (_). Characters to avoid include: 

! @ # $ % ^ &* ( ) ' " [ ] { } \ | ~ ` <> . ? / 

Selection-expressions describe the class of atoms you are referencing. Most of them require identifiers to complete the specification. For example, the selector resi references biopolymer residues by sequence number, and the identifier gives the number. The selector name references atoms according to the names described in the PDB, and the identifier gives the name (ca for alpha carbons, cb for beta carbons, etc). A handful of selection-expressions don't require identifiers, but most do. 

PyMOL uses several logical operators to increase the generality or specificity of selection-expressions. Logical combinations of selectors can get complex, so PyMOL accepts short forms and macros that express them with a minimum of keystrokes. This section describes named-selections, and then gives the syntax for making selections in a progression from simple one-word selectors to complex combinations of selectors, using macros and short forms. 

Retrieved from "[https://pymolwiki.org/index.php?title=Selection-expressions&oldid=8445](https://pymolwiki.org/index.php?title=Selection-expressions&oldid=8445)"


---

