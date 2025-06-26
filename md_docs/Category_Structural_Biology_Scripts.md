# Category: Structural Biology Scripts

## AAindex

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/aaindex.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/aaindex.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.aaindex](http://pymol.org/psicoredirect.php?psico.aaindex)  
  
## Contents

  * 1 Introduction
  * 2 Python Example
  * 3 PyMOL Example
  * 4 See Also



## Introduction

[![](/images/4/44/AAindexExample.png)](/index.php/File:AAindexExample.png)

[](/index.php/File:AAindexExample.png "Enlarge")

Hydrophobicity coloring with KYTJ820101

AAindex is a database of numerical indices representing various physicochemical and biochemical properties of amino acids and pairs of amino acids. See <http://www.genome.jp/aaindex/>

This script is a python parser for the AAindex flat files which will be downloaded from <ftp://ftp.genome.jp/pub/db/community/aaindex/>

The script provides two PyMOL commands (but can also be used without PyMOL). 

  * aaindex2b: Loads numerical indices from aaindex1 as b-factors into your structure
  * pmf: Potential of Mean Force (aaindex3)



## Python Example

Consider the script is called aaindex.py, it is placed somewhere in your PYTHONPATH and the aaindex flatfiles are found in the current directory. 
    
    
    import aaindex
    aaindex.init(path='.')
    
    aaindex.grep('volume')
    x = aaindex.get('KRIW790103')
    print(x)
    print(x.get('A'))
    
    aaindex.init(index='2')
    aaindex.grep('blosum')
    x = aaindex.get('HENS920102')
    print(x.get('A', 'K'))
    

## PyMOL Example

Solvent-accessible surface coloring by [Hydropathy index (Kyte-Doolittle, 1982)](https://www.genome.jp/dbget-bin/www_bget?aaindex:KYTJ820101)
    
    
    run aaindex.py
    aaindex2b KYTJ820101
    spectrum b, white yellow forest
    set surface_solvent
    show surface
    

## See Also

  * Protscale from [rTools](http://www.rubor.de/pymol_extensions_de.html) does a similar job in coloring by amino acid properties



Retrieved from "[https://pymolwiki.org/index.php?title=AAindex&oldid=13884](https://pymolwiki.org/index.php?title=AAindex&oldid=13884)"


---

## Angle between domains

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [orientation.py](https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/orientation.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.orientation](http://pymol.org/psicoredirect.php?psico.orientation)  
  
**angle_between_domains** calculates the rotation and displacement that would happen when [aligning](/index.php/Align "Align") two atom selections. The typical use case is for two conformations of a multi-domain structure, which is first aligned on domain 1 (e.g. chain A), and then **angle_between_domains** is calculated for domain 2 (e.g. chain B). 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See also



## Usage
    
    
    angle_between_domains selection1, selection2 [, method
        [, state1 [, state2 [, visualize ]]]]
    

## Arguments

  * **selection1** = str: atom selection in conformation 1
  * **selection2** = str: atom selection in conformation 2
  * **method** = str: alignment command like "align", "super" or "cealign" {default: align}



## Example
    
    
    # import the "angle_between_domains" command
    run https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/orientation.py
    
    # two conformations of a two-chain structure
    fetch 1ake, s1, async=0
    fetch 4ake, s2, async=0
    
    # align on chain A
    align s1 & chain A, s2 & chain A
    
    # measure rotation and displacement of chain B
    angle_between_domains s1 & chain B, s2 & chain B
    

For 1ake and 4ake, this reports: 
    
    
    Angle: 15.55 deg, Displacement: 56.08 angstrom
    

## See also

  * [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices")
  * [align](/index.php/Align "Align"), [super](/index.php/Super "Super"), [cealign](/index.php/Cealign "Cealign")
  * [pymol-users 01 Nov 2012](https://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10995.html)
  * [pymol-users 13 Jan 2015](https://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg13015.html)
  * [elbow_angle](/index.php/Elbow_angle "Elbow angle")
  * [RotationAxis](/index.php/RotationAxis "RotationAxis")



Retrieved from "[https://pymolwiki.org/index.php?title=Angle_between_domains&oldid=12487](https://pymolwiki.org/index.php?title=Angle_between_domains&oldid=12487)"


---

## AngleBetweenHelices

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/anglebetweenhelices.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/anglebetweenhelices.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.orientation](http://pymol.org/psicoredirect.php?psico.orientation)  
  
  * [![Orientation of a helix](/images/a/ad/Helix_orientation1.png)](/index.php/File:Helix_orientation1.png "Orientation of a helix")

Orientation of a helix 

  * [![Orientation of two helices, their angle that separates the two is 145.08 degrees as reported by the script](/images/5/5b/Helix_orientation2.png)](/index.php/File:Helix_orientation2.png "Orientation of two helices, their angle that separates the two is 145.08 degrees as reported by the script")

Orientation of two helices, their angle that separates the two is 145.08 degrees as reported by the script 




Calculate angle between alpha-helices or beta-sheets. There are four different methods implemented to fit a helix, two of them also work for sheets or loops. 

# Commands
    
    
    angle_between_helices selection1, selection2 [, method [, visualize ]]
    
    
    
    helix_orientation selection [, visualize [, sigma_cutoff ]]
    
    
    
    helix_orientation_hbond selection [, visualize [, cutoff ]]
    
    
    
    loop_orientation selection [, visualize ]
    
    
    
    cafit_orientation selection [, visualize ]
    

# Example
    
    
    import anglebetweenhelices
    
    fetch 2x19, async=0
    select hel1, /2x19//B/23-36/
    select hel2, /2x19//B/40-54/
     
    # just calculate/visualize orientation of single alpha-helix
    helix_orientation_hbond hel1
     
    # get angle between two helices
    angle_between_helices hel1, hel2
    angle_between_helices hel1, hel2, method=1
    angle_between_helices hel1, hel2, method=2
     
    # get angle between beta-sheets
    select sheet1, A/47-54/
    select sheet6, A/146-149/
    angle_between_helices sheet1, sheet6, method=loop_orientation
    angle_between_helices sheet1, sheet6, method=cafit_orientation
    

Output: 
    
    
    PyMOL>angle_between_helices hel1, hel2, method=cafit_orientation
    Using method: cafit_orientation
    Angle: 145.08 deg
    

## See Also

  * [angle_between_domains](/index.php/Angle_between_domains "Angle between domains")



Retrieved from "[https://pymolwiki.org/index.php?title=AngleBetweenHelices&oldid=13885](https://pymolwiki.org/index.php?title=AngleBetweenHelices&oldid=13885)"


---

## AutoMultiFit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# OVERVIEW

AutoMultiFit will fit all given chain(s) in all given state(s) to any other given chain(s) in any other given state(s) for some given selection within a homo-oligomer. So, if you have an MD trajectory or NMR ensemble of a tertramer, for example, over 20 states this program could calculate the RMS values for residues 10-30 in chain A's state 17, against all of chain C. Or all states in all chains against each other. 

See the notes in the code for usage. 

# The Code
    
    
    #
    # autoMultiFit.py -- Given a homo-oligomer in N-states, fit any combination of chains across states
    #
    # AUTHOR: Jason Vertrees
    # DATE  : 2009-11-02
    #
    import pymol
    from pymol import cmd
    
    def autoMultiFit(sel, fromChain="*", fromState="*", toChain="*", toState="*", verbose=None):
        """
        FUNCTION
          Given a homo-oligomer in N-states (say from MD), fit any
          combination of chains/selections across states. See Usage, Examples and Notes.
    
        USAGE
          autoMultiFit sel, [fromChain[, fromState[, toChain[, toState[, verbose ]]]]]
    
          sel
            The selection that must exist in all chains (what you're fitting)
    
          fromChain
            Fit from this chain, leave blank or use '*' for all chains;
    
          fromState
            Fit from this state, leave blank or use '*' for all states;
    
          toChain
            Fit to this chain only, use '*' or leave blank for all chains;
    
          toState
            Fit to this state only, use '*' or leave blank for all states
    
          verbose
            If turned on, this prints the fit for all pairs; this could be VERY slow.
    
        RETURNS
          A hash of hashes of hashes of hashes. :-)  Access the results like this:
          
          fitValues = autoMultiFit( mySelection )
          #fitValues[fromChain][fromState][toChain][toState]
          # to get the result of the fit from chain A, state 4 to chain C, state 17 type,
          fitValues['A']['4']['C']['17']
          
        EXAMPLES
          * Fit ALL chains to each other across ALL STATES for PDB 2kb1
          autoMultiFit 2kb1
    
          * Fit chain A's residue 22-34 against all states
          autoMutliFit 2kb1 and i. 22-34, fromChain=A
    
          * Fit chain A, B and C's residues 8-17, states 10 and 11, against chain D's state 8
          myVals = autoMultiFit("2kb1 and i. 8-17", fromChain="a b c", fromState="10 11", toChain="d", toState=8)
          myVals['B']['10']['D']['8']
    
          NOTES
            * The return value uses UPPERCASE and STRINGS for ALL hash keys, so if you used chain 'a' it will
              be 'A' in the map; if you used 'chainXX' it will be 'CHAINXX' in the map.
            * Continuing from the above note, all statese are accessible by their strings, so '17' gets state
              17 from the hash, not just the number 17.
            * This can be very slow: it's O( nFromChains * nFromStates * nToChains * nToStates )
            * fromChain, toChain, fromStates and toStates can all be single values, like "8", they can be a SPACE
              separated list of values, like "2 4 5", or they can be "*" which means ALL.
        """
    
        fromChain = processChain(fromChain,sel)
        fromState = processState(fromState, sel)
        toChain   = processChain(toChain, sel)
        toState   = processState(toState, sel)
    
        res = {}
    
        for frC in fromChain:
            res[frC] = {}
            for frS in fromState:
                res[frC][str(frS)] = {}
                cmd.create( "__tmpA", sel + " and c. " + frC, frS, 1 )
                for toC in toChain:
                    res[frC][str(frS)][toC] = {}
                    for toS in toState:
                        cmd.create("__tmpB", sel + " and c. " + toC, toS, 1 )
                        curVal = cmd.pair_fit( "__tmpA", "__tmpB", quiet=1 )
                        res[frC][str(frS)][toC][str(toS)] = curVal
                        if verbose!=None:
                            print "Pair Fitted: from (chain: %s, state: %s) to (chain: %s, states %s) with RMSD %s" % (frC, frS, toC, toS, curVal)
                        cmd.delete("__tmpB")
                cmd.delete("__tmpA")
        return res
    
    def processChain(ch, sel):
        """Turn a string-based space separated list into an array of chains"""
        if ch == "*":
            # just in case
            return map(str.upper, cmd.get_chains(sel))
        else:
            return map(str.upper, processTextList(ch))
    
    def processState(st, sel):
        """Tur a string-based space separated list into an array of states"""
        if st == "*":
            # assumes that if there is a state X, there is a state, X-1
            return range(cmd.count_states(sel))
        else:
            return map(int, processTextList(st))
    
    def processTextList(t):
        """ make a space-separated list into a real Python list, eg.
        input: a b c d
        output: [ 'a', 'b', 'c', 'd' ]
        """
        t = str(t).split()
        return t
    
        
    if __name__ == "__main__":
        assert processTextList( "a b c d" ) == ['a', 'b', 'c', 'd']
        assert processTextList( " 1 45 s s s s j p k") == [ '1', '45', 's', 's', 's', 's', 'j', 'p', 'k' ]
    
        assert processChain( "a b c dd", "blankSel" ) == ["a", "b", "c", "dd"]
        assert processState( "1 2 3 44 5", "blankSel" ) == [ 1, 2, 3, 44, 5 ]
    
    
    cmd.extend( "autoMultiFit", autoMultiFit )
    

Retrieved from "[https://pymolwiki.org/index.php?title=AutoMultiFit&oldid=8056](https://pymolwiki.org/index.php?title=AutoMultiFit&oldid=8056)"


---

## Average b

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 average_b.py
    * 1.1 usage
    * 1.2 author
    * 1.3 code



# average_b.py

  * calculate the average B-factor of a selection.



## usage

  * copy code to text file and save as average_b.py
  * type "run average_b.py" in PyMOL
  * then type "average_b (selection)"



## author

Gregor Hagelueken 

## code
    
    
    from pymol import cmd,stored
    def average_b(selection):
    	stored.tempfactor = 0
    	stored.atomnumber = 0
    	cmd.iterate(selection, "stored.tempfactor = stored.tempfactor + b")
    	cmd.iterate(selection, "stored.atomnumber = stored.atomnumber + 1")
    	print("Your selection: %s" % selection)
    	print("sum of B factors: %s" % stored.tempfactor)
    	print("number of atoms: %s" % stored.atomnumber)
    	averagetempfactor = stored.tempfactor / stored.atomnumber
    	print("average B of '%s': %s" % (selection, averagetempfactor))
    cmd.extend("average_b", average_b)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Average_b&oldid=13555](https://pymolwiki.org/index.php?title=Average_b&oldid=13555)"


---

## BbPlane

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/bbPlane.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/bbPlane.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Blaine Bell](/index.php?title=User:Bell&action=edit&redlink=1 "User:Bell \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO plane between the backbone atoms of two neighboring residues. This is to show the planarity of the atoms. The image style this is meant to represent can be found many places, like "Introduction to Protein Structure" by Branden and Tooze (2nd ed. pp. 8). 

  * [![Close up of planar atoms](/images/8/8d/BbPlane3.png)](/index.php/File:BbPlane3.png "Close up of planar atoms")

Close up of planar atoms 

  * [![A few more](/images/8/89/BbPlane1.png)](/index.php/File:BbPlane1.png "A few more")

A few more 

  * [![Global view](/images/b/b8/BbPlane2.png)](/index.php/File:BbPlane2.png "Global view")

Global view 




# Examples
    
    
    # download the source and save as bbPlane.py
    run bbPlane.py
    fetch 1cll
    # make planes for residues 4-9
    bbPlane i. 4-10
    

Retrieved from "[https://pymolwiki.org/index.php?title=BbPlane&oldid=13888](https://pymolwiki.org/index.php?title=BbPlane&oldid=13888)"


---

## BiologicalUnit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Note:** PyMOL 1.8 can load biological units from mmCIF files with the **[assembly](/index.php/Assembly "Assembly")** setting. 

This script can be used to re-create biological units for proteins. (This was created as a workaround of PyMOL's semi-functioning [Symexp](/index.php/Symexp "Symexp") command.) It's also a fun script to play with for learning about symmetry. 

  * [![single unit, before running this program](/images/a/a9/Before.png)](/index.php/File:Before.png "single unit, before running this program")

single unit, before running this program 

  * [![after the expansion of the 60 units](/images/c/c6/After.png)](/index.php/File:After.png "after the expansion of the 60 units")

after the expansion of the 60 units 

  * [![Example creating the biological unit for PDB 1RMV created from fiber diffraction.](/images/7/75/1rmv_fiber.png)](/index.php/File:1rmv_fiber.png "Example creating the biological unit for PDB 1RMV created from fiber diffraction.")

Example creating the biological unit for PDB 1RMV created from fiber diffraction. 




## Contents

  * 1 Usage
  * 2 The Code
  * 3 Notes
  * 4 See Also



# Usage
    
    
    load /path/to/some/pdbFile.pdb
    symMat = readSymmetry("/path/to/some/pdbFile.pdb","pdbFile")
    biologicalUnit("mates", "pdbFile", symMat)
    

# The Code
    
    
    #
    # Jason Vertrees <Jason-dot-Vertrees-at-schrodinger_dot_com>, 2010.
    #
    import pymol
    from pymol import cmd
    
    def readSymmetry(inFile, verbose=None):
      """
      This function will read "inFile" and glean the
      symmetry operations, if any, from it.
      
      PARAMS
        inFile
          (string) path to PDB file
          
        verbose
          (boolean) if verbose is not None, print more
          
      RETURNS
        matrix
          Array of lists.  One 16-element list per symmetry operation.  Feed this matrix
          into manualSymExp in order to make the other symmetry mates in the biological unit
      """
      # a remark-350 lines has:
      # REMARK 350 BIOMTn TAGn X Y Z Tx
      REM, TAG, BIOMT, OPNO, X, Y, Z, TX = range(8)
      
      thePDB = open(inFile, 'rb').readlines()
      
      matrices = []
      curTrans = -1
      
      # The transformation is,
      # output = U*input + Tx
      
      for l in thePDB:
        tokens = l.split()
        if len(tokens)!=8:
          continue
        if tokens[REM]=="REMARK" and tokens[TAG]=="350" and tokens[BIOMT].startswith("BIOMT"):
          if tokens[OPNO]!=curTrans:
            # new transformation matrix
            matrices.append([])
          
          matrices[-1].append( map( lambda s: float(s), tokens[X:]))
          curTrans = tokens[OPNO]
    
      if verbose!=None:
        print "Found %s symmetry operators in %s." % (len(matrices), inFile)
      return matrices
    
    
    def biologicalUnit(prefix, objSel, matrices ):
      """
      Manually expands the object in "objSel" by the symmetry operations provided in "matrices" and
      prefixes the new objects with "prefix".
      
      PARAMS
        prefix
          (string) prefix name for new objects
        
        objSel
          (string) name of object to expand
          
        matrices
          (list of 16-element lists) array of matrices from readSymmetry
          
        RETUNRS
          None
      
        SIDE EFFECTS
          Creates N new obects each rotated and translated according to the symmetry operators, where N
          equals len(matrices).
      """
      for m in matrices:
        n = cmd.get_unused_name(prefix)
        cmd.create(n, objSel)
        s1 = "%s + (x*%s + y*%s + z*%s)" % (m[0][3], m[0][0], m[0][1], m[0][2])
        s2 = "%s + (x*%s + y*%s + z*%s)" % (m[1][3], m[1][0], m[1][1], m[1][2])
        s3 = "%s + (x*%s + y*%s + z*%s)" % (m[2][3], m[2][0], m[2][1], m[2][2])
        cmd.alter_state(1, n, "(x,y,z) = (%s, %s, %s)" % (s1, s2, s3) )
    

# Notes

This is slow compared to [Symexp](/index.php/Symexp "Symexp"); use the above for learning, playing and when [Symexp](/index.php/Symexp "Symexp") doesn't work as advertised. 

# See Also

  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat") (alternative implementation)
  * [Symexp](/index.php/Symexp "Symexp")
  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [PDB Tutorial Biol. Units](http://www.rcsb.org/pdb/static.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html)
  * [Wikipedia article](http://en.wikipedia.org/wiki/Fiber_diffraction)
  * [assembly](/index.php/Assembly "Assembly")



Retrieved from "[https://pymolwiki.org/index.php?title=BiologicalUnit&oldid=12433](https://pymolwiki.org/index.php?title=BiologicalUnit&oldid=12433)"


---

## BiologicalUnit/Quat

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Note:** PyMOL 1.8 can load biological units from mmCIF files with the **[assembly](/index.php/Assembly "Assembly")** setting. 

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [(biomolecule) psico.xtal (biomolecule)](http://pymol.org/psicoredirect.php?psico.xtal)  
  
This script does more or less the same as [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit"). It creates the **quat** ernary structure (BIOMOLECULE 1 assembly) from the REMARK 350 header. 

This script is convenient to use because it searches automatically for the PDB file in the current directory, in [fetch_path](/index.php/Fetch_Path "Fetch Path") and (if available) in your local PDB mirror. 

_Also available from[psico](/index.php/Psico "Psico"), but the command is called **biomolecule** instead of **quat**. _

# Example Usage
    
    
    fetch 3bw1, type=pdb
    quat 3bw1
    as cartoon
    

# The Code
    
    
    '''
    (c) 2010-2011 Thomas Holder, MPI for Developmental Biology
    
    Module for reading REMARK records from PDB files and in particular
    generate quaterny structure from REMARK 350.
    '''
    
    import sys, os
    from pymol import cmd, stored
    
    local_mirror_divided = '/mnt/bio/db/pdb.divided'
    
    def pdbremarks(filename):
        '''
        Read REMARK lines from PDB file. Return dictionary with remarkNum as key
        and list of lines as value.
        '''
        remarks = dict()
        if not isinstance(filename, basestring):
            f = filename
        elif filename[-3:] == '.gz':
            import gzip
            f = gzip.open(filename)
        else:
            f = open(filename)
        for line in f:
            recname = line[0:6]
            if recname == 'REMARK':
                num = int(line[7:10])
                lstring = line[11:]
                remarks.setdefault(num, []).append(lstring)
        return remarks
    
    def quat350(rem350):
        '''
        Get transformation matrices for biomolecule 1 from REMARK 350.
        '''
        biomt = dict()
        chains = tuple()
        seenbiomolecule = False
        for line in rem350:
            if line.startswith('BIOMOLECULE:'):
                if seenbiomolecule:
                    break
                seenbiomolecule = True
            elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
                chains = tuple(chain.strip() for chain in line[30:].split(','))
            elif line.startswith('                   AND CHAINS:'):
                chains += tuple(chain.strip() for chain in line[30:].split(','))
            elif line.startswith('  BIOMT'):
                row = int(line[7])
                num = int(line[8:12])
                vec = line[12:].split()
                vec = map(float, vec)
                biomt.setdefault(chains, dict()).setdefault(num, []).extend(vec)
        return biomt
    
    def quat(name=None, filename=None, prefix=None, quiet=0):
        '''
    DESCRIPTION
    
        Read REMARK 350 from `filename` and create biological unit
        (quaternary structure)
    
    USAGE
    
        quat [name [, filename [, prefix]]]
    
    ARGUMENTS
    
        name = string: name of object and basename of PDB file, if
        filename is not given {default: first loaded object}
    
        filename = string: file path {default: <name>.pdb}
    
        prefix = string: prefix for new objects {default: <name>}
    
    EXAMPLE
    
        fetch 1rmv, type=pdb
        quat 1rmv
        '''
        quiet = int(quiet)
        if name is None:
            name = cmd.get_object_list()[0]
        if prefix is None:
            prefix = name
        if filename is None:
            candidates = [
                '%s.pdb' % (name),
                '%s/%s.pdb' % (cmd.get('fetch_path'), name),
                '%s/%s/pdb%s.ent.gz' % (local_mirror_divided, name[1:3], name),
            ]
            for filename in candidates:
                if os.path.exists(filename):
                    break
            else:
                print 'please provide filename'
                return
            if not quiet:
                print 'loading from %s' % (filename)
        remarks = pdbremarks(filename)
        if 350 not in remarks:
            print 'There is no REMARK 350 in', filename
            return
        quat = quat350(remarks[350])
        for chains in quat:
            matrices = quat[chains]
            for num in matrices:
                mat = matrices[num][0:12]
                mat.extend([0,0,0,1])
                copy = '%s_%d' % (prefix, num)
                if not quiet:
                    print 'creating %s' % (copy)
                cmd.create(copy, '/%s//%s' % (name, '+'.join(chains)))
                cmd.alter(copy, 'segi="%d"' % (num))
                cmd.transform_object(copy, mat)
        cmd.disable(name)
        cmd.group('%s_quat' % (prefix), '%s_*' % (prefix))
    
    cmd.extend('quat', quat)
    
    # vi:expandtab:smarttab
    

# See Also

  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [Supercell](/index.php/Supercell "Supercell")
  * [assembly](/index.php/Assembly "Assembly")



Retrieved from "[https://pymolwiki.org/index.php?title=BiologicalUnit/Quat&oldid=12743](https://pymolwiki.org/index.php?title=BiologicalUnit/Quat&oldid=12743)"


---

## Bondpack

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | <https://github.com/rasbt/BondPack>  
Author(s)  | Sebastian Raschka   
License  | GNU GENERAL PUBLIC LICENSE   
  
## Contents

  * 1 BondPack
  * 2 Introduction
  * 3 HydroBond
  * 4 BondVis
  * 5 BondLab
  * 6 Github



## BondPack

A collection of PyMOL plugins to visualize atomic bonds.  


  


## Introduction

PyMOL is without any doubt a great tool to visualize protein and ligand molecules.  
However, drawing interactions between atoms can be often quite cumbersome when done manually.  
For the sake of convenience, I developed three plugins for PyMOL that will make our life as protein biologists a little bit easier.  
All three PyMOL plugins can be installed and used separately; they don't depend on each other, but rather complement each other.  
At the end of this article, you will find brief instructions on how to install plugins in PyMOL - a very quick and simple process.  


## HydroBond

HydroBond visualizes all potential polar contacts between protein and ligand molecules within a user-specified distance.  
The underlying function is based on the different atom types, such as hydrogen bond acceptor and donor atoms,  
and thus it is required to have hydrogen atoms present in the structure.   
If your structure file doesn't contain hydrogen atoms already, you can add them directly in PyMOL as shown in the image below.  
  
[![Add hydrogens.png](/images/e/e0/Add_hydrogens.png)](/index.php/File:Add_hydrogens.png)   
  
HydroBond is related to PyMOLs "[A]->find->polar contacts" command, however,  
it doesn't consider geometry criteria and angle thresholds,  
but is rather based on atom types.  
When you select HydroBond from the "Plugin" menu, you will be prompted to enter the name of the protein object,  
the ligand object, and a distance cutoff as shown in the figure below.  
If HydroBond was able to detect hydrogen bond and acceptor atoms within the  
specified distance, potential interactions will be visualized as yellow dashed lines.  
  
[![Hydrobond action.png](/images/8/8b/Hydrobond_action.png)](/index.php/File:Hydrobond_action.png)   
  


## BondVis

The BondVis plugin lets you visualize interactions between any pair of atoms you specified.  
Often I find it helpful for my understanding (and for verification) to visualize the bonds between certain atoms  
that were assigned in docking or any other prediction software.  
Most software will provide you with information about the atoms that were "connected" to perform the analysis.  
If you run BondVis from the "Plugin" menu, it will ask you to select a "bond info file."  
  
[![Bondinfo.png](/images/3/33/Bondinfo.png)](/index.php/File:Bondinfo.png)   
  
This is just a simple text file that contains the atom numbers of connected atoms in pairs.  
Those can be separated by spaces, tabs, or commas. An example file with bond information could look like this:  


`
    
    
    1174		1357
    1175		1358
    1176		1359
    

`

  
When you selected a "bond info" file, BondVis will connect all atom pairs by yellow dashed lines  
and print out the connected atoms in the command field.  
  
[![Bondvis.png](/images/d/db/Bondvis.png)](/index.php/File:Bondvis.png)   
  


## BondLab

If you are not happy with the looks of the lines that were drawn to visualize connections between atoms,  
you will like to use BondLab.  
This plugin offers a simple interface that allows you to change the diameter,  
gap width, and color of the lines.   
  
[![Bondlab.png](/images/b/b3/Bondlab.png)](/index.php/File:Bondlab.png)  
The following video demonstrates the different features of BondLab:  
<http://youtu.be/14UZctxtK3w>   


## Github

If you are interested, you can follow the BondPack Github repository  
<https://github.com/rasbt/BondPack> for updates 

Retrieved from "[https://pymolwiki.org/index.php?title=Bondpack&oldid=11553](https://pymolwiki.org/index.php?title=Bondpack&oldid=11553)"


---

## Cart to frac

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will convert the real space orthonormal cartesian coordinates of a selection to the fractional coordinates given the loaded cell symmetry. 

  * Thanks to [Wikipedia](http://en.wikipedia.org/wiki/Fractional_coordinates).



# Example
    
    
    # load the script
    
    run cart_to_frac.py
    
    # fetch a protein and test
    
    fetch 1rx1, async=0
    
    # get the coordinates for the organic small
    # molecule as fractional coordinates
    
    print cart_to_frac("org")
    

# The Code
    
    
    from pymol import cmd
    
    def cart_to_frac(objSel,quiet=0,_self=cmd):
        """
        Returns an array of fractional atomic coordinates
        for a given object or selection.
    
        PARAMS
          objSel -- any object or selection
    
          quiet -- suppress output (default, no)
    
          _self -- core CMD object; or none
    
        RETURNS
          Python list of fractional coordinates
    
        NOTES/EXAMPLES
          cart_to_frac org
    
          x = cart_to_frac("solvent", quiet=1)
        """
        import numpy
        from numpy import cos, sin, sqrt
    
        a2r = numpy.pi / 180.0
    
        # get the model and coordinates
    
        m = _self.get_model(objSel)
        cart_coord = numpy.matrix(m.get_coord_list())
    
        # get the symmetry information
        try:
            a,b,c,alpha,beta,gamma,gp = _self.get_symmetry(objSel)
        except:
            print "Error-Failed to get symmetry. Please ensure you have a"
            print "valid object with proper crystal symmetry loaded."
            return None
    
        # convert to radians
    
        alpha = a2r * alpha
        beta  = a2r * beta
        gamma = a2r * gamma
        
        # (scaled) volume of the cell
        
        v = sqrt(1 -cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma) + 2*cos(alpha)*cos(beta)*cos(gamma))
    
        tmat = numpy.matrix( [
          [ 1.0 / a, -cos(gamma)/(a*sin(gamma)), (cos(alpha)*cos(gamma)-cos(beta)) / (a*v*sin(gamma))  ],
          [ 0.0,     1.0 / (b*sin(gamma)),         (cos(beta) *cos(gamma)-cos(alpha))/ (b*v*sin(gamma))  ],
          [ 0.0,     0.0,                        sin(gamma) / (c*v)                                    ] ]
          )
    
        r = cart_coord * tmat.T
    
        if not quiet:
            for (x,y,z) in r.tolist():
                print '%8.5f %8.5f %8.5f' % (x,y,z)
    
        # return the Nx3 results
        return r
        
    cmd.extend("cart_to_frac", cart_to_frac)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cart_to_frac&oldid=12442](https://pymolwiki.org/index.php?title=Cart_to_frac&oldid=12442)"


---

## Ccp4 contact

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ccp4_contact.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ccp4_contact.py)  
Author(s)  | [Gerhard Reitmayr and Dalia Daujotyte](/index.php?title=User:Dalyte&action=edit&redlink=1 "User:Dalyte \(page does not exist\)")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Example 1
  * 4 Getting a CONTACT file
    * 4.1 Install CCP4 - for Linux
    * 4.2 Use of CONTACT - for Linux



## Overview

[![](/images/6/64/HhaExample.png)](/index.php/File:HhaExample.png)

[](/index.php/File:HhaExample.png "Enlarge")

Interface residues (at cutoff <4A) in the 2c7r.pdb were found using NCONT, but similar results can be obtained using this script and CONTACT output. Usage of ccp4_contact script in PyMOL allows easy selection of residues and atoms listed in the output file. Interacting protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

The script selects residues and atoms from the list of the contacts found by CONTACT from CCP4 Program Suite (CONTACT analyses contacts between subsets of atoms in a PDB file). First, we run CONTACT on our pdb file to find interface residues. Then by using the ccp4_contact script in PyMOL we separately select residues and atoms listed in the output file. This generates two selections (atoms and residues) for each interacting chain, allowing quick manipulation of (sometimes) extensive lists in CONTACT log file. 

## Usage

ccp4_contact( contactsfile, selName1 = "source", selName2 = "target" ) 

  


## Example 1

First use CONTACT to find interface residues/atoms in the pdb file. Once you have the log file proceed to PyMOL. Make sure you import the ccp4_contact script first. 
    
    
    fetch 2c7r
    ccp4_contact 2c7r.contact, selName1=prot, selName2=dna
    

[![](/images/8/81/HhaI20example.png)](/index.php/File:HhaI20example.png)

[](/index.php/File:HhaI20example.png "Enlarge")

Quick and easy selection of interacting residues and atoms listed in the CONTACT log file. Protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

Download: [examples/ccp4_contact_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_contact_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_contact_1.pml>" highlight="python" />

## Getting a CONTACT file

### Install CCP4 - for Linux

Goto: <http://www.ccp4.ac.uk/download.php>   
Click: automated Downloads Pages   
Select: Linux, generic linux (x86)   
Select: Customized installation   
Select: Only CCP4 Program Suite, Executables -> Continue   
No additional packages -> Continue   
Download   


Extract for example to: **/home/YOU/Software/CCP** 4   
Then run:   

    
    
    $ /home/YOU/Software/CCP4/install.sh
    

write yes, read agreement, push y to agree license   
For sourcing scripts, say yes.   
See the changes to your environmental virables:   

    
    
    $ less ~/.bashrc
    

### Use of CONTACT - for Linux

See here for the CONTACT program and options: <http://www.ccp4.ac.uk/html/contact.html>   
Locate the pdb, and now run in terminal:   

    
    
    $ contact XYZIN 2c7r.pdb >> 2c7r.contact << eof   (#press enter)
    > MODE ISUB  (#press enter)
    > ATYPE NON-CARBON  (#press enter)
    > eof        (#press enter, and now the program runs, and shell saves to 2c7r.contact)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ccp4_contact&oldid=13889](https://pymolwiki.org/index.php?title=Ccp4_contact&oldid=13889)"


---

## Ccp4 ncont

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ccp4_ncont.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ccp4_ncont.py)  
Author(s)  | [Gerhard Reitmayr and Dalia Daujotyte](/index.php?title=User:Dalyte&action=edit&redlink=1 "User:Dalyte \(page does not exist\)")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 Getting a NCONT file
    * 4.1 Install CCP4 - for Linux
    * 4.2 Use of NCONT - for Linux



## Overview

[![](/images/6/64/HhaExample.png)](/index.php/File:HhaExample.png)

[](/index.php/File:HhaExample.png "Enlarge")

Interface residues (at cutoff <4A) in the 2c7r.pdb were found using NCONT. Usage of ccp4_ncont script in PyMOL allows easy selection of residues and atoms listed in ncont.log file. Interacting protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

The script selects residues and atoms from the list of the contacts found by NCONT from CCP4 Program Suite (NCONT analyses contacts between subsets of atoms in a PDB file). First, we run NCONT on our pdb file to find interface residues. Then by using the ccp4_ncont script in PyMOL we separately select residues and atoms listed in a ncont.log file. This generates two selections (atoms and residues) for each interacting chain, allowing quick manipulation of (sometimes) extensive lists in NCONT log file. 

This script works best for intermolecular contacts (when NCONT target and source selections don't overlap). If crystal contacts (NCONT parameter cell = 1 or 2) are included then additional coding is required to distinguish inter from intramolecular contacts. 

## Usage

ccp4_ncont( contactsfile, selName1 = "source", selName2 = "target" ) 

  


## Examples

First use NCONT to find interface residues/atoms in the pdb file. Once you have ncont.log file proceed to PyMOL. Make sure you import the ccp4_ncont script first. 
    
    
    fetch 2c7r
    ccp4_ncont 2c7r.ncont, selName1=prot, selName2=dna
    

[![](/images/8/81/HhaI20example.png)](/index.php/File:HhaI20example.png)

[](/index.php/File:HhaI20example.png "Enlarge")

Quick and easy selection of interacting residues and atoms listed in the NCONT log file. Protein and DNA residues are colored in red and slate, respectively. Atoms in contact are shown in dots.

Download: [examples/ccp4_ncont_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_ncont_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_ncont_1.pml>" highlight="python" />

## Getting a NCONT file

### Install CCP4 - for Linux

Goto: <http://www.ccp4.ac.uk/download.php>   
Click: automated Downloads Pages   
Select: Linux, generic linux (x86)   
Select: Customized installation   
Select: Only CCP4 Program Suite, Executables -> Continue   
No additional packages -> Continue   
Download   


Extract for example to: **/home/YOU/Software/CCP** 4   
Then run:   

    
    
    $ /home/YOU/Software/CCP4/install.sh
    

write yes, read agreement, push y to agree license   
For sourcing scripts, say yes.   
See the changes to your environmental virables:   

    
    
    $ less ~/.bashrc
    

### Use of NCONT - for Linux

See here for the NCONT program and options:   
<http://www.ccp4.ac.uk/html/ncont.html>   
<http://www.ccp4.ac.uk/html/pdbcur.html#atom_selection>   
Locate the pdb, and now run in terminal:   

    
    
    $ ncont XYZIN 2c7r.pdb >> 2c7r.ncont << eof   (#press enter)
    > source A    (#press enter)
    > target C,D  (#press enter)
    > eof         (#press enter, and now the program runs, and shell saves to 2c7r.ncont)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ccp4_ncont&oldid=13890](https://pymolwiki.org/index.php?title=Ccp4_ncont&oldid=13890)"


---

## Ccp4 pisa

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/ccp4_pisa.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/ccp4_pisa.py)  
Author(s)  | [Gerhard Reitmayr and Dalia Daujotyte](/index.php?title=User:Dalyte&action=edit&redlink=1 "User:Dalyte \(page does not exist\)")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Overview

The script selects atoms from the list of the contacts found by PISA. First, we run PISA on our pdb file to find the interfaces. Then by using the ccp4_pisa script in PyMOL we separately select atoms for all interface types and individual interfaces. This generates many selections, two for each interface, allowing quick manipulation of (sometimes) extensive lists in PISA log file. 

## Usage

ccp4_pisa( pisafile ) 

  


## Example 1

The script parses the XML output files from the PISA service or command line tool. A short description of how to download the XML output files is available here <http://www.ebi.ac.uk/msd-srv/prot_int/pi_download.html>. 

(For example, the following URL downloads the interfaces in 2c7r.pdb <http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?2c7r>) 

Make sure you import the ccp4_pisa script first. 

Download: [examples/ccp4_pisa_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_pisa_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/ccp4_pisa_1.pml>" highlight="python" />

Retrieved from "[https://pymolwiki.org/index.php?title=Ccp4_pisa&oldid=13891](https://pymolwiki.org/index.php?title=Ccp4_pisa&oldid=13891)"


---

## Centroid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/centroid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/centroid.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
    * 3.1 See Also



## Overview

Centroid is a small script that returns the value of the geometric center (or centroid) of your selection. It also can translate the object of that selection to the origin. 

## Syntax
    
    
    centroid (selection=PyMOLSelection), [center=boolean]
    

## Examples
    
    
    # get the centroid of the polymer
    import centroid
    fetch 4ins, async=0
    centroid polymer
    
    # move some 'ligand' to the origin
    centroid ligand, center=1
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [Center_Of_Mass](/index.php/Center_Of_Mass "Center Of Mass")



Retrieved from "[https://pymolwiki.org/index.php?title=Centroid&oldid=13893](https://pymolwiki.org/index.php?title=Centroid&oldid=13893)"


---

## Color by conservation

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/color_by_conservation.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/color_by_conservation.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | Free   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview

This script reads an alignment object and colors the protein objects in the alignment by the sequence conservation found in the alignment. 

[![Cbcon.png](/images/9/94/Cbcon.png)](/index.php/File:Cbcon.png)

# Example Usage
    
    
    reinitialize
    import color_by_conservation
     
    # get some kinases
    fetch 1opk 3dtc 3p86 2eva 3efw, async=0
     
    # turn on the sequence viewer
    set seq_view
    
    # align them into the "algn" object
    for x in cmd.get_names(): cmd.align(x, "3efw and c. A", object="algn")
     
    # color
    color_by_conservation aln=algn, as_putty=1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_by_conservation&oldid=13895](https://pymolwiki.org/index.php?title=Color_by_conservation&oldid=13895)"


---

## Color By Mutations

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script creates an alignment of two proteins and superimposes them. Aligned residues that are different in the two (i.e. mutations) are highlighted and colored according to their difference in the BLOSUM90 matrix. It is meant to be used for similar proteins, e.g. close homologs or point mutants, to visualize their differences. 

Example 
    
    
    color_by_mutation protein1, protein2
    

Example: rat trypsin and cow trypsin colored by color_by_mutation. 

[![Color by mutation.png](/images/0/06/Color_by_mutation.png)](/index.php/File:Color_by_mutation.png)

  

    
    
    '''
    created by Christoph Malisi.
    
    Creates an alignment of two proteins and superimposes them. 
    Aligned residues that are different in the two (i.e. mutations) are highlighted and 
    colored according to their difference in the BLOSUM90 matrix. 
    Is meant to be used for similar proteins, e.g. close homologs or point mutants, 
    to visualize their differences.
    
    '''
    
    from pymol import cmd
    
    aa_3l = {'ALA':0, 'ARG':1, 'ASN':2, 'ASP':3, 'CYS':4, 'GLN':5, 'GLU':6, 'GLY':7, 'HIS':8, 'ILE':9, 'LEU':10, 'LYS':11,
            'MET':12, 'PHE':13, 'PRO':14, 'SER':15, 'THR':16, 'TRP':17, 'TYR':18, 'VAL':19, 'B':20, 'Z':21, 'X':22, '*':23}
    
    #            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    blosum90 = [[ 5, -2, -2, -3, -1, -1, -1,  0, -2, -2, -2, -1, -2, -3, -1,  1,  0, -4, -3, -1, -2, -1, -1, -6], 
                [-2,  6, -1, -3, -5,  1, -1, -3,  0, -4, -3,  2, -2, -4, -3, -1, -2, -4, -3, -3, -2,  0, -2, -6], 
                [-2, -1,  7,  1, -4,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -5, -3, -4,  4, -1, -2, -6], 
                [-3, -3,  1,  7, -5, -1,  1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5,  4,  0, -2, -6], 
                [-1, -5, -4, -5,  9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -5, -3, -6], 
                [-1,  1,  0, -1, -4,  7,  2, -3,  1, -4, -3,  1,  0, -4, -2, -1, -1, -3, -3, -3, -1,  4, -1, -6], 
                [-1, -1, -1,  1, -6,  2,  6, -3, -1, -4, -4,  0, -3, -5, -2, -1, -1, -5, -4, -3,  0,  4, -2, -6], 
                [ 0, -3, -1, -2, -4, -3, -3,  6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -3, -2, -6], 
                [-2,  0,  0, -2, -5,  1, -1, -3,  8, -4, -4, -1, -3, -2, -3, -2, -2, -3,  1, -4, -1,  0, -2, -6], 
                [-2, -4, -4, -5, -2, -4, -4, -5, -4,  5,  1, -4,  1, -1, -4, -3, -1, -4, -2,  3, -5, -4, -2, -6], 
                [-2, -3, -4, -5, -2, -3, -4, -5, -4,  1,  5, -3,  2,  0, -4, -3, -2, -3, -2,  0, -5, -4, -2, -6], 
                [-1,  2,  0, -1, -4,  1,  0, -2, -1, -4, -3,  6, -2, -4, -2, -1, -1, -5, -3, -3, -1,  1, -1, -6], 
                [-2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7, -1, -3, -2, -1, -2, -2,  0, -4, -2, -1, -6], 
                [-3, -4, -4, -5, -3, -4, -5, -5, -2, -1,  0, -4, -1,  7, -4, -3, -3,  0,  3, -2, -4, -4, -2, -6], 
                [-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8, -2, -2, -5, -4, -3, -3, -2, -2, -6], 
                [ 1, -1,  0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2,  5,  1, -4, -3, -2,  0, -1, -1, -6], 
                [ 0, -2,  0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2,  1,  6, -4, -2, -1, -1, -1, -1, -6], 
                [-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2,  0, -5, -4, -4, 11,  2, -3, -6, -4, -3, -6], 
                [-3, -3, -3, -4, -4, -3, -4, -5,  1, -2, -2, -3, -2,  3, -4, -3, -2,  2,  8, -3, -4, -3, -2, -6], 
                [-1, -3, -4, -5, -2, -3, -3, -5, -4,  3,  0, -3,  0, -2, -3, -2, -1, -3, -3,  5, -4, -3, -2, -6], 
                [-2, -2,  4,  4, -4, -1,  0, -2, -1, -5, -5, -1, -4, -4, -3,  0, -1, -6, -4, -4,  4,  0, -2, -6], 
                [-1,  0, -1,  0, -5,  4,  4, -3,  0, -4, -4,  1, -2, -4, -2, -1, -1, -4, -3, -3,  0,  4, -1, -6], 
                [-1, -2, -2, -2, -3, -1, -2, -2, -2, -2, -2, -1, -1, -2, -2, -1, -1, -3, -2, -2, -2, -1, -2, -6], 
                [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1]] 
    
    def getBlosum90ColorName(aa1, aa2):
        '''returns a rgb color name of a color that represents the similarity of the two residues according to
           the BLOSUM90 matrix. the color is on a spectrum from blue to red, where blue is very similar, and 
           red very disimilar.'''
        # return red for residues that are not part of the 20 amino acids
        if aa1 not in aa_3l or aa2 not in aa_3l:
            return 'red'
        
        # if the two are the same, return blue
        if aa1 == aa2:
            return 'blue'
        i1 = aa_3l[aa1]
        i2 = aa_3l[aa2]
        b = blosum90[i1][i2]
    
        # 3 is the highest score for non-identical substitutions, so substract 4 to get into range [-10, -1]
        b = abs(b - 4)
    
        # map to (0, 1]:
        b = 1. - (b / 10.0)
    
        # red = [1.0, 0.0, 0.0], blue = [0.0, 0.0, 1.0]
        colvec = [(0., 0., 1.), (1.,0.,0.)]
        bcolor = (1.-b, 0., b)
        col_name = '0x%02x%02x%02x' % tuple(int(b * 0xFF) for b in bcolor)
        return col_name
    
    def color_by_mutation(obj1, obj2, waters=0, labels=0):
        '''
        DESCRIPTION
        
            Creates an alignment of two proteins and superimposes them. 
            Aligned residues that are different in the two (i.e. mutations) are highlighted and 
            colored according to their difference in the BLOSUM90 matrix. 
            Is meant to be used for similar proteins, e.g. close homologs or point mutants, 
            to visualize their differences.      
        
        USAGE
        
            color_by_mutation selection1, selection2 [,waters [,labels ]]
            
        ARGUMENTS
        
            obj1: object or selection
            
            obj2: object or selection    
            
            waters: bool (0 or 1). If 1, waters are included in the view, colored
                    differently for the both input structures.
                    default = 0
    
            labels: bool (0 or 1). If 1, the possibly mutated sidechains are 
                    labeled by their chain, name and id
                    default = 0
            
        EXAMPLE
            
            color_by_mutation protein1, protein2
            
        SEE ALSO
    
            super
        '''
        from pymol import stored, CmdException
    
        if cmd.count_atoms(obj1) == 0:
            print('%s is empty'%obj1)
            return
        if cmd.count_atoms(obj2) == 0:
            print('%s is empty'%obj2)
            return
        waters = int(waters)
        labels = int(labels)
        
        # align the two proteins
        aln = '__aln'
        
        # first, an alignment with 0 cycles (no atoms are rejected, which maximized the number of aligned residues)
        # for some mutations in the same protein this works fine). This is essentially done to get a
        # sequence alignment
        cmd.super(obj1, obj2, object=aln, cycles=0)
        
        # superimpose the the object using the default parameters to get a slightly better superimposition,
        # i.e. get the best structural alignment    
        cmd.super(obj1, obj2)
    
        stored.resn1, stored.resn2 = [], []
        stored.resi1, stored.resi2 = [], []
        stored.chain1, stored.chain2 = [], []
        
        # store residue ids, residue names and chains of aligned residues
        cmd.iterate(obj1 + ' and name CA and ' + aln, 'stored.resn1.append(resn)')
        cmd.iterate(obj2 + ' and name CA and ' + aln, 'stored.resn2.append(resn)')
    
        cmd.iterate(obj1 + ' and name CA and ' + aln, 'stored.resi1.append(resi)')
        cmd.iterate(obj2 + ' and name CA and ' + aln, 'stored.resi2.append(resi)')
    
        cmd.iterate(obj1 + ' and name CA and ' + aln, 'stored.chain1.append(chain)')
        cmd.iterate(obj2 + ' and name CA and ' + aln, 'stored.chain2.append(chain)')
        
        
        mutant_selection = '' 
        non_mutant_selection = 'none or '
        colors = []
    
        # loop over the aligned residues
        for n1, n2, i1, i2, c1, c2 in zip(stored.resn1, stored.resn2,
                                          stored.resi1, stored.resi2,
                                          stored.chain1, stored.chain2):
            # take care of 'empty' chain names
            if c1 == '':
                c1 = '""'
            if c2 == '':
                c2 = '""'
            if n1 == n2:
                non_mutant_selection += '((%s and resi %s and chain %s) or (%s and resi %s and chain %s)) or '%(obj1, i1, c1, obj2, i2, c2 )            
            else:
                mutant_selection += '((%s and resi %s and chain %s) or (%s and resi %s and chain %s)) or '%(obj1, i1, c1, obj2, i2, c2 )
                # get the similarity (according to the blosum matrix) of the two residues and
                c = getBlosum90ColorName(n1, n2)
                colors.append((c, '%s and resi %s and chain %s and elem C'%(obj2, i2, c2)))
    
        if mutant_selection == '':
            raise CmdException('No mutations found')
    
        # create selections        
        cmd.select('mutations', mutant_selection[:-4])
        cmd.select('non_mutations', non_mutant_selection[:-4])
        cmd.select('not_aligned', '(%s or %s) and not mutations and not non_mutations'%(obj1, obj2))
        
        # create the view and coloring
        cmd.hide('everything', '%s or %s'%(obj1, obj2))
        cmd.show('cartoon', '%s or %s'%(obj1, obj2))
        cmd.show('lines', '(%s or %s) and ((non_mutations or not_aligned) and not name c+o+n)'%(obj1, obj2))
        cmd.show('sticks', '(%s or %s) and mutations and not name c+o+n'%(obj1, obj2))
        cmd.color('gray', 'elem C and not_aligned')
        cmd.color('white', 'elem C and non_mutations')
        cmd.color('blue', 'elem C and mutations and %s'%obj1)
        for (col, sel) in colors:
            cmd.color(col, sel)
    
        cmd.hide('everything', '(hydro) and (%s or %s)'%(obj1, obj2))        
        cmd.center('%s or %s'%(obj1, obj2))
        if labels:
            cmd.label('mutations and name CA','"(%s-%s-%s)"%(chain, resi, resn)')    
        if waters:
            cmd.set('sphere_scale', '0.1')
            cmd.show('spheres', 'resn HOH and (%s or %s)'%(obj1, obj2))
            cmd.color('red', 'resn HOH and %s'%obj1)
            cmd.color('salmon', 'resn HOH and %s'%obj2)
        print('''
                 Mutations are highlighted in blue and red.
                 All mutated sidechains of %s are colored blue, the corresponding ones from %s are
                 colored on a spectrum from blue to red according to how similar the two amino acids are
                 (as measured by the BLOSUM90 substitution matrix).
                 Aligned regions without mutations are colored white.
                 Regions not used for the alignment are gray.
                 NOTE: There could be mutations in the gray regions that were not detected.'''%(obj1, obj2))
        cmd.delete(aln)
        cmd.deselect()
        
    
        
    cmd.extend("color_by_mutation", color_by_mutation)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_By_Mutations&oldid=13306](https://pymolwiki.org/index.php?title=Color_By_Mutations&oldid=13306)"


---

## Colorblindfriendly

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/colorblindfriendly.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorblindfriendly.py)  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  


## Contents

  * 1 Introduction
  * 2 Colors
  * 3 Usage
    * 3.1 Import as a module
    * 3.2 Run the latest version of the script from Github
    * 3.3 Download the script and run locally
  * 4 Requirements



## Introduction

Certain colors are indistinguishable to people with the various forms of color blindness, and therefore are better not used in figures intended for public viewing. 

This script generates a palette of named colors for PyMOL that are unambiguous both to colorblind and non-colorblind individuals. 

The colors listed here are defined according to recommendations by Okabe and Ito previously available at [J*FLY (archived)](https://web.archive.org/web/20170608112249/http://jfly.iam.u-tokyo.ac.jp/color/#pallet). This website is a good reference to consult when making all kinds of figures, not just those made using PyMOL. 

## Colors

These are the 0-255 RGB values from the J*FLY page that are used in the script, with the defined color names and alternate names. 

| name | R | G | B | alternate names   
---|---|---|---|---|---  
| cb_black | 0 | 0 | 0 |   
| cb_orange | 230 | 159 | 0 |   
| cb_sky_blue | 86 | 180 | 233 | cb_skyblue, cb_light_blue, cb_lightblue   
| cb_bluish_green | 0 | 158 | 115 | cb_bluishgreen, cb_green   
| cb_yellow | 240 | 228 | 66 |   
| cb_blue | 0 | 114 | 178 |   
| cb_vermillion | 213 | 94 | 0 | cb_red, cb_red_orange, cb_redorange   
| cb_reddish_purple | 204 | 121 | 167 | cb_rose, cb_violet, cb_magenta   
  
## Usage

### Import as a module

[![](/images/c/ce/Colorblindfriendly_menu.png)](/index.php/File:Colorblindfriendly_menu.png)

[](/index.php/File:Colorblindfriendly_menu.png "Enlarge")

Screenshot of the `cb_colors` color menu in the OpenGL GUI in PyMOL 2.0.

After importing the module, call the `set_colors()` function to add the colors to PyMOL's color palette. Then, use these color names just like any other named color, using the `[color](/index.php/Color "Color")` command. 
    
    
    import colorblindfriendly as cbf
    
    # Add the new colors
    cbf.set_colors()
    color cb_red, myObject
    

The colors can also be made to replace the built-in colors (i.e. they are created both with and without the "`cb_`" prefix.). Do this by passing the `replace` keyword argument. 
    
    
    # Replace built-in colors with cbf ones
    cbf.set_colors(replace=True)
    color yellow, myOtherObject   # actually cb_yellow
    

One can also add an entry to the color menu in the right-side OpenGL GUI. So clicking on [C], there will now be a `cb_colors` menu item, which expands to give all the color blind-friendly colors, except black, which is available in the `grays` menu. 
    
    
    # Add a `cb_colors` menu to the OpenGL GUI (see screenshot)
    # This will also add the colors if `set_colors()` hasn't yet been run.
    cbf.add_menu()
    

### Run the latest version of the script from Github

In a PyMOL session, run at the command line: 
    
    
    run <https://github.com/Pymol-Scripts/Pymol-script-repo/raw/master/colorblindfriendly.py>
    

This will add all the colors as well as the OpenGL menu. 

  


### Download the script and run locally

Save the script from the link in the box at the top right to your computer. 

In a PyMOL session, run at the command line: 
    
    
    run /path/to/colorblindfriendly.py
    

This will add all the colors as well as the OpenGL menu. 

## Requirements

The `cb_colors` GUI menu (generated by the `add_menu()` function) requires PyMOL 1.6.0 and later. 

Retrieved from "[https://pymolwiki.org/index.php?title=Colorblindfriendly&oldid=13896](https://pymolwiki.org/index.php?title=Colorblindfriendly&oldid=13896)"


---

## Colorbydisplacement

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/colorbydisplacement.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorbydisplacement.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Bug in code
    * 2.1 Examples
  * 3 Example 1



## Introduction

This script allows you to color two structures by distance displacement between an Open and Closed form of a protein, as calculated by PyMol's internal distance command. The pairwise distance is calculated between all-atoms. The distance displacement values are stored as B-factors of these residues, which are colored by a _rainbow_ color spectrum, with blue specifying minimum and red indicating maximum. 

Do keep in mind, all original B-factors values are overwritten! 

There exist one version.   
ColorByDisplacement**All** is between All atoms in residues and is quite slow => 3-5 mins for a run. Ideal for sticks representation. 

**You have to specify which residues should be used in the alignment procedure, or it will take all residues as standard**

V.2 is implemented the 2011.01.06 - Due to a bug in coloring. 

## Bug in code

A bug in the boolean operator of the spectrum command has been found. This versions work for version 1.3 Educational product.   
For other versions of pymol, try to change (comment/uncomment) the **cmd.spectrum** line. The other spectrum line works for Open-Source PyMOL 1.2r3pre, Incentive product 

### Examples
    
    
    ColorByDisplacementAll O5NT, C5NT, super1=resi 26-355, super2=resi 26-355, doColor=t, doAlign=t
    
    ColorByDisplacementAll O5NT, C5NT, super1=resi 26-355, super2=resi 26-355, doColor=t, doAlign=t, AlignedWhite='no'
    

  * [![ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.](/images/a/a2/ColorByDisplacement-All-1.png)](/index.php/File:ColorByDisplacement-All-1.png "ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.")

ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement. 

  * [![ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.](/images/6/69/ColorByDisplacement-All-2.png)](/index.php/File:ColorByDisplacement-All-2.png "ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement.")

ColorByDisplacementAll used on 1HP1 and 1HPU aligned and colored by distance displacement. 




Dark blue is low displacement, higher displacements are in orange/yellow/red.   
Residues used for alignment is colored white. Can be turned off in top of algorithm. Residues not in both pdb files is colored black 

## Example 1

Download: [examples/colorbydisplacement_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/colorbydisplacement_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/colorbydisplacement_1.pml>" highlight="python" />

Retrieved from "[https://pymolwiki.org/index.php?title=Colorbydisplacement&oldid=13897](https://pymolwiki.org/index.php?title=Colorbydisplacement&oldid=13897)"


---

## ColorByRMSD

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/colorbyrmsd.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/colorbyrmsd.py)  
Author(s)  | [Shivender Shandilya](/index.php/User:Shiven "User:Shiven"), [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Arguments
  * 4 Examples



## Introduction

This script allows you to color two structures by Root Mean Square Deviation (RMSD). The distances between aligned C-alpha atom pairs are stored as B-factors of these residues, which are colored by a color spectrum, with blue specifying the minimum pairwise RMSD and red indicating the maximum. Unaligned residues are colored gray. 

## Usage
    
    
    colorbyrmsd mobile, target [, doAlign [, doPretty [, guide [, method ]]]]
    

## Arguments

  * **mobile** = string: atom selection for mobile atoms
  * **target** = string: atom selection for target atoms
  * **doAlign** = 0 or 1: Superpose selections before calculating distances {default: 1}
  * **doPretty** = 0 or 1: Show nice representation and colors {default: 1}
  * **guide** = 0 or 1: Only use C-alpha atoms {default: 1}
  * **method** = align or super: Method to match atoms {default: super}



## Examples
    
    
    # example #1
    colorbyrmsd 1cbs, 1hmt, doAlign=1, doPretty=1
    # example #2
    colorbyrmsd 1eaz, 1fao, doAlign=1, doPretty=1
    

  * [![1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.](/images/5/57/ColorByRMSD_1cbs_1hmt.png)](/index.php/File:ColorByRMSD_1cbs_1hmt.png "1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.")

1cbs and 1hmt aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white. 

  * [![1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.](/images/b/b6/ColorByRMSD_1eaz_1fao.png)](/index.php/File:ColorByRMSD_1eaz_1fao.png "1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white.")

1eaz and 1fao aligned and colored by RMSD. Dark blue is good alignment, higher deviations are in orange/yellow/red. Residues not used for alignment are colored white. 




Retrieved from "[https://pymolwiki.org/index.php?title=ColorByRMSD&oldid=13898](https://pymolwiki.org/index.php?title=ColorByRMSD&oldid=13898)"


---

## Consistent View/ Map Inspect

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**DESCRIPTION**

This is a toolkit for rapidly inspecting multiple maps and models. The right & left arrow keys navigate sequentially through the amino acid chain, but alternating between two structures (could be NCS-related chains or models solved in different space groups). Each view is rendered in a consistent orientation (the default is centered on alpha carbon, with nitrogen right, carbonyl left & beta carbon up). The view can be customized. It is necessary to edit the script to define the behavior for the problem at hand. 

**INSTALLATION**

Three components are needed: 

  * A Python module, support.py, giving matrix algebra & Kabsch-style structure alignment (source code given below). This module must be in the current working directory, and PyMol started from the command line, e.g., on MacOSX:


    
    
    /Applications/MacPyMol.app/Contents/MacOS/MacPyMOL
    

  * Atomic coordinates must be loaded first, and if appropriate, maps and isomesh levels. An example .pml file is shown here, consisting of two structural models solved in different space groups:


    
    
    load sand2b_001_2mFo-DFc_map.xplor, map2
    isomesh msh2,map2,2.0
    color blue, msh2
    
    load sand3_rotfix_001_2mFo-DFc_map.xplor, map3
    isomesh msh3,map3,2.0
    color marine, msh3
    
    load sand2b_001.pdb, mod02
    load sand3_rotfix_001.pdb, mod03
    set antialias, 2
    set line_width, 3
    

  * A final pml is then loaded (@consistent.pml) to define the viewing sequence for the problem at hand. An example is shown here. (The python / python end construct will not work in PyMOL versions below 1.0.)


    
    
    python
    
    from support import matrix,residue,std_residue,kabsch_rotation
    from support import view_matrix, new_set_view
    
    #  This section is the generic preset view--DO NOT EDIT #########################
    #  Centered on CA, CB up, N left, C=O right
    preset_view = view_matrix((
         0.295077115,    0.440619737,   -0.847830832,
        -0.408159286,    0.860427082,    0.305112839,
         0.863915384,    0.256014407,    0.433732271,
         0.,             0.,           -30.,
         0.,             0.,             0.,
        26.,            33.,             0.
    ))
    preset_view.difference_to_this_CA = matrix.col((0.,0.,0.,))
    
    preset_std_residue = std_residue()
    preset_std_residue.N = matrix.col((16.1469993591, 33.0660018921, -20.0760002136))
    preset_std_residue.CA = matrix.col((15.4309997559, 34.2599983215, -20.4640007019))
    preset_std_residue.C = matrix.col((15.2889995575, 34.1189994812, -21.9699993134))
    preset_std_residue.CB = matrix.col((16.2469997406, 35.516998291, -20.1380004883))
    #  END GENERIC VIEW #############################################################
    
    def view_preset():
      cmd.current_std_view = preset_view
      cmd.current_std_residue = preset_std_residue
    
    view_preset()
    
    #  EDIT THIS SECTION TO REFLECT THE PROBLEM AT HAND #############################
    #  In this example, there are two objects: structure mod02 & mod03
    #  We will be looking at chain B in each structure, each having 262 residues,
    #  for a total of 524 amino acids to view. 
    #  For each of 524 views, we select the appropriate map and isomesh
    
    cmd.no_of_structures = 2
    cmd.global_ptr = 1 * cmd.no_of_structures
    cmd.last_residue = preset_std_residue
    
    def seqwalk():
      chain_id = "B"
    
      if cmd.global_ptr<524:
      
        resid = cmd.global_ptr/2
        this_object = cmd.global_ptr%2
        if this_object==0:
          atoms = cmd.get_model("object mod03 and chain %s"%chain_id)
          obj_id = "mod03"
          map = "map3"
          mesh = "msh3"
        else:
          atoms = cmd.get_model("object mod02 and chain %s"%chain_id)
          obj_id = "mod02"
          map = "map2"
          mesh = "msh2"
        
        this = residue(chain_id,resid,atoms)
        cmd.last_residue = this #remember this residue in case reset command is issued
        print "Centering on /%s//%s/%s %d/CA; waiting for right mouse key"%(
          obj_id,chain_id,this.residue_type,resid,)
        cmd.center("/%s//%s/%d/CA"%(obj_id,chain_id,resid))
        cmd.isomesh(
         name=mesh,map=map,level= 2.0, selection="object %s and chain %s and resid %d"%(obj_id,chain_id,resid),buffer=8)
    #  No more edits below this line
        
        kr = kabsch_rotation(cmd.current_std_residue.sites(),this.sites())
        view_matrix = new_set_view(cmd.current_std_view,kr,this,verbose=False).view_matrix()
        cmd.set_view(view_matrix)
        cmd.global_ptr=cmd.global_ptr+1
    
    #  END OF PROBLEM-SPECIFIC SECTION TO EDIT  #####################################
    
    def seqwalk2():
      cmd.global_ptr-=2
      seqwalk()
        
    cmd.set_key("right",seqwalk)
    cmd.set_key("left",seqwalk2)
    
    def view_reset():
      cmd.current_std_view = view_matrix(cmd.get_view()).set_diff_to_CA(cmd.last_residue)
      cmd.current_std_residue = cmd.last_residue
    
    def view_goto(number):
      cmd.global_ptr = cmd.no_of_structures*number
    
    cmd.extend("view_reset",view_reset)  # After customizing the view with GUI mouse & clicks, make the
                                         # view persistent by typing "view reset"
    cmd.extend("view_preset",view_preset)# Forget the user-customized view, go back to the generic view
                                         # defined at the beginning of the program
    cmd.extend("view_goto",view_goto)    # Example: view_goto(36) --resets at residue 36
    
    python end
    

**USAGE**

Once the code is set up as described above, you must mouse-click in PyMOL's 3D display window to enable the right & left arrow keys. These keys will navigate through the amino acid sequence, displaying each residue in a predefined orientation. Maps, if defined, will be redrawn for each view. 

Three commands can be issued at the prompt: 

  * view_goto(N) -- will go to residue N. The view will not actually change until you mouse click in the 3D display window and hit the right arrow key.


  * view_reset -- will accept the current view (in relation to the present alpha carbon) as the default. This is extremely useful for preparing views for presentation, which compare sidechains or ligands in two homologous structures. For example, if you are interested in a Mg++ ion adjacent to residue 56, you will first use the arrow keys to center on the residue 56 alpha carbon. Then recenter on the Mg++ ion, and rotate to the exact view desired. Typing view_reset will then produce the same view of the Mg++ ion in both structures.


  * view_preset -- revert back to the generic view centered on the alpha carbon.



**METHODOLOGY**

Four atoms of the current residue (N, CA, CB, and C) are used for a Kabsch-style alignment against a reference orientation. For glycines, a hypothetical beta-carbon is modelled (but not shown) based on tetrahedral geometry. A known limitation of the approach is that the alignment is very local, i.e., neighboring residues may not precisely align between structures. 

Adaptation of PyMOL's set_view command to display aligned views of separate structures was suggested by Herb Klei. 

**SUPPORTING MODULE**

The following Python module, support.py, must be placed in the current working directory. Code is based on CCTBX (cctbx.sf.net) and is released under the cctbx license. 
    
    
    import math
    
    #####################################################################
    # This code has been adapted from cctbx.sf.net, file scitbx/matrix.py
    flex = None
    numpy = None
    
    class rec(object):
    
      container_type = tuple
    
      def __init__(self, elems, n):
        assert len(n) == 2
        if (not isinstance(elems, self.container_type)):
          elems = self.container_type(elems)
        assert len(elems) == n[0] * n[1]
        self.elems = elems
        self.n = tuple(n)
    
      def n_rows(self):
        return self.n[0]
    
      def n_columns(self):
        return self.n[1]
    
      def __neg__(self):
        return rec([-e for e in self.elems], self.n)
    
      def __add__(self, other):
        assert self.n == other.n
        a = self.elems
        b = other.elems
        return rec([a[i] + b[i] for i in xrange(len(a))], self.n)
    
      def __sub__(self, other):
        assert self.n == other.n
        a = self.elems
        b = other.elems
        return rec([a[i] - b[i] for i in xrange(len(a))], self.n)
    
      def __mul__(self, other):
        if (not hasattr(other, "elems")):
          if (not isinstance(other, (list, tuple))):
            return rec([x * other for x in self.elems], self.n)
          other = col(other)
        a = self.elems
        ar = self.n_rows()
        ac = self.n_columns()
        b = other.elems
        if (other.n_rows() != ac):
          raise RuntimeError(
            "Incompatible matrices:\n"
            "  self.n:  %s\n"
            "  other.n: %s" % (str(self.n), str(other.n)))
        bc = other.n_columns()
        if (ac == 0):
          # Roy Featherstone, Springer, New York, 2007, p. 53 footnote
          return rec((0,)*(ar*bc), (ar,bc))
        result = []
        for i in xrange(ar):
          for k in xrange(bc):
            s = 0
            for j in xrange(ac):
              s += a[i * ac + j] * b[j * bc + k]
            result.append(s)
        if (ar == bc):
          return sqr(result)
        return rec(result, (ar, bc))
    
      def __rmul__(self, other):
        "scalar * matrix"
        if (isinstance(other, rec)): # work around odd Python 2.2 feature
          return other.__mul__(self)
        return self * other
    
      def transpose_multiply(self, other=None):
        a = self.elems
        ar = self.n_rows()
        ac = self.n_columns()
        if (other is None):
          result = [0] * (ac * ac)
          jac = 0
          for j in xrange(ar):
            ik = 0
            for i in xrange(ac):
              for k in xrange(ac):
                result[ik] += a[jac + i] * a[jac + k]
                ik += 1
            jac += ac
          return sqr(result)
        b = other.elems
        assert other.n_rows() == ar, "Incompatible matrices."
        bc = other.n_columns()
        result = [0] * (ac * bc)
        jac = 0
        jbc = 0
        for j in xrange(ar):
          ik = 0
          for i in xrange(ac):
            for k in xrange(bc):
              result[ik] += a[jac + i] * b[jbc + k]
              ik += 1
          jac += ac
          jbc += bc
        if (ac == bc):
          return sqr(result)
        return rec(result, (ac, bc))
    
      def __div__(self, other):
        return rec([e/other for e in self.elems], self.n)
    
      def __truediv__(self, other):
        return rec([e/other for e in self.elems], self.n)
    
      def __floordiv__(self, other):
        return rec([e//other for e in self.elems], self.n)
    
      def __mod__(self, other):
        return rec([ e % other for e in self.elems], self.n)
    
      def __call__(self, ir, ic):
        return self.elems[ir * self.n_columns() + ic]
    
      def __len__(self):
        return len(self.elems)
    
      def __getitem__(self, i):
        return self.elems[i]
    
      def as_float(self):
        return rec([float(e) for e in self.elems], self.n)
    
      def as_int(self, rounding=True):
        if rounding:
          return rec([int(round(e)) for e in self.elems], self.n)
        else:
          return rec([int(e) for e in self.elems], self.n)
    
      def each_abs(self):
        return rec([abs(e) for e in self.elems], self.n)
    
      def min(self):
        result = None
        for e in self.elems:
          if (result is None or result > e):
            result = e
        return result
    
      def max(self):
        result = None
        for e in self.elems:
          if (result is None or result < e):
            result = e
        return result
    
      def min_index(self):
        result = None
        for i in xrange(len(self.elems)):
          if (result is None or self.elems[result] > self.elems[i]):
            result = i
        return result
    
      def max_index(self):
        result = None
        for i in xrange(len(self.elems)):
          if (result is None or self.elems[result] < self.elems[i]):
            result = i
        return result
    
      def sum(self):
        result = 0
        for e in self.elems:
          result += e
        return result
    
      def product(self):
        result = 1
        for e in self.elems:
          result *= e
        return result
    
      def trace(self):
        assert self.n_rows() == self.n_columns()
        n = self.n_rows()
        result = 0
        for i in xrange(n):
          result += self.elems[i*n+i]
        return result
    
      def norm_sq(self):
        result = 0
        for e in self.elems:
          result += e*e
        return result
    
      def round(self, digits):
        return rec([ round(x, digits) for x in self.elems ], self.n)
    
      def __abs__(self):
        assert self.n_rows() == 1 or self.n_columns() == 1
        return math.sqrt(self.norm_sq())
    
      length_sq = norm_sq # for compatibility with scitbx/vec3.h
      length = __abs__
    
      def normalize(self):
        return self / abs(self)
    
      def dot(self, other=None):
        result = 0
        a = self.elems
        if (other is None):
          for i in xrange(len(a)):
            v = a[i]
            result += v * v
        else:
          assert len(self.elems) == len(other.elems)
          b = other.elems
          for i in xrange(len(a)):
            result += a[i] * b[i]
        return result
    
      def cross(self, other):
        assert self.n in ((3,1), (1,3))
        assert self.n == other.n
        a = self.elems
        b = other.elems
        return rec((
          a[1] * b[2] - b[1] * a[2],
          a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]), self.n)
    
      def is_r3_rotation_matrix_rms(self):
        if (self.n != (3,3)): raise RuntimeError("Not a 3x3 matrix.")
        rtr = self.transpose_multiply()
        return (rtr - identity(n=3)).norm_sq()**0.5
    
      def is_r3_rotation_matrix(self, rms_tolerance=1e-8):
        return self.is_r3_rotation_matrix_rms() < rms_tolerance
    
      def unit_quaternion_as_r3_rotation_matrix(self):
        assert self.n in [(1,4), (4,1)]
        q0,q1,q2,q3 = self.elems
        return sqr((
          2*(q0*q0+q1*q1)-1, 2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2),
          2*(q1*q2+q0*q3),   2*(q0*q0+q2*q2)-1, 2*(q2*q3-q0*q1),
          2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   2*(q0*q0+q3*q3)-1))
    
      def r3_rotation_matrix_as_unit_quaternion(self):
        # Based on work by:
        #   Shepperd (1978), J. Guidance and Control, 1, 223-224.
        #   Sam Buss, http://math.ucsd.edu/~sbuss/MathCG
        #   Robert Hanson, jmol/Jmol/src/org/jmol/util/Quaternion.java
        if (self.n != (3,3)): raise RuntimeError("Not a 3x3 matrix.")
        m00,m01,m02,m10,m11,m12,m20,m21,m22 = self.elems
        trace = m00 + m11 + m22
        if (trace >= 0.5):
          w = (1 + trace)**0.5
          d = w + w
          w *= 0.5
          x = (m21 - m12) / d
          y = (m02 - m20) / d
          z = (m10 - m01) / d
        else:
          if (m00 > m11):
            if (m00 > m22): mx = 0
            else:           mx = 2
          elif (m11 > m22): mx = 1
          else:             mx = 2
          invalid_cutoff = 0.8 # not critical; true value is closer to 0.83
          invalid_message = "Not a r3_rotation matrix."
          if (mx == 0):
            x_sq = 1 + m00 - m11 - m22
            if (x_sq < invalid_cutoff): raise RuntimeError(invalid_message)
            x = x_sq**0.5
            d = x + x
            x *= 0.5
            w = (m21 - m12) / d
            y = (m10 + m01) / d
            z = (m20 + m02) / d
          elif (mx == 1):
            y_sq = 1 + m11 - m00 - m22
            if (y_sq < invalid_cutoff): raise RuntimeError(invalid_message)
            y = y_sq**0.5
            d = y + y
            y *= 0.5
            w = (m02 - m20) / d
            x = (m10 + m01) / d
            z = (m21 + m12) / d
          else:
            z_sq = 1 + m22 - m00 - m11
            if (z_sq < invalid_cutoff): raise RuntimeError(invalid_message)
            z = z_sq**0.5
            d = z + z
            z *= 0.5
            w = (m10 - m01) / d
            x = (m20 + m02) / d
            y = (m21 + m12) / d
        return col((w, x, y, z))
    
      def unit_quaternion_product(self, other):
        assert self.n in [(1,4), (4,1)]
        assert other.n in [(1,4), (4,1)]
        q0,q1,q2,q3 = self.elems
        o0,o1,o2,o3 = other.elems
        return col((
          q0*o0 - q1*o1 - q2*o2 - q3*o3,
          q0*o1 + q1*o0 + q2*o3 - q3*o2,
          q0*o2 - q1*o3 + q2*o0 + q3*o1,
          q0*o3 + q1*o2 - q2*o1 + q3*o0))
    
      def axis_and_angle_as_unit_quaternion(self, angle, deg=False):
        assert self.n in ((3,1), (1,3))
        if (deg): angle *= math.pi/180
        h = angle * 0.5
        c, s = math.cos(h), math.sin(h)
        u,v,w = self.normalize().elems
        return col((c, u*s, v*s, w*s))
    
      def axis_and_angle_as_r3_rotation_matrix(self, angle, deg=False):
        uq = self.axis_and_angle_as_unit_quaternion(angle=angle, deg=deg)
        return uq.unit_quaternion_as_r3_rotation_matrix()
    
      def rt_for_rotation_around_axis_through(self, point, angle, deg=False):
        assert self.n in ((3,1), (1,3))
        assert point.n in ((3,1), (1,3))
        r = (point - self).axis_and_angle_as_r3_rotation_matrix(
          angle=angle, deg=deg)
        return rt((r, self-r*self))
    
      def ortho(self):
        assert self.n in ((3,1), (1,3))
        x, y, z = self.elems
        a, b, c = abs(x), abs(y), abs(z)
        if c <= a and c <= b:
          return col((-y, x, 0))
        if b <= a and b <= c:
          return col((-z, 0, x))
        return col((0, -z, y))
    
      def rotate_around_origin(self, axis, angle, deg=False):
        assert self.n in ((3,1), (1,3))
        assert axis.n == self.n
        if deg: angle *= math.pi/180
        n = axis.normalize()
        x = self
        c, s = math.cos(angle), math.sin(angle)
        return x*c + n*n.dot(x)*(1-c) + n.cross(x)*s
    
      def rotate(self, axis, angle, deg=False):
        import warnings
        warnings.warn(
          message=
            "The .rotate() method has been renamed to .rotate_around_origin()"
            " for clarity. Please update the code calling this method.",
          category=DeprecationWarning,
          stacklevel=2)
        return self.rotate_around_origin(axis=axis, angle=angle, deg=deg)
    
      def vector_to_001_rotation(self,
            sin_angle_is_zero_threshold=1.e-10,
            is_normal_vector_threshold=1.e-10):
        assert self.n in ((3,1), (1,3))
        x,y,c = self.elems
        xxyy = x*x + y*y
        if (abs(xxyy + c*c - 1) > is_normal_vector_threshold):
          raise RuntimeError("self is not a normal vector.")
        s = (xxyy)**0.5
        if (s < sin_angle_is_zero_threshold):
          if (c > 0):
            return sqr((1,0,0,0,1,0,0,0,1))
          return sqr((1,0,0,0,-1,0,0,0,-1))
        us = y
        vs = -x
        u = us / s
        v = vs / s
        oc = 1-c
        return sqr((c + u*u*oc, u*v*oc, vs, u*v*oc, c + v*v*oc, -us, -vs, us, c))
    
      def outer_product(self, other=None):
        if (other is None): other = self
        assert self.n[0] == 1 or self.n[1] == 1
        assert other.n[0] == 1 or other.n[1] == 1
        result = []
        for a in self.elems:
          for b in other.elems:
            result.append(a*b)
        return rec(result, (len(self.elems), len(other.elems)))
    
      def cos_angle(self, other, value_if_undefined=None):
        self_norm_sq = self.norm_sq()
        if (self_norm_sq == 0): return value_if_undefined
        other_norm_sq = other.norm_sq()
        if (other_norm_sq == 0): return value_if_undefined
        d = self_norm_sq * other_norm_sq
        if (d == 0): return value_if_undefined
        return self.dot(other) / math.sqrt(d)
    
      def angle(self, other, value_if_undefined=None, deg=False):
        cos_angle = self.cos_angle(other=other)
        if (cos_angle is None): return value_if_undefined
        result = math.acos(max(-1,min(1,cos_angle)))
        if (deg): result *= 180/math.pi
        return result
    
      def accute_angle(self, other, value_if_undefined=None, deg=False):
        cos_angle = self.cos_angle(other=other)
        if (cos_angle is None): return value_if_undefined
        if (cos_angle < 0): cos_angle *= -1
        result = math.acos(min(1,cos_angle))
        if (deg): result *= 180/math.pi
        return result
    
      def is_square(self):
        return self.n[0] == self.n[1]
    
      def determinant(self):
        assert self.is_square()
        m = self.elems
        n = self.n[0]
        if (n == 1):
          return m[0]
        if (n == 2):
          return m[0]*m[3] - m[1]*m[2]
        if (n == 3):
          return   m[0] * (m[4] * m[8] - m[5] * m[7]) \
                 - m[1] * (m[3] * m[8] - m[5] * m[6]) \
                 + m[2] * (m[3] * m[7] - m[4] * m[6])
        if (flex is not None):
          m = flex.double(m)
          m.resize(flex.grid(self.n))
          return m.matrix_determinant_via_lu()
        return determinant_via_lu(m=self)
    
      def co_factor_matrix_transposed(self):
        n = self.n
        if (n == (0,0)):
          return rec(elems=(), n=n)
        if (n == (1,1)):
          return rec(elems=(1,), n=n)
        m = self.elems
        if (n == (2,2)):
          return rec(elems=(m[3], -m[1], -m[2], m[0]), n=n)
        if (n == (3,3)):
          return rec(elems=(
             m[4] * m[8] - m[5] * m[7],
            -m[1] * m[8] + m[2] * m[7],
             m[1] * m[5] - m[2] * m[4],
            -m[3] * m[8] + m[5] * m[6],
             m[0] * m[8] - m[2] * m[6],
            -m[0] * m[5] + m[2] * m[3],
             m[3] * m[7] - m[4] * m[6],
            -m[0] * m[7] + m[1] * m[6],
             m[0] * m[4] - m[1] * m[3]), n=n)
        assert self.is_square()
        raise RuntimeError("Not implemented.")
    
      def inverse(self):
        assert self.is_square()
        n = self.n
        if (n[0] < 4):
          determinant = self.determinant()
          assert determinant != 0
          return self.co_factor_matrix_transposed() / determinant
        if (flex is not None):
          m = flex.double(self.elems)
          m.resize(flex.grid(n))
          m.matrix_inversion_in_place()
          return rec(elems=m, n=n)
        if (numpy is not None):
          m = numpy.asarray(self.elems)
          m.shape = n
          m = numpy.ravel(numpy.linalg.inv(m))
          return rec(elems=m, n=n)
        return inverse_via_lu(m=self)
    
      def transpose(self):
        elems = []
        for j in xrange(self.n_columns()):
          for i in xrange(self.n_rows()):
            elems.append(self(i,j))
        return rec(elems, (self.n_columns(), self.n_rows()))
    
      def _mathematica_or_matlab_form(self,
            outer_open, outer_close,
            inner_open, inner_close, inner_close_follow,
            label,
            one_row_per_line,
            format,
            prefix):
        nr = self.n_rows()
        nc = self.n_columns()
        s = prefix
        indent = prefix
        if (label):
          s += label + "="
          indent += " " * (len(label) + 1)
        s += outer_open
        if (nc != 0):
          for ir in xrange(nr):
            s += inner_open
            for ic in xrange(nc):
              if (format is None):
                s += str(self(ir, ic))
              else:
                s += format % self(ir, ic)
              if (ic+1 != nc): s += ", "
              elif (ir+1 != nr or len(inner_open) != 0): s += inner_close
            if (ir+1 != nr):
              s += inner_close_follow
              if (one_row_per_line):
                s += "\n"
                s += indent
              s += " "
        return s + outer_close
    
      def mathematica_form(self,
            label="",
            one_row_per_line=False,
            format=None,
            prefix="",
            matrix_form=False):
        result = self._mathematica_or_matlab_form(
          outer_open="{", outer_close="}",
          inner_open="{", inner_close="}", inner_close_follow=",",
          label=label,
          one_row_per_line=one_row_per_line,
          format=format,
          prefix=prefix)
        if matrix_form: result += "//MatrixForm"
        result = result.replace('e', '*^')
        return result
    
      def matlab_form(self,
            label="",
            one_row_per_line=False,
            format=None,
            prefix=""):
        return self._mathematica_or_matlab_form(
          outer_open="[", outer_close="]",
          inner_open="", inner_close=";", inner_close_follow="",
          label=label,
          one_row_per_line=one_row_per_line,
          format=format,
          prefix=prefix)
    
      def __repr__(self):
        n0, n1 = self.n
        e = self.elems
        if (len(e) <= 3):
          e = str(e)
        else:
          e = "(%s, ..., %s)" % (str(e[0]), str(e[-1]))
        return "matrix.rec(elems=%s, n=(%d,%d))" % (e, n0, n1)
    
      def __str__(self):
        return self.mathematica_form(one_row_per_line=True)
    
      def as_list_of_lists(self):
        result = []
        nr,nc = self.n
        for ir in xrange(nr):
          result.append(list(self.elems[ir*nc:(ir+1)*nc]))
        return result
    
      def as_sym_mat3(self):
        assert self.n == (3,3)
        m = self.elems
        return (m[0],m[4],m[8],
                (m[1]+m[3])/2.,
                (m[2]+m[6])/2.,
                (m[5]+m[7])/2.)
    
      def as_mat3(self):
        assert self.n == (3,3)
        return self.elems
    
      def as_flex_double_matrix(self):
        assert flex is not None
        result = flex.double(self.elems)
        result.reshape(flex.grid(self.n))
        return result
    
      def as_flex_int_matrix(self):
        assert flex is not None
        result = flex.int(self.elems)
        result.reshape(flex.grid(self.n))
        return result
    
      def extract_block(self, stop, start=(0,0), step=(1,1)):
        assert 0 <= stop[0] <= self.n[0]
        assert 0 <= stop[1] <= self.n[1]
        i_rows = range(start[0], stop[0], step[0])
        i_colums = range(start[1], stop[1], step[1])
        result = []
        for ir in i_rows:
          for ic in i_colums:
            result.append(self(ir,ic))
        return rec(result, (len(i_rows),len(i_colums)))
    
      def __eq__(self, other):
        if self is other: return True
        if other is None: return False
        if issubclass(type(other), rec):
          return self.elems == other.elems
        for ir in xrange(self.n_rows()):
          for ic in xrange(self.n_columns()):
            if self(ir,ic) != other[ir,ic]: return False
        return True
    
      def resolve_partitions(self):
        nr,nc = self.n
        result_nr = 0
        for ir in xrange(nr):
          part_nr = 0
          for ic in xrange(nc):
            part = self(ir,ic)
            assert isinstance(part, rec)
            if (ic == 0): part_nr = part.n[0]
            else: assert part.n[0] == part_nr
          result_nr += part_nr
        result_nc = 0
        for ic in xrange(nc):
          part_nc = 0
          for ir in xrange(nr):
            part = self(ir,ic)
            if (ir == 0): part_nc = part.n[1]
            else: assert part.n[1] == part_nc
          result_nc += part_nc
        result_elems = [0] * (result_nr * result_nc)
        result_ir = 0
        for ir in xrange(nr):
          result_ic = 0
          for ic in xrange(nc):
            part = self(ir,ic)
            part_nr,part_nc = part.n
            i_part = 0
            for part_ir in xrange(part_nr):
              i_result = (result_ir + part_ir) * result_nc + result_ic
              for part_ic in xrange(part_nc):
                result_elems[i_result + part_ic] = part[i_part]
                i_part += 1
            result_ic += part_nc
          assert result_ic == result_nc
          result_ir += part_nr
        assert result_ir == result_nr
        return rec(elems=result_elems, n=(result_nr, result_nc))
    
    class mutable_rec(rec):
      container_type = list
    
      def __setitem__(self, i, x):
        self.elems[i] = x
    
    class row_mixin(object):
    
      def __init__(self, elems):
        super(row_mixin, self).__init__(elems, (1, len(elems)))
    
    class row(row_mixin, rec): pass
    class mutable_row(row_mixin, mutable_rec): pass
    
    class col_mixin(object):
    
      def __init__(self, elems):
        super(col_mixin, self).__init__(elems, (len(elems), 1))
    
      def random(cls, n, a, b):
        uniform = random.uniform
        return cls([ uniform(a,b) for i in xrange(n) ])
      random = classmethod(random)
    
    class col(col_mixin, rec): pass
    class mutable_col(col_mixin, mutable_rec): pass
    
    class sqr(rec):
    
      def __init__(self, elems):
        l = len(elems)
        n = int(l**(.5) + 0.5)
        assert l == n * n
        rec.__init__(self, elems, (n,n))
    
    ################################################################################
    # This code has been adapted from cctbx.sf.net, file scitbx/matrix/eigensystem.h
    class real_symmetric:
      def __init__(self,symmat3):
        m = symmat3
        self.relative_epsilon = 1.e-10
        self.absolute_epsilon = 0
        
        self.a = [symmat3[0],symmat3[3],symmat3[1],symmat3[4],symmat3[5],symmat3[2]]
        self.vectors_ = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
        self.values_ = [0.,0.,0.];
        self.min_abs_pivot_ = self.real_symmetric_given_lower_triangle(
           self.a,
           3,
           self.vectors_,
           self.values_,
           self.relative_epsilon,
           self.absolute_epsilon);
      def real_symmetric_given_lower_triangle(self,
           a,  # size of memory pointed to by a must be n*(n+1)/2
           n,  
           eigenvectors,
           eigenvalues,
           relative_epsilon,
           absolute_epsilon):
          assert (relative_epsilon >= 0);
          assert (absolute_epsilon >= 0);
          if (n == 0): return 0;
          # The matrix that will hold the results is initially = I.
    
          for x in xrange(n*n):
            eigenvectors[x] = 0.0;
    
          for x in xrange(0,(n*n),(n+1)):
            eigenvectors[x] = 1.0;
    
          # Setup variables
          #std::size_t il, ilq, ilr, im, imq, imr, ind, iq;
          #std::size_t j, k, km, l, ll, lm, lq, m, mm, mq;
          #FloatType am, anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;
          # Initial and final norms (anorm & anrmx).
          anorm=0.0;
          iq=0;
          for i in xrange(n):
            for j in xrange(i+1):
              if (j!=i): anorm+=a[iq]*a[iq];
              iq+=1
          
          anorm=math.sqrt(2.0*anorm);
          anrmx=relative_epsilon*anorm/n;
          if (anrmx < absolute_epsilon): anrmx = absolute_epsilon;
          if (anorm>0.0):
            # Compute threshold and initialise flag.
            thr=anorm;
            while (thr>anrmx): # Compare threshold with final norm
              thr/=n;
              ind=1;
              while (ind):
                ind=0;
                l=0;
                while (l != n-1): # Test for l beyond penultimate column
                  lq=l*(l+1)/2;
                  ll=l+lq;
                  m=l+1;
                  ilq=n*l;
                  while (m != n): # Test for m beyond last column
                    # Compute sin & cos.
                    mq=m*(m+1)/2;
                    lm=l+mq;
                    if (a[lm]*a[lm]>thr*thr):
                      ind=1;
                      mm=m+mq;
                      x=0.5*(a[ll]-a[mm]);
                      denominator=math.sqrt(a[lm]*a[lm]+x*x);
                      assert (denominator != 0);
                      y=-a[lm]/denominator;
                      if (x<0.0): y=-y;
                      sinx=y/math.sqrt(2.0*(1.0+(math.sqrt(1.0-y*y))));
                      sinx2=sinx*sinx;
                      cosx=math.sqrt(1.0-sinx2);
                      cosx2=cosx*cosx;
                      sincs=sinx*cosx;
                      # Rotate l & m columns.
                      imq=n*m;
                      for i in xrange(n):
                        iq=i*(i+1)/2;
                        if (i!=l and i!=m):
                          if (i<m): im=i+mq;
                          else:     im=m+iq;
                          if (i<l): il=i+lq;
                          else:     il=l+iq;
                          x=a[il]*cosx-a[im]*sinx;
                          a[im]=a[il]*sinx+a[im]*cosx;
                          a[il]=x;
                        
                        ilr=ilq+i;
                        imr=imq+i;
                        x = (eigenvectors[ilr]*cosx) \
                          - (eigenvectors[imr]*sinx);
                        eigenvectors[imr] = (eigenvectors[ilr]*sinx) \
                                            + (eigenvectors[imr]*cosx);
                        eigenvectors[ilr] = x;
                      
                      x=2.0*a[lm]*sincs;
                      y=a[ll]*cosx2+a[mm]*sinx2-x;
                      x=a[ll]*sinx2+a[mm]*cosx2+x;
                      a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2);
                      a[ll]=y;
                      a[mm]=x;
                    
                    m+=1;
                  
                  l+=1;
                
              
            
          
          # Sort eigenvalues & eigenvectors in order of descending eigenvalue.
          k=0;
          for i in xrange(n-1):
            im=i;
            km=k;
            am=a[k];
            l=0;
            for j in xrange(n):
              if (j>i and a[l]>am):
                im=j;
                km=l;
                am=a[l];
              
              l+=j+2;
            
            if (im!=i):
              a[km]=a[k];
              a[k]=am;
              l=n*i;
              m=n*im;
              for jj in xrange(n):
                am=eigenvectors[l];
                eigenvectors[l] = eigenvectors[m];
    	    l+=1;
                eigenvectors[m] = am;
    	    m+=1
              
            
            k+=i+2;
          
          # place sorted eigenvalues into the matrix_vector structure
          k = 0
          for j in xrange(n):
            eigenvalues[j]=a[k];
            k+=j+2;
          
          return anrmx;
    
      def vectors(self):
        return matrix.sqr(self.vectors_)
      def values(self):
        return matrix.col(self.values_)
    
    
    class module: pass
    
    eigensystem = module()
    eigensystem.real_symmetric = real_symmetric
    
    matrix = module()
    matrix.sqr = sqr
    matrix.col = col
    matrix.rec = rec
    
    ######################################################################
    #This code is adapted from cctbx.sf.net, file scitbx/math/superpose.py
    def kabsch_rotation(reference_sites, other_sites):
      """
    Kabsch, W. (1976). Acta Cryst. A32, 922-923.
    A solution for the best rotation to relate two sets of vectors
    
    Based on a prototype by Erik McKee and Reetal K. Pai.
    
      """
      assert len(reference_sites) == len(other_sites)
      sts = matrix.sqr(other_sites.transpose_multiply(reference_sites))
      sym_mat3_input = (sts * sts.transpose()).as_sym_mat3()
      eigs = eigensystem.real_symmetric(sym_mat3_input)
      vals = list(eigs.values())
      vecs = list(eigs.vectors())
      a3 = list(matrix.col(vecs[:3]).cross(matrix.col(vecs[3:6])))
      a = matrix.sqr(list(vecs[:6])+a3)
      b = list(a * sts)
      for i in xrange(3):
        d = math.sqrt(math.fabs(vals[i]))
        if (d > 0):
          for j in xrange(3):
            b[i*3+j] /= d
      b3 = list(matrix.col(b[:3]).cross(matrix.col(b[3:6])))
      b = matrix.sqr(b[:6]+b3)
      return b.transpose() * a
    
    #######################################
    #This code is new for this application:
    class residue:
      def __init__(self,chain,resid,atoms):
        self.C = None
        self.CA = None
        self.CB = None
        self.N = None
        for atom in atoms.atom:
          if chain != str(atom.chain): continue
          if str(resid) != str(atom.resi): continue
          if atom.name=="N":
    	self.N = matrix.col((atom.coord[0],atom.coord[1],atom.coord[2]))
          if atom.name=="CA":
            self.residue_type = atom.resn
    	self.CA = matrix.col((atom.coord[0],atom.coord[1],atom.coord[2]))
          if atom.name=="C":
    	self.C = matrix.col((atom.coord[0],atom.coord[1],atom.coord[2]))
          if atom.name=="CB":
    	self.CB = matrix.col((atom.coord[0],atom.coord[1],atom.coord[2]))
        assert self.N != None
        assert self.CA!= None
        assert self.C != None
        if self.CB==None:  #generate a CB position for Glycine
          self.CB = self.constructed_CB()
          
      def constructed_CB(self):
        #refer to the documentation for geometrical construction
        unit_N = (self.N - self.CA).normalize()
        assert abs(unit_N.length() - 1.0) < 0.0001
        unit_C = (self.C - self.CA).normalize()
        on_bisector = (unit_N + unit_C).normalize()
        unit_rotation_axis = (unit_N - unit_C).normalize()
        expected_angle_bisector_to_CB = math.acos (-1./math.sqrt(3.))
        #Use Euler-Rodrigues formula for rotation
        unit_on_CB = self.rotation_formula(on_bisector,unit_rotation_axis,
                                           expected_angle_bisector_to_CB)
        O_to_CB = 1.53 * unit_on_CB
        return self.CA + O_to_CB
        
      def rotation_formula(self,vector,axis,theta):
        return math.cos(theta)*vector + \
               math.sin(theta)*(axis.cross(vector)) + \
    	   (1.-math.cos(theta))*((axis.dot(vector))*axis)
        
      def sites(self):
        diff_N = self.N - self.CA
        diff_C = self.C - self.CA
        diff_CB= self.CB- self.CA
        all = []
        for site in [diff_N,diff_C,diff_CB]:
          for coord in [0,1,2]:
            all.append(site[coord])
        return matrix.rec(all,(3,3))
    
    class std_residue(residue):
      def __init__(self): pass
    
    class view_matrix:
      def __init__(self,getview_output):
        #3x3 rotation matrix which transforms model to camera space
        self.rotmat = matrix.sqr(getview_output[0:9])
        #camera position in model space and relative to the origin of rotation
        self.camera_position = matrix.col(getview_output[9:12])
        #origin of rotation in model space
        self.origin_of_rotation = matrix.col(getview_output[12:15])
        #front plane distance from camera
        self.front_plane = getview_output[15]
        self.back_plane = getview_output[16]
        self.orthoscopic_flag = getview_output[17]
        
      def set_diff_to_CA(self,residue):
        self.difference_to_this_CA = self.origin_of_rotation - residue.CA
        return self
    
    
    class new_set_view:
      def __init__(self,std_view,kr,residue, verbose=True):
        R_prime = (kr.inverse() * std_view.rotmat)
        
        if verbose: print "delta",std_view.origin_of_rotation - residue.CA
        test = residue.CA + (kr.inverse() * std_view.difference_to_this_CA)
        
        if verbose:
          print "set_view ( ",
          for x in R_prime.elems:  print "%10.4f,"%x,
          for x in std_view.camera_position.elems:  print "%10.4f,"%x,
          for x in residue.CA.elems:  print "%10.7f,"%x,
          for x in [std_view.front_plane,std_view.back_plane,std_view.orthoscopic_flag]:
            print "%10.4f,"%x,
          print ")"    
        
        self.result = list(R_prime.elems)+\
                      list(std_view.camera_position.elems)+\
    		  list(test.elems)+\
    		  [std_view.front_plane,std_view.back_plane,std_view.orthoscopic_flag]
    
      def view_matrix(self): return self.result
    

_Author:_ **Nick Sauter**

Retrieved from "[https://pymolwiki.org/index.php?title=Consistent_View/_Map_Inspect&oldid=8166](https://pymolwiki.org/index.php?title=Consistent_View/_Map_Inspect&oldid=8166)"


---

## Contact Surface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script calculates individual or global contact areas between a receptor molecule and a (multimodel) bundle of docked ligand structures. The exact contact surface area values (in Angstrom^2) are printed to the screen and also appended to a file called contactareas.txt. If only a single global contact surface is calculated, a selection named "contact" is created that includes all receptor atoms within 3.9A of any ligand atom to illustrate the _approximate_ contact surface. 

  
The parameters are: 

**receptor** _(string)_

    

    The name of the selection/object representing the receptor protein

**ligand** _(string)_

    

    The name of the selection/object representing the ligand
    Note that this may be another protein!

**states** _(integer)_ , default:0 

    

    Calculate contact surface between the receptor and the first n states of the ligand.
    If states = 0, the script calculates a global contact surface which takes all possible ligand states into account.

  


# Usage
    
    
    contact_surface receptor, ligand, [states=0]
    

  


# The Code
    
    
    #contact_surface v.3.0
    #Copyleft Martin Christen, 2013
    
    from pymol import cmd,stored
    def contact_surface(receptor,ligand,states=0):
    
    	"""
    	AUTHOR
    	Martin Christen
    	
    	DESCRIPTION
    	This script calculates individual or global contact surfaces between a
    	receptor molecule and a bundle of docked ligand structures (which have
    	to be loaded into PyMOL as a multimodel object).
    	
    	The exact contact surface area values (in Angstrom^2) are printed to
    	the screen and also appended to a file called contactareas.txt
    	
    	If only a single global contact surface is calculated, a selection
    	named "contact" is created that includes all receptor atoms within
    	3.9A of any ligand atom.
    	
    	USAGE
    	contact_surface receptor, ligand, [states=0]
    	
    	PARAMETERS
    	
    	receptor (string)
    	The name of the selection/object representing the receptor protein
    	
    	ligand (string)
    	The name of the selection/object representing the ligand.
    	Note that this may be another protein!
    	
    	states (integer)
    	Calculate contact surface between the receptor and the first n states
    	of the ligand. If states = 0 (default), the script calculates a global
    	contact surface which takes  all possible ligand states into account.
    	"""
    	# sanity check the number of states
    	states = abs(int(states))
    	
    	# make sure all atoms within an object occlude one another
    	cmd.flag('ignore','none')
    	
    	# use solvent-accessible surface with high sampling density
    	cmd.set('dot_solvent','1')
    	cmd.set('dot_density','3')
    	
    	#if the 'states' parameter = 0 create a superposition of all ligand states
    	if states == 0:
    		cmd.split_states(ligand)
    		cmd.group('ligandtemp',ligand+"_*")
    		cmd.create(ligand+"_all",'ligandtemp')
    		cmd.delete('ligandtemp')
    		
    		#create complex
    		cmd.create('complextemp',ligand+"_all "+receptor)
    	
    		#measure area
    		ligand_area=cmd.get_area(ligand+"_all")
    		receptor_area=cmd.get_area(receptor)
    		complex_area=cmd.get_area('complextemp')
    		#normalize since the area is counted TWICE (once on receptor and once on ligand)
    		contact_area=((ligand_area + receptor_area) - complex_area) / 2
    		#delete complex
    		cmd.delete('complextemp')
    		
    		#create the contact surface
    		cmd.select('contact',"("+receptor+" and ("+ligand+"_all around 3.9))")
    		
    		#print contact surface area
    		f=open('contactareas.txt','a')
    		print "%s - %s : " % (receptor,ligand),
    		print >>f, "%-s\t%-s\t" % (receptor,ligand),
    		print >>f, "%-s" % (contact_area)
    		print contact_area
    		f.close()
    		print "The GLOBAL contact area between "+receptor+ " and "+ligand+" is (A^2):"
    		print ((ligand_area + receptor_area) - complex_area) / 2
    	
    	#If 'states' <> 0 calculate the contact areas to the first 'states' ligand states.
    	#No individual contact surface objects are created to avoid overloading PyMOL.
    	else:
    		#create an object for each ligand state
    		cmd.split_states(ligand)
    		
    		#sanity check: do not exceed that maximum number of states
    		if states > cmd.count_states(ligand):
    			states = cmd.count_states(ligand)
    		
    		#calculate contact surface area
    		print "The contact areas between "+receptor+" and "+ligand+" [states 1 - "+str(states)+"] are (A^2):"
    		#start looping
    		for s in range(1,states+1):
    			#create complex
    			cmd.create("tmp",ligand,s,1)
    			cmd.create('complextemp',"tmp "+receptor)
    			#measure areas
    			ligand_area=cmd.get_area('tmp')
    			receptor_area=cmd.get_area(receptor)
    			complex_area=cmd.get_area('complextemp')
    			#normalize since the area is counted TWICE (once on receptor and once on ligand)
    			contact_area=((ligand_area + receptor_area) - complex_area)/2
    			#delete temporary files
    			cmd.delete('tmp')
    			cmd.delete(ligand+"_*")
    			cmd.delete('complextemp')
    			#print contact surface area
    			f=open('contactareas.txt','a')
    			print "%s - %s_%-5s: " % (receptor,ligand,s),
    			print >>f, "%-s\t%-s_%-5s\t" % (receptor,ligand,s),
    			print >>f, "%-s" % (contact_area)
    			print contact_area
    			f.close()
    
    cmd.extend("contact_surface",contact_surface)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Contact_Surface&oldid=10991](https://pymolwiki.org/index.php?title=Contact_Surface&oldid=10991)"


---

## CreateSecondaryStructure

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 DESCRIPTION
  * 2 SETUP
  * 3 NOTES / STATUS
  * 4 USAGE
  * 5 EXAMPLES
  * 6 SCRIPTS (CreateSecondaryStructures.py)



### DESCRIPTION

To build a peptide sequence. 

### SETUP

run CreateSecondaryStructure.py 

### NOTES / STATUS

  * Tested on Pymolv1.0, Windows platform
  * Many bugs to be found



### USAGE

  * seqInfo = getTableFromCsvFile("seqInfo.csv")


    
    
    with in the file something like:
    MET,-54,-47
    PRO,-54,-47
    

  * seqInfo = [['MET',-57,-47],['PRO',-57,-47]]
  * createPeptide(seqInfo)



### EXAMPLES
    
    
    seqInfo = [['MET',-57,-47],['PRO',-57,-47]]
    createPeptide(seqInfo)
    

### SCRIPTS (CreateSecondaryStructures.py)

CreateSecondaryStructures.py 
    
    
    ##############################################
    #   Original Author:  Dan Kulp
    #   Date  :  9/8/2005
    #    MOdified by Jurgen F. Doreleijers
    #    For Hamid Eghbalnia               
    #
    #############################################
    # Call in window like : 
    # @C:\Documents and Settings\jurgen.WHELK.000\workspace\Wattos\python\Wattos\Utils\CreateSecondaryStructures.py
    # Next line is a pymol directive
    python
    import urllib
    
    
    # Well I guess one can build a protein with it but the vdw contacts would be horrible.
    # Peptide needs to be at least 2 residues.
    def createPeptide(seqInfo):
        cmd.delete("all")
        # Creates residue TWO
        editor.attach_amino_acid('pk1',seqInfo[1][0]) 
        # Creates residue ONE
        createSS('resi 2', sequence=seqInfo[0][0],terminal='N')
        print "found sequence info for number of residues: ", len(seqInfo)
        for i in range(2,len(seqInfo) ):
            # resn is the residue number of the new residue
            resn = i + 1
            print "Adding residue: ", resn,   seqInfo[i][0]
            # Note that the previous residue is numbered i. 
            resi = 'resi '+`i`
            createSS(resi, sequence=seqInfo[i][0])
        for i in range( len(seqInfo) ):
            resi = 'resi '+`i+1`
    #        print "Setting backbone angles for residue: ", (i+1),   seqInfo[i][0],seqInfo[i][1],seqInfo[i][2]
            set_phipsi(resi,seqInfo[i][1],seqInfo[i][2])
        
    # Create generic secondary structure, based off a selection
    def createSS(sel, sequence='ALA',repeat=1,terminal='C'):
    
        # Set selection
        selection = "%s and name %s" % (sel,terminal)
    
        # Pick atom for editing - interestingly only need to do this for the first addition
        cmd.edit(selection,None,None,None,pkresi=0,pkbond=0)
    
        # Array of residues
        seq = sequence.split(",")
    
        # Get residue numbering .. potential bug here if number is inconsistent.. (Only works at c-terminal)
        resi = int(cmd.get_model(sel).atom[0].resi) + 1
        # Loop and build new residues
        for i in range(1,repeat+1):
            for s in seq:
    #            print "residue[%i]: %s %s" % (i,s,terminal)
                editor.attach_amino_acid('pk1',s)
    
        # Remove extra OXT carboxylate atom (OXT1, OXT2 ?) .. fix as needed
        if terminal == 'C':
            cmd.remove("%s and name OXT" % sel)
        
        
    def set_phipsi(sel,phi,psi):
        # Get atoms from selection
        atoms = cmd.get_model("byres ("+sel+")")
    
        # Loop through atoms in selection        
        for at in atoms.atom:
            if at.name == "N":
                # Check for a null chain id (some PDBs contain this) 
                unit_select = ""
                if not at.chain == "":
                   unit_select = "chain "+str(at.chain)+" and "
        
                try:
                    # Define residue selections     
                    residue_def_prev = unit_select+'resi '+str(int(at.resi)-1)
                    residue_def      = unit_select+'resi '+str(at.resi)        
    #                print "residue_def_prev: [%s]" % residue_def_prev
    #                print "residue_def     : [%s]" % residue_def
                    if at.resn == "PRO":
                        print "Skipping setting phi for PRO"
                    else:
                        old_phi = cmd.get_dihedral(residue_def_prev+' and name C',residue_def+' and name N', residue_def+' and name CA',residue_def+' and name C')        
                        cmd.set_dihedral(          residue_def_prev+' and name C',residue_def+' and name N', residue_def+' and name CA',residue_def+' and name C'      ,phi)
                        print "Changed residue %4s %4s phi: from %6.1f to %6.1f" % (at.resn, at.resi, old_phi, float(phi))        
                except:
                    
                    print "Note skipping set of phi because of error; this is normal for a N-terminal residue"
                try:
                    residue_def      = unit_select+'resi '+str(at.resi)
                    residue_def_next = unit_select+'resi '+str(int(at.resi)+1)
    #                print "residue_def     : [%s]" % residue_def
    #                print "residue_def_next: [%s]" % residue_def_next
                    old_psi = cmd.get_dihedral(residue_def     +' and name N',residue_def+' and name CA',residue_def+' and name C', residue_def_next+' and name N')
                    cmd.set_dihedral(          residue_def     +' and name N',residue_def+' and name CA',residue_def+' and name C', residue_def_next+' and name N',psi)
                    print "Changed residue %4s %4s psi: from %6.1f to %6.1f" % (at.resn, at.resi, old_psi, float(psi))        
                except:
                    print "Note skipping set of psi; this is normal for a C terminal residue"
    
    def getTableFromCsvFile(urlLocation):
      result = []
      r1 = urllib.urlopen(urlLocation)
      data = r1.read()
      r1.close()  
      dataLines = data.split("\n")   
      for dataLine in dataLines:
        if dataLine:
            result.append( dataLine.split(',') )     
      return result
    
    # next line is a pymol directive.
    python end
    
    os.chdir("C:\Documents and Settings\jurgen.WHELK.000\workspace\Wattos\python\Wattos\Utils")
    seqInfo = getTableFromCsvFile("seqInfo.csv")
    createPeptide(seqInfo)
    

Retrieved from "[https://pymolwiki.org/index.php?title=CreateSecondaryStructure&oldid=7667](https://pymolwiki.org/index.php?title=CreateSecondaryStructure&oldid=7667)"


---

## Cyspka

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/cyspka.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cyspka.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 Overview
    * 1.2 Algorithm development
    * 1.3 Correction to article
  * 2 Example of use
  * 3 References



## Introduction

This script is an experimental surface cysteine pKa predictor.  
The script is solely based on the work by: 

_Predicting Reactivities of Protein Surface Cysteines as Part of a Strategy for Selective Multiple Labeling_.  
Maik H. Jacob, Dan Amir, Vladimir Ratner, Eugene Gussakowsky, and Elisha Haas.  
**Biochemistry**. Vol 44, p. 13664-13672, [doi:10.1021/bi051205t](http://dx.doi.org/10.1021/bi051205t)

Questions to the article should be send to [Maik Jacob](http://www.jacobs-university.de/node/1816). 

[![Cyspka.png](/images/2/2b/Cyspka.png)](/index.php/File:Cyspka.png)

[](/index.php/File:Cyspka.png "Enlarge")

### Overview

The authors Jacob _et al_. were able to describe a computational algorithm that could predict the reactivity of surface cysteines. The algorithm was based on reaction rates with Ellmans reagent, Riddles _et al_., on 26 single cysteine mutants of adenylate kinase. The authors could predict the reactivity of the cysteines with a pearson correlation coefficient of 0.92. The algorithm was based on predicting the pKa values of cysteines by a calculation of electrostatic interactions to the backbone and sidechains of the protein and a energetic solvation effect from the number of atom neighbours. The algorithm is different from other pKa algorithms, since it calculates a Boltzmann energy distribution for the rotational states of cysteine. The reaction rate with Ellman's reagent was set proportional to the fraction of negatively charged cysteines, Bulaj _et al_. 

The authors Ratner _et al_., used the prediction of cysteine reactivity to selectively label several two cysteine mutants of adenylate kinase. A double mutant was selected with a high and a low reactive cysteine. A dye was first reacted with the most reactive cysteine and purified. Subsequent the second dye was attached to the low reactive cysteine under unfolding conditions. The strategy is interesting since in meet the double challenge of both site-specificity and double labeling of proteins. 

The prediction algorithm was tested against a popular, (see Sanchez _et al_., Sundaramoorthy _et al_.), web-based pKa prediction program [PROPKA](http://propka.ki.ku.dk/). [ A script](/index.php/Propka "Propka") was developed to to send a structure from PyMOL to the server and fetch and display the calculated pKa values in PyMOL. The adenylate kinase protein was virtual mutated at the 26 positions described by Jacob et al. and showed a lower pearson correlation coefficient of 0.7. In collaboration with Dr. Maik Jacob, the computational algorithm was developed as a general cysteine prediction algorithm for PyMOL. The script is here published at the pymolwiki, but the phase of validating the algorithm against other experimental pKa values has not been finished. The validation phase was hindered by the limited amount of available cysteine pKa data. For example revealed a search in the "[PPD a database of protein ionization constants](http://www.ddg-pharmfac.net/ppd/PPD/pKasearch.htm)" only 12 canditates, where several was dimerized and only 2 candidate cysteines were surface exposed. 

### Algorithm development

The algorithm is based on electrostatic calculations, where some parameters have been fine-tuned.   
The distance from the sulphur atom (SG) of the cysteine to the nearest backbone amide groups and residues with a partial charge, is considered in the electrostatic model.   
The model is including a evalution of Boltzman distribution of the rotation of the SG atom around the CA->CB bond. 

Twenty-six mutants of Escherichia coli adenylate kinase (4AKE) were produced, each containing a single cysteine at the protein surface, and the rates of the reaction with [Ellman's reagent](http://en.wikipedia.org/wiki/Ellman's_reagent) were measured. The reaction rate was set proportional to the pKa, to fine-tune the parameters in the electro static model. 

### Correction to article

There is a type error in equation 6. There is missing a minus "-". The equations should read:  
W M C , S C ( i ) = − ( ∑ W M C ( i ) + ∑ W S C ( i ) ) {\displaystyle W_{MC,SC(i)}=-\left(\sum W_{MC(i)}+\sum W_{SC(i)}\right)} ![{\\displaystyle W_{MC,SC\(i\)}=-\\left\(\\sum W_{MC\(i\)}+\\sum W_{SC\(i\)}\\right\)}](https://wikimedia.org/api/rest_v1/media/math/render/svg/7d3970750d5609b67ea693104acdba3f8e6910f5)

## Example of use

[Escherichia coli adenylate kinase.](http://www.proteopedia.org/wiki/index.php/4ake)   

    
    
    reinitialize
    import cyspka
    
    fetch 4AKE, async=0
    create 4AKE-A, /4AKE//A and not resn HOH
    delete 4AKE
    hide everything
    show cartoon, 4AKE-A
    cyspka 4AKE-A, A, 18
    
    ### You can loop over several residues. 
    loopcyspka 4AKE-A, A, residue=18.25.41-42
    
    ### OR for the original 26 residues. Takes a long time, so not to many at the time.
    #loopcyspka 4AKE-A, A, residue=18.25.41-42.55.73.90.113.162.188-189.203.28.58.75.102.138.142.148.154.169.214.3.24.86.109
    

## References

  * _Ellman's Reagent: 5,5'-Dithiobis(2-nitrobenzoic Acid) a Reexamination_. Peter W. Riddles, Robert L. Blakeley, and Burt Zerner. **Analytical Biochemistry**. Vol 94, p. 75-81, 1979
  * _Ionization Reactivity Relationships for Cysteine Thiols in Polypeptides_. Grzegorz Bulaj, Tanja Kortemme, and David P. Goldenberg. **Biochemistry**. Vol 37, p. 8965-8972, 1998
  * _A General Strategy for Site-Specific Double Labelling of Globular Proteins for Kinetic FRET Studies_. V. Ratner, E. Kahana, M. Eichler, and E. Haas. **Bioconjugate Chem.**. Vol 13, p. 1163-1170, 2002
  * _Prediction of reversibly oxidized protein cysteine thiols using protein structure properties_. Ricardo Sanchez, Megan Riddle, Jongwook Woo, and Jamil Momand. **Protein Science**. Vol 17, p. 473-481, 2008
  * _Predicting protein homocysteinylation targets based on dihedral strain energy and pKa of cysteines_. Elayanambi Sundaramoorthy, Souvik Maiti, Samir K. Brahmachari, and Shantanu Sengupta. **Proteins**. December, p. 1475-1483, 2007 [doi:10.1002/prot.21846](http://dx.doi.org/10.1002/prot.21846)



Retrieved from "[https://pymolwiki.org/index.php?title=Cyspka&oldid=13900](https://pymolwiki.org/index.php?title=Cyspka&oldid=13900)"


---

## Displacementmap

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/displacementmap.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/displacementmap.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Overview
  * 2 Bugs
  * 3 Example
  * 4 Output
  * 5 Text output
  * 6 Example 1
  * 7 References



## Overview

DisplacementMap is made for easy investigations of suitable positions for site-directed mutagenesis of amino residues into cysteines and FRET/EPR pair labelling. 

A Open and Closed form of a protein should be loaded. New objects should be created for the selected asymmetric unit. Parts of the protein should be aligned, leaving the flexible part in two different positions. 

The input is the objects, Open (molecule1) and Closed (molecule2).   
Further is the criteria for selecting which atom the distance should be calculated between. Standard is atom='CA' (atom).   
Then one selects the Förster distance R0 (mindist). This is the minimum distance between the residues. This depends of the selection of the FRET pair and protein at hand. But usually in the range 40 - 80 Angstrom is suitable.   
Then one defines the minimum displacement that is accepted. Usually R0/2 (mindelta).   
The script will find the 5 best (listlength=5) positive and negative distance displacement between the two objects. 

It parses the results back to Pymol, that is standard set to show it as sticks (showsticks='yes').   
If one is looking for a particular residue range for the FRET pair, this can be specified with two input. resi1=24.45-47.86 resi2=100-105.107 resi1 is "from" and resi2 is "to". Individual residues are split by a ".", and ranges are defined with "-".  
In the end, it makes a large data-matrix with all the distances. It also produces a gnuplot file, for easy visualisation. Just drag the .plt file for win gnuplot command window and it plots your datamatrix. 

## Bugs

If the criterion is set to low, the memory gets flooded in the data-matrix file, making the file unreadable. No solutions found yet. 

## Example
    
    
    dispmap(molecule1="NIL", molecule2="NIL", mindist=30.0, mindelta=15.0, resi1=str(0), resi2=str(0), atom='CA', listlength=5, showsticks='yes'):
    

Use of functions 
    
    
    import displacementmap
    dispmap Open5NT, Closed5NT, 40.0, 15.0, resi1=206, resi2=1-512.515-550 
    dispmap Open5NT, Closed5NT, 30.0, 15.0, resi2=1-512.515-550, atom=CA, listlength=10
    

## Output

Suggestions are created in pymol, and gnuplot file is created for easy visualisation of pair data-matrix and the general backbone displacement. 

  * [![O5NT-C5NT-CA-dist.png](/images/1/1a/O5NT-C5NT-CA-dist.png)](/index.php/File:O5NT-C5NT-CA-dist.png)

  * [![O5NT-C5NT-CA-dist-menu.png](/images/c/cd/O5NT-C5NT-CA-dist-menu.png)](/index.php/File:O5NT-C5NT-CA-dist-menu.png)

  * [![O5NT-C5NT-CA-dist-gnuplot.png](/images/6/66/O5NT-C5NT-CA-dist-gnuplot.png)](/index.php/File:O5NT-C5NT-CA-dist-gnuplot.png)

  * [![O5NT-C5NT-CA-dist-backbone.png](/images/c/ce/O5NT-C5NT-CA-dist-backbone.png)](/index.php/File:O5NT-C5NT-CA-dist-backbone.png)




## Text output

In the data-matrix.txt file, you find the best suggestions 
    
    
    # Input 1: Open5NT  and Input 2: Closed5NT
    # Find for: CA  with min. residue-residue dist: 30.0 Angstrom
    # Looking for min. displacement dist: 15.0 Angstrom
    # I give nr# suggestions: 5, and do I show sticks in pymol?: yes
    # I look for suggestions in the range: ([0]=>means all residues)
    # for Input 1: ['0'] and for Input 2: ['0']
    # Mutation info is BLOSUM62 log-odds likelihood score and PAM250 is probability in % for evolutionary distance
    ###########################################################################################################
    # Max Negative and positive distances                                       #       Mutation info         #
    ###########################################################################################################
    # Obj.1   Obj.2   Delta   Op-Op Cl-Cl # Obj.1   Obj.2   Delta   Op-Op Cl-Cl # Res.1  Res.2 # Res.1  Res.2 #
    # Res.1   Res.2   -Dist   Dist  Dist  # Res.1   Res.2   +Dist   Dist  Dist  # B62/PAM250%  # B62/PAM250%  #
    ###########################################################################################################
    # PRO241  ASP456  -25.7   59.1   33.4 # PRO274  PRO513   26.8   31.2   58.0 # -3/ 2 -3/ 1  # -3/ 2 -3/ 2  #
    # LYS197  ASP456  -25.6   57.3   31.7 # THR273  PRO513   26.1   31.6   57.7 # -1/ 1 -3/ 1  # -1/ 2 -3/ 2  #
    # PRO513  ASP456  -25.4   32.4    7.0 # PRO274  GLY514   24.8   32.9   57.6 # -3/ 2 -3/ 1  # -3/ 2 -3/ 2  #
    # LEU198  ASP456  -25.3   59.0   33.7 # PRO274  LYS512   24.7   30.3   55.0 # -1/ 1 -3/ 1  # -3/ 2 -1/ 1  #
    # GLN201  ASP456  -25.2   62.8   37.6 # ASN311  PRO513   24.7   35.6   60.3 # -3/ 1 -3/ 1  # -3/ 1 -3/ 2  #
    

The script also automatically make the gnuplot plot file (.plt), with all the defined variables, for easy visualisation of the data-matrix.txt and the backbone displacement. 

## Example 1

Download: [examples/displacementmap_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/displacementmap_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/displacementmap_1.pml>" highlight="python" />

## References

For EPR considerations   
_Conformation of T4 Lysozyme in Solution. Hinge-Bending Motion and the Substrate-Induced Conformational Transition Studied by Site-Directed Spin Labeling_   
Hassane S. Mchaourab, Kyoung Joon Oh, Celia J. Fang, and Wayne L. Hubbell   
_Biochemistry_ **1997** , 36, 307-316 

_Probing Single-Molecule T4 Lysozyme Conformational Dynamics by Intramolecular Fluorescence Energy Transfer_   
Yu Chen, Dehong Hu, Erich R. Vorpagel, and H. Peter Lu   
J. Phys. Chem. B **2003** , 107, 7947-7956 

For FRET pair selection and considerations   
_Fluorescent probes and bioconjugation chemistries for single-molecule fluorescence analysis of biomolecules_   
Achillefs N. Kapanidisa and Shimon Weiss   
_Journal of chemical physics_ VOLUME 117, Number 24 22 **December 2002**

**For inspiration to DisplacementMap. Fig: 6, Difference-distance matrix for the difference in CA-CA distances.**   
_Structure of a Hinge-bending Bacteriophage T4 Lysozyme Mutant, Ile3 - > Pro_   
M. M. Dixon, H. Nicholsont, L. Shewchuk W. A. Baase and B. W. Matthews1   
J. Mol. Biol. (**1992**) 227. 917-933 

Retrieved from "[https://pymolwiki.org/index.php?title=Displacementmap&oldid=13901](https://pymolwiki.org/index.php?title=Displacementmap&oldid=13901)"


---

## DistancesRH

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [figshare](http://dx.doi.org/10.6084/m9.figshare.1031580)  
Author(s)  | [Pietro Gatti-Lafranconi](/index.php/User:PietroGattiLafranconi "User:PietroGattiLafranconi")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Examples
  * 6 Notes, Further Developments and Requests
  * 7 Updates
  * 8 The Code



## Introduction

Given an object, a reference amino acid and a radius, this scripts calculates pairwise distances between atom(s) in the reference and all atoms that fall within the given range. 

The script is intended to provide a range of useful tools while requiring minimal coding skills (i.e. it is not optimised for efficiency). 

By changing input parameters it is possible to: 

  * choose the output (none, on screen, on file)
  * measure distances from one single atom
  * display distance objects in pymol
  * restrict distances between hydrogen-bond forming atoms only



  
The script also runs a number of controls to verify that inputs and outputs exists/are suitable. 

  


## Usage
    
    
    distancesRH obj, ref, dist, [chain, [output S/P/N, [show Y/N, [Hbonds Y/N, [aname]]]]]
    

  


## Required Arguments

  * **obj** = object name
  * **ref** = reference amino acid id
  * **dist** = radius in Angstroms



  


## Optional Arguments

  * _chain_ = selected chain
  * _output_ = accepts S (screen), P (print) or N (none), default=S
  * _show_ = shows distances in pymol window, default=N
  * _Hbonds_ = restricts results to atoms that can form Hbonds, default=N
  * _aname_ = selects a specific atom in reference object and to calculate distances from one individual atom



  


## Examples

**example #1**
    
    
    distancesRH 2WHH, 25, 3.2, show=Y
    

Will show distances between all atoms of residue 25 in 2WHH and any atom within 3.2 Å and return the result on screen (figure 1, left). 
    
    
    distancesRH 2WHH, 25, 3.2, show=Y, aname=OD1, output=P
    

Same result as before, but only distances from atom OD1 are measured and shown. Results will be printed in 'distancesRH_25-OD1.txt' (figure 1, right). 

[![example #1](/images/3/37/DistancesRH_ex1.png)](/index.php/File:DistancesRH_ex1.png "example #1")

  
**example #2**
    
    
    distancesRH 1EFA, 901, 3.5, show=Y
    

Will show distances between all atoms in resi 901 of 1EFA and all atoms within 3.5 Å (figure 2, left). 
    
    
    distancesRH 1EFA, 901, 3.5, show=Y, Hbonds=Y
    

Distances will now be calculated between atoms that can form H bonds only (figure 2, right). 

[![example #2](/images/b/b6/DistancesRH_ex2.png)](/index.php/File:DistancesRH_ex2.png "example #2")

  
**example #3**
    
    
    distancesRH 1EFA, 18, 4.7, chain=B, show=Y, aname=NE2
    

Distances between atom NE2 of residue 19 in chain B of 1EFA and all atoms within 4.7 Å (all of which are in the DNA fragment of chain D, figure 3). 

[![example #3](/images/7/71/DistancesRH_ex3.png)](/index.php/File:DistancesRH_ex3.png "example #3")

  


## Notes, Further Developments and Requests

  * This script works with visible objects (handling of multiple states will hopefully be implemented soon)
  * As it is intended to be user-friendly more than fast, this script is not suitable to calculate distances between large objects. Use [Quick_dist](/index.php/Quick_dist "Quick dist") instead.



  


## Updates

  * **11/01/2013** \- edit #1: _distances measured between different chains_
  * **13/01/2013** \- edit #2: _print output optimsed_



  


## The Code

Copy the following text and save it as distancesRH.py 
    
    
    """
    DESCRIPTION
    Given an object, a reference amino acid and a radius, this scripts calculates pairwise
    distances between all atoms in the reference and all atoms that fall within the given range.
    Distances are either returned on screen (default) or printed on file.
    Distances can be calculated from a single atom only.
    Optionally, only atoms that can form H-bonds can be returned (but check distances!)
    
    More information at: PymolWiki
    http://pymolwiki.org/index.php/DistancesRH
    
    AUTHOR
    Pietro Gatti-Lafranconi, 2012
    Please inform me if you use/improve/like/dislike/publish with this script.
    CC BY-NC-SA
    """
    
    from pymol import cmd, stored, math
    	
    def distancesRH (obj, ref, dist, chain=0, output="S", show="N", Hbonds="N", aname="*"):
    	"""	
    	usage: distancesRH obj, ref, dist, [chain, [output S/P/N, [show Y/N, [Hbonds Y/N, [aname]]]]]
    	
    	obj: object name, ref:reference aa, dist: radius in A, chain: selected chain
    	output: accepts S (screen), P (print) or N (none), default=S
    	show: shows distances in pymol window, default=N
    	Hbonds: restricts results to atoms that can form Hbonds, default=N
    	aname: name of a specific atom in ref (to calculate distances from one atom only)
    		
    	example: distancesRH 2WHH, 25, 3.2, [B, [S, [N, [Y, [CB]]]]]
    	"""
    	print ""
    	#delete previous selections and displayed distances
    	cmd.delete ("dist*")
    	cmd.delete ("Atoms")
    	cmd.delete ("Residues")
    	cmd.delete ("Reference")
    	cmd.delete ("Hbonds")
    	cmd.delete ("RefAtom")
    	#creates and names selections	
    	if chain==0: cen=obj
    	else: cen=obj+" and chain "+chain
    	refA=" and resi "+ref+" and name "+aname
    	cmd.select("Reference", "byres "+cen+refA)
    	m1=cmd.get_model ("Reference")
    	if aname!="*":
    		cmd.select("RefAtom", cen+refA)
    		m1=cmd.get_model ("RefAtom")
    	cmd.select("Atoms", cen+refA+" around "+str(dist)+" and not resi "+ref+" and "+obj)
    	cmd.show("nb_spheres", "(Atoms or Reference) and (resn HOH or HET)")
    	cmd.select("Residues", "byres Atoms")
    	m2=cmd.get_model("Atoms and visible")	
    	if Hbonds=="Y":
    		cmd.select("__Hbonds1", "Reference and not symbol c")
    		if aname!="*":
    			cmd.select("__Hbonds1", "RefAtom and not symbol c")
    		cmd.select("__Hbonds2", "Atoms and not symbol c")
    		cmd.select("Hbonds", "__Hbonds1 or __Hbonds2")
    		cmd.show("lines", "Hbonds")
    		m1=cmd.get_model ("__Hbonds1")
    		m2=cmd.get_model ("__Hbonds2 and visible")	
    	
    	#controllers	
    	if (not (output=="S" or output=="N" or output=="P")) or (not (show=="Y" or show=="N")) or (not (Hbonds=="Y" or Hbonds=="N")):
    		print "Malformed input, please try again"
    		return
    	abort=1
    	abort2=0
    	mc=cmd.get_model(cen)
    	for c0 in range(len(mc.atom)):
    		if mc.atom[c0].resi==ref: abort=0
    	mc2=cmd.get_model(cen+refA)
    	if aname!="*":
    		abort2=1
    		for c00 in range (len(mc2.atom)):	
    			if mc2.atom[c00].name==aname: abort2=0
    	if Hbonds=="Y":
    		if aname[0]=="C":
    			print "Warning: no atoms in the selection can form H-bonds"
    			return
    	if abort==1: print "Warning: no such residue in the reference object"
    	if abort2==1: 
    		if abort==0: print "Warning: no such atom in the reference object"
    	if abort==1 or abort2==1: return
    
    	#measures distances
    	distance2=[]
    	distances=[]
    	screenOP=""	
    	fileOP=""
    	for c1 in range(len(m1.atom)):
    		for c2 in range(len(m2.atom)):
    			if math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord,m2.atom[c2].coord))))<float(dist):
    				distance=cmd.distance (obj+"//"+m1.atom[c1].chain+"/"+m1.atom[c1].resi+"/"+m1.atom[c1].name, obj+"//"+m2.atom[c2].chain+"/"+m2.atom[c2].resi+"/"+m2.atom[c2].name)
    				distances.append(distance)
    				screenOP+="%s/%s/%s/%s to %s/%s/%s/%s: %.3f\n" % (m1.atom[c1].chain,m1.atom[c1].resn,m1.atom[c1].resi,m1.atom[c1].name,m2.atom[c2].chain,m2.atom[c2].resn,m2.atom[c2].resi,m2.atom[c2].name, distance)
    				fileOP+="%s/%s/%s/%s\t\t\t%s/%s/%s/%s\t\t\t\t%.5f\n"%  (m1.atom[c1].chain,m1.atom[c1].resn,m1.atom[c1].resi,m1.atom[c1].name,m2.atom[c2].chain,m2.atom[c2].resn,m2.atom[c2].resi,m2.atom[c2].name, distance)	
    	#last controller
    	if len(distances)==0:
    		if abort2==0:
    			print "Warning: no atoms within the selected range"
    		return
    
    	#data outputs
    	if output=="N": print "Data recorded, but not printed"
    	if output=="S":
    		print "Distance of atoms in "+obj+"/"+ref+"/"+aname+" to all atoms within "+dist+" Angstroms"
    		if Hbonds=="Y": print "Distances restricted to potential H-bond pairs"
    		print screenOP	
    	if output=="P":
    		if aname=="*": aname="ALL"
    		print 'data printed in "distances_rH_'+ref+'-'+aname+'.txt" file'
    		f=open('distances_rH_'+ref+'-'+aname+'.txt','w')
    		f.write("Pairwise distances betweet atoms in %s and all atoms within %s Angstroms\n" % (obj+"/"+ref+"/"+aname, dist))
    		if Hbonds=="Y": f.write("Distances restricted to potential H-bond pairs\n")
    		f.write("Atom in reference\tAtom in destination\t\tdistance\n")
    		f.write(fileOP)
    		f.close()
    	print "All done! Good luck with your data"
    
    	#graphical output
    	if show=="N": cmd.delete ("dist*")
    	if show=="Y":
    		cmd.show("lines", "Reference")
    		cmd.show("lines", "Residues")
    				
    cmd.extend("distancesRH", distancesRH);
    

Retrieved from "[https://pymolwiki.org/index.php?title=DistancesRH&oldid=11567](https://pymolwiki.org/index.php?title=DistancesRH&oldid=11567)"


---

## Draw Protein Dimensions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/Draw_Protein_Dimensions.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/Draw_Protein_Dimensions.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw the dimensions of a protein based on an Inertia Axis Aligned Bounding Box (IABB). 

  
The idea behind this script is to calculate an approximate minimal bounding box (MBB) to extract the cell dimensions of a protein. To calculate the MBB is not trivial and usually the Axis Aligned Bounding Box (AABB) does not show up the real dimensions of the protein. This script calculates the inertia tensor of the object, extract the eigenvectors and use them to rotate the molecule (using as rotation matrix the transpose of the eigenvector matrix). The result is a molecule oriented with the inertia axis aligned with the cartesian axis. A new Bounding Box is calculated, which is called Inertia Axis Aligned Bounding Box (IABB), with a volume always lower than the AABB volume, and in many cases may correspond with the MBB. 

As always with these type of things, you have to use at your own risk. I did not try all the possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado (at) gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run path-to-the-script/Draw_Protein_Dimensions.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_Protein_Dimensions _selection_**

This will draw the cell dimensions of your selection based on a IABB box. It also generates the IABB box and the inertia axis, you just need to do "show cgo" to display them. 

You could also try: 

**draw_BB _selection_**

This will draw the AABB and IABB boxes with their cell dimensions and show up their volumes, you can compare them. 

  


# Examples
    
    
    # download the source and save as Draw_Protein_Dimensions.py
    run Draw_Protein_Dimensions.py
    fetch 2vak
    # calculate the dimensions of the full ASU
    draw_Protein_Dimensions 2vak
    

  * [![Dimensions of 2vak ASU](/images/7/77/2vak_ASU.png)](/index.php/File:2vak_ASU.png "Dimensions of 2vak ASU")

Dimensions of 2vak ASU 



    
    
    # you can extract only one chain and calculates the dimensions
    draw_Protein_Dimensions obj01
    

  * [![Dimensions of one protomer \(chain A - 2vak\)](/images/e/ed/2vak_A.png)](/index.php/File:2vak_A.png "Dimensions of one protomer \(chain A - 2vak\)")

Dimensions of one protomer (chain A - 2vak) 



    
    
    # You can also draw the Bounding boxes (AABB and IABB) to compare them.
    fetch 2vak
    draw_BB 2vak
    

  


  * [![Axis-aligned bounding box](/images/4/46/2vak_AABB.png)](/index.php/File:2vak_AABB.png "Axis-aligned bounding box")

Axis-aligned bounding box 

  * [![Inertia-axis-aligned bounding box](/images/4/41/2vak_IABB.png)](/index.php/File:2vak_IABB.png "Inertia-axis-aligned bounding box")

Inertia-axis-aligned bounding box 




Retrieved from "[https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902](https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902)"


---

## Dssp

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.editing](http://pymol.org/psicoredirect.php?psico.editing)  
  
dssp is a wrapper for the popular [DSSP](http://swift.cmbi.ru.nl/gv/dssp/) program, which computes secondary structure. The command updates PyMOL's **ss** atom property. 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Example
  * 5 See Also



## Installation

The dssp command is available from the [psico](/index.php/Psico "Psico") package and requires the [dssp](http://swift.cmbi.ru.nl/gv/dssp/) binary (or mkdssp, dsspcmbi, dssp-2). 

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c salilab -c speleo3 dssp 
    

## Usage
    
    
    dssp [ selection [, exe [, raw [, state ]]]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **exe** = str: path to dssp binary {default: search $PATH for dsspcmbi, dssp, dssp-2, mkdssp}
  * **raw** = str: [atom property](/index.php/Iterate#Exposed_Variables "Iterate") to assign the dssp secondary structure type to. This type will also be translated to PyMOL's H/S/L types and assiged to the **ss** property {default: }
  * **state** = int: object state {default: -1 (current state)}



## Example
    
    
    import psico.editing
    fetch 1ubq, async=0
    
    dssp all, raw=custom
    
    color gray
    color red, custom H
    color orange, custom G
    color yellow, custom E
    color wheat, custom B
    color forest, custom T
    color green, custom S
    set cartoon_discrete_colors, 1
    

## See Also

  * [DSSP Stride](/index.php/DSSP_Stride "DSSP Stride") (plugin with GUI)
  * [dss](/index.php/Dss "Dss")



Retrieved from "[https://pymolwiki.org/index.php?title=Dssp&oldid=12651](https://pymolwiki.org/index.php?title=Dssp&oldid=12651)"


---

## Dssr block

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [script/dssr_block.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/script/dssr_block.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![Trna.png](/images/d/dd/Trna.png)](/index.php/File:Trna.png)

[](/index.php/File:Trna.png "Enlarge")

This script adds the dssr_block command, which is a simple wrapper for the **DSSR** program and can create "block" shaped cartoons for nucleic acid bases and base pairs. 

_Requires the**x3dna-dssr** executable, obtainable from <http://x3dna.org>_

Instead of installing **x3dna-dssr** and the plugin locally, you can also use the web server at <http://skmatic.x3dna.org/>

See also the blog post by DSSR author Xiang-Jun Lu: <http://x3dna.org/highlights/dssr-base-blocks-in-pymol-interactively>

## Contents

  * 1 Recommended Setup
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 See Also



## Recommended Setup

Place the **x3dna-dssr** executable into **$PATH**. Examples: 

  * Linux/Mac: **/usr/bin/x3dna-dssr**
  * Windows: **C:\Program Files\PyMOL\PyMOL\x3dna-dssr.exe**



The script can be [run](/index.php/Run "Run"), imported as a Python module, or installed with the [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager"). 

## Usage
    
    
    dssr_block [ selection [, state [, block_file [, block_depth
           [, block_color [, name [, exe ]]]]]]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **state** = int: object state (0 for all states) {default: -1, current state}
  * **block_file** = face|edge|wc|equal|minor|gray: Corresponds to the **\--block-file** option (see [DSSR manual](http://x3dna.org/)). Values can be combined, e.g. "wc-minor". {default: face}
  * **block_depth** = float: thickness of rectangular blocks {default: 0.5}
  * **block_color** = str: Corresponds to the **\--block-color** option (new in DSSR v1.5.2) {default: }
  * **name** = str: name of new CGO object {default: dssr_block##}
  * **exe** = str: path to "x3dna-dssr" executable {default: x3dna-dssr}



## Examples

Combining DSSR block representation with regular PyMOL cartoons: 
    
    
    fetch 1ehz, async=0
    as cartoon
    set cartoon_ladder_radius, 0.1
    set cartoon_ladder_color, gray
    set cartoon_nucleic_acid_mode, 1
    dssr_block
    

Joined base-pair blocks (block_file=wc): 
    
    
    fetch 1ehz, async=0
    dssr_block block_file=wc
    

Multi-state Example: 
    
    
    fetch 2n2d, async=0
    dssr_block 2n2d, 0
    set all_states
    

Custom coloring: 
    
    
    fetch 1msy, async=0
    dssr_block block_color=N red | minor 0.9 | major yellow
    

## See Also

  * [3DNA](/index.php/3DNA "3DNA") (old, obsoleted by x3dna-dssr)
  * [Overview of nucleic acid cartoons](/index.php/Overview_of_nucleic_acid_cartoons "Overview of nucleic acid cartoons")
  * <http://skmatic.x3dna.org/>



Retrieved from "[https://pymolwiki.org/index.php?title=Dssr_block&oldid=13940](https://pymolwiki.org/index.php?title=Dssr_block&oldid=13940)"


---

## DynoPlot

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/dynoplot.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/dynoplot.py)  
Author(s)  | [Dan Kulp](/index.php/User:Tmwsiy "User:Tmwsiy")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 IMAGES
    * 1.2 SETUP
    * 1.3 NOTES / STATUS
    * 1.4 USAGE
    * 1.5 EXAMPLES



## Introduction

This script was setup to do generic plotting, that is given a set of data and axis labels it would create a plot. Initially, I had it setup to draw the plot directly in the PyMol window (allowing for both 2D and 3D style plots), but because I couldn't figure out how to billboard CGO objects (Warren told me at the time that it couldn't be done) I took a different approach. The plot now exists in it's own window and can only do 2D plots. It is however interactive. I only have here a Rama.(phi,psi) plot, but the code can be easily extended to other types of data. For instance, I had this working for an energy vs distance data that I had generated by another script. 

This script will create a Phi vs Psi(Ramachandran) plot of the selection given. The plot will display data points which can be dragged around Phi,Psi space with the corresponding residue's Phi,Psi angles changing in the structure (PyMol window). 

### IMAGES

  * [![Initial Ramachandran plot of 1ENV](/images/e/e0/RamaPlotInitComposite.png)](/index.php/File:RamaPlotInitComposite.png "Initial Ramachandran plot of 1ENV")

Initial Ramachandran plot of 1ENV 

  * [![Modified pose and plot of 1ENV](/images/c/c7/RamaPlotBentComposite.png)](/index.php/File:RamaPlotBentComposite.png "Modified pose and plot of 1ENV")

Modified pose and plot of 1ENV 




### SETUP

Install from the plugins menu with _Plugin > Manage Plugins > Install ..._ or just [run](/index.php/Run "Run") the script. 

### NOTES / STATUS

  * Tested on Linux, PyMol version 1.4
  * Left, Right mouse buttons do different things; Right = identify data point, Left = drag data point around
  * Post comments/questions or send them to: dwkulp@mail.med.upenn.edu



### USAGE
    
    
    rama [ sel [, name [, symbols [, filename ]]]]
    

### EXAMPLES
    
    
    fetch 1ENV, async=0 # (download it or use the PDB loader plugin)
    select sel01, resi 129-136
    rama sel01
    rock   # the object needs to be moving in order for the angles to be updated.
    

Don't create callback object, use symbols by secondary structure and dump canvas as postscript file: 
    
    
    fetch 2x19, async=0
    color yellow, chain A
    color forest, chain B
    rama polymer, none, ss, /tmp/canvasdump.ps
    rama ss H,    none, aa, /tmp/canvasdump_helix.ps
    rama ss S,    none, aa, /tmp/canvasdump_sheet.ps
    

Retrieved from "[https://pymolwiki.org/index.php?title=DynoPlot&oldid=10435](https://pymolwiki.org/index.php?title=DynoPlot&oldid=10435)"


---

## Elbow angle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/elbow_angle.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/elbow_angle.py)  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  


## Contents

  * 1 Introduction
  * 2 Syntax
  * 3 Examples
    * 3.1 Basic usage
    * 3.2 Included Example
  * 4 Implementation



## Introduction

This script allows you to calculate the elbow angle of an antibody Fab fragment object and optionally draw a graphical representation of the vectors used to calculate the elbow angle. 

  


  


## Syntax
    
    
    elbow_angle object [, light=L [, heavy=H [, limit_l=107 [, limit_h=113 [, draw=0 ]]]]]
    

The first argument should be a PyMOL object. You can calculate the elbow angles of multiple Fabs from the same object (one at a time) by specifying the chains manually for each Fab. 

The `light` and `heavy` arguments are for PDB chain IDs, and `limit_l` and `limit_h` are the last residue of the light and heavy chain variable domains, respectively. (For Kabat-numbered structures, these limits will be 107 and 113, respectively. For structures with different numbering schemes, the limits can be estimated by visual inspection of the PDB file.) Setting `draw=1` will draw the "dumbbells" seen in the images below. 

## Examples

### Basic usage
    
    
    # load an antibody Fab fragment from the PDB
    fetch 3ghe, async=0
    # calculate the elbow angle and draw the vectors
    elbow_angle 3ghe, draw=1
    

This will produce something like the first image below. 

Note that if you don't specify the light/heavy chain IDs or the variable domain limit residue numbers, the default values (L/H and 107/113, respectively) will be used. If your antibody is not Kabat or Chothia numbered, or has different chain names, this will result in an incorrect value or, in the case of wrong chain IDs, may cause the script to fail entirely due to an empty selection. Have a look at Dr. Andrew Martin's [Abnum](http://www.bioinf.org.uk/abs/abnum/) for more information on antibody numbering. 

  * [![Fab fragment 3ghe shown with draw=1.](/images/1/11/Elbow_angle_3ghe.png)](/index.php/File:Elbow_angle_3ghe.png "Fab fragment 3ghe shown with draw=1.")

Fab fragment 3ghe shown with draw=1. 

  * [![5 PDB examples from Stanfield, et al., JMB 2006, shown in the same orientations as in Figure 1 of that paper.](/images/9/95/Stanfield.png)](/index.php/File:Stanfield.png "5 PDB examples from Stanfield, et al., JMB 2006, shown in the same orientations as in Figure 1 of that paper.")

5 PDB examples from Stanfield, et al., JMB 2006, shown in the same orientations as in Figure 1 of that paper. 




The black "dumbbells" pass through the centers of mass of the variable and constant domains of each Fab (as determined using [com.py](/index.php/Com "Com")). The green and red dumbbells denote the residues used to split the variable and constant domains, with a green ball for the light chain, and a red ball for the heavy chain. 

### Included Example

There is also an example .pml file in the Pymol-script-repo/examples directory which can be run by the following: 
    
    
    import ex
    ex elbow_angle
    

This will run the following .pml file: 

Type  | [PyMOL Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [examples/elbow_angle.pml](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/examples/elbow_angle.pml)  
Author(s)  | [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/elbow_angle.pml>" highlight="python" />

and will produce the second image above. 

## Implementation

The elbow angle is the angle between the pseudo-twofold axes between the light and heavy chain variable and constant domains, respectively. The rotation matrix to superpose VL onto VH and CL onto CH are calculated and the vectors corresponding to the rotation axes are determined (using Christoph Gohlke's [transformations.py](/index.php/Transformations "Transformations")). The elbow angle is the obtuse angle obtained by taking the arccos of the dot product of the two vectors. For consistency, the elbow angle is reported as less than 180° when the cross product of the Variable and Constant domain rotation axis vectors ( **V** ⨯ **C** ) is a vector pointing the same direction as that from the limit_h residue alpha carbon atom to the limit_l alpha carbon atom. 

Retrieved from "[https://pymolwiki.org/index.php?title=Elbow_angle&oldid=13905](https://pymolwiki.org/index.php?title=Elbow_angle&oldid=13905)"


---

## FindSurfaceCharge

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/findSurfaceCharge.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findSurfaceCharge.py)  
Author(s)  | [Teddy Warner](/index.php/User:TeddyWarner "User:TeddyWarner")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Drawing upon the [findSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") script, the findSurfaceCharge script will identify and output a list of all charged residues on the surface of a selectionand calculates the ionization state of a protein at a given pH. The charge can be calculated for either a folded or denatured protein. This function is intended to be used to give buffer conditions for mass spectrometry. 

# Usage
    
    
       findSurfaceCharge [selection, [pH, [folded ,[cutoff]]]]
    

# Arguments

  * **selection** = str: The object or selection for which to find exposed residues {default: all}
  * **pH** = float: The pH to calculate the surface charge at {default: 7.0}
  * **folded** = bool: Whether the program should calculate the charge of a folded (True) or denatured (False) protein.
  * **cutoff** = float: The cutoff in square Angstroms that defines exposed or not. Those atoms with > cutoff Ang^2 exposed will be considered _exposed_ {default: 2.5 Ang^2}



# Examples

  * [![Result of 4FIX.pdb at pH 7.0.](/images/2/2a/SurfaceCharge.PNG)](/index.php/File:SurfaceCharge.PNG "Result of 4FIX.pdb at pH 7.0.")

Result of 4FIX.pdb at pH 7.0. 



    
    
    run findSurfaceResiduesListCharged.py
    fetch 4FIX
    
    findSurfaceResiduesListCharged
    
    # see how pH changes the protein surface charge:
    findSurfaceCharge("4fix",7.0)
        Exposed charged residues: 
            ERREDRKEE...
        The expected surface charge of 4fix at pH 7.0 is: +3.24
    
    findSurfaceCharge("4fix",7.0,False)
        Charged residues: 
            HHHHHHRHERREDRKEE...
        The expected charge of denatured 4fix at pH 7.0 is: +0.13
    
    findSurfaceCharge("4fix",10.0)
        Charged residues: ...
        The expected surface charge of 4fix at pH 10.0 is: -3.86
    

# See Also

  * [findSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [get_area](/index.php/Get_Area "Get Area")



Retrieved from "[https://pymolwiki.org/index.php?title=FindSurfaceCharge&oldid=13951](https://pymolwiki.org/index.php?title=FindSurfaceCharge&oldid=13951)"


---

## FindSurfaceResidues

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/findSurfaceResidues.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findSurfaceResidues.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The findSurfaceResidues script will select (and color if requested) surface residues and atoms on an object or selection. See the options below. 

Each time, the script will create two new selections called, **exposed_res_XYZ** and **exposed_atm_XYZ** where _XYZ_ is some random number. This is done so that no other selections/objects are overwritten. 

# Usage
    
    
    findSurfaceResidues [ selection=all [, cutoff=2.5 [, doShow=0 ]]]
    

# Arguments

  * **selection** = str: The object or selection for which to find exposed residues {default: all}
  * **cutoff** = float: The cutoff in square Angstroms that defines exposed or not. Those atoms with > cutoff Ang^2 exposed will be considered _exposed_ {default: 2.5 Ang^2}
  * **doShow** = 0/1: Change the visualization to highlight the exposed residues vs interior {default: 0}



# Examples

  * [![Result of $TUT/1hpv.pdb at 2.5 Ang cutoff.](/images/f/fe/FindExRes.png)](/index.php/File:FindExRes.png "Result of $TUT/1hpv.pdb at 2.5 Ang cutoff.")

Result of $TUT/1hpv.pdb at 2.5 Ang cutoff. 

  * [![Example coloring of surface residues](/images/8/83/Surface_residues_ex.png)](/index.php/File:Surface_residues_ex.png "Example coloring of surface residues")

Example coloring of surface residues 



    
    
    run findSurfaceResidues.py
    fetch 1hpv, async=0
    
    findSurfaceResidues
    
    # now show the exposed
    findSurfaceResidues doShow=1
    
    # watch how the visualization changes:
    findSurfaceResidues doShow=1, cutoff=0.5
    findSurfaceResidues doShow=1, cutoff=1.0
    findSurfaceResidues doShow=1, cutoff=1.5
    findSurfaceResidues doShow=1, cutoff=2.0
    findSurfaceResidues doShow=1, cutoff=2.5
    findSurfaceResidues doShow=1, cutoff=3.0
    

# See Also

  * [get_sasa_relative](/index.php/Get_sasa_relative "Get sasa relative")
  * [get_area](/index.php/Get_Area "Get Area")



Retrieved from "[https://pymolwiki.org/index.php?title=FindSurfaceResidues&oldid=13953](https://pymolwiki.org/index.php?title=FindSurfaceResidues&oldid=13953)"


---

## Helicity check

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**helicity_check** shows the evolution of O—N distances over an amino acid sequence 

See for further info: 

["Models for the 3(10)-helix/coil, pi-helix/coil, and alpha-helix/3(10)-helix/coil transitions in isolated peptides"](http://www.proteinscience.org/cgi/content/abstract/5/8/1687). _Protein Sci_ ROHL and DOIG 5 (8) 1687. 

Uses: 

  * in the pymol console:


    
    
     >run pymol_helicity_check.py
       ----> select some consecutive amino acids
              - this is nicely done with the Display->Sequence tool
     >helicity_check()
    

  * installing helicity_check


    
    
     copy pymol_helicity_check.py in $PYMOL_INSTALL_DIR/modules/pmg_tk/startup
     launch Pymol: you now have a new option in the Plugin menu
    
    

helicity_check uses [gnuplot](http://www.gnuplot.info) to display its results. As a consequence gnuplot needs to be installed. 

This plugin was tested on linux only, it might need some modifications to run on other OSes. (hints: launching gnuplot and path to dumpfile) 

  

    
    
    # pymol_helicity_check.py
    # Copyright (c) 2006-2007 Julien Lefeuvre <lefeuvrejulien@yahoo.fr>
    #
    
    """
    Pymol plugin for checking helicity type
    
    helicity_check() takes as input a selection ('sele' by default)
    of at least 5 amino acids and computes the distances between
    O(i) - N(i+3)
    O(i) - N(i+4)
    O(i) - N(i+5)
    See for further info:
    Protein Sci ROHL and DOIG 5 (8) 1687
    'Models for the 3(10)-helix/coil, pi-helix/coil,
    and alpha-helix/3(10)-helix/coil transitions in isolated peptides.'
    
    uses:
    *in the pymol console:
      >run pymol_helicity_check.py
        ----> select some consecutive amino acids
               - this is nicely done with the Display->Sequence tool
      >helicity_check()
    *installing helicity_check
      copy pymol_helicity_check.py in $PYMOL_INSTALL_DIR/modules/pmg_tk/startup
      launch Pymol: you now have a new option in the Plugin menu
    
    helicity_check uses gnuplot (http://www.gnuplot.info) to display its results
    As a consequence gnuplot needs to be installed.
    
    This plugin was tested on linux only, it my need some modifications to run on
    other OSes (hints: launching gnuplot and path to dumpfile)
    """
    
    __author__ =    "Julien Lefeuvre <lefeuvrejulien@yahoo.fr>"
    __version__ =   "1.0"
    __date__ =      "2007-04-02"
    __copyright__ = "Copyright (c) 2007 %s. All rights reserved." % __author__
    __licence__ =   "BSD"
    
    from pymol import cmd
    from math import sqrt
    import sys
    import os
    import subprocess
    import time
    
    def __init__(self):
        """init function in order to have a nice menu option in Pymol"""
        self.menuBar.addmenuitem('Plugin', 'command', 'Helicity Check',
                 label='Helicity Check', command = lambda: helicity_check())
    
    
    class Residue(object):
    
        def __init__(self):
            self.name=None
            self.index=None
            self.Ocoord=None
            self.Ncoord=None
    
    
    def calc_distON(Ocoord,Ncoord):
        """return the distance between 2 atoms given their coordinates"""
        sum = 0
        for o, n in zip(Ocoord, Ncoord):
            sum += (o - n)**2
        return sqrt(sum)
    
    
    def helicity_check(selection='sele'):
        """calcultate distance O[res i]-N[res i+3]
                               O[res i]-N[res i+4]
                               O[res i]-N[res i+5]
        """
        seq_model = cmd.get_model(selection) #get info from selection
        res_lim = seq_model.get_residues()
    
        if len(res_lim)<5:
            sys.stderr.write("\nPlease select at least 5 residues\n")
            return
    
        atom_list = seq_model.atom
        res_data=[]
    
        for start,end in res_lim:   #extract the data we are interested in
            res=Residue()
            for atom in atom_list[start:end]:
                if atom.name == 'N':
                    res.name = atom.resn
                    res.index = atom.resi
                    res.Ncoord = atom.coord
                elif atom.name == 'O':
                    res.Ocoord = atom.coord
            if res.Ocoord and res.Ncoord and res.name and res.index:
                res_data.append(res)
            else:
                sys.stderr.write("\nPlease select complete protein residues\n")
                return
    
        res_list = [int(res.index) for res in res_data]
    
        if res_list != range(res_list[0], res_list[-1]+1):
            sys.stderr.write("\nPlease select a unbrocken residue sequence\n")
            return
    
        distON3 = []
        distON4 = []
        distON5 = []
        distONs = [distON3, distON4, distON5]
    
        for i,res in enumerate(res_data[:-5]): #distances calculations
            resis = res_data[i+3:i+6]
            for resi, distONi in zip(resis, distONs):
                distONi.append(calc_distON(res.Ocoord, resi.Ncoord))
    
        dump = os.tmpnam()+'.dat'
        dumpfile = file(dump, 'w')
    
        sys.stdout.write('\n#Distances O(i)---N(i+n)\n'
               '#ResNum , d(O(i)-N(i+3)) , d(O(i)-N(i+4)) , d(O(i)-N(i+4))\n')
        for i, d3, d4, d5 in zip(res_list, distON3, distON4, distON5):
            #writing console output
            sys.stdout.write(
                  '  %i ,      %f ,       %f ,       %f \n'%(i, d3, d4, d5))
            #writing data to a dump file for use by gnuplot
            dumpfile.write(
                  '  %i       %f        %f        %f \n'%(i, d3, d4, d5))
        dumpfile.flush()
    
        #launch a gnuplot window to show the distances
        gnuplotcmd = subprocess.Popen(['/usr/bin/gnuplot'], shell=True,
                                   stdin=subprocess.PIPE)
        gnuplotcmd.stdin.write('set autoscale\n')
        gnuplotcmd.stdin.write("plot "
             "'%s' using 1:2 title 'd(O(i)-N(i+3))' with lines, "
             "'%s' using 1:3 title 'd(O(i)-N(i+4))' with lines, "
             "'%s' using 1:4 title 'd(O(i)-N(i+5))' with lines\n'"
                              % (dump, dump, dump))
        time.sleep(3)
        dumpfile.close()
        os.remove(dump)
    

## Download

[Helicity_check-1.0.tar.zip](/images/f/fb/Helicity_check-1.0.tar.zip "Helicity check-1.0.tar.zip")

Retrieved from "[https://pymolwiki.org/index.php?title=Helicity_check&oldid=12200](https://pymolwiki.org/index.php?title=Helicity_check&oldid=12200)"


---

## HighlightAlignedSS

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![Highlight ss1.png](/images/2/22/Highlight_ss1.png)](/index.php/File:Highlight_ss1.png)

[](/index.php/File:Highlight_ss1.png "Enlarge")

This script will align and color the paired secondary structures of the two proteins a similar rainbow color. 

  

    
    
    from pymol import cmd, util
    
    def highlight_aligned_ss(obj1,obj2,transform=1,quiet=1):
        """
    DESCRIPTION
     
        Aligns two structures and colors their matching
        secondary structure elements with a matching
        rainbow colorscheme.
     
    USAGE
     
        highlight_aligned_ss obj1, obj2
    
        If transform=0 then the proteins are not
        moved after alignment.
    
    EXAMPLES
     
        highlight_aligned_ss 1cll, 1ggz
        highlight_aligned_ss 1rlw, 1byn and state 1
     
    SEE ALSO
     
        align
    
        JV 3-2-11
        """
    
    
        if not cmd.count_atoms(obj1):
            print "Error: Object 1 needs at least a few atoms to align."
            return None
        if not cmd.count_atoms(obj2):
            print "Error: Object 2 needs at least a few atoms to align."
            return None
    
        # align them
        uAln = cmd.get_unused_name("aln")
        cmd.align(obj1,obj2,object=uAln,transform=int(transform))
        cmd.hide("cgo", uAln)
    
        # select atoms of similar SS
        uSimSS = cmd.get_unused_name("similar_ss_")
        cmd.select(uSimSS, "((%s or %s) in %s) and (ss 'S' or ss 'H')" %
                   (obj1,obj2,uAln))
    
        # color by rainbow; these could be 
        # customized by function parameters
        util.rainbow(uSimSS + " and " + obj1)
        util.rainbow(uSimSS + " and " + obj2)
    
        # now color everything else grey
        cmd.color("grey70", "not (%s)" % uSimSS)
    
        # could also be an option to
        # update the representation
        # as cartoon
    
        # hide indicators
        cmd.select("none")
    
    cmd.extend("highlight_aligned_ss", highlight_aligned_ss)
    

Retrieved from "[https://pymolwiki.org/index.php?title=HighlightAlignedSS&oldid=11246](https://pymolwiki.org/index.php?title=HighlightAlignedSS&oldid=11246)"


---

## InterfaceResidues

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/0/09/1h4gi.png)](/index.php/File:1h4gi.png)

[](/index.php/File:1h4gi.png "Enlarge")

The red residues were detected as being in the interface of the two proteins. This is PDB 1H4G from the Example below.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
    * 3.1 Simple Example
    * 3.2 More Complex Example
  * 4 The Code



# Overview

This script finds interface residues between two proteins or chains, using the following concept. First, we take the area of the complex. Then, we split the complex into two pieces, one for each chain. Next, we calculate the chain-only surface area. Lastly, we take the difference between the comeplex-based areas and the chain-only-based areas. If that value is greater than your supplied cutoff, then we call it an interface residue. 

# Usage
    
    
    interfaceResidue complexName[, cA=firstChainName[, cB=secondChainName[, cutoff=dAsaCutoff[, selName=selectionNameToReturn ]]]]
    

where, 

**complexName**

    

    The name of the complex. cA and cB must be within in this complex

**cA**

    

    The name of the first chain to investigate

**cB**

    

    The name of the 2nd chain to investigate

**cutoff**

    

    The dASA cutoff, in sqaure Angstroms

**selName**

    

    Name of the selection to return.

For each residue in the newly defined interface, this script returns the model (cA or cB) name, the residue number and the change in area. To get that information use interfaceResidues in the api form as: 
    
    
    myInterfaceResidues = interfaceResidue(complexName[, cA=firstChainName[, cB=secondChainName[, cutoff=dAsaCutoff[, selName=selectionNameToReturn ]]]])
    

and the result, 
    
    
    myInterfaceResidues

will have all the value for you. 

# Examples

For these two examples, make sure you've run the script first. 

## Simple Example

This just finds the residues b/t chain A and chain B 
    
    
    fetch 1h4g, async=0
    interfaceResidues 1h4g
    

## More Complex Example
    
    
    fetch 1qox, async=0
    foundResidues = interfaceResidues("1qox", cA="c. I", cB="c. J", cutoff=0.75, selName="foundIn1QOX")
    

[![](/images/8/8d/1qoxi.png)](/index.php/File:1qoxi.png)

[](/index.php/File:1qoxi.png "Enlarge")

Interface residues between chains I and J in 1QOX.

# The Code
    
    
    from pymol import cmd, stored
    
    def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    	"""
    	interfaceResidues -- finds 'interface' residues between two chains in a complex.
    	
    	PARAMS
    		cmpx
    			The complex containing cA and cB
    		
    		cA
    			The first chain in which we search for residues at an interface
    			with cB
    		
    		cB
    			The second chain in which we search for residues at an interface
    			with cA
    		
    		cutoff
    			The difference in area OVER which residues are considered
    			interface residues.  Residues whose dASA from the complex to
    			a single chain is greater than this cutoff are kept.  Zero
    			keeps all residues.
    			
    		selName
    			The name of the selection to return.
    			
    	RETURNS
    		* A selection of interface residues is created and named
    			depending on what you passed into selName
    		* An array of values is returned where each value is:
    			( modelName, residueNumber, dASA )
    			
    	NOTES
    		If you have two chains that are not from the same PDB that you want
    		to complex together, use the create command like:
    			create myComplex, pdb1WithChainA or pdb2withChainX
    		then pass myComplex to this script like:
    			interfaceResidues myComlpex, c. A, c. X
    			
    		This script calculates the area of the complex as a whole.  Then,
    		it separates the two chains that you pass in through the arguments
    		cA and cB, alone.  Once it has this, it calculates the difference
    		and any residues ABOVE the cutoff are called interface residues.
    			
    	AUTHOR:
    		Jason Vertrees, 2009.		
    	"""
    	# Save user's settings, before setting dot_solvent
    	oldDS = cmd.get("dot_solvent")
    	cmd.set("dot_solvent", 1)
    	
    	# set some string names for temporary objects/selections
    	tempC, selName1 = "tempComplex", selName+"1"
    	chA, chB = "chA", "chB"
    	
    	# operate on a new object & turn off the original
    	cmd.create(tempC, cmpx)
    	cmd.disable(cmpx)
    	
    	# remove cruft and inrrelevant chains
    	cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
    	
    	# get the area of the complete complex
    	cmd.get_area(tempC, load_b=1)
    	# copy the areas from the loaded b to the q, field.
    	cmd.alter(tempC, 'q=b')
    	
    	# extract the two chains and calc. the new area
    	# note: the q fields are copied to the new objects
    	# chA and chB
    	cmd.extract(chA, tempC + " and (" + cA + ")")
    	cmd.extract(chB, tempC + " and (" + cB + ")")
    	cmd.get_area(chA, load_b=1)
    	cmd.get_area(chB, load_b=1)
    	
    	# update the chain-only objects w/the difference
    	cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
    	
    	# The calculations are done.  Now, all we need to
    	# do is to determine which residues are over the cutoff
    	# and save them.
    	stored.r, rVal, seen = [], [], []
    	cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')
    
    	cmd.enable(cmpx)
    	cmd.select(selName1, 'none')
    	for (model,resi,diff) in stored.r:
    		key=resi+"-"+model
    		if abs(diff)>=float(cutoff):
    			if key in seen: continue
    			else: seen.append(key)
    			rVal.append( (model,resi,diff) )
    			# expand the selection here; I chose to iterate over stored.r instead of
    			# creating one large selection b/c if there are too many residues PyMOL
    			# might crash on a very large selection.  This is pretty much guaranteed
    			# not to kill PyMOL; but, it might take a little longer to run.
    			cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))
    
    	# this is how you transfer a selection to another object.
    	cmd.select(selName, cmpx + " in " + selName1)
    	# clean up after ourselves
    	cmd.delete(selName1)
    	cmd.delete(chA)
    	cmd.delete(chB)
    	cmd.delete(tempC)
    	# show the selection
    	cmd.enable(selName)
    	
    	# reset users settings
    	cmd.set("dot_solvent", oldDS)
    	
    	return rVal
    
    cmd.extend("interfaceResidues", interfaceResidues)
    

Retrieved from "[https://pymolwiki.org/index.php?title=InterfaceResidues&oldid=12214](https://pymolwiki.org/index.php?title=InterfaceResidues&oldid=12214)"


---

## Intra xfit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.fitting](http://pymol.org/psicoredirect.php?psico.fitting)  
  
[![Intra xfit-1nmr.png](/images/2/2f/Intra_xfit-1nmr.png)](/index.php/File:Intra_xfit-1nmr.png)

[](/index.php/File:Intra_xfit-1nmr.png "Enlarge")

intra_xfit does a weighted superposition of a multi-state object, like an NMR-ensemble. The weights are estimated with maximum likelihood. 

This typically gives much better looking ensemble-fittings than the built-in [intra_fit](/index.php/Intra_fit "Intra fit") command, which does not do any outlier rejection or weighting. 

To superpose two different objects, use [xfit](/index.php/Xfit "Xfit"). 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Example
  * 5 See Also



## Installation

intra_xfit is available from the [psico](/index.php/Psico "Psico") package and requires [CSB](https://github.com/csb-toolbox/CSB). 

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c speleo3 csb
    

## Usage
    
    
    intra_xfit selection [, load_b [, cycles [, guide [, seed ]]]]
    

## Arguments

  * **selection** = string: atom selection from a single multi-state object
  * **load_b** = 0 or 1: save -log(weights) into B-factor column {default: 0}
  * **cycles** = int: number of weight refinement cycles {default: 20}
  * **guide** = 0 or 1: use only CA-atoms (protein) or C4' (nucleic acid) {default: 1}
  * **seed** = 0 or 1: use initial weights from current positions {default: 0}



## Example
    
    
    import psico.fitting
    fetch 1nmr, async=0
    
    intra_xfit 1nmr, load_b=1
    
    spectrum b, blue_white_red, guide
    

## See Also

  * [xfit](/index.php/Xfit "Xfit")
  * [intra_fit](/index.php/Intra_fit "Intra fit")



Retrieved from "[https://pymolwiki.org/index.php?title=Intra_xfit&oldid=12660](https://pymolwiki.org/index.php?title=Intra_xfit&oldid=12660)"


---

## Iterate sses

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A slightly more complex version of [ss](/index.php/Ss "Ss") that allows the user to pass in a function to act on the sse list. Of course, this requires the user to know about the internals of the sse class, but since this code is all open I doubt that matters.. 
    
    
    def iterate_sses(selection, action):
    
        class SSE(object):
    
            def __init__(self, start, typ, sseNumber):
                self.start, self.typ = start, typ
                self.end = -1
                self.sseNumber = sseNumber
    
            def __repr__(self):
                return "%s-%s %s" % (self.start, self.end, self.typ)
    
        stored.pairs = []
        cmd.iterate(selection, "stored.pairs.append((resi, ss))")
        num, currentType = stored.pairs[0]
    
        sseNumber = 1
        sses = [SSE(num, currentType, sseNumber)]
        currentSSE = sses[0]
        for resi, ssType in stored.pairs:
            if ssType == currentType:
                currentSSE.end = resi
            else:
                sseNumber += 1
                sses.append(SSE(resi, ssType, sseNumber))
                currentSSE = sses[-1]
                currentType = ssType
    
        for sse in sses:
            action(sse)
    
    cmd.extend("iterate_sses", iterate_sses)
    

As an example, here is a function that makes a series of selections, one for each sse, called "H1", "E2", and so on. Use it like: "iterate_sses('my_protein', doSelect)". 
    
    
    def doSelect(sse):
        cmd.select("%s%s" % (sse.typ, sse.sseNumber), "i. %s-%s" % (sse.start, sse.end))
    

Retrieved from "[https://pymolwiki.org/index.php?title=Iterate_sses&oldid=6366](https://pymolwiki.org/index.php?title=Iterate_sses&oldid=6366)"


---

## LigAlign

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 IMAGES
  * 2 DESCRIPTION
  * 3 Additional Information
  * 4 See Also



### IMAGES

  * [![Ligand-based alignment of 1XBB and 1OPJ](/images/7/7c/1xbb_proteins.png)](/index.php/File:1xbb_proteins.png "Ligand-based alignment of 1XBB and 1OPJ")

Ligand-based alignment of 1XBB and 1OPJ 

  * [![Atom-to-atom mapping of ligands](/images/b/b4/Mapping_example1.png)](/index.php/File:Mapping_example1.png "Atom-to-atom mapping of ligands")

Atom-to-atom mapping of ligands 

  * [![Alignment of 1V07 and 1HBI using heme](/images/c/c1/Heme18.png)](/index.php/File:Heme18.png "Alignment of 1V07 and 1HBI using heme")

Alignment of 1V07 and 1HBI using heme 




### DESCRIPTION

LigAlign is a tool to compare protein active-sites and investigate ligand binding. The active-site alignment is guided by the orientation of bound ligands in the protein active sites. LigAlign supports analysis of flexible ligand via automatic fragment-based alignment: first computing a natural fragmentation of the query ligand, aligning each fragment of the query independently against the baseline, and then permitting easy visualization of each active site subcavity.   
  


We use protein-ligand complexes to compare the active sites of several proteins which interact with a chosen ligand. Beginning with a user-specified protein-ligand structure, LigAlign gathers experimental structures of other proteins bound to the ligand from the Protein Data Bank. The tool then aligns the ligands bound in the structures to minimize the ligand-to-ligand RMSD. This transformation also aligns the active sites. Finally, the user can examine the aligned active sites to identify structural patterns, such as conserved steric hindrance or hydrophobicity.  
  


However, a flexible ligand can bend itself into different active sites, where the active site subcavities have different relative positions or orientations. Therefore, when comparing two active sites, a rigid RMSD-minimizing transform on docked flexible ligands may fail to align the correct portions of the active site. However, each subcavity should still exhibit chemical or geometric complementarity to the piece of the ligand which it binds. Please see [the website](http://compbio.cs.toronto.edu/ligalign) for more information.  
  


LigAlign simplifies a number of protein analysis tasks. For example, LigAlign will align similar but distinct ligands which, in the context of structure-based drug discovery, permits the comparison of the docking of different ligands. Alternatively, if the user only specifies one protein-ligand complex, LigAlign will find chemically similar ligands automatically via [the Protein-Small Molecule Database](http://compbio.cs.toronto.edu/psmdb). Finally, LigAlign improves workflow by automatically fetching necessary data from the [Protein Data Bank](http://www.rcsb.org).  
  


LigAlign is contributed by [Abraham Heifets](http://www.cs.toronto.edu/~aheifets) and [Ryan Lilien](http://www.cs.toronto.edu/~lilien) at the University of Toronto.  


### Additional Information

  * **[Downloads](http://compbio.cs.toronto.edu/ligalign/downloads.html)**
  * **[LigAlign Website](http://compbio.cs.toronto.edu/ligalign)**
  * **[Tutorial](http://compbio.cs.toronto.edu/ligalign/README.html)**



## See Also

  * [mcsalign](/index.php/Mcsalign "Mcsalign") (psico)



Retrieved from "[https://pymolwiki.org/index.php?title=LigAlign&oldid=12656](https://pymolwiki.org/index.php?title=LigAlign&oldid=12656)"


---

## Load aln

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.importing](http://pymol.org/psicoredirect.php?psico.importing)  
  
load_aln loads a pairwise sequence alignment file as an alignment object into PyMOL. 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Example
  * 5 See Also



## Installation

load_aln is available from the [psico](/index.php/Psico "Psico") package and requires [biopython](http://biopython.org/). 

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install biopython
    

## Usage
    
    
    load_aln filename [, object [, mobile [, target
        [, mobile_id [, target_id [, format [, transform ]]]]]]]
    

## Arguments

  * **filename** = str: alignment file
  * **object** = str: name of the object {default: filename prefix}
  * **mobile, target** = str: atom selections {default: ids from alignment file}
  * **mobile_id, target_id** = str: ids from alignment file {default: first two}
  * **format** = str: file format, see <http://biopython.org/wiki/AlignIO> {default: guess from first line in file}
  * **transform** = 0/1: superpose mobile on target (using [fit](/index.php/Fit "Fit")) {default: 0}



## Example

Alignment file (alignment.faa): 
    
    
    >seq1
    ACDEFG----HIKLMN
    >seq2
    ACNEYGGGGGHVRLMN
    

PyMOL script: 
    
    
    # create objects
    fab ACDEFGHIKLMN, m1
    fab ACNEYGGGGGHVRLMN, m2
    
    # load alignment
    import psico.importing
    load_aln alignment.faa, mobile=m1, target=m2
    
    # show sequence viewer
    set seq_view
    

Use the alignment object to superpose with [xfit](/index.php/Xfit "Xfit"): 
    
    
    import psico.fitting
    xfit m1, m2, match=alignment, cycles=100
    rebuild
    

## See Also

  * [align#Alignment Objects](/index.php/Align#Alignment_Objects "Align")



Retrieved from "[https://pymolwiki.org/index.php?title=Load_aln&oldid=12741](https://pymolwiki.org/index.php?title=Load_aln&oldid=12741)"


---

## Load new B-factors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [figshare](http://dx.doi.org/10.6084/m9.figshare.1176991)  
Author(s)  | [Pietro Gatti-Lafranconi](/index.php/User:PietroGattiLafranconi "User:PietroGattiLafranconi")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Limitations
  * 6 Examples
  * 7 The Code



## Introduction

Quick and simple script to replace B-factor values in a PDB structure. 

New B-factors are provided in an external .txt file, one per line. 

By default, the script will also redraw the selected molecule as cartoon putty and colour by B-factor 

  


## Usage
    
    
    loadBfacts mol, [startaa, [source, [visual Y/N]]]
    

  


## Required Arguments

  * **mol** = any object selection (within one single object though)



  


## Optional Arguments

  * **startaa** = number of amino acid the first value in 'source' refers to (default=1)
  * **source** = name of the file containing new B-factor values (default=newBfactors.txt)
  * **visual** = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)



## Limitations

For its very nature, this script is not suitable for complex cases in which only some B-factors are to be replaced, or on complex selections. In such cases, [data2bfactor.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/) is a better choice. 

B-factors have to be entered one per line, with no labels or amino acid ID. They will replace values of a continuous stretch of amino acids of equal length (or shorter) starting from the position provided with 'startaa'. 

'newBfactors.txt' looks like 
    
    
    2
    0
    0
    3
    7
    ...
    

  


## Examples
    
    
    PyMOL>loadbfacts 1LVM
    

[![example #1](/images/8/86/LoadBfacts1.png)](/index.php/File:LoadBfacts1.png "example #1")

  

    
    
    PyMOL>loadbfacts 1LVM and chain a, startaa=135
    

[![example #2](/images/2/2b/LoadBfacts02.png)](/index.php/File:LoadBfacts02.png "example #2")

  


## The Code
    
    
    from pymol import cmd, stored, math
    	
    def loadBfacts (mol,startaa=1,source="newBfactors.txt", visual="Y"):
    	"""
    	Replaces B-factors with a list of values contained in a plain txt file
    	
    	usage: loadBfacts mol, [startaa, [source, [visual]]]
     
    	mol = any object selection (within one single object though)
    	startaa = number of first amino acid in 'new B-factors' file (default=1)
    	source = name of the file containing new B-factor values (default=newBfactors.txt)
    	visual = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)
     
    	example: loadBfacts 1LVM and chain A
    	"""
    	obj=cmd.get_object_list(mol)[0]
    	cmd.alter(mol,"b=-1.0")
    	inFile = open(source, 'r')
    	counter=int(startaa)
    	bfacts=[]
    	for line in inFile.readlines():	
    		bfact=float(line)
    		bfacts.append(bfact)
    		cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%bfact)
    		counter=counter+1
    	if visual=="Y":
    		cmd.show_as("cartoon",mol)
    		cmd.cartoon("putty", mol)
    		cmd.set("cartoon_putty_scale_min", min(bfacts),obj)
    		cmd.set("cartoon_putty_scale_max", max(bfacts),obj)
    		cmd.set("cartoon_putty_transform", 0,obj)
    		cmd.set("cartoon_putty_radius", 0.2,obj)
    		cmd.spectrum("b","rainbow", "%s and n. CA " %mol)
    		cmd.ramp_new("count", obj, [min(bfacts), max(bfacts)], "rainbow")
    		cmd.recolor()
    
    cmd.extend("loadBfacts", loadBfacts);
    

Retrieved from "[https://pymolwiki.org/index.php?title=Load_new_B-factors&oldid=11711](https://pymolwiki.org/index.php?title=Load_new_B-factors&oldid=11711)"


---

## Mcsalign

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.mcsalign](http://pymol.org/psicoredirect.php?psico.mcsalign)  
  
mcsalign aligns two small-molecule selections based on Maximum-Common-Substructure. 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Example
  * 5 See Also



## Installation

mcsalign is available from the [psico](/index.php/Psico "Psico") package and requires [rdkit](http://www.rdkit.org/) and [csb](https://github.com/csb-toolbox/CSB). 

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c rdkit rdkit
    conda install -c speleo3 csb
    

## Usage
    
    
    mcsalign mobile, target [, mobile_state [, target_state
        [, cycles [, timeout [, method ]]]]]
    

## Arguments

  * **mobile** = str: atom selection of mobile object
  * **target** = str: atom selection of target object
  * **mobile_state** = int: object state of mobile selection {default: -1 = current state}
  * **target_state** = int: object state of target selection {default: -1 = current state}
  * **cycles** = int: number of weight-refinement iterations for weighted RMS fitting {default: 5}
  * **timeout** = int: MCS search timeout in seconds {default: 10}
  * **method** = indigo or rdkit {default: check availability}



## Example

Align Cytochrome C and Hemoglobin based on their Heme moieties: 
    
    
    fetch 3zcf 4n8t, async=0
    zoom /4n8t//A/HEM, animate=2, buffer=3
    
    import psico.mcsalign
    mcsalign /3zcf//A/HEC, /4n8t//A/HEM
    

Align a set of small molecules to a reference 
    
    
    fetch ala ile arg trp met
    extra_fit *, ala, method=mcsalign
    

## See Also

  * [align](/index.php/Align "Align")
  * [pair_fit](/index.php/Pair_fit "Pair fit")
  * [LigAlign](/index.php/LigAlign "LigAlign")



Retrieved from "[https://pymolwiki.org/index.php?title=Mcsalign&oldid=13302](https://pymolwiki.org/index.php?title=Mcsalign&oldid=13302)"


---

## Measure Distance

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    # This script writes the distance from 
    # atom mol1///25/ha to atom mol1///26/ha
    # out to the file "dist.txt"
    # Simply change your selections to see different distances.
    
    # import PyMOL's command namespace
    from pymol import cmd
    
    # open dist.txt for writing
    f=open('dist.txt','w')
    
    # calculate the distance and store it in dst
    dst=cmd.distance('tmp','mol1///25/ha','mol1///26/ha')
    
    # write the formatted value of the distance (dst)
    # to the output file
    f.write("%8.3f\n"%dst)
    
    # close the output file.
    f.close()
    

## See Also

[Distance](/index.php/Distance "Distance")

Retrieved from "[https://pymolwiki.org/index.php?title=Measure_Distance&oldid=6368](https://pymolwiki.org/index.php?title=Measure_Distance&oldid=6368)"


---

## Motif

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Small script to show backbone motifs. Use like "motif 10-12" to show a segment from residues 10 to 12. 

The (optional) chain argument is the chain letter, but the output lists the chain along with the residue, so it's not essential to use this. 

When run, the output looks like this (for pdbid 1arb): 
    
    
    PyMOL>motif 10-12
     VAL-10 : (  -62,  -31 ) AR
     VAL-11 : (  -72,   -8 ) AR
     CYS-12 : (  -70,  170 ) BR
    

The "AR", "BR", "AL", "BL" are conformation symbols, and cover quite large regions of the Ramachandran plot. to change them, change the typemap. 
    
    
    # very rough bounds!
    typemap = { (-180,0, 90, 180) : "BR", (-150, -30, -60, 60) : "AR", (0, 180, 90, 180) : "BL", (30, 150, -60, 60) : "AL" }
    
    def determinetype(phipsi):
        phi, psi = phipsi
        for bound in typemap:
            if bound[0] < phi < bound[1] and bound[2] < psi < bound[3]:
                return typemap[bound]
        return "?"
    
    def my_phi_psi(selection):
        r = cmd.get_phipsi(selection)
    
        if r is not None:
            keys = r.keys()
            keys.sort()
    
            cmd.feedback('push')
            cmd.feedback('disable','executive','actions')
            for key in keys:
                phipsiType = determinetype(r[key])
                argtuple = r[key] + (phipsiType,)
                cmd.iterate("(%s`%d)" % key, "print ' %-5s " + ("( %4.0f, %4.0f ) %s" % argtuple) + "'%(resn+'-'+resi+' '+chain+':')")
            cmd.feedback('pop')
    
    def motif(residueRange, chain=None):
        """
        Use like "motif 10-12" to show a backbone segment from residues 10 to 12.
        """
        name = residueRange
        if chain is None:
            selection = "name ca+n+c+o+h and not hetatm and resi %s" % residueRange
        else:
            selection = "name ca+n+c+o+h and not hetatm and resi %s and chain %s" % (residueRange, chain)
        cmd.select(name, selection)
        cmd.show("sticks", name)
        cmd.zoom(name)
        cmd.disable(name)
        cmd.label("name ca and %s" % selection, "resi")
        my_phi_psi(name)
    
    cmd.extend("motif", motif)
    cmd.extend("my_phi_psi", my_phi_psi)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Motif&oldid=6370](https://pymolwiki.org/index.php?title=Motif&oldid=6370)"


---

## Msms surface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.msms](http://pymol.org/psicoredirect.php?psico.msms)  
  
msms_surface calculates a molecular surface with [MSMS](http://mgltools.scripps.edu/packages/MSMS) and loads it into PyMOL as a CGO. 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 See Also



## Installation

msms_surface is available from the [psico](/index.php/Psico "Psico") package and requires the [msms](http://mgltools.scripps.edu/downloads#msms) binary. 

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c bioconda -c speleo3 msms
    

## Usage
    
    
    msms_surface [ selection [, state [, density [, name
        [, atomcolors [, exe [, preserve ]]]]]]]
    

## Arguments

  * **selection** = str: atom selection {default: polymer}
  * **state** = int: object state {default: 1}
  * **density** = float: MSMS surface point density {default: 3}
  * **name** = str: name of CGO object to create
  * **atomcolors** = 0/1: color surface by atom colors {default: 0}
  * **exe** = str: path to msms binary {default: find in $PATH}
  * **preserve** = 0/1: don't delete msms generated files {default: 0}



## Examples

Use atom colors: 
    
    
    import psico.msms
    fetch 1ubq, async=0
    msms_surface atomcolors=1
    

Use a solid color: 
    
    
    msms_surface name=bluesurf
    color blue, bluesurf
    

## See Also

  * [MSMS](/index.php/MSMS "MSMS") (plugin with GUI)
  * [surface](/index.php/Surface "Surface")



Retrieved from "[https://pymolwiki.org/index.php?title=Msms_surface&oldid=12650](https://pymolwiki.org/index.php?title=Msms_surface&oldid=12650)"


---

## Pairwise distances

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [figshare](http://dx.doi.org/10.6084/m9.figshare.1031600)  
Author(s)  | [Pietro Gatti-Lafranconi](/index.php/User:PietroGattiLafranconi "User:PietroGattiLafranconi")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Examples
  * 6 The Code



## Introduction

Given any two selections, this script calculates and returns the pairwise distances between all atoms that fall within a defined distance. 

Can be used to measure distances within the same chain, between different chains or different objects. 

Distances can be restricted to sidechain atoms only and the outputs either displayed on screen or printed on file. 

  


## Usage
    
    
    pairwise_dist sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]
    

  


## Required Arguments

  * **sel1** = first selection
  * **sel2** = second selection
  * **max_dist** = max distance in Angstroms



  


## Optional Arguments

  * **output** = accepts Screen/Print/None (default N)
  * **sidechain** = limits (Y) results to sidechain atoms (default N)
  * **show** = shows (Y) individual distances in pymol menu (default=N)



  


## Examples

**example #1**
    
    
    PyMOL>pairwise_dist 1efa and chain D, 1efa and chain B, 3, output=S, show=Y
     
    1efa/D/DC/13/OP1 to 1efa/B/TYR/47/OH: 2.765
    1efa/D/DC/13/OP2 to 1efa/B/LEU/6/N: 2.983
    1efa/D/DC/13/OP2 to 1efa/B/LEU/6/CB: 2.928
    1efa/D/DT/14/O4' to 1efa/B/ALA/57/CB: 2.827
    1efa/D/DT/14/OP1 to 1efa/B/ASN/25/OD1: 2.858
    1efa/D/DT/14/OP1 to 1efa/B/GLN/54/NE2: 2.996
    1efa/D/DT/14/OP2 to 1efa/B/SER/21/OG: 2.517
    1efa/D/DC/15/N4 to 1efa/B/GLN/18/NE2: 2.723
    1efa/D/DA/16/N6 to 1efa/B/GLN/18/NE2: 2.931
     
    Number of distances calculated: 9
    

[![example #1](/images/d/df/Pairwise1.png)](/index.php/File:Pairwise1.png "example #1")

  
**example #2**
    
    
    PyMOL>pairwise_dist 2w8s and chain a, 2W8s and chain b, 3, sidechain=Y, output=S, show=Y
     
    2W8S/A/GLN/432/OE1 to 2W8S/B/ARG/503/NH2: 2.758
    2W8S/A/ASP/434/OD1 to 2W8S/B/SER/493/OG: 2.444
    2W8S/A/TYR/447/OH to 2W8S/B/GLN/485/NE2: 2.878
    2W8S/A/HIS/449/NE2 to 2W8S/B/SER/489/OG: 2.686
    2W8S/A/GLN/485/OE1 to 2W8S/B/TYR/447/OH: 2.971
    2W8S/A/SER/489/OG to 2W8S/B/HIS/449/NE2: 2.913
    2W8S/A/SER/493/OG to 2W8S/B/ASP/434/OD1: 2.491
    2W8S/A/ARG/503/NH2 to 2W8S/B/GLN/432/OE1: 2.653
     
    Number of distances calculated: 8
    

[![example #2](/images/3/38/Pairwise2.png)](/index.php/File:Pairwise2.png "example #2")

  
**example #3**
    
    
    PyMOL>pairwise_dist 1FOS and chain H and resi 290:300, 1FOS and chain G, 4, sidechain=Y, show=Y, output=P
     
    Results saved in IntAtoms_4.txt
    Number of distances calculated: 24
    

[![example #3](/images/d/d2/Pairwise3.png)](/index.php/File:Pairwise3.png "example #3")

  


## The Code

Copy the following text and save it as pairwisedistances.py 
    
    
    from __future__ import print_function
    from pymol import cmd, stored, math
    
    def pairwise_dist(sel1, sel2, max_dist, output="N", sidechain="N", show="N"):
    	"""
    	usage: pairwise_dist sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]
    	sel1 and sel2 can be any to pre-existing or newly defined selections
    	max_dist: maximum distance in Angstrom between atoms in the two selections
    	--optional settings:
    	output: accepts Screen/Print/None (default N)
    	sidechain: limits (Y) results to sidechain atoms (default N)
    	show: shows (Y) individual distances in pymol menu (default=N)
    	"""
    	print("")
    	cmd.delete("dist*")
    	extra=""
    	if sidechain=="Y":
    		extra=" and not name c+o+n"
    	
    	#builds models
    	m1 = cmd.get_model(sel2+" around "+str(max_dist)+" and "+sel1+extra)
    	m1o = cmd.get_object_list(sel1)
    	m2 = cmd.get_model(sel1+" around "+str(max_dist)+" and "+sel2+extra)
    	m2o = cmd.get_object_list(sel2)
    
    	#defines selections
    	cmd.select("__tsel1a", sel1+" around "+str(max_dist)+" and "+sel2+extra)
    	cmd.select("__tsel1", "__tsel1a and "+sel2+extra)
    	cmd.select("__tsel2a", sel2+" around "+str(max_dist)+" and "+sel1+extra)
    	cmd.select("__tsel2", "__tsel2a and "+sel1+extra)
    	cmd.select("IntAtoms_"+max_dist, "__tsel1 or __tsel2")
    	cmd.select("IntRes_"+max_dist, "byres IntAtoms_"+max_dist)
     
    	#controlers-1
    	if len(m1o)==0: 
    		print("warning, '"+sel1+extra+"' does not contain any atoms.")
    		return
    	if len(m2o)==0: 
    		print("warning, '"+sel2+extra+"' does not contain any atoms.")
    		return
    	
    	#measures distances
    	s=""
    	counter=0
    	for c1 in range(len(m1.atom)):
    		for c2 in range(len(m2.atom)):
    			distance=math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord,m2.atom[c2].coord))))
    			if distance<float(max_dist):
    				s+="%s/%s/%s/%s/%s to %s/%s/%s/%s/%s: %.3f\n" % (m1o[0],m1.atom[c1].chain,m1.atom[c1].resn,m1.atom[c1].resi,m1.atom[c1].name,m2o[0],m2.atom[c2].chain,m2.atom[c2].resn,m2.atom[c2].resi,m2.atom[c2].name, distance)
    				counter+=1
    				if show=="Y": cmd.distance (m1o[0]+" and "+m1.atom[c1].chain+"/"+m1.atom[c1].resi+"/"+m1.atom[c1].name, m2o[0]+" and "+m2.atom[c2].chain+"/"+m2.atom[c2].resi+"/"+m2.atom[c2].name)
    
    	#controler-2
    	if counter==0: 
    		print("warning, no distances were measured! Check your selections/max_dist value")
    		return
    	
    	#outputs
    	if output=="S":
    		print(s)
    	if output=="P":
    		f=open('IntAtoms_'+max_dist+'.txt','w')
    		f.write("Number of distances calculated: %s\n" % (counter))
    		f.write(s)
    		f.close()
    		print("Results saved in IntAtoms_%s.txt" % max_dist)
    	print("Number of distances calculated: %s" % (counter))
    	cmd.hide("lines", "IntRes_*")
    	if show=="Y": cmd.show("lines","IntRes_"+max_dist)
    	cmd.deselect()
      
    cmd.extend("pairwise_dist", pairwise_dist)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pairwise_distances&oldid=12898](https://pymolwiki.org/index.php?title=Pairwise_distances&oldid=12898)"


---

## PDB plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 PDB PyMOL plugin
    * 1.1 Installation
      * 1.1.1 PyMOL 1.6
      * 1.1.2 PyMOL 1.7 or later
    * 1.2 PDB Analysis - All
    * 1.3 PDB Analysis - Molecules
    * 1.4 PDB Analysis - Domains
    * 1.5 PDB Analysis - Validation
    * 1.6 PDB Analysis - Assemblies
    * 1.7 View assemblies for your own mmCIF file
    * 1.8 Acknowledgements



## PDB PyMOL plugin

PDBe’s PyMOL plugin provides an easy way to visualise PDB data and annotations in PyMOL. This includes displaying individual molecules, Pfam, SCOP and CATH domains and viewing geometric validation for PDB entries based on PDBe's API. 

### Installation

The PDB PyMOL plugin works on PyMOL 1.6 or later 

#### PyMOL 1.6

Download the plugin from the following URL: 

<http://www.ebi.ac.uk/pdbe/pdb-component-library/pymol_plugin/PDB_plugin.py>

#### PyMOL 1.7 or later

Open the “plugin manager” Enter the following url in the “Install from PyMOL Wiki or URL 

<http://www.ebi.ac.uk/pdbe/pdb-component-library/pymol_plugin/PDB_plugin.py>

This will add the following items to the Plugins menu: 

PDB Analysis - All 

PDB Analysis - Molecules 

PDB Analysis - Domains 

PDB Analysis - Validation 

PDB Analysis - Assemblies 

  
Screenshot of the menu 

  
[![Screenshot of the PyMOL plugin menu with the PDB plugin installed](/images/e/e1/PDB_plugin.png)](/index.php/File:PDB_plugin.png "Screenshot of the PyMOL plugin menu with the PDB plugin installed")

  
To use each component click the menu item in the Plugins menu and then enter and PDB ID into the dialogue box. 

### PDB Analysis - All

This displays chemically distinct molecules, domains which make up a PDB entry. It also displays the assemblies within a PDB entry. 

Chemically distinct molecules which make up a PDB entry, including proteins, nucleic acids and ligands, are highlighted. Each molecule is given a unique colour and selected as a different object. Polymers (protein, DNA and RNA chains) are shown as cartoon or ribbon representation and ligands are shown as spheres. Pfam, SCOP and CATH domains within a PDB entry are highlighted. Each domain is coloured in a different colour and selected as a different object. This enables each domain to be turned on and off with PyMOL. Assemblies annotated within a PDB entry are highlighted. This requires pymol 1.76 or later. 

All PDB annotations show in PDB entry 3unb 

  
[![All PDB annotations shown on PDB entry 3unb](/images/2/26/Pdbe_3unb_all_annotation.png)](/index.php/File:Pdbe_3unb_all_annotation.png "All PDB annotations shown on PDB entry 3unb")

  


### PDB Analysis - Molecules

This highlights the chemically distinct molecules which make up a PDB entry, including proteins, nucleic acids and ligands. Each molecule is given a unique colour and selected as a different object. Polymers (protein, DNA and RNA chains) are shown as cartoon or ribbon representation and ligands are shown as spheres. 

In the below example the distinct molecules in PDB entry 3l2p are shown. 

  
[![Distinct molecules in PDB entry 3l2p](/images/b/b8/Pdbe_3l2p_entity.png)](/index.php/File:Pdbe_3l2p_entity.png "Distinct molecules in PDB entry 3l2p")

  


### PDB Analysis - Domains

This highlights the Pfam, Rfam, SCOP and CATH domains within a PDB entry. Each domain is coloured in a different colour and selected as a different object. This enables each domain to be turned on and off with PyMOL. The Domains Plugin also highlights chemically distinct molecules and domains are overlaid on these molecules. 

The CATH domains in PDB entry 3b43 are highlighted in the figure below. 

  
[![CATH domains shown on PDB entry 3b43](/images/0/0f/PDBe_3b43_domains_2.png)](/index.php/File:PDBe_3b43_domains_2.png "CATH domains shown on PDB entry 3b43")

### PDB Analysis - Validation

This overlays geometric validation on a PDB entry. Geometry outliers are coloured as per the PDB validation report: Residues with no geometric outliers are shown in green Residues with one geometric outliers are shown in yellow Residues with two geometric outliers are shown in orange Residues with three or more geometric outliers are shown in red. 

Geometric validation is shown on PDB entry 2gc2 

  
[![Geometric validation shown on PDB entry 2gc2](/images/7/78/Pdbe_2gc2_validation.png)](/index.php/File:Pdbe_2gc2_validation.png "Geometric validation shown on PDB entry 2gc2")

### PDB Analysis - Assemblies

This highlights the assemblies within a PDB entry. This requires pymol 1.76 or later. 

Assemblies are shown within PDB entry 5j96, with each assembly as an object in pymol and coloured in a distinct colour. 

  
[![Assemblies shown for PDB entry 5j96](/images/8/85/PDB_5j96_assemblies.png)](/index.php/File:PDB_5j96_assemblies.png "Assemblies shown for PDB entry 5j96")

  


### View assemblies for your own mmCIF file

To view annotated assemblies in your own mmCIF file - for example after PDB annotation - run pymol with the following command: 

pymol -r PDB_plugin.py -- mmCIF_file=YOUR_MMCIF_FILE 

where YOUR_MMCIF_FILE is the filename you want display the assemblies for. 

  


### Acknowledgements

This plugin was developed by PDBe [![PDBe logo](/images/7/75/Pdbe_logo.png)](http://pdbe.org "PDBe logo") who are a founding member of the wwPDB [![wwPDB logo](/images/d/de/Wwpdb-logo.png)](http://www.wwpdb.org/ "wwPDB logo")

Retrieved from "[https://pymolwiki.org/index.php?title=PDB_plugin&oldid=12772](https://pymolwiki.org/index.php?title=PDB_plugin&oldid=12772)"


---

## Propka

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/propka.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/propka.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Acknowledgement
  * 2 Introduction
  * 3 Dependency of python module: mechanize
    * 3.1 Examples
  * 4 Use of script
    * 4.1 Input paramaters
  * 5 Examples
    * 5.1 Example 1 - Mutagenesis analysis - part 1
    * 5.2 Example 2 - Mutagenesis analysis - part 2
    * 5.3 Example 3 - Scan a range of proteins
  * 6 ScriptVersion
  * 7 Changelog
    * 7.1 Known bugs



## Acknowledgement

propka.py contact and relies on the result from the [propka](http://propka.ki.ku.dk) server 

The PROPKA method is developed by the [Jensen Research Group](http://propka.ki.ku.dk/~luca/wiki/index.php/JansPage) , Department of Chemistry, University of Copenhagen.   


## Introduction

This script can fetch the pka values for a protein from the [propka](http://propka.ki.ku.dk/) server. The "propka" function downloads the results and processes them.   
It also automatically writes a pymol command file and let pymol execute it. This command file make pka atoms, rename them, label them and color them according to the pka value. 

If you put the mechanize folder and the propka.py script somewhere in your pymol search path, then getting the pka values is made super easy. By the way, did you know, that you don't have to prepare the .pdb file by adding/removing hydrogens? The [propka](http://propka.ki.ku.dk/) server uses its own internal hydrogen placement algorithm. 
    
    
    import propka
    fetch 4ins, async=0
    propka
    

If there is no web connection, it is possible to process a result file from a previous run or from a downloaded propka webpage result. This can be a handsome feature in a teaching/seminar situation, since it speeds up the pymol result or that an available web connection can be doubtful. Just point to the .pka file: Remember the dot "." which means "current directory". 
    
    
    import propka
    load 4ins.pdb
    propka pkafile=./Results_propka/4ins"LOGTIME".pka, resi=18.25-30, resn=cys
    

The last possibility, is just to ask for the pka values of a recognized PDB id. This is done with the "getpropka" function. 
    
    
    import propka
    getpropka source=ID, PDBID=4ake, logtime=_, showresult=yes
    

## Dependency of python module: mechanize

The script needs mechanize to run. The module is included in the project, [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro). 

  * On windows, it is not easy to make additional modules available for pymol. So put in into your working folder.
  * The easy manual way:


  1. Go to: <http://wwwsearch.sourceforge.net/mechanize/download.html>
  2. Download mechanize-0.2.5.zip. <http://pypi.python.org/packages/source/m/mechanize/mechanize-0.2.5.zip>
  3. Extract to .\mechanize-0.2.5 then move the in-side folder "mechanize" to your folder with propka.py. The rest of .\mechanize-0.2.5 you don't need.


  * You can also see other places where you could put the "mechanize" folder. Write this in pymol to see the paths where pymol is searching for "mechanize"
  * import sys; print(sys.path)



### Examples

Read about the proteins here:  
<http://www.proteopedia.org/wiki/index.php/4ins>   
<http://www.proteopedia.org/wiki/index.php/1hp1>   

    
    
    import propka
    
    fetch 4ins, async=0
    propka        OR
    propka 4ins   OR
    propka 4ins, resi=19.20, resn=ASP.TYR, logtime=_, verbose=yes
    
    import propka
    fetch 1hp1, async=0
    propka molecule=1hp1, chain=A, resi=305-308.513, resn=CYS, logtime=_
    
    import propka
    getpropka source=ID, PDBID=4ins, logtime=_, server_wait=3.0, verbose=yes, showresult=yes
    

  * [![propka used on 1HP1.](/images/9/94/Propka1HP1.png)](/index.php/File:Propka1HP1.png "propka used on 1HP1.")

propka used on 1HP1. 

  * [![propka used on 1HP1, zoom on ATP ligand.](/images/2/29/Propka1HP1-ATP.png)](/index.php/File:Propka1HP1-ATP.png "propka used on 1HP1, zoom on ATP ligand.")

propka used on 1HP1, zoom on ATP ligand. 

  * [![propka used on 1HP1, zoom on ZN ligand metal center.](/images/a/af/Propka1HP1-ZN.png)](/index.php/File:Propka1HP1-ZN.png "propka used on 1HP1, zoom on ZN ligand metal center.")

propka used on 1HP1, zoom on ZN ligand metal center. 



  * [![propka used on 4INS.](/images/7/75/Propka4INS.png)](/index.php/File:Propka4INS.png "propka used on 4INS.")

propka used on 4INS. 

  * [![The easy menu does it easy to click on/off ligands and bonds](/images/f/f7/Pymolpropkamenu.png)](/index.php/File:Pymolpropkamenu.png "The easy menu does it easy to click on/off ligands and bonds")

The easy menu does it easy to click on/off ligands and bonds 

  * [![An appending logfile ./Results_propka/_Results.log saves the input commands and the result for the residues specified with resi= and resn=](/images/a/a6/Resultslog.png)](/index.php/File:Resultslog.png "An appending logfile ./Results_propka/_Results.log saves the input commands and the result for the residues specified with resi= and resn=")

An appending logfile ./Results_propka/_Results.log saves the input commands and the result for the residues specified with resi= and resn= 




pka atoms are created and renamed for their pka value. That makes it easy to "click" the atom in pymol and instantly see the pka value. 

The atoms b value are also altered to the pka value, and the atoms are then spectrum colored from pka=0-14. 

The pka value of 99.9 represent a di-sulphide bond, and is colored gold and the sphere size is set a little bigger. 

If one wants to see the specified result, the logfile ./Results_propka/_Results.log saves the link to the propka server. Here one can see in an interactive Jmol appp, the interactions to the pka residues. 

## Use of script
    
    
    cd /home/tlinnet/test
    
    import propka
    
    ### The fastest method is just to write propka. Then the last pymol molecule is assumed and send to server. verbose=yes makes the script gossip mode.
    fetch 4ins, async=0
    propka
    
    ### Larger protein
    fetch 1hp1, async=0
    propka logtime=_, resi=5-10.20-30, resn=CYS.ATP.TRP, verbose=yes
    
    ### Fetch 4ins from web. async make sure, we dont execute script before molecule is loaded. The resi and resn prints the interesting results right to command line.
    fetch 4ins, async=0
    propka chain=*, resi=5-10.20-30, resn=ASP.CYS, logtime=_
    
    ### If there is no web connection, one can process a local .pka file. Either from a previous run or from a downloaded propka webpage result.
    ### Then run and point to .pka file with: pkafile=./Results_propka/pkafile.pka Remember the dot "." in the start, to make it start in the current directory.
    load 4ins.pdb
    propka pkafile=./Results_propka/4ins_.pka, resi=18.25-30, resn=cys,
    
    ### Some more examples. This molecule has 550 residues, so takes a longer time. We select to run the last molecule, by writing: molecule=1hp1
    fetch 4ins, async=0
    fetch 1hp1, async=0
    propka molecule=1hp1, chain=A, resi=300-308.513, resn=CYS.ATP.TRP, logtime=_, verbose=no, showresult=no
    propka molecule=1hp1, pkafile=./Results_propka/1hp1_.pka, verbose=yes
    

### Input paramaters
    
    
    ############################################Input parameters: propka############################################
    ############# The order of input and changable things:
    propka(molecule="NIL",chain="*",resi="0",resn="NIL",method="upload",logtime=time.strftime("%m%d",time.localtime()),server_wait=3.0,version="v3.1",verbose="no",showresult="no",pkafile="NIL")
    # method : method=upload is default. This sends .pdb file and request result from propka server.
    ## method=file will only process a manual .pka file, and write a pymol command file. No use of mechanize.
    ## If one points to an local .pka file, then method is auto-changed to method=file. This is handsome in off-line environment, ex. teaching or seminar.
    # pkafile: Write the path to .pka file. Ex: pkafile=./Results_propka/4ins_.pka
    # molecule : name of the molecule. Ending of file is assumed to be .pdb
    # chain : which chains are saved to file, before molecule file is send to server. Separate with "." Ex: chain=A.b
    # resi : Select by residue number, which residues should be printed to screen and saved to the log file: /Results_propka/_Results.log.
    ## Separate with "." or make ranges with "-". Ex: resi=35.40-50
    # resn : Select by residue name, which residues should be printed to screen and saved to the log file: /Results_propka/_Results.log.
    ## Separate with "." Ex: resn=cys.tyr
    # logtime : Each execution give a set of files with the job id=logtime. If logtime is not provided, the current time is used.
    ## Normal it usefull to set it empty. Ex: logtime=_
    # verbose : Verbose is switch, to turn on messages for the mechanize section. This is handsome to see how mechanize works, and for error searching.
    # showresult : Switch, to turn on all results in pymol command window. Ex: showresult=yes
    # server_wait=10.0 is default. This defines how long time between asking the server for a result. Set no lower than 3 seconds.
    # version=v3.1 is default. This is what version of propka which would be used.
    ## Possible: 'v3.1','v3.0','v2.0'. If a newer version is available than the current v3.1, a error message is raised to make user update the script.
    ############################################Input parameters: getpropka############################################
    ############# The order of input and changable things:
    getpropka(PDB="NIL",chain="*",resi="0",resn="NIL",source="upload",PDBID="",logtime=time.strftime("%Y%m%d%H%M%S",time.localtime()),server_wait=3.0,version="v3.1",verbose="no",showresult="no")
    # PDB: points the path to a .pdb file. This is auto-set from propka function.
    # source : source=upload is default and is set at the propka webpage.
    # source=ID, PDBID=4ake , one can print to the command line, the pka value for any official pdb ID. No files are displayed in pymol.
    # PDBID: is used as the 4 number/letter pdb code, when invoking source=ID.
    

## Examples

### Example 1 - Mutagenesis analysis - part 1

This script was developed with the intention of making analysis of possible mutants easier. For example, the reactivity of Cysteines in FRET maleimide labelling is determined by the fraction of the Cysteine residue which is negatively charged (C-). This fraction is related to its pKa value and the pH of the buffer: f(C-)=1/(10(pK-pH)+1). So, one would be interested in having the lowest possible pKa value as possible. Ideally lower than the pH of the buffer. To analyse where to make the best mutant in your protein, you could do the following for several residues. We do the mutagenesis in the command line, since we then could loop over the residues in the protein.  
So, in a loop with defined residues, this could look like the following code. Note, now we are quite happy for the result log file, since it collects the pka for the mutants. 

Download: [examples/propka_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_1.pml>" highlight="python" />

### Example 2 - Mutagenesis analysis - part 2

To only loop over surface residues, you might want to find these with the script [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues"). 

Download: [examples/propka_2.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_2.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_2.pml>" highlight="python" />

A little warning though. You need to be carefull about the rotamer you're choosing. It can happen that the first rotamer ends up being in physically non-reasonable contact distance to other residues, so atoms become overlayed. Also, the mutagenesis wizard can have the funny habit of sometimes not adding hydrogens to terminal -C or -N after mutating a residue. 

### Example 3 - Scan a range of proteins

Could be done with this script 

Download: [examples/propka_3.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_3.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/propka_3.pml>" highlight="python" />

## ScriptVersion

Current_Version=20111202 

## Changelog

  * 20111202 
    1. This code has been put under version control. In the project, [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro).


  * 20110823 
    1. Fixed some issues with selection algebra of Ligands.
    2. Now colors Ligands automatically to purple scheme.


  * 20110822 
    1. Made the naming scheme consistent, so one can work with multiple proteins, and the grouping still works.
    2. Bonds to N-terminal and C-terminal did not show up. Fixed.
    3. If one just write "propka", the "last" molecule in the pymol object list is now assumed, instead of the first. This makes mutagenesis analysis easier.
    4. The pka difference from assumed standard values are now also displayed. Standard values are set to: pkadictio = {'ASP':3.9, 'GLU':4.3, 'ARG':12.0, 'LYS':10.5, 'HIS':6.0, 'CYS':8.3, 'TYR':10.1}
    5. The menu size is made bigger, so it can fit the long names for the bonding partners. There's also a slider that you can drag using the mouse 

    \--it's the tiny tab between the command line and the movie controls. Drag that left or right.


  * 20110821 
    1. The resn function was made wrong, only taking the last item of list. Fixed.
    2. resn= can now also be print to screen/log for ligands like ATP. Ex: resn=CYS.ATP.TRP
    3. Ligands are now also written to the "stripped-file". "./Results_propka/PDB.stripped"
    4. Bonds are now also made for Ligands.


  * 20110820 
    1. If one points to a result .pka file, then "method" is automatically set to "method=file".
    2. Added the ability for pka values for ligands
    3. Bonds are now generated for the pka atoms.
    4. The color scheme is changed from "rainbow" to "red_white_blue". This is easier to interpret.
    * CC: COULOMBIC INTERACTION. Color is red.
    * SH: SIDECHAIN HYDROGEN BOND. Color is brightorange.
    * BH: BACKBONE HYDROGEN BOND. Color is lightorange.


  * 20110817 
    1. If just invoking with "propka", it will select the first molecule. And now it is possible also to write. "propka all".
    2. Removed the "raise UserWarning" if the script is oudated. Only a warning message is printed.


  * 20110816 
    1. Made the execution of the pymol command script silent, by only using cmd.API. This will raise a Warning, which can be ignored.
    2. Built-in a Script version control, to inform the user if the propka script is updated on this page
    3. The alternate attribute for the labeling atoms are reset. It was found pymol objected altering names for atoms which had alternate positions/values.
    4. Reorganized the input order, which means that: molecule=all is default.



### Known bugs

  * Bonds can be multiplied from amino acids to Ligands like ZN or ATP. Assume the shortest bond to be correct. 
    1. ZN: This is caused, since ZN has no ID number and when there are several in the same chain. This can be reproduced for 1HP1, bond: A41ZNCC. Here it shows to bonds, where only the shortes is correct. No fix possible.
    2. ATP: This is caused, since propka sometimes eats the ending of the atom name. In the .pdb file is written O3', which propka represents like O3. Therefore a wildcard "*" is generel inserted, which can cause multiple bonds to the ATP molecule. This can be reproduced for 1HP1, bond: A504ATPSH. Here multiple bonds are made to all the O3*,O2* atoms. No fix possible.
  * The propka server sometimes use another naming scheme for Ligands. These Ligands will not be pka labelled or shown in pymol. The results will still downloaded to "/.Results_propka/file.pka" 
    1. This can be reproduced with 1BJ6. Here the propka webpage predicts the pka values for the ligands: AN7,GN1,CN3,CO2,GO2, but the .pdb file names these ligands: DA,DA,DC,DG.
  * Alternative configurations of a Ligand is at the moment a problem, and will not be shown. For example 1HXB and the ligand ROC. 

    The script extraxt AROC and BROC from pymol, which does match with ROC in the .pka file. Try save the protein with only one configuration of the Ligand.
    * >In pymol: 

    fetch 1hxb, async=0
    create 1hxbA, 1hxb and not alt B
    save 1hxbA.pdb, 1hxbA
    * >Quit pymol
    * >Manually replace with text editor in .pdb file "AROC" with " ROC" Remember the space " "!
    * >Then in pymol 

    load 1hxbA.pdb
    import propka
    propka verbose=yes



Retrieved from "[https://pymolwiki.org/index.php?title=Propka&oldid=13919](https://pymolwiki.org/index.php?title=Propka&oldid=13919)"


---

## Resicolor plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/resicolor_plugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/resicolor_plugin.py)  
Author(s)  | [Philaltist](/index.php/User:Philaltist "User:Philaltist")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Description

A small plugin to color proteins according to residue type. This functionality is also available as a script ([resicolor](/index.php/Resicolor "Resicolor")). 

## Example of use

If you follow the instructions for the [Pymol-script-repo](http://www.pymolwiki.org/index.php/Git_intro), it is available from the Plugin menu. 

Retrieved from "[https://pymolwiki.org/index.php?title=Resicolor_plugin&oldid=10440](https://pymolwiki.org/index.php?title=Resicolor_plugin&oldid=10440)"


---

## RmsdByResidue

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [Zhenting Gao](/index.php/User:Zhentg "User:Zhentg")  
License  |   
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Examples
  * 6 The Code



## Introduction

RMSD between two structures of the same protein. The concept is similar as RMSF between two structures. 

  * programming 
    * Platform 
      * PyMOL 
        * rms_cur
    * Feature 
      * Output 
        * RMSD of all atoms of each residues pairs
        * Least RMSD of all atoms of each residues pairs 
          * symmetry of Phe, Tyr, His, Asp, Glu, Gln, Asn, Arg, Leu and Val needs to be considered 
            * switch the atom name and then calculate the RMSD again
            * Selected least RMSD of a residue pair for report 
              * RMSD of backbone atoms of each residues pairs
              * RMSD of C alpha atoms of each residues pairs
          * With defined residues pairs 
            * Residue pair can be limited to within binding site
    * Workflow 
      * Read reference and target pdb files 
        * two structures should be superposed before using this function
  * Note 
    * Python
    * PyMOL 
      * Clean attributes 
        * otherwise rms_cur will fail
      * How to get residue name? 
        * residue name, residue index and etc. can only be read from an atom
  * Reference 
    * <https://sourceforge.net/p/pymol/mailman/message/32710889/>



## Usage

  * Open PyMOL
  * Load PDB files
  * run this Python script inside PyMOL
  * call the function 
    * rmsdByRes pdb1, resSelection, pdb2



## Required Arguments

Text 

## Optional Arguments

Text 

## Examples

Text 

## The Code
    
    
    #load library
    from pymol import cmd, stored
    import sys
    import os
    
    #function to judge if a file exists
    def is_non_zero_file(fpath):
        return os.path.isfile(fpath) and os.path.getsize(fpath) > 0
    
    #function to switch atom names within a residue
    def switchName(residueSelection, atomName1, atomName2):
      """
      switch the name of two atoms
      """
      cmd.alter(residueSelection+" and name "+atomName1, 'name="gzt"')
      cmd.alter(residueSelection+" and name "+atomName2, 'name="'+atomName1+'"')
      cmd.alter(residueSelection+" and name gzt", 'name="'+atomName2+'"')
    
    
    #function to change atom names of some residues with symetric sidechain
    def flipAtomName(targetResidueSelection):
      """
      switch the atom names of specific residues
      """
    # Create flipped residue
      cmd.create("flippedRes",targetResidueSelection+" and not alt B")
      targetResidueCa=cmd.get_model("flippedRes and name CA")
      for g in targetResidueCa.atom:
    #   print g.resn
        if g.resn=='ARG':
         switchName("flippedRes", "NH1", "NH2")
        elif g.resn=='HIS':
         switchName("flippedRes", "ND1", "CD2")
         switchName("flippedRes", "CE1", "NE2")
        elif g.resn=='ASP':
         switchName("flippedRes", "OD1", "OD2")
        elif g.resn=='PHE':
         switchName("flippedRes", "CD1", "CD2")
         switchName("flippedRes", "CE1", "CE2")
        elif g.resn=='GLN':
         switchName("flippedRes", "OE1", "NE2")
        elif g.resn=='GLU':
         switchName("flippedRes", "OE1", "OE2")
        elif g.resn=='LEU':
         switchName("flippedRes", "CD1", "CD2")
        elif g.resn=='ASN':
         switchName("flippedRes", "OD1", "ND2")
        elif g.resn=='TYR':
         switchName("flippedRes", "CD1", "CD2")
         switchName("flippedRes", "CE1", "CE2")
        elif g.resn=='VAL':
          switchName("flippedRes", "CG1", "CG2")
      cmd.sort()
    # cmd.label("flippedRes","name")
      return "flippedRes"
    
      
    
    #main function
    def rmsdByRes(referenceProteinChain,sel, targetProteinChain):
      """
      Update
        Zhenting Gao on 7/28/2016
    
      USAGE
    
        rmsf referenceProteinChain, targetProteinChain, selection [,byres=0], [reference_state=1]
    
        Calculate the RMSD for each residue pairs from two chains of the same protein from two crystal structures.
    
      Workflow
        Read reference and target pdb files
        Align two structures
            sel target, proA and chain A
                #define target protein chain
            sel refrence, proB and chain A
                #define reference protein chain
            align target, reference
                #automatical alignment
        Clean attributes
            otherwise rms_cur will fail
    
    
      """
    # Create temporary objects, exclude alternative conformation B
      cmd.create("ref_gzt", referenceProteinChain+" and polymer and not alt B")
      cmd.alter("ref_gzt", "chain='A'")
      cmd.alter("ref_gzt", "segi=''")
      cmd.create("target_gzt", targetProteinChain+" and polymer and not alt B")
      cmd.alter("target_gzt", "chain='A'")
      cmd.alter("target_gzt", "segi=''")
    #  cmd.align("target_gzt","ref_gzt",object="align")
    # parameters
      outputText=""
      res2Check=['HIS','ASP','ARG','PHE','GLN','GLU','LEU','ASN','TYR','VAL']
           
    # select alpha carbon of selected residues in reference structure
      calpha=cmd.get_model(sel+" and name CA and not alt B")
      
      for g in calpha.atom:
    #  print g.resi+g.resn
       if cmd.count_atoms("ref_gzt and polymer and resi "+g.resi)==cmd.count_atoms("target_gzt and polymer and resi "+g.resi):
        rmsdRes=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi,"target_gzt and polymer and resi "+g.resi)
        rmsdResCa=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi+" and name ca","target_gzt and polymer and resi "+g.resi+" and name ca")
        rmsdResBackbone=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi+" and name ca+n+c+o","target_gzt and polymer and resi "+g.resi+" and name ca+n+c+o")
    #  calculate minimum rmsd
        rmsdResMin=rmsdRes
        if g.resn in res2Check:
          flippedRes=flipAtomName("target_gzt and polymer and resi "+g.resi)
          rmsdFlippedRes=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi,flippedRes)
          if rmsdFlippedRes<rmsdRes:
           rmsdResMin=rmsdFlippedRes
    
    #    print cmd.count_atoms("ref_gzt and polymer and resi "+g.resi),cmd.count_atoms("target_gzt and polymer and resi "+g.resi)
        outputText+="%s,%s,%s,%.3f,%.3f,%.3f,%.3f\n" % (targetProteinChain,g.resn,g.resi,rmsdRes,rmsdResCa,rmsdResBackbone,rmsdResMin)
      
      print outputText
    # Destroy temporary objects
      cmd.delete("ref_gzt target_gzt align res_gzt "+flippedRes)
      
    # Save data into csv
      outputFile='rmsdByRes_'+sel+'.csv'
      f=open(outputFile,'a')
      if not is_non_zero_file(outputFile):
       f.write("targe,residueName,residueId,allAtomRMSD,rmsdResCa,rmsdResBackbone,allAtomRMSDMin\n")
      f.write(outputText)
      f.close()
      print "Results saved in "+outputFile
      
    cmd.extend("rmsdByRes",rmsdByRes)
    

Retrieved from "[https://pymolwiki.org/index.php?title=RmsdByResidue&oldid=12811](https://pymolwiki.org/index.php?title=RmsdByResidue&oldid=12811)"


---

## Rotamer Toggle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 DESCRIPTION
  * 2 IMAGES
  * 3 SETUP
  * 4 NOTES / STATUS
  * 5 USAGE
  * 6 EXAMPLES
  * 7 REFERENCES
  * 8 SCRIPTS (Rotamers.py ; MyMenu.py)



### DESCRIPTION

Backbone-Dependent Rotamer library (Dunbrack, Cohen ; see ref) is imported into pymol giving access to this information. There are a number of different ways to use the data, I've only implemented a few as well as added extra functions that seemed useful. 

  * Rotamer Menu - an added menu into menu.py, which displays the most common rotamers for the given(clicked) residue; you can also set the residue any of the common rotamers as well
  * colorRotamers - color rotamers by closest matching rotamer angles from database; i.e. color by how common each rotamer of selection, blue - red (least to most common).
  * set_rotamer - routine called by above menu, but can be called manually to set a specific residues side-chain angles
  * set_phipsi - set all phi,psi angles of given selection to given angles (useful for creating secondary structures)
  * createRotamerPDBs - create pdb for each rotamer of given selection ; filter by rotamer-probability



### IMAGES

  * [![Rotamer Menu for a GLN residue](/images/d/d1/RotamerMenu.png)](/index.php/File:RotamerMenu.png "Rotamer Menu for a GLN residue")

Rotamer Menu for a GLN residue 

  * [![Rotamer Comparison of crystal structure and most common for GLU; just as an example](/images/5/52/GLURotamerComparison5.png)](/index.php/File:GLURotamerComparison5.png "Rotamer Comparison of crystal structure and most common for GLU; just as an example")

Rotamer Comparison of crystal structure and most common for GLU; just as an example 




Print out while selecting most common rotamer from above-left image (GLN residue): 
    
    
     Given GLN:40 PHI,PSI (-171.626373291,-96.0500335693) : bin (-170,-100)
     CHIs: [179.18069458007812, 72.539344787597656, -47.217315673828125]
     Setting Chi1 to -176.9
     Setting Chi2 to 177.4
     Setting Chi3 to 0.7
    

### SETUP

run "rotamers.py" and use functions from commandline. 

or 

To setup a rotamer menu inside the residue menu (default windows pymol installation): 

  * copy rotamers.py to C:/Program Files/DeLano Scientific/PyMol/modules/pymol/rotamers.py
  * copy mymenu.py to C:/Program Files/DeLano Scientific/PyMol/modules/pymol/menu.py (WARNING : overwrites default menu.py - use at your own risk)
  * copy bbdep02.May.sortlib to C:/Program Files/DeLano Scientific/PyMol/modules/pymol/bbdep02.May.sortlib (or newer version of sorted bbdep)



This is only one possible way to do this, I am sure there are many others. I'm not going to post the bbdep, but there is a link in the References section to Dunbrack's download page (get the "sorted" lib) 

### NOTES / STATUS

  * Tested on Pymolv0.97, Windows platform, Red Hat Linux 9.0 and Fedora Core 4. Will test v0.98 and MacOSX later on.
  * The way it's setup now, when you import rotamers , it will automatically read-in the rotamer database; this may not be what you want.
  * Post problems in the discussion page, on 'my talk' page or just email me : dwkulp@mail.med.upenn.edu



TASKS TODO: 

  * Rotamer Movie, using mset, etc create movie to watch cycle through rotamers
  * Code could be organized a bit better; due to time constraints this is good for now..



TASKS DONE: 

  * Store crystal structure in rotamer menu, so you can go back to original orientation



### USAGE
    
    
    colorRotamers selection
    set_rotamer selection, chi1_angle [,chi2_angle] [,chi3_angle] [,chi4_angle]
    set_phipsi selection phi_angle, psi_angle
    createRotamerPBDs selection [,ncutoff] [,pcutoff] [,prefix]
    

### EXAMPLES
    
    
      colorRotamers chain A
      set_rotamer resi 40, -60,-40   (only set chi1,chi2 angles)
      set_phipsi resi 10-40, -60,-60 (create an alpha-helical-like section)
      createRotamerPDBs resi 10-12, ncutoff=3 (create 9 PDBs; each with one of the 3 most probable rotamers for resi 10,11,12)
      createRotamerPDBs resi 14, pcutoff=0.4  (create a pdb file for each rotamer of residue 14 with probablity > 0.4)
    

### REFERENCES

Dunbrack and Cohen. Protein Science 1997 

[Dunbrack Lab Page (Contains backbone-dependent library)](http://dunbrack.fccc.edu/bbdep/index.php)

### SCRIPTS (Rotamers.py ; MyMenu.py)

Rotamers.py 
    
    
    ##################################################################
    # File:          Rotamers.py
    # Author:        Dan Kulp
    # Creation Date: 6/8/05
    # Contact:       dwkulp@mail.med.upenn.edu
    #
    # Notes:
    #     Incorporation of Rotamer library
    #     readRotLib() - fills rotdat; 
    #        indexed by "RES:PHI_BIN:PSI_BIN".
    #
    #     Three main functions:
    #     1. colorRotamers - colors according
    #          to rotamer probablitity
    #     2. getBins(sel)
    #           phi,psi bin for rotamer
    #     3. set_rotamer - set a side-chain 
    #           to a specific rotamer
    #
    #     To setup a rotamer menu in the 
    #   right click, under "Residue"
    #        1. cp mymenu.py modules/pymol/menu.py
    #        2. cp rotamers.py modules/pymol/rotamers.py (update ROTLIB)
    #
    # Requirements:
    #  set ROTLIB to path for rotamer library
    # Reference: 
    #  Dunbrack and Cohen. Protein Science 1997
    ####################################################################
     
    import colorsys,sys
    import editing
    import os
    import cmd
    import math
     
    # Path for library
    ROTLIB=os.environ['PYMOL_PATH']+"/modules/pymol/bbdep02.May.sortlib"
     
    # Place for library in memory..
    rotdat = {}
     
    def readRotLib():
        # Column indexes in rotamer library..
        RES  = 0
        PHI  = 1
        PSI  = 2
        PROB = 8
        CHI1 = 9
        CHI2 = 10
        CHI3 = 11
        CHI4 = 12
     
        if os.path.exists(ROTLIB):
                    print "File exists: "+ROTLIB
                    input = open(ROTLIB, 'r')
                    for line in input:
     
                        # Parse by whitespace (I believe format is white space and not fixed-width columns)
                        dat = line.split()
     
                        # Add to rotamer library in memory : 
                        #   key format       RES:PHI_BIN:PSI_BIN
                        #   value format     PROB, CHI1, CHI2, CHI3, CHI4
                        key=dat[RES]+":"+dat[PHI]+":"+dat[PSI]
                        if key in rotdat:
                            rotdat[key].append([ dat[PROB], dat[CHI1], dat[CHI2], dat[CHI3], dat[CHI4] ])
                        else:
                            rotdat[key] = [ [ dat[PROB], dat[CHI1], dat[CHI2], dat[CHI3], dat[CHI4] ] ]
     
     
        else:
            print "Couldn't find Rotamer library"
     
     
    # Atoms for each side-chain angle for each residue
    CHIS = {}
    CHIS["ARG"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD" ],
                    ["CB","CG","CD","NE" ],
                    ["CG","CD","NE","CZ" ]
                  ]
     
    CHIS["ASN"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","OD2" ]
                  ]
     
    CHIS["ASP"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","OD1" ]
                  ]
    CHIS["CYS"] = [ ["N","CA","CB","SG" ]
                  ]
     
    CHIS["GLN"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD" ],
                    ["CB","CG","CD","OE1"]
                  ]
     
    CHIS["GLU"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD" ],
                    ["CB","CG","CD","OE1"]
                  ]
     
    CHIS["HIS"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","ND1"]
                  ]
     
    CHIS["ILE"] = [ ["N","CA","CB","CG1" ],
                    ["CA","CB","CG1","CD1" ]
                  ]
     
    CHIS["LEU"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD1" ]
                  ]
     
    CHIS["LYS"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD" ],
                    ["CB","CG","CD","CE"],
                    ["CG","CD","CE","NZ"]
                  ]
     
    CHIS["MET"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","SD" ],
                    ["CB","CG","SD","CE"]
                  ]
     
    CHIS["PHE"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD1" ]
                  ]
     
    CHIS["PRO"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD" ]
                  ]
     
    CHIS["SER"] = [ ["N","CA","CB","OG" ]
                  ]
     
    CHIS["THR"] = [ ["N","CA","CB","OG1" ]
                  ]
     
    CHIS["TRP"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD1"]
                  ]
     
    CHIS["TYR"] = [ ["N","CA","CB","CG" ],
                    ["CA","CB","CG","CD1" ]
                  ]
     
    CHIS["VAL"] = [ ["N","CA","CB","CG1" ]
                  ]
     
    # Color Rotamer by side-chain angle position
    #  'bin' side-chain angles into closest
    def colorRotamers(sel):
        doRotamers(sel)
     
    # Utility function, to set phi,psi angles for a given selection
    # Note: Cartoon, Ribbon functionality will not display correctly after this
    def set_phipsi(sel, phi,psi):
        doRotamers(sel,angles=[phi,psi],type="set")
     
    # Set a rotamer, based on a selection, a restype and chi angles
    def set_rotamer(sel, chi1, chi2=0,chi3=0,chi4=0):
        at = cmd.get_model("byres ("+sel+")").atom[0]
     
        list = [chi1,chi2,chi3,chi4]
        for i in range(len(CHIS[at.resn])):
            print "Setting Chi"+str(i+1)+" to "+str(list[i])
            editing.set_dihedral(sel + ' and name '+CHIS[at.resn][i][0],
                                 sel + ' and name '+CHIS[at.resn][i][1],
                                 sel + ' and name '+CHIS[at.resn][i][2],
                                 sel + ' and name '+CHIS[at.resn][i][3], str(list[i]))
     
        # Remove some objects that got created
        cmd.delete("pk1")
        cmd.delete("pk2")
        cmd.delete("pkmol")
     
    # Get Phi,Psi bins for given selection
    # WARNING:  assume selection is single residue (will only return first residue bins)
    def getBins(sel):
        return doRotamers(sel, type="bins")
     
    # Color Ramp...
    def rot_color(vals): 
            nbins = 10
            vals.sort(key=lambda x:x[1])
    #       print "End sort: "+str(len(vals))+" : "+str(nbins)
     
     
            # Coloring scheme...
            j = 0
            rgb = [0.0,0.0,0.0]
            sel_str = ""
            for i in range(len(vals)):
                    if int(len(vals)/nbins) == 0 or i % int(len(vals)/nbins) == 0:
                          hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), 1.0, 1.0)
     
                          #convert to rgb and append to color list
                          rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
                          if j < nbins-1:
                                  j += 1
     
                    cmd.set_color("RotProbColor"+str(i), rgb)
                    cmd.color("RotProbColor"+str(i), str(vals[i][0]))
     
     
    # Main function
    def doRotamers(sel,angles=[], type="color"):
     
            # Read in Rotamer library if not already done
            if len(rotdat) == 0:
                    readRotLib()
     
            # Set up some variables..
            residues = ['dummy']  # Keep track of residues already done
            probs = []            # probability of each residue conformation
            phi = 0               # phi,psi angles of current residue
            psi = 0
     
            # Get atoms from selection
            atoms = cmd.get_model("byres ("+sel+")")
     
            # Loop through atoms in selection
            for at in atoms.atom:
                try:
                   # Don't process Glycines or Alanines
                   if not (at.resn == 'GLY' or at.resn == 'ALA'):
                    if at.chain+":"+at.resn+":"+at.resi not in residues:
                        residues.append(at.chain+":"+at.resn+":"+at.resi)
     
                        # Check for a null chain id (some PDBs contain this) 
                        unit_select = ""
                        if at.chain != "":
                            unit_select = "chain "+str(at.chain)+" and "
     
                        # Define selections for residue i-1, i and i+1
                        residue_def = unit_select+'resi '+str(at.resi)
                        residue_def_prev = unit_select+'resi '+str(int(at.resi)-1)
                        residue_def_next = unit_select+'resi '+str(int(at.resi)+1)
     
                        # Compute phi/psi angle
     
                        phi = cmd.get_dihedral(residue_def_prev+' and name C',residue_def+' and name N',residue_def+' and name CA',residue_def+' and name C')
                        psi = cmd.get_dihedral(residue_def+' and name N',residue_def+' and name CA',residue_def+' and name C',residue_def_next+' and name N')
                        if type == "set":
                                print "Changing "+at.resn+str(at.resi)+" from "+str(phi)+","+str(psi)+" to "+str(angles[0])+","+str(angles[1])
                                cmd.set_dihedral(residue_def_prev+' and name C',residue_def+' and name N',residue_def+' and name CA',residue_def+' and name C',angles[0])
                                cmd.set_dihedral(residue_def+' and name N',residue_def+' and name CA',residue_def+' and name C',residue_def_next+' and name N', angles[1])
                                continue
     
                        # Find correct 10x10 degree bin                                     
                        phi_digit = abs(int(phi)) - abs(int(phi/10)*10)
                        psi_digit = abs(int(psi)) - abs(int(psi/10)*10)
     
                        # Remember sign of phi,psi angles
                        phi_sign = 1
                        if phi < 0:    phi_sign = -1
     
                        psi_sign = 1
                        if psi < 0:    psi_sign = -1
     
                        # Compute phi,psi bins
                        phi_bin = int(math.floor(abs(phi/10))*10*phi_sign)
                        if phi_digit >= 5:    phi_bin = int(math.ceil(abs(phi/10))*10*phi_sign)
     
                        psi_bin = int(math.floor(abs(psi/10))*10*psi_sign)
                        if psi_digit >= 5:    psi_bin = int(math.ceil(abs(psi/10))*10*psi_sign)
     
                        print "Given "+at.resn+":"+at.resi+" PHI,PSI ("+str(phi)+","+str(psi)+") : bin ("+str(phi_bin)+","+str(psi_bin)+")"
     
     
                        # Get current chi angle measurements
                        chi = []
                        for i in range(len(CHIS[at.resn])):
                           chi.append(cmd.get_dihedral(residue_def + ' and name '+CHIS[at.resn][i][0],
                                                         residue_def + ' and name '+CHIS[at.resn][i][1],
                                                         residue_def + ' and name '+CHIS[at.resn][i][2],
                                                         residue_def + ' and name '+CHIS[at.resn][i][3]))
                        print "CHIs: "+str(chi)
                        if type == 'bins':
                             return [at.resn, phi_bin,psi_bin]
     
                        # Compute probabilities for given chi angles
                        prob = 0
                        prob_box = 22                   
                        for item in range(len(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)])):
                            print "Rotamer from db: "+str(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item])
                            if chi[0]:
                                if chi[0] >= float(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][1]) - (prob_box/2) and \
                                    chi[0] <= float(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][1]) + (prob_box/2):
                                    if len(chi) == 1:
                                            prob = rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][0]
                                            break
                                    if chi[1] >= float(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][2]) - (prob_box/2) and \
                                     float(chi[1] <= rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][2]) + (prob_box/2):
                                            if len(chi) == 2:
                                                prob = rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][0]
                                                break
                                            if chi[2] >= float(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][3]) - (prob_box/2) and \
                                               float(chi[2] <= rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][3]) + (prob_box/2):
                                                if len(chi) == 3:
                                                    prob = rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][0]
                                                    break
                                                if chi[3] >= float(rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][4]) - (prob_box/2) and \
                                                   float(chi[3] <= rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][4]) + (prob_box/2):
                                                    prob = rotdat[at.resn+":"+str(phi_bin)+":"+str(psi_bin)][item][0]
                                                    break
     
     
                        print "PROB OF ROTAMER: "+str(prob)
                        print "---------------------------"
                        probs.append([residue_def, prob])
     
                except:
    #               probs.append([residue_def, -1])
                    print "Exception found"
                    continue
     
            # Color according to rotamer probability
            rot_color(probs)
     
     
     
     
    #  Create PDB files containing most probable rotamers
    def createRotamerPDBs(sel,ncutoff=10,pcutoff=0,prefix="ROTAMER"):
     
            # Get atoms from selection
            atoms = cmd.get_model("byres ("+sel+")")
     
            # Set up some variables..
            residues = ['dummy']  # Keep track of residues already done
     
            # Loop through atoms in selection
            for at in atoms.atom:
                    if at.resn in ('GLY','ALA') or "%s:%s:%s" % (at.chain,at.resn,at.resi) in residues:
                            continue
     
                    # Add to residue list (keep track of which ones we've done)
                    residues.append("%s:%s:%s" % (at.chain,at.resn,at.resi))
     
                    # Check for a null chain id (some PDBs contain this)
                    unit_select = ""
                    if not at.chain == "":
                        unit_select = "chain "+str(at.chain)+" and "
     
                    # Define selections for residue 
                    residue_def = unit_select+'resi '+str(at.resi)
     
                    # Get bin (phi,psi) definitions for this residue
                    bin = doRotamers(residue_def, type='bins')
     
                    # Store crystal angle
                    crystal_angles = [0.0,0.0,0.0,0.0]
                    for angle in range(3):
                            try:
                                    crystal_angles[angle] = bin[3][angle]
                            except IndexError:
                                    break
     
                    # Retreive list of rotamers for this phi,psi bin + residue type
                    match_rotamers = rotdat["%s:%s:%s" % (bin[0],str(bin[1]),str(bin[2]))]
     
                    count = 0
                    for item in range(len(match_rotamers)):
     
                            # Store probablity
                            prob = match_rotamers[item][0]
     
                            # Check cutoffs
                            if float(prob) <= float(pcutoff):
                                    continue
     
                            if float(count) >= float(ncutoff):
                                    break
     
                            # Increment count
                            count += 1
     
                            # Output to screen ...
                            print "Residue %s%s, rotamer %i, prob %s" % (str(at.resn),str(at.resi),int(item),str(prob))
     
                            # Set to new rotamer
                            set_rotamer(residue_def,match_rotamers[item][1],match_rotamers[item][2],match_rotamers[item][3],match_rotamers[item][4])                                                                                                
     
                            # Store in PDB file
                            cmd.save("%s_%s%s_%i_%s.pdb" % (prefix,str(at.resn),str(at.resi),int(item),str(prob)))
     
                            # Reset crystal angle
                            set_rotamer(residue_def,crystal_angles[0],crystal_angles[1],crystal_angles[2],crystal_angles[3])
     
    # Uncommenting this is nice because it loads rotamer library upon startup
    #  however, it slows the PyMOL loading process a lot
    #  instead I've put this call into the menuing code..
    # readRotLib()
     
    cmd.extend('set_phipsi',set_phipsi)
    cmd.extend('set_rotamer',set_rotamer)
    cmd.extend('colorRotamers',colorRotamers)
    cmd.extend('createRotamerPDBs',createRotamerPDBs)
    

MyMenu.py 
    
    
    Since menu.py is copyrighted I can't post my edited version, but you can create it very simply by adding these two peices of code
    

1\. In the "pick_option(title,s,object=0)" function of menu.py add the following code after the first "result =" statement 
    
    
    # Edit dwkulp 6/11/05 , add a rotamer menu to residue menu
       if title == 'Residue':
    	result.extend([[ 1, 'rotamers'    , rotamer_menu(s)]])
    

2\. At the end of the file add this: 
    
    
    ###############################################
    # Dan Kulp
    # Added Rotamer list to residue menu..
    # rotamer.py must be importable (I placed it in 
    # the same directory as menu.py)
    ###############################################
    
    import rotamers
    
    
    def rotamer_menu(s):
    	# Check for rotamer library being loaded
    	if not rotamers.rotdat:
                 rotamers.readRotLib()
    #	     return [ [2, "Must run rotamers.py first",'']]
    
    	# Check for valid rotamer residue..
    	res = cmd.get_model("byres ("+s+")").atom[0].resn
            if not res in rotamers.CHIS.keys():
    	    return [ [2, "Residue: "+res+" not known sidechain or does not have rotamers", '']]
    
    	# Get PHI,PSI bins for rotamer (also prints out current phi,psi, chi1,chi2,chi3,chi4)
    	bins = rotamers.doRotamers(s,type='bins')
    
    	# Add a title to the menu
    	result = [ [2, bins[0]+' Rotamers in bin ('+str(bins[1])+','+str(bins[2])+')','' ], [1, ':::PROB,CHI1,CHI2,CHI3,CHI4:::','']]
    
            # Grab the entries for this residue and phi,psi bins
    	match_rotamers = rotamers.rotdat[bins[0]+":"+str(bins[1])+":"+str(bins[2])]
    
    	# Set max number of rotamers to display (probably should be somewhere 'higher up' in the code)
    	max_rotamers = min(10, len(match_rotamers))
    
    	# Create menu entry for each possible rotamer
            for item in range(max_rotamers):
                 result.append( [ 1, str(match_rotamers[item]), 'rotamers.set_rotamer("'+s+'","'\
    										    +str(match_rotamers[item][1])+'","'\
    										    +str(match_rotamers[item][2])+'","'\
    										    +str(match_rotamers[item][3])+'","'\
    										    +str(match_rotamers[item][4])+'")'])
    	return result
    

Retrieved from "[https://pymolwiki.org/index.php?title=Rotamer_Toggle&oldid=8412](https://pymolwiki.org/index.php?title=Rotamer_Toggle&oldid=8412)"


---

## RotationAxis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/draw_rotation_axis.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/draw_rotation_axis.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO cylinder representing a rotation axis for a given transformation. It's very useful for drawing the axes of rotational symmetry in an oligomeric assembly. 

The idea is to align two molecules/domains/chains/selections (using cmd.super) and extract the trasformation (TTT) matrix (T). The direction of the rotation axis, and a point are obtained from T and are used to create a cgo object representing the axis. The script generates two measures: one in the graphical screen (the rotation axis value, and the norm of the translation vector along the rotation axis) and some basic information in the command-line (the transformation matrix, the rotation angle, distance between centroids, and some pml lines that you can use to reproduce the axis...) 

As always with these type of things, you have to use it at your own risk. I did not try all possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado@gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run _path-to-the-script_ /draw_rotation_axis.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_axis('selection1', 'selection2')**

Please, pay attention to the apostrophes around the selections, you MUST use them. Also works with chains: 

**draw_axis('chain A', 'chain B')**

or objects: 

**draw_axis('obj01', 'obj02')**

Also, you can play a bit with the length, width and color of the axis. 

**draw_axis('selection1', 'selection2', scale_factor, width, r1, g1, b1, r2, g2, b2)**

**scale_factor** = to control the length of the axis, the default is 20 

**width** = to control the width of the axis. Default is 0.6 

**r1, g1, b1** = first RGB color code. Default is 1, 1, 1 

**r2, g2, b2** = second RGB color code to create the gradient. Default is 1, 0, 0. 

To create a single color axis, just made r1,g1,b1=r2,g2,b2 

# Examples
    
    
    # download the source and save as draw_rotation_axis.py
    run draw_rotation_axis.py
    fetch 2vak
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/d/df/2vak4.png)](/index.php/File:2vak4.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/3/37/2vak2.png)](/index.php/File:2vak2.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Another example, calculate the rotation axis of an homodimer
    run draw_rotation_axis.py
    fetch 3rfg
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/5/55/3rfg1.png)](/index.php/File:3rfg1.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/d/d6/3rfg4.png)](/index.php/File:3rfg4.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Clearly, each of the domains are made up with motifs with internal symmetry
    draw_axis('resi 200-236 and chain A', 'resi 238-274 and chain A', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # Also, you can create first the selections and use them to calculate the axis
    sele selection1, resi 200-236 and chain A
    sele selection2, resi 238-274 and chain A
    draw_axis('selection1', 'selection2', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # will produce the same result.
    

  * [![Internal rotation axis](/images/6/60/3rfg2.png)](/index.php/File:3rfg2.png "Internal rotation axis")

Internal rotation axis 

  * [![Some basic information is printed in the screen](/images/0/03/3rfg3.png)](/index.php/File:3rfg3.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 




Retrieved from "[https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923](https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923)"


---

## Show contacts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/show_contacts.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/show_contacts.py)  
Author(s)  | [David Ryan Koes](/index.php?title=User:DavidKoes&action=edit&redlink=1 "User:DavidKoes \(page does not exist\)")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments



## Introduction

PyMOL Plugin for displaying polar contacts. Good hydrogen bonds (as determined by PyMOL) are shown in yellow. Electrostatic clashes (donor-donor or acceptor-acceptor) are shown in red. Close (<4.0 A) but not ideal contacts are shown in purple. Cutoffs are configurable. Exports the command `contacts` which takes two selections and an optional name for the generated contacts group. Alternatively, the selections can be chosen using a dialog box accessible from the Plugins menu. 

  


## Usage
    
    
    show_contacts sel1, sel2, result="contacts", cutoff=3.6, bigcutoff=4.0
    

## Required Arguments

  * **sel1** = first selection
  * **sel2** = second selection



  


## Optional Arguments

  * **result** = name of created group containing contacts
  * **cutoff** = cutoff for good contacts (yellow dashes)
  * **bigcutoff** = cutoff for suboptimal contacts (purple dashes)



Retrieved from "[https://pymolwiki.org/index.php?title=Show_contacts&oldid=13265](https://pymolwiki.org/index.php?title=Show_contacts&oldid=13265)"


---

## Nmr cnstr

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/nmr_cnstr.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/nmr_cnstr.py)  
Author(s)  | [Evangelos Papadopoulos](/index.php?title=User:Evangelos&action=edit&redlink=1 "User:Evangelos \(page does not exist\)")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  
This script will display the NMR constrains used for a structure calculation atop a structure. Usage: Save this as "NMRcnstr.py" load your protein in PyMOL, and run the script. type upl('fname') or cns('fname') where fname is the filename with the NMR constrains you want to display. It is still a very preliminary version. 

If you generated the structure by CYANA type: 

cyana> read final.pdb 

to input the structure in cyana then: 

cyana>pseudo=1 

before exporting the structure again by: 

cyana> write final.pdb 

this way the structure will contain the appropriate pseudoatoms nomeclature. 

Welcome to contact me if you need some help to set it up. 

[![Cns.png](/images/b/ba/Cns.png)](/index.php/File:Cns.png) [![Upl.png](/images/8/8f/Upl.png)](/index.php/File:Upl.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Nmr_cnstr&oldid=13915](https://pymolwiki.org/index.php?title=Nmr_cnstr&oldid=13915)"


---

## ShowLigandWaters

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [Gianluca Tomasello](/index.php?title=User:GianlucaTomasello&action=edit&redlink=1 "User:GianlucaTomasello \(page does not exist\)")  
License  | [CC BY-NC-SA](http://creativecommons.org/licenses/by-nc-sa/3.0)  
  
  


## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Examples
  * 5 The Code



## Introduction

This script detects waters molecules within a specified distance from the ligand.  
Detected water molecules are shown as spheres.  
Distances between water molecules and O or N atoms of ligand (potential H-bonds) are shown by dotted lines.  
An output file containing a list of distance between waters and ligand atoms and the number of interactions is written in the main PyMol folder. 

  


## Usage

waters [ligand name, distance] 

  


## Required Arguments

  * **ligand name** = the ligand residue name
  * **distance** = max distance in Angstroms



  


## Examples

**example #1 on PDB structure (1C0L)**
    
    
    PyMOL>waters FAD, 2.8
    

Output in GUI of PyMol: 

Total number of water molecules at 2.8 A from ligand FAD: 3 Number of water molecules that interact with ligand: 3 Number of interactions between water molecules and ligand: 4 

output file: waters.txt 

HOH residue 2002 -- and -- ['FAD', '1363', 'O1P'] ---> 2.611289 A HOH residue 2002 -- and -- ['FAD', '1363', 'O5B'] ---> 3.592691 A HOH residue 2008 -- and -- ['FAD', '1363', 'O1A'] ---> 2.678604 A HOH residue 2009 -- and -- ['FAD', '1363', 'O2P'] ---> 2.643039 A 

* * *

Total number of water molecules at 2.8 A from ligand FAD: 3 Number of water molecules that interact with ligand: 3 Number of interactions between water molecules and ligand: 4 

[![example #1](/images/e/e7/Example_1_FAD_1C0L.png)](/index.php/File:Example_1_FAD_1C0L.png "example #1")

  


## The Code

Copy the following text and save it as waters.py 
    
    
    # -*- coding: cp1252 -*-
    """
    This script detects waters molecules within a specified distance from the ligand.
    Water molecules are shown.
    Distance between water molecules and O or N atoms of ligand are shown.
    
    Author: Gianluca Tomasello, Gianluca Molla
    University of Insubria, Varese, Italy
    09/02/2013
    gianluca.molla@uninsubria.it
    
    Usage
    waters, [ligand name, distance]
    
    Parameters
    
    ligand : the ligand residue name 
    
    distance : a float number that specify the maximum distance from the ligand to consider the water molecule
               
    Output
    -A file is produced containing a list of distance between waters and ligand atoms and the number of interactions
    -A graphic output show the water molecules interacting whith the ligand atoms by showing the distances between them
    """
    
    from pymol import cmd,stored
    #define function: waters
    def waters(ligand, distance):
        stored.HOH = [] 
        cmd.set('label_color','white') 
        cmd.delete ("sele")
        cmd.hide ("everything")
        cmd.show_as ("sticks", "resn %s" % ligand)
        #iterate all water molecules nearby the ligand
        cmd.iterate('resn %s around %s and resn HOH' % (ligand,distance),'stored.HOH.append(resi)') 
        f = open("waters.txt","a")
        count=0
        count_int=0
        
        for i in range(0,len(stored.HOH)):        
            cmd.distance('dist_HOH_FAD', 'resi ' + stored.HOH[i], '(resn %s and n. O*+N*) w. 3.6 of resi %s'% (ligand, stored.HOH[i]))
            stored.name = []
            #iterate all ligand atoms within a predetermined distance from the water molecule
            cmd.iterate('(resn %s and n. O*+N*) w. 3.6 of resi %s'% (ligand, stored.HOH[i]),'stored.name.append([resn,resi,name])') 
    
            if stored.name:           
                count = count+1
                count_int = count_int+len(stored.name)
    
            for j in range(0,len(stored.name)):            
                cmd.select('base', 'resi ' + stored.HOH[i])
                cmd.select('var','resn '+ligand+ ' and n. ' + stored.name[j][2]) 
                #calculate the distance between a specific atom and the water molecule
                dist = cmd.get_distance('base','var')            
                f.write('HOH residue %s -- and -- %s  --->  %f A\n'%(stored.HOH[i],stored.name[j], dist))  
    
        cmd.select ("waters","resn %s around %s and resn HOH" % (ligand,distance))
        cmd.show_as ("spheres", "waters")
        cmd.zoom ("visible")
        num_atm = cmd.count_atoms ("waters")
        print ("Total number of water molecules at %s A from ligand %s: %s \n" % (distance,ligand,num_atm))
        print ("Number of water molecules that interact with ligand: %d\n" % (count))
        print ("Number of interactions between water molecules and ligand: %d\n" % count_int)
        f.write('-------------------------------------------------------------\n')
        f.write("Total number of water molecules at %s A from ligand %s: %s \n" % (distance,ligand,num_atm))
        f.write("Number of water molecules that interact with ligand: %d\n" % count)
        f.write("Number of interactions between water molecules and ligand: %d\n\n\n\n" % count_int)
        f.close()
        cmd.delete ("waters")
        
    cmd.extend("waters",waters)
    

Retrieved from "[https://pymolwiki.org/index.php?title=ShowLigandWaters&oldid=11415](https://pymolwiki.org/index.php?title=ShowLigandWaters&oldid=11415)"


---

## Sidechaincenters

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/e/ef/SidechaincentersExample.png)](/index.php/File:SidechaincentersExample.png)

[](/index.php/File:SidechaincentersExample.png "Enlarge")

sidechain centers

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.creating](http://pymol.org/psicoredirect.php?psico.creating)  
  
Pseudo single-atom representation of sidechains. Usefull for [pair potential calculation](/index.php/AAindex "AAindex") for example. 

# Example
    
    
    fetch 2x19
    sidechaincenters scc, 2x19
    

# The Script
    
    
    '''
    (c) 2010 Thomas Holder
    '''
    
    from pymol import cmd
    from chempy import Atom, cpv, models
    
    sidechaincenteratoms = {
        'GLY': ('CA',),
        'ALA': ('CB',),
        'VAL': ('CG1', 'CG2'),
        'ILE': ('CD1',),
        'LEU': ('CD1', 'CD2'),
        'SER': ('OG',),
        'THR': ('OG1', 'CG2'),
        'ASP': ('OD1', 'OD2'),
        'ASN': ('OD1', 'ND2'),
        'GLU': ('OE1', 'OE2'),
        'GLN': ('OE1', 'NE2'),
        'LYS': ('NZ',),
        'ARG': ('NE', 'NH1', 'NH2'),
        'CYS': ('SG',),
        'MET': ('SD',),
        'MSE': ('SE',),
        'PHE': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
        'TYR': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'),
        'TRP': ('CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3'),
        'HIS': ('CG', 'ND1', 'CD2', 'CE1', 'NE2'),
        'PRO': ('CB', 'CG', 'CD'),
    }
    
    def sidechaincenters(object='scc', selection='all', name='PS1', method='bahar1996'):
        '''
    DESCRIPTION
    
        Creates an object with sidechain representing pseudoatoms for each residue
        in selection.
    
        Sidechain interaction centers as defined by Bahar and Jernigan 1996
        http://www.ncbi.nlm.nih.gov/pubmed/9080182
    
    USAGE
    
        sidechaincenters object [, selection]
    
    ARGUMENTS
    
        object = string: name of object to create
    
        selection = string: atoms to consider {default: (all)}
    
        name = string: atom name of pseudoatoms {default: PS1}
    
    SEE ALSO
    
        sidechaincentroids, pseudoatom
        '''
        atmap = dict()
        if method == 'bahar1996':
            modelAll = cmd.get_model('(%s) and resn %s' % (selection, '+'.join(sidechaincenteratoms)))
            for at in modelAll.atom:
                if at.name in sidechaincenteratoms[at.resn]:
                    atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)
        elif method == 'centroid':
            modelAll = cmd.get_model('(%s) and not (hydro or name C+N+O)' % selection)
            for at in modelAll.atom:
                atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)
        else:
            print('Error: unknown method:', method)
            return
    
        model = models.Indexed()
        for centeratoms in atmap.values():
            center = cpv.get_null()
            for at in centeratoms:
                center = cpv.add(center, at.coord)
            center = cpv.scale(center, 1./len(centeratoms))
            atom = Atom()
            atom.coord = center
            atom.index = model.nAtom + 1
            atom.name = name
            for key in ['resn','chain','resi','resi_number','hetatm','ss','segi']:
                setattr(atom, key, getattr(at, key))
            model.add_atom(atom)
        model.update_index()
        if object in cmd.get_object_list():
            cmd.delete(object)
        cmd.load_model(model, object)
        return model
    
    def sidechaincentroids(object='scc', selection='all', name='PS1'):
        '''
    DESCRIPTION
    
        Sidechain centroids. Works like "sidechaincenters", but the
        pseudoatom is the centroid of all atoms except hydrogens and backbone atoms
        (N, C and O).
    
    NOTE
    
        If you want to exclude C-alpha atoms from sidechains, modify the selection
        like in this example:
    
        sidechaincentroids newobject, all and (not name CA or resn GLY)
    
    SEE ALSO
    
        sidechaincenters
        '''
        return sidechaincenters(object, selection, name, method='centroid')
    
    cmd.extend('sidechaincenters', sidechaincenters)
    cmd.extend('sidechaincentroids', sidechaincentroids)
    
    cmd.auto_arg[1].update({
        'sidechaincenters'     : [ cmd.selection_sc           , 'selection'       , ''   ],
        'sidechaincentroids'   : [ cmd.selection_sc           , 'selection'       , ''   ],
    })
    
    # vi: ts=4:sw=4:smarttab:expandtab
    

Retrieved from "[https://pymolwiki.org/index.php?title=Sidechaincenters&oldid=13468](https://pymolwiki.org/index.php?title=Sidechaincenters&oldid=13468)"


---

## Ss

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A command to list a summary of the secondary structure for a selection. Use like "ss my_protein" where "my_protein" is the name of the chain or structure in view. 
    
    
    from pymol import cmd
    from pymol import stored
    
    def ss(selection):
    
        class SSE(object):
    
            def __init__(self, start, typ):
                self.start, self.typ = start, typ
                self.end = -1
    
            def __repr__(self):
                return "%s-%s %s" % (self.start, self.end, self.typ)
    
        stored.pairs = []
        cmd.iterate("%s and n. ca" % selection, "stored.pairs.append((resi, ss))")
        num, currentType = stored.pairs[0]
    
        sses = [SSE(num, currentType)]
        currentSSE = sses[0]
        for resi, ssType in stored.pairs:
            if ssType == currentType:
                currentSSE.end = resi
            else:
                sses.append(SSE(resi, ssType))
                currentSSE = sses[-1]
                currentType = ssType
    
        for sse in sses:
            print sse
    
    cmd.extend("ss", ss)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ss&oldid=7116](https://pymolwiki.org/index.php?title=Ss&oldid=7116)"


---

## Sst

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.editing](http://pymol.org/psicoredirect.php?psico.editing)  
  
sst is a wrapper for the [SST secondary structural assignment web service](http://lcb.infotech.monash.edu.au/sstweb2/). The command updates PyMOL's **ss** atom property. 

SST does minimum message length inference of secondary structure based on C-alpha atom positions. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    sst [ selection [, raw [, state ]]]
    

## Arguments

  * **selection** = str: atom selection {default: all}
  * **raw** = str: [atom property](/index.php/Iterate#Exposed_Variables "Iterate") to assign the sst secondary structure type to. This type will also be translated to PyMOL's H/S/L types and assiged to the **ss** property {default: }
  * **state** = int: object state {default: -1 (current state)}



## Example
    
    
    import psico.editing
    fetch 1ubq, async=0
    
    sst all, raw=custom
    as cartoon
    
    color gray
    spectrum custom, selection=guide
    set cartoon_discrete_colors, 1
    label guide, custom
    

## See Also

  * [dssp](/index.php/Dssp "Dssp")
  * [dss](/index.php/Dss "Dss")



Retrieved from "[https://pymolwiki.org/index.php?title=Sst&oldid=12677](https://pymolwiki.org/index.php?title=Sst&oldid=12677)"


---

## Transformations

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/transformations.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/transformations.py)  
Author(s)  | Christoph Gohlke   
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Introduction

This script provides a library of functions for manipulation of matrices and vectors, useful in comparing alignments of various objects. It was written by Christoph Gohlke and is available along with many other scripts from his [website](http://www.lfd.uci.edu/~gohlke/). 

## Usage

The script must be imported in any other script that uses it. The easiest way to do this is to use the [Pymol-script-repo](/index.php/Git "Git"). 
    
    
    import transformations
    

## See also

The following scripts use transformations.py: 

  * [elbow_angle](/index.php/Elbow_angle "Elbow angle")



Retrieved from "[https://pymolwiki.org/index.php?title=Transformations&oldid=13933](https://pymolwiki.org/index.php?title=Transformations&oldid=13933)"


---

## Translate And Measure

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

To use, you would call it like : 
    
    
    print translateAndMeasure("molA", "molB", [1,0,0], 4)
    

which would print "overlap" if any of the atoms in molA or molB were within 4 Angstrom after translating by 1 along X. 

Of course, this could be improved to report exactly _which_ atoms were overlapping, or to make distance objects (using cmd.distance) to show them. 
    
    
    def translateAndMeasure(selection, other, translationVector, cutoff):
        cmd.translate(translationVector, selection)
        return checkDistances(selection, other, cutoff)
    
    def checkDistances(moleculeA, moleculeB, cutoff):
        ids_A = getIds(moleculeA)
        ids_B = getIds(moleculeB)
        for idA in ids_A:
            for idB in idsB:
                d = distance(moleculeA, idA, moleculeB, idB)
                if d > cutoff: return "overlap"
        return "no overlap"
    
    def distance(a, idA, b, idB):
        atomA = "%s and id %s" % (a, idA)
        atomB = "%s and id %s" % (b, idB)
        return cmd.get_distance(atomA, atomB)
    
    def getIds(selection):
        my_dict = { 'my_list' : [] }
        cmd.iterate(selection, "my_list.append(ID)", space=my_dict)
        return my_dict['my_list']
    

Retrieved from "[https://pymolwiki.org/index.php?title=Translate_And_Measure&oldid=6369](https://pymolwiki.org/index.php?title=Translate_And_Measure&oldid=6369)"


---

## WriteSS

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Instructions
  * 3 The Code
  * 4 Example



## Overview

This script will write the secondary structural elements for each PDB file in some specified directory, for each alpha carbon in the protein, to an output text file. Residues without secondary structure definition will get a "." in the string. 

This is untested code and has little chance of working for anyone but me. Although, it does work for me. 

## Instructions

  1. Copy the code to your machine
  2. Find **FIXME** in the code below and change "files" to either a glob like the one already there (for a whole directory) or a list with one element (for just one file).
  3. change the output filename. By default it's **PDBCODE.ss**.



## The Code
    
    
    import os
    import os.path
    import glob
    
    from pymol import cmd
    from pymol import stored
    from pymol import selector
    
    files = glob.glob("/tmp/thy_model/*.pdb")
    
    for file in files:
            pdbName = os.path.basename(file).split(".")[0]
            cmd.load(file, pdbName)
            outFile = open(pdbName + '.ss', 'wb')
            stored.ss = ""
            cmd.iterate( '(n. CA)', 'stored.ss = stored.ss + ("%1s"%ss)')
            for c in stored.ss:
                    if c  == " ":
                            outFile.write('.')
                    else:
                            outFile.write(c)
            cmd.delete(pdbName)
            outFile.close()
    

## Example
    
    
    python ss.pym  # my directory has 1d7p.pdb in it
    

output for 1D7P.pdb: 
    
    
    ...............HHH.......SSS..SSSHHHHH..................................................................SS..SS..................................................
    

Retrieved from "[https://pymolwiki.org/index.php?title=WriteSS&oldid=7668](https://pymolwiki.org/index.php?title=WriteSS&oldid=7668)"


---

## Xfit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.fitting](http://pymol.org/psicoredirect.php?psico.fitting)  
  
xfit does a weighted superposition of two structures. The weights are estimated with maximum likelihood. 

This typically gives better looking superpositions than the [align](/index.php/Align "Align") command, which does simple outlier rejection (equivalent to weight=0 on the rejected atoms, and weight=1 on the included atoms). 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Example
  * 5 See Also



## Installation

xfit is available from the [psico](/index.php/Psico "Psico") package and requires [CSB](https://github.com/csb-toolbox/CSB). 

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c speleo3 csb
    

## Usage
    
    
    xfit mobile, target [, mobile_state [, target_state [, load_b
        [, cycles [, match [, guide [, seed ]]]]]]]
    

## Arguments

  * **mobile** = string: atom selection
  * **target** = string: atom selection
  * **mobile_state** = int: object state of mobile selection {default: current}
  * **target_state** = int: object state of target selection {default: current}
  * **load_b** = 0 or 1: save -log(weights) into B-factor column {default: 0}
  * **cycles** = int: number of weight refinement cycles {default: 10}
  * **match** = align or super: alignment method {default: align}
  * **guide** = 0 or 1: use only CA-atoms (protein) or C4' (nucleic acid) {default: 1}
  * **seed** = 0 or 1: use initial weights from current positions {default: 0}



## Example
    
    
    import psico.fitting
    fetch 4akeA 1akeA, async=0
    
    xfit 4akeA, 1akeA, load_b=1
    
    spectrum b, blue_white_red, 4akeA & guide
    

## See Also

  * [intra_xfit](/index.php/Intra_xfit "Intra xfit")
  * [align](/index.php/Align "Align")



Retrieved from "[https://pymolwiki.org/index.php?title=Xfit&oldid=12658](https://pymolwiki.org/index.php?title=Xfit&oldid=12658)"


---

