# Category: Script Library

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

## AlphaToAll

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_**alphatoall** (all lower case) is a core command since PyMOL 1.7.2_

# Overview

AlphaToAll will expand the alpha-carbon property that you specify to all atoms in the residues. For example, if you set the b-factor column to some different value for all alpha carbons, this script will propagate those values into all the atoms as well. 

**In the example below, notice the color of the sticks.**

  * [![The alpha carbon here has a color property that does not match the rest of the residue.](/images/b/b7/Alphac0.png)](/index.php/File:Alphac0.png "The alpha carbon here has a color property that does not match the rest of the residue.")

The alpha carbon here has a color property that does not match the rest of the residue. 

  * [![This script will fix that, expanding the color attribute to the rest of the atoms.](/images/5/5c/Alphac1.png)](/index.php/File:Alphac1.png "This script will fix that, expanding the color attribute to the rest of the atoms.")

This script will fix that, expanding the color attribute to the rest of the atoms. 




# Usage
    
    
    # expand propertyName in objName on all the alpha carbons
    # to every atom, for each residue
    alphaToAll objName, propertyName
    
    # Example #1
    # For all objects, for each residue, set every atom's b-column value
    # to match that residue's alpha carbon's b-column value.
    alphaToAll *, b
    
    # Example #2
    # For selection 1qox and c. A expand the alpha carbon
    # colors to all atoms in the selection
    alphaToAll 1qox and c. A, color
    

# The Code
    
    
    from pymol import cmd, CmdException
    
    def alphaToAll(sel, col="b",forceRebuild=False):
    	"""
    	alphaToAll -- expand any property of the alpha carbons to all atoms in the residue
    	
    	PARAMS
    		sel
    			The selection or object (include "*" for all objects) to operate on.  This will
    			read the alpha carbon's "COL" property and expand that down to all atoms in the
    			given residues in sel.
    			
    		col
    			Any valid PyMOL property.  For example, 'b' or 'q' or 'color', etc.
    			DEFAULT: b, for B-factor as we most often overwrite that column
    			
    		forceRebuild
    			If a color changes, this will rebuild the objects for us.  For speed, we
    			DEFAULT this to False.
    			
    	RETURNS
    		None, just epxands the properties as dicsussed above.
    		
    	NOTES
    		This script is useful for instances where we update the b-factor column for the
    		alpha carbons and then want to expand that color to all atoms for smoother looking
    		images.
    		
    	AUTHOR:
    		Jason Vertrees, 2009.
    			
    	"""
    	col = '(' + ','.join(col.split()) + ')'
    	space = {'props': dict()}
    	cmd.iterate('byca (%s)' % (sel), 'props[model,segi,chain,resi] = ' + col, space=space)
    	cmd.alter(sel, col + ' = props.get((model,segi,chain,resi), ' + col + ')', space=space)
    	if forceRebuild != False:
    		cmd.rebuild(sel)
    
    cmd.extend("alphaToAll", alphaToAll)
    

Retrieved from "[https://pymolwiki.org/index.php?title=AlphaToAll&oldid=12886](https://pymolwiki.org/index.php?title=AlphaToAll&oldid=12886)"


---

## Annotate v

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/annotate_v.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/annotate_v.py)  
Author(s)  | [Xin Yu](/index.php/User:XinYu "User:XinYu")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Annotation Schemes
  * 3 Dependencies and Limitations
  * 4 How to use
  * 5 Contact



## Description

_Version: 1.0_

A simple-to-use script for annotation of VH or VL sequence of an antibody. It creates a Pymol object for each FR or CDR region. It utilizes the REST API interface of _Abnum_ (<http://www.bioinf.org.uk/abs/abnum/>) from Dr Andrew Martin's group at UCL 

## Annotation Schemes

Currently supports **Kabat, Chothia, Contact, IMGT** 1

1Definitions for Kabat, Chothia, Contact, and IIMGT are same as listed in the table (<http://www.bioinf.org.uk/abs/info.html#kabatnum>), except that for IMGT, H-CDR2 is defined as **H51-H57** in this script, as opposed to of **H51-H56** in the table. This slight revision generates result that matches that from _IMGT website_ (<http://www.imgt.org/>) 

## Dependencies and Limitations

1\. Import **request** module 

2\. Relies on internet connection to _Abnum_

3\. Incomplete VH or VL sequence might not be annotated 

## How to use

1\. Create a selection object for the V region, (such as **VH** and **VL** shown here). Only 1 V-region (VH OR VL) can be annotated each time. Alternatively, simply select the residues into the default **< sele>** group. 

  


[![](/images/8/82/Annotate_v_pic1.png)](/index.php/File:Annotate_v_pic1.png)

[](/index.php/File:Annotate_v_pic1.png "Enlarge")

Create a selection group for input seq, or use the default sele

2\. Copy the entire script and create a `.py` file. Run script. 

3\. The general syntax is: 
    
    
    annotate_v("selection_group_name", "scheme")
    

**selection_group_name:** name of the selection group where you have the input sequence. Must be single-letter amino acid coded. Either VH or VL but not both 

**scheme:** currently supports Kabat, Chothia, Contact, IMGT. Must be lowercase 

For example: 
    
    
    annotate_v("VH", "kabat") #input sequence from VH selection group, annotate using "kabat"
    annotate_v("sele", "imgt") #input sequence from the default sele group, annotate using "imgt". You can also try "contact", "chothia"
    

  
4\. In the output window, the script will print the FR and CDR regions (if found). It also automatically creates a selection group for each of the FR and CDR region. 

[![](/images/6/63/Annotate_v_pic2.png)](/index.php/File:Annotate_v_pic2.png)

[](/index.php/File:Annotate_v_pic2.png "Enlarge")

example output from command lines

[![](/images/7/72/Annotate_v_pic3.png)](/index.php/File:Annotate_v_pic3.png)

[](/index.php/File:Annotate_v_pic3.png "Enlarge")

new FR and CDR selection groups are created

  


## Contact

For bugs & questions, emai @ xinyu18018@gmail.com. Enjoy! 

Retrieved from "[https://pymolwiki.org/index.php?title=Annotate_v&oldid=13886](https://pymolwiki.org/index.php?title=Annotate_v&oldid=13886)"


---

## Apropos

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.helping](http://pymol.org/psicoredirect.php?psico.helping)  
  
**DESCRIPTION**

"apropos" searches through the documentation of all currently defined commands and lists those commands for which the keyword is either contained in the documentation or matches the command name itself. 

If an appropriate **DESCRIPTION** section is provided in the documentation of the command, the first 80 characters are listed as a summary. 

**INSTALLATION**

1\. copy the script bellow to a folder and name it _apropos.py_

2\. run it from within pymol with the **run** command: 
    
    
    run apropos.py
    

3\. see example bellow on how to use it 

**USAGE**
    
    
    apropos [keyword or regexp]
    

**EXAMPLE**

_apropos fit_
    
    
    ###EXACT MATCH FOR: fit ==> try 'help fit' at the prompt.
    
    ###The following commands are NOT documented.
    
          vdw_fit
    
    ###The following commands are documented.  'help command' 
    
              fit : "fit" superimposes the model in the first selection on to the model
        intra_fit : "intra_fit" fits all states of an object to an atom selection
              rms : "rms" computes a RMS fit between two atom selections, but does not
         pair_fit : "pair_fit" fits a set of atom pairs between two models.  Each atom
    intra_rms_cur : "intra_rms_cur" calculates rms values for all states of an object
         commands : >>>>>>>> Ooopsie, no DESCRIPTION found for this command!!! <<<<<<
             zoom : "zoom" scales and translates the window and the origin to cover the
        intra_rms : "intra_rms" calculates rms fit values for all states of an object
            align : "align" performs a sequence alignment followed by a structural
          rms_cur : "rms_cur" computes the RMS difference between two atom
          fitting : "fitting" allows the superpositioning of object1 onto object2 using
    

**SEE ALSO**
    
    
       grepset(www.pymolwiki.org), Python re module
    

**apropos.py**
    
    
    #
    # apropos.py 
    # Author: Ezequiel Panepucci
    # Date: 2006-07-20
    #
    from pymol import cmd
    import re
    
    def apropos(regexp=''):
        '''
    DESCRIPTION
            "apropos" searches through the documentation of all currently 
            defined commands and lists those commands for which the keyword
            is either contained in the documentation or matches the command
            name itself.
            
            If an appropriate "DESCRIPTION" section is provided in the documentation
            of the command, the first 80 characters are listed as a summary.
    
    USAGE
            apropos [keyword or regexp]
    
    EXAMPLE
            apropos fit
    
    ###EXACT MATCH FOR: fit ==> try 'help fit' at the prompt.
    
    ###The following commands are NOT documented.
    
          vdw_fit
    
    ###The following commands are documented.  'help command' 
    
              fit : "fit" superimposes the model in the first selection on to the model
        intra_fit : "intra_fit" fits all states of an object to an atom selection
              rms : "rms" computes a RMS fit between two atom selections, but does not
         pair_fit : "pair_fit" fits a set of atom pairs between two models.  Each atom
    intra_rms_cur : "intra_rms_cur" calculates rms values for all states of an object
         commands : >>>>>>>> Ooopsie, no DESCRIPTION found for this command!!! <<<<<<
             zoom : "zoom" scales and translates the window and the origin to cover the
        intra_rms : "intra_rms" calculates rms fit values for all states of an object
            align : "align" performs a sequence alignment followed by a structural
          rms_cur : "rms_cur" computes the RMS difference between two atom
          fitting : "fitting" allows the superpositioning of object1 onto object2 using
    
    SEE ALSO
        grepset(www.pymolwiki.org), Python re module
        '''
        cmd.set("text","1",quiet=1)
    
        count=0
        docre = re.compile(regexp, re.MULTILINE | re.IGNORECASE)
        cmdre = re.compile(regexp, re.IGNORECASE)
    
        matches_with_help = []
        matches_without_help = []
    
        maxcclen=0
        for cc in cmd.keyword:
            if cc == regexp:
                print '\n###EXACT MATCH FOR: %s ==> try \'help %s\' at the prompt.' % (cc,cc)
    
            doc = cmd.keyword[cc][0].__doc__
    
            if doc == None:
                if re.search(regexp, cc, re.IGNORECASE):
                    count += 1
                    matches_without_help.append(cc)
                continue
    
            if re.search(regexp, doc, re.MULTILINE | re.IGNORECASE):
                count += 1
                if len(cc) > maxcclen:
                    maxcclen = len(cc)
    
                docmatches = re.match(r"""^\s+DESCRIPTION\s+(.{0,80})\S*""", doc, re.IGNORECASE)
                if docmatches == None:
                    desc = '>>>>>>>> Ooopsie, no DESCRIPTION found for this command!!! <<<<<<'
                else:
                    desc = docmatches.group(1)
                matches_with_help.append( (cc, desc ) )
    
                
        if len(matches_without_help) > 0:
            print '\n###The following commands are NOT documented.\n'
            for cc in matches_without_help:
                print '%*s' % (maxcclen, cc)
    
        if len(matches_with_help) > 0:
            print '\n###The following commands are documented.  \'help command\' \n'
            for cc,desc in matches_with_help:
                print '%*s : %s' % (maxcclen, cc,desc)
    
    cmd.extend('apropos',apropos)
    

_Author: Ezequiel (Zac) Panepucci_

Retrieved from "[https://pymolwiki.org/index.php?title=Apropos&oldid=12885](https://pymolwiki.org/index.php?title=Apropos&oldid=12885)"


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

## Axes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Axes with text labels
    
    
    # axes.py
    from pymol.cgo import *
    from pymol import cmd
    from pymol.vfont import plain
    
    # create the axes object, draw axes with cylinders coloured red, green,
    #blue for X, Y and Z
    
    obj = [
       CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
       CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
       CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
       ]
    
    # add labels to axes object (requires pymol version 0.8 or greater, I
    # believe
    
    cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    
    # then we load it into PyMOL
    cmd.load_cgo(obj,'axes')
    

## Axes with nice cones

[![Axes demo.png](/images/e/e9/Axes_demo.png)](/index.php/File:Axes_demo.png)

This script draws a simple cartesian coordinate system. 
    
    
    from pymol.cgo import *
    from pymol import cmd
    
    w = 0.06 # cylinder width 
    l = 0.75 # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter
    
    obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
           CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
           CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    
    cmd.load_cgo(obj, 'axes')
    

## Axes which always stay in the lower left corner
    
    
    from pymol import cmd
    from chempy import cpv
    
    class PutCenterCallback(object):
        prev_v = None
    
        def __init__(self, name, corner=0):
            self.name = name
            self.corner = corner
            self.cb_name = cmd.get_unused_name('_cb')
    
        def load(self):
            cmd.load_callback(self, self.cb_name)
    
        def __call__(self):
            if self.name not in cmd.get_names('objects'):
                import threading
                threading.Thread(None, cmd.delete, args=(self.cb_name,)).start()
                return
    
            v = cmd.get_view()
            if v == self.prev_v:
                return
            self.prev_v = v
    
            t = v[12:15]
    
            if self.corner:
                vp = cmd.get_viewport()
                R_mc = [v[0:3], v[3:6], v[6:9]]
                off_c = [0.15 * v[11] * vp[0] / vp[1], 0.15 * v[11], 0.0]
                if self.corner in [2,3]:
                    off_c[0] *= -1
                if self.corner in [3,4]:
                    off_c[1] *= -1
                off_m = cpv.transform(R_mc, off_c)
                t = cpv.add(t, off_m)
    
            z = -v[11] / 30.0
            m = [z, 0, 0, 0, 0, z, 0, 0, 0, 0, z, 0, t[0] / z, t[1] / z, t[2] / z, 1]
            cmd.set_object_ttt(self.name, m)
    
    def axes(name='axes'):
        '''
    DESCRIPTION
    
        Puts coordinate axes to the lower left corner of the viewport.
        '''
        from pymol import cgo
    
        cmd.set('auto_zoom', 0)
    
        w = 0.06 # cylinder width
        l = 0.75 # cylinder length
        h = 0.25 # cone hight
        d = w * 1.618 # cone base diameter
    
        obj = [cgo.CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
               cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
               cgo.CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
               cgo.CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
               cgo.CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    
        PutCenterCallback(name, 1).load()
        cmd.load_cgo(obj, name)
    
    cmd.extend('axes', axes)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Axes&oldid=12542](https://pymolwiki.org/index.php?title=Axes&oldid=12542)"


---

## B2transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/b2transparency.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/b2transparency.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
b2transparency can set transparency settings like [transparency](/index.php/Transparency "Transparency") or [sphere_transparency](/index.php/Sphere_transparency "Sphere transparency") for each atom scaled by b-factor (or any other numeric atomic property). 

[![B2.png](/images/d/d1/B2.png)](/index.php/File:B2.png)

## Usage
    
    
    b2transparency [ selection [, setting [, minimum [, maximum [, var ]]]]]
    

## Example
    
    
    # some dummy molecule with increasing b-factor along the chain
    fab AAAAAAAAAAAAA
    alter all, b=resv
    
    b2transparency all
    
    as sticks
    show surface
    
    set surface_color, gray
    

## See Also

  * [spectrum](/index.php/Spectrum "Spectrum")



Retrieved from "[https://pymolwiki.org/index.php?title=B2transparency&oldid=13887](https://pymolwiki.org/index.php?title=B2transparency&oldid=13887)"


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

## Bounding Box

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    import math
    from pymol import querying
    from pymol.cgo import *
    from pymol import cmd
    
    #NOTE!! : These packages (numarray, Scientific) must be present in pymol's 
    #'$PYMOLDIR/py23/lib/python2.3/site-packages/' directory!!
    # OR..if you are using Mac PyMOL 0.99beta19 or later then the system installs will be used
    from numarray import *
    from numarray.ma import average
    from numarray import linear_algebra as la
    
    from Scientific.Geometry import Vector, Tensor, Transformation
    
    def boundingBox(selectionName, colourRGB=[1,1,1]):
            """
            The main function to call : 
    
                    run "box.py"
                    boundingBox("peptide")
    
            Should make a box around "peptide" (assuming it exists!). For a different colour use:
    
                    boundingBox("peptide", colourRGB=[1, 0, 1])
    
            Or whatever. The box should be a cgo called 'peptide-box' (for this example).
            """
            model = querying.get_model(selectionName)
            coords = model.get_coord_list()
    
            #find the least square plane through the coordinates
            eigenvalues, eigenvectors, centroid = findPlaneWithEigenVectors(coords)
            normal = eigenvectors[eigenvalues.argmin()]
            eigenvectors.remove(normal)
    
            #orient the axes correctly
            x, y, normal = orientAxes(eigenvectors[0], eigenvectors[1], normal)
    
            #determine the dimensions and the structure's orientation wrt the origin
            minDimensions, rotation = findMinDimensionsAndRotation(coords, centroid, x, y, normal)
    
            #'create' the box(IE: make the corners) and 'draw' it (IE: make a cgo)
            box = makeBox(minDimensions, rotation, centroid)
            drawBox(box, selectionName, colourRGB)
    
    def findPlaneWithEigenVectors(coords):
            centroid = average(coords)
            coords -= centroid
            B = matrixmultiply(transpose(coords), coords)
            eigenvalues, eigenvectors = la.eigenvectors(B)
            #return eigenvalues, [Vector(e) for e in eigenvectors], Vector(centroid)
            return eigenvalues, [Vector([i for i in e]) for e in eigenvectors], Vector(centroid) #not sure why I had to make this change!
    
    def orientAxes(x, y, z):
            XcrossY = x.cross(y)
            #ZXY = around(math.degrees(z.angle(XcrossY)))
            ZXY = int(around(math.degrees(z.angle(XcrossY)))) #again, a bit of a hack!
            if (ZXY == 180): x, y = y, x
            return x, y, z
    
    def makeBox(dimensions, rotation, centroid):
            x, y, z = dimensions
            v = [[]] * 8
    
            #make a cuboid with the lower corner on the origin
            v[0] = [0, 0, 0]                # [0, 0, 0]
            v[1] = [x, 0, 0]                # [1, 0, 0]
            v[2] = [x, y, 0]                # [1, 1, 0]
            v[3] = [x, 0, z]                # [1, 0, 1]
            v[4] = [0, y, 0]                # [0, 1, 0]
            v[5] = [0, 0, z]                # [0, 0, 1]
            v[6] = [0, y, z]                # [0, 1, 1]
            v[7] = [x, y, z]                # [1, 1, 1]
    
            #move to the origin, THEN move to the centroid of the points, then rotate
            translationToOrigin = Transformation.Translation(-Vector(x/2, y/2, z/2))
            translationToCentroid = Transformation.Translation(centroid)
            transform = translationToCentroid * rotation * translationToOrigin
    
            #use the Transformation to transform the corners of the box
            v = [transform(Vector(i)) for i in v]
    
            bot =  [v[0], v[1], v[2], v[4]] # O, x, xy, y
            top =  [v[7], v[3], v[5], v[6]] # xyz, xz, z, yz
            minL = [v[0], v[4], v[6], v[5]] # O, y, yz, z
            minR = [v[0], v[5], v[3], v[1]] # O, z, xz, x
            maxL = [v[4], v[2], v[7], v[6]] # y, xy, xyz, yz
            maxR = [v[3], v[1], v[2], v[7]] # xz, x, xy, xyz
            box = [bot, minR, minL, maxR, maxL, top]
    
            return box
    
    def drawBox(box, name, colourRGB):
            boxObj = []
            for side in box:
                    boxObj.append(BEGIN)
                    boxObj.append(LINE_STRIP)
                    boxObj.append(COLOR)
                    boxObj.extend(colourRGB)
                    for point in side:
                            boxObj.append(VERTEX)
                            boxObj.extend(point)
                    boxObj.append(END)
    
            cmd.set('auto_zoom', 0)
            cmd.load_cgo(boxObj, "%s-box" % name)
            cmd.set('auto_zoom', 1)
    
    def findMinDimensionsAndRotation(coords, centroid, x, y, z):
            O = Vector(0, 0, 0)
            X = Vector(1, 0, 0)
            Y = Vector(0, 1, 0)
            Z = Vector(0, 0, 1)
    
            #Create a Transformation t = |x, y, z| . |X, Y, Z| ^ -1
            mfinal = array([array(X), array(Y), array(Z)])
            morig = array([array(x), array(y), array(z)])
            rotmatrix = matrixmultiply(morig, transpose(mfinal))
            tTotal = Transformation.Rotation(Tensor(rotmatrix))
    
            #Transform the coordinates and find the min, max dimensions
            coordArray = array([array(tTotal(Vector(c))) for c in coords])
            minDimensions = [max(coordArray[:,i]) - min(coordArray[:,i]) for i in range(3)]
    
            #Now, compose the inverse rotation used to move the bounding box to the right place
            tForward = Transformation.Rotation(Tensor(matrixmultiply(mfinal, transpose(morig))))
    
            return minDimensions, tForward
    

Retrieved from "[https://pymolwiki.org/index.php?title=Bounding_Box&oldid=6395](https://pymolwiki.org/index.php?title=Bounding_Box&oldid=6395)"


---

## CalcArea

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This code will calculate the area of a given selection. If you have a bunch of disparate selections, or one selection made up of various objects then you can use this to find the area. 

# The Code
    
    
    def calcArea(sel, ASA=0, density=3, quiet=1):
    	"""
    DESCRIPTION
    
        Calculate the area of an object or selection (as if it were one complete object).
    
    ARGUMENTS
    
        sel = string: atom selection
    
        ASA = 0/1: solvent excluded or solvent accessible surface area
    
        density = int: sampling quality {default: 3}
    	"""
    	tName = "__temp__for__doArea__"
    	cmd.create(tName, sel)
    	cmd.flag("ignore", tName, "clear")
    	cmd.set("dot_solvent", int(ASA), tName)
    	cmd.set("dot_density", int(density), tName)
    	theArea = cmd.get_area(tName)
    	cmd.delete(tName)
    	if not int(quiet):
    		print 'Area (solvent %s): %.3f Angstroms^2' % (['excluded',
    			'accessible'][int(ASA)], theArea)
    	return theArea
    
    cmd.extend('calcArea', calcArea)
    

# See Also

  * [Get_Area](/index.php/Get_Area "Get Area")
  * [Set](/index.php/Set "Set")
  * [Flag](/index.php/Flag "Flag")
  * [Create](/index.php/Create "Create")
  * [Delete](/index.php/Delete "Delete")



Retrieved from "[https://pymolwiki.org/index.php?title=CalcArea&oldid=9125](https://pymolwiki.org/index.php?title=CalcArea&oldid=9125)"


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

## Cealign plugin

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Go directly to[DOWNLOAD](/index.php/Cealign#Version_0.8-RBS "Cealign")**

Note: CEAlign is now built into PyMOL as a native command. See the open-source project page. 

This page is the home page of the open-source CEAlign PyMOL plugin. The CE algorithm is a fast and accurate protein structure alignment algorithm, pioneered by Drs. Shindyalov and Bourne (See References). 

## Contents

  * 1 Introduction
  * 2 Comparison to PyMol
    * 2.1 Fit vs. optAlign
      * 2.1.1 Take Home messages
      * 2.1.2 Discussion
  * 3 Examples
    * 3.1 Usage
      * 3.1.1 Syntax
        * 3.1.1.1 Examples
        * 3.1.1.2 Multiple Structure Alignments
    * 3.2 Results
  * 4 Installation
    * 4.1 Mac OS X (10.5, 10.6)
    * 4.2 Windows systems
      * 4.2.1 CEAlign 0.9
        * 4.2.1.1 Requirements
        * 4.2.1.2 Directions
      * 4.2.2 CEAlign 0.8
        * 4.2.2.1 Requirements
        * 4.2.2.2 Directions
    * 4.3 Gentoo Linux
    * 4.4 *nix systems
      * 4.4.1 Requirements
      * 4.4.2 Directions
        * 4.4.2.1 Pre-compiled Hackish Install
  * 5 The Code
    * 5.1 Version 0.8-RBS
    * 5.2 Beta Version 0.9
  * 6 Coming Soon
  * 7 Updates
    * 7.1 2008-03-25
    * 7.2 2007-04-14
  * 8 Troubleshooting
    * 8.1 Unicode Issues in Python/Numpy
    * 8.2 LinAlg Module Not Found
    * 8.3 CCEAlign & NumPy Modules Not Found
    * 8.4 The Function SimpAlign() is not found
    * 8.5 Short Alignments Don't Work
    * 8.6 It Worked A Second Ago!
    * 8.7 file is not of required architecture
  * 9 References
  * 10 License



## Introduction

There are a few changes from the original CE publication (See Notes). The source code is implemented in C (and another in C++) with the rotations finally done by Numpy in Python (or C++ in version 0.9). Because the computationally complex portion of the code is written in C, it's quick. That is, on my machines --- relatively fast 64-bit machines --- I can align two 400+ amino acid structures in about 0.300 s with the C++ implementation. 

This plugs into PyMol very easily. See [the code](/index.php/Cealign#The_Code "Cealign") and [examples](/index.php/Cealign#Examples "Cealign") for installation and usage. 

## Comparison to PyMol

**Why should you use this?**

PyMOL's structure alignment algorithm is fast and robust. However, its first step is to perform a sequence alignment of the two selections. Thus, proteins in the **twilight zone** or those having a low sequence identity, may not align well. Because CE is a structure-based alignment, this is not a problem. Consider the following example. The two images below demonstrate the difference superimposing [1C0M](http://www.rcsb.org/pdb/explore/explore.do?structureId=1C0m) chain B onto [1BCO](http://www.rcsb.org/pdb/explore/explore.do?structureId=1BCO). The first image below shows the results from PyMol's `align` command: an alignment of **221 atoms** (not residues) to an RMSD of **15.7 Angstroms**. The second image is the result of CEAlign, which used alpha carbons of **152 residues** with an RMSD of **4.96 Angstroms**. 

  * [![PyMol's results \(221 atoms; 15.7 Ang. \)](/images/a/aa/Pymol_align.png)](/index.php/File:Pymol_align.png "PyMol's results \(221 atoms; 15.7 Ang. \)")

PyMol's results (221 atoms; 15.7 Ang. ) 

  * [![Cealign's results \(152 aligned; 4.96 Ang.\)](/images/d/d8/Cealign_ex1.png)](/index.php/File:Cealign_ex1.png "Cealign's results \(152 aligned; 4.96 Ang.\)")

Cealign's results (152 aligned; 4.96 Ang.) 




### Fit vs. optAlign

#### Take Home messages

  * [fit](/index.php/Fit "Fit") and [optAlign](/index.php/OptAlign "OptAlign") perform nearly equally as well
  * if you need an algorithm with an appropriate reference, use [optAlign](/index.php/OptAlign "OptAlign") (references at bottom of page).
  * [fit](/index.php/Fit "Fit") is faster -- if you're aligning many structures, use it over [optAlign](/index.php/OptAlign "OptAlign")



#### Discussion

[optAlign](/index.php/OptAlign "OptAlign") is a function within the [Cealign](/index.php/Cealign "Cealign") package that performs the optimal superposition of two objects of equal length. [optAlign](/index.php/OptAlign "OptAlign") follows the Kabsch algorithm which is a closed form, and provably optimal solution to the problem. [fit](/index.php/Fit "Fit") on the other hand uses the Jacobi rotations to iteratively arrive at the solution of optimal superposition. The difference in error between [optAilgn](/index.php?title=OptAilgn&action=edit&redlink=1 "OptAilgn \(page does not exist\)") and [fit](/index.php/Fit "Fit") seems to be a non-issue (see below) as they both arrive at equivalent solutions for the rotation matrix. The two algorithms are undertake different approaches to orthogonally diagonalizing the correlation matrix. 

PyMOL's [fit](/index.php/Fit "Fit") is fast and works well. If you have to use something with a known reference then check out the "optAlign" function from the qkabsch.py file that comes with this [Cealign](/index.php/Cealign "Cealign") package. If not, you can just use [fit](/index.php/Fit "Fit") and avoid installing new software. :-) 

optAlign is slower than fit. I just tested both on a sample NMR ensemble; and, while not an extensive validation of "fit" it shows that (1) fit is faster; and (2) fit gets the same exact RMSD as "optAlign" (when optAlign is told to use all atoms, not just CA). To make optAlign use all atoms and not just the alpha-carbon backbones, comment out (that is, put a "#" at the start of) lines 183 and 184 in qkabsch.py, where it says "CUT HERE." 
    
    
    fetch 1nmr
    split_states 1nmr
    delete 1nmr
    
    # compare fit and optAlign RMSDs
    for x in cmd.get_names(): print cmd.fit("1nmr_0001", x)
    for x in cmd.get_names(): optAlign(x, "1nmr_0001")
    
    
    
    # results from fit
    0.0
    4.50344991684
    5.33588504791
    5.78613853455
    7.25597000122
    6.67145586014
    3.25131297112
    3.36766290665
    6.74802017212
    5.1579709053
    5.96959495544
    6.68093347549
    4.13217163086
    5.51539039612
    6.24266338348
    6.03838825226
    5.01363992691
    5.33336305618
    6.87617444992
    7.797062397
    
    #results from optAlign
    RMSD=0.000000
    RMSD=4.503450
    RMSD=5.335886
    RMSD=5.786138
    RMSD=7.255970
    RMSD=6.671456
    RMSD=3.251313
    RMSD=3.367663
    RMSD=6.748021
    RMSD=5.157971
    RMSD=5.969595
    RMSD=6.680934
    RMSD=4.132172
    RMSD=5.515390
    RMSD=6.242664
    RMSD=6.038388
    RMSD=5.013640
    RMSD=5.333363
    RMSD=6.876174
    RMSD=7.797062
    

## Examples

### Usage

#### Syntax

CEAlign has the semantic, and syntactic formalism of 
    
    
    cealign MASTER, TARGET
    

where a post-condition of the algorithm is that the coordinates of the **MASTER** protein are unchanged. This allows for easier multi-protein alignments. For example, 
    
    
    cealign 1AUE, 1BZ4
    cealign 1AUE, 1B68
    cealign 1AUE, 1A7V
    cealign 1AUE, 1CPR
    

will superimpose all the TARGETS onto the MASTER. 

##### Examples
    
    
    cealign 1cll and i. 42-55, 1ggz and c. A
    cealign 1kao, 1ctq
    cealign 1fao, 1eaz
    

##### Multiple Structure Alignments

Use the **alignto** command, now provided with cealign. Just type, 
    
    
    alignto PROT
    

to align all your proteins in PyMOL to the one called, **PROT**. 

### Results

See **Changes** for updates. But, overall, the results here are great. 

  * Note: PyMOL v1.5.0.4 (svn revision 4001) has updates that improve some alignments slightly. These improved results are shown here.


  * [![EASY: 1FAO vs. 1EAZ; 96 residues, 1.09 Ang](/images/e/e2/V7_1fao_1eaz.png)](/index.php/File:V7_1fao_1eaz.png "EASY: 1FAO vs. 1EAZ; 96 residues, 1.09 Ang")

EASY: 1FAO vs. 1EAZ; 96 residues, 1.09 Ang 

  * [![EASY: 1CBS vs. 1HMT; 128 residues, 2.01 Ang](/images/2/25/V7_1cbs_1hmt.png)](/index.php/File:V7_1cbs_1hmt.png "EASY: 1CBS vs. 1HMT; 128 residues, 2.01 Ang")

EASY: 1CBS vs. 1HMT; 128 residues, 2.01 Ang 

  * [![MODERATE: 1A15 vs 1B50; 56 residues, 2.54 Ang.](/images/8/81/V7_1a15_1b50.png)](/index.php/File:V7_1a15_1b50.png "MODERATE: 1A15 vs 1B50; 56 residues, 2.54 Ang.")

MODERATE: 1A15 vs 1B50; 56 residues, 2.54 Ang. 

  * [![EASY: 1OAN vs. 1S6N \(state 1\); 96 residues aligned to 2.66 Ang. RMSD.](/images/5/58/V7_1oan_1s6n.png)](/index.php/File:V7_1oan_1s6n.png "EASY: 1OAN vs. 1S6N \(state 1\); 96 residues aligned to 2.66 Ang. RMSD.")

EASY: 1OAN vs. 1S6N (state 1); 96 residues aligned to 2.66 Ang. RMSD. 

  * [![HARD: 1RLW to 1BYN; 104 residues; 2.21 Ang.](/images/3/3a/V7_1rlw_1byn.png)](/index.php/File:V7_1rlw_1byn.png "HARD: 1RLW to 1BYN; 104 residues; 2.21 Ang.")

HARD: 1RLW to 1BYN; 104 residues; 2.21 Ang. 

  * [![HARD: 1TEN vs. 3HHR; 80 residues, 2.96 Ang.](/images/e/e3/V7_1ten_3hhr.png)](/index.php/File:V7_1ten_3hhr.png "HARD: 1TEN vs. 3HHR; 80 residues, 2.96 Ang.")

HARD: 1TEN vs. 3HHR; 80 residues, 2.96 Ang. 

  * [![HARD: 2SIM vs. 1NSB; 272 residues, 4.92 Ang.](/images/a/ad/V7_2sim_1nsb.png)](/index.php/File:V7_2sim_1nsb.png "HARD: 2SIM vs. 1NSB; 272 residues, 4.92 Ang.")

HARD: 2SIM vs. 1NSB; 272 residues, 4.92 Ang. 

  * [![HARD: 1CEW vs. 1MOL; 80 residues, 3.67 Ang.](/images/a/a6/V7_1cew_1mol.png)](/index.php/File:V7_1cew_1mol.png "HARD: 1CEW vs. 1MOL; 80 residues, 3.67 Ang.")

HARD: 1CEW vs. 1MOL; 80 residues, 3.67 Ang. 




## Installation

### Mac OS X (10.5, 10.6)

[![](/images/9/90/Cealign_mac_os_x.png)](/index.php/File:Cealign_mac_os_x.png)

[](/index.php/File:Cealign_mac_os_x.png "Enlarge")

CEAlign running on Mac OS X (10.5)

  * Install PyMOL under fink.
  * Download and install cealign (download instructions below)


    
    
    sudo /sw/bin/python setup.py install
    

  * In PyMOL, run the two scripts needed for cealign: "cealign.py" and "qkabsch.py". These are located in the cealign directory you previously downloaded.
  * Voila!
  * Note that the above python version must match the same version that is used by PyMOL. If you are using the pre-compiled version of MacPyMOL, the above instructions won't work.
  * Note: if you get an error about **-Wno-long-double** then your gcc is mismatched. I fixed this by pointing the symbolic link _/usr/bin/gcc_ from _/usr/bin/gcc-4.2_ to _/usr/bin/gcc-4.0_. Or, in code,


    
    
    # These command are commented out to stop people from copy/pasting b/c 
    # these are possibly dangerous for your system.  Ensure that /usr/bin/gcc
    # is a symbolic link and not a real binary.  If so, I used the following
    # to fix the -Wno-long-double error.
    # sudo rm /usr/bin/gcc
    # sudo ln -s /usr/bin/gcc-4.0 /usr/bin/gcc
    

### Windows systems

#### CEAlign 0.9

This is a Win32 build of CEAlign 0.9 [[1]](http://pymolwiki.org/index.php/Cealign#Beta_Version_0.9)

##### Requirements

  * Christoph Gohlke's latest **unofficial** PyMol build: <http://www.lfd.uci.edu/~gohlke/#pythonlibs>
  * "Python 2.6.2 Windows installer" from python.org: <http://www.python.org/download/>
  * **CEAlign09Win32.zip** from: <http://users.umassmed.edu/shivender.shandilya/pymol/CEAlign09Win32.zip>



##### Directions

  1. Download the **CEAlign09Win32.zip** file
  2. Unzip the downloaded file and follow the directions as per the included README.txt
  3. Enjoy the _awesomeness_ that is CEAlign!



#### CEAlign 0.8

This is a quick and dirty method to use CEAlign 0.8 on Win32 system with the **official** Pymol builds... 

##### Requirements

  * Latest PyMol, installed on your system
  * Numpy for python 2.4 -- quick download of just what's needed: <http://users.umassmed.edu/shivender.shandilya/pymol/cealign08/numpy.zip>



[Note: If this file is corrupt, you may download the latest 'Numpy for Python 2.4' directly from SourceForge.net 

  * Pre-compiled ccealign.pyd python module: <http://users.umassmed.edu/Shivender.Shandilya/pymol/cealign08/ccealign.zip>
  * Modified pymolrc: <http://users.umassmed.edu/Shivender.Shandilya/pymol/cealign08/pymolrc>
  * cealign.py and qkabsch.py from the Cealign-0.8-RBS package: download below



##### Directions

  1. Unzip the numpy.zip file, which will give you a folder named **numpy**
  2. Move this entire folder to: C:\Program Files\DeLano Scientific\PyMOL\modules\ (or the corresponding location on your system)
  3. Unzip ccealign.zip, which will give you a file called **ccealign.pyd**
  4. Move this pyd file to: C:\Program Files\DeLano Scientific\PyMOL\py24\DLLs\ (or the corresponding location on your system)
  5. Copy the downloaded **pymolrc** file to: C:\Program Files\DeLano Scientific\PyMOL\ (or the corresponding location on your system)
  6. Extract and copy the files cealign.py and qkabsch.py from the Cealign-0.8-RBS package to: C:\Program Files\DeLano Scientific\PyMOL\py24\Lib\ (or the corresponding location on your system)
  7. Run PyMol and load some molecules
  8. Run this command in Pymol: **cealign molecule1, molecule2**
  9. Enjoy!



### Gentoo Linux

Add the science overlay via 
    
    
    layman -a sci
    

and emerge the cealign plugin 
    
    
    emerge pymol-plugins-cealign
    

### *nix systems

#### Requirements

  * C compiler
  * Python 2.4+ with distutils
  * Numpy 
    * for User-compiled PyMOL: 
          
          python setup.py install
          

    * for the precompiled version of PyMOL 
          
          python setup.py install --prefix "" --root /DIR_TO/pymol/ext/
          




#### Directions

  1. uncompress the distribution file **cealign-VERSION.tgz**
  2. cd cealign-VERSION
  3. sudo python setup.py install # if you installed by PyMOL by hand 
     1. python setup.py install --prefix "" --root /DIR/TO/pymol/ext/ # if you are using the precompiled binary download
  4. insert "run DIR_TO_CEALIGN/cealign.py" and "run DIR_TO_CEALIGN/qkabsch.py" into your **.pymolrc** file, or just run the two Python scripts by hand.
  5. load some molecules
  6. run, **cealign molecule1, molecule2**
  7. enjoy



##### Pre-compiled Hackish Install

For those people that prefer to use the pre-compiled version of PyMOL, here are the basics for your install. **This is a poor method of installing Cealign. I suggest users compile and install their own PyMOL.** The final goal is to get 

  1. **ccealign.so** module into **PYMOL/ext/lib/python2.4/site-packages**
  2. numpy installed (get the numpy directory into (or linked into) **PYMOL/ext/lib/python2.4/site-packages**
  3. and be able to run cealign.py and qkabsch.py from PyMOL.



If you can do the above three steps, **cealign** should run from the pre-compiled PyMOL. 

In more detail, on a completely fictitious machine --- that is, I created the following commands from a fake machine and I don't expect a copy/paste of this to work **anywhere** , but the commands should be helpful enough to those who need it: 
    
    
    # NOTES:
    # This is fake code: don't copy/paste it.
    #
    # PYMOL='dir to precompiled PyMOL install'
    # CEALIGN='dir where you will unpack cealign'
    # replace lib with lib64 for x86-64
    # install numpy
    apt-get install numpy
    
    # link numpy to PyMOL
    ln -s /usr/local/lib/python2.4/site-packages/numpy PYMOL/ext/lib/python2.4/site-packages
    
    # download and install Cealign
    wget http://www.pymolwiki.org/images/e/ed/Cealign-0.6.tar.bz2
    tar -jxvf Cealign-0.6.tar.bz2
    cd cealign-0.6
    sudo python setup.py build
    cp build/lib-XYZ-linux/ccealign.so PYMOL/ext/lib/python2.4/site-packages
    
    # run pymol and try it out
    pymol
    run CEALIGN/cealign.py
    run CEALIGN/qkabsch.py
    fetch 1cew 1mol, async=0
    cealign 1c, 1m
    

## The Code

Please unpack and read the documentation. All comments/questions should be directed to Jason Vertrees (javertre _at_ utmb ...dot... edu). 

**LATEST IS v0.8-RBS**. (Dedicated to Bryan Sutton for allowing me to use his computer for testing.) 

### Version 0.8-RBS

  * **Download:[CE Align v0.8-RBS](/images/5/58/Cealign-0.8-RBS.tar.bz2 "Cealign-0.8-RBS.tar.bz2") (bz2)**
  * **Download:[CE Align v0.8-RBS](/images/5/59/Cealign-0.8-RBS.zip "Cealign-0.8-RBS.zip") (zip)**



### Beta Version 0.9

Use at your own peril. Please report any problems or inconsistent alignments to this discussion page, or to me directly; my email address all over this page. 

**Improvements/Changes** : 

  * All C++ 
    * So, faster
    * comes with the dependencies built in
  * No numpy



**Download:[CE Align v0.9](/images/0/03/Cealign-0.9.zip "Cealign-0.9.zip") (zip)**

## Coming Soon

  * Windows binary
  * Linux Binaries (32bit, x86-64)
  * Better instructions for precompiled distributions
  * Optimization



## Updates

### 2008-03-25

Pure C++ code released. See the beta version above. 

### 2007-04-14

v0.8-RBS source updated. Found the bug that had been plaguing 32-bit machines. This should be the last release for a little while. 

Also, I provide the option of aligning based solely upon RMSD or upon the better CE-Score. See the **References** for information on the **CE Score**. 

## Troubleshooting

Post your problems/solutions here. 

### Unicode Issues in Python/Numpy

**Problem** : Running/Installing cealign gives 
    
    
    Traceback (most recent call last):
      File "/home/byron/software/pymol_1.00b17/pymol/modules/pymol/parser.py",
    line 308, in parse
      File "/home/byron/software/pymol_1.00b17/pymol/modules/pymol/parsing.py",
    line 410, in run_file
      File "qkabsch.py", line 86, in ?
        import numpy
      File "/usr/lib/python2.4/site-packages/numpy/__init__.py", line 36, in ?
        import core
      File "/usr/lib/python2.4/site-packages/numpy/core/__init__.py", line 5, in ?
        import multiarray
    ImportError: /home/byron/software/pymol/ext/lib/python2.4/site-packages/numpy/core/multiarray.so:
    undefined symbol: _PyUnicodeUCS4_IsWhitespace
    

where the important line is 
    
    
    undefined symbol: _PyUnicodeUCS4_IsWhitespace
    

This problem indicates that your Numpy Unicode is using a different byte-size for unicode characters than is the Python distribution your PyMOL is running from. For example, this can happen if you use the pre-built PyMOL and some other pre-built Numpy package. 

  


**Solution** : Hand-install Numpy. 

  


### LinAlg Module Not Found

**Problem** : Running CE Align gives the following error message: 
    
    
    run qkabsch.py
    Traceback (most recent call last):
    File "/usr/lib/python2.4/site-packages/pymol/parser.py", line 285, in parse
    parsing.run_file(exp_path(args[nest][0]),pymol_names,pymol_names)
    File "/usr/lib/python2.4/site-packages/pymol/parsing.py", line 407, in run_file
    execfile(file,global_ns,local_ns)
    File "qkabsch.py", line 86, in ?
    import numpy
    File "/usr/lib/python2.4/site-packages/numpy/__init__.py", line 40, in ?
    import linalg
    ImportError: No module named linalg
    

  


**Solution** : You do not have the linear algebra module installed (or Python can't find it) on your machine. One workaround is to install [Scientific Python](http://www.scipy.org/). (on debian/ubuntu this can be done by: sudo apt-get install python-scipy) Another is to reinstall the Numpy package from source, ensuring that you have the necessary requirements for the linear algebra module (linpack, lapack, fft, etc.). 

### CCEAlign & NumPy Modules Not Found

**Problem** : Running CE Align gives the following error message: 
    
    
    PyMOL>run cealign.py
    Traceback (most recent call last):
      File "/home/local/warren/MacPyMOL060530/build/Deployment/MacPyMOL.app/pymol/modules/pymol/parser.py", line 297, in parse
      File "/home/local/warren/MacPyMOL060530/build/Deployment/MacPyMOL.app/pymol/modules/pymol/parsing.py", line 408, in run_file
      File "/usr/local/pymol/scripts/cealign-0.1/cealign.py", line 59, in ?
        from ccealign import ccealign
    ImportError: No module named ccealign
    run qkabsch.py
    Traceback (most recent call last):
    File "/home/local/warren/MacPyMOL060530/build/Deployment/MacPyMOL.app/pymol/modules/pymol/parser.py", line 297, in parse
    File "/home/local/warren/MacPyMOL060530/build/Deployment/MacPyMOL.app/pymol/modules/pymol/parsing.py", line 408, in run_file
    File "qkabsch.py", line 86, in ?
    import numpy
    ImportError: No module named numpy
    

  


**Solution** : This problem occurs under [Apple Mac OS X](http://www.apple.com/macosx) if (a) the Apple's python executable on your machine (/usr/bin/python, currently version 2.3.5) is superseded by [Fink](http://fink.sourceforge.net/)'s python executable (/sw/bin/python, currently version 2.5) and (b) you are using [precompiled versions of PyMOL](http://delsci.com/rel/099/#MacOSX) (MacPyMOL, PyMOLX11Hybrid or PyMOL for Mac OS X/X11). These executables ignore Fink's python and instead use Apple's - so, in order to run CE Align, one must install NumPy (as well as CE Align itself) using Apple's python. To do so, first download the [Numpy source code archive](http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103) (currently version 1.0.1), unpack it, change directory to numpy-1.0.1 and specify the full path to Apple's python executable during installation: `sudo /usr/bin/python setup.py install | tee install.log`. Then, donwload the [CE Align source code archive](http://www.pymolwiki.org/index.php/Cealign#The_Code) (currently version 0.2), unpack it, change directory to cealign-0.2 and finally install CE Align as follows: `sudo /usr/bin/python setup.py install | tee install.log`. [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)") 05:11, 25 January 2007 (CST). 

### The Function SimpAlign() is not found

**Problem** : Running CE Align gives the following error message: 
    
    
    PyMOL>cealign 1CLL,1GGZ
    Traceback (most recent call last):
      File "C:\Program Files (x86)\DeLano Scientific\PyMOL/modules\pymol\parser.py", line 203, in parse
        result=apply(kw[nest][0],args[nest],kw_args[nest])
      File "py24/Lib/cealign.py", line 177, in cealign
        curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )
    NameError: global name 'simpAlign' is not defined
    

I am running PyMOL v. 0.99rc6 on Win XP Professional x64 edition version 2003 sp2 and have followed the windows install procedure as described above. 

**Answer** : This simply means that PyMOL couldn't find the simplAlign function. To let PyMOL know about this, you must run the following commands before running [cealign](/index.php/Cealign "Cealign"): 
    
    
    run /your/path/to/cealign/qkabsch.py
    run /your/path/to/cealign/cealign.py
    

but most people that use cealign would just put these two lines in their **.pymolrc** file. 

### Short Alignments Don't Work

If you are trying to align fewer than 16 residues then use [align](/index.php/Align "Align"), [super](/index.php/Super "Super"), or [optAlign](/index.php/OptAlign "OptAlign"). CE uses a window size of 8; and to build a path of more than one window, you need 2*8=16 residues. I will insert some code to re-route small alignments to one of the aforementioned alignment algorithms. 

### It Worked A Second Ago!

[![](/images/0/01/Rewind.png)](/index.php/File:Rewind.png)

[](/index.php/File:Rewind.png "Enlarge")

Showing the rewind button to rewind to state 1.

If you were using cealign (or alignto) and now the commands don't work -- that is, they return an RMSD, but don't actually superimpose the objects, then you have a simple problem dealing with states. Most likely the cause of this oddness was (1) when you issued "cealign prot1, prot2" one of them was actually an ensemble of states or (2) you are trying to align to proteins with only one state, but are not looking at state one (because the last protein you were considering had more than one state and you quit editing that protein on a state that's not state 1). To fix this, use the rewind button to get the proteins back into state 1 & reissue the cealign/alignto command. 

### file is not of required architecture

This error happens on a Mac when you compile one bit of code with gcc-4.0/g++-4.0 and then try to make a library with code compiled from gcc-4.2/g++-4.2. If you recent installed Snow Leopard (Mac OS X 10.6) then this might bother you when you try to install Cealign or even PyMOL. To get around this, ensure that you're building all components with the same gcc/g++ executable. Here's how I did it, 
    
    
    # sudo rm /usr/bin/gcc /usr/bin/g++
    # sudo ln -s /usr/bin/gcc-4.0 /usr/bin/gcc
    # sudo ln -s /usr/bin/g++-4.0 /usr/bin/g++
    

I commented out those lines to stop people from blindly copy/pasting possible harmful lines. Please ensure that your /usr/bin/gcc and /usr/bin/g++ are actually symbolic links, otherwise you could be doing bad things to your computer. In my case, I only relinked gcc and not g++, hence the error. 

## References

Text taken from PubMed and formatted for the wiki. The first reference is the most important for this code. 

  1. Shindyalov IN, Bourne PE. **Protein structure alignment by incremental combinatorial extension (CE) of the optimal path.** _Protein Eng._ 1998 Sep;11(9):739-47. PMID: 9796821 [PubMed - indexed for MEDLINE]
  2. Jia Y, Dewey TG, Shindyalov IN, Bourne PE. **A new scoring function and associated statistical significance for structure alignment by CE.** _J Comput Biol._ 2004;11(5):787-99. PMID: 15700402 [PubMed - indexed for MEDLINE]
  3. Pekurovsky D, Shindyalov IN, Bourne PE. **A case study of high-throughput biological data processing on parallel platforms.** _Bioinformatics._ 2004 Aug 12;20(12):1940-7. Epub 2004 Mar 25. PMID: 15044237 [PubMed - indexed for MEDLINE]
  4. Shindyalov IN, Bourne PE. **An alternative view of protein fold space.** _Proteins._ 2000 Feb 15;38(3):247-60. PMID: 10713986 [PubMed - indexed for MEDLINE]



## License

The CEAlign and all its subprograms that I wrote, are released under the open source Free BSD License (BSDL). 

Retrieved from "[https://pymolwiki.org/index.php?title=Cealign_plugin&oldid=11715](https://pymolwiki.org/index.php?title=Cealign_plugin&oldid=11715)"


---

## Center of mass

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/center_of_mass.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/center_of_mass.py)  
Author(s)  | [Henschel](/index.php?title=User:Henschel&action=edit&redlink=1 "User:Henschel \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Examples
  * 4 See Also



### Description

This script calculates the true center-of-mass (COM) or the center-of-geometry (COG) for a given selection and returns the x, y, z values in the form of a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") (rather than a CGO sphere). The benefit of using a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") is that it can be selected and used in calculations. In addition, this script also iteratively goes through all states of a selection if more than one state exists and appends the corresponding COM/COG values as states into the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). The script itself is quite simple yet robust enough to be applied in different settings. As well, the calculation of the COM/COG is handled independently from the formation of the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") and can be called as an independent function where applicable. 

### Usage
    
    
    com selection [,state=None [,mass=None [,object=None]]]
    

  


### Examples
    
    
    import center_of_mass
    fetch 1c3y, finish=1, multiplex=0
    
    com 1c3y, state=1
    #Create a pseudoatom representing the 1c3y COG and store it as "1c3y_COM"
    #The "1c3y_COM" object will contain 1 state only
    
    com 1c3y, state=1, object=COG
    #Create a pseudoatom representing the 1c3y COG and store it as "COG"
    #The "COG" object will contain 1 state only
    
    com 1c3y, state=1, object=COM, mass=1
    #Create a pseudoatom representing the 1c3y COM and store it as "COM"
    #The "COM" object will contain 1 state only
    
    com 1c3y, object=COM, mass=1
    #Create a single pseudoatom containing the COM for each state found in 1c3y and store it as "COM"
    #The "COM" object will contain MULTIPLE states!
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [COM](/index.php/COM "COM")
  * [get_extent](/index.php/Get_extent "Get extent")
  * [get_position](/index.php/Get_position "Get position")



Retrieved from "[https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892](https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892)"


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

## Cgo arrow

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/cgo_arrow.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_arrow.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO arrow between two picked atoms. 

## Usage
    
    
    cgo_arrow [ atom1 [, atom2 [, radius [, gap [, hlength [, hradius [, color [, name ]]]]]]]]
    

## Example
    
    
    run cgo_arrow.py
    
    fetch 1rx1, async=0
    preset.pretty("*")
    
    cgo_arrow [34.9, 68.4, 19.1], A/164/O3X, gap=1.0
    

[![Cgo arrow example.png](/images/1/17/Cgo_arrow_example.png)](/index.php/File:Cgo_arrow_example.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_arrow&oldid=13894](https://pymolwiki.org/index.php?title=Cgo_arrow&oldid=13894)"


---

## Cgo grid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/cgo_grid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_grid.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![cgo_grid creates flowing mesh objects](/images/2/27/Cgo_grid.gif)](/index.php/File:Cgo_grid.gif "cgo_grid creates flowing mesh objects")

## Contents

  * 1 About cgo_grid
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 SEE ALSO



## About cgo_grid

**cgo_grid** will generate a flowing mesh object using the points provided or the current view. By default is will generate a flowing membrane. The shape is affected substantially by the arguments!  


  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage

**cgo_grid** has many arguments, but not all need to be set necessarily (see arguments or examples). 
    
    
    cgo_grid [ pos1 [, pos2 [, pos3 [, length_x [, length_z [, npoints_x [, npoints_z
    [, nwaves_x [, nwaves_z [, offset_x [, offset_z [, gain_x [, gain_z
    [, thickness [, color [, nstates [, startframe [, endframe
    [, mode [, view [, name [, quiet ]]]]]]]]]]]]]]]]]]]]]]
    

|   
---|---  
  
  


## Arguments
    
    
        pos1 = single atom selection (='pk1') or list of 3 floats {default: [0,0,0]}
    
        pos2 = single atom selection (='pk2') or list of 3 floats {default: [1,0,0]}
    
        pos3 = single atom selection (='pk3') or list of 3 floats {default: [0,0,1]}
    
        --> the plane is defined by pos1 (origin) and vectors to pos2 and pos3, respectively
    
        length_x = <float>: length of membrane {default: 30}
        length_z = <float>: length of membrane {default: ''} # same as length_x
    
        npoints_x = <int>: number of points(lines) along x-direction
                    {default: ''} #will be set to give a ~1 unit grid
        npoints_z = <int>: number of points(lines) along z-direction
                    {default: ''} #will be set to give a ~1 unit grid
                    {minimum: 1 # automatic}
    
        nwaves_x =   <float>: number of complete sin waves along object x-axis
                     {default: 2}
        nwaves_z =  <float>: number of complete sin waves along object z-axis
                    {default: ''} # same as nwaves_x
                    define separately to adjust number of waves in each direction
    
        offset_x = <float> phase delay (in degrees) of sin wave in x-axis
                 can be set to affect shape and starting amplitude {default: 0}
        offset_z = <float> phase delay (in degrees) of sin wave in z-axis
                 can be set to affect shape and starting amplitude
                 {default: ''} # same as  offset_x
        offset_x and offset_z can be used together to phase
        otherwise identical objects
    
        gain_x = <float>: multiplication factor for y-amplitude for x-direction
                 {default: 1}
        gain_z = <float>: multiplication factor for y-amplitude for z-direction
                 {default: ''} #=gain_x
    
        thickness = <float>: line thickness {default: 2}
    
        color = color name <string> (e.g. 'skyblue') OR
                rgb-value list of 3 floats (e.g. [1.0,1.0,1.0]) OR
                {default: ''} // opposite of background
                input illegal values for random coloring
    
        nstates =  <int>: number of states; {default: 60}
                   this setting will define how many states
                   the object will have (per wave) and how fluent and fast the
                   animation will be.
                   Higher values will promote 'fluent' transitions,
                   but decrease flow speed.
                       Note: Frame animation cycles thought the states one at a time
                       and needs to be set accordingly. Can also be used to phase
                       otherwise identical objects.
                   Set to 1 for static object {automatic minimum}
    
        startframe: specify starting frame <int> or set (='') to use current frame
                    set to 'append' to extend movie from the last frame {default: 1}
          endframe: specify end frame <int> or set (='') to use last frame
                    if 'append' is used for startframe,
                    endframe becomes the number of frames to be appended instead
                    {default: 1}
                    Note: if start- and endframe are the same, movie animation will
                    be skipped, the object will be loaded and can be used afterwards
    
        mode: defines positioning {default: 0}:
        0: pos1 is center
        1: pos1 is corner
    
        view {default: 0}:
        '0': off/ uses provided points to create CGO
        '1': overrides atom selections and uses current orienatation for positioning
             - pos1 = origin/center
             - pos2 = origin +1 in camera y
             - pos3 = origin +1 in camera z
    
        name: <string> name of cgo object {default: ''} / automatic
    
        quiet: <boolean> toggles output
    

| Concept sketch:  
[![cgo_grid concept sketch](/images/d/d1/Cgo_grid.png)](/index.php/File:Cgo_grid.png "cgo_grid concept sketch")  
---|---  
  
  


## Examples

The behaviour or shape of the cgo_grid object are substantially influenced by the arguments 
    
    
    delete all
    set_view (\
         0.263772517,   -0.113038681,    0.957937598,\
        -0.040853567,    0.990910411,    0.128179103,\
        -0.963716805,   -0.072944686,    0.256756991,\
         0.000000000,    0.000000000, -131.816467285,\
         0.000000000,    0.000000000,    0.000000000,\
       -50.008331299,  353.641235352,  -20.000000000 )
    
    #membrane
    cgo_grid color=blue
    
    #swimming worm, random color
    cgo_grid \
    pos1=[0,-5,0], pos2=[1,-5,1], pos3=[0,-5,1],\
    length_x=15,\
    npoints_z=1,\
    gain_x=2,\
    gain_z=0,\
    thickness=20,\
    color=3,\
    mode=1,\
    name="worm"
    
    #Moving Ladder
    cgo_grid \
    length_x=15,\
    pos1=[0,10,0], pos2=[0,10,1], pos3=[0,9,0],\
    npoints_x=2, npoints_z=30,\
    name="ladder"
    
    #Roof
    cgo_grid \
    nstates=1,\
    npoints_x=15,\
    npoints_z=15,\
    gain_x=20,\
    gain_z=20,\
    nwaves_x=0.5,\
    thickness=5,\
    color=cyan,\
    name="roof"
    
    #Boxes
    cgo_grid \
    pos1=[0,-10,0], pos2=[1,-10,0], pos3=[0,-10,1],\
    nstates=1,\
    npoints_x=50,\
    npoints_z=50,\
    nwaves_x=0,\
    color=[0.00 , 0.53 , 0.22],\
    thickness=5,\
    name="bottom"
    
    cgo_grid \
    nstates=1,\
    npoints_x=2,\
    npoints_z=2,\
    nwaves_x=0,\
    color=gray60,\
    thickness=10,\
    name="top"
    
    cgo_grid \
    pos1=[-15,-10,15], pos2=[-14,-10,15], pos3=[-15,-9,15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=-2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="front"
    
    cgo_grid \
    pos1=[-15,-10,-15], pos2=[-14,-10,-15], pos3=[-15,-9,-15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="back"
    
    set ray_trace_frames, 0
    set movie_loop,1
    mplay
    
    # play around with the ARGUMENTS! :-)
    

|   
---|---  
  
  


## SEE ALSO

  * [CgoCircle](/index.php/CgoCircle "CgoCircle")



Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938](https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938)"


---

## CGO Text

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/5/5f/Cgo_text.png)](/index.php/File:Cgo_text.png)

[](/index.php/File:Cgo_text.png "Enlarge")

CGO Text Example from following example
    
    
    # draw text using cgo
    from pymol import cmd
    from pymol.cgo import *
    from pymol.vfont import plain
    
    cgo = []
    
    axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]
    
    pos = [0.0,0.0,0.0]
    wire_text(cgo,plain,pos,'Hello World',axes)
    
    pos = [0.0,-3.0,0.0]
    cyl_text(cgo,plain,pos,'Hello Universe',0.10,axes=axes)
    
    cmd.set("cgo_line_radius",0.03)
    cmd.load_cgo(cgo,'txt')
    cmd.zoom("all",2.0)
    

Retrieved from "[https://pymolwiki.org/index.php?title=CGO_Text&oldid=10878](https://pymolwiki.org/index.php?title=CGO_Text&oldid=10878)"


---

## CgoCircle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will create a CGO circle with the origin at the specified X,Y,Z coordinates. Also, you can specify the radius and the colors. See the examples. 

If you want to draw a circle around an object or selection, use **circleSelection**. If you want pure flexibility over your circle then use **cgoCircle**. 

  
**There are two functions here:**

**cgoCircle**

    

    — creates a CGO circle at some user-specified location

**circleSelection**

    

    —creates a circle around the named object or selection.

  


  * [![Drawn circle.](/images/5/52/Circle1.png)](/index.php/File:Circle1.png "Drawn circle.")

Drawn circle. 

  * [![CGO circle.](/images/8/87/Circle2.png)](/index.php/File:Circle2.png "CGO circle.")

CGO circle. 

  * [![Circles of specified radius.](/images/4/47/Circle3.png)](/index.php/File:Circle3.png "Circles of specified radius.")

Circles of specified radius. 

  * [![Circle with specified width.](/images/d/da/CircleR.png)](/index.php/File:CircleR.png "Circle with specified width.")

Circle with specified width. 

  * [![Circle with a line width of 150. Pores anyone?](/images/d/d3/CircleR2.png)](/index.php/File:CircleR2.png "Circle with a line width of 150. Pores anyone?")

Circle with a line width of 150. Pores anyone? 




# Usage
    
    
    import math
    import pymol
    from pymol.cgo import *
    import random
    
    def cgoCircle(x, y, z, r=8.0, cr=1.0, cg=0.4, cb=0.8, w=2.0):
      """
      Create a CGO circle
    
      PARAMS
            x, y, z
              X, Y and Z coordinates of the origin
    
            r
              Radius of the circle
    
            cr, cg, cb
              Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].
    
            w
              Line width of the circle
    
      RETURNS
            the CGO object (it also loads it into PyMOL, too).
    
      """
      x = float(x)
      y = float(y)
      z = float(z)
      r = abs(float(r))
      cr = abs(float(cr))
      cg = abs(float(cg))
      cb = abs(float(cb))
      w = float(w)
    
      obj = [ BEGIN, LINES, COLOR, cr, cg, cb ]
      for i in range(180):
            obj.append( VERTEX )
            obj.append(r*math.cos(i) + x )
            obj.append(r*math.sin(i) + y )
            obj.append(z)
            obj.append( VERTEX )
            obj.append(r*math.cos(i+0.1) + x )
            obj.append(r*math.sin(i+0.1) + y )
            obj.append(z)
      obj.append(END)
     
      cName = cmd.get_unused_name("circle_")
      cmd.load_cgo( obj, cName )
      cmd.set("cgo_line_width", w, cName )
      return obj
    
    
    def circleSelection( selName, r=None, cr=1.0, cg=0.4, cb=0.8, w=2.0 ):
      """
      circleSelection -- draws a cgo circle around a given selection or object
    
      PARAMS
            selName
              Name of the thing to encircle.
    
            r
              Radius of circle.
              DEFAULT: This cript automatically defines the radius for you.  If
              you select one atom and the resultant circle is too small, then
              you can override the script's calculation of r and specify your own.
    
            cr, cg, cb
              red, green and blue coloring, each a value in the range [0.0, 1.0]
    
      RETURNS
            The circle object.
    
      """
      ((minX, minY, minZ), (maxX, maxY, maxZ)) = cmd.get_extent(selName)
    
      if r==None:
            r = max( [maxX-minX, maxY-minY, maxZ-minZ] )
    
      stored.coords = []
      cmd.iterate_state(1, selName, "stored.coords.append([x,y,z])")
      l = len(stored.coords)
    
      centerX = sum(map(lambda x: x[0], stored.coords)) / l
      centerY = sum(map(lambda x: x[1], stored.coords)) / l
      centerZ = sum(map(lambda x: x[2], stored.coords)) / l
    
      return cgoCircle( centerX, centerY, centerZ, r, cr, cg, cb, w )
    
    
    cmd.extend( "cgoCircle", cgoCircle )
    cmd.extend( "circleSelection", circleSelection )
    

# Updates

  * Line width option
  * better circle naming



Retrieved from "[https://pymolwiki.org/index.php?title=CgoCircle&oldid=10881](https://pymolwiki.org/index.php?title=CgoCircle&oldid=10881)"


---

## Check Key

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [check_key.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/check_key.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 DESCRIPTION
  * 2 USAGE
  * 3 PyMOL API
  * 4 SEE ALSO



### DESCRIPTION

A simple script used to check if a given key is valid for for the [Set Key](/index.php/Set_Key "Set Key") command. This is useful when the user would like to use a keyboard key as shortcut/hotkey in their script but need to check if the key is a valid one. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    check_key (keystroke)
    

If the key specified is a valid one as defined in [Set Key](/index.php/Set_Key "Set Key") then it returns the value of the keystroke. Otherwise, it returns None. 

For an example that calls this script, see [Jump](/index.php/Jump "Jump"). 

### PyMOL API
    
    
    from pymol import cmd
    import re
    
    allowed_keys=re.compile('(F1|F2|left|right|pgup|pgdn|home|insert|(CTRL-[A-Z])|ALT-[A-Z0-9])')
    
    def check_key (keystroke):
    
        """
        Author Sean M. Law
        University of Michigan
        seanlaw_(at)_umich_dot_edu
        """
    
        keystroke=keystroke.strip('\"\'')
        out=allowed_keys.search(keystroke)
    
        if (out != None):
            return keystroke         
        else:
            print "Error: Invalid key \""+keystroke+"\""
            print "Valid Keys: F1, F2, left, right, pdup, pgdn, home, insert"
            print "            CTRL-A to CTRL-Z"
            print "            ALT-0 to ALT-9, ALT-A to ALT-Z"
            return
    cmd.extend("check_key", check_key)
    

### SEE ALSO

[Button](/index.php/Button "Button") [Set Key](/index.php/Set_Key "Set Key")

Retrieved from "[https://pymolwiki.org/index.php?title=Check_Key&oldid=11147](https://pymolwiki.org/index.php?title=Check_Key&oldid=11147)"


---

## Cluster Count

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 Example
  * 4 The Code



# Overview

This script calculates statistics on the B-values for all atoms in the selected object, prints the information on screen and appends it to a file called "cluster_count.txt". 

  


# Usage
    
    
    cluster_count object
    

  


# Example
    
    
    cluster_count 1ubq
    
    #output on screen:
    
    Number of atoms in ' 1ubq ':  602
    Minimum and Maximum B-values:  2.0 42.75
    Average B-value:  13.4131063029
    Standard deviation of the B-values:  8.70767140923
    This data will be appended to cluster_count.txt
    
    #output in file cluster_count.txt; the format is:
    #objectname N minB maxB aveB stdevB
    1ubq            602    2.000   42.750   13.413    8.708
    

  


# The Code
    
    
    # Script: cluster_count.py
    # Copyleft 2010 Martin Christen
    
    from pymol import cmd,stored
    def cluster_count(selection):
    
            """
     AUTHOR
     
     Martin Christen
     
     DESCRIPTION
     
     This script calculates statistics on the B-values for all atoms in
     the selected object, prints the information on screen and appends
     it to the file "cluster_count.txt".
     
     Output format on screen:
     ------------------------
     Number of atoms in 'selection': 0
     Minimum and Maximum B-values:  0.0
     Average B-value :  0.0
     Standard deviation of the B-values (best): 0.0 (0.0)
     This data will be appended to cluster_count.txt
    
     Output format in cluster_count.txt:
     -----------------------------------
     selection N minB maxB aveB stdevB
    
     EXAMPLE
    
     cluster_count 1ubq
    
     Number of atoms in ' 1ubq ':  602
     Minimum and Maximum B-values:  2.0 42.75
     Average B-value:  13.4131063029
     Standard deviation of the B-values:  8.70767140923
     This data will be appended to cluster_count.txt
    
     (in cluster_count.txt:)
     1ubq            602    2.000   42.750   13.413    8.708
    
     USAGE
    
     cluster_count selection
    
            """
    	
    # get list of B-factors from selection
    	m = cmd.get_model(selection)
    	sel = []
    	b_list = []
    	dummy = []
    	for i in range(len(m.atom)):
    		b_list.append(m.atom[i].b)
    
    #determine min and max
    	try: max_b = max(b_list)
    	except ValueError: max_b=0
    	try: min_b = min(b_list)
    	except ValueError: min_b=0
    
    #determine average
    	try: average_b= float(sum(b_list)) / len(b_list)
    	except ZeroDivisionError: average_b=0
    
    #determine standard deviation
    	for i in range(len(m.atom)):
    		if m.atom[i]>average_b:
    			dummy.append((m.atom[i].b-average_b)**2)
    		if m.atom[i]<average_b:
    			dummy.append((average_b-m.atom[i].b)**2)
    	try: stdev_b= (sum(dummy) / (len(m.atom)-1))**(1/2.0)
    	except ZeroDivisionError: stdev_b=0
    
    #print values on screen
    	print "Number of atoms in '", selection,"': ", len(b_list)
    	print "Minimum and Maximum B-values: ", min_b, max_b
    	print "Average B-value: ", average_b
    	print "Standard deviation of the B-values: ", stdev_b
    	print "This data will be appended to cluster_count.txt"
    
    #write information to cluster_count.txt
    	f=open('cluster_count.txt','a')
    	print >>f, '%-10s %8d %8.3f %8.3f %8.3f %8.3f' % (selection, len(m.atom), min_b, max_b, average_b, stdev_b)
    	f.close()
    cmd.extend("cluster_count",cluster_count)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cluster_Count&oldid=10986](https://pymolwiki.org/index.php?title=Cluster_Count&oldid=10986)"


---

## CMPyMOL

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Software   
---|---  
Download  | <https://github.com/emptyewer/CMPyMOL>  
Author(s)  | [Venkatramanan Krishnamani](/index.php/User:Venky "User:Venky")  
License  | The MIT License (MIT)   
  
CMPyMOL is an add-on software to molecular visualization program PyMOL. It combines the protein 3D visualization capabilities of PyMOL and the protein's 2D contact map with an interactive interface for scientific analysis. This software is freely distributed under the MIT license for Linux and Mac OS X platforms. 

## Contents

  * 1 Website
  * 2 Version History
  * 3 Prerequisites
    * 3.1 Mac OS X
    * 3.2 Linux
  * 4 Installation
  * 5 Usage
  * 6 Software
    * 6.1 Inputs
    * 6.2 Overlays
    * 6.3 Plots
    * 6.4 Word of Caution
  * 7 Requests and Disclaimer
  * 8 Main Window
  * 9 Pairwise Aminoacid Heatmap
  * 10 License



## Website

<http://emptyewer.github.io/CMPyMOL/>

## Version History

  * Added support for reading multi-model PDB files. (Multi-Model PDB file for NMR structure or Trajectory from MD simulations.)
  * Supports displaying variance of contact points from a series of contact maps generated from multi-model PDB files.
  * CMPyMOL stores the calculated contact maps, heat maps and contact density information in a local SQLite database for fast and easy subsequent access.
  * Cleaner GUI.
  * Parallelized the code for contact map calculation for multi-frame PDB files.



## Prerequisites

  * Python 2.7
  * Python module dependency 
    * [wxpython](http://www.wxpython.org)
    * [matplotlib](http://matplotlib.org)
    * [python imaging library (PIL)](http://www.pythonware.com/products/pil/)
    * [numpy](http://www.numpy.org)
  * PyMOL. (It is recommended that the user add the PyMOL installation directory to the $PATH environment variable.)
  * Stride secondary structure assignment tool. This program can be downloaded from <http://webclu.bio.wzw.tum.de/stride/> and compiled into a stand-alone executable. It is recommended that the Stride executable or its installation directory is added to the $PATH environment variable. NOTE: If this executable is not detected in the $PATH variable, the secondary structure calculation will be disabled in CMPyMOL.



### Mac OS X

  * Users can install the python libraries using "easy_install" or "pip". It is recommended that the user use Enthought Canopy python distribution and management package downloaded from <https://www.enthought.com/products/canopy/>. This package includes a robust python library management software and a python IDE.
  * PyMOL 1.5.x can be installed using MacPorts <http://www.macports.org>. NOTE: This automatically adds the executable into the $PATH environment variable.



### Linux

  * The python dependencies and PyMOL can be installed using apt-get (aptitude) or a similar package management system.



## Installation

There is no need for installation of the script. Optionally, a standalone executable can be complied using "pyinstaller" or "py2exe" or "py2app" package, depending on the users operating system. 

## Usage
    
    
    python /<path to CMPyMOL directory>/CMPyMOL.py
    

  * This command will automatically invoke the PyMOL executable and the user is led through the rest of the program with a series of pop-up windows.



## Software

Clicking (left) and dragging a selection of contact points on the displayed contact map will highlight the corresponding residues in the PyMOL window (as red and blue colored atoms in spheres representation). In addition, several structural/biochemical properties can be overlaid on top of the contact map. The contact-map data can also be plotted in other representations. The calculated contact-map, heat-map and contact density information is stored in a local SQL database. Any subsequent access of the same PDB with matching parameters will be read from the database for fast access. The code for calculating contact map for trajectory files is parallelized for efficiency. 

### Inputs

  * Single-frame PDB files (local)
  * Multi-frame PDB trajectory files (local)
  * Multi-frame trajectory files should have the following format.


    
    
    MODEL X
    .
    .
    .
    ATOM ...
    ATOM ...
    .
    .
    .
    ENDMDL
    

NOTE: The PDB can include REMARKS, CRYST and other standard PDB information entries. The MODEL line is essential for the software to work properly (ENDMDL is optional). 

### Overlays

  * Secondary structure of the protein is overlaid as translucent strips over the contact map. This button won't be active if secondary structure calculation program stride is not found in the system path ($PATH). (Button: Secondary Structure)
  * Contact points where a Charge-Charge interaction occurs are highlighted. (Button: Charged Interactions)
  * Residues that interact via hydrophobic interaction are highlighted. (Button: Hydrophobic Interactions)
  * Contact regions that have a B-factor that is higher than a certain cutoff are highlighted (Button: B-factor). The b-factor cutoff can be varied using a slider (Slider).
  * Highlights a contact point/region where the pair of selected residues are in contact (selected by checking the checkboxes). If only one aminoacid is selected from the list, interaction site of the selected aminoacid with another one of the same type is highlighted. (List of checkboxes for each aminoacid)



### Plots

  * Pairwise Heat Map - Plots a 20x20 matrix of pairwise aminoacid interaction count.
  * Contacts Histogram - Plots the number of contacts around a given residue. Selecting a particular bar highlights the corresponding residue in the PyMOL window.
  * Variance Contact Map - For Multi-frame PDB files (trajectory), this button toggles the displays the variance contact map starting from the initial frame until the current frame. This view can be used to identifying the dynamic regions in a protein.



### Word of Caution

When using a multi-frame PDB file, the contact-map for the next frame(s) are being pre-calculated in the background (depending on the number of free CPU cores available). Clicking on "Next Frame" in rapid succession may lead to undesired results and/or a crash. 

In the event of a crash, delete the database that is created in the working directory and relaunch the program. 

## Requests and Disclaimer

Users are welcome to send me an email to request the addition of a specific feature or to report a bug. 

## Main Window

[![CMPyMOL.PNG](/images/d/dd/CMPyMOL.PNG)](/index.php/File:CMPyMOL.PNG)

[](/index.php/File:CMPyMOL.PNG "Enlarge")

## Pairwise Aminoacid Heatmap

[![Heatmap-CMPyMOL.PNG](/images/3/39/Heatmap-CMPyMOL.PNG)](/index.php/File:Heatmap-CMPyMOL.PNG)

[](/index.php/File:Heatmap-CMPyMOL.PNG "Enlarge")

## License
    
    
    #  The MIT License (MIT)
    # =======================
    # 
    # The PyMOL Plugin source code in this file is copyrighted, but you are
    # free to use and copy it as long as you don't change or remove any of
    # the copyright notices.
    # 
    # -----------------------------------------------------------------------------------
    # CMPyMOL
    # Copyright (C) 2013 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy of this 
    # software and associated documentation files (the "Software"), to deal in the Software 
    # without restriction, including without limitation the rights to use, copy, modify, 
    # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
    # permit persons to whom the Software is furnished to do so, subject to the following 
    # conditions:
    #
    # The above copyright notice and this permission notice shall be included in all copies 
    # or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
    # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    # OTHER DEALINGS IN THE SOFTWARE.
    

Retrieved from "[https://pymolwiki.org/index.php?title=CMPyMOL&oldid=12359](https://pymolwiki.org/index.php?title=CMPyMOL&oldid=12359)"


---

## CollapseSel

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

CollapseSel is a small utility function to compress strings for selections. So, if you have a selection with residues 1+2+3+4+20+21+22+100-120 this will return, 1-4+20-22+100-120. 

**CollapseIDs** is a small utility function to compress strings for an array of IDs. This does NOT have the logic for detecting duplicates, CollapseSel does. 

# Example
    
    
    run /dir/to/collapseSel.py
    fetch 1cll
    select EE, resn GLU
    print collapseSel("EE")
    #
    # and PyMOL should output:
    #
    # 6-7+11+14+31+45+47+54+67+82-84+87+104+114+119-120+123+127+139-140
    #
    

# The Code
    
    
    import pymol
    from pymol import stored
    
    def collapseIDs(ids):
    	"""
    	Helper function to make a smart list of IDs: eg turns 1+2+3+4+5 into 1-5.
    	"""
    	rVal = []
    	if len(ids)==0:
    		return ""
    
    	scanning=False
    	anchor = 0
    	start = 0
    	# 1-5 7-10 12 21-23
    	for cur in range(0,len(ids)-1):
    		if ids[cur]+1 == ids[cur+1]:
    			if scanning:
    				scanning=True
    				continue
    			else:
    				scanning=True
    				start = cur
    		else:
    			if scanning:
    				scanning=False
    				rVal.append(str(ids[start]) + "-" + str(ids[cur]))
    				start = cur
    			else:
    				scanning=False
    				rVal.append(str(ids[cur]))
    	if scanning:
    		rVal.append( str(ids[start]) + "-" + str(ids[cur+1]))
    	else:
    		rVal.append(str(ids[-1]))
    	return rVal
    
    def collapseSel(sel=None,lType="resi"):
            """
            collapseSel -- given a valid PyMOL selection and list type, return a collapsed
                    list of numbers corresponding to the lType.  For example, to compactly
                    print the residue identifiers for all the waters, try:
                            select theWaters, resn HOH
                            print collapseSel( "theWaters" )
    
                    This will convert: 1+2+3+4+10+11+12 into 1-4+10-12
    
            PARAMS
                    sel
                            The selection to collapse
                    lType
                            The identifier type: 'resi', 'ID', 'id', or any numerical property.
    
            RETURNS
                    a string of collapsed IDs.
            """
            if sel==None:
                    return ""
    
            stored.s = set()
            cmd.iterate( sel, "stored.s.add(int(float(%s)))" % lType)
            l = list(stored.s)
            l.sort()
            return "+".join(collapseIDs(list(l)))
    
    cmd.extend("collapseSel", collapseSel)
    

Retrieved from "[https://pymolwiki.org/index.php?title=CollapseSel&oldid=7213](https://pymolwiki.org/index.php?title=CollapseSel&oldid=7213)"


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

## Color cbcobj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [ Oscar Conchillo-Solé](/index.php?title=User:OscarConchillo-Sol%C3%A9&action=edit&redlink=1 "User:OscarConchillo-Solé \(page does not exist\)")  
License  |   
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 Code



# Overview

color_cbcobj 

Color all chains in different objects with a different color 

color_cbcobj(selection='(all)', first_color=0, quiet=0, _self=cmd): 

uses colors and order as in _color_cycle defined here: <https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/util.py>

quiet=0 means verbose (not quiet) 

# Usage
    
    
    color_cbcobj selection, [first_color=0..n, [quiet=0/1 ]]
    

where: 

  * **selection** can be any existing or newly-defined selection
  * **first_color** first color tu use from "_color_cycle" (first is 0)
  * **quiet** (default 0) print colored chains and models quiet=0 means verbose (not quiet)



# Examples
    
    
    PyMOL>color_cbcobj all and elem C,2,1
    

colors all chains in all objects differently, only carbon atoms, starting by color 154 (lightmagenta), does not show what has been colored 
    
    
    PyMOL>color_cbcobj all and elem C,2
     util.cbcobj: color 154 and model 6Y6C_BC', (chain B)
     Executive: Colored 670 atoms.
     util.cbcobj: color 6 and model 6Y6C_BC', (chain C)
     Executive: Colored 1131 atoms.
     util.cbcobj: color 9 and model 6Y6C_AD', (chain A)
     Executive: Colored 664 atoms.
     util.cbcobj: color 29 and model 6Y6C_AD', (chain D)
     Executive: Colored 1135 atoms.
    

# Code

Copy the following text and save it as color_cbcobj.py 
    
    
    from pymol import cmd, CmdException
    
    _color_cycle = [
        26   , # /* carbon */
        5    , # /* cyan */
        154  , # /* lightmagenta */
        6    , # /* yellow */
        9    , # /* salmon */
        29   , # /* hydrogen */
        11   , # /* slate */
        13   , # /* orange */
        10   , # /* lime */
        5262 , # /* deepteal */
        12   , # /* hotpink */
        36   , # /* yelloworange */
        5271 , # /* violetpurple */
        124  , # /* grey70 */
        17   , # /* marine */
        18   , # /* olive */
        5270 , # /* smudge */
        20   , # /* teal */
        5272 , # /* dirtyviolet */
        52   , # /* wheat */
        5258 , # /* deepsalmon */
        5274 , # /* lightpink */
        5257 , # /* aquamarine */
        5256 , # /* paleyellow */
        15   , # /* limegreen */
        5277 , # /* skyblue */
        5279 , # /* warmpink */
        5276 , # /* limon */
        53   , # /* violet */
        5278 , # /* bluewhite */
        5275 , # /* greencyan */
        5269 , # /* sand */
        22   , # /* forest */
        5266 , # /* lightteal */
        5280 , # /* darksalmon */
        5267 , # /* splitpea */
        5268 , # /* raspberry */
        104  , # /* grey50 */
        23   , # /* deepblue */
        51   , # /* brown */
        ]
    
    _color_cycle_len = len(_color_cycle)
    
    def color_cbcobj(selection='(all)', first_color=0, quiet=0, _self=cmd):
        '''
    DESCRIPTION
    
        Color all chains in different objects a different color
    
    SEE ALSO
    
        util.cbc, https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/util.py
        split_chains.py, https://pymolwiki.org/index.php/Split_chains
        '''
        quiet = int(quiet)
        count = int(first_color)
        models = cmd.get_object_list('(' + selection + ')')
        for model in models:
            for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
                if len(chain.split()) != 1:
                    chain = '"' + chain + '"'
                color = _color_cycle[count % _color_cycle_len]
                if not quiet:
                    print(" util.cbcobj: color %s and model %s', (chain %s)" % (color, model, chain))
                #_self.color(color, "(chain %s and (%s))" % (chain, model), quiet=quiet)
                _self.color(color, "(chain %s and (%s) and (%s))" % (chain, model, selection), quiet)
                count+=1
    
    cmd.extend('color_cbcobj', color_cbcobj)
    
    #color_cbcobj all and elem C,1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_cbcobj&oldid=13496](https://pymolwiki.org/index.php?title=Color_cbcobj&oldid=13496)"


---

## Color h

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script colors the selection passed in based on the hydrophobicity scale as defined by: 

    

    Source: <http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html>
    Amino acid scale: Normalized consensus hydrophobicity scale
    Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
    Reference: J. Mol. Biol. 179:125-142 (1984)

Or check out hydrophobicity coloring, with rTools from Kristian Rother. <http://www.rubor.de/pymol_extensions_de.html>

# The Code
    
    
    # color_h
    # -------
     
    # PyMOL command to color protein molecules according to the Eisenberg hydrophobicity scale
     
    #
    # Source: http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html
    # Amino acid scale: Normalized consensus hydrophobicity scale
    # Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
    # Reference: J. Mol. Biol. 179:125-142 (1984)
    #
    # Amino acid scale values:
    #
    # Ala:  0.620
    # Arg: -2.530
    # Asn: -0.780
    # Asp: -0.900
    # Cys:  0.290
    # Gln: -0.850
    # Glu: -0.740
    # Gly:  0.480
    # His: -0.400
    # Ile:  1.380
    # Leu:  1.060
    # Lys: -1.500
    # Met:  0.640
    # Phe:  1.190
    # Pro:  0.120
    # Ser: -0.180
    # Thr: -0.050
    # Trp:  0.810
    # Tyr:  0.260
    # Val:  1.080
    #
    # Usage:
    # color_h (selection)
    #
    from pymol import cmd
    
    def color_h(selection='all'):
            s = str(selection)
            print(s)
            cmd.set_color('color_ile',[0.996,0.062,0.062])
            cmd.set_color('color_phe',[0.996,0.109,0.109])
            cmd.set_color('color_val',[0.992,0.156,0.156])
            cmd.set_color('color_leu',[0.992,0.207,0.207])
            cmd.set_color('color_trp',[0.992,0.254,0.254])
            cmd.set_color('color_met',[0.988,0.301,0.301])
            cmd.set_color('color_ala',[0.988,0.348,0.348])
            cmd.set_color('color_gly',[0.984,0.394,0.394])
            cmd.set_color('color_cys',[0.984,0.445,0.445])
            cmd.set_color('color_tyr',[0.984,0.492,0.492])
            cmd.set_color('color_pro',[0.980,0.539,0.539])
            cmd.set_color('color_thr',[0.980,0.586,0.586])
            cmd.set_color('color_ser',[0.980,0.637,0.637])
            cmd.set_color('color_his',[0.977,0.684,0.684])
            cmd.set_color('color_glu',[0.977,0.730,0.730])
            cmd.set_color('color_asn',[0.973,0.777,0.777])
            cmd.set_color('color_gln',[0.973,0.824,0.824])
            cmd.set_color('color_asp',[0.973,0.875,0.875])
            cmd.set_color('color_lys',[0.899,0.922,0.922])
            cmd.set_color('color_arg',[0.899,0.969,0.969])
            cmd.color("color_ile","("+s+" and resn ile)")
            cmd.color("color_phe","("+s+" and resn phe)")
            cmd.color("color_val","("+s+" and resn val)")
            cmd.color("color_leu","("+s+" and resn leu)")
            cmd.color("color_trp","("+s+" and resn trp)")
            cmd.color("color_met","("+s+" and resn met)")
            cmd.color("color_ala","("+s+" and resn ala)")
            cmd.color("color_gly","("+s+" and resn gly)")
            cmd.color("color_cys","("+s+" and resn cys)")
            cmd.color("color_tyr","("+s+" and resn tyr)")
            cmd.color("color_pro","("+s+" and resn pro)")
            cmd.color("color_thr","("+s+" and resn thr)")
            cmd.color("color_ser","("+s+" and resn ser)")
            cmd.color("color_his","("+s+" and resn his)")
            cmd.color("color_glu","("+s+" and resn glu)")
            cmd.color("color_asn","("+s+" and resn asn)")
            cmd.color("color_gln","("+s+" and resn gln)")
            cmd.color("color_asp","("+s+" and resn asp)")
            cmd.color("color_lys","("+s+" and resn lys)")
            cmd.color("color_arg","("+s+" and resn arg)")
    cmd.extend('color_h',color_h)
    
    def color_h2(selection='all'):
            s = str(selection)
            print(s)
            cmd.set_color("color_ile2",[0.938,1,0.938])
            cmd.set_color("color_phe2",[0.891,1,0.891])
            cmd.set_color("color_val2",[0.844,1,0.844])
            cmd.set_color("color_leu2",[0.793,1,0.793])
            cmd.set_color("color_trp2",[0.746,1,0.746])
            cmd.set_color("color_met2",[0.699,1,0.699])
            cmd.set_color("color_ala2",[0.652,1,0.652])
            cmd.set_color("color_gly2",[0.606,1,0.606])
            cmd.set_color("color_cys2",[0.555,1,0.555])
            cmd.set_color("color_tyr2",[0.508,1,0.508])
            cmd.set_color("color_pro2",[0.461,1,0.461])
            cmd.set_color("color_thr2",[0.414,1,0.414])
            cmd.set_color("color_ser2",[0.363,1,0.363])
            cmd.set_color("color_his2",[0.316,1,0.316])
            cmd.set_color("color_glu2",[0.27,1,0.27])
            cmd.set_color("color_asn2",[0.223,1,0.223])
            cmd.set_color("color_gln2",[0.176,1,0.176])
            cmd.set_color("color_asp2",[0.125,1,0.125])
            cmd.set_color("color_lys2",[0.078,1,0.078])
            cmd.set_color("color_arg2",[0.031,1,0.031])
            cmd.color("color_ile2","("+s+" and resn ile)")
            cmd.color("color_phe2","("+s+" and resn phe)")
            cmd.color("color_val2","("+s+" and resn val)")
            cmd.color("color_leu2","("+s+" and resn leu)")
            cmd.color("color_trp2","("+s+" and resn trp)")
            cmd.color("color_met2","("+s+" and resn met)")
            cmd.color("color_ala2","("+s+" and resn ala)")
            cmd.color("color_gly2","("+s+" and resn gly)")
            cmd.color("color_cys2","("+s+" and resn cys)")
            cmd.color("color_tyr2","("+s+" and resn tyr)")
            cmd.color("color_pro2","("+s+" and resn pro)")
            cmd.color("color_thr2","("+s+" and resn thr)")
            cmd.color("color_ser2","("+s+" and resn ser)")
            cmd.color("color_his2","("+s+" and resn his)")
            cmd.color("color_glu2","("+s+" and resn glu)")
            cmd.color("color_asn2","("+s+" and resn asn)")
            cmd.color("color_gln2","("+s+" and resn gln)")
            cmd.color("color_asp2","("+s+" and resn asp)")
            cmd.color("color_lys2","("+s+" and resn lys)")
            cmd.color("color_arg2","("+s+" and resn arg)")
    cmd.extend('color_h2',color_h2)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_h&oldid=13186](https://pymolwiki.org/index.php?title=Color_h&oldid=13186)"


---

## Color Objects

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | Gareth Stockwell   
License  |   
      
    
    #####################################################################
    #
    # Colour by object
    #
    #####################################################################
     
    from pymol import cmd
    
    def color_obj(rainbow=0):
     
            """
     
    AUTHOR 
     
            Gareth Stockwell
     
    USAGE
     
            color_obj(rainbow=0)
     
            This function colours each object currently in the PyMOL heirarchy
            with a different colour.  Colours used are either the 22 named
            colours used by PyMOL (in which case the 23rd object, if it exists,
            gets the same colour as the first), or are the colours of the rainbow
     
    SEE ALSO
    
            util.color_objs()
            """
     
            # Process arguments
            rainbow = int(rainbow)
     
            # Get names of all PyMOL objects
            obj_list = cmd.get_names('objects')
     
            if rainbow:
     
               print "\nColouring objects as rainbow\n"
     
               nobj = len(obj_list)
     
               # Create colours starting at blue(240) to red(0), using intervals
               # of 240/(nobj-1)
               for j in range(nobj):
                  hsv = (240-j*240/(nobj-1), 1, 1)
                  # Convert to RGB
                  rgb = hsv_to_rgb(hsv)
                  # Define the new colour
                  cmd.set_color("col" + str(j), rgb)
                  print obj_list[j], rgb
                  # Colour the object
                  cmd.color("col" + str(j), obj_list[j])
     
            else:
     
               print "\nColouring objects using PyMOL defined colours\n"
     
               # List of available colours
               colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
               'wheat', 'white', 'grey' ]
               ncolours = len(colours)
     
               # Loop over objects
               i = 0
               for obj in obj_list:
                  print "  ", obj, colours[i]
                  cmd.color(colours[i], obj)
                  i = i+1
                  if(i == ncolours):
                     i = 0
     
     
    # HSV to RGB routine taken from Robert L. Campbell's color_b.py script
    #   See http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
    # Original algorithm from: http://www.cs.rit.edu/~ncs/color/t_convert.html
    def hsv_to_rgb(hsv):
     
            h = float(hsv[0])
            s = float(hsv[1])
            v = float(hsv[2])
     
            if( s == 0 ) :
                    #achromatic (grey)
                    r = g = b = v
     
            else:
                    # sector 0 to 5
                    h = h/60.            
                    i = int(h)
                    f = h - i                       # factorial part of h
                    #print h,i,f
                    p = v * ( 1 - s )
                    q = v * ( 1 - s * f )
                    t = v * ( 1 - s * ( 1 - f ) )
     
                    if i == 0:
                            (r,g,b) = (v,t,p)
                    elif i == 1:
                            (r,g,b) = (q,v,p)
                    elif i == 2:
                            (r,g,b) = (p,v,t)
                    elif i == 3:
                            (r,g,b) = (p,q,v)
                    elif i == 4:
                            (r,g,b) = (t,p,v)
                    elif i == 5:
                            (r,g,b) = (v,p,q)
                    else:
                            (r,g,b) = (v,v,v)
                            print "error, i not equal 1-5"
     
            return [r,g,b]
     
     
     
    # Add color_obj to the PyMOL command list 
    cmd.extend("color_obj",color_obj)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Color_Objects&oldid=10269](https://pymolwiki.org/index.php?title=Color_Objects&oldid=10269)"


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

## ConnectedCloud

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 The Code



# Overview

ConnectedCloud grows a cloud of objects that are "connected" through single-linkage clustering. Two objects are "connected" if they are within the user-defined cutoff, radius. It starts with a selection, and grows until no more objects are withing _radius_ Angstroms of the cloud. 

Let's say we want to find the cloud of docking guesses (blue lines) that are connected to the lowest energy docking pose (red sticks). I call, 
    
    
    connectedCloud lig1, lig*, radius=1, doPrint=True
    

and get 

  * [![Find the cloud growing from the red sticks. Consider all blue lines for inclusion.](/images/1/18/ConnectedCloud.png)](/index.php/File:ConnectedCloud.png "Find the cloud growing from the red sticks. Consider all blue lines for inclusion.")

Find the cloud growing from the red sticks. Consider all blue lines for inclusion. 

  * [![Those in the cloud are shown in red. Notice the small blue cloud above the red was not selected. At no point were those atoms close enough to the currently growing cloud, to be considered.](/images/6/62/ConnectedCloud2.png)](/index.php/File:ConnectedCloud2.png "Those in the cloud are shown in red. Notice the small blue cloud above the red was not selected. At no point were those atoms close enough to the currently growing cloud, to be considered.")

Those in the cloud are shown in red. Notice the small blue cloud above the red was not selected. At no point were those atoms close enough to the currently growing cloud, to be considered. 




# Syntax
    
    
    connectedCloud origSel, subSel, [radius=1, [doPrint=True]]
    

  * origSel -- the original selection to grow from
  * subSel -- only consider these objects when growing
  * radius -- defines "connectedness". Anything <= radius Angstroms from the growing cloud and in subSel is considered
  * doPrint -- print some misc. info



# Examples
    
    
    # Find all objects in the set "lig*" (all ligands from a docking run, for example) that are
    # in the connected cloud of "connectedness" set to 1 Ang.  Print the number of atoms in the
    # selection and the surface area of the selection (if it were an object).
    connectedCloud lig427, lig*, radius=1, doPrint=True
    
    # select everything within 0.85 of lig2467 and don't print anything
    connectedCloud lig2467, *, radius=0.85, doPrint=False
    

# The Code
    
    
    # -*- coding: utf-8 -*-
    def connectedCloud( origSel, subSel, radius=1, doPrint=True ):
    	"""
    	connectedCloud -- find all atoms, by object in subSel connected to
    		origSel by radius Ang.  In other words, start from origSel
    		and find all objects <=radius Ang. from it & add that to the
    		result.  Once nothing more is added, you have your cloud of
    		points.
    
    	origSel,
    		the selection around which we start expanding.
    
    	subSel,
    		the subset of objects to consider
    
    	radius,
    		the maximum radius in which to consider; anything >radius
    		Ang. from any objet in the current selection is ignored (unless
    		it is added later by some other path).
    
    	doPrint
    		print some values before returning the selection name
    
    	Author:
    		Jason Vertrees, 2009.
    	
    	"""
    	cName = "__tempCloud__"
    	cmd.select(cName, origSel)
    
    	# we keep track of cloud size based on the number of atoms in the selection; oldCount
    	# is the previous iteration and curCount is after we add everything for this iteration.
    	oldCount=0
    	curCount=cmd.count_atoms( cName )
    
    	# flag to stop at 1000 iterations. if something goes wrong we don't want to hijack the computer.
    	flag=0
    
    	# while something was added in this iteration
    	while oldCount != curCount:
    		flag += 1
    		# grow current cloud
    		cmd.select( cName, "bo. " + subSel + " within " + radius + " of " + cName )
    		# update atom counts in the cloud
    		oldCount=curCount
    		curCount=cmd.count_atoms(cName)
    		
    		if ( flag > 1000 ):
    			print "Warning: I gave up.  I grew the cloud 1000 times and it kept growing."
    			print "Warning: If this is valid (weird), then increases the 'flag' in the source"
    			print "Warning: code to something >1000."
    			break
    
    	# create the object for getting its surface area
    	cmd.create("connectedCloud", cName)
    
    	# print cloud info
    	if doPrint:
    		print "Number of atoms in cloud:", str(cmd.count_atoms(cName))
    		print "Surface area of cloud:", str(cmd.get_area("connectedCloud"))
    
    	# delete our work 
    	cmd.delete("connectedCloud")
    	# rename and return
    	cmd.set_name(cName, "connectedCloud")
    	return "connectedCloud"
    
    cmd.extend("connectedCloud", connectedCloud)
    

Retrieved from "[https://pymolwiki.org/index.php?title=ConnectedCloud&oldid=9962](https://pymolwiki.org/index.php?title=ConnectedCloud&oldid=9962)"


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

## Count molecules in selection

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will return the number of distinct molecular objects in a given selection. If no selection is given, the script assumes "sele" the named mouse selection. 

# Example
    
    
    # run the script
    
    run count_molecules_in_sel.py
    
    # grab a protein from the PDB
    
    fetch 2f56, async=0
    
    # select all the urea molecules
    
    select resn URE
    
    # should count 11 molecules
    
    count_molecules_in_selection 
    
    # change the selection and try again
    
    select polymer
    
    # should return 3, as there's 3 chains
    
    count_molecules_in_selection
    

# The Code
    
    
    import pymol
    from pymol import cmd
    
    def count_mols_in_sel(sel="sele"):
        """
        Returns the number of distinct molecules in a given selection.
        """
    
        sel_copy = "__selcopy"
    
        cmd.select(sel_copy, sel)
    
        num_objs = 0
    
        atoms_in_sel = cmd.count_atoms(sel_copy)
    
        while atoms_in_sel > 0:
    
            num_objs += 1
    
            cmd.select(sel_copy, "%s and not (bm. first %s)" % (sel_copy, sel_copy))
    
            atoms_in_sel = cmd.count_atoms(sel_copy)
    
        print "There are %d distinct molecules in the selection '%s'." % (num_objs, sel)
    
        return num_objs
    
    
    cmd.extend("count_molecules_in_selection", count_mols_in_sel)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Count_molecules_in_selection&oldid=10860](https://pymolwiki.org/index.php?title=Count_molecules_in_selection&oldid=10860)"


---

## CreateAtom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A script to create an atom C at a point some distance d from a pair of atoms (A, B), along the line of the bond A-B. The main function takes a modelName (usually, the name of the file loaded, like "1pqr" or "peptide"), a distance, and some parameters to identify the atoms A, B. 

Use like: 
    
    
    createAtomAlongBond("gly", 3, 23, "H", 23, "N", "O")

. 
    
    
    import cmd
    from chempy import models, cpv
    
    """
    Create an atom at a distance 'distance' along the bond between atomA and atomB
    """
    def createAtomAlongBond(modelName, distance, resiA, atomNameA, resiB, atomNameB, atomNameC):
        model = cmd.get_model(modelName)
        p1 = getAtomCoords(model, str(resiA), atomNameA)
        p2 = getAtomCoords(model, str(resiB), atomNameB)
        if p1 is None:
            print "atom not found!", modelName, resiA, atomNameA
        elif p2 is None:
            print "atom not found!", modelName, resiB, atomNameB
        else:
            p3 = calculateNewPoint(p1, p2, distance)
    
            # the details of the new atom
            atomDetails = {}
            atomDetails['residueName'] = "HOH"
            atomDetails['residueNumber'] = "1"
            atomDetails['symbol'] = "O"
            atomDetails['name'] = atomNameC
            atomDetails['coords'] = p3
            
            # make an atom with index n+1 and chain "X"
            newAtom = makeAtom(model.nAtom + 1, atomDetails, "X")
            model.add_atom(newAtom)
            model.update_index()
            cmd.load_model(model, "newpeptide")
    
    def getAtomCoords(model, resi, atomName):
        for a in model.atom:
            if a.resi == resi and a.name == atomName:
                return a.coord
        return None
    
    def calculateNewPoint(p1, p2, distance):
        v1 = cpv.normalize(cpv.sub(p1, p2))
        return cpv.add(p1, cpv.scale(v1, distance))
    
    def makeAtom(index, atomDetails, chain):
        atom = chempy.Atom()
        atom.index = index
        atom.name = atomDetails['name']
        atom.symbol = atomDetails['symbol']
        atom.resn = atomDetails['residueName']
        atom.chain = chain
        atom.resi = atomDetails['residueNumber']
        atom.resi_number = int(atomDetails['residueNumber'])
        atom.coord = atomDetails['coords']
        atom.hetatm = False
        return atom
    
    cmd.extend("createAtomAlongBond", createAtomAlongBond)
    

Retrieved from "[https://pymolwiki.org/index.php?title=CreateAtom&oldid=6387](https://pymolwiki.org/index.php?title=CreateAtom&oldid=6387)"


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

## Cubes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/cubes.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cubes.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the **cubes** and **tetrahedra** commands which create CGOs with cube or tetrahedron representations for all atoms in selection. 

## Usage
    
    
    cubes [ selection [, name [, state [, scale [, atomcolors ]]]]]
    
    
    
    tetrahedra [ selection [, name [, state [, scale [, atomcolors ]]]]]
    

## Example
    
    
    fetch 1rx1, async=0
    cubes solvent & b < 50
    tetrahedra solvent & b > 50
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cubes&oldid=13899](https://pymolwiki.org/index.php?title=Cubes&oldid=13899)"


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

## Distancetoatom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/distancetoatom.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/distancetoatom.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar") and [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![get distance to atoms within cutoff](/images/c/c2/Distancetoatom_0.png)](/index.php/File:Distancetoatom_0.png "get distance to atoms within cutoff")

## Contents

  * 1 About distancetoatom
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
  * 4 SEE ALSO



## About distancetoatom

**distancetoatom** prints all distances between a specified atom, coordinate or group selection center and all atoms within cutoff distance that are part of the selection.  
All coordinates and distances can be saved in a csv-style text file report and/or can be stored in a (custom) atom property, if defined. 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage
    
    
        distancetoatom [ origin [, cutoff [, filename [, selection [, state [, property_name [, coordinates [, decimals [, sort [, quiet ]]]]]]]]]]
    

|   
---|---  
  
### Arguments

Arguments for distancetoatom   
---  
Keyword  | Default  | What it does   
origin  | pk1  | "origin" defines the coordinates for the origin and can be: 1\. a list with coordinates [x,y,z]  
2\. a single atom selection string  
3\. a multi-atom selection string (center will be used)   
cutoff  | 10  | "cutoff" defines the maximum distance, atoms further away are not considered   
filename  | None  | "filename" is the name of a report file that will be created (optional). set to e.g. 'report.txt' to create a report. This file is CSV style and can be imported into EXCEL.  
(omit or set to "", None, 0 or False to disable)   
selection  | all  | only atoms within "selection" (and cutoff distance) will be used for calculation   
state  | 0  | "state" is the state that will be used for calculations use 0 (omit), to automatically use the current state   
property_name  | p.dist  | "property_name" defines a (custom) property in which the calculated distance will be stored. This can be used e.g. for labeling or spectrum coloring etc. (cf. examples)  
NB! older PyMOL versions only support b or q and the distance will overwrite this property if set.   
coordinates  | 0  | "coordinates" toggles whether, besides distance, the atom coordinates will be included in the report.   
decimals  | 3  | "decimals" is the number of decimal places calculated. Note that PDB files do not support a higher resolution than 3.   
sort  | 1  | "sort" defines how the output will be sorted: 1: ascending (default)  
0: no sorting (by selection)  
-1: descending   
quiet  | 1  | toggle verbosity   
  
  


## Examples
    
    
    # make function available to PyMOL
    import distancetoatom
    
    ##############################
    # Example 1: print all distances
    frag LYS
    edit name CA
    distancetoatom
    
    ##############################
    # Example 2: print all distances to file
    fetch 1cll, async=0
    distancetoatom origin=/1cll//A/CA`150/CA, filename=distances.txt, selection=elem O, property_name=b
    # the file is CSV-format and can be imported, e.g. to EXCEL
    
    # format
    hide everything
    select byres (resi 150 expand 5 and not resn hoh) 
    show_as sticks, sele
    util.cbaw sele
    show_as spheres, resi 150
    zoom sele
    
    ##########
    # Label by stored distance (property_name)
    label sele and elem O, "%.2f"%b
    
    
    ##############################
    # Example 3: color an object by distance to a coordinate
    
    fetch 1cll, async=0
    distancetoatom [0,0,0], cutoff=1000, property_name=b
    # the distance is stored in the b-factor in this case
    # newer PyMOL versions support custom properties, e.g. "p.dist"
    spectrum b, rainbow_rev
    

|  [![Example 2](/images/5/5c/Distancetoatom_2.png)](/index.php/File:Distancetoatom_2.png "Example 2")  
  
[![Example 3](/images/3/3a/Distancetoatom_1.png)](/index.php/File:Distancetoatom_1.png "Example 3")  
---|---  
  
  


## SEE ALSO

  * [Distance](/index.php/Distance "Distance")
  * [Get Distance](/index.php/Get_Distance "Get Distance")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")



Retrieved from "[https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939](https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939)"


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

## DrawGridBox

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/drawgridbox.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/drawgridbox.py)  
Author(s)  | [Cunliang Geng](/index.php?title=User:Clgeng&action=edit&redlink=1 "User:Clgeng \(page does not exist\)")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the [drawgridbox](/index.php?title=Drawgridbox&action=edit&redlink=1 "Drawgridbox \(page does not exist\)") command, which draw grid boxes around a selection. 

## Usage
    
    
     drawgridbox [ selection,  [nx, [ny, [nz, [padding, [lw, [r, [g, [b]]]]]]]]]
    

  


## Prameters

  * **selection** = str: atom selection {default: all}
  * **nx** = int: number of grids on axis X {default: 10}
  * **ny** = int: number of grids on axis Y {default: 10}
  * **nz** = int: number of grids on axis Z {default: 10}
  * **padding** = float: default 0.0
  * **lw** = float: line width {default: 2.0}
  * **r** = float: red color component, valid range is [0.0, 1.0] {default 1.0}
  * **g** = float: green color component, valid range is [0.0, 1.0] {default 1.0}
  * **b** = float: blue color component, valid range is [0.0, 1.0] {default 1.0}



  


## Examples
    
    
    load drawgridbox.py
    fetch 1cbh
    show surface
    drawgridbox 1cbh, nx=5, ny=5, nz=5,  lw=1, g=0, b=0
    

[![Drawgridbox.png](/images/5/5a/Drawgridbox.png)](/index.php/File:Drawgridbox.png)

  

    
    
    drawgridbox 1cbh, nx=5, ny=5, nz=5, padding=5, lw=1, r=0, g=0
    

[![Drawgridbox padding.png](/images/3/32/Drawgridbox_padding.png)](/index.php/File:Drawgridbox_padding.png)

Retrieved from "[https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904](https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904)"


---

## DrawBoundingBox

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Example
  * 3 Installation
  * 4 See Also



# Overview

Draws a bounding box around a given selection. 

  * [![Example of a bounding box](/images/5/51/DrawMinBB.png)](/index.php/File:DrawMinBB.png "Example of a bounding box")

Example of a bounding box 

  * [![Two bounding boxes, one with 0 padding, the other with 10 Ang padding.](/images/6/63/Bb_with_padding.png)](/index.php/File:Bb_with_padding.png "Two bounding boxes, one with 0 padding, the other with 10 Ang padding.")

Two bounding boxes, one with 0 padding, the other with 10 Ang padding. 




# Example
    
    
    run ~/drawBoundingBox.py
    fetch 1jsd
    drawBoundingBox 1jsd, r=0.33, g=0.80
    
    # example from above w/padding, draw it light blue
    drawBoundingBox padding=10, r=0.5, g=0.8, b=1.0
    

# Installation

Copy the source code to your computer, and execute it by issuing the command "run /path/to/drawBoundingBox.py" in PyMOL. 
    
    
    # -*- coding: utf-8 -*-                                                                                     
    from pymol.cgo import *                                                                                     
    from pymol import cmd                                                                                       
    from random import randint                                                                                  
    
    #############################################################################
    #                                                                            
    # drawBoundingBox.py -- Draws a box surrounding a selection 
    #
    #                                                                            
    # AUTHOR: Jason Vertrees                                                   
    # DATE  : 2/20/2009                                                          
    # NOTES : See comments below.                                                
    #                                                                            
    #############################################################################
    def drawBoundingBox(selection="(all)", padding=0.0, linewidth=2.0, r=1.0, g=1.0, b=1.0):     
            """                                                                  
            DESCRIPTION                                                          
                    Given selection, draw the bounding box around it.          
    
            USAGE:
                    drawBoundingBox [selection, [padding, [linewidth, [r, [g, b]]]]]
    
            PARAMETERS:
                    selection,              the selection to enboxen.  :-)
                                            defaults to (all)
       
                    padding,                defaults to 0
    
                    linewidth,              width of box lines
                                            defaults to 2.0
    
                    r,                      red color component, valid range is [0.0, 1.0]
                                            defaults to 1.0                               
    
                    g,                      green color component, valid range is [0.0, 1.0]
                                            defaults to 1.0                                 
    
                    b,                      blue color component, valid range is [0.0, 1.0]
                                            defaults to 1.0                                
    
            RETURNS
                    string, the name of the CGO box
    
            NOTES
                    * This function creates a randomly named CGO box that minimally spans the protein. The
                    user can specify the width of the lines, the padding and also the color.                            
            """                                                                                                    
    
            ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)
    
            print "Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)
    
            minX = minX - float(padding)
            minY = minY - float(padding)
            minZ = minZ - float(padding)
            maxX = maxX + float(padding)
            maxY = maxY + float(padding)
            maxZ = maxZ + float(padding)
    
            if padding != 0:
                     print "Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)
    
            boundingBox = [
                    LINEWIDTH, float(linewidth),
    
                    BEGIN, LINES,
                    COLOR, float(r), float(g), float(b),
    
                    VERTEX, minX, minY, minZ,       #1
                    VERTEX, minX, minY, maxZ,       #2
    
                    VERTEX, minX, maxY, minZ,       #3
                    VERTEX, minX, maxY, maxZ,       #4
    
                    VERTEX, maxX, minY, minZ,       #5
                    VERTEX, maxX, minY, maxZ,       #6
    
                    VERTEX, maxX, maxY, minZ,       #7
                    VERTEX, maxX, maxY, maxZ,       #8
    
    
                    VERTEX, minX, minY, minZ,       #1
                    VERTEX, maxX, minY, minZ,       #5
    
                    VERTEX, minX, maxY, minZ,       #3
                    VERTEX, maxX, maxY, minZ,       #7
    
                    VERTEX, minX, maxY, maxZ,       #4
                    VERTEX, maxX, maxY, maxZ,       #8
    
                    VERTEX, minX, minY, maxZ,       #2
                    VERTEX, maxX, minY, maxZ,       #6
    
    
                    VERTEX, minX, minY, minZ,       #1
                    VERTEX, minX, maxY, minZ,       #3
    
                    VERTEX, maxX, minY, minZ,       #5
                    VERTEX, maxX, maxY, minZ,       #7
    
                    VERTEX, minX, minY, maxZ,       #2
                    VERTEX, minX, maxY, maxZ,       #4
    
                    VERTEX, maxX, minY, maxZ,       #6
                    VERTEX, maxX, maxY, maxZ,       #8
    
                    END
            ]
    
            boxName = "box_" + str(randint(0,10000))
            while boxName in cmd.get_names():
                    boxName = "box_" + str(randint(0,10000))
    
            cmd.load_cgo(boundingBox,boxName)
            return boxName
    
    cmd.extend ("drawBoundingBox", drawBoundingBox)
    

# See Also

[Bounding_Box](/index.php/Bounding_Box "Bounding Box")

Retrieved from "[https://pymolwiki.org/index.php?title=DrawBoundingBox&oldid=10935](https://pymolwiki.org/index.php?title=DrawBoundingBox&oldid=10935)"


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

## Dump2CGO

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

Dumps a PyMOL object to a CGO object. 

# The Code
    
    
    from pymol import cmd
    from pymol.cgo import *
    
    def dump2surfaceCGO():
        CGOobj = []
        dumpedFile = open("dump.tmp").read()
        for block in dumpedFile.split('\n\n'):
            CGOobj.append(BEGIN)
            CGOobj.append(TRIANGLES)
    
            for line in block.split('\n'):
                if line == '':
                    continue
    
                vals = line.split()
                CGOobj.append(NORMAL)
                CGOobj.append(float(vals[3]))
                CGOobj.append(float(vals[4]))
                CGOobj.append(float(vals[5]))
                CGOobj.append(VERTEX)
                CGOobj.append(float(vals[0]))
                CGOobj.append(float(vals[1]))
                CGOobj.append(float(vals[2]))
    
            CGOobj.append(END)
        return CGOobj
    
    def dump2meshCGO():
        CGOobj = []
        dumpedFile = open("dump.tmp").read()
        for block in dumpedFile.split('\n\n'):
            CGOobj.append(BEGIN)
            CGOobj.append(LINE_STRIP)
    
            for line in block.split('\n'):
                if line == '':
                    continue
    
                CGOobj.append(VERTEX)
                vals = line.split()
    
                CGOobj.append(float(vals[0]))
                CGOobj.append(float(vals[1]))
                CGOobj.append(float(vals[2]))
    
            CGOobj.append(END)
        return CGOobj
    
    def getType(objname):
        session = cmd.get_session()['names']
        for obj in session:
            if obj == None:
                continue
            if obj[0] != objname:
                continue
            return obj[4]
        return -1
    
    
    def dump2CGO(obj):
        cmd.dump("dump.tmp", obj)
        type = getType(obj)
        cgo = []
        if (type == 3): # Mesh
            cgo = dump2meshCGO()
        elif (type == 7): #Surface
            cgo = dump2surfaceCGO()
        else:
            print("Unknown type")
            return
    
        cmd.load_cgo(cgo, "CGO " + obj)
    
    cmd.extend('dump2CGO', dump2CGO)
    cmd.auto_arg[0]['dump2CGO'] = [cmd.object_sc, 'object', '']
    

Retrieved from "[https://pymolwiki.org/index.php?title=Dump2CGO&oldid=12733](https://pymolwiki.org/index.php?title=Dump2CGO&oldid=12733)"


---

## Dynamic mesh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/dynamic_mesh.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/dynamic_mesh.py)  
Author(s)  | [Takanori Nakane](/index.php?title=User:TakanoriNakane&action=edit&redlink=1 "User:TakanoriNakane \(page does not exist\)")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**dynamic_mesh** displays [isomesh](/index.php/Isomesh "Isomesh") around the center of the view. When the view is moved, the isomesh will be updated automatically. You can also change contour level by PageDown/PageUp keys. This script is intended to implement interface similar to Coot for examing electron density maps. 

Note: PyMOL's [Density Wizard](/index.php?title=Density_Wizard&action=edit&redlink=1 "Density Wizard \(page does not exist\)") (Menu > Wizard > Density) provides similar functionality. It is implemented using [wizard](/index.php/Wizard "Wizard") framework, while this uses [CallBack](/index.php?title=CallBack&action=edit&redlink=1 "CallBack \(page does not exist\)") object. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    dynamic_mesh map_name [, level [, radius [, name [ sym_source ]]]]
    

## Arguments

  * map_name = string: name of volumetric object(map) to display
  * level = float: contour level of isomesh {default: 1.0}
  * radius = float: radius of isomesh around the center of the view {default: 8}
  * name = string: name of mesh object {default: dynamic_mesh}
  * sym_source = string: name of object from which symmetry information is derived {default: map_name}



## Example
    
    
    fetch 1HWK, async=1
    fetch 1HWK, 1hwk_map, type=2fofc, async=1
    run dynamic_mesh.py
    dynamic_mesh 1hwk_map, sym_source=1hwk
    show sticks, resn 117
    show ribbon
    zoom chain A and resn 117
    

Note: On PyMOL <= 1.4, you have to download the electron density map from the Uppsala Electron Density Server manually. 

## See Also

  * [isomesh](/index.php/Isomesh "Isomesh")
  * [get_position](/index.php/Get_position "Get position")
  * [Density Wizard](/index.php?title=Density_Wizard&action=edit&redlink=1 "Density Wizard \(page does not exist\)")
  * [map_auto_expand_sym](/index.php/Map_auto_expand_sym "Map auto expand sym")
  * Coot <http://lmb.bioch.ox.ac.uk/coot/>



Retrieved from "[https://pymolwiki.org/index.php?title=Dynamic_mesh&oldid=13903](https://pymolwiki.org/index.php?title=Dynamic_mesh&oldid=13903)"


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

## Ellipsoid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script provides methods that create [cgos](/index.php?title=Cgo&action=edit&redlink=1 "Cgo \(page does not exist\)") as triangles. It uses code that is ported from [this c++ code](http://www.gamedev.net/reference/articles/article1172.asp) and seems to be correct! 

Here is the script. The last four lines show this in use, by making an ellipse and a toroid and loading them into pymol. This is done most easily by something like "cmd.load_cgo(makeEllipsoid(1, 1, 1, 2, 3, 4), 'ellipsoid')" which makes an ellipsoid at x, y, z = 1, 1, 1 and dimensions 2, 3, 4 and called 'ellipsoid'. 
    
    
    from pymol.cgo import BEGIN, COLOR, TRIANGLES, VERTEX, NORMAL, END
    from pymol import cmd
    
    def signOfFloat(f):
            if f < 0: return -1
            if f > 0: return 1
            return 0
    
    def sqC(v, n):
            return signOfFloat(math.cos(v)) *  math.pow(math.fabs(math.cos(v)), n)
    
    def sqCT(v, n, alpha):
            return alpha + sqC(v, n)
    
    def sqS(v, n):
            return signOfFloat(math.sin(v)) * math.pow(math.fabs(math.sin(v)), n)
    
    def sqEllipsoid(x, y, z, a1, a2, a3, u, v, n, e):
            x = a1 * sqC(u, n) * sqC(v, e) + x
            y = a2 * sqC(u, n) * sqS(v, e) + y
            z = a3 * sqS(u, n) + z
            nx = sqC(u, 2 - n) * sqC(v, 2 - e) / a1
            ny = sqC(u, 2 - n) * sqS(v, 2 - e) / a2
            nz = sqS(u, 2 - n) / a3
            return x, y, z, nx, ny, nz
    
    def sqToroid(x, y, z, a1, a2, a3, u, v, n, e, alpha):
            a1prime = 1.0 / (a1 + alpha)
            a2prime = 1.0 / (a2 + alpha)
            a3prime = 1.0 / (a3 + alpha)
            x = a1prime * sqCT(u, e, alpha) * sqC(v, n)
            y = a2prime * sqCT(u, e, alpha) * sqS(v, n)
            z = a3prime * sqS(u, e)
            nx = sqC(u, 2 - e) * sqC(v, 2 - n) / a1prime
            ny = sqC(u, 2 - e) * sqS(v, 2 - n) / a2prime
            nz = sqS(u, 2 - e) / a3prime
            return x, y, z, nx, ny, nz
    
    def makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, n, e, u1, u2, v1, v2, u_segs, v_segs, color=[0.5, 0.5, 0.5]):
    
            r, g, b = color
    
            # Calculate delta variables */
            dU = (u2 - u1) / u_segs
            dV = (v2 - v1) / v_segs
    
            o = [ BEGIN, TRIANGLES ]
    
            U = u1
            for Y in range(0, u_segs):
                    # Initialize variables for loop */
                    V = v1
                    for X in range(0, v_segs):
                            # VERTEX #1 */
                            x1, y1, z1, n1x, n1y, n1z = sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, e)
                            x2, y2, z2, n2x, n2y, n2z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, e)
                            x3, y3, z3, n3x, n3y, n3z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e)
                            x4, y4, z4, n4x, n4y, n4z = sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, e)
    
                            o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
    
                            # Update variables for next loop */
                            V += dV
                    # Update variables for next loop */
                    U += dU
            o.append(END)
            return o
    
    def makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, n, e, u1, u2, v1, v2, u_segs, v_segs, color=[0.5, 0.5, 0.5]):
    
            r, g, b = color
    
            # Calculate delta variables */
            dU = (u2 - u1) / u_segs
            dV = (v2 - v1) / v_segs
    
            o = [ BEGIN, TRIANGLES ]
    
            U = u1
            for Y in range(0, u_segs):
                    # Initialize variables for loop */
                    V = v1
                    for X in range(0, v_segs):
                            # VERTEX #1 */
                            x1, y1, z1, n1x, n1y, n1z = sqToroid(x, y, z, a1, a2, a3, U, V, n, e, alpha)
                            x2, y2, z2, n2x, n2y, n2z = sqToroid(x, y, z, a1, a2, a3, U + dU, V, n, e, alpha)
                            x3, y3, z3, n3x, n3y, n3z = sqToroid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e, alpha)
                            x4, y4, z4, n4x, n4y, n4z = sqToroid(x, y, z, a1, a2, a3, U, V + dV, n, e, alpha)
    
                            o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
    
                            # Update variables for next loop */
                            V += dV
                    # Update variables for next loop */
                    U += dU
            o.append(END)
            return o
    
    def makeEllipsoid(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeCylinder(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeSpindle(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeDoublePyramid(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 2.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makePillow(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 0.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)
    
    def makeRoundCube(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.2, 0.2, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeToroid(x, y, z, a1, a2, a3, alpha):
                    return makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, 1.0, 1.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)
    
    x, y, z, rx, ry, rz = 1, 1, 1, 1, 2, 3
    cmd.load_cgo(makeEllipsoid(x, y, z, rx, ry, rz), 'ellipsoid-cgo')
    x, y, z, rx, ry, rz = 1, 1, 1, 8, 2, 2
    cmd.load_cgo(makeToroid(x, y, z, rx, ry, rz, 3), 'toroid-cgo')
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoid&oldid=6396](https://pymolwiki.org/index.php?title=Ellipsoid&oldid=6396)"


---

## Expand To Surface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Example
  * 3 The Code
  * 4 See Also



# Overview

This little script will color the surface of your protein the same color as the b-factors from the alpha carbons of each amino acid. This script is useful when you want to color a protein by the b-factor column, but all you have in the b-factor for each alpha carbon. If you try to color a surface based only on the alpha carbon b-factors, then the surface isn't rendered as you typically want. 

**Note:** This script, once run in PyMOL, is called **e2s**. See examples, below. 

# Example
    
    
    # in PyMOL type the following
    
    # fetch 1cll
    fetch 1cll
    
    # load your b-factor column; for "1cll" you can simulate fake data
    # by writing 144 random integers to /tmp/bb.  If you have some real
    # data, then put that there.  1 row for each alpha carbon.
    f = open('/tmp/bb', 'r').readlines();
    
    # alter the alpha carbon coords
    alter *,b=0
    alter n. CA, b=f.pop(0)
    
    # color the surface
    e2s 1cll
    

You can make fake data in a shell by quickly executing: 
    
    
    # in a BASH shell, type the following
    
    for ((i=0; i<144; i++)); do echo $RANDOM.$RANDOM >> /tmp/bb; done
    

# The Code
    
    
    ##########################################################
    #
    # Expand_to_surface (e2s): Expands alpha-carbon-only 
    # b-factor coloring to an entire surface of the protein
    # of your choice.
    #
    # AUTHOR: Jason Vertrees -- Python code; Warren DeLano,
    #         the original code.
    #
    # COPYRIGHT: BSDL.  Feel free to use it.
    #
    # DATE: 2007-12-03
    #
    ##########################################################
    from pymol import cmd
    def e2s(sel):
            """
            e2s: Expand to surface
    
            e2s will color a surface based upon the b-factors
            from the alpha carbons in your selection.
    
            usage: e2s protName
            """
    
            cmd.create("ca_obj", sel + " and n. CA" )
            cmd.ramp_new("ramp_obj", "ca_obj", [0, 10], [-1, -1, 0] );
            cmd.set("surface_color", "ramp_obj", sel )
    
    cmd.extend("e2s", e2s);
    

# See Also

  * [alphaToAll](/index.php/AlphaToAll "AlphaToAll")
  * [surface](/index.php/Surface "Surface")
  * [isomesh](/index.php/Isomesh "Isomesh")
  * [isosurface](/index.php/Isosurface "Isosurface")
  * [ramp_new](/index.php/Ramp_new "Ramp new")
  * [map_new](/index.php/Map_new "Map new")



Retrieved from "[https://pymolwiki.org/index.php?title=Expand_To_Surface&oldid=8335](https://pymolwiki.org/index.php?title=Expand_To_Surface&oldid=8335)"


---

## Extra fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

extra_fit aligns multiple objects to a reference object. 

_New in PyMOL 1.7.2_

It can use any of PyMOL's pairwise alignment methods ([align](/index.php/Align "Align"), [super](/index.php/Super "Super"), [cealign](/index.php/Cealign "Cealign"), [fit](/index.php/Fit "Fit")...). More precisely it can use any function which takes arguments **mobile** and **target** , so it will for example also work with [tmalign](/index.php/Tmalign "Tmalign"). Additional keyword arguments are passed to the used method, so you can for example adjust outlier cutoff or create an alignment object. 

There are several similar commands/scripts for this job, like "A > align > all to this" from the PyMOL panel, the "[alignto](/index.php?title=Alignto&action=edit&redlink=1 "Alignto \(page does not exist\)")" command, [align_all.py and super_all.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/) scripts from Robert Campbell. 

## Usage
    
    
    extra_fit [ selection [, reference [, method ]]]
    

## Example

This will align 4 structures on CA atoms using the [super](/index.php/Super "Super") method. It will also create an alignment object. 
    
    
    fetch 1e4y 1ake 4ake 3hpq, async=0
    remove not chain A
    
    extra_fit name CA, 1ake, super, object=aln_super
    

Same, but with tmalign method (see [TMalign](/index.php/TMalign "TMalign")) 
    
    
    import tmalign
    extra_fit name CA, 1ake, tmalign, object=aln_tmalign
    

## See Also

  * [align_all.py and super_all.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/)



Retrieved from "[https://pymolwiki.org/index.php?title=Extra_fit&oldid=12909](https://pymolwiki.org/index.php?title=Extra_fit&oldid=12909)"


---

## FetchLocal

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Note
  * 3 Syntax
  * 4 The Code
  * 5 See Also



## Overview

Try fetching local copy of PDB files before going on internet to fetch a remote copy. 

Sometimes you have a local copy of the PDB structure files, or [biological unit](http://www.rcsb.org/pdb/static.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html) files. They are, however, usually organized into sub-directories (for a good reason). Setting [Fetch_Path](/index.php/Fetch_Path "Fetch Path") to the top directory will not make PyMOL to use local files before go on internet to fetch a remote copy. This script extends the functionality of the PyMOL command [Fetch](/index.php/Fetch "Fetch") by searching local sub-directories for target files first. 

The search order of the script: 1. local copy of PDB; 2. [Fetch_Path](/index.php/Fetch_Path "Fetch Path"); 3. remote repository ([Fetch_Host](/index.php/Fetch_Host "Fetch Host")). 

## Note

The script only copes with `.pdb` and `.pdb1` (or `.pdb2`, `.pdb3`, ...) file types. 

## Syntax

See [Fetch](/index.php/Fetch "Fetch")

## The Code
    
    
    """ 2012_09_18, Hongbo Zhu <hongbo.zhu.cn gmail>
        DESCRIPTION: Look for pdb* on local disk before fetching remotely.
    """
    import os
    import cmd
    
    #####################
    # user configuration
    localpdbdir = '/your/dir/to/PDB/divided'
    localbudir  = '/your/dir/to/PDB_BU/divided' # set to '' (empty string) if N/A
    
    def fetchlocal(code, name='', state=0, finish=1, discrete=-1,
                   multiplex=-2, zoom=-1, type='pdb', async=-1, path=None,
                   file=None, quiet=1):  
        """ Default parameters are the same as in function definition for fetch() 
            in file ./modules/pymol/importing.py
        """
        if type.startswith('pdb'): # pdb files, or biological unit files .pdb1 .pdb2 etc.
            if type == 'pdb': # pdb files
                localdir = localpdbdir
            else:
                localdir = localbudir
            for c in string.split(str(code)):
                subdir = c.lower()[1:3]
                if type == 'pdb': 
                    if os.path.isfile('%s/%s/pdb%s.ent.gz' % (localdir,subdir,c.lower())):
                        localfn = '%s/%s/pdb%s.ent.gz' % (localdir,subdir,c.lower())
                    elif os.path.isfile('%s/%s/pdb%s.ent' % (localdir,subdir,c.lower())):
                        localfn = '%s/%s/pdb%s.ent.gz' % (localdir,subdir,c.lower())
                    else:
                        localfn = None
                else:             
                    if os.path.isfile('%s/%s/%s.%s.gz' % (localdir,subdir,c.lower(),type)):
                        localfn =  '%s/%s/%s.%s.gz' % (localdir,subdir,c.lower(),type)
                    elif os.path.isfile('%s/%s/%s.%s' % (localdir,subdir,c.lower(),type)):
                        localfn = '%s/%s/%s.%s' % (localdir,subdir,c.lower(),type)
                    else:
                        localfn = None
    
                if localfn is not None:
                    print 'Load local file instead of fetching from internet: ', localfn 
                    cmd.load(filename=localfn, object=name, format='pdb', state=state, 
                             finish=finish, discrete=discrete, multiplex=multiplex, 
                             zoom=zoom, quiet=quiet)
                else:  # otherwise hand it over to pymol fetch function
                    cmd.fetch(c,name,state,finish,discrete,multiplex,zoom,type,async,path,file,quiet)
        else:
            cmd.fetch(code,name,state,finish,discrete,multiplex,zoom,type,async,path,file,quiet)
    
    
    cmd.extend('fetchlocal', fetchlocal)
    

# See Also

[Fetch](/index.php/Fetch "Fetch"), [Fetch_Path](/index.php/Fetch_Path "Fetch Path"), [Fetch_Host](/index.php/Fetch_Host "Fetch Host"), [Psico](/index.php/Psico "Psico")

Retrieved from "[https://pymolwiki.org/index.php?title=FetchLocal&oldid=12456](https://pymolwiki.org/index.php?title=FetchLocal&oldid=12456)"


---

## FilterByMol

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This script filters through all the PDBs in the parent dir (you can easily the the directory it scans). For each molecule, it saves **just** the ligands/heteroatoms (excluding the waters). This gives you a simple way to filter through a database of proteins looking only at their ligands. 

This script, as noted below, works on the objects at the level of a **molecule**. While we can [iterate](/index.php/Iterate "Iterate") over atom number (ID), residue number (resi), etc we do not have any such "MOLID". So, we provide this simple workaround. You might need this file because if you have a residue (like #111 from 3BEP) that consists of a molecule and an atom then there's no other way to save the separate pieces (of molecule/atom) into two (or more files). As you can see in the following listing, if we iterate over the hetero atoms (and not waters) in 3BEP we get, 
    
    
    PyMOL>iterate bymol het, print resi, resn, ID, chain, segi, alt
    111 5CY 6473 C  
    111 5CY 6474 C  
    111 5CY 6476 C  
    111 5CY 6477 C  
    111 5CY 6478 C  
    111 5CY 6479 C  
    111 5CY 6480 C  
    111 5CY 6481 C  
    111 5CY 6482 C  
    111 5CY 6483 C  
    111 5CY 6484 C  
    111 5CY 6485 C  
    111 5CY 6486 C  
    111 5CY 6487 C  
    111 5CY 6488 C  
    111 5CY 6489 C  
    111 5CY 6490 C
    

which does not allow us to separate the two pieces. 

## The Code
    
    
    python
    
    #
    # This simple script will filter through all PDBs in a directory, and for each one
    # save all the ligands/heterotoms (that aren't waters) to their own file.  This
    # script operates at the level of molecules, not residues, atoms, etc.  Thus, if
    # you have a ligand that PyMOL is treating as ONE residue, but is actually two
    # separate molecules, or a molecule and an atom, then you will get multiple files.
    #
    
    from glob import glob
    from os import path
    from pymol import stored
    
    theFiles = glob("../*.pdb");
    
    for f in theFiles:
        # load the file
        cmd.load(f);
        # remove the protein and waters
        cmd.remove("polymer or resn HOH");
    
        cmd.select("input", "all")
        cmd.select("processed", "none")
        mol_cnt = 0
    
        while cmd.count_atoms("input"):
            # filter through the selections, updating the lists
            cmd.select("current","bymolecule first input")
            cmd.select("processed","processed or current")
            cmd.select("input","input and not current")
    
            # prepare the output parameters
            curOut = path.basename(f).split(".")[0] + "_" + str(mol_cnt).zfill(5) + "_het.pdb"
            curSel = "current"
            
            # save the file
            cmd.save( curOut, curSel );
            print "Saved " + curSel + " to " + curOut
            
            mol_cnt = mol_cnt + 1;
    
        # remove all to move to next molecule
        cmd.delete("*");        
    
    python end
    

Retrieved from "[https://pymolwiki.org/index.php?title=FilterByMol&oldid=7665](https://pymolwiki.org/index.php?title=FilterByMol&oldid=7665)"


---

## Find buried waters

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script finds and turns a selection of waters determined to be buried (no solvent accessibility) in the original selection. If you prefer to define your own cutoff of "accessible" use the 2nd parameter. 

## Usage
    
    
    find_buried_waters [ sele [, cutoff [, state ]]]
    

## Script
    
    
    from pymol import cmd
    
    def find_buried_waters(sele='all', cutoff=-1, state=1, quiet=1, _self=cmd):
        '''
    DESCRIPTION
    
        Finds and turns a selection of waters determined to be buried (no solvent
        accessibility) in the original selection.
    
    ARGUMENTS
    
        sele = string: atom selection {default: all}
    
        cutoff = float: threshold on what one considers an "exposed"
        atom (in A**2) {default: surface_residue_cutoff}
    
        state = int: object state {default: 1}
        '''
        cutoff, state, quiet = float(cutoff), int(state), int(quiet)
    
        tmpObj=_self.get_unused_name("__tmp")
        _self.create(tmpObj, sele, state, 1, zoom=0)
    
        _self.set("dot_solvent", 1, tmpObj);
        _self.get_area(tmpObj, state=1, load_b=1)
    
        if cutoff < 0:
            cutoff = _self.get("surface_residue_cutoff")
        _self.remove(tmpObj + " and not solvent")
        _self.remove(tmpObj + " and b > %s" % cutoff)
    
        exposed = set()
        _self.iterate(tmpObj, "exposed.add((chain,resv))", space=locals())
    
        selName = _self.get_unused_name("buried")
        _self.select(selName, "(%s) in %s" % (sele, tmpObj))
    
        # clean up
        _self.delete(tmpObj)
    
        if not quiet:
            print ' Found %d buried water atoms' % (len(exposed))
    
        return sorted(exposed)
    
    cmd.extend('find_buried_waters', find_buried_waters)
    

## See Also

  * [get_area](/index.php/Get_Area "Get Area")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")



Retrieved from "[https://pymolwiki.org/index.php?title=Find_buried_waters&oldid=10740](https://pymolwiki.org/index.php?title=Find_buried_waters&oldid=10740)"


---

## Find symbol

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script looks for the given symbol in a given module(s). This can be handy if you forgot which module something is in or just want to query for a substring in a set of modules. 

# Usage
    
    
    # search pymol and pymol.cmd for any symbol starting with 'mov'
    
    fs("mov")
    ['cmd.get_movie_length',
     'cmd.get_movie_locked',
     'cmd.get_movie_playing',
     'cmd.mmove',
     'cmd.move',
     'cmd.movie',
     'cmd.moving',
     'cmd.remove',
     'cmd.remove_picked']
    
    # Search PyMOL's CMD module for something called align
    
    fs("align", "cmd")
    ['cmd.align', 'cmd.alignto', 'cmd.cealign', 'cmd.get_raw_alignment']
    

# The Code
    
    
    import pymol
    import inspect
    import pprint
    
    def fs(needle,haystack=["pymol","cmd"]):
        """
        This script will find the 'needle' in the 'haystack,' where the former is
        a symbol to be found in the latter, which is a module.
        """
        
        if type("") == type(haystack):
    
            haystack = [haystack,]
    
        for mod in haystack:
    
            found_list = map(lambda x: "%s.%s" % (mod,x), [name for name,obj in inspect.getmembers(eval(mod)) if needle in name])
    
        pprint.pprint(found_list)
    
        return found_list
    

Retrieved from "[https://pymolwiki.org/index.php?title=Find_symbol&oldid=11434](https://pymolwiki.org/index.php?title=Find_symbol&oldid=11434)"


---

## FindObjectsNearby

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This script returns the names of all objects (in state 1) that are within a certain cutoff distance from your provided target. For example, if you load an MD trajectory and are watching the binding pocket, you can select the ligand and then find out which loaded files are within 5 Ang, say, of that ligand atom. 

## The Code
    
    
    # -*- coding: utf-8 -*-
    #
    # findObjectsNearby.pml -- Get object names within a given radius of the target
    #
    def findObjectsNearby( target="", radius=2.0, doPrint=False, fileName="" ):
    	"""
    	DESCRIPTION:
    		finds all PyMOL object names within [radius] Angstroms of [target].
    	
    	PARAMETERS:
    		target,		the base selection/object around which to search for things
    		radius,		the radius of the sphere centered around target within which to search
    		doPrint,	boolean, if True print the results to console, False no printing
    		fileName,	string, if not blank the list called [fileName] is written to disk
    	
    	RETURNS:
    		[list] of results or None if bad input.
    		
    	NOTES:
    		* This function ALWAYS returns a list of results or None if the user submits malformed input.
    		* The user may opt to print the list to console and also save it to a file.
    
            EXAMPLE:
    		* Find all objects nearby 1j01 and heteroatoms that aren't water.  (Contrived example)
    		findObjectsNearby 1j01 and (het not resn HOH), 5.5, doPrint=True
    
    	AUTHOR:
    		Jason Vertrees, 2009.  Simplified with the help of Warren DeLano.
    	
    	"""
    	if ( len(target)==0 ):
    		print "Error: please provide a target selection."
    		return None
    	
    	stored.objs = {}
    	cmd.iterate_state(1, "(" + target + ") around " + str(radius), "stored.objs[model]=1" )
    
    	# save to file?
    	if len(fileName) != 0:
    		try:
    			outFile = open(fileName, 'wb')
    			for x in stored.objs.keys(): f.write( x + "\n" )
    			outFile.close()
    		except IOError:
    			print "Error: couldn't open/write to output file, ", fileName
    	# print?
    	if doPrint:
    		print stored.objs.keys()
    
    	return stored.objs.keys()
    
    cmd.extend("findObjectsNearby", findObjectsNearby)
    

## Troubleshooting

If the script isn't working like you think it should, then make sure that you're in **state #1**. 

Retrieved from "[https://pymolwiki.org/index.php?title=FindObjectsNearby&oldid=9957](https://pymolwiki.org/index.php?title=FindObjectsNearby&oldid=9957)"


---

## Findseq

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/findseq.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/findseq.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
# Overview & Motivation

Anyone ever give you a protein and then say, find the sequence "FLVEW"? Well, this script will find any string or regular expression in a given object and return that selection for you. Here's an example, 
    
    
    reinitialize
    import findseq
    # fetch two sugar-binding PDB
    fetch 1tvn, async=0
    # Now, find FATEW in 1tvn, similarly
    findseq FATEW, 1tvn
    # lower-case works, too
    findseq fatew, 1tvn
    # how about a regular expression?
    findseq F.*W, 1tvn
     
    # Find the regular expression:
    #  ..H[TA]LVWH
    # in the few proteins loaded.
    # I then showed them as sticks and colored them to highlight matched AAs
    for x in cmd.get_names(): findseq.findseq("..H[TA]LVWH", x, "sele_"+x, firstOnly=1)
    

[![](/images/c/c3/SeqFinder.png)](/index.php/File:SeqFinder.png)

[](/index.php/File:SeqFinder.png "Enlarge")

Red residues were those matching the regular expression '..H[TA]LVWH'.

# Usage

I built this to be rather flexible. You call it as: 
    
    
    findseq needle, haystack[, selName[, het[, firstOnly ]]]
    

where the options are: 

    

    **needle** the sequence of amino acids to find. Should be a string of one letter amino acid abbreviations. Can also be a string-style regular expression (eg. FW.*QQ).
    **haystack** the PyMOL object or selection in which to search
    **selName** the name of the returned selection. If you leave this blank, it'll be _foundSeqXYZ_ where XYZ is some random integer (eg. foundSeq1435); if you supply _sele_ then the usual PyMOL _(sele)_ is used; and, finally, if it's anything else, then that will be used verbatim. Defaults to _foundSeqXYZ_ so as not to overwrite any selections you might have in _sele_.
    **het** 0/1 -- if 0 then heteroatoms are not considered; if 1 then they are; defaults to 0
    **firstOnly** 0/1 -- if 0 then all matches are selected and returned; if 1 then only the first is returned

# See Also

select_pepseq@[Psico](/index.php/Psico "Psico")

Retrieved from "[https://pymolwiki.org/index.php?title=Findseq&oldid=13907](https://pymolwiki.org/index.php?title=Findseq&oldid=13907)"


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

## Flatten obj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/flatten_obj.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/flatten_obj.py)  
Author(s)  | [Spencer Bliven](/index.php/User:Sbliven "User:Sbliven")  
License  | Public Domain   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Arguments
  * 4 Notes
  * 5 Examples
  * 6 See Also



# Description

The `flatten_obj` python script combines multiple objects or states into a single object, renaming chains where required. 

This is particularly useful for dealing with biological assemblies, which are loaded as multi-state objects when fetched using `[fetch](/index.php/Fetch "Fetch") PDBID, type=pdb1`. It can also be used as a quick way to combine multiple objects without causing collisions between chain identifiers. 

The command re-letters chains to avoid collisions. Older versions of PyMOL restrict the chain id to a single character, so the script will fail for assemblies with >62 chains. With more recent versions, this problem is solved with multi-character chain IDs. Several options are available for how re-lettering should occur. 

# Usage
    
    
       flatten_obj name, selection[, state[, rename[, quiet[, chain_map]]]]
    

# Arguments

  * name = a unique name for the flattened object {default: flat}


  * selection = the set of objects to include in the flattening. The selection will be expanded to include all atoms of objects. {default: all}


  * state = the source state to select. Use 0 or -1 to flatten all states {default: 0}


  * rename = The scheme to use for renaming chains: {default: 0} 
    * (0) preserve chains IDs where possible, rename other chains alphabetically
    * (1) rename all chains alphabetically
    * (2) rename chains using the original chain letter, object name, and state


  * quiet = If set to 0, print some additional information about progress and chain renaming {default: 1}


  * chain_map = An attribute name for the 'stored' scratch object. If specified, `stored.<chain_map>` will be populated with a dictionary mapping the new chain names to a tuple giving the originated object, state, and chainID. {default: ""}



# Notes

Like the select command, if name is omitted then the default object name ("flat") is used as the name argument. 

Chain renaming is tricky. PDB files originally limited chains to single letter identifiers containing [A-Za-z0-9]. When this was found to be limiting, multi-letter chains (ideally < 4 chars) were allowed. This is supported as of PyMOL 1.7. Earlier versions do not accept rename=2, and will raise an exception when flattening a structure with more than 62 chains. 

# Examples
    
    
       flatten_obj flat, nmrObj
       flatten_obj ( obj1 or obj2 )
    

# See Also

  * [split_states](/index.php/Split_states "Split states")



Retrieved from "[https://pymolwiki.org/index.php?title=Flatten_obj&oldid=13954](https://pymolwiki.org/index.php?title=Flatten_obj&oldid=13954)"


---

## FocalBlur

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [focal_blur.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/focal_blur.py)  
Author(s)  | [Jarl Underhaug](/index.php?title=User:Jarl.Underhaug&action=edit&redlink=1 "User:Jarl.Underhaug \(page does not exist\)")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![StylizedFocalBlur.png](/images/a/af/StylizedFocalBlur.png)](/index.php/File:StylizedFocalBlur.png)

  


## Contents

  * 1 Description
    * 1.1 Usage
    * 1.2 Notes
    * 1.3 Bugs
    * 1.4 Examples
  * 2 Script



## Description

This script creates fancy figures by introducing a focal blur to the image. The object at the origin will be in focus. 

### Usage

Load the script using the [run](/index.php/Run "Run") command. Execute the script using PyMOL syntax: 
    
    
    FocalBlur aperture=2.0,samples=20,ray=1
    

or using python syntax: 
    
    
    FocalBlur(aperture=2.0,samples=20,ray=1)
    

  
For additional options, see the script comments. 

### Notes

  * When using raytracing, the image creation will take _n_ times longer than normal, where _n_ is the number of samples.
  * The aperture is related to the aperture on a camera, and is inversely proportional to the f-number.



### Bugs

  * FocalBlur uses the Python Image Library (PIL), a necessary components of PIL is missing in the Windows version of PyMOL
  * There is a bug when not using ray tracing with the free version of PyMOL



### Examples

  * [![FocalBlur aperture=1,samples=100,ray=1](/images/4/47/FocalBlur_a1.0_r1.png)](/index.php/File:FocalBlur_a1.0_r1.png "FocalBlur aperture=1,samples=100,ray=1")

FocalBlur aperture=1,samples=100,ray=1 

  * [![FocalBlur aperture=2,samples=100,ray=1](/images/f/f0/FocalBlur_a2.0_r1.png)](/index.php/File:FocalBlur_a2.0_r1.png "FocalBlur aperture=2,samples=100,ray=1")

FocalBlur aperture=2,samples=100,ray=1 

  * [![FocalBlur aperture=4,samples=400,ray=1](/images/5/5a/FocalBlur_a4.0_r1.png)](/index.php/File:FocalBlur_a4.0_r1.png "FocalBlur aperture=4,samples=400,ray=1")

FocalBlur aperture=4,samples=400,ray=1 

  * [![FocalBlur aperture=4,samples=400,ray=0](/images/e/e8/FocalBlur_a4.0_r0.png)](/index.php/File:FocalBlur_a4.0_r0.png "FocalBlur aperture=4,samples=400,ray=0")

FocalBlur aperture=4,samples=400,ray=0 

  * [![Focal blur ex6.png](/images/7/7a/Focal_blur_ex6.png)](/index.php/File:Focal_blur_ex6.png)

  * [![Focal blur ex ap3.png](/images/d/d6/Focal_blur_ex_ap3.png)](/index.php/File:Focal_blur_ex_ap3.png)

  * [![Focal blur ex ap3 mode1.png](/images/1/1e/Focal_blur_ex_ap3_mode1.png)](/index.php/File:Focal_blur_ex_ap3_mode1.png)

  * [![Focal blur ex ap3 mode2.png](/images/8/86/Focal_blur_ex_ap3_mode2.png)](/index.php/File:Focal_blur_ex_ap3_mode2.png)

  * [![Focal blur ex ap3 mode3.png](/images/f/fd/Focal_blur_ex_ap3_mode3.png)](/index.php/File:Focal_blur_ex_ap3_mode3.png)




## Script

Load the script using the [run](/index.php/Run "Run") command 
    
    
    from pymol import cmd
    from tempfile import mkdtemp
    from shutil import rmtree
    from math import sin,cos,pi,sqrt
    from PIL import Image
     
    def FocalBlur(aperture=2.0,samples=10,ray=0,width=0,height=0):
        '''
    DESCRIPTION
     
        Creates fancy figures by introducing a focal blur to the image. The object
        at the origin will be in focus. 
     
    AUTHOR
     
        Jarl Underhaug
        University of Bergen
        jarl_dot_underhaug_at_gmail_dot_com
    
        Updates by Jason Vertrees and Thomas Holder
     
    USAGE
     
        FocalBlur aperture=float, samples=int, ray=0/1, width=int, height=int
     
    EXAMPELS
     
        FocalBlur aperture=1, samples=100
        FocalBlur aperture=2, samples=100, ray=1, width=600, height=400
        '''
    
        # Formalize the parameter types
        ray = (ray in ("True", "true", 1, "1"))
        aperture, samples = float(aperture), int(samples)
        width, height = int(width), int(height)
     
        # Create a temporary directory
        tmpdir = mkdtemp()
    
        # Get the orientation of the protein and the light
        light = cmd.get('light')[1:-1]
        light = [float(s) for s in light.split(',')]
        view = cmd.get_view()
     
        # Rotate the protein and the light in order to create the blur
        for frame in range(samples):
            # Angles to rotate protein and light
            # Populate angles as Fermat's spiral
            theta = frame * pi * 110.0/144.0
            radius = 0.5 * aperture * sqrt(frame/float(samples-1))
            x = cos(theta) * radius
            y = sin(theta) * radius
            xr = x/180.0*pi
            yr = y/180.0*pi
     
            # Rotate the protein
            cmd.turn('x',x)
            cmd.turn('y',y)
     
            # Rotate the light
            ly = light[1]*cos(xr)-light[2]*sin(xr)
            lz = light[2]*cos(xr)+light[1]*sin(xr)
            lx = light[0]*cos(yr)+lz*sin(yr)
            lz = lz*cos(yr)-lx*sin(yr)
            cmd.set('light',[lx,ly,lz])
     
            curFile = "%s/frame-%04d.png" % (tmpdir,frame)
            print "Created frame %i/%i (%0.0f%%)" % (frame+1,samples,100*(frame+1)/samples)
    
            # Save the image to temporary directory
    	if ray:
                    cmd.ray(width,height)
                    cmd.png(curFile)
    	else:
            	cmd.png(curFile,quiet=1)
            
            # Create the average/blured image
            try:
                avg = Image.blend(avg,Image.open(curFile),1.0/(frame+1))
            except:
                avg = Image.open(curFile)
            
            # Return the protein and the light to the original orientation
            cmd.set('light',light)
            cmd.set_view(view)
     
        # Load the blured image
        avg.save('%s/avg.png' % (tmpdir))
        cmd.load('%s/avg.png' % (tmpdir))
     
        # Delete the temporary files
        rmtree(tmpdir)
    
    cmd.extend('FocalBlur', FocalBlur)
    

Retrieved from "[https://pymolwiki.org/index.php?title=FocalBlur&oldid=11652](https://pymolwiki.org/index.php?title=FocalBlur&oldid=11652)"


---

## Format bonds

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/format_bonds.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/format_bonds.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The script **format_bonds** will automatically format bonds in amino acids. 

## Contents

  * 1 Usage
  * 2 Examples
  * 3 Notes
  * 4 SEE ALSO



## Usage
    
    
    format_bonds [ selection [, bonds ]]
    

## Examples

  * [![format_bonds bonds=1](/images/1/10/PHE_valence_0.png)](/index.php/File:PHE_valence_0.png "format_bonds bonds=1")

format_bonds bonds=1 

  * [![format_bonds bonds=2](/images/9/9f/PHE_valence_1_mode_1.png)](/index.php/File:PHE_valence_1_mode_1.png "format_bonds bonds=2")

format_bonds bonds=2 

  * [![format_bonds](/images/c/ce/PHE_delocalized.png)](/index.php/File:PHE_delocalized.png "format_bonds")

format_bonds 



    
    
    import format_bonds
    
    frag PHE
    format_bonds
    
    format_bonds bonds=2
    

  


## Notes

  * Remember to correctly configure plugin import (see: [Git intro](/index.php/Git_intro "Git intro"))
  * **format_bonds** will introduce delocalized bonds by default or when _bonds_ is larger than 2.  

  * Setting _bonds=1_ will simply disable valence display (globally!)
  * The _selection_ argument is 'all' by default and can be used to restrict editing to selected residues.  

  * Note that **format_bonds** will also format acidic residues, the C-terminus, arginine and nitro groups.
  * Tip: press the TAB key after entering format_bonds to get argument suggestions



  


# SEE ALSO

[Git intro](/index.php/Git_intro "Git intro"), [Valence](/index.php/Valence "Valence"), [Modeling and Editing Structures](/index.php/Modeling_and_Editing_Structures "Modeling and Editing Structures")

Retrieved from "[https://pymolwiki.org/index.php?title=Format_bonds&oldid=13941](https://pymolwiki.org/index.php?title=Format_bonds&oldid=13941)"


---

## Forster distance calculator

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/forster_distance_calculator.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/forster_distance_calculator.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Spectre input
  * 3 Getting spectres
  * 4 Input values
  * 5 How the script works
  * 6 How to run the script



## Introduction

This script can be handsome, if one is working with Förster resonance energy transfer in proteins. The script calculates the [Förster distance](http://en.wikipedia.org/wiki/F%C3%B6rster_resonance_energy_transfer): R0, from two dyes excitation/emission spectres.   
This script is very handsome, if one want's to pick two dyes from different companies, and the Förster distance is not provided in a table. 

This script does no calculation of proteins in pymol, but is made as a "hack/shortcut" for a python script.   
We use the python part of pymol to do the calculations, so a student would not need to install python at home, but simply pymol.   


## Spectre input

**There should be provided the path to four datafiles:**   
Excitation and Emission spectre for the Donor dye: D_Exi="path.txt" D_Emi="path.txt   
Excitation and Emission spectre for the Acceptor dye: A_Exi="path.txt" A_Emi="path.txt   


Each of the files shall be a two column file. Separated by space. Numbers are "." dot noted and not by "," comma. One can for example prepare the files by search-and-replace in a small text editor.   
The first column is the wavelength in nanometers "nm" or centimetres "cm". Second column is the arbitrary units of excitation/emission. **For example:**
    
    
    Absorption
    Wavelength "Alexa Fluor 488 H2O (AB)" 
    
    300.00 0.1949100000 
    301.00 0.1991200000 
    302.00 0.2045100000 
    303.00 0.2078800000 
    ....
    

## Getting spectres

For example provides the company [ATTO-TEC](http://www.atto-tec.com/index.php?id=128&L=1&language=en) spectre for their dyes in excel format [Spectre(excel-file)](https://www.atto-tec.com/index.php?id=113&no_cache=1&L=1). These can easily be copied from excel to a flat two column file. 

Some of the most cited dyes are the Alexa Fluor dyes from the company [invitrogen](http://www.invitrogen.com/site/us/en/home/References/Molecular-Probes-The-Handbook/Fluorophores-and-Their-Amine-Reactive-Derivatives/Alexa-Fluor-Dyes-Spanning-the-Visible-and-Infrared-Spectrum.html). Invitrogen does not provide spectres for their dyes in dataformat, but as flat picture files. 

Luckily, a group on University of Arizona have traced several spectre of dyes from the literature with a graphics program. They have made these spectre easily public at <http://www.spectra.arizona.edu/>. I highly recommend this homepage. With these Spectra, the script can calculate the Forster Distance for different dyes from different companies.   


Download one spectrum at the time by "deselecting" in the right side of the graph window. Then get the datafile with the button in the lower right corner. 

## Input values

The calculation needs input for: 

  * **Qd** : Fluorescence quantum yield of the donor in the absence of the acceptor. 

    This is normally provided from the manufacture, and this information is normally also available in the [database](http://www.spectra.arizona.edu/).
    It is recognized in the literature, that the value for **Qd** , changes on behalf on which position its located on the protein.
    Changes in 5-10 % is normal, but much larger values can be expected if the dye make hydrophobic interactions with the protein. This is dye-nature dependant.
    This value can only be determined precisely with experiments for the protein at hand.
    But anyway, for a "starting" selection for a FRET pair, the manufacture **Qd** should be "ok".
  * **Kappa2 = 2.0/3.0** : Dipole orientation factor. 

    In the literature, I often see this equal 2/3=0,667 This corresponds for free diffusing dye.
    For a dye attached to a protein, and restricted in one direction, the value is equal Kappa2=0.476
    But since we take the 6th root of the value, the difference ends up only being 5% in relative difference for the R0 value.
    Kappa2=2/3 can be a valid fast approximation, for the purpose of ordering dyes. But be careful and check the literature for discussions.
  * **n = 1.33** : 

    water=1.33, protein=1.4, n(2M GuHCl)=1.375
    Again, we take the sixth root, so the differences are "not so important".



## How the script works

  * The script integrate the Donor emission spectre, to divide with it. 

    This is done by simple numerical integration by the [Rectangle method](http://en.wikipedia.org/wiki/Rectangle_method)
  * All datapoints for the Acceptor excitation spectre are scaled, so that **e**(molar extinction coefficient) fits the wavelength where it is determined.
  * For the multiplication of the spectre, the script test that both data points exist before multiplication. 

    This can be troublesome, if there is 1-2 nm between each data point. Or if they are "un-even". Like Donor has 100.2 100.4 100.4 and Acceptor has 100.1 100.3 100.5
    But most spectre have much better resolution. Mostly 0.2 nm, and "even spaced"
    So this should work quite good as well.



The output is a text file, with all the datapoints. 

  1. Column: Donor: The input Emission wavelength.
  2. Column: Donor: The input Emission, arbitrary units of excitation.
  3. Column: Donor: The input Emission, arbitrary units of excitation, divided by the integrated area.
  4. Column: Acceptor: The input excitation wavelength.
  5. Column: Acceptor: The input excitation **e**(molar extinction coefficient).
  6. Column: Acceptor: The input excitation **e**(molar extinction coefficient), scaled correctly.
  7. Column: The calculated overlap for this datapoint.
  8. Column: The summed values for calculated overlap, until this datapoint.



  
It also automatically generates a [gnuplot](http://www.gnuplot.info/) .plt file, for super fast making graphs of the spectres. In linux, gnuplot is normally part of the installation. Just write: gnuplot FILENAME.plt If you are in windows, [download gnuplot](http://sourceforge.net/projects/gnuplot/files/) and open the .plt file. 

Gnuplot makes three graphs in .eps format. They can be converted to pdf with the program: epstopdf or eps2pdf. They are for example part of LaTeX: C:\Program Files (x86)\MiKTeX 2.9\miktex\bin or you can [download it here](http://tinyurl.com/eps2pdf). The format of .eps is choses, since gnuplot then allows for math symbols in the graphs. 

  * [![All graphs plotted together, but not scaled correctly.](/images/7/77/1-ALEXA488Emi-ALEXA633Exi-overlap-all-spectre.png)](/index.php/File:1-ALEXA488Emi-ALEXA633Exi-overlap-all-spectre.png "All graphs plotted together, but not scaled correctly.")

All graphs plotted together, but not scaled correctly. 

  * [![The Donor emission \(y1\) and Acceptor excitation \(y2\) scaled correctly.](/images/7/7d/2-ALEXA488Emi-ALEXA633Exi-overlap-normalized-spectre.png)](/index.php/File:2-ALEXA488Emi-ALEXA633Exi-overlap-normalized-spectre.png "The Donor emission \(y1\) and Acceptor excitation \(y2\) scaled correctly.")

The Donor emission (y1) and Acceptor excitation (y2) scaled correctly. 

  * [![The Rectangle integral datapoint \(y1\) and the summed integration \(y2\).](/images/0/09/3-ALEXA488Emi-ALEXA633Exi-overlap-integral.png)](/index.php/File:3-ALEXA488Emi-ALEXA633Exi-overlap-integral.png "The Rectangle integral datapoint \(y1\) and the summed integration \(y2\).")

The Rectangle integral datapoint (y1) and the summed integration (y2). 




## How to run the script

Make a pymol .pml file like this. Write in the required information, and execute/run the script with pymol. Then open the .plt file afterwards with gnuplot. 
    
    
    ## Change to your directory
    cd /homes/YOU/Documents/Atto-dyes/Spectre/ALEXA488-ALEXA633
    import forster_distance_calculator
    forster D_Exi=ALEXA488Exi.txt, D_Emi=ALEXA488Emi.txt, A_Exi=ALEXA633Exi.txt, A_Emi=ALEXA633Emi.txt, A_e_Max_Y=159000, A_e_Max_X=621, Qd=0.92
    

Retrieved from "[https://pymolwiki.org/index.php?title=Forster_distance_calculator&oldid=13908](https://pymolwiki.org/index.php?title=Forster_distance_calculator&oldid=13908)"


---

## Frame slider

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/frame_slider.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/frame_slider.py)  
Author(s)  | Matthew Baumgartner   
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  
The Frame slider plugin provides a simple slider bar that allows you to quickly skip through frames in PyMOL. It also has a text field that you can type the desired frame in and it will automatically go to that frame. 

[![Frame slider.png](/images/a/ad/Frame_slider.png)](/index.php/File:Frame_slider.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Frame_slider&oldid=13909](https://pymolwiki.org/index.php?title=Frame_slider&oldid=13909)"


---

## Get colors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/get_colors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/get_colors.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
  * **get_colors** contains two functions that can be useful when working with colors 
    * _get_colors_ : returns all available PyMOL colors
    * _get_random_color_ : returns a random available PyMOL color


  * Note that basic colors can be accessed manually without this script from the PyMOL menu under **Setting** \--> **Colors...**
  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



[![1LSD colored randomly by residue](/images/1/13/1LSD_random_colors.png)](/index.php/File:1LSD_random_colors.png "1LSD colored randomly by residue")

  


## Usage
    
    
    get_colors [ selection [, quiet ]]
    get_random_color [ selection [, quiet ]]
    

## Examples
    
    
    # basic example
    get_colors # basic colors
    get colors all # larger range with intermediates
    
    
    
    #get disco funky
    import get_colors
    from get_colors import get_colors
    from get_colors import get_random_color
    
    cmd.delete('all')
    cmd.fetch('1LSD', async=0) # :P
    cmd.hide('everything')
    cmd.show_as('sticks','not hetatm')
    cmd.orient()
    
    python # start a python block
    from pymol import stored
    stored.atom_list=[]
    cmd.iterate('name CA','stored.atom_list.append([model, resi])')
    resi_list=["model %s and resi %s"%(value[0],value[1]) for value in stored.atom_list]
    for resi in resi_list: cmd.color(get_random_color(),resi)
    python end # end python block
    cmd.ray()
    

|   
---|---  
  
  


## SEE ALSO

  * [Color](/index.php/Color "Color")
  * [Get Color Indices](/index.php/Get_Color_Indices "Get Color Indices")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_colors&oldid=13942](https://pymolwiki.org/index.php?title=Get_colors&oldid=13942)"


---

## Get Coordinates I

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

There are several ways to extract or load atomic coordinates in PyMOL using the python API. 

## Contents

  * 1 Extract coordinates using get_coords
  * 2 Extract coordinates using get_coordset
  * 3 Extract coordinates using get_model
  * 4 Extract coordinates using iterate_state
  * 5 Load coordinates using alter_state
  * 6 Load coordinates using update
  * 7 Load coordinates using load_coords
  * 8 Load coordinates using load_coordset
  * 9 Load coordinates using get_coordset



## Extract coordinates using [get_coords](/index.php/Get_coords "Get coords")

_New in PyMOL 1.7.4_. This is the fastest method to extract coordinates from a selection. It considers the object rotation matrix. 
    
    
    xyz = cmd.get_coords('sele', 1)
    

## Extract coordinates using [get_coordset](/index.php/Get_coordset "Get coordset")

_New in PyMOL 1.7.4_. Operates on the object-state level, not on selections. Does **not** consider the object rotation matrix. Retrieves coordinates in original order (e.g. PDB file atom order), not in sorted atom order. Faster than [get_coords](/index.php/Get_coords "Get coords"). 
    
    
    xyz = cmd.get_coordset(objectname, 1)
    

## Extract coordinates using [get_model](/index.php/Get_Model "Get Model")

Before 1.7.4, this was the fastest method to extract coordinates. It considers the object rotation matrix. 
    
    
    xyz = cmd.get_model('sele', 1).get_coord_list()
    

## Extract coordinates using [iterate_state](/index.php/Iterate_state "Iterate state")

This is much slower than the first method. It does **not** consider the object rotation matrix. 
    
    
    xyz = []
    cmd.iterate_state(1, 'sele', 'xyz.append([x,y,z])', space=locals(), atomic=0)
    

## Load coordinates using [alter_state](/index.php/Alter_State "Alter State")

This is the most convenient way to load coordinates and works equivalent to **iterate_state**. 
    
    
    xyz = [...] # some Nx3 list with coordinates
    xyz_iter = iter(xyz)
    cmd.alter_state(1, 'sele', '(x,y,z) = xyz_iter.next()', space=locals())
    

## Load coordinates using [update](/index.php/Update "Update")

This example gets a copy of the coordinates in Python, rotates the object about the Z axis, and then updates the coordinates in the original object. 
    
    
    model = cmd.get_model('pept')
    for a in model.atom:
        a.coord = [ -a.coord[1], a.coord[0], a.coord[2]]
    
    cmd.load_model(model, "_tmp")
    cmd.update("pept", "_tmp")
    cmd.delete("_tmp")
    

## Load coordinates using [load_coords](/index.php/Load_coords "Load coords")

_Changed in PyMOL 1.7.3_. Update selection coordinates from a Nx3 float array. 
    
    
    cmd.load_coords(xyz, selection)
    

## Load coordinates using [load_coordset](/index.php?title=Load_coordset&action=edit&redlink=1 "Load coordset \(page does not exist\)")

_New in PyMOL 1.7.4_. Update object state coordinates from a Nx3 float array. Can also append a state. Order of coordinates equivalent to [get_coordset](/index.php/Get_coordset "Get coordset"). 
    
    
    cmd.load_coordset(xyz, objectname, 1)
    

## Load coordinates using [get_coordset](/index.php/Get_coordset "Get coordset")

_New in PyMOL 1.7.4_. The [get_coordset](/index.php/Get_coordset "Get coordset") function can also return a memory view instead of a copy. This allows modifying the coordinates in place. 
    
    
    xyz = cmd.get_coordset(objectname, 1, copy=0)
    xyz[:] = newxyz
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Coordinates_I&oldid=12261](https://pymolwiki.org/index.php?title=Get_Coordinates_I&oldid=12261)"


---

## Get Coordinates II

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    from pymol import cmd
    from pymol import stored
    
    stored.xyz = []
    cmd.iterate_state(1,"pept","stored.xyz.append([x,y,z])")
    
    # at this point, stored.xyz is a native Python array holding
    # the coordinates, which you can modify as required
    
    stored.xyz = map(lambda v:[-v[1],v[0],v[2]],stored.xyz)
    
    # and now you can update the internal coordinate sets
    
    cmd.alter_state(1,"pept","(x,y,z)=stored.xyz.pop(0)")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_Coordinates_II&oldid=6350](https://pymolwiki.org/index.php?title=Get_Coordinates_II&oldid=6350)"


---

## Get raw distances

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/get_raw_distances.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/get_raw_distances.py)  
Author(s)  | [Takanori Nakane](/index.php?title=User:TakanoriNakane&action=edit&redlink=1 "User:TakanoriNakane \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.querying](http://pymol.org/psicoredirect.php?psico.querying)  
  
get_raw_distances can dump distance objects, created with [distance](/index.php/Distance "Distance"). 

This script also provides the command **select_distances** , which selects atoms from distance objects. 

_Warning: the atoms are hashed by coordinates; this could cause issues if coordinates are altered after distance objects have been created (see also[dynamic_measures](/index.php?title=Dynamic_measures&action=edit&redlink=1 "Dynamic measures \(page does not exist\)") setting)._

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Example Atom Selection
  * 5 See Also



## Usage
    
    
    get_raw_distances [ names [, state [, selection ]]]
    
    
    
    select_distances [ names [, name [, state [, selection [, cutoff ]]]]]
    

## Arguments

  * **names** = string: names of distance objects (no wildcards!) {default: all measurement objects}
  * **state** = integer: object state {default: 1}
  * **selection** = string: atom selection {default: all}



## Example
    
    
    fetch 2xwu, async=0
    
    # interface polar contacts
    distance iface_hbonds, chain A, chain B, mode=2
    
    # dump (model,index) information
    get_raw_distances iface_hbonds
    

## Example Atom Selection
    
    
    # select atoms
    select_distances iface_hbonds, hbsele
    
    # nice representation
    set cartoon_side_chain_helper
    set sphere_scale, 0.5
    as cartoon, 2xwu
    show sticks, byres hbsele
    show spheres, hbsele
    

## See Also

  * [distance](/index.php/Distance "Distance")
  * [dynamic_measures](/index.php?title=Dynamic_measures&action=edit&redlink=1 "Dynamic measures \(page does not exist\)")
  * [find_pairs](/index.php/Find_pairs "Find pairs")
  * [get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_raw_distances&oldid=13910](https://pymolwiki.org/index.php?title=Get_raw_distances&oldid=13910)"


---

## GetNamesInSel

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script returns the list of object names that exist in any given selection. For example, if you have 100 objects loaded (viz. ligand poses from small molecule docking) and select all atoms w/in some cutoff, then you have a selection. What objects are actually in that selection? This script will tell you. 

# The Code
    
    
    from pymol import cmd, stored
    def getNamesInSel(sel, sorted=0, reversed=0, quiet=1):
        """
        PARAMETERS
            sel,
                The selection, object or group to iterate over
            sorted (boolean),
                Should the list be sorted?
            reversed (boolean)
                Should the list be reversed before returned?  (Combined
                with the above, you can return a decreasing-sorted list
                of names
     
        RETURNS
            list[] of strings, representing the object names desired.
    
        """
        rList = cmd.get_object_list('(' + sel + ')')
    
        # if you want the list reversed or sorted,
        # uncomment the following lines
        if int(sorted):
            rList.sort()
        if int(reversed):
            rList.reverse()
    
        if not int(quiet):
            print ' getNamesInSel: ', rList
    
        return rList
    
    cmd.extend("getNamesInSel", getNamesInSel)
    

# See Also

  * [Get_Names](/index.php/Get_Names "Get Names")
  * [get_object_list](/index.php/Get_object_list "Get object list")



Retrieved from "[https://pymolwiki.org/index.php?title=GetNamesInSel&oldid=9417](https://pymolwiki.org/index.php?title=GetNamesInSel&oldid=9417)"


---

## Grepsel

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Create selections matching motifs, using python regular expression syntax. 

This is very similar to the [FindSeq](/index.php/FindSeq "FindSeq") script. 

## The Code
    
    
    #Create named selections using regular expressions for the protein sequence
    
    import pymol
    import re
    cmd = pymol.cmd
    
    aa = { 'ASP' : 'D' , 'GLU' : 'E' , 'GLN' : 'Q' , 'ASN' : 'N' , 'SER' : 'S' ,
           'THR' : 'T' , 'CYS' : 'C' , 'HIS' : 'H' , 'ARG' : 'R' , 'LYS' : 'K' ,
           'MET' : 'M' , 'ALA' : 'A' , 'ILE' : 'I' , 'LEU' : 'L' , 'VAL' : 'V' ,
           'GLY' : 'G' , 'PRO' : 'P' , 'TRP' : 'W' , 'PHE' : 'F' , 'TYR' : 'Y' ,
           'SCY' : 'U' , 'ASX' : 'B' , 'GLX' : 'Z' , 'XXX' : 'X'}
    
    #made this before the sequence view option, probably another way to do it now
    
    def seqoneint(model):
       pymol.stored.seq = []
       cmd.iterate("%s and name ca"%model,"stored.seq.append(resn)")
       seq = ""
       for x in pymol.stored.seq:
          if aa.has_key(x):
             res = aa[x]
             seq = seq+res
          else:
             seq = seq + '-'
       return seq
    
    
    
    def grepsel(selection="(all)",stretch="",prefix="",combined="0",single="1"):
       '''
    DESCRIPTION
    
        Create selections matching motifs, using python regular expression syntax.
        Motif is automatically converted to uppercase. Motif selections are labelled
        as "prefix_motif_###", where ### is the index for the first residue of the
        match. Prefix defaults to selection name. combined = 1 creates one selection
        for all occurences. single = 1 creates one selection for each occurance
        (the default).
        
    USAGE
    
        grepsel selection, motif, [prefix, [combined, [single ]]]
    
    EXAMPLES
    
        Create selections for all motifs matching "ESS" (selection_ESS_###,...):
        grepsel selection, ess
    
        Create selections for the PxGY motif with prefix m (m_P.CY_###,...):
        grepsel selection, p.gy, m
        '''
     
       if selection == "(all)":
          selection = "all"
       if prefix == "":
          prefix=selection
    
       stretch = stretch.upper() 
       seq = seqoneint(selection)
       pymol.stored.resi = []
       pymol.stored.chain = []
       cmd.iterate("%s and name ca"%selection,"stored.resi.append(resi);stored.chain.append(chain)")
       motif = re.compile(stretch)
       occurrences = motif.finditer(seq)
       stretchmod = stretch.replace("+","\+")
       stretchmod = stretchmod.replace("?","\?")
    
       print stretchmod
       if combined == "1":
          cmd.select("%s_%s"%(prefix,stretch), "none")
    
    
       for find in occurrences:      
    
          mb = pymol.stored.resi[find.start()]
          me = pymol.stored.resi[find.end()-1]
    
          ch = pymol.stored.chain[find.start()]
          cmd.select("%s_%s_%s%s"%(prefix,stretch,me,ch), "chain %s and (i; %s-%s)"%(ch,int(mb),int(me)))
          if combined == "1":
             cmd.select("%s_%s"%(prefix,stretch),"\"%s_%s\" | (%s and chain %s and (i; %s-%s))"%(prefix,stretchmod,selection,ch,int(mb),int(me)))
    
       cmd.select("none")
       cmd.delete("sel*")
       
    cmd.extend("grepsel",grepsel)
    

## See Also

  * [FindSeq](/index.php/FindSeq "FindSeq")
  * [pepseq](/index.php/Selection_Algebra "Selection Algebra") selection operator
  * [seq_select](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/) by Robert Campbell



Retrieved from "[https://pymolwiki.org/index.php?title=Grepsel&oldid=9486](https://pymolwiki.org/index.php?title=Grepsel&oldid=9486)"


---

## Grepset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/grepset.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/grepset.py)  
Author(s)  | [Ezequiel Panepucci](/index.php?title=User:Zac&action=edit&redlink=1 "User:Zac \(page does not exist\)")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.helping](http://pymol.org/psicoredirect.php?psico.helping)  
  
Use this little script to explore PyMOL's myriad settings. 

Usefull for newbies and those with not so good memory skills... 

To use: 

  1. put the script in a file called **grepset.py**
  2. from within PyMOL execute **run grepset.py**. If you have the Pymol-script library configured, then **import grepset**
  3. try it out, see examples below.



Example 1: **grepset light**
    
    
    cartoon_highlight_color        default
    dot_lighting                   on
    light                          [ -0.40000, -0.40000, -1.00000 ]
    mesh_lighting                  off
    two_sided_lighting             off
    5 settings matched
    

Example 2: **grepset trans**
    
    
    cartoon_transparency           0.00000
    ray_transparency_contrast      1.00000
    ray_transparency_shadows       1
    ray_transparency_spec_cut      0.90000
    ray_transparency_specular      0.40000
    sphere_transparency            0.00000
    stick_transparency             0.00000
    transparency                   0.00000
    transparency_mode              2
    transparency_picking_mode      2
    10 settings matched
    

Example 3: **grepset ^trans**
    
    
    transparency                   0.00000
    transparency_mode              2
    transparency_picking_mode      2
    3 settings matched
    

Retrieved from "[https://pymolwiki.org/index.php?title=Grepset&oldid=13911](https://pymolwiki.org/index.php?title=Grepset&oldid=13911)"


---

## Hbplus

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [scripts/hbplus.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/hbplus.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
hbplus is a PyMOL wrapper for the [HBPLUS](http://www.biochem.ucl.ac.uk/bsm/hbplus/home.html) program. It creates distance objects for the hydrogen bonds found by HBPLUS. 

## Usage
    
    
    hbplus [ selection [, exe [, prefix [, state ]]]]
    

## See Also

  * [distance](/index.php/Distance "Distance")



Retrieved from "[https://pymolwiki.org/index.php?title=Hbplus&oldid=13912](https://pymolwiki.org/index.php?title=Hbplus&oldid=13912)"


---

## Colorama

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/colorama.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/colorama.py)  
Author(s)  | [Gregor Hagelueken](/index.php?title=User:Gha&action=edit&redlink=1 "User:Gha \(page does not exist\)")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 COLORAMA
    * 1.1 Program features
    * 1.2 Screenshot
    * 1.3 Usage
      * 1.3.1 Single color
      * 1.3.2 Color gradients
    * 1.4 Contact



# COLORAMA

COLORAMA is a PyMOL plugin which allows to color objects using adjustable scale bars. 

## Program features

In RGB mode, each color R, G, B is represented by one scale bar which can be manually adjusted while the selected object is colored in real time. For convenience, it is as well possible to switch to the HSV color system. 

Additionally, a color gradient with user-defined start- and end-colors can be created for the selected molecule. 

## Screenshot

[![COLORAMA-screenshot.jpg](/images/9/95/COLORAMA-screenshot.jpg)](/index.php/File:COLORAMA-screenshot.jpg)

## Usage

This plugin is included in the project [ Pymol-script-repo](/index.php/Git_intro "Git intro"). 

Manual install the plugin by copying the code below into an empty text file (e.g. "colorama.py") located in the \Pymol\modules\pmg_tk\startup directory. After PyMOL has been started, the program can be launched from the PLUGINS menu. COLORAMA has not been tested with PyMOL versions older than 1.0. 

### Single color

  1. Type the name of a PyMOL object to be colored into the COLORAMA entry field.
  2. Push the "Set" button.
  3. The scales are adjusted to the current color which is additionally visualized in a field left to the scales.
  4. If one of the scales is moved, the color of the selected object will change in real-time.
  5. Pushing the RGB or HSV buttons on the left allows to switch between both color systems.



### Color gradients

  1. After an object has been selected, push the "G" button (gradient).
  2. Select the start color by pushing "C1" and adjusting it using the scales.
  3. Select the end color "C2" in the same way.
  4. To create the gradient, push "Set gradient".



A new object will be created which is called "dummy-OLD_OBJECT". The B-factor column of this object is overwritten and now contains the number of each residue. The original object is left unchanged. The gradient mode can be left by pushing "M" (monochrome). This part of the program uses a modified version of the [color_b script](https://github.com/zigeuner/robert_campbell_pymol_scripts/tree/master/work_pymol) by Robert L. Campbell & James Stroud. 

## Contact

Gregor Hagelueken, hagelueken'at'pc.uni-bonn.de 

Retrieved from "[https://pymolwiki.org/index.php?title=Colorama&oldid=13845](https://pymolwiki.org/index.php?title=Colorama&oldid=13845)"


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

## ImmersiveViz

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 ImmersiveViz
  * 2 Availability
    * 2.1 Additional Information
    * 2.2 PyMol Script
    * 2.3 Head Tracking



## ImmersiveViz

ImmersiveViz was developed as a component to monitor head tracking and rotate the molecule displayed in such a manner to provide an immersive experience. The ImmersiveViz (MolViz) software integrates two forms of head tracking: Wiimote _infrared (IR) based_ (active tracking) and _webcam based_ (passive tracking). By rotating the molecule in a direction opposite to the motion of the user's head we provide a 3D experience; to the user, it appears as if they are 'peeking' around the side of the object.  
  


Our system contains a head tracking thread which communicates the users position to a PyMol script via a socket. The PyMol script unpacks the message and updates the world respectively. All rotation is done around the virtual origin in PyMol and zoom is also considered.  
  


We are also developing a Wiimote-based interface whereby the Wii remote can be used as an high degree-of-freedom input device (i.e. a 3d mouse).  
**Comments, questions, and feedback** can be sent to [molviz (at) googlegroups (dot) com](http://groups.google.com/group/molviz)  
  


**Contributed by[Christian Muise](http://www.haz.ca/) and [Ryan Lilien](http://www.cs.toronto.edu/~lilien) at the University of Toronto.**  


## Availability

### Additional Information

  * **Project Information:** <http://molviz.cs.toronto.edu/molviz>
  * **Project Discussion:** <http://groups.google.com/group/molviz>



### PyMol Script

  * **Project Page:** <http://code.google.com/p/immersive-viz/>
  * **Source:** [here](http://code.google.com/p/immersive-viz/source/browse/trunk/MolViz.py) and [here](http://code.google.com/p/immersive-viz/source/browse)
  * **Instructions:** [here](http://code.google.com/p/immersive-viz/) and [here](http://code.google.com/p/immersive-viz/wiki/UserManual)



### Head Tracking

  * **Project Page:** <http://code.google.com/p/htdp/>
  * **Source:** <http://code.google.com/p/htdp/source/browse>
  * **Instructions:** <http://code.google.com/p/htdp/wiki/Users>
  * **[RUCAP UM-5](/index.php/RUCAP_UM-5 "RUCAP UM-5") tracker:** <http://rucap.ru/en/um5/about>



Retrieved from "[https://pymolwiki.org/index.php?title=ImmersiveViz&oldid=11233](https://pymolwiki.org/index.php?title=ImmersiveViz&oldid=11233)"


---

## Inertia tensor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/inertia_tensor.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/inertia_tensor.py)  
Author(s)  | [Mateusz Maciejewski](http://www.mattmaciejewski.com)  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![It.png](/images/d/dd/It.png)](/index.php/File:It.png)

This script will draw the eigenvectors of the inertia tensor of the selected part of the molecule. 

## Usage
    
    
    tensor selection, name, model
    

## Examples

To draw the inertia tensor of the first model in PDB file 4B2R do: 
    
    
    fetch 4b2r, async=0
    hide lines
    show cartoon
    run inertia_tensor.py
    tensor 4b2r & n. n+ca+c & i. 569-623, tensor_dom11, 1
    

## Reference

A figure generated using this script can be seen in the following reference: 

  * _Estimation of interdomain flexibility of N-terminus of factor H using residual dipolar couplings_. Maciejewski M, Tjandra N, Barlow PN. **Biochemistry**. 50(38) 2011, p. 8138-49, fig. 1 [doi:10.1021/bi200575b](http://dx.doi.org/10.1021/bi200575b)



Retrieved from "[https://pymolwiki.org/index.php?title=Inertia_tensor&oldid=13943](https://pymolwiki.org/index.php?title=Inertia_tensor&oldid=13943)"


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

## Jump

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [jump.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/jump.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Example
  * 4 PyMOL API
  * 5 See Also



### Description

In the case of movies or simulation trajectories, it is often useful to quickly jump between different frames but not necessarily in any order. This script allows the user to specify a list of frames and flipping through the set of frames is controlled by a user-defined keystroke (limited to those available through the [Set Key](/index.php/Set_Key "Set Key") command). This script is a nice tool for identifying different protein states from a simulation and can be used in conjunction with the [Modevectors](/index.php/Modevectors "Modevectors") script in a very useful way by first characterizing the directionality of the movements (using [Modevectors](/index.php/Modevectors "Modevectors")) and then visualizing the change directly using this script. 

Note: This script also requires the [Check Key](/index.php/Check_Key "Check Key") script! 

  


### Usage

load the script using the [run](/index.php/Run "Run") command 
    
    
    jump [frame1 [,...frameN [,sleep=0.1 [,x=1 [,key='pgdn']]]]]
    

### Example
    
    
    jump 1, 10, 100
    
    #This causes the movie to jump from frame 1 to frame 10 and finally to frame 100 each time the 'pgdn' key is pressed
    #Note that there is no limit to the number of frames that can be specified
    
    jump 1, 10, 100, 1000, key='F1'
    #This sets the hotkey to F1.  When pressed, it will jump through the frames 1, 10, 100, 1000 respectively
    
    jump 1, 10, 100, sleep=0.01
    #This does the same thing as the first example except that the pause between frames is 10x faster
    #One could consider the "sleep" option to be similar to the frame rate
    
    jump 1, 10, 100, x=10
    #This does the same thing as the first example except that each press of the hotkey will 
    #cause the script to run through the specified frames 10 times before stopping.
    #Note that you will not be able to control or exit the script during the 10 iterations!
    

### PyMOL API
    
    
    from pymol import cmd
    import time
    
    def jump (*args, **kwargs):
        
        """
        Author Sean M. Law
        University of Michigan
        seanlaw_(at)_umich_dot_edu
        """
    
        keystroke='pgdn'
        
        for key in kwargs:
            if (key == "key"):
                keystroke=kwargs["key"]
                keystroke=check_key(keystroke)
                if (keystroke == None):
                    return
            else:
                continue
    
        cmd.set_key(keystroke, jumpit, args, kwargs)
    
        return
    
    cmd.extend("jump", jump)
    
    def jumpit (*args, **kwargs):
    
        x=1
        sleep=0.1
    
        for key in kwargs:
            if (key == "sleep"):
                sleep=kwargs["sleep"]
                sleep=float(sleep)
            elif (key == "x"):
                x=kwargs["x"]
                x=int(x)
            else:
                continue
    
        args=list(args)
        
        for i in range (x):
            for j in range (len(args)):
                cmd.frame(int(args[j]))
                cmd.refresh()
                time.sleep(sleep)
    
        return
    cmd.extend("jumpit", jumpit)
    

  


### See Also

[Check Key](/index.php/Check_Key "Check Key")

Retrieved from "[https://pymolwiki.org/index.php?title=Jump&oldid=11148](https://pymolwiki.org/index.php?title=Jump&oldid=11148)"


---

## Kabsch

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Note:** PyMOL has built-in commands to do RMSD fitting. This script is typically not needed. 

In particular, `optAlign (sele1), (sele2)` is identical to `fit (sele1), (sele2), matchmaker=-1`. 

See also: [fit](/index.php/Fit "Fit"), [align](/index.php/Align "Align"). 

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [GPL](http://opensource.org/licenses/GPL-2.0)  
  
## Contents

  * 1 Intro
  * 2 To use
    * 2.1 Caveats
  * 3 Notes
    * 3.1 Examples
  * 4 The Code
    * 4.1 The Old Code
  * 5 References



## Intro

The Kabsch algorithm uses linear and vector algebra to find the optimal rotation and translation of two sets of points in N-dimensional space as to minimize the RMSD between them. The following program is a Python implementation of the Kabsch algorithm. 

This program when called will align the two selections, optimally, convert the proteins in the selection to ribbons and change the color of the selections to show the matched alignments. 

**WHAT THIS DOESN'T DO** : This program does NOT provide a pairwise alignment of two structures from scratch. You have to tell it what the equivalent items are. See [Cealign](/index.php/Cealign "Cealign"). 

**NOTE:** This has **NOT** been tested on any other machine than mine, by me. It works on all PyMols 0.97 and newer (haven't tested earlier versions) and I use Python 2.3's Numeric Version 23.3. 

**NOTE:** I have added new Kabsch code. The new code uses SVD, and fixed an old bug. For ease of use, try the new code (requires numpy, though). 

## To use

  1. Save this script to Kabsch.py
  2. Open PyMol
  3. Load the alignment script: **run Kabsch.py** (The command **optAlign** is now defined in PyMol.)
  4. Load your proteins
  5. Align the proper segments (see below examples)



To align two equivalent sets of residues do: 
    
    
    optAlign SEL1 and n. CA and i. a-b, SEL2 and n. CA and i. c-d
    

where 

  * **SEL1** is the first protein
  * **a-b** is the range of residues to align in the first protein
  * **SEL2** is the second protein
  * **c-d** is the range of residues to align in the second protein



### Caveats

  * Ensure that you're equivalencing **N** atoms to **N** atoms (run [Count_Atoms](/index.php/Count_Atoms "Count Atoms") over your two selections to ensure they are the same length).
  * Sometimes PyMol doesn't seem to superimpose them right the first time. Hit the up-arrow and rerun the program if this happens. It always superimposes the correctly the second time. I think it has something to do with the orientation. I'll fix this when I find the error.
  * The RMSD is only between the equivalent atoms. Use PyMol's [Rms_Cur](/index.php/Rms_Cur "Rms Cur") if you want a full RMSD.
  * Make sure your atom selections are numbered correctly. Many times PDB files start residue numbers at something other than 0 or 1. To ensure you're aligning the right things, do 
        
        set seq_view,1
        

to turn on the sequence viewer and double check your residue numbers. If a protein has residue one numbered as something other than one, say 2064, simply run 
        
        alter (SEL), resi=str(int(resi)-2064)
        

and then 
        
        sort
        

where **SEL** is the name of the protein and 2064 is the offset to adjust by. Your protein will now work as needed. See [Alter](/index.php/Alter "Alter"). This capability is also provided in a script; See [zero_residues](/index.php/Zero_residues "Zero residues").  




## Notes

  1. Windows users are having problems running the script. Python tells them first off "TypeError: Can't convert rank-0 arrays to Python scalars." The fix to that breaks some code in Numeric -- which I don't maintain.
  

  2. However, to make this work, you can change the code in **Numeric.py** supplied with Pymol, located in the folder "<Pymol Home>\modules\Numeric\" (for example: "C:\Program Files\DeLano Scientific\PyMOL\modules\Numeric").  
  
Essentially, you need to search for the line: 
         
         if axis2 < 0: axis2 = axis1 + n  # (should be around line 250)
         

and replace it with: 
         
         if axis2 < 0: axis2 = axis2 + n
         




### Examples
    
    
    optAlign 1cll and n. CA and i. 4-20+30-60, 1ggz and n. CA and i. 4-20+30-60
    
    
    
    optAlign 1kao and n. CA and i. 20-50, 1ctq and n. CA and i. 20-50
    

  * [![1cll and 1ggz loaded](/images/b/be/OptAlign1.png)](/index.php/File:OptAlign1.png "1cll and 1ggz loaded")

1cll and 1ggz loaded 

  * [![1cll and 1ggz aligned to residues 5-50+55-80 shown in red](/images/0/0d/OptAlign2.png)](/index.php/File:OptAlign2.png "1cll and 1ggz aligned to residues 5-50+55-80 shown in red")

1cll and 1ggz aligned to residues 5-50+55-80 shown in red 




Kabsch can also align hetero-atoms: 
    
    
    load 1cll.pdb
    load 1ggz.pdb
    optAlign 1cll and e. CA, 1ggz and e. CA
    

The above aligns the 4 Calciums in each structure. 

## The Code
    
    
    #!python
     
    ##############################################################################
    #
    # @SUMMARY: -- QKabsch.py.  A python implementation of the optimal superposition
    #     of two sets of vectors as proposed by Kabsch 1976 & 1978.
    #
    # @AUTHOR: Jason Vertrees
    # @COPYRIGHT: Jason Vertrees (C), 2005-2007
    # @LICENSE: Released under GPL:
    # This program is free software; you can redistribute it and/or modify
    #    it under the terms of the GNU General Public License as published by
    #    the Free Software Foundation; either version 2 of the License, or
    #    (at your option) any later version.
    # This program is distributed in the hope that it will be useful, but WITHOUT
    # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    # FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License along with
    # this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
    # Street, Fifth Floor, Boston, MA 02110-1301, USA 
    #
    # DATE  : 2007-01-01
    # REV   : 2
    # REQUIREMENTS: numpy
    #
    #############################################################################
    from array import *
     
    # system stuff
    import os
    import copy
     
    # pretty printing
    import pprint
     
    # for importing as a plugin into PyMol
    from pymol import cmd
    from pymol import stored
    from pymol import selector
     
    # using numpy for linear algebra
    import numpy
     
    def optAlign( sel1, sel2 ):
    	"""
    	optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
    	Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
    	Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA
     
    	Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
    	PyMol based upon your selections.
     
    	By default, this program will optimally align the ALPHA CARBONS of the selections provided.
    	To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
     
    	@param sel1: First PyMol selection with N-atoms
    	@param sel2: Second PyMol selection with N-atoms
    	"""
    	cmd.reset()
     
    	# make the lists for holding coordinates
    	# partial lists
    	stored.sel1 = []
    	stored.sel2 = []
    	# full lists
    	stored.mol1 = []
    	stored.mol2 = []
     
    	# -- CUT HERE
    	sel1 += " and N. CA"
    	sel2 += " and N. CA"
    	# -- CUT HERE
     
    	# Get the selected coordinates.  We
    	# align these coords.
    	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
    	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
     
    	# get molecule name
    	mol1 = cmd.identify(sel1,1)[0][0]
    	mol2 = cmd.identify(sel2,1)[0][0]
     
    	# Get all molecule coords.  We do this because
    	# we have to rotate the whole molcule, not just
    	# the aligned selection
    	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
    	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
     
    	# check for consistency
    	assert len(stored.sel1) == len(stored.sel2)
    	L = len(stored.sel1)
    	assert L > 0
     
    	# must alway center the two proteins to avoid
    	# affine transformations.  Center the two proteins
    	# to their selections.
    	COM1 = numpy.sum(stored.sel1,axis=0) / float(L)
    	COM2 = numpy.sum(stored.sel2,axis=0) / float(L)
    	stored.sel1 -= COM1
    	stored.sel2 -= COM2
     
    	# Initial residual, see Kabsch.
    	E0 = numpy.sum( numpy.sum(stored.sel1 * stored.sel1,axis=0),axis=0) + numpy.sum( numpy.sum(stored.sel2 * stored.sel2,axis=0),axis=0)
     
    	#
    	# This beautiful step provides the answer.  V and Wt are the orthonormal
    	# bases that when multiplied by each other give us the rotation matrix, U.
    	# S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    	V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(stored.sel2), stored.sel1))
     
    	# we already have our solution, in the results from SVD.
    	# we just need to check for reflections and then produce
    	# the rotation.  V and Wt are orthonormal, so their det's
    	# are +/-1.
    	reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
     
    	if reflect == -1.0:
    		S[-1] = -S[-1]
    		V[:,-1] = -V[:,-1]
     
    	RMSD = E0 - (2.0 * sum(S))
    	RMSD = numpy.sqrt(abs(RMSD / L))
     
    	#U is simply V*Wt
    	U = numpy.dot(V, Wt)
     
    	# rotate and translate the molecule
    	stored.sel2 = numpy.dot((stored.mol2 - COM2), U)
    	stored.sel2 = stored.sel2.tolist()
    	# center the molecule
    	stored.sel1 = stored.mol1 - COM1
    	stored.sel1 = stored.sel1.tolist()
     
    	# let PyMol know about the changes to the coordinates
    	cmd.alter_state(1,mol1,"(x,y,z)=stored.sel1.pop(0)")
    	cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
     
    	print("RMSD=%f" % RMSD)
     
    	# make the alignment OBVIOUS
    	cmd.hide('everything')
    	cmd.show('ribbon', sel1 + ' or ' + sel2)
    	cmd.color('gray70', mol1 )
    	cmd.color('paleyellow', mol2 )
    	cmd.color('red', 'visible')
    	cmd.show('ribbon', 'not visible')
    	cmd.center('visible')
    	cmd.orient()
    	cmd.zoom('visible')
     
    cmd.extend("optAlign", optAlign)
    

### The Old Code
    
    
    #!python
    
    ##############################################################################
    #
    # @SUMMARY: -- Kabsch.py.  A python implementation of the optimal superposition
    #     of two sets of vectors as proposed by Kabsch 1976 & 1978.
    #
    # @AUTHOR: Jason Vertrees
    # @COPYRIGHT: Jason Vertrees (C), 2005-2007
    # @LICENSE: Released under GPL:
    # This program is free software; you can redistribute it and/or modify
    #    it under the terms of the GNU General Public License as published by
    #    the Free Software Foundation; either version 2 of the License, or
    #    (at your option) any later version.
    # This program is distributed in the hope that it will be useful, but WITHOUT
    # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    # FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License along with
    # this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
    # Street, Fifth Floor, Boston, MA 02110-1301, USA 
    #
    # DATE  : 2005-04-07
    # REV   : 2
    # NOTES: Updated RMSD, notes, cleaned up the code a little.
    #
    #############################################################################
    # math imports
    import math
    import Numeric
    import LinearAlgebra
    import Matrix
    
    from array import *
    
    # system stuff
    import os
    import copy
    
    # pretty printing
    import pprint
    
    # for importing as a plugin into PyMol
    #import tkSimpleDialog
    #import tkMessageBox
    from pymol import cmd
    from pymol import stored
    from pymol import selector
    
    class kabsch:
    	"""
    	Kabsch alignment of two set of vectors to produce and optimal alignment.
    
    	Steps
    	=====	
    		1.  Calculate the center of mass for each protein.  Then, move the protein
    		to its center of mass.  We choose as a convention, to use the origin as 
    		the center of mass of the first set of coordinates.  This will allow us to
    		return one translation vector, instead of two.
    		
    		Update: Since any rotation around a point that's not the origin, is in fact
    		an affine rotation.  So, to beat that, we translate both to the center.
    		
    		NAME: superpose(c1, c2)
    		POST: T is now defined.
    		
    		2.  Calculate the matrix, R, by (eq7, 1976).  r_{i,j} = sum_n w_n * y_{ni} * x_{nj},
    		where y_{ni} is the ith component of the vector y_n.
    		NAME: calcR
    		POST: R is now defined
    		
    		3.  Calculate RtR (R-transpose * R).
    		NAME: calcRtR
    		POST: RtR is now defined
    		
    		4.  Calculate the corresponding eigenpairs for RtR, and sort them accordingly:
    		m1 >= m2 >= m3; set v3 = v1 x v2 to ensure a RHS
    		NAME: calcEigenPairs
    		POST: The eigen pairs are calculated, sorted such that m1 >= m2 >= m3 and
    		v3 = v1 x v2.
    		
    		5.  Calculate R*v_k and normalize the first two vectors to obtain b_1, b_2 and set
    		b_3 = b_1 x b_2.  This also takes care of m1 > m2 = 0. 
    		NAME: calcBVectors
    		POST: b-Vectors are defined and returned
    		
    		6.  Calculate U=(u_{ij})=(sum n b_{ki} * a_{kj}) to obtain the best rotation.  Set
    		sigma_3 = -1 if b3(Ra3) < 0 else sigma_3 = +1.
    		NAME: calcU
    		POST: U is defined
    		
    		7.  Calculate the RMSD.  The residual error is then
    		The E = E0 - sqrt(m1) - sqrt(m2) - sigma_3(sqrt(m3)).
    		NAME: calcRMSD
    		POST: RMSD is computed.
    	
    	 @note: This should be a static method that takes three parameters
    		
    		1. The first protein's coordinates, f.  This program will 
    		accept coordinates in the form an array of 3D vectors/lists
    		
    		2. The second protein's coordinates, g.
    		
    		3. The array of integer pairs representing the pairs to align.
    		Coordinates should be formatted as as array of 2D vectors/lists.
    	
    
    	
    	"""
    	
    	def __init__(self):
    		"""
    		Constructor.  @see kabsch.align.
    		
    		"""
    
    		#
    		# Instance Variables:  All three of these will be updated
    		# every time the user calls ~.align.  Just to warn ya'.
    		#		
    		# U, the rotation matrix
    		self.U = []
    		# T, the translation vector
    		self.T = []
    		# R, the RMSD
    		self.R = -1.0
    
    		#self.menuBar.addmenuitem('Plugin', 'command', 'Kabsch Align', label = "Align Selections to Optial RMSD", command = lamda s=self: fetchPDBDialog(s))
    		
    				
    	def align(self, c1, c2, pairs):
    		"""
    		Finds the best alignment of c1 and c2's pairs resulting in
    		the smallest possible RMSD.
    
    				
    		@note:
    			- All weights in this first version are set to 1.  Kabsch allows,
    			differential weighting.  In the future, I may extend to this option,
    			and then pairs may become 3D (r1, r2, weight) or I may add another
    			parameter.
    		
    			- Helper functions will soon be provided such that the user may
    			just use this package alone to compute the rotation.
    		
    		@param c1: coordinats of the first vectors, as an array of 3D vectors.
    		@type  c1: Python list
    		   
    		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
    		@type  c2: Python list
    		   
    		@param pairs: the list of pairs as an array of 2D pairs.
    		@type  pairs: Python list of 2D lists.
    		
    		@return: U, the rotation matrix that gives the optimal rotation between the two proteins.
    			
    			T, the translation vector given to align the two objects centers of mass.
    			
    			R, the RMSD of the final rotation.
    		"""
    		
    		#
    		# First we move the center of mass of one protein, to the
    		# center of mass of the other.  This removes any translation
    		# between the two.
    		#
    		T1, T2, c1, c2 = self.superpose(c1, c2)
    		# Calculate the initial RMSD
    		E0 = self.calcE0(c1, c2)
    		# Calculate R via eq. 7.
    		R = self.calcR(c1, c2)
    		# Calculate R(transpose)*R
    		RtR = self.calcRtR(R)
    		# Determined the eigenpairs for the matrix RtR.
    		eValues, eVectors = self.calcEigenPairs(RtR)
    		# Determine the bVectors as required
    		bVectors = self.calcBVectors(R, eVectors)
    		# Calculate the roation matrix
    		U = self.calcU(eVectors, bVectors)
    		# Calculate the final RMSD using U.
    		RMSD = self.calcRMSD(E0, eValues, eVectors, bVectors, R, len(c1))
    		
    		return U, T1, T2, RMSD, c1, c2
    		
    		
    	def superpose(self, c1, c2 ):
    		"""
    		Calculate the center of mass for each protein.  Then, move the protein
    		to its center of mass.  We choose as a convention, to use the origin as 
    		the center of mass of the first set of coordinates.  This will allow us to
    		return one translation vector, instead of two.
    		(CORRECT)
    		
    		@precondition: c1 and c2 are well defined lists of N-dimensional points with length > 0.
    		@postcondition: T is now defined.
    		
    		@param c1: coordinats of the first vectors, as an array of 3D vectors.
    		@type  c1: Python list
    		   
    		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
    		@type  c2: Python list
    		
    		@return: T the translation vector.
    		
    			c2 one list of coordinates that of c2 translated to the COM of c1.
    		
    		"""
    		
    		# make sure we don't get bad data
    		if len(c1) != len(c2):
    			print "Two different length selections, with lengths, %d and %d." % (len(c1), len(c2))
    			print "This algorithm must be used with selections of the same length."	
    			print "In PyMol, type 'count_atoms sel1' where sel1 are your selections to find out their lengths."
    			print "Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122"
            		print "Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA"
    
    			
    		assert len(c1) == len(c2) != 0
    
    		L = len(c1)
    		
    		#
    		# Centers of Mass
    		#
    		c1COM = Numeric.zeros((3,1), Numeric.Float64)
    		c2COM = Numeric.zeros((3,1), Numeric.Float64)
    
    		# calculate the CsOM		
    		for i in range(L):
    			for j in range(3):
    				c1COM[j] += c1[i][j]
    				c2COM[j] += c2[i][j]
    		
    		T1 = - c1COM / L
    		T2 = - c2COM / L
    
    		# move everything back to the origin.
    		for i in range(L):
    			for j in range(3):
    				c1[i][j] += T1[j]
    				c2[i][j] += T2[j]
    				
    		return T1, T2, c1, c2
    
    					
    	def calcR( self, c1, c2 ):
    		"""
    		Calculate the matrix, R, by (eq7, 1976).  M{r_{i,j} = sum_n w_n * y_{ni} * x_{nj}},
    		where M{y_{ni}} is the ith component of the vector M{y_n}.
    		(CORRECT)
    		
    		@param c1: coordinates of the first vectors, as an array of 3D vectors.
    		@type  c1: Python list
    		   
    		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
    		@type  c2: Python list
    
    		@postcondition: R is now defined.
    		
    		@return: R
    		@rtype : 3x3 matrix
    				
    		"""
    		
    		# Create the 3x3 matrix
    		R = Numeric.zeros((3,3), Numeric.Float64)
    		L = len(c1)
    		
    		for k in range(L):
    			for i in range(3):
    				for j in range(3):
    					R[i][j] += c2[k][i] * c1[k][j]
    		
    		# return R the 3x3 PSD Matrix.
    		return R
    		
    	
    	def calcRtR( self, R ):
    		"""
    		Calculate RtR (R-transpose * R).
    		(CORRECT)
    		
    		@param R: Matrix
    		@type  R: 3x3 Matrix
    		
    		@precondition: R is a the well formed matrix as per Kabsch.
    		@postcondition: RtR is now defined
    		
    		@return: M{R^tR}
    		@rtype : 3x3 matrix
    		
    		"""
    		
    		RtR = Numeric.matrixmultiply(Numeric.transpose(R), R)
    		
    		return RtR
    	
    	
    	def calcEigenPairs( self, RtR ):
    		"""
    		Calculate the corresponding eigenpairs for RtR, and sort them accordingly:
    		M{m1 >= m2 >= m3}; set M{v3 = v1 x v2} to ensure a RHS
    		(CORRECT)
    		
    		@postcondition: The eigen pairs are calculated, sorted such that M{m1 >= m2 >= m3} and
    		M{v3 = v1 x v2}.
    		
    		@param RtR: 3x3 Matrix of M{R^t * R}.
    		@type  RtR: 3x3 Matrix
    		@return : Eigenpairs for the RtR matrix.
    		@rtype  : List of stuff
    		
    		"""
    		
    		eVal, eVec = LinearAlgebra.eigenvectors(RtR)
    
    		# This is cool.  We sort it using Numeric.sort(eVal)
    		# then we reverse it using nifty-crazy ass notation [::-1].
    		eVal2 = Numeric.sort(eVal)[::-1]
    		eVec2 = [[],[],[]] #Numeric.zeros((3,3), Numeric.Float64)
    				
    		# Map the vectors to their appropriate owners		
    		if eVal2[0] == eVal[0]:
    			eVec2[0] = eVec[0]
    			if eVal2[1] == eVal[1]:
    				eVec2[1] = eVec[1]
    				eVec2[2] = eVec[2]
    			else:
    				eVec2[1] = eVec[2]
    				eVec2[2] = eVec[1]
    		elif eVal2[0] == eVal[1]:
    			eVec2[0] = eVec[1]
    			if eVal2[1] == eVal[0]:
    				eVec2[1] = eVec[0]
    				eVec2[2] = eVec[2]
    			else:
    				eVec2[1] = eVec[2]
    				eVec2[2] = eVec[0]
    		elif eVal2[0] == eVal[2]:
    			eVec2[0] = eVec[2]
    			if eVal2[1] == eVal[1]:
    				eVec2[1] = eVec[1]
    				eVec2[2] = eVec[0]
    			else:
    				eVec2[1] = eVec[0]
    				eVec2[2] = eVec[1]
    
    		eVec2[2][0] = eVec2[0][1]*eVec2[1][2] - eVec2[0][2]*eVec2[1][1]
    		eVec2[2][1] = eVec2[0][2]*eVec2[1][0] - eVec2[0][0]*eVec2[1][2]
    		eVec2[2][2] = eVec2[0][0]*eVec2[1][1] - eVec2[0][1]*eVec2[1][0]
    		
    		return [eVal2, eVec2]
    		
    	
    	def calcBVectors( self, R, eVectors ):
    		"""
    		Calculate M{R*a_k} and normalize the first two vectors to obtain M{b_1, b_2} and set
    		M{b_3 = b_1 x b_2}.  This also takes care of {m2 > m3 = 0}. 
    		(CORRECT)
    		
    		@postcondition: b-Vectors are defined and returned
    		
    		@return: The three B-vectors
    		@rtype: List of 3D vectors (Python LOL).
    		"""
    		
    		bVectors = Numeric.zeros((3,3), Numeric.Float64)
    
    		for i in range(3):
    			bVectors[i] = Numeric.matrixmultiply(R, eVectors[i])
    
    		bVectors[0] = bVectors[0] / Numeric.sqrt(Numeric.add.reduce(bVectors[0]**2))
    		bVectors[1] = bVectors[1] / Numeric.sqrt(Numeric.add.reduce(bVectors[1]**2))
    		bVectors[2] = bVectors[2] / Numeric.sqrt(Numeric.add.reduce(bVectors[2]**2))
    		
    		bVectors[2][0] = bVectors[0][1]*bVectors[1][2] - bVectors[0][2]*bVectors[1][1]
    		bVectors[2][1] = bVectors[0][2]*bVectors[1][0] - bVectors[0][0]*bVectors[1][2]
    		bVectors[2][2] = bVectors[0][0]*bVectors[1][1] - bVectors[0][1]*bVectors[1][0]
    		
    		return bVectors
    		
    		
    		
    	def calcU( self, eVectors, bVectors ):
    		"""
    		Calculate M{U=(u_{ij})=(sum n b_{ki} * a_{kj})} to obtain the best rotation.  Set
    		M{sigma_3 = -1 if b3(Ra3) < 0 else sigma_3 = +1}.
    		(CORRECT)
    		
    		@postcondition: U is defined
    		
    		@param eVectors: Eigenvectors for the system.
    		@type  eVectors: Eigenvectors
    		
    		@param bVectors: BVectors as described by Kabsch.
    		@type  bVectors: BVectors
    		
    		@return: U the rotation matrix.
    		@rtype  :3x3 matrix.
    		
    		"""
    		
    		U = Numeric.zeros((3,3), Numeric.Float64)
    		
    		for k in range(3):
    			for i in range(3):
    				for j in range(3):
    					U[i][j] += Numeric.matrixmultiply(bVectors[k][i], eVectors[k][j])
    		
    		return U
    	
    		
    	def calcE0( self, c1, c2 ):
    		"""
    		Calculates the initial RMSD, which Kacbsch called E0.
    		(CORRECT)
    		
    		@param c1: coordinats of the first vectors, as an array of 3D vectors.
    		@type  c1: Python list
    		   
    		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
    		@type  c2: Python list
    		
    		@return: E0 the initial RMSD.
    		@rtype : float.
    				
    		"""
    		
    		E0 = 0.0
    		
    		L = len(c1)
    		for i in range( L ):
    			for j in range(3):
    				E0 += 0.5*( c1[i][j]*c1[i][j]+c2[i][j]*c2[i][j])
    		
    		return E0
    	
    	def calcRMSD( self, E0, eValues, eVectors, bVectors, R, N):
    		"""
    		Calculate the RMSD.  The residual error is then
    		The M{E = E0 - sqrt(m1) - sqrt(m2) - sigma_3(sqrt(m3))}.
    		
    		@param E0: Initial RMSD as calculated in L{calcE0}.
    		@type  E0: float.
    		
    		@param eVectors: Eigenvectors as calculated from L{calcEigenPairs}
    		@type eVectors: vectors, dammit!
    		
    		@param bVectors: B-vectors as calc. from L{calcBVectors}
    		@type  bVectors: More vectors.
    		
    		@param R: The matrix R, from L{calcR}.
    		@type  R: 3x3 matrix.		
    
    		@param N: Number of equivalenced points
    		@type  N: integer
    		
    		@postcondition: RMSD is computed.
    		@return: The RMSD.
    		
    		"""
    		sigma3 = 0
    		if Numeric.matrixmultiply(bVectors[2], Numeric.matrixmultiply( R, eVectors[2])) < 0:
    			sigma3 = -1
    		else:
    			sigma3 = 1
    
    		E = math.sqrt( 2*(E0 - math.sqrt(eValues[0]) - math.sqrt(eValues[1]) - sigma3*(math.sqrt(eValues[2]))) / N)
    		
    		return E
    		
    		
    	def calcSimpleRMSD( self, c1, c2 ):
    		"""
    		Calculates the usual concept of RMSD between two set of points.  The CalcRMSD above
    		sticks to Kabsch's alignment method protocol and calculates something much different.
    		@see kabsch.calcRMSD
    		
    		@param c1: List of points #1
    		@type  c1: LOL
    		
    		@param c2: List of points #2
    		@type  c2: LOL
    		
    		@return: RMSD between the two
    		
    		"""
    		
    		RMSD = 0.0
    		for i in range(len(c1)):
    			for j in range(3):
    				RMSD += (c2[i][j]-c1[i][j])**2
    				
    		RMSD = RMSD / len(c1)
    		RMSD = Numeric.sqrt(RMSD)
    		return RMSD
    		
    		
    	#####################################################################
    	#
    	# UTIL Functions
    	def rotatePoints(self, U, c2):
    		"""
    		Rotate all points in c2 based on the rotation matrix U.
    		
    		@param U: 3x3 Rotation matrix
    		@type  U: 3x3 matrix
    		
    		@param c2: List of points to rotate
    		@type  c2: List of 3D vectors
    		
    		@return: List of rotated points
    		
    		"""
    	
    		L = len(c2)
    
    		for n in range(L):
    			c2[n][0] = c2[n][0] * U[0][0] + c2[n][1] * U[1][0] + c2[n][2] * U[2][0]
    			c2[n][1] = c2[n][0] * U[0][1] + c2[n][1] * U[1][1] + c2[n][2] * U[2][1]
    			c2[n][2] = c2[n][0] * U[0][2] + c2[n][1] * U[1][2] + c2[n][2] * U[2][2]
    		
    		return c2
    	
    	def writeU( self, U, fileName ):
    		"""
    		Convenience function.  Writes U to disk.
    		
    		"""
    		
    		if len(fileName) == 0:
    			fileName = "./U"
    			
    		outFile = open( fileName, "wb")
    		for i in range(3):
    			for j in range(3):
    				outFile.write( str(U[i][j]).ljust(20) )
    			outFile.write("\n")
    		outFile.close()		
    
    	
    
    def optAlign( sel1, sel2 ):
    	"""
    	optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
    	Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
    	Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA
    	
    	Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
    	PyMol based upon your selections.
    
    	By default, this program will optimally align the ALPHA CARBONS of the selections provided.
    	To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
    	
    	@param sel1: First PyMol selection with N-atoms
    	@param sel2: Second PyMol selection with N-atoms
    	"""
    	cmd.reset()
    
    	# make the lists for holding coordinates
    	# partial lists
    	stored.sel1 = []
    	stored.sel2 = []
    	# full lists
    	stored.mol1 = []
    	stored.mol2 = []
    
    	# now put the coordinates into a list
    	# partials
    
    	# -- REMOVE ALPHA CARBONS
    	sel1 += " and N. CA"
    	sel2 += " and N. CA"
    	# -- REMOVE ALPHA CARBONS
    
    	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
    	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
    	# full molecule
    	mol1 = cmd.identify(sel1,1)[0][0]
    	mol2 = cmd.identify(sel2,1)[0][0]
    	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
    	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
    
    	K = kabsch()
    	U, T1, T2, RMSD, c1, c2 = K.align(stored.sel1, stored.sel2, [])
    
    	stored.mol2 = map(lambda v:[T2[0]+((v[0]*U[0][0])+(v[1]*U[1][0])+(v[2]*U[2][0])),T2[1]+((v[0]*U[0][1])+(v[1]*U[1][1])+(v[2]*U[2][1])),T2[2]+((v[0]*U[0][2])+(v[1]*U[1][2])+(v[2]*U[2][2]))],stored.mol2)
    	#stored.mol1 = map(lambda v:[ v[0]+T1[0], v[1]+T1[1], v[2]+T1[2] ], stored.mol1)
    	stored.mol1 = map(lambda v:[ v[0]+T1[0], v[1]+T1[1], v[2]+T1[2] ], stored.mol1)
    
    	cmd.alter_state(1,mol1,"(x,y,z)=stored.mol1.pop(0)")
    	cmd.alter_state(1,mol2,"(x,y,z)=stored.mol2.pop(0)")
    	cmd.alter( 'all',"segi=''")
    	cmd.alter('all', "chain=''")
    	print "RMSD=%f" % cmd.rms_cur(sel1, sel2)
    	print "MY RMSD=%f" % RMSD
    	cmd.hide('everything')
    	cmd.show('ribbon', sel1 + ' or ' + sel2)
    	cmd.color('gray70', mol1 )
    	cmd.color('paleyellow', mol2 )
    	cmd.color('red', 'visible')
    	cmd.show('ribbon', 'not visible')
    	cmd.center('visible')
    	cmd.orient()
    	cmd.zoom('visible')
    
    cmd.extend("optAlign", optAlign)
    

## References

**[Kabsch, 1976]** Kabsch, W. (1976). A solution for the best rotation to relate two sets of vectors. _Acta. Crystal_ , 32A:922-923. 

**[Kabsch, 1978]** Kabsch, W. (1978). A discussion of the solution for the best rotation to related two sets of vectors. _Acta. Crystal_ , 34A:827-828. 

Retrieved from "[https://pymolwiki.org/index.php?title=Kabsch&oldid=12910](https://pymolwiki.org/index.php?title=Kabsch&oldid=12910)"


---

## Key Wait

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #use "spawn spawn_demo.py, local" to invoke this
    #python script from within PyMOL
    
    wait=""
    i=0
    while i < 12 and wait!="x":
      cmd.turn("y", 30)
      print "Press enter key to continue or x + enter to terminate"
      wait=raw_input()
      i=i+1
    print "Done"
    

_Note: this approach only works when you have a console window open._

Retrieved from "[https://pymolwiki.org/index.php?title=Key_Wait&oldid=11390](https://pymolwiki.org/index.php?title=Key_Wait&oldid=11390)"


---

## LatticeGenerator

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is a simple script to generate a repeating atom lattice. It just writes atoms in PDB format, and PyMOL will automatically [connect adjacent atoms](/index.php/Connect_mode "Connect mode") when loading the PDB file. You may adjust the `c` and/or `x_shift`, `y_shift` arrays to obtain different geometries, and the number of iterations in the loops to adjust the lattice size. 

See also the [thread on pymol-users mailing list](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg11139.html). 
    
    
    from numpy import array
    
    # coordinates of the repeating unit (from cyclohexane)
    c = array([[ 0.064,   2.851,  -1.085 ],
               [ 0.260,   1.969,   0.159 ]])
    x_shift = array([ 1.67441517, -0.91605961,  1.66504574])
    y_shift = array([-0.69477826, -0.40578592,  2.40198410])
    
    # template string for PDB hetatom line
    s = 'HETATM %4d  C03 UNK     1    %8.3f%8.3f%8.3f  0.00  0.00           C  '
    
    out = open('lattice.pdb', 'w')
    
    i = 0
    for x in range(10):
        for y in range(10):
            for v in (c + (x-y//2) * x_shift + y * y_shift):
                i += 1
                out.write(s % (i, v[0], v[1], v[2]) + "\n")
    
    out.close()
    
    try:
        from pymol import cmd
        cmd.load('lattice.pdb')
    except ImportError:
        print('Please load lattice.pdb with PyMOL')
    

Retrieved from "[https://pymolwiki.org/index.php?title=LatticeGenerator&oldid=12720](https://pymolwiki.org/index.php?title=LatticeGenerator&oldid=12720)"


---

## Launching From a Script

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The recommended way to run PyMOL-Python scripts is by using PyMOL as the interpreter. This is supported by all versions of PyMOL, including the pre-compiled bundles provided by Schrödinger. 

Example from a shell: 
    
    
    shell> pymol -cq script.py
    

With arguments (_sys.argv_ becomes _["script.py", "foo", "bar"]_): 
    
    
    shell> pymol -cq script.py -- foo bar
    

Example from a running PyMOL instance: 
    
    
    PyMOL> run script.py
    

For advanced users, the following PyMOL versions also allow to run PyMOL from an existing Python process: 

  * [PyMOL 2.0](https://pymol.org/2/#download) based on Anaconda (using Anaconda's python, which is included in bundles provided by Schrödinger)
  * Open-Source PyMOL
  * Schrödinger-provided "[Mac alternative X11-only build](http://pymol.org/download)" of the 1.8.x series



After importing the **pymol** module, PyMOL's event loop has to be started with a call to **pymol.finish_launching()** _(not supported on macOS)_. 

## Contents

  * 1 Usage
  * 2 Example 1
  * 3 Example 2
  * 4 STDOUT
  * 5 Independent PyMOL Instances
  * 6 Independent PyMOL Instances (Context Manager)
  * 7 See Also



## Usage

With PyMOL 2.1, calling any pymol.cmd function will automatically start a backend process without the GUI in the main thread. "finish_launching" should not be necessary, and will launch PyMOL in a new thread with an event loop, which will cause 100% CPU usage (at least with "-c"). 
    
    
    from pymol import cmd
    cmd.fragment('ala')
    cmd.zoom()
    cmd.png('/tmp/test.png', 300, 200)
    

Since PyMOL 1.7.4, the following form is sufficient: 

_This is not supported on macOS (see[bug report](https://github.com/schrodinger/pymol-open-source/issues/28))_
    
    
    import pymol
    pymol.finish_launching(['pymol', '-q'])
    

Before 1.7.4, either "-K" was needed as an additional argument, or arguments had to be assigned to __main__.pymol_argv or pymol.pymol_argv. 

## Example 1

Here is an example script that launches PyMol for stereo viewing on a [VisBox](http://www.visbox.com/boxMain.html). It runs PyMol fullscreen stereo, and disables the internal gui. The environment (PYTHON_PATH and PYMOL_PATH) must already be set up for this example to work (see Example 2 below for how to setup within the script). 
    
    
    #!/usr/bin/env python
     
    # Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
    import __main__
    __main__.pymol_argv = [ 'pymol', '-qei' ]
     
    import pymol
     
    # Call the function below before using any PyMOL modules.
    pymol.finish_launching()  # not supported on macOS
     
    from pymol import cmd
    cmd.stereo('walleye')
    cmd.set('stereo_shift', 0.23)
    cmd.set('stereo_angle', 1.0)
    

## Example 2

This script launches PyMOL without any GUI for scripting only. It enables tab-completion on the python command line and does the PyMOL environment setup (you need to adjust the **moddir** variable!). _Hint: You may save this as "pymol-cli" executable._
    
    
    #!/usr/bin/python2.6 -i
    
    import sys, os
    
    # autocompletion
    import readline
    import rlcompleter
    readline.parse_and_bind('tab: complete')
    
    # pymol environment
    moddir='/opt/pymol-svn/modules'
    sys.path.insert(0, moddir)
    os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')
    
    # pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
    import pymol
    pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
    pymol.finish_launching()
    cmd = pymol.cmd
    

## STDOUT

PyMOL captures **sys.stdout** and **sys.stderr** , to control it with it's own [feedback](/index.php/Feedback "Feedback") mechanism. To prevent that, save and restore both streams, e.g.: 
    
    
    import sys
    stdout = sys.stdout
    stderr = sys.stderr
    pymol.finish_launching(['pymol', '-xiq'])  # not supported on macOS
    sys.stdout = stdout
    sys.stderr = stderr
    

## Independent PyMOL Instances

It's possible to have multiple independent instances. 
    
    
    import pymol2
    p1 = pymol2.PyMOL()
    p1.start()
    
    p2 = pymol2.PyMOL()
    p2.start()
    
    p1.cmd.fragment('ala')
    p1.cmd.zoom()
    p1.cmd.png('/tmp/ala.png', 1000, 800, dpi=150, ray=1)
    
    p2.cmd.fragment('ser')
    p2.cmd.zoom()
    p2.cmd.png('/tmp/ser.png', 1000, 800, dpi=150, ray=1)
    
    p1.stop()
    p2.stop()
    

## Independent PyMOL Instances (Context Manager)

PyMOL 2.2 adds context manager support. 
    
    
    import pymol2
    with pymol2.PyMOL() as p1:
        p1.cmd.fragment('ala')
        p1.cmd.zoom()
        p1.cmd.png('/tmp/ala.png', 1000, 800, dpi=150, ray=1)
    

## See Also

  * [Command Line Options](/index.php/Command_Line_Options "Command Line Options")
  * [Jupyter](/index.php/Jupyter "Jupyter")



Retrieved from "[https://pymolwiki.org/index.php?title=Launching_From_a_Script&oldid=13294](https://pymolwiki.org/index.php?title=Launching_From_a_Script&oldid=13294)"


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

## List Colors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    # This is how to do it from the PyMOL command line or .pml script:
    #
    iterate all, print color
    
    #! /usr/bin/python
    #
    # and this in a Python script
    #
    import pymol
    pymol.color_list = []
    cmd.iterate('all', 'pymol.color_list.append(color)')
    print pymol.color_list
    

Retrieved from "[https://pymolwiki.org/index.php?title=List_Colors&oldid=6345](https://pymolwiki.org/index.php?title=List_Colors&oldid=6345)"


---

## List Secondary Structures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    import pymol
    pymol.stored_ss = []
    cmd.iterate('all', 'pymol.stored_ss.append(string.ljust(ss,1))')
    print string.join(pymol.stored_ss)
    

Retrieved from "[https://pymolwiki.org/index.php?title=List_Secondary_Structures&oldid=6346](https://pymolwiki.org/index.php?title=List_Secondary_Structures&oldid=6346)"


---

## List Selection

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    # Using PyMOL commands:
    #
    list=[]
    iterate (name ca),list.append((resi,resn))
    print list
    
    [('ASP', '1'), ('CYS', '2'), ('ALA', '3'), ('TRP', '4'), ('HIS', '5'), ('LEU',
     '6'), ('GLY', '7'), ('GLU', '8'), ('LEU', '9'), ('VAL', '10'), ('TRP', '11'), 
    ('CYS', '12'), ('THR', '13')]
    
    #!/usr/bin/python
    #
    # as Python script:
    #
    from pymol import cmd,stored
    stored.list=[]
    cmd.iterate("(name ca)","stored.list.append((resi,resn))")
    print stored.list
    
    [('1', 'ASP'), ('2', 'CYS'), ('3', 'ALA'), ('4', 'TRP'), ('5', 'HIS'), ('6', '
    LEU'), ('7', 'GLY'), ('8', 'GLU'), ('9', 'LEU'), ('10', 'VAL'), ('11', 'TRP'), 
    ('12', 'CYS'), ('13', 'THR')]
    
    # Note from Warren: 
    #
    # The above approach uses a the global pymol variable "stored"
    # In recent versions, "cmd.iterate" has been extended to take a dictionary
    # as an argument so that you no longer have to use a global variable.
    # Avoiding globals helps prevent conflicts between scripts.
    #
    from pymol import cmd
    my_dict = { 'my_list' : [] }
    cmd.iterate("(name ca)","my_list.append((resi,resn))",space=my_dict)
    print my_dict['my_list']
    
    [('1', 'ASP'), ('2', 'CYS'), ('3', 'ALA'), ('4', 'TRP'), ('5', 'HIS'), ('6', '
    LEU'), ('7', 'GLY'), ('8', 'GLU'), ('9', 'LEU'), ('10', 'VAL'), ('11', 'TRP'), 
    ('12', 'CYS'), ('13', 'THR')]
    

Retrieved from "[https://pymolwiki.org/index.php?title=List_Selection&oldid=6344](https://pymolwiki.org/index.php?title=List_Selection&oldid=6344)"


---

## ListSelection2

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [figshare](http://dx.doi.org/10.6084/m9.figshare.1031599)  
Author(s)  | [Pietro Gatti-Lafranconi](/index.php/User:PietroGattiLafranconi "User:PietroGattiLafranconi")  
License  | [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/)  
  
## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 Code



# Overview

Alternative script to [List_Selection](/index.php/List_Selection "List Selection") that will: 

  * list multiple residues once
  * include water and hetatoms
  * account for different chains/objects
  * produce a screen/printed output



  


# Usage
    
    
    listselection selection, [output=N/S/P, [HOH=Y/N ]]
    

where: 

  * **selection** can be any existing or newly-defined selection
  * **output** controls if the list is hidden (default), printed on screen (S) or saved in a file (P)
  * **HOH** (default Y) allows to exclude (N) water residues from the list



  


# Examples
    
    
    PyMOL>listselection 1efa, output=P
    Residues in '1efa': 1043
    Results saved in listselection_1efa.txt
    

  

    
    
    PyMOL>listselection sele, output=S, HOH=N
    Residues in 'sele, without HOH': 7
    1LBH/A/ARG/86
    1LBH/A/ALA/87
    1LBH/A/ASP/88
    1LBH/A/GLN/89
    1LBH/A/LEU/90
    1LBH/A/GLY/91
    1LBH/A/ALA/92
    

  

    
    
    PyMOL>listselection 1efa and resn LEU, output=S
    Residues in '1efa and resn LEU': 108
    1efa/A/LEU/6
    1efa/A/LEU/45
    1efa/A/LEU/56
    [...]
    1efa/C/LEU/323
    1efa/C/LEU/330
    

  


# Code

Copy the following text and save it as listselection.py 
    
    
    from pymol import cmd, stored
     
    def listselection (selection, output="N", HOH="Y"):
    	"""	
    	usage: listselection selection, [output=N/S/P, [HOH=Y/N ]]
    	
    	More information at: PymolWiki: http://http://pymolwiki.org/index.php/ListSelection2
    	AUTHOR: Pietro Gatti-Lafranconi, 2013
    	Please inform me if you use/improve/like/dislike/publish with this script.
    	CC BY-NC-SA
    	"""
    	printedselection=""
    	extra=""
    	counter=0
    	sel=selection
    	objs=cmd.get_object_list(sel)
    
    	if HOH=="N":
    		sel=selection+" and not resn HOH"
    		extra=", without HOH"
    	
    	for a in range(len(objs)):
    		m1=cmd.get_model(sel+" and "+objs[a])
    		for x in range(len(m1.atom)):
    			if m1.atom[x-1].resi!=m1.atom[x].resi:
    				printedselection+="%s/%s/%s/%s\n" % (objs[a], m1.atom[x].chain, m1.atom[x].resn, m1.atom[x].resi)
    				counter+=1
    				
    	print "Residues in '%s%s': %s" % (selection, extra, counter)
    	if output=="S": print printedselection
    	if output=="P":
    		f=open('listselection_'+selection+'.txt','w')
    		f.write("Residues in '%s%s': %s\n" % (selection, extra, counter))
    		f.write(printedselection)
    		f.close()
    		print "Results saved in listselection_%s.txt" % selection
    		
    cmd.extend('listselection',listselection)
    

Retrieved from "[https://pymolwiki.org/index.php?title=ListSelection2&oldid=11569](https://pymolwiki.org/index.php?title=ListSelection2&oldid=11569)"


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

## Load img stack

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/load_img_stack.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/load_img_stack.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the load_img_stack command, which loads a set of images as a map object. The images can either be individual files, or contained in a multi-page TIFF file. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    load_img_stack pattern [, name [, grid [, channel [, normalize [, extent ]]]]]
    

## Arguments

  * **pattern** = str: image filename or pattern
  * **name** = str: map object name to create
  * **grid** = float: grid spacing in Angstrom {default: 1.0}
  * **channel** = int: color channel for RGB images {default: 0}
  * **normalize** = 0 or 1: normalize data {default: 1}
  * **extent** = 3-float: (a,b,c) edge lengths in Angstrom, overwrites "grid" arguments if given {default: }



## Example
    
    
    load_img_stack img*.png, name=map, extent=(10.0, 10.0, 5.0)
    volume vol, map, rainbow
    

## See Also

  * [tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")
  * [volume](/index.php/Volume "Volume")



Retrieved from "[https://pymolwiki.org/index.php?title=Load_img_stack&oldid=13944](https://pymolwiki.org/index.php?title=Load_img_stack&oldid=13944)"


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

## LoadDir

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Load all files of the suffix **suff** from the directory **dirName** , where **suff** and **dirName** are function parameters. 

_This script is superseded by the[loadall](/index.php/Loadall "Loadall") command, which was added in PyMOL 1.7.2._

## Install

  1. copy the source below to the a file called "loadDir.pml" somewhere on your computer
  2. load the file with "run /your/path/toLoadDir/loadDir.pml"
  3. run loadDir. See examples below.



## Examples
    
    
    # load the script
    run ~/loadDir.pml
    
    # load all SD files from /tmp
    loadDir /tmp, sdf
    loadDir /tmp, .sdf
    loadDir /tmp, *.sdf
    # even stupid stuff works; hopefully as one would want.
    # load all PDBs from /tmp
    loadDir /tmp, foo.woo.pdb
    
    # load all the PDBs in all the directories under ./binders/ERK
    loadDir ./binders/ERK/*, .pdb
    
    # load the PDBs into groups: now we can load all the files in the tree under
    # ./ERK into the group "ERK" and the files from ./SYK into the group "SYK"
    loadDir ./binders/ERK/*, .pdb, group=ERKb
    loadDir ./binders/SYK/*, .pdb, group=SYKb
    

## The Code
    
    
    from glob import glob
    from os.path import sep, basename
    
    def loadDir(dirName=".", suff="pdb", group=None):
            """
            Loads all files with the suffix suff (the input parameter) from the directory dirName).
    
            dirName:        directory path
            suff:           file suffix.  Should be simply "pdb" or "sdf" or similar.  Will accept the
                            wildcard and dot in case the user doesn't read this.  So, "*.pdb", ".pdb",
                            and "pdb" should work.  The suffix can be anything valid that PyMOL knows
                            how to natively load.
            group:          groupName to add the files to.
    
            example:
                    # load all the PDBs in the current directory
                    loadDir
    
                    # load all SD files from /tmp
                    loadDir /tmp, "sdf"
    
            notes:
                    make sure you call this script w/o quotes around your parameters:
                            loadDir ., .pdb
                    as opposed to
                            loadDir ".", "*.pdb"
                    Use the former.
            """
    
            g = dirName + sep + "*." + suff.split(".")[-1]
    
            for c in glob( g ):
                    cmd.load(c)
    
                    if ( group != None ):
                            cmd.group( group, basename(c).split(".")[0], "add" )
    
    cmd.extend("loadDir", loadDir)
    

Retrieved from "[https://pymolwiki.org/index.php?title=LoadDir&oldid=12540](https://pymolwiki.org/index.php?title=LoadDir&oldid=12540)"


---

## Make Figures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 Examples
  * 4 The Code



# Overview

This script will aid you in making publication quality figures for the currently displayed scene. It understands a variety of preset "modes" and sizes. 

  


# Usage
    
    
    make_figure filename [, mode] [, size (default=900 pixels)] [,opaque]
    

The parameters are: 

**filename** : The name of the resulting image file. 

    

    _The extension**.png** is added automatically._

  


_NOTE_ : if no further arguments are given, the script generates quick figures for general use: 

    

    

  * figures are not ray-traced
  * TWO figures at 300 x 300 px, 72 dpi
  * views are rotated 180° about y (front and back view)
  * _front_quick and _back_quick are appended to the filename



  


**mode** : Type of figure desired. Possible values are as follows: 

  


    

    _**single**_ = "single view" 

  * one ray-traced 300 dpi figure of the current view



  


    

    _**fb**_ = "front and back view" 

  * TWO ray-traced 300 dpi figures
  * views are rotated 180° about y
  * _front and _back are appended to the filename



  


    

    _**sides**_ = "four side views" 

  * FOUR ray-traced 300 dpi figures
  * view is incrementally rotated by 90° about y
  * _1, _2, _3, and _4 are appended to the filename



  


    

    _**stereo**_ = "stereoview" 

  * TWO ray-traced 300 dpi figures
  * views are shifted by +/- 3 degrees
  * image dimensions are fixed at 750 x 750 pixels (size arguments are ignored)
  * _L and _R are appended to the filename
  * the output files are meant to be combined side by side to generate a stereo image



  


**size** : Size of the figure(s) in pixels or in "panels" 

    

    

  * DEFAULT = 900 pixels if a mode is specified
  * if size is 12 or smaller, the script interprets it as a number of "panels" to make.
  * panel sizes in pixels are hardcoded in the source but can easily be modified.



  


**opaque** : Create an opaque background. 

    

    

  * By default the figure background is 100% transparent.



  


# Examples
    
    
    #quick 'n' dirty front and back views of the scene
    #300x300 pixels at 72 dpi, transparent background
    # filenames will be output_front_quick.png and output_back_quick.png
    make_figure output
    
    # one ray-traced PNG file 975x975 pixels at 300dpi, with opaque background
    # filename will be output.png
    make_figure output, single, 975, opaque
    
    # two panels (1350x1350 px each at 300dpi) of the "front" and "back" view on transparent background
    # filenames will be output_front.png and output_back.png
    make_figure output, fb, 2
    
    # four panels (900x900 px each) where the view is incrementally rotated by 90° about y, transparent background
    # filenames will be output_1.png, output_2.png, output_3.png and output_4.png
    make_figure output, sides,4
    
    #stereoview of the current scene with opaque background
    #size is fixed to 2x 750x750px
    # filenames will be output_L.png and output_R.png
    make_figure output, stereo, opaque
    

  


# The Code
    
    
    #make_figure v.3.0
    #Copyleft Martin Christen, 2010
    
    from pymol import cmd
    def make_figure(output='', mode='',size=900,opaque='transparent'):
    
    	"""
    
    AUTHOR
    
    	Martin Christen
    
    
    DESCRIPTION
    
    	"make_figure" creates publication-quality figures of the current scene.
    	It understands several predefined "modes" and sizes.
    
    
    USAGE
    
    	make_figure filename [, mode] [, size (default=900 pixels)] [,opaque]
    
    
    ARGUMENTS
    
    	mode = string: type of desired figure (single, fb, sides, stereo or -nothing-)
    	size = integer: size of the figure in pixels OR # panels (if <= 12)
    	opaque = specify an opaque background.
    	         By default, the script makes the background transparent.
    EXAMPLES
    
    	make_figure output
    	make_figure output, single, 975, opaque
    	make_figure output, fb, 2
    	make_figure output, sides,4
    	make_figure output, stereo
    
    
    NOTES
    
    	"single" mode makes a single 300 dpi figure
    	
    	"fb" mode makes TWO 300 dpi figure
    	("front" and "back", rotating by 180 degrees about y)
    	
    	"sides" mode makes FOUR 300 dpi figures
    	("front" "left" "right" and back, rotating by 90 degrees clockwise about y)
    	
    	"stereo" generates two 300 dpi, 750 px figures
    	("L" and "R", to be combined as a stereo image)
    	If you specify the stereo mode, the size argument is IGNORED.
    		
    	If no mode argument is given, the script generates quick figures
    	for general	use: TWO figures (front and back) at 300 x 300 px, 72 dpi.
    	
    	Size is interpreted as pixels, except if the number is ridiculously small
    	(<=12),	in which case the script as "number of panels" to make.
    
    	Edit the script manually to define corresponding values.
    	
    	"""
    
    	#define sizes here (in pixels)
    	panel1 = 1800
    	panel2 = 1350
    	panel3 = 900
    	panel4 = 900
    	panel5 = 750
    	panel6 = 750
    	panel7 = 675
    	panel8 = 675
    	panel9 = 600
    	panel10 = 585
    	panel11 = 585
    	panel12 = 585
    	
    	#verify size is an integer number and convert to pixels
    	size = int(size)
    	if size > 12:
    		pixels = size
    	
    	elif size == 1:
    		pixels = panel1
    	
    	elif size == 2:
    		pixels = panel2
    	
    	elif size == 3:
    		pixels = panel3
    	
    	elif size == 4:
    		pixels = panel4
    	
    	elif size == 5:
    		pixels = panel5
    	
    	elif size == 6:
    		pixels = panel6
    	
    	elif size == 7:
    		pixels = panel7
    	
    	elif size == 8:
    		pixels = panel8
    	
    	elif size == 9:
    		pixels = panel9
    	
    	elif size == 10:
    		pixels = panel10
    	
    	elif size == 11:
    		pixels = panel11
    	
    	elif size == 3:
    		pixels = panel12
    
    	#change background
    	cmd.unset('opaque_background')
    	if opaque == 'opaque':
    		cmd.set('opaque_background')
    	
    	#apply mode
    	if output == '':
    		print 'no output filename defined\n'
    		print 'try: \'make_figure filename\''
    		return -1
    		# abort if no output file name given
    
    	if mode =='':
    		cmd.set('surface_quality',1)
    		cmd.set('opaque_background')
    		cmd.png(output+"_back_quick",300,300,dpi=72)
    		cmd.turn('y',180)
    		cmd.png(output+"_front_quick",300,300,dpi=72)
    		cmd.turn('y',180)
    		cmd.set('surface_quality',0)
    		# make front and back figures for quick mode
    
    	elif mode == 'single':
    		cmd.set('surface_quality',1)
    		cmd.set('ray_shadow',0)
    		cmd.ray(pixels, pixels)
    		cmd.png(output, dpi=300)
    		cmd.set('surface_quality',0)
    		# make a figure for single mode
    		
    	elif mode == 'fb':
    		cmd.set('surface_quality',1)
    		cmd.set('ray_shadow',0)
    		cmd.ray(pixels, pixels)
    		cmd.png(output+"_front", dpi=300)
    		cmd.turn('y',180)
    		cmd.ray(pixels, pixels)
    		cmd.png(output+"_back", dpi=300)
    		cmd.turn('y',180)
    		cmd.set('surface_quality',0)
    		# make front and back figures for single mode
    
    	elif mode == 'sides':
    		cmd.set('surface_quality',1)
    		cmd.set('ray_shadow',0)
    		cmd.ray(pixels, pixels)
    		cmd.png(output+"_1", dpi=300)
    		cmd.turn('y',90)
    		cmd.ray(pixels, pixels)
    		cmd.png(output+"_2", dpi=300)
    		cmd.turn('y',90)
    		cmd.ray(pixels, pixels)
    		cmd.png(output+"_3", dpi=300)
    		cmd.turn('y',90)
    		cmd.ray(pixels, pixels)
    		cmd.png(output+"_4", dpi=300)
    		cmd.turn('y',90)
    		cmd.set('surface_quality',0)
    		# make front and back figures for single mode
    		
    	elif mode == 'stereo':
    		cmd.set('surface_quality',1)
    		cmd.set('ray_shadow',0)
    		cmd.ray(750, 750, angle=-3)
    		cmd.png(output+"_R", dpi=300)
    		cmd.ray(750, 750, angle=3)
    		cmd.png(output+"_L", dpi=300)
    		# make stereo figure (for more control use stereo_ray)
    		
    cmd.extend('make_figure',make_figure)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Make_Figures&oldid=10985](https://pymolwiki.org/index.php?title=Make_Figures&oldid=10985)"


---

## MakeVinaCommand

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Usage
  * 3 The Code
  * 4 Example
  * 5 See Also



# Overview

[Vina](http://vina.scripps.edu/) is a new, very fast, molecular docking program, written by Oleg Trott. In order to run Vina, you should have [MGLTools1.5.4](http://mgltools.scripps.edu/downloads) installed (you need prepare_receptor4.py and prepare_ligand4.py). But, you don't need MGLTools for this script to run. 

To run Vina it needs to know the center of the protein and also its extents (size in XYZ directions). We can easily do this with PyMOL. 

Requirements: 

  * [COM](/index.php/COM "COM") \-- simple script to get the center of mass of a selection



Todo: 

  * Commandline options
  * Robustness to vagaries of input
  * Usage & Help



# Usage
    
    
    pymol -cq ./makeVinaCommandFromProtein.py proteinFile ligandFile
    

For high-throughput screening, where we compare each ligand to each protein, I typically wrap the script in a for-loop like: 
    
    
    # foreach protein in my proteins directory
    for protein in proteinsDir/*; do
    
      # foreach ligand in my ligands directory
      for ligand in ligandsDir/*; do;
    
        # make a Vina command to run the protein vs the ligand.
        # note the backticks to run the output
        `pymol -cq ./makeVinaCommandFromProtein.py $protein $ligand | grep vina`;
    
      done;
    done;
    

# The Code
    
    
    # -*- coding: utf-8 -*-
    #
    # makeVinaCommandFromProtein.py -- automatically make a valid Vina command from a PDB file and a ligand (name)
    #
    # Author: Jason Vertrees
    # Date  : 2/2009
    #
    from pymol import cmd
    from sys import argv
    from os import path
    
    # try to keep PyMOL quiet
    cmd.feedback("disable","all","actions")
    cmd.feedback("disable","all","results")
    
    # prepare some output names
    protName= path.basename(argv[-2])
    ligName = path.basename(argv[-1])
    outName = protName.split(".")[0] + "." + ligName.split(".")[0] + ".docked.pdbqt"
    logName = protName.split(".")[0] + "." + ligName.split(".")[0] + ".log"
    
    # very unsafe commandline checking; needs updating
    cmd.load(argv[-2], protName)
    cmd.load(argv[-1], ligName)
    
    # remove the ligand before calculating the center of mass
    cmd.delete(ligName)
    
    # load center of mass script
    cmd.do("run /home/path/to/COM.py")
    
    # calculate the center of mass and extents
    (comX, comY, comZ) = COM(protName)
    ((maxX, maxY, maxZ), (minX, minY, minZ)) = cmd.get_extent(protName)
    
    # print the command line
    print "vina --receptor "+protName"qt --ligand "+ligName+"qt --center_x ", str(comX), " --center_y ", str(comY),\
    " --center_z ", str(comZ), " --size_x ", str(abs(maxX-minX)), " --size_y ", str(abs(maxY-minY)), " --size_z ", \
    str(abs(maxZ-minZ)), " --all ", outName , " --exhaustiveness 200 --log ", logName, " \n"
    

# Example
    
    
    # execute the script
    > pymol -cq ./makeVinaCommandFromProtein.py 1ia6_gh09.pdb glucose.pdb | grep vina
    
    # the output
    vina --receptor  1ia6_gh09.pdbqt  --ligand  glucose_dimer.pdbqt --center_x  1.86797851457  --center_y  17.7951449088  --center_z  65.2250072289  --size_x  55.9499988556  --size_y  49.7459993362  --size_z  58.1769981384  --all  1ia6_gh09.glucose_dimer.docked.pdbqt  --exhaustiveness 100 --log  1ia6_gh09.glucose_dimer.log
    

# See Also

  * [COM](/index.php/COM "COM") \-- you need this script
  * [Vina](http://vina.scripps.edu/) \-- the docking software



Retrieved from "[https://pymolwiki.org/index.php?title=MakeVinaCommand&oldid=8421](https://pymolwiki.org/index.php?title=MakeVinaCommand&oldid=8421)"


---

## Mark center

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [mark_center.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/mark_center.py)  
Author(s)  | [User:Inchoate](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Mark_center is a small callback that puts a cross at the center of space, where the camera is pointing. 

You could also use this to setup arbitrary CGOs (points, lines, etc) in 3D space. 
    
    
    from pymol import cmd
    
    def crosshair_put_center():
       t = cmd.get_position()
       m = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0] + t + [1]
       cmd.set_object_ttt('crosshair', m)
    
    cmd.load_callback(crosshair_put_center, '_crosshair_cb')
    
    cmd.pseudoatom('crosshair', pos=(0,0,0))
    cmd.show_as('nonbonded', 'crosshair')
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mark_center&oldid=12541](https://pymolwiki.org/index.php?title=Mark_center&oldid=12541)"


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

## Modevectors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/modevectors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/modevectors.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
### DESCRIPTION

[![](/images/3/34/Modevectors.png)](/index.php/File:Modevectors.png)

[](/index.php/File:Modevectors.png "Enlarge")

Example showing modevectors in action. (See Examples below).

**modevectors.py** is a PyMol script that was originally written to visualize results obtained from Normal Mode Analysis (NMA) by drawing arrows or vectors from a starting structure to a final structure. However, this script can be used to visualize the direction of motion between two specified states (e.g. targeted MD, comparing open and closed structures, etc). The strength of this script is that the arrows are highly customizable with numerous options to choose from (see script below). It is important to note that the two states MUST be identical except for the X, Y, and Z coordinates. That is, you cannot draw vectors between homologous structures. The default settings sets the background to white and displays the arrows along with the first object frame (in cartoon representation). 

  
Note: Ray tracing these CGO arrows appears to be a memory intensive process so please save a session before ray tracing and use "pymol -qc" where necessary. If anybody can come up with a method to improve this please e-mail me (see address in script). 

Update: A new way of drawing the arrows has been implemented for PyMOL v.1.1 and has been updated in this code (which automatically detects the version that you are using). The new method uses the built in cone primitive rather than concentric circles decreasing in size. Thus, the new method is much less memory intensive and is a significant improvement over the old method. Additionally, you may want to turn off [ray shadow](/index.php/Ray_shadow "Ray shadow") depending on the [scene](/index.php/Scene "Scene") that you are viewing as shadows can be confusing to others when seen in publications. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    modevectors first_obj_frame, last_obj_frame [,first_state=1 [,last_state=1 [,outname=modevectors \
      [,head=1.0 [,tail=0.3 \[,head_length=1.5 [,headrgb=(1.0,1.0,1.0) [,tailrgb=(1.0,1.0,1.0) [,cutoff=4.0 
      [,skip=0 [,cut=0.5 [,atom=CA [,stat=show [,factor=1.0 [,notail=0]]]]]]]]]]]]]]
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "modevectors" (which can be changed through the options). 

  


### EXAMPLES
    
    
    modevectors 1E3M, 1W7A
    modevectors 1E3M, 1W7A, outname="arrows"
    modevectors 1E3M, 1W7A, atom="P"
    

Copy/paste the following code to see an example of modevectors. This uses a multistate protein and the arrows are connected between the first and last states. 
    
    
    import modevectors
    # fetch the PDBs from pdb.org
    fetch 1c3y, finish=1, multiplex=0, async=0
    # separate the first and last states of the NMR ensemble to individual objects
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
    # run the modevectors code
    modevectors 1c3y_0001, 1c3y_0023
    # just setup a nice representation
    as cartoon, 1c3y_0001 or 1c3y_0023
    show cgo, modevectors
    color marine, 1c3y_0001
    color purpleblue, 1c3y_0023
    

The following set of examples will illustrate the power of each optional argument. Each example should be compared to the default figure in the table below. 

[![](/images/3/35/Mv_default.png)](/index.php/File:Mv_default.png) [](/index.php/File:Mv_default.png "Enlarge")Fig.1 - Default Settings | [![](/images/b/b0/Mv_fig2.png)](/index.php/File:Mv_fig2.png) [](/index.php/File:Mv_fig2.png "Enlarge")Fig.2 - Arrow direction drawn from CA to CA | [![](/images/4/4f/Mv_fig3.png)](/index.php/File:Mv_fig3.png) [](/index.php/File:Mv_fig3.png "Enlarge")Fig.3 - Arrow head radius  
---|---|---  
[![](/images/4/4f/Mv_fig4.png)](/index.php/File:Mv_fig4.png) [](/index.php/File:Mv_fig4.png "Enlarge")Fig.4 - Arrow tail radius | [![](/images/f/f1/Mv_fig5.png)](/index.php/File:Mv_fig5.png) [](/index.php/File:Mv_fig5.png "Enlarge")Fig.5 - Arrow head length | [![](/images/e/e4/Mv_fig6.png)](/index.php/File:Mv_fig6.png) [](/index.php/File:Mv_fig6.png "Enlarge")Fig.6 - Arrow head RGB  
[![](/images/f/fb/Mv_fig7.png)](/index.php/File:Mv_fig7.png) [](/index.php/File:Mv_fig7.png "Enlarge")Fig.7 - Arrow tail RGB | [![](/images/b/bc/Mv_fig8.png)](/index.php/File:Mv_fig8.png) [](/index.php/File:Mv_fig8.png "Enlarge")Fig.8 - Arrow length cutoff | [![](/images/7/7f/Mv_fig9.png)](/index.php/File:Mv_fig9.png) [](/index.php/File:Mv_fig9.png "Enlarge")Fig.9 - Skip every other arrow  
[![](/images/d/d7/Mv_fig10.png)](/index.php/File:Mv_fig10.png) [](/index.php/File:Mv_fig10.png "Enlarge")Fig.10 - Subtract value from vector length | [![](/images/4/42/Mv_fig11.png)](/index.php/File:Mv_fig11.png) [](/index.php/File:Mv_fig11.png "Enlarge")Fig.11 - Draw arrows from CB atoms | [![](/images/b/b1/Mv_fig12.png)](/index.php/File:Mv_fig12.png) [](/index.php/File:Mv_fig12.png "Enlarge")Fig.12 - Shorten by 50%  
| [![](/images/a/a3/Mv_fig13.png)](/index.php/File:Mv_fig13.png) [](/index.php/File:Mv_fig13.png "Enlarge")Fig.13 - Multiple options being used (see final example below)  
      
    
    reinitialize
    import modevectors
    
    fetch 1c3y, async=0
    
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
     
    #This is the default setting (Fig.1)
    modevectors 1c3y_0001, 1c3y_0023
     
    #This is shows that the direction of the arrows is drawn from the 1c3y_0001 towards 1c3y_0023 (Fig.2)
    show cartoon, 1c3y_0023
    color red, 1c3y_0023
     
    #This controls the base radius of the cone/arrow head in angstroms (Fig.3)
    modevectors 1c3y_0001, 1c3y_0023, head=2.5
     
    #This controls the radius of the cylinders/arrow tails in angstroms (Fig.4)
    modevectors 1c3y_0001, 1c3y_0023, tail=0.75
     
    #This controls the length of the cone/arrow head in angstroms (Fig.5)
    #Note that this option does NOT increase the vector length but simply changes the tail length
    modevectors 1c3y_0001, 1c3y_0023, head_length=3.0
     
    #This controls the colour of the cone/arrow head using RGB values (Fig.6)
    modevectors 1c3y_0001, 1c3y_0023, headrgb=(1.0,0.0,0.0)
     
    #This controls the colour of the cylinder/arrow tails using RGB values (Fig.7)
    modevectors 1c3y_0001, 1c3y_0023, tailrgb=(1.0,0.0,0.0)
     
    #This controls the which vectors to show based on a specific cutoff value in angstroms.  Vector lengths that are less 
    #than the cutoff value will not be displayed (Fig.8)
    modevectors 1c3y_0001, 1c3y_0023, cutoff=30.0
     
    #This controls how many vectors to skip (integer value) and is useful when there are too many vectors showing.  
    #Skip=1 will show every other arrow (Fig.9)
    modevectors 1c3y_0001, 1c3y_0023, skip=1
     
    #This controls how much to cut from each vector (in angstroms).  Note that arrows will point in the opposite directions 
    #when too much is cutoff (resulting in a negative vector length) (Fig.10) and should be used with caution!
    modevectors 1c3y_0001, 1c3y_0023, cut=15.0
     
    #This controls which atom to draw a vector from (Fig.11).  Note that this is case-sensitive and is really only useful 
    #when atom=CA or when atom=P (for DNA backbones)
    modevectors 1c3y_0001, 1c3y_0023, atom=CB
     
    #This controls how much to multiple the length of each vector by (percentage increase/decrease) (Fig.12)
    #This example halves the length of each vector (50%)
    modevectors 1c3y_0001, 1c3y_0023, factor=0.5
     
    #This hides the statistics which count atoms skipped, atoms counted (number of arrows showing), and number of atoms 
    #that did not meet the cutoff and are not shown
    modevectors 1c3y_0001, 1c3y_0023, stat=hide
     
    #This example shows multiple options being used (Fig.13)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, tail=1.0, head_length=2.0, headrgb=(1.0,0.0,0.0), tailrgb=(0.0,0.0,1.0),
    cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8
     
    #Finally, this example hides all arrow tails and only uses arrow heads via the notail option(No Figure)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8, notail=1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Modevectors&oldid=13914](https://pymolwiki.org/index.php?title=Modevectors&oldid=13914)"


---

## Monitor file continuously

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script can be used to continuously check the modification timestamp on a file (any format, although this example assumes it's a PDB file), and re-loads it whenever the timestamp changes. As written it is intended to be started from the command line, but this is not a requirement. 

# The Code
    
    
    from pymol import cmd
    import threading
    import time
    import os
    import sys
    
    class pymol_file_monitor (object) :
      def __init__ (self,
          file_name, 
          time_wait=1) : # time in seconds between mtime check
        self.file_name = file_name
        self.time_wait = time_wait
        self.watch = True # this can be toggled elsewhere to stop updating
        self.mtime = 0
        t = threading.Thread(target=self.check_file)
        t.setDaemon(1)
        t.start()
        print "Watching file %s" % file_name
    
      def check_file (self) :
        while (self.watch) :
          if (os.path.exists(self.file_name)) :
            print "checking..."
            mtime = os.path.getmtime(self.file_name)
            if (mtime > self.mtime) :
              self.mtime = mtime
              print "Re-loading %s" % self.file_name
              cmd.load(self.file_name, state=1)
          time.sleep(self.time_wait)
    
    if (__name__ == "pymol") :
      monitor = pymol_file_monitor("status.pdb")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Monitor_file_continuously&oldid=10469](https://pymolwiki.org/index.php?title=Monitor_file_continuously&oldid=10469)"


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

## Mouse modes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Add the code below to the file _mouse_modes.py_ and run it from within pymol. 

After you run it, click on the mouse mode indicator to cycle it at least one time in order to get the new bindings. 

Replace `['three_button_viewing']` by `['three_button_editing']` to edit the **3-Button Editing** mode. 

## mouse_modes.py
    
    
    from pymol.controlling import ring_dict,mode_name_dict,mode_dict
    
    # redefine the three_button_viewing mode
    mode_name_dict['three_button_viewing'] = 'My 3-But View'
    mode_dict['three_button_viewing'] =  [ ('l','none','rota'),
                          ('m','none','move'),
                          ('r','none','movz'),
                          ('l','shft','+Box'),
                          ('m','shft','-Box'),
                          ('r','shft','clip'),                 
                          ('l','ctrl','+/-'),
                          ('m','ctrl','pkat'),
                          ('r','ctrl','pk1'),                 
                          ('l','ctsh','Sele'),
                          ('m','ctsh','orig'),
                          ('r','ctsh','menu'),
                          ('w','none','slab'),
                          ('w','shft','movs'),
                          ('w','ctrl','mvsz'),
                          ('w','ctsh','movz'),
                          ('double_left','none','menu'),
                          ('double_middle','none','none'),
                          ('double_right','none', 'pk1'),
                          ('single_left','none','+/-'),
                          ('single_middle','none','cent'),
                          ('single_right','none', 'pkat'),
                          ]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mouse_modes&oldid=8361](https://pymolwiki.org/index.php?title=Mouse_modes&oldid=8361)"


---

## Movie color fade

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/movie_color_fade.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/movie_color_fade.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![movie_color_fade example](/images/0/05/Movie_color_fade_example_1.gif)](/index.php/File:Movie_color_fade_example_1.gif "movie_color_fade example")

**movie_color_fade** is like [movie_fade](/index.php/Movie_fade "Movie fade"), but will fade colors in a movie. Simply specify the arguments (see below). 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Examples
  * 3 SEE ALSO



## Usage
    
    
     movie_color_fade [ startframe [, startcolor [, endframe [, endcolor [, selection ]]]]]
     help movie_color_fade
    

|   
---|---  
  
### Arguments

  * startframe, endframe = beginning and end movie frame for fading 
    * if omitted, the current or last frame will be respectively set automatically
  * startcolor, endcolor = coloring at start and end
  * selection: target selection



## Examples

This example yields the movie to the top right 
    
    
    # make function available to PyMOL
    import movie_color_fade
    
    # object prep
    bg black
    delete all
    fetch 1hpv, async=0
    as cartoon
    orient
    color yellow, chain A
    color skyblue, chain B
    # movie prep
    mset 1x120
    
    # color fading
    movie_color_fade 1, yellow, 60, skyblue, chain A
    movie_color_fade 60, skyblue, 120, yellow, chain A
    movie_color_fade 1, skyblue, 60, yellow, chain B
    movie_color_fade 60, yellow, 120, skyblue, chain B
    

|   
---|---  
  
  


## SEE ALSO

  * [Movie_fade](/index.php/Movie_fade "Movie fade")
  * [mappend](/index.php/Mappend "Mappend")
  * [Color](/index.php/Color "Color")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_color_fade&oldid=13945](https://pymolwiki.org/index.php?title=Movie_color_fade&oldid=13945)"


---

## Movie fade

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/movie_fade.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/movie_fade.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will help fade in and out settings in a movie. Just specify the setting, it's initial value at an initial frame and it's ending value and frame. 

## Usage
    
    
    movie_fade setting, startFrame, startVal, endFrame, endVal [, selection ]
    

## Examples

To fade in sticks from fully transparent to fully opaque across 60 frames do: 
    
    
    mset 1x60
    movie_fade stick_transparency, 1, 1., 60, 0.
    

More complex example which involves camera motion: 
    
    
    fetch 1rx1, async=0
    as cartoon
    show surface
    mset 1x80
    movie.roll
    movie_fade transparency,  1, 0., 40, 1.
    movie_fade transparency, 41, 1., 80, 0.
    

**Example file:**

[![](/images/9/98/Movie_fade_example.gif)](/index.php/File:Movie_fade_example.gif)

[](/index.php/File:Movie_fade_example.gif "Enlarge")

Result
    
    
    import movie_fade
    
    fetch 1hpv, async=0
    orient
    
    #format
    bg black
    show_as cartoon
    color marine, ss s
    color red, ss h
    color white, ss l+""
    show sticks, resn 478
    util.cbao resn 478
    set surface_color, white, 1hpv
    show surface
    
    #movie
    mset 1x120
    movie_fade transparency,1, 0, 60, 1, 1hpv
    movie_fade transparency,61, 1, 120, 0, 1hpv
    

## See Also

  * [mdo](/index.php/Mdo "Mdo")
  * [mappend](/index.php/Mappend "Mappend")
  * [set](/index.php/Set "Set")
  * [movie_color_fade](/index.php/Movie_color_fade "Movie color fade")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_fade&oldid=13946](https://pymolwiki.org/index.php?title=Movie_fade&oldid=13946)"


---

## Movit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Hi everyone ! 

I wanted to see roughly if one protein could dock on another one, and could not find any mouse setting allowing to do that. Here is a script that does just that. 

USAGE : movit <selection>

The selection can now be translated/rotated with the following keyboard shortcuts: 

q,a : x translation 

w,s : y translation 

e,d : z translation 

  
r,f : x rotation 

t,g : y rotation 

y,h : z rotation 

  


Use : run movit.py during startup. Then, 'movit _selection'_. 

It is not very practical, does anybody know how to do that with the mouse instead of just one step at a time like here ? 
    
    
    from pymol import cmd
    
    def movit(selection="all"):
        """USAGE : movit <selection>
                   The selection can now be translated/rotated
                   with the following keyboard shortcuts:
    
                   q,a : x translation
                   w,s : y translation
                   e,d : z translation
    
                   r,f : x rotation
                   t,g : y rotation
                   y,h : z rotation
        """
    
        line="translate [1,0,0],"+selection
        cmd.alias ('q',line)
        line="translate [0,1,0],"+selection
        cmd.alias ('w',line)
        line="translate [0,0,1],"+selection
        cmd.alias ('e',line)
    
        line="translate [-1,0,0],"+selection
        cmd.alias ('a',line)
        line="translate [0,-1,0],"+selection
        cmd.alias ('s',line)
        line="translate [0,0,-1],"+selection
        cmd.alias ('d',line)
    
        line="rotate x,5,"+selection
        cmd.alias ('r',line)
        line="rotate y,5,"+selection    
        cmd.alias ('t',line)
        line="rotate z,5,"+selection
        cmd.alias ('y',line)
    
        line="rotate x,-5,"+selection
        cmd.alias ('f',line)
        line="rotate y,-5,"+selection
        cmd.alias ('g',line)
        line="rotate z,-5,"+selection
        cmd.alias ('h',line)
    
    cmd.extend ("movit",movit)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Movit&oldid=7119](https://pymolwiki.org/index.php?title=Movit&oldid=7119)"


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

## ObjectByArrows

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

When sifting through many similar structures for days on end, mousing around the object pane can become strenuous. This script binds the left and right arrows on your keyboard to moving up and down the object list. Note that this script supports objects, but is untested for groups. 
    
    
    from pymol import cmd
    
    def move_down():
            enabled_objs = cmd.get_names("object",enabled_only=1)
            all_objs = cmd.get_names("object",enabled_only=0)
            for obj in enabled_objs:
                    cmd.disable(obj)
                    last_obj=obj
                    for i in range(0,len(all_objs)):
                            if all_objs[i] == obj:
                                    if i+1 >= len(all_objs):
                                            cmd.enable( all_objs[0] )
                                    else:
                                            cmd.enable( all_objs[i+1] )
            cmd.orient
    def move_up():
            enabled_objs = cmd.get_names("object",enabled_only=1)
            all_objs = cmd.get_names("object",enabled_only=0)
            for obj in enabled_objs:
                    cmd.disable(obj)
                    last_obj=obj
                    for i in range(0,len(all_objs)):
                            if all_objs[i] == obj:
                                    if i-1 < 0:
                                            cmd.enable( all_objs[-1] )
                                    else:
                                            cmd.enable( all_objs[i-1] )
            cmd.orient
    
    cmd.set_key('left', move_up)
    cmd.set_key('right', move_down)
    

Retrieved from "[https://pymolwiki.org/index.php?title=ObjectByArrows&oldid=6332](https://pymolwiki.org/index.php?title=ObjectByArrows&oldid=6332)"


---

## ObjectFocus

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

A simple script to walk over objects in the object menu panel. 

# The Code
    
    
    import pymol
    from pymol import cmd
    
    UP, DOWN = -1, 1
    
    in_focus = 0
    
    def object_focus(direction):
        """
        A simple script that remaps the PgUp and PgDn keys to walk through
        objects in the object list.
        """
    
        global in_focus, UP, DOWN
    
        cmd.wizard()
    
        names = cmd.get_names("public_objects")
    
        in_focus += direction
    
        if in_focus<0:
            in_focus = 0
        if in_focus > len(names)-1:
            in_focus = len(names)-1
    
        cur_obj = names[in_focus]
    
        cmd.orient(cur_obj,animate=1)
        cmd.wizard("message", "Object: %s" % (cur_obj))
    
    
    cmd.set_key("PgUp", object_focus, [UP])
    cmd.set_key("PgDn", object_focus, [DOWN])
    

Retrieved from "[https://pymolwiki.org/index.php?title=ObjectFocus&oldid=10848](https://pymolwiki.org/index.php?title=ObjectFocus&oldid=10848)"


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

## PDB Web Services Script

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

A short snippet of code utilizing the new PDB Web Services. 

Notes: 

  1. See [PDB Web Services](http://www.rcsb.org/robohelp/webservices/summary.htm)
  2. You need [SOAP for Python](http://pywebsvcs.sourceforge.net/) installed
  3. The sequence of a chain fetched from the PDB does not always equal the sequence of a chain fetched from the PDB Web Services. I believe the PDB Web Services has the complete biological chain, not just what's in the structure.
  4. The offset is a hack; it uses a trick in PyMOL and may not be robust at all. Use at your own caution.



## The Code
    
    
    import pymol
    from pymol import stored, cmd
    import SOAPpy
    
    # define a mapping from three letters to one
    # Thanks to whomever posted this on the PyMOLWiki:
    #   http://www.pymolwiki.org/index.php/Label
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}
    
    # fetch some molecule; yes, 1foo exists in the PDB
    cmd.fetch("1foo", async=0)
    # get the sequence from PyMOL
    stored.pyMOLSeq = []
    cmd.iterate( "1foo and c. A and n. CA" , "stored.pyMOLSeq.append( one_letter[resn] )" )
    
    # open up the PDB web services and make a connection
    stored.pyMOLSeq = ''.join( stored.pyMOLSeq )
    server = SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")
    
    # fetch the secondary structure assignment from the PDB &
    # split it up into an array of characters
    stored.ss = list( server.getKabschSander("1foo", "A") )
    
    # get the sequence
    pdbSeq = server.getSequenceForStructureAndChain("1foo", "A")
    # danger: hackishly foolish, but worked for 1foo.
    offset = pdbSeq.find( stored.pyMOLSeq )
    
    # make the assignment in PyMOL
    cmd.alter( "1foo and c. A and n. CA and i. " + str(offset) + "-", "ss=stored.ss.pop(0)" )
    
    # show as cartoons
    cmd.as("cartoon")
    

Retrieved from "[https://pymolwiki.org/index.php?title=PDB_Web_Services_Script&oldid=7993](https://pymolwiki.org/index.php?title=PDB_Web_Services_Script&oldid=7993)"


---

## Pdbsurvey

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Hi everyone ! 

I have created a script because I was tired of browsing the pdb for endless searches of structures relevant to my field. This scripts takes as input a text file in which you copied your favourite keywords, and the number of days you want to search back (by default, it will look at the structures added the last 50 days). It generates a report text file that contains the pdb id and name of the relevant structures that have been recently added. All you need to do is add 'run pdbsurvey.py' to your startup file, and create a text file called 'keywords.txt' with your keywords separated by an end-of-line character. Then you're ready to go. Just hit 'pdbsurvey' from within your PyMol instance, and the program returns the report file. The pdb ftp server is updated weekly. 
    
    
    from pymol import cmd
    
    def pdbsurvey(days=50):
    
        """USAGE : pdbsurvey (<days>)
        Surveys the updates added to the PDB (ftp.rcsb.org) in the last
        50 days (or otherwise specified when calling this function) for
        entries that contain the words specified in the file
        keywords.txt.
        """
        print days
    
        import ftplib
        import time
        import os
    
    
    
        def todaymerge():
            """Puts today's date in a pdb format string.
            """
            date=time.localtime()
            fyear="%i" %(date[0])
            fmonth="%i" %(date[1])
            if date[1]<10:
                fmonth="0"+"%i" %(date[1])
            fday="%i" %(date[2])
            if date[2]<10:
                fday="0"+"%i" %(date[2])
            dateS=fyear+fmonth+fday
            return dateS
    
        def file2list(filename):
            """Low-level routine to brainlessly implement
            file.read().
            """
            fq=open(filename,'rb')
            linesS=fq.read()
            fq.close()
            LIST=linesS.splitlines()
            return LIST
    
        def connect2pdb():
            """Opens an anonymous socket to ftp://ftp.rcsb.org
            """
            f=ftplib.FTP()
            f.connect ('ftp.rcsb.org')
            f.login ()
            print "Remote connection established","\n"
            return f
    
        def decrementdate(dateS):
            """given a string date (pdb format yyyymmdd)
            this routine returns a string of the day before
            (sadly assuming that every month has 31 days, but
            no big deal here...).
            """
            #decompose dateS into components
            yearS=dateS[0]+dateS[1]+dateS[2]+dateS[3]
            monthS=dateS[4]+dateS[5]
            dayS=dateS[6]+dateS[7]
    
            #convert each into integers
            yearI=int(yearS)
            monthI=int(monthS)
            dayI=int(dayS)
    
            #actual stuff
            dayI=dayI-1
            if dayI==0:
                dayI=31
                monthI-=1
                if monthI==0:
                    monthI=12
                    yearI-=1
            dayS="%i" %(dayI)
            monthS="%i" %(monthI)
            yearS="%i" %(yearI)
            if dayI<10:
                dayS="0"+dayS
            if monthI<10:
                monthS="0"+monthS
            #and finally...
            dateS=yearS+monthS+dayS
            return dateS
    
        def findlastdir(dateS,f,days):
            """Puts the names of the "recent" directories in the
            list named "directoriesL".
            """
            directoriesL=[]
            for p in range(days):
                dateS=decrementdate(dateS)
                attempt="/pub/pdb/data/status/"+dateS
                try:
                    f.cwd(attempt)
                    directoriesL.append(attempt)
                except:
                    pass
            return directoriesL
    
        def compilinfile(directoriesL,f):
            """lists all structures in the added.pdb files
            contained in the directories specified in directoriesL
            """
            command="RETR added.pdb"
            handle=open("donotedit.dat","wrb")
            for k in directoriesL:
                f.cwd(k)
                print "Currently in directory ",f.pwd()
                f.retrbinary(command,handle.write)
            handle.close()
            return len(directoriesL)
    
        def listparser():
            """Extracts the pdbids from donotedit.dat file,
            and stacks them into the list pdbidsL
            """
            linesL=file2list("donotedit.dat")
            pdbidsL=[]
            for iter in linesL:
                pdbidsL.append(iter[57:61])
            pdbidsL.sort()
            return pdbidsL
    
        def currentrelease(f):
            """Stores the content of cmpd_res.idx file
            This file contains the equivalencies pdbid<->title
            for all current entries of the PDB.
            """
            command="RETR cmpd_res.idx"
            f.cwd("/pub/pdb/derived_data/index/")
            print "Currently in directory ",f.pwd()
            fq=open("dictionnary.dat",'wrb')
            f.retrbinary(command,fq.write)
            fq.close()
            dictL=file2list("dictionnary.dat")
            return dictL
    
        def extract(pdbidsL,dictL):
            """Populates dictionnaryD with pdb entries found in the
            latest releases.
            """
            dictionnaryD={}
            problemL=[]
            extractL=[dictionnaryD,problemL]
            for i in dictL:
                tempS=i[0:4].lower()
                for ii in pdbidsL:
                    if ii == tempS:
                        title=i[14:216]
                        extractL[0][ii]=title
            if len(extractL[0]) != len(pdbidsL):
                print "Dimension mismatch, seeking troublemaker..."
                for i in pdbidsL:
                    if not any(i == ii for ii in extractL[0]):
                        extractL[1].append(i)
            return extractL
    
        def disconnectpdb(f):
            """Diconnects the current ftp session
            """
            f.quit()
            print "Remote connection terminated","\n"
            return f
    
        def releventries(dictionnaryD):
            """Generates a cleaned dictionary with only entries
            that have one or more keywords specified in the local
            user-defined keywords.txt file
            """
            keywL=file2list("keywords.txt")
            relevdicD={}
            for i in keywL:
                for elem,temp in dictionnaryD.items():
                    if i in temp:
                        relevdicD[elem]=temp
            return relevdicD
    
        def diskcleanup(filelist=["donotedit.dat","dictionnary.dat"]):
            """Lo-level disk cleanup to free up memory without the user
            """
            for filename in filelist:
                command='DEL '+filename
                os.system(command)
            return "clean"
    
    
    
    
        print "Welcome in the auto-PDB updater !"
    
        print "Survey of updates made since",days,"days ago."
        
        print "Acquisition of local time..."
        dateS=todaymerge()                                                 #Initializes dateS
        print "today is ",dateS
        print "Connecting to remote ftp server..."
        f=connect2pdb()                                                    #Connect anonymously to ftp.rcsb.org
    
        print "Acquisition of latest added remote directories..."
        directoriesL=findlastdir(dateS,f,days)                             #Lists recent directories in directoriesL
        if len(directoriesL)==0:
            print "No updates have been found since",days,"ago. Starting over with 50 days ago."
            directoriesL=findlastdir(dateS,f,50)
    
        print "Acquisition of latest addedremote files..."
        updatesnumberI=compilinfile(directoriesL,f)                        #Concatenates the corresponding added.pdb into donotedit.dat
    
        print "Parsing of latest entries..."
        pdbidsL=listparser()                                               #Recent names now present in the pdbidsL list (one name per element)
    
        print "Acquisition of the current pdb distribution..."
        dictL=currentrelease(f)                                            #Populates dictL with the current entries of the PDB
    
        print "Parsing of the current pdb distribution into [code,title] tuples..."
        extractL=extract(pdbidsL,dictL)                                    #generates the dictionnary of latest releases key:PDBid ; definition:pdbtitle
    
        print "Disconnection from the remote ftp server..."
        f=disconnectpdb(f)                                                 #Closes the ftp instance
    
        print "Extraction of the relevant entries..."
        relevdicD=releventries(extractL[0])                               #Generates a subset of dictionnary D with criterion being "has keywords contained in keywords.txt in its title"
    
        print "Cleaning program-generated temporary files..."
        clean=diskcleanup()                                                #Cleans the mess generated by the program
    
        reportL=[]
        reportL.append("\n")
        reportL.append("###############REPORT########################################\n")
        reportL.append("\n")
        lendictS="%i" %(len(dictL))
        chmilblik = 'The current pdb version (as of '+dateS+") has "+lendictS+" entries.\n"
        reportL.append(chmilblik)
        line="The most recent directory is : "+directoriesL[0]+".\n"
        reportL.append(line)
        updatesnumberS="%i" %(updatesnumberI)
        entriesnumber="%i" %(len(extractL[0]))
        line="The "+updatesnumberS+" last updates ("+entriesnumber+" entries) have been examined.\n"
        reportL.append(line)
        diclengthS="%i" %(len(relevdicD))
        line=diclengthS+" are relevant to you :\n"
        reportL.append(line)
        for i,x in relevdicD.items():
            entry=i+" : "+x+"\n"
            reportL.append(entry)
        problemS=""
        for i in extractL[1]:
            problemS=i+";"+problemS
        problemS="["+problemS
        problemS=problemS.strip(";")
        problemS=problemS+"]"
        lineS="The entries "+problemS+" raised problems,"
        reportL.append(lineS)
        reportL.append("they should be examined manually.")
        reportL.append("\n")
        reportL.append("###############END OF REPORT#################################\n")
        report=open("report.aut","w")
        for elem in reportL:
            print elem
            report.writelines(elem+'\n')
        report.close()
        command2='start keywords.txt'
        command3='start report.aut'
        os.system(command2)
        os.system(command3)
    
    cmd.extend("pdbsurvey",pdbsurvey)
    

Thank you for any feedback, ways of improving it,... 

Retrieved from "[https://pymolwiki.org/index.php?title=Pdbsurvey&oldid=8010](https://pymolwiki.org/index.php?title=Pdbsurvey&oldid=8010)"


---

## Perp maker

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/perp_maker.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/perp_maker.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The **perp_maker.py** script creates perpendicular planes. 

Nothing to do with cops. Given a simple PyMol scene, attempts to create a CGO background triangle perpendicular to the vector created - which is parallel to the line segment drawn through the camera point and current center of mass - as obtained by "get_position," or "get_view." 

## Usage

  * Load your scene
  * Orient the scene as you wish
  * [Run](/index.php/Running_Scripts "Running Scripts") the script



Could it be any simpler? 

You can rotate and move the plane using the editing features (A > drag, Shift+Drag with right, middle or left mouse button). 

Retrieved from "[https://pymolwiki.org/index.php?title=Perp_maker&oldid=13916](https://pymolwiki.org/index.php?title=Perp_maker&oldid=13916)"


---

## Plane Wizard

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Introduction
  * 2 Original
  * 3 Modification
    * 3.1 Gallery
    * 3.2 plane.py
    * 3.3 Examples



## Introduction

This wizard has a simple purpose - to draw a cgo plane that passes through three points picked by the user. Most of the wizard itself was copied from the measure wizard. 

To use, put it in the same directory as the other wizards. This is not quality code, and there may be bugs, but it seems to work okay. 

## Original
    
    
    import pymol
    from pymol import cmd
    from pymol.wizard import Wizard
    from chempy import cpv
    from pymol.cgo import *
    
    def makePrimitive(cgo, name):
        az = cmd.get('auto_zoom', quiet=1)
        cmd.set('auto_zoom', 0, quiet=1)
        cmd.load_cgo(cgo, name)
        cmd.set('auto_zoom', az, quiet=1)
    
    def point(p):
        x, y, z = p
        return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]
    
    def line(p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]
    
    def plane(corner1, corner2, corner3, corner4, normal):
        planeObj = []
        planeObj.extend(point(corner1))
        planeObj.extend(point(corner2))
        planeObj.extend(point(corner3))
        planeObj.extend(point(corner4))
        planeObj.extend(line(corner1, corner2))
        planeObj.extend(line(corner2, corner3))
        planeObj.extend(line(corner3, corner4))
        planeObj.extend(line(corner4, corner1))
    
        planeObj.extend([COLOR, 0.8, 0.8, 0.8])
        planeObj.extend([BEGIN, TRIANGLE_STRIP])
        planeObj.append(NORMAL)
        planeObj.extend(normal)
        for corner in [corner1, corner2, corner3, corner4, corner1]:
            planeObj.append(VERTEX)
            planeObj.extend(corner)
        planeObj.append(END)
        return planeObj
    
    def planeFromPoints(point1, point2, point3, facetSize):
        v1 = cpv.normalize(cpv.sub(point2, point1))
        v2 = cpv.normalize(cpv.sub(point3, point1))
        normal = cpv.cross_product(v1, v2)
        v2 = cpv.cross_product(normal, v1)
        x = cpv.scale(v1, facetSize)
        y = cpv.scale(v2, facetSize)
        center = point2
        corner1 = cpv.add(cpv.add(center, x), y)
        corner2 = cpv.sub(cpv.add(center, x), y)
        corner3 = cpv.sub(cpv.sub(center, x), y)
        corner4 = cpv.add(cpv.sub(center, x), y)
        return plane(corner1, corner2, corner3, corner4, normal)
    
    
    class PlaneWizard(Wizard):
    
        def __init__(self):
            Wizard.__init__(self)
    
            # some attributes to do with picking
            self.pick_count = 0
            self.object_count = 0
            self.object_prefix = "pw"
    
            # the plane facet size (the 'radius' of the section of plane we show)
            self.facetSize = 5
    
            self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
            cmd.set("mouse_selection_mode",0) # set selection mode to atomic
            cmd.deselect()
    
        def reset(self):
            cmd.delete(self.object_prefix + "*")
            cmd.delete("sele*")
            cmd.delete("_indicate*")
            cmd.unpick()
            cmd.refresh_wizard()
    
        def delete_all(self):
            cmd.delete("plane*")
    
        def cleanup(self):
            cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
            self.reset()
            self.delete_all()
    
        def get_prompt(self):
            self.prompt = None
            if self.pick_count == 0:
                self.prompt = [ 'Please click on the first atom...']
            elif self.pick_count == 1:
                self.prompt = [ 'Please click on the second atom...' ]
            elif self.pick_count == 2:
                self.prompt = [ 'Please click on the third atom...' ]
            return self.prompt
    
        def do_select(self, name):
            # "edit" only this atom, and not others with the object prefix
            try:
                cmd.edit("%s and not %s*" % (name, self.object_prefix))
                self.do_pick(0)
            except pymol.CmdException, pmce:
                print pmce
    
        def pickNextAtom(self, atom_name):
            # transfer the click selection to a named selection
            cmd.select(atom_name, "(pk1)")
    
            # delete the click selection
            cmd.unpick()
    
            # using the magic of indicate, highlight stuff
            indicate_selection = "_indicate" + self.object_prefix
            cmd.select(indicate_selection, atom_name)
            cmd.enable(indicate_selection)
    
            self.pick_count += 1
            self.error = None
    
            # necessary to force update of the prompt
            cmd.refresh_wizard()
    
        def do_pick(self, picked_bond):
    
            # this shouldn't actually happen if going through the "do_select"
            if picked_bond:
                self.error = "Error: please select bonds, not atoms"
                print self.error
                return
    
            atom_name = self.object_prefix + str(self.pick_count)
            if self.pick_count < 2:
                self.pickNextAtom(atom_name)
            else:
                self.pickNextAtom(atom_name)
    
                point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
                point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
                point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
                plane = planeFromPoints(point1, point2, point3, self.facetSize)
    
                planeName = "plane-%02d" % self.object_count
                self.object_count += 1
                makePrimitive(plane, planeName)
                cmd.show("cgo", "plane*")
    
                self.pick_count = 0
                self.reset()
    
        def get_panel(self):
            return [
                [ 1, 'Plane Wizard',''],
                [ 2, 'Reset','cmd.get_wizard().reset()'],
                [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
                [ 2, 'Done','cmd.set_wizard()'],
            ]
    
    # create an instance
    
    wiz = PlaneWizard()
    
    # make this the active wizard
    
    cmd.set_wizard(wiz)
    

## Modification

Small modifications to the same code as above. 

### Gallery

  * [![Drawing a plane in a protein](/images/2/28/Plane_img1.png)](/index.php/File:Plane_img1.png "Drawing a plane in a protein")

Drawing a plane in a protein 

  * [![Drawing 6 planes in a box](/images/d/d3/Plane_img2.png)](/index.php/File:Plane_img2.png "Drawing 6 planes in a box")

Drawing 6 planes in a box 

  * [![Drawing 6 planes in a box](/images/2/25/Plane_img3.png)](/index.php/File:Plane_img3.png "Drawing 6 planes in a box")

Drawing 6 planes in a box 




### plane.py

Make a **plane.py** file in the same directory where you are working 
    
    
    '''
    Described at PyMOL wiki:
    https://pymolwiki.org/index.php/Plane_Wizard
    
    Authors : Troels Schwarz-Linnet
    Date    : Dec 2016
    Modified: From previous contributors. 
    '''
    
    import pymol
    from pymol import cmd
    from pymol.wizard import Wizard
    from chempy import cpv
    from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, TRIANGLE_STRIP, NORMAL, VERTEX, END, ALPHA
    
    def makePrimitive(cgo, name):
        az = cmd.get('auto_zoom', quiet=1)
        cmd.set('auto_zoom', 0, quiet=1)
        cmd.load_cgo(cgo, name)
        cmd.set('auto_zoom', az, quiet=1)
    
    def point(p):
        x, y, z = p
        return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]
    
    
    def line(p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]
    
    
    def plane(corner1, corner2, corner3, corner4, normal, settings):
        planeObj = []
        planeObj.extend(point(corner1))
        planeObj.extend(point(corner2))
        planeObj.extend(point(corner3))
        planeObj.extend(point(corner4))
        planeObj.extend(line(corner1, corner2))
        planeObj.extend(line(corner2, corner3))
        planeObj.extend(line(corner3, corner4))
        planeObj.extend(line(corner4, corner1))
    
        # Make settings
        if 'ALPHA' in settings:
            planeObj.extend([ALPHA, settings['ALPHA']])
    
        if 'COLOR' in settings:
            planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
        else:
            planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish
    
        planeObj.extend([BEGIN, TRIANGLE_STRIP])
        planeObj.append(NORMAL)
    
        if 'INVERT' in settings:
            if settings['INVERT']==True:
                planeObj.extend(cpv.negate(normal))
            else:
                planeObj.extend(normal)
        else:
            planeObj.extend(normal)
    
    
        for corner in [corner1, corner2, corner3, corner4, corner1]:
            planeObj.append(VERTEX)
            planeObj.extend(corner)
        planeObj.append(END)
        return planeObj
    
    
    def planeFromPoints(p1, p2, p3, vm1=None, vm2=None, center=True, settings={}):
        v1 = cpv.sub(p1, p2)
        v2 = cpv.sub(p3, p2)
        normal = cpv.cross_product(v1, v2)
    
        if 'translate' in settings:
            vtran = cpv.scale(cpv.normalize(normal), settings['translate'])
            p1_t = cpv.sub(p1, vtran)
            p2_t = cpv.sub(p2, vtran)
            p3_t = cpv.sub(p3, vtran)
            print("New coordinates are:")
            print_info("New", p1_t, p2_t, p3_t)
            print("New coordinates are for normalized plane:")
            v1_t = cpv.normalize(cpv.sub(p1_t, p2_t))
            v2_t = cpv.normalize(cpv.sub(p3_t, p2_t))
            normal_t = cpv.normalize(cpv.cross_product(v1_t, v2_t))
            v2_t = cpv.normalize(cpv.cross_product(normal_t, v1_t))
            p1_t2 = cpv.add(v1_t, p2_t)
            p3_t2 = cpv.add(v2_t, p2_t)
            print_info("Newnormal", p1_t2, p2_t, p3_t2)
    
        if vm1!=None:
            v1 = cpv.scale(cpv.normalize(v1), vm1)
        if vm2!=None:
            v2 = cpv.scale(cpv.normalize(v2), vm2)
    
        centrum = p2
        if center:
            corner1 = cpv.add(cpv.add(centrum, v1), v2)
            corner2 = cpv.sub(cpv.add(centrum, v1), v2)
            corner3 = cpv.sub(cpv.sub(centrum, v1), v2)
            corner4 = cpv.add(cpv.sub(centrum, v1), v2)
        else:
            corner1 = cpv.add(cpv.add(centrum, v1), v2)
            corner2 = cpv.add(centrum, v1)
            corner3 = centrum
            corner4 = cpv.add(centrum, v2)
    
        return plane(corner1, corner2, corner3, corner4, normal, settings)
    
    
    def print_info(name, coor1, coor2, coor3):
        cs1 = (map(float, [ '%.2f' % elem for elem in coor1 ]) )
        cs2 = (map(float, [ '%.2f' % elem for elem in coor2 ]) )
        cs3 = (map(float, [ '%.2f' % elem for elem in coor3 ]) )
        print("You can also use the function calls with these coordinates")
        print("plane.make_plane_points(name='%s', l1=%s, l2=%s, l3=%s)"%(name, cs1, cs2, cs3))
    
    
    def make_plane(name,a1='(pk1)',a2='(pk2)',a3='(pk3)', vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
        """
        DESCRIPTION
        Create a CGO plane from three atomic coordinates
    
        USAGE
        make_plane name, a1, a2, a3
    
        where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
        """
        # get coordinates from atom selections
        coor1 = cmd.get_model(a1).get_coord_list()[0]
        coor2 = cmd.get_model(a2).get_coord_list()[0]
        coor3 = cmd.get_model(a3).get_coord_list()[0]
    
        # Help with alternative
        print_info(name, coor1, coor2, coor3)
    
        # Get the plane
        plane = planeFromPoints(p1=coor1, p2=coor2, p3=coor3, vm1=vm1, vm2=vm2, center=center, settings=settings)
        makePrimitive(plane, name)
        #cmd.show("cgo", "plane*")
    
        if makepseudo:
            cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=coor1)
            cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=coor2)
            cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=coor3)
    
    # Extend function to be called inside pymol
    cmd.extend("make_plane", make_plane)
    
    def make_plane_points(name,l1=None,l2=None,l3=None, vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
        """
        DESCRIPTION
        Create a CGO plane from three atomic coordinates
    
        USAGE
        make_plane name, l1, l2, l3
    
        where each xys is a list with floats of x,y,z coordinates
        """
        if l1==None or l2==None or l3==None:
            print("Please provide a list of xyz floats for each 3 positions")
            return
        if type(l1) is not list or type(l2) is not list or type(l3) is not list:
            print(type(l1),type(l2),type(l3))
            print("Please provide 3 list of xyz floats for each 3 positions")
            return
    
        plane = planeFromPoints(p1=l1, p2=l2, p3=l3, vm1=vm1, vm2=vm2, center=center, settings=settings)
        makePrimitive(plane, name)
    
        if makepseudo:
            cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=l1)
            cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=l2)
            cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=l3)
    
    # Extend function to be called inside pymol
    cmd.extend("make_plane_points", make_plane_points)
    
    class PlaneWizard(Wizard):
        def __init__(self):
            Wizard.__init__(self)
    
            # some attributes to do with picking
            self.pick_count = 0
            self.object_count = 0
            self.object_prefix = "pw"
    
            self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
            cmd.set("mouse_selection_mode",0) # set selection mode to atomic
            cmd.deselect()
    
        def reset(self):
            cmd.delete(self.object_prefix + "*")
            cmd.delete("sele*")
            cmd.delete("_indicate*")
            cmd.unpick()
            cmd.refresh_wizard()
    
        def delete_all(self):
            cmd.delete("plane*")
    
        def cleanup(self):
            cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
            self.reset()
            self.delete_all()
    
        def get_prompt(self):
            self.prompt = None
            if self.pick_count == 0:
                self.prompt = [ 'Please click on the first atom...']
            elif self.pick_count == 1:
                self.prompt = [ 'Please click on the second atom...' ]
            elif self.pick_count == 2:
                self.prompt = [ 'Please click on the third atom...' ]
            return self.prompt
    
        def do_select(self, name):
            # "edit" only this atom, and not others with the object prefix
            try:
                cmd.edit("%s and not %s*" % (name, self.object_prefix))
                self.do_pick(0)
            except pymol.CmdException, pmce:
                print pmce
    
        def pickNextAtom(self, atom_name):
            # transfer the click selection to a named selection
            cmd.select(atom_name, "(pk1)")
    
            # delete the click selection
            cmd.unpick()
    
            # using the magic of indicate, highlight stuff
            indicate_selection = "_indicate" + self.object_prefix
            cmd.select(indicate_selection, atom_name)
            cmd.enable(indicate_selection)
    
            self.pick_count += 1
            self.error = None
    
            # necessary to force update of the prompt
            cmd.refresh_wizard()
    
        def do_pick(self, picked_bond):
    
            # this shouldn't actually happen if going through the "do_select"
            if picked_bond:
                self.error = "Error: please select bonds, not atoms"
                print self.error
                return
    
            atom_name = self.object_prefix + str(self.pick_count)
            if self.pick_count < 2:
                self.pickNextAtom(atom_name)
            else:
                self.pickNextAtom(atom_name)
    
                point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
                point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
                point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
                plane = planeFromPoints(point1, point2, point3)
    
                planeName = "plane-%02d" % self.object_count
    
                print_info(planeName, point1, point2, point3)
    
                self.object_count += 1
                makePrimitive(plane, planeName)
                cmd.show("cgo", "plane*")
    
                self.pick_count = 0
                self.reset()
    
        def get_panel(self):
            return [
                [ 1, 'Plane Wizard',''],
                [ 2, 'Reset','cmd.get_wizard().reset()'],
                [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
                [ 2, 'Done','cmd.set_wizard()'],
            ]
    

### Examples

Plane in a protein 
    
    
    reinitialize
    
    fetch 1ubq, async=0
    show_as cartoon, all 
    
    import plane
    # Start the wizard, and do manual picking
    cmd.set_wizard(plane.PlaneWizard())
    

Or by selection 
    
    
    reinitialize
    
    fetch 1ubq, async=0
    show_as cartoon, all 
    
    import plane
    
    # Set alpha level and color
    #violetpurple: 0.55, 0.25, 0.60 #yellow: 1.0, 1.0, 0.0, #blue: 0.0, 0.0, 1.0 #orange: 1.0, 0.5, 0.0 #forest: 0.2, 0.6, 0.2 #red: 1.0, 0.0, 0.0
    dict = {'ALPHA':0.6, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':True}
    plane.make_plane(name='test', a1='/1ubq//A/24/CA', a2='/1ubq//A/29/CA', a3='/1ubq//A/40/CA', center=False, settings=dict)
    

Or by atom coordinates 
    
    
    reinitialize
    
    fetch 1ubq, async=0
    show_as cartoon, all 
    
    import plane
    
    # Or  make from atom coordinates
    #plane.make_plane_points(name='test', l1=[35.03, 21.72, 17.07], l2=[37.47, 27.39, 10.67], l3=[37.74, 31.64, 23.71])
    
    # Define plane, 10 angstrom in length
    #plane.make_plane_points(name='p1', l1=[0.0, 10.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 0.0, 10.0], center=False, makepseudo=False)
    

Or make a color cube 
    
    
    reinitialize
    #violetpurple: 0.55, 0.25, 0.60 #yellow: 1.0, 1.0, 0.0, #blue: 0.0, 0.0, 1.0 #orange: 1.0, 0.5, 0.0 #forest: 0.2, 0.6, 0.2 #red: 1.0, 0.0, 0.0
    
    # YZ Plane, #purple
    dict = {'ALPHA':0.4, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':True}
    plane.make_plane_points(name='p1', l1=[0.0, 10.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 0.0, 10.0], center=False, makepseudo=False, settings=dict)
    # YZ Plane, shifted in X, #yellow
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 1.0, 0.0]}
    plane.make_plane_points(name='p6', l1=[10.0, 10.0, 0.0], l2=[10.0, 0.0, 0.0], l3=[10.0, 0.0, 10.0], center=False, makepseudo=False, settings=dict)
    
    # XZ Plane, blue
    dict = {'ALPHA':0.4, 'COLOR':[0.0, 0.0, 1.0]}
    plane.make_plane_points(name='p2', l1=[10.0, 0.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 0.0, 10.0], center=False, makepseudo=False, settings=dict)
    # XZ Plane, shifted in Y, #orange
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.5, 0.0], 'INVERT':True}
    plane.make_plane_points(name='p5', l1=[10.0, 10.0, 0.0], l2=[0.0, 10.0, 0.0], l3=[0.0, 10.0, 10.0], center=False, makepseudo=False, settings=dict)
    
    # XY Plane, forest
    dict = {'ALPHA':0.4, 'COLOR':[0.2, 0.6, 0.2], 'INVERT':True}
    plane.make_plane_points(name='p4', l1=[10.0, 0.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 10.0, 0.0], center=False, makepseudo=False, settings=dict)
    # XY Plane, shifted in Z, red
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.0, 0.0]}
    plane.make_plane_points(name='p3', l1=[10.0, 0.0, 10.0], l2=[0.0, 0.0, 10.0], l3=[0.0, 10.0, 10.0], center=False, makepseudo=False, settings=dict)
    
    zoom all
    

Or make a color cube, by initial fixed plane 
    
    
    import plane
    
    python
    from chempy import cpv
    p1 = [23.76, -47.69, 45.23]
    p2 = [34.96, -18.57, -1.25]
    p3 = [90.76, -4.31, 21.69]
    v1 = cpv.sub(p1, p2)
    v2 = cpv.sub(p3, p2)
    normal = cpv.cross_product(v1, v2)
    normal_norm = cpv.normalize(normal)
    v3 = cpv.scale(normal_norm, 40)
    p4 = cpv.add(p2, v3)
    cmd.pseudoatom("mine", color="tv_blue", pos=p4)
    python end
    
    # XY Plane, forest
    dict = {'ALPHA':0.4, 'COLOR':[0.2, 0.6, 0.2], 'INVERT':True}
    plane.make_plane_points(name='p4', l1=p1, l2=p2, l3=p3, center=False, makepseudo=False, settings=dict)
    # XY Plane, shifted in Z, red
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.0, 0.0]}
    plane.make_plane_points(name='p3', l1=cpv.add(p1, v3), l2=cpv.add(p2, v3), l3=cpv.add(p3, v3), center=False, makepseudo=False, settings=dict)
    
    # XZ Plane, blue
    dict = {'ALPHA':0.4, 'COLOR':[0.0, 0.0, 1.0]}
    plane.make_plane_points(name='p2', l1=p1, l2=p2, l3=cpv.add(cpv.sub(p3, v2), v3), center=False, makepseudo=False, settings=dict)
    # XZ Plane, shifted in Y, #orange
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.5, 0.0], 'INVERT':True}
    plane.make_plane_points(name='p5', l1=cpv.add(p1, v2), l2=cpv.add(p2, v2), l3=cpv.add(p3, v3), center=False, makepseudo=False, settings=dict)
    
    # YZ Plane, #purple
    dict = {'ALPHA':0.4, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':True}
    plane.make_plane_points(name='p1', l1=cpv.add(cpv.sub(p1, v1), v3), l2=p2, l3=p3, center=False, makepseudo=False, settings=dict)
    ## YZ Plane, shifted in X, #yellow
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 1.0, 0.0]}
    plane.make_plane_points(name='p6', l1=cpv.add(p1, v3), l2=cpv.add(p2, v1), l3=cpv.add(p3, v1), center=False, makepseudo=False, settings=dict)
    
    zoom all
    

Retrieved from "[https://pymolwiki.org/index.php?title=Plane_Wizard&oldid=12500](https://pymolwiki.org/index.php?title=Plane_Wizard&oldid=12500)"


---

## Plot noe

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/plot_noe.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/plot_noe.py)  
Author(s)  | [Justin Lorieau](http://www.lorieau.com/) and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
plot_noe reads a XPLOR NOE distance restraints file and shows them as distance objects in PyMOL. 

This is a modified version of [Justin's original code](http://lorieau.com/software/biophysics-software/40-plot-xplor-noes-in-pymol.html). You may visit his site for more description and examples. 

## Usage
    
    
    plot_noe filename [, selection [, line_color [, line_width ]]]
    

## See Also

  * [distance](/index.php/Distance "Distance")



Retrieved from "[https://pymolwiki.org/index.php?title=Plot_noe&oldid=13917](https://pymolwiki.org/index.php?title=Plot_noe&oldid=13917)"


---

## Pml2py

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script converts a pml script to a python script. 

See [pymol-users mailing list](http://sourceforge.net/mailarchive/message.php?msg_id=27331786) (Subject: Convert pml script to Pymol Python script, Fri, 8 Apr 2011). 
    
    
    from __future__ import print_function
    
    import sys
    from pymol import cmd, parsing
    
    def pml2py(filename, out=sys.stdout):
        '''
    DESCRIPTION
    
        Convert a pml script to python syntax.
    
    USAGE
    
        pml2py infile [, outfile]
        '''
        def quote(args):
            args = iter(args)
            for arg in args:
                if '=' not in arg:
                    prefix = ''
                else:
                    prefix, arg = arg.split('=', 1)
                    prefix += '='
                    arg = arg.lstrip()
                yield prefix + repr(arg)
    
        class stackiter:
            def __init__(self, collection):
                self.iterator = iter(collection)
                self.stack = []
            def __iter__(self):
                return self
            def next(self):
                if len(self.stack):
                    return self.stack.pop()
                return next(self.iterator)
            __next__ = next
            def push(self, v):
                self.stack.append(v)
    
        basestring = (bytes, unicode) if (bytes is str) else (bytes, str)
        if isinstance(out, basestring):
            out = open(out, 'w')
    
        print('''
    # automatically converted from "%s"
    import pymol
    from pymol import *
    ''' % (filename), file=out)
    
        handle = stackiter(open(filename, 'rU'))
        for line in handle:
            while line.endswith('\\\n'):
                line = line[:-2] + next(handle)
    
            a = line.split(None, 1)
            if len(a) > 1 and a[0] == '_':
                line = a[1]
                a = line.split(None, 1)
    
            try:
                name = a[0]
                if name.startswith('/'):
                    line = line.lstrip()[1:]
                    raise
                name = cmd.kwhash.shortcut.get(name, name)
                kw = cmd.keyword[name]
                assert kw[4] != parsing.PYTHON
                func = kw[0]
            except:
                out.write(line)
                continue
            
            # PyMOL stuff without named python function
            if func.__name__ == '<lambda>' or name.startswith('@'):
                print('cmd.do(%s)' % (repr(line.strip())), file=out)
                continue
    
            # code blocks
            if name == 'python':
                for line in handle:
                    if line.split() == ['python', 'end']:
                        break
                    out.write(line)
                continue
            if name == 'skip':
                for line in handle:
                    if line.split() == ['skip', 'end']:
                        break
                continue
    
            # split args
            tok = ','
            if kw[4] == parsing.MOVIE:
                tok = kw[3]
                split_mx = 1
            else:
                split_mx = kw[4] - parsing.LITERAL
            if len(a) > 1:
                a[1] = parsing.split(a[1], '#', 1)[0] # strip of comment
                if split_mx < 0:
                    # strip of additional commands
                    a[1:] = parsing.split(a[1], ';', 1)
                    if len(a) > 2:
                        handle.push(a[2])
                if split_mx == 0:
                    args = [a[1]]
                else:
                    args = [i.strip() for i in parsing.split(a[1], tok, max(0, split_mx))]
            else:
                args = []
    
            # old syntax: set property=value
            if kw[4] == parsing.LEGACY and '=' in args[0]:
                args = [i.strip() for i in args[0].split('=', 1)] + args[1:]
    
            # register alias
            if name == 'alias':
                cmd.alias(*args)
    
            # use 'cmd' module if possible
            try:
                test = getattr(cmd, func.__name__)
                assert func == test
                module = 'cmd'
            except:
                module = func.__module__
    
            print('%s.%s(%s)' % (module, func.__name__, ', '.join(quote(args))), file=out)
    
    cmd.extend('pml2py', pml2py)
    
    # vi:expandtab:smarttab
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pml2py&oldid=13247](https://pymolwiki.org/index.php?title=Pml2py&oldid=13247)"


---

## PoseView

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/poseview.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/poseview.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**poseview** is a PyMOL wrapper for the [PoseView](http://www.biosolveit.de/poseview/) program. PoseView generates 2D structure-diagrams of protein-ligand complexes. 

## Setup

You need the poseview executable and a valid license file. Example for a launcher script which should be saved to **/usr/bin/poseview** so that PyMOL can find it: 
    
    
    #!/bin/sh
    POSEVIEW_PATH="/opt/poseview-1.1.2-Linux-x64"
    export BIOSOLVE_LICENSE_FILE="$POSEVIEW_PATH/poseview.lic"
    exec "$POSEVIEW_PATH/poseview" "$@"
    

## Usage
    
    
    poseview [ ligand [, protein [, width [, height [, filename [, exe [, state ]]]]]]]
    

## Example
    
    
    fetch 1a00, async=0
    poseview chain C and resn HEM
    

[![Poseviewscreenshot.png](/images/c/c7/Poseviewscreenshot.png)](/index.php/File:Poseviewscreenshot.png)

Retrieved from "[https://pymolwiki.org/index.php?title=PoseView&oldid=13918](https://pymolwiki.org/index.php?title=PoseView&oldid=13918)"


---

## PowerMate Dial OS X

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Link to W. G. Scott's web page](http://scottlab.ucsc.edu/~wgscott/xtal/wiki/index.php/Using_powermate_dials_with_PyMOL) (off site). Use of the script below is explained in that link. Briefly, the PowerMate driver must be installed and configured as explained in detail in that link before this script will work. 

Pymol lets you define some keys (the function keys, ALT-[A-Z,0-9], left, right, up, down, etc) as pythonic pymol commands, like those below. I chose 2 degree increments for rotation and 0.5 for translations because these values give a visually smooth response, but change this to suit your own taste. Type 
    
    
    help set_key 
    

in pymol for more information. 

The code below dynamically assigns functions to the four programmed keys by running python scripts in pymol to reset the key bindings. Clicking the Powermate dial to activate the toggle function (assigned to F5) allows you to toggle cyclically through dial functionality options, thus giving you three sets of rotations and translations. 

In other words, turn the dial right or left to rotate in the plus y or minus y directions. Push down while turning right or left to get the plus and minus y translations. Click the dial (press down and release rapidly) to toggle to an analogous dial set for z, and then click again for x, and then again for y. 
    
    
    from pymol import cmd
    
    # Define aliases for mapping in [x,y,z] rotations and translations into a single Powermate
    # dial.  Toggling between the three is possible if you then assign these to special keys.
    
    # Functions for x,y,z rotations and translations using Powermate dial
    # Program F1 and F2 for Rotate Right and Rotate Left
    # Program F3 and F4 for Click & Rotate Right and Click & Rotate Left
    # Program F5 for  Click  (to toggle between dialsets)
    
    # dialset = 2
    
    def dialx(): \
        global dialset \
        dialset = 1 \
        cmd.set_key ('F1', cmd.turn,('x',-2.0)) \
        cmd.set_key ('F2', cmd.turn,('x',2.0)) \
        cmd.set_key ('F3', cmd.move,('x',-0.5)) \
        cmd.set_key ('F4', cmd.move,('x',0.5)) \
        print "dialset ", dialset, " [ X ]\n" \
        return dialset
    
    def dialy(): \
        global dialset \
        dialset = 2 \
        cmd.set_key ('F1', cmd.turn,('y',-2.0)) \
        cmd.set_key ('F2', cmd.turn,('y',2.0)) \
        cmd.set_key ('F3', cmd.move,('y',-0.5)) \
        cmd.set_key ('F4', cmd.move,('y',0.5)) \
        print "dialset ", dialset, " [ Y ]\n" \
        return dialset
    
    
    def dialz(): \
        global dialset \
        dialset = 3 \
        cmd.set_key ('F1', cmd.turn,('z',-2.0)) \
        cmd.set_key ('F2', cmd.turn,('z',2.0)) \
        cmd.set_key ('F3', cmd.move,('z',-0.5)) \
        cmd.set_key ('F4', cmd.move,('z',0.5)) \
        print "dialset ", dialset, " [ Z ]\n" \
        return dialset
    
    def toggle_dial(): \
        if dialset == 1 : \
            print "Changing to y" \
            dialy() \
        elif dialset == 2 : \
            print "Changing to z" \
            dialz() \
        elif dialset == 3 : \
            print "Changing to x" \
            dialx() \
        else: print "Dial assignment isn't working"
    
    
    cmd.set_key ('F5', toggle_dial)
    
    # Start default dial state for rotate y  (arbitrary choice)
    
    dialy()
    

Retrieved from "[https://pymolwiki.org/index.php?title=PowerMate_Dial_OS_X&oldid=11416](https://pymolwiki.org/index.php?title=PowerMate_Dial_OS_X&oldid=11416)"


---

## Process All Files In Directory

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Explanation

For a given directory with PDB files in it, the following code will output, for each PDB, the bound disulfide bond lengths like this: 
    
    
    1alk.pdb
      A CYS 168 SG
      A CYS 178 SG
        1.975
      A CYS 286 SG
      A CYS 336 SG
        1.995
      B CYS 168 SG
      B CYS 178 SG
        1.996
      B CYS 286 SG
      B CYS 336 SG
        2.032
    1btu.pdb
       CYS 42 SG
       CYS 58 SG
        2.039
       CYS 136 SG
       CYS 201 SG
        2.031
       CYS 168 SG
       CYS 182 SG
        2.001
       CYS 191 SG
       CYS 220 SG
        2.019
    ...
    

# Bound Disulfides
    
    
    from pymol import cmd
    from glob import glob
    
    for file in glob("*.pdb"):
        print file
        cmd.load(file,'prot')
        for a in cmd.index("elem s and bound_to elem s"):
            if cmd.select("s1","%s`%d"%a) and \
               cmd.select("s2","elem s and bound_to %s`%d"%a):
                if cmd.select("(s1|s2) and not ?skip"):
                    cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
                    print '   ',round(cmd.dist("tmp","s1","s2"),3)
                    cmd.select("skip","s1|s2|?skip")
        cmd.delete("all")
    

# All Sulfur Distances

Note that the above is for bonded sulfurs in disulfides. For all intra-cysteine gamma sulfur distances, you'd want to do something more like: 
    
    
    1alk.pdb
      A CYS 168 SG
      A CYS 178 SG
        1.975
      A CYS 168 SG
      A CYS 286 SG
        35.845
      A CYS 168 SG
      A CYS 336 SG
        35.029
      A CYS 168 SG
      B CYS 168 SG
        63.64
      A CYS 168 SG
      B CYS 178 SG
        63.775
      A CYS 168 SG
      B CYS 286 SG
        39.02
      A CYS 168 SG
      B CYS 336 SG
        39.314
    1btu.pdb
       CYS 42 SG
       CYS 58 SG
        2.039
       CYS 42 SG
       CYS 136 SG
    
    
    
    from pymol import cmd
    from glob import glob
    
    for file in glob("*.pdb"):
        print file
        cmd.load(file,'prot')
        for a in cmd.index("CYS/SG"):
            for b in cmd.index("CYS/SG"):
                if a[1]<b[1]:
                    cmd.select("s1","%s`%d"%a)
                    cmd.select("s2","%s`%d"%b)
                    if cmd.select("(s1|s2) and not ?skip"):
                        cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
                        print '   ',round(cmd.dist("tmp","s1","s2"),3)
                        cmd.select("skip","s1|s2|?skip")
        cmd.delete("all")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Process_All_Files_In_Directory&oldid=8031](https://pymolwiki.org/index.php?title=Process_All_Files_In_Directory&oldid=8031)"


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

## Psico

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | <https://github.com/speleo3/pymol-psico>  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3") and Steffen Schmidt   
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
  
psico is a python module which extends PyMOL with many commands. It was developed as an in-house script collection, but grew to a quite comprehensive PyMOL extension. Some of the provided commands are already in the PyMOLWiki [Script Library](/index.php/Category:Script_Library "Category:Script Library"), but it contains also many new stuff. 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Initialization
  * 4 Command List



## Installation

Using PyMOL 2.0 from the [schrodinger anaconda channel](https://anaconda.org/schrodinger/): 
    
    
    conda install -c schrodinger pymol
    conda install -c speleo3 pymol-psico
    

Alternatively, copy the **psico** folder to your PYTHONPATH or use the **setup.py** based installation. 

## Usage

To import all commands with one line: 
    
    
    import psico.fullinit
    

A reference document with all commands can be generated with: 
    
    
    import psico.fullinit
    import psico.helping
    psico.helping.write_html_ref('psico-commands.html')
    

## Initialization

By importing **psico.fullinit** , the PyMOL core commands [save](/index.php/Save "Save") and [fetch](/index.php/Fetch "Fetch") become overloaded. In addition, all psico commands will be added to the **pymol.cmd** namespace. The save command will have enhanced PDB output, namely saving of secondary structure and crystal information, if available. The fetch command will also fetch 5-letter codes (code + chain) as well as CATH and SCOP identifiers (SCOP requires setting up a MySQL database). To skip the overloading, initialize psico like this: 
    
    
    import psico
    psico.init(save=0, fetch=0, pymolapi=0)
    

## Command List

This list might change. If you install psico, it's recommended to generate a reference document as described above to see which commands are actually available. 

A "outdated" note in the table means that the code in psico is newer (somehow enhanced) than on the corresponding wiki page. 

Command | Already in PyMOLWiki... | Added to PyMOL   
---|---|---  
add_missing_atoms |   
aaindex2b | [AAindex](/index.php/AAindex "AAindex") (outdated)   
alignwithanymethod | [TMalign](/index.php/TMalign "TMalign")  
alphatoall | [AlphaToAll](/index.php/AlphaToAll "AlphaToAll") (outdated)   
[angle_between_domains](/index.php/Angle_between_domains "Angle between domains") |   
angle_between_helices | [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices") (outdated)   
apbs_surface |   
api_info |  | 1.7.2 (as [api](/index.php?title=Api&action=edit&redlink=1 "Api \(page does not exist\)"))   
apropos | [apropos](/index.php/Apropos "Apropos")  
atmtypenumbers |   
biomolecule | [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat") (outdated)   
cafit_orientation | [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices")  
cbm | similar: [Color Objects](/index.php/Color_Objects "Color Objects")  
cbs |   
centerofmass | similar: [center_of_mass](/index.php/Center_of_mass "Center of mass") | [1.7.2](/index.php/Centerofmass "Centerofmass")  
closest_keyframe |   
collapse_resi | similar: [CollapseSel](/index.php/CollapseSel "CollapseSel")  
consurfdb |   
corina |   
delaunay | similar: [PyDet](/index.php/PyDet "PyDet") plugin   
diff |   
dss_promotif |   
[dssp](/index.php/Dssp "Dssp") | similar: [dssp stride](/index.php/Dssp_stride "Dssp stride") plugin   
dyndom |   
extra_fit | [extra_fit](/index.php/Extra_fit "Extra fit") | 1.7.2   
fasta | see also: [Aa_codes](/index.php/Aa_codes "Aa codes")  
fetch |   
frames2states |   
gdt_ts |   
get_keyframes |   
get_raw_distances | [get_raw_distances](/index.php/Get_raw_distances "Get raw distances")  
get_sasa |   
get_sasa_ball |   
get_sasa_mmtk |   
grepset | [grepset](/index.php/Grepset "Grepset")  
gyradius | [radius of gyration](/index.php/Radius_of_gyration "Radius of gyration") (outdated)   
helix_orientation | [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices") (outdated)   
hydropathy2b |   
intra_promix |   
[intra_theseus](/index.php/Theseus "Theseus") |   
[intra_xfit](/index.php/Intra_xfit "Intra xfit") |   
iterate_plot |   
join_states |  | [1.7.2](/index.php/Join_states "Join states")  
load_3d |   
[load_aln](/index.php/Load_aln "Load aln") |   
load_consurf |   
load_gro |   
load_msms_surface |   
load_mtz_cctbx |   
load_traj_crd |   
loadall |  | [1.7.2](/index.php/Loadall "Loadall")  
local_rms |   
loop_orientation | [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices")  
map_new_apbs |   
matrix_to_ttt |   
[mcsalign](/index.php/Mcsalign "Mcsalign") |   
morpheasy |   
morpheasy_linear |   
mse2met |  | [1.6.0](/index.php?title=Mse2met&action=edit&redlink=1 "Mse2met \(page does not exist\)")  
[msms_surface](/index.php/Msms_surface "Msms surface") |   
mutate |   
mutate_all |   
nice |   
normalmodes_mmtk |   
normalmodes_pdbmat |   
normalmodes_prody |   
paper_png |   
paper_settings |   
pca_plot |   
pdb2pqr |   
peptide_rebuild |   
peptide_rebuild_modeller |   
pir |   
plane_orientation |   
polyala |   
promix |   
prosmart |   
ramp_levels |   
remove_alt | similar: [removealt](/index.php/Removealt "Removealt")  
rms_plot |   
rmsf2b |   
save |   
save_colored_fasta |   
save_movie_mpeg1 |   
save_pdb |   
save_settings | [save_settings](/index.php/Save_settings "Save settings")  
[save_traj](/index.php/Save_traj "Save traj") | [save2traj](/index.php/Save2traj "Save2traj") (outdated)   
save_xyzr |   
sculpt_relax |   
select_distances | [get_raw_distances](/index.php/Get_raw_distances "Get raw distances")  
select_pepseq |   
select_sspick |   
set_phipsi |   
set_sequence |   
[sidechaincenters](/index.php/Sidechaincenters "Sidechaincenters") |   
snp_ncbi |   
snp_uniprot |   
spectrum_states | [spectrum_states](/index.php/Spectrum_states "Spectrum states")  
spectrumany | [spectrumany](/index.php/Spectrumany "Spectrumany") | 1.6.0 (as [spectrum](/index.php/Spectrum "Spectrum"))   
split_chains | [split_chains](/index.php/Split_chains "Split chains") | 1.6.0   
split_molecules | similar: [split_object](/index.php/Split_object "Split object")  
[sst](/index.php/Sst "Sst") |   
stride | similar: [dssp stride](/index.php/Dssp_stride "Dssp stride") plugin   
stub2ala |   
supercell | [supercell](/index.php/Supercell "Supercell")  
symdiff |   
symexpcell | [supercell](/index.php/Supercell "Supercell")  
[theseus](/index.php/Theseus "Theseus") |   
tmalign | [TMalign](/index.php/TMalign "TMalign")  
update_identifiers |   
[xfit](/index.php/Xfit "Xfit") |   
  
Retrieved from "[https://pymolwiki.org/index.php?title=Psico&oldid=13470](https://pymolwiki.org/index.php?title=Psico&oldid=13470)"


---

## Pucker

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [pucker.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/pucker.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 DESCRIPTION
  * 2 USAGE
  * 3 EXAMPLES
  * 4 PYMOL API



### DESCRIPTION

**pucker.py** is a PyMol script that returns the sugar pucker information (phase, amplitude, pucker) for a given selection. 

This script uses its own dihedral calculation scheme rather than the get_dihedral command. Thus, it is lightning fast! 

If a selection does not contain any ribose sugars then an error message is returned. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    pucker selection
    #Calculates the sugar information for the given selection
    
    pucker selection, state=1
    #Calculates the sugar information for the given selection in state 1
    
    pucker selection, state=0 
    #Iterates over all states and calculates the sugar information for the given selection
    

### EXAMPLES
    
    
    #fetch 1BNA.pdb
    fetch 1bna
    
    #Select DNA only
    #Otherwise, you will get an error for water not having sugars
    select DNA, not solvent
    
    #Execute pucker command
    pucker DNA
    
    #The output should look like this
     Phase     Amp    Pucker  Residue
    161.22   55.32  C2'-endo  A 1
    139.52   41.67   C1'-exo  A 2
     92.82   38.31  O4'-endo  A 3
    166.35   48.47  C2'-endo  A 4
    128.57   46.30   C1'-exo  A 5
    126.92   49.75   C1'-exo  A 6
    101.30   47.32  O4'-endo  A 7
    115.62   49.23   C1'-exo  A 8
    140.44   46.37   C1'-exo  A 9
    145.79   53.36  C2'-endo  A 10
    147.47   47.04  C2'-endo  A 11
    113.80   51.69   C1'-exo  A 12
     
     Phase     Amp    Pucker  Residue
    153.24   43.15  C2'-endo  B 13
    128.49   45.01   C1'-exo  B 14
     67.74   43.84   C4'-exo  B 15
    149.33   41.01  C2'-endo  B 16
    169.27   42.30  C2'-endo  B 17
    147.03   42.30  C2'-endo  B 18
    116.29   47.52   C1'-exo  B 19
    129.62   49.92   C1'-exo  B 20
    113.61   42.93   C1'-exo  B 21
    156.34   50.98  C2'-endo  B 22
    116.89   44.21   C1'-exo  B 23
     34.70   45.55  C3'-endo  B 24
    

  


### PYMOL API
    
    
    from pymol.cgo import *    # get constants
    from math import *
    from pymol import cmd
    
    def pucker(selection,state=1):
    
    	# Author: Sean Law
    	# Institute: University of Michigan
    	# E-mail: seanlaw@umich.edu
    
    	#Comparison to output from 3DNA is identical
    	#Note that the 3DNA output has the second chain
    	#printed out in reverse order
    	state=int(state)
    	if state == 0:
    		for state in range(1,cmd.count_states()+1):
    			sel=cmd.get_model(selection,state)
    			if state == 1:
    				print " " #Blank line to separate chain output
    				print "%5s  %6s  %6s  %8s  Residue" % ("State","Phase","Amp", "Pucker")
    			get_pucker(sel,iterate=state)
    	else:
    		sel=cmd.get_model(selection,state)
    		get_pucker(sel,iterate=0)
    	return
    
    def get_pucker (sel,iterate=0):
    	iterate=int(iterate)
    	first=1
    	old=" "
    	oldchain=" "
    	residue={}
    	theta={}
    	count=0
    	for atom in sel.atom:
    		new=atom.chain+" "+str(atom.resi)
    		newchain=atom.chain+" "+atom.segi
    		if oldchain != newchain and first:
    			if iterate == 0:
    				print " " #Blank line to separate chain output
    				print "%6s  %6s  %8s  Residue" % ("Phase", "Amp", "Pucker")
    		if new != old and not first:
    			#Check that all 5 atoms exist
    			if len(residue) == 15:
    				#Construct planes
    				get_dihedrals(residue,theta)
    				#Calculate pucker
    				info = pseudo(residue,theta)
    				print info+"  "+old
    			else:
    				print "There is no sugar in this residue"
    			if oldchain != newchain:
    				print " " #Blank line to separate chain output
    				print "%6s  %6s  %8s  Residue" % ("Phase", "Amp", "Pucker")
    			#Clear values
    			residue={}
    			dihedrals={}
    			theta={}
    			#Store new value
    			store_atom(atom,residue)
    		else:
    			store_atom(atom,residue)
    		first=0
    		old=new
    		oldchain=newchain
    	#Final Residue
    	#Calculate dihedrals for final residue
    	if len(residue) == 15:
    		#Construct planes
    		get_dihedrals(residue,theta)
    		#Calculate pucker for final residue
    		info = pseudo(residue,theta)
    		if iterate == 0:
    			print info+"  "+old
    		else:
    			print "%5d  %s" % (iterate,info+"  "+old)
    	else:
    		print "There is no sugar in this residue"
    	return
    
    def sele_exists(sele):
    	return sele in cmd.get_names("selections");
    
    def pseudo(residue,theta):
    	other=2*(sin(math.radians(36.0))+sin(math.radians(72.0)))
    	
    	#phase=atan2((theta4+theta1)-(theta3+theta0),2*theta2*(sin(math.radians(36.0))+sin(math.radians(72.0))))
    	phase=atan2((theta['4']+theta['1'])-(theta['3']+theta['0']),theta['2']*other)
    	amplitude=theta['2']/cos(phase)
    	phase=math.degrees(phase)
    	if phase < 0:
    		phase+=360 # 0 <= Phase < 360
    	#Determine pucker
    	if phase < 36:
    		pucker = "C3'-endo"
    	elif phase < 72:
    		pucker = "C4'-exo"
    	elif phase < 108:
    		pucker = "O4'-endo"
    	elif phase < 144:
    		pucker = "C1'-exo"
    	elif phase < 180:
    		pucker = "C2'-endo"
    	elif phase < 216:
    		pucker = "C3'-exo"
    	elif phase < 252:
    		pucker = "C4'-endo"
    	elif phase < 288:
    		pucker = "O4'-exo"
    	elif phase < 324:
    		pucker = "C1'-endo"
    	elif phase < 360:
    		pucker = "C2'-exo"
    	else:
    		pucker = "Phase is out of range"
    	info = "%6.2f  %6.2f  %8s" % (phase, amplitude, pucker)
    	return info
    	
    
    def store_atom(atom,residue):
    	if atom.name == "O4'" or atom.name == "O4*":
    		residue['O4*X'] = atom.coord[0]
    		residue['O4*Y'] = atom.coord[1]
    		residue['O4*Z'] = atom.coord[2]
    	elif atom.name == "C1'" or atom.name == "C1*":
    		residue['C1*X'] = atom.coord[0]
    		residue['C1*Y'] = atom.coord[1]
    		residue['C1*Z'] = atom.coord[2]
    	elif atom.name == "C2'" or atom.name == "C2*":
    		residue['C2*X'] = atom.coord[0]
    		residue['C2*Y'] = atom.coord[1]
    		residue['C2*Z'] = atom.coord[2]
    	elif atom.name == "C3'" or atom.name == "C3*":
    		residue['C3*X'] = atom.coord[0]
    		residue['C3*Y'] = atom.coord[1]
    		residue['C3*Z'] = atom.coord[2]
    	elif atom.name == "C4'" or atom.name == "C4*":
    		residue['C4*X'] = atom.coord[0]
    		residue['C4*Y'] = atom.coord[1]
    		residue['C4*Z'] = atom.coord[2]
    	return
    
    def get_dihedrals(residue,theta):
    
    	C = []
    	ribose = residue.keys()
    	ribose.sort()
    	
    	shift_up(ribose,6)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['0']=calcdihedral(C)
    	
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['1']=calcdihedral(C)
    
    
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['2']=calcdihedral(C)
    
    
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['3']=calcdihedral(C)
    
    	C = []
    	shift_down(ribose,3)
    	for i in range (0,12):
    		C.append(residue[ribose[i]])
    	theta['4']=calcdihedral(C)
    	
    	return
    
    def shift_up(list,value):
    	for i in range (0,value):
    		list.insert(0,list.pop())
    	return
    
    def shift_down(list,value):
    	for i in range (0,value):
    		list.append(list.pop(0))
    	return
    
    def calcdihedral(C):
    	
    	DX12=C[0]-C[3]
    	DY12=C[1]-C[4]
    	DZ12=C[2]-C[5]
    	
    	DX23=C[3]-C[6]
    	DY23=C[4]-C[7]
    	DZ23=C[5]-C[8]
    	
    	DX43=C[9]-C[6];
    	DY43=C[10]-C[7];
    	DZ43=C[11]-C[8];
    
    	#Cross product to get normal
    	
    	PX1=DY12*DZ23-DY23*DZ12;
    	PY1=DZ12*DX23-DZ23*DX12;
    	PZ1=DX12*DY23-DX23*DY12;
    
    	NP1=sqrt(PX1*PX1+PY1*PY1+PZ1*PZ1);
    
    	PX1=PX1/NP1
    	PY1=PY1/NP1
    	PZ1=PZ1/NP1
    
    	PX2=DY43*DZ23-DY23*DZ43;
    	PY2=DZ43*DX23-DZ23*DX43;
    	PZ2=DX43*DY23-DX23*DY43;
    
    	NP2=sqrt(PX2*PX2+PY2*PY2+PZ2*PZ2);
    	
    	PX2=PX2/NP2
    	PY2=PY2/NP2
    	PZ2=PZ2/NP2
    
    	DP12=PX1*PX2+PY1*PY2+PZ1*PZ2
    
    	TS=1.0-DP12*DP12
    
    	if TS < 0:
    		TS=0
    	else:
    		TS=sqrt(TS)
    	
    	ANGLE=math.pi/2.0-atan2(DP12,TS)
    
    	PX3=PY1*PZ2-PY2*PZ1
    	PY3=PZ1*PX2-PZ2*PX1
    	PZ3=PX1*PY2-PX2*PY1
    
    	DP233=PX3*DX23+PY3*DY23+PZ3*DZ23
    
    	if DP233 > 0:
    		ANGLE=-ANGLE
    
    	ANGLE=math.degrees(ANGLE)
    
    	if ANGLE > 180:
    		ANGLE-=360
    	if ANGLE < -180:
    		ANGLE+=360
    
    	return ANGLE
    
    cmd.extend("pucker",pucker)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pucker&oldid=11145](https://pymolwiki.org/index.php?title=Pucker&oldid=11145)"


---

## Pymol2glmol

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/pymol2glmol.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/pymol2glmol.py)  
Author(s)  | Takanori Nakane   
License  | [LGPL3](http://opensource.org/licenses/LGPL-3.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
  * 2 Example of use
  * 3 System requirements
    * 3.1 Fix explained for Mint 12, gnome shell 3, chrome 16



## Introduction

A script to export a scene in pymol to GLmol. GLmol is a molecular viewer for Web browsers written in WebGL/Javascript. 

With **pymol2glmol** , you can publish your pymol scene to a Web page. Visitors can rotate, zoom the molecule on the page. 

Compared to export of polygon coordinates (VRML or Object3D), the published web page contain only atomic coordinates so that the file size is much smaller and visitors can even change representation. 

[Examples and script can be downloaded from my web page.](http://webglmol.sourceforge.jp/pymol_exporter/index.html)

This script uses **cmd.get_session** to extract which representations is enabled on each part of the molecule. I think this technique is useful for many purposes, for example, writing exporters, copying representations between aligned structures, etc. 

Comments and suggestions are welcome. 

Best regards,   
Takanori Nakane 

[![Pymol2glmol.png](/images/d/d8/Pymol2glmol.png)](/index.php/File:Pymol2glmol.png)

[](/index.php/File:Pymol2glmol.png "Enlarge")

## Example of use
    
    
    reinitialize
    import pymol2glmol
    
    fetch 1TQN, async=0
    preset.pretty_solv("1TQN")
    select heme, organic
    
    pymol2glmol 1TQN
    import webbrowser
    webbrowser.open("1TQN.html")
    

## System requirements

GLmol runs on newer versions of Firefox, Chrome, Safari or Opera.   
[See here how to activate in windows.](http://www.windows7hacker.com/index.php/2010/02/how-to-enable-webgl-in-firefox-chrome-and-safari/) Internet Explorer is not supported because IE doesn't implement WebGL. 

GLmol also runs on Sony Ericsson's Android devices which support WebGL.   
Support for Firefox Mobile is currently underway.   
Reportedly, GLmol also runs on WebGL enabled safari in iOS. 

If you see only black screen and you are using 

  * Internet Explorer: sorry. IE doesn't support WebGL.
  * Firefox (version 4 or later): [try force enable WebGL.](https://wiki.mozilla.org/Blocklisting/Blocked_Graphics_Drivers#How_to_force-enable_blocked_graphics_features)
  * Chrome: [try force enable WebGL.](http://www.google.com/support/forum/p/Chrome/thread?tid=4b9244822aa2f2e0&hl=en)
  * Safari: [enable WebGL.](https://discussions.apple.com/thread/3300585?start=0&tstart=0)



### Fix explained for Mint 12, gnome shell 3, chrome 16

Install menu editor, and run it 
    
    
    sudo apt-get install alacarte
    sudo alacarte
    

Locate internet menu. Edit the "Google Chrome" entry from->to 
    
    
    /opt/google/chrome/google-chrome %U
    /opt/google/chrome/google-chrome --enable-webgl %U
    

Then find out what handles your default opening of **text/html** , and edit that .desktop file.  
If you use the repository version of chrome, then search for **chromium-browser.desktop** instead. 
    
    
    less ~/.local/share/applications/mimeapps.list | grep 'text/html'
    sudo find / -name 'google-chrome.desktop' | xargs sudo gedit
    

Then locate the line **Exec=** and insert in line **\--enable-webgl** : 
    
    
    Exec=/opt/google/chrome/google-chrome --enable-webgl %U
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pymol2glmol&oldid=13920](https://pymolwiki.org/index.php?title=Pymol2glmol&oldid=13920)"


---

## Pytms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/pytms.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/pytms.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [GPL-2.0](http://opensource.org/licenses/GPL-2.0)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![PyTMs example: PTMs in red](/images/d/dd/Pytms.gif)](/index.php/File:Pytms.gif "PyTMs example: PTMs in red")

**PyTMs** is a PyMOL plugin that enables to easily introduce a set of common post-translational modifications (PTMs) into protein models.  
Latest version: **1.2 (November 2015)**  


## Contents

  * 1 Installation
  * 2 Using the GUI
  * 3 PyMOL API functions
  * 4 Selections and custom property selectors
  * 5 PyTMs in python scripts
  * 6 Basic residue-based optimization, animation, vdW clashes, surface selections
  * 7 Hints
  * 8 Update notes
  * 9 References
  * 10 SEE ALSO



## Installation

  * For instruction on setting up plugins see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



PyTMS can be [downloaded separately](https://github.com/andwar/Pymol-script-repo/blob/master/plugins/pytms.py), or together with the [pymol-script-repository](https://github.com/Pymol-Scripts/Pymol-script-repo)  
The **pytms.py** file should be placed in the **plugins** folder. 

## Using the GUI

In PyMOL, navigate to: _P_ lugin --> PyTMs  
[![PyTMs GUI](/images/6/6a/Pytms_menu.png)](/index.php/File:Pytms_menu.png "PyTMs GUI")  
The supported PTM can be selected from the left panel. The options will appear on the right. For help, click the help button.  
Note that available objects/selection can be selected using the buttons or be entered into the field. Target atoms are automatically sub-selected and do not need to be specified. 

  


## PyMOL API functions

List of functions   
---  
PTM  | PyTMs API command  | p.PTM property selector *  | Target amino acid(s)  | Modifyable N-terminus? **   
acetylation  | acetylate  | acetylation  | Lysine  | Yes, biologically relvant   
carbamylation  | carbamylate  | carbamylation  | Lysine  | Yes, implemented   
citrullination  | citrullinate  | citrullination  | Arginine  | No   
cysteine oxidations  | oxidize_cys  | oxidation_cys  | Cysteine /(Seleno-)  | No   
malondialdehyde adducts  | mda_modify  | MDA  | Lysine  | Yes, biologically relevant   
methionine oxidation  | oxidize_met  | oxidation_met  | Methionine /(Seleno-)  | No   
methylation  | methylate  | methylation  | Lysine  | Yes, implemented   
nitration  | nitrate  | nitration  | Tyrosine (Tryptophan)  | No   
nitrosylation  | nitrosylate  | nitrosylation  | Cysteine  | No   
proline hydroxylation  | hydroxy_pro  | hydroxylation_pro  | Proline  | No   
phosphorylation  | phosphorylate  | phosphorylation  | Serine, Threonine, Tyrosine  | No   
  
  * (*) requires incentive PyMOL version 1.7, or higher
  * (**) N-terminal modifiactions excluding Proline
  * Note that each function has an **individual help description** , callable by entering, e.g.:


    
    
    help citrullinate
    

  


## Selections and custom property selectors

PTM amino acids have altered names, and altered/added atoms are identifyable by unique names.  
The nomenclature is adapted mostly from [RCSB](http://www.rcsb.org/pdb/home/home.do). 
    
    
    Examples
    select resn CIR # citrullines
    select resn NIY # nitro-tyrosines
    # etc. ...
    

  * Note that more information can be found in the associated help of individual functions
  * Hint: using the [label](/index.php/Label "Label") command is a convenient way to get hold of the names of atoms or residues



  
Incentive PyMOL version (>1.7) support custom property selectors.  
PyTMs supports these by intoducting the property **p.PTM** selector (cf. table above), which will select the altered PTM atoms:  

    
    
    # To select a PTM, use e.g.:
    select p.PTM in nitration #selects nitro groups
    select p.PTM in * #selects any PTM
    

  


## PyTMs in python scripts

In simple API scripts, the API commands can be used directly.  
For use in python scripts PyTMs can also be imported, e.g.: 
    
    
    import pmg_tk.startup.pytms
    pmg_tk.startup.pytms.mda_modify()
    

  


## Basic residue-based optimization, animation, vdW clashes, surface selections

Some functions (MDA & phosphorylation) support basic structure optimization that is used to position PTMs in more favorable locations. Note that this optimization is based solely on sterical vdW strain and currently ignores charge. Due to interation of testing, there may be a significant calculation time associated with this procedure. The optimization can be followed in the console and graphically in the PyMOL window (only in case the GUI is used). There is an option to animate the states from original to final after calculation for reference (cf. example).  
Note that the presence of hydrogens will significantly affect both the calculation during optimization and the display of vdW clashes!  
  
All functions support the display of vdW clashes after modification to allow assessment of potential steric interactions/displacements.  
This is essentially an integrated version of the original [ show_bumps script](/index.php/Show_bumps "Show bumps") by Thomas Holder.  
Note that a dedicated **display vdW strain** function is available to perform this in retrospect, or on unmodified models.   
PyTMs also enables the user to sub-select surface accessible residues. This corresponds to an integrated version of [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") by Jason Vertrees. This option is enables by providing the required cutoff (in A^2).  | **Animation of optimization with steric vdW clashes:**  
[![pytms example: optimization, animation, clashes](/images/0/08/Pytms_optimization_animation.gif)](/index.php/File:Pytms_optimization_animation.gif "pytms example: optimization, animation, clashes")  
---|---  
  
  


  


## Hints

  * console output appears first in the PyMOL console prior to the API interface
  * when using the GUI menu, the modification process can be followed in the display window



  


## Update notes
    
    
    '''
        0.90 pre release testing
        1.00 first release
            * minor changes and typo correction
            * fixed non-incentive users receiving
              error messages related to the 'alter' command
            * new feature: integrated surface selection
            * new feature: integrated display of vdW clashes for all PTMs
            * new function: display of vdW clashes indpendent of modification
        1.1 Minor fixes
            * changed boolean processing of 'delocalized' keyword
        1.2 Additional function and minor improvements
            * fixed crashed related to random selections
            * new function: nitrosylate (Cysteine S-Nitrosylation)
            * updated usage descriptions
            * the user interface window is now resizable
            * changes in nitrate-function:
                * the torsion angle can now be defined for the nitro group and has
                  a default of ~22.352 (slight angle);
                  this feature is currently only supported for Tyrosines!
                * selecting a surface cutoff in conjuction with random placement
                  of CE1 or CE2 for Nitro-tyrosines will now select the most
                  accessible CE atom rather than a random one
            * font size can now be adjusted from the Main menu
    '''
    

## References

The original citation can be found on [PubMED](http://www.ncbi.nlm.nih.gov/pubmed/25431162) A current (non-repository) version of PyTMs can be downloaded here: [pytms.py](https://github.com/andwar/Pymol-script-repo/blob/master/plugins/pytms.py)  
In case of suggestions, questions or feedback please feel free to [ contact me](/index.php/User:Andwar "User:Andwar"). 

  


## SEE ALSO

  * [Git intro](/index.php/Git_intro "Git intro")
  * [Plugin_manager](/index.php/Plugin_manager "Plugin manager")
  * [Plugins](/index.php/Plugins "Plugins")



Retrieved from "[https://pymolwiki.org/index.php?title=Pytms&oldid=12231](https://pymolwiki.org/index.php?title=Pytms&oldid=12231)"


---

## Quick dist

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**quick_dist** will quickly find and write out the pairwise distances between two selections in PyMOL. On my machine this will write about 80,000 distance measures/second. 

  

    
    
    #
    # quick_dist.py
    #
    def quick_dist(s1, s2, inFile=None):
      import math
      m1 = cmd.get_model(s1)
      m2 = cmd.get_model(s2)
    
      if inFile!=None:
        f = open(inFile, 'w')
        f.write("ATOM1\tATOM2\tDIST\n")
    
      s=""
    
      for c1 in range(len(m1.atom)):
        for c2 in range(len(m2.atom)):
          s+="%s/%d\t%s/%d\t%3.2f\n" % (m1.atom[c1].name, m1.atom[c1].index,\
            m2.atom[c2].name, m2.atom[c2].index,\
            math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(a1.coord,a2.coord)))))
    
      if inFile!=None:
        f.write(s)
        f.close()
      else:
        print s
    
    
    cmd.extend("quick_dist", quick_dist)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Quick_dist&oldid=8896](https://pymolwiki.org/index.php?title=Quick_dist&oldid=8896)"


---

## Quickdisplays

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/quickdisplays.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/quickdisplays.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The module **quickdisplays** is a package containing several functions for quick standard displays:  


List of functions   
---  
Function  | What it does  | Usage   
disp_list  | prints a list of functions  | disp_list   
disp_ss  | secondary structure  | disp_ss [ selection [, colors [, only ]]]   
disp_ball_stick  | balls and sticks  | disp_ball_stick [ selection [, hydrogens [, only ]]]   
disp_mesh  | surface mesh  | disp_mesh [ selection [, color_m [, hydrogens [, only [, limits ]]]]]   
disp_surf  | surface  | disp_surf [ selection [, color_s [, transparency [, hydrogens [, solvent [, ramp_above [, only [, limits ]]]]]]]]   
disp_putty  | putty b-factor sausage  | disp_putty [ selection [, limits [, only ]]]   
  
  * Note that each function has an **individual help description** , callable by entering, e.g.:


    
    
    help disp_surf
    

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



  


Examples   
---  
[![Quickdisplays disp surf.png](/images/3/3d/Quickdisplays_disp_surf.png)](/index.php/File:Quickdisplays_disp_surf.png)**disp_surf**  
 _surface_ | [![Quickdisplays disp ss.png](/images/a/a8/Quickdisplays_disp_ss.png)](/index.php/File:Quickdisplays_disp_ss.png)**disp_ss**  
 _secondary structure_ | [![Quickdisplays disp putty.png](/images/5/51/Quickdisplays_disp_putty.png)](/index.php/File:Quickdisplays_disp_putty.png)**disp_putty**  
 _putty b-factor sausage_ |  The _combined displays_ example was produced like this: 
    
    
    import quickdisplays
    fetch 1hpv, async=0
    
    disp_putty all, limits=10
    disp_surf color_s=lightblue, transparency=0.8
    disp_ss chain B
    disp_ball_stick hetatm
    util.cbao hetatm
    set mesh_mode, 1 # make sure hetams are meshed
    set mesh_cutoff, 4.5
    disp_mesh resn 478, color_m=green
      
  
[![Quickdisplays disp mesh.png](/images/0/07/Quickdisplays_disp_mesh.png)](/index.php/File:Quickdisplays_disp_mesh.png)**disp_mesh**  
 _mesh_ | [![Quickdisplays disp ball stick.png](/images/4/4b/Quickdisplays_disp_ball_stick.png)](/index.php/File:Quickdisplays_disp_ball_stick.png)**disp_ball_stick**  
 _balls and sticks_ | [![Quickdisplays combo.png](/images/f/fc/Quickdisplays_combo.png)](/index.php/File:Quickdisplays_combo.png)_combined displays_  
see example   
  
## Some notes on the arguments

  * **selection** can be used to restrict the display to certain objects or selections, the default is 'all'
  * **only** can be _True_ or _False_ and will toggle whether the display will be added (_show_) or replace others (_show_as_)
  * **disp_mesh** and **disp_surf** support the color **putty** , which will color the object/selection by b-factor
  * **limits** defines the b-factor color range: 
    * a list entry will define absolute values **[min;max]**
    * a float value will calculate the corresponding percentiles **±value%** , _default=5_
  * setting **hydrogens** will add hydrogen to the model; if set to _=False_ , hydrogens are removed, if omitted the state of hydogens will be 'as is'
  * **colors** in **disp_ss** is flexible: 
    * set e.g. three colors, e.g. 'red green blue' for sheets, helices and loops, respectively (cf. example)
    * alternatively certain 'util.' functions can be used (cf. example)
    * setting _False_ once or for individual colors will suppress coloring



  


## Applied example

Enter the following lines one-by-one and observe how the display will be affected 
    
    
    # import function to PyMOL
    import quickdisplays
    
    # prep
    fetch 1hpv, async=0
    bg black
    orient
    disp_list # list of functions
    
    help disp_ss # check out help
    disp_ss
    disp_ss only=T
    disp_ss colors=smudge orange grey
    disp_ss colors=chartreuse False False # False suppresses coloring
    disp_ss colors=util.chainbow # You can use util.cbc, util.rainbow, ...
    
    help disp_ball_stick # check out help
    disp_ball_stick hetatm
    disp_ball_stick all, only=T
    util.cbaw
    disp_ball_stick hydrogens=1
    disp_ball_stick hydrogens=-1
    disp_stick_ball # will this work?
    
    help disp_mesh # check out help
    disp_mesh
    disp_mesh color_m=green, only=False
    disp_mesh color_m=green, only=T
    disp_ball_stick hetatm
    set mesh_skip, 2 # omits lines in mesh
    disp_mesh color_m=default, hydrogens=1, only=T
    disp_mesh color_m=putty, hydrogens=0, limits=20
    disp_mesh color_m=putty, limits=30
    disp_mesh color_m=putty, limits=[15,50]
    
    help disp_putty # check out help
    disp_putty chain A, only=False, limits=[15,50] #only default is True
    disp_putty all, limits=0
    disp_putty all, limits=10
    
    help disp_surf # check out help
    disp_surf color_s=white, solvent=1
    disp_surf color_s=white# solvent=0
    disp_surf color_s=lightblue, transparency=0.8
    

  


## SEE ALSO

[Git intro](/index.php/Git_intro "Git intro"), [Displaying Biochemical Properties](/index.php/Displaying_Biochemical_Properties "Displaying Biochemical Properties")

Retrieved from "[https://pymolwiki.org/index.php?title=Quickdisplays&oldid=13947](https://pymolwiki.org/index.php?title=Quickdisplays&oldid=13947)"


---

## Radius of gyration

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.querying#gyradius](http://pymol.org/psicoredirect.php?psico.querying#gyradius)  
  
This script does calcualate the [radius of gyration](http://en.wikipedia.org/wiki/Radius_of_gyration) of a molecule. 

Thanks to Tsjerk Wassenaar for posting this script on the [pymol-users mailing list](http://sourceforge.net/mailarchive/message.php?msg_id=27288491)! 

_The function was added to[psico](/index.php/Psico "Psico") as "gyradius"._
    
    
    from pymol import cmd
    import math
    
    def rgyrate(selection='(all)', quiet=1):
        '''
    DESCRIPTION
    
        Radius of gyration
    
    USAGE
    
        rgyrate [ selection ]
        '''
        try:
            from itertools import izip
        except ImportError:
            izip = zip
        quiet = int(quiet)
        model = cmd.get_model(selection).atom
        x = [i.coord for i in model]
        mass = [i.get_mass() for i in model]
        xm = [(m*i,m*j,m*k) for (i,j,k),m in izip(x,mass)]
        tmass = sum(mass)
        rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in izip(x,xm))
        mm = sum((sum(i)/tmass)**2 for i in izip(*xm))
        rg = math.sqrt(rr/tmass - mm)
        if not quiet:
            print("Radius of gyration: %.2f" % (rg))
        return rg
    
    cmd.extend("rgyrate", rgyrate)
    
    # vi:expandtab
    

### Example

[![Radius of gyration Example.png](/images/2/2b/Radius_of_gyration_Example.png)](/index.php/File:Radius_of_gyration_Example.png)

[](/index.php/File:Radius_of_gyration_Example.png "Enlarge")

This example uses the radius of gyration and the [center of mass](/index.php/Center_Of_Mass "Center Of Mass") to display a semitransparent sphere. 
    
    
    from pymol import cmd, stored, util
    import centerOfMass, radiusOfGyration
    
    cmd.set('sphere_transparency', 0.5)
    
    cmd.fetch('2xwu')
    cmd.hide('everything')
    cmd.show('cartoon', 'chain A')
    
    r = radiusOfGyration.rgyrate('chain A and polymer')
    centerOfMass.com('chain A and polymer', object='com', vdw=r)
    
    util.cbc()
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass")
  * [Center Of Mass](/index.php/Center_Of_Mass "Center Of Mass")
  * [COM](/index.php/COM "COM")



Retrieved from "[https://pymolwiki.org/index.php?title=Radius_of_gyration&oldid=13172](https://pymolwiki.org/index.php?title=Radius_of_gyration&oldid=13172)"


---

## Rasmolify

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Here it is! Long awaited, less tested; 

## Contents

  * 1 Install
  * 2 Usage
  * 3 Related
  * 4 TODO
  * 5 Code



## Install

Linux
    In your ~/.pymolrc set something like the following 
    
    
     run ~/pymolscripts/rasmolify.py 

Finally, make a directory called ~/pymolscripts and copy the code below into a file called rasmolify.py - That should do the trick. You may also like to add a line that reads 
    
    
    set virtual_trackball, off

in your ~/.pymolrc

Windows
    ???

## Usage

Think 'rasmol' 

## Related

  * <http://arcib.dowling.edu/sbevsl/>



## TODO

  * Check if a 'selection' exists, and limit commands to that selection (map the concept of a rasmol 'selection' onto the concept of a pymol selection).
  * Implement 'scaling' units for display functions
  * Fix the mouse behaviour?
  * Add a rasmol GUI!



## Code
    
    
    ## This is just a quick hack. For something more meaty see;
    ## http://arcib.dowling.edu/sbevsl/
     
    ## Version 0.0.00-000/1
     
     
    ## Turn off the virtual_trackball
    cmd.set("virtual_trackball", "off")
     
     
    ## spacefill
    def spacefill(p1=''):
        if(p1=='off'):
            cmd.hide("spheres")
        elif(p1==''):
            cmd.show("spheres")
        else:
            print("feh!")
    cmd.extend("spacefill", spacefill)
     
    ## cartoon
    def cartoon(p1=''):
        if(p1=='off'):
            cmd.hide("cartoon")
        elif(p1==''):
            cmd.show("cartoon")
        else:
            print("feh!")
    cmd.extend("cartoon", cartoon)
     
    ## wireframe
    def wireframe(p1=''):
        if(p1=='off'):
            cmd.hide("lines")
        elif(p1==''):
            cmd.show("lines")
        else:
            print("feh!")
    cmd.extend("wireframe", wireframe)
     
     
    ## exit
    def exit():
        cmd.quit()
    cmd.extend("exit", exit)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Rasmolify&oldid=6375](https://pymolwiki.org/index.php?title=Rasmolify&oldid=6375)"


---

## Read PDB-String

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    delete all
    cmd.read_pdbstr("""HEADER    CREATED BY CONVERTPROSPECT              27-JAN-02   2tnf  \
    REMARK   1                                                                       \
    ATOM      1  N   PRO A   9       1.895  67.213 -38.182  1.00  0.00           N   \
    ATOM      2  CA  PRO A   9       1.703  68.680 -38.402  1.00  0.00           C   \
    ....
    ATOM   1153  C   GLY A 157       6.927  59.108 -38.901  1.00  6.00           C   \
    ATOM   1154  O   GLY A 157       6.700  59.292 -37.676  1.00  6.00           O   \
    TER    1155      GLY A 157                                                       \
    MASTER                                                                           \
    END                                                                              \
    ""","2tnfa")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Read_PDB-String&oldid=6362](https://pymolwiki.org/index.php?title=Read_PDB-String&oldid=6362)"


---

## Removealt

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/removealt.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/removealt.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | Free   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[removeAlt](/index.php/RemoveAlt "RemoveAlt") removes all atoms from **obj** that have alternate locations but aren't altloc **keep**. 

## Usage
    
    
    removealt [ obj [, keep ]]
    

## Example
    
    
    import removealt
    
    fetch 1hxb, async=0
    
    select ligA, alt a
    select ligB, alt b
    
    count_atoms ligA    # 49 atoms
    count_atoms ligB    # 49 atoms
    
    removealt 1hxb, b
    
    count_atoms ligA    # 0 atoms
    count_atoms ligB    # 49 atoms
    

## See Also

[Handling Alternate Locations](/index.php/Property_Selectors#Alternate_Locations "Property Selectors") and [alter](/index.php/Alter "Alter")

Retrieved from "[https://pymolwiki.org/index.php?title=Removealt&oldid=13921](https://pymolwiki.org/index.php?title=Removealt&oldid=13921)"


---

## Renumber

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/renumber.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/renumber.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
renumber sets new residue numbers (resi) for a polymer based on connectivity. 

## Example

This examples takes a pdb structure with insertion codes and sets a new, linear numbering based on [Q8N2U3_HUMAN](http://www.uniprot.org/uniprot/Q8N2U3_HUMAN). 
    
    
    fetch 1h4w, async=0
    
    # move everything which is not polymer to another chain
    alter not polymer, chain="B"
    
    # renumber polymer, first 27 residues of Q8N2U3_HUMAN missing.
    renumber chain A, 28
    

This example fixes numbering after concatenating two chains with [fuse](/index.php/Fuse "Fuse"). Note that the cartoon representation and the [sequence viewer](/index.php/Seq_view "Seq view") need [sorting](/index.php/Sort "Sort") to display correctly. 
    
    
    fab ACDEFG, chain1
    fab HIKLMN, chain2
    
    disable chain1
    as cartoon
    
    fuse last (chain1 and name C), first (chain2 and name N)
    renumber chain2
    sort chain2
    

## Troubleshooting

For large molecules it might be necessary to increase the [python recursion limit](http://docs.python.org/library/sys.html#sys.getrecursionlimit): 
    
    
    import sys
    sys.setrecursionlimit(10**5)
    

## See Also

  * [alter](/index.php/Alter "Alter")
  * [rename](/index.php/Rename "Rename")
  * [pdb_retain_ids](/index.php/Pdb_retain_ids "Pdb retain ids")



Retrieved from "[https://pymolwiki.org/index.php?title=Renumber&oldid=13922](https://pymolwiki.org/index.php?title=Renumber&oldid=13922)"


---

## Resicolor

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Description

Here is a small script to color proteins according to residue type. Call it from your .pymolrc or from within pymol: "run resicolor.py". Then, issuing "resicolor" will apply the scheme. Selections are supported: "resicolor chain A" will apply the scheme only to chain A. This functionality is also available as a plugin ([Resicolor_plugin](/index.php/Resicolor_plugin "Resicolor plugin")). 
    
    
    #script contributed by Philippe Garteiser; garteiserp@omrf.ouhsc.edu
    from pymol import cmd
    
    def resicolor(selection='all'):
        
        '''USAGE: resicolor <selection>
        colors all or the given selection with arbitrary
        coloring scheme.
        '''
        cmd.select ('calcium','resn ca or resn cal')
        cmd.select ('acid','resn asp or resn glu or resn cgu')
        cmd.select ('basic','resn arg or resn lys or resn his')
        cmd.select ('nonpolar','resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
        cmd.select ('polar','resn ser or resn thr or resn asn or resn gln or resn tyr')
        cmd.select ('cys','resn cys or resn cyx')
        cmd.select ('backbone','name ca or name n or name c or name o')
        cmd.select ('none')
    
        print selection
        code={'acid'    :  'red'    ,
              'basic'   :  'blue'   ,
              'nonpolar':  'orange' ,
              'polar'   :  'green'  ,
              'cys'     :  'yellow'}
        cmd.select ('none')
        for elem in code:
            line='color '+code[elem]+','+elem+'&'+selection
            print line
            cmd.do (line)
        word='color white,backbone &'+selection
        print word
        cmd.do (word)                  #Used to be in code, but looks like
                                       #dictionnaries are accessed at random
        cmd.hide ('everything','resn HOH')
    
    cmd.extend ('resicolor',resicolor)
    

## Download

[Resicolor-0.1.tar.zip‎](/images/f/ff/Resicolor-0.1.tar.zip "Resicolor-0.1.tar.zip")

Retrieved from "[https://pymolwiki.org/index.php?title=Resicolor&oldid=12198](https://pymolwiki.org/index.php?title=Resicolor&oldid=12198)"


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

## Rotkit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/rotkit.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/rotkit.py)  
Author(s)  | [Troels E. Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Introduction
    * 1.1 Functions available in PyMOL
  * 2 Example of use
    * 2.1 Example 1 - Make a rotation of domain
    * 2.2 Example 2 - Simulate dye freedom
    * 2.3 Example 3 - Create distance distribution histogram
    * 2.4 Example 4 - A tutorial file



## Introduction

This script-kit is a collection of small script to be able to precisely to put a molecule (like a dye) where you want in relation to a protein. 

You can also create rotational states of a domain or simulate a dye freedom. 

It simply makes the [ PyMOL TTT matrixes](/index.php/Object_Matrix "Object Matrix"), in a easy and user friendly way. The calls to the functions available in PyMOL, takes care of all the conversion of input and such. 

If you are interested in this, you might also want to check out the PyMOL [Chempy](/index.php/Chempy "Chempy") module that is included in PyMOL. It provides handy vector and matrix functions. 

### Functions available in PyMOL

  * rotateline(Pos1,Pos2,degangle,molecule): 

    "Pos1->Pos2" define a line whereabout "molecule" will be rotated "degangle" degrees
    **rotateline Pos1=P513C_CA, Pos2=P513C_CB, degangle=5, molecule=Atto590**
    **rotateline Pos1=dyeatom87, Pos2=dyeatom85, degangle=10, molecule=Atto590**
  * mutate(molecule,chain,resi,target="CYS",mutframe="1"): 

    Mutate a /molecule//chain/resi into a target, and selecting most probable frame 1
    **mutate 1HP1, chain=A, resi=515, target=CYS, mutframe=1**
  * toline(Pos1,Pos2,atom,molecule,dist=1): 

    Translate molecule atom, 1 angstrom away in the same direction Pos1->Pos2 specify
    **toline Pos1=P513C_CA, Pos2=P513C_CB, atom=dyeatom87, molecule=Atto590, dist=3**



**Available through rotkit.functionname**

  * printMat(matrix): 

    prints the TTT matrix in a readable format. (4X4)
  * getxyz(Sel): 

    output is a list [x,y,z] in float. The input can be a list, a string(list) or a selection.
  * vector(Sel1,Sel2): 

    Finds the vector between points. Gets the xyz list from getxyz, so input can be anything.
  * vectorstr(vector): 

    turn a vector in list format into string. No real function actually.
  * transmat(vector,dist=1): 

    Makes a TTT translation matrix for according to the input vector. The vector is multiplied with dist.
  * unitvector(vector): 

    Make a vector a unitvector.
  * radangle(angle): 

    Convert degree to radians. Not that all input are assumed to be in degrees, and are converted automatically.
  * rotmat(angle,vectornorm,pointcoord): 

    This function is the most important. That makes the TTT matrix that rotates a molecule around a normalized vector, which goes through a coordinate point.
  * crossprod(Vector1, Vector2): 

    Makes a crossproduct between two vectors
  * crosspoint(Pos1, crossprod): 

    Returns the endpoint for the Position plus the crossproduct vector. Suitable if one would like to rotate around a crossvector.


  * [![Place your dye how you want. You always get same ending position.](/images/1/16/Rotkitex1.png)](/index.php/File:Rotkitex1.png "Place your dye how you want. You always get same ending position.")

Place your dye how you want. You always get same ending position. 

  * [![Mutate a residue through command line, and quickly put your dye there.](/images/d/de/Rotkitex2.png)](/index.php/File:Rotkitex2.png "Mutate a residue through command line, and quickly put your dye there.")

Mutate a residue through command line, and quickly put your dye there. 




## Example of use

### Example 1 - Make a rotation of domain

Download: [examples/rotkit_1.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_1.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_1.pml>" highlight="python" />

### Example 2 - Simulate dye freedom

Download: [examples/rotkit_2.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_2.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_2.pml>" highlight="python" />

### Example 3 - Create distance distribution histogram

Download: [examples/rotkit_3.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_3.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_3.pml>" highlight="python" />

### Example 4 - A tutorial file

To understand how the functions works, read through the tutorial. Hash/Unhash "##" each step at the time to see the effect. To be able to follow the tutorial, you need the dye molecule, which is loaded from the Pymol-script-repository. 

Download: [examples/rotkit_4.pml](https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_4.pml)  
---  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
<include src="<https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/examples/rotkit_4.pml>" highlight="python" />

Retrieved from "[https://pymolwiki.org/index.php?title=Rotkit&oldid=13924](https://pymolwiki.org/index.php?title=Rotkit&oldid=13924)"


---

## RUCAP UM-5

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 RUCAP UM-5 tracker
  * 2 Operating principle
  * 3 5 degrees of freedom
  * 4 Getting UM-5 data of the users position into socket by your script rucap.py
  * 5 How use MolVizGui.py and WxPython in PyMol with Cygwin
  * 6 Availability
    * 6.1 Additional Information
    * 6.2 PyMol Script of ImmersiveViz



## RUCAP UM-5 tracker

RUCAP UM-5 tracker is a wireless joystick for view control in PyMol. The code based on project [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz").  
  


The script, which is reproduced below, transmits data to [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz"). [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz") runs on top of PyMol and adjusts the PyMol screen view in accordance with the position of your eyes in the space in front of the monitor.  
  


Our system contains a head tracking thread which communicates the users position to a PyMol script via a socket. The PyMol script unpacks the message and updates the world respectively. All rotation is done around the virtual origin in PyMol and zoom is also considered.  
  


## Operating principle

Tracker RUCAP UM-5 in real time detects the absolute position and orientation of Your head relative to computer monitor. Antenna RU constantly traces the position of emitter CAP and transmits the data to RUCAP UM-5 Service program rucap.dll for further processing. Information from ultrasonic sensors is processed into concrete commands for the PyMol or computer program which runs RUCAP UM-5 profile at the moment.  
  


## 5 degrees of freedom

There are 6 types of movement in three-dimensional space: linear move along X, Y, Z axes and also rotations of head around each of the axes. RUCAP UM-5 supports 5 degrees of freedom: all types of movement and head rotations left-right (yaw) and up-down (pitch). These are all types of view control that are used in computer games that's why RUCAP UM-5 tracker gives You complete view control in PyMol. 

## Getting UM-5 data of the users position into socket by your script rucap.py

This device has got only OS Windows support on current time. 

First, you need runing PyMol with [MolViz.py](https://code.google.com/p/immersive-viz/source/browse/trunk/MolViz.py). And optional with [MolVizGui.py](https://code.google.com/p/immersive-viz/source/browse/trunk/MolVizGeneratedGui.py), that make fine settings of interface RUCAP (That is [UserManual](https://code.google.com/p/immersive-viz/wiki/UserManual) for GUI). Any model is good, but for example, [1jnx.pdb](https://code.google.com/p/immersive-viz/source/browse/trunk/1jnx.pdb) you will take from the same [repo](https://code.google.com/p/immersive-viz/source/browse/trunk/). 
    
    
        ./pymol -l MolViz.py 1jnx.pdb | python MolVizGui.py
    

MolViz.py runing the [SocketServer.py](https://code.google.com/p/immersive-viz/source/browse/trunk/SocketServer.py). 

For convert data to MolViz format you need run next socket script client **rucap.py** : 
    
    
    # -*- coding: cp1251
    # rucup.py - wraper for RUCAP driver (rucap.dll - http://rucap.ru/um5/support#sdk)
    # need code from https://code.google.com/p/immersive-viz/source/browse/trunk/
    # Take data from UM-5 and put to the socket on port 4440.
    #----- The server used to listen for head tracking
    #		self.server = Server(self, 4440)
    #----- The server used to listen for the GUI to update the values / display info
    #		self.gui = Server(self, 4441)
    # 
    
    from ctypes import *
    import socket
    
    HOST = "localhost"                 # remote host-computer (localhost)
    PORT = 4440              # port on remote host
    
    RUCAP_VERSION = 3
    RUCAP_RESULT_OK = 0;
    RUCAP_RESULT_FAIL = 1;
    RUCAP_RESULT_TRACKERAPP_OFF = 2;
    RUCAP_RESULT_UNPLUGGED = 3;
    RUCAP_RESULT_DEVICE_NOT_CREATED = 4;
    RUCAP_RESULT_BAD_VERSION = 5;
    
    # load dll
    a = CDLL("rucap.dll")
    # call function, that returned int
    a.RucapCreateDevice.restype = c_int
    print "RucapCreateDevice:"
    print a.RucapCreateDevice(c_int(RUCAP_VERSION))
    
    # initialisation of device
    if (a.RucapCreateDevice(c_int(RUCAP_VERSION)) == RUCAP_RESULT_OK): print "Can`t create ruCap device\n"
    
    a.RucapGetHeadPos.restype = c_int
    x = c_float(0.0)
    y = c_float(0.0)
    z = c_float(0.0)
    
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((HOST, PORT))
    
    # main loop
    while (1) :
        if a.RucapGetHeadPos(byref(x),byref(y),byref(z)):
            print "No Position!"
        else :
            #print '%f %f %f\n' % ( x.value, y.value, z.value)
            print x.value, y.value, z.value
            sock.send('X%.15f;Y%.15f;Z%.15f;\n' % ( x.value, y.value, z.value))
            #result = sock.recv(1024)
            #print "Recive:", result
    
    sock.close()
    

## How use MolVizGui.py and WxPython in PyMol with Cygwin

Istruction about compile WxPython for Cygwin on OS Windows - [WxPythonCygwin](http://gnuradio.org/redmine/projects/gnuradio/wiki/WxPythonCygwin). 

## Availability

### Additional Information

  * **RUCAP UM-5 tracker:** <http://rucap.ru/en/um5/about>



### PyMol Script of [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz")

  * **Project Page:** <http://code.google.com/p/immersive-viz/>
  * **Source:** [here](http://code.google.com/p/immersive-viz/source/browse/trunk/MolViz.py) and [here](http://code.google.com/p/immersive-viz/source/browse)
  * **Instructions:** [here](http://code.google.com/p/immersive-viz/) and [here](http://code.google.com/p/immersive-viz/wiki/UserManual)



Retrieved from "[https://pymolwiki.org/index.php?title=RUCAP_UM-5&oldid=11279](https://pymolwiki.org/index.php?title=RUCAP_UM-5&oldid=11279)"


---

## Save Mopac

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/save_mopac.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/save_mopac.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
**save_mopac** attempts to save a PDB file in the MOPAC file format. 

## Usage
    
    
    save_mopac filename [, selection [, zero [, state ]]]
    

## See Also

  * [thread on the pymol-users mailing list](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10678.html)
  * [openmopac.net](http://openmopac.net/)



Retrieved from "[https://pymolwiki.org/index.php?title=Save_Mopac&oldid=13948](https://pymolwiki.org/index.php?title=Save_Mopac&oldid=13948)"


---

## Save pdb with anisou

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/save_pdb_with_anisou.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/save_pdb_with_anisou.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Save in PDB format including ANISOU records. 

**Warning:** PyMOL does not rotate anisotropic temperature factors if you rotate an object (for example with [align](/index.php/Align "Align") or [rotate](/index.php/Rotate "Rotate")). See also related bugs on Sourceforge. 

## Usage
    
    
    save_pdb_with_anisou filename [, selection [, state ]]
    

## See Also

  * [save](/index.php/Save "Save")
  * [ellipsoids](/index.php/Ellipsoids "Ellipsoids")
  * On Sourceforge Tracker: 
    * [#3371398](https://sourceforge.net/tracker/?func=detail&aid=3371398&group_id=4546&atid=354546)
    * [#2112541](https://sourceforge.net/tracker/?func=detail&aid=2112541&group_id=4546&atid=104546)



Retrieved from "[https://pymolwiki.org/index.php?title=Save_pdb_with_anisou&oldid=13925](https://pymolwiki.org/index.php?title=Save_pdb_with_anisou&oldid=13925)"


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

## Save settings

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/save_settings.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/save_settings.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
save_settings is a prototype implementation to save all settings with non-default values to a file. 

By default settings will be stored to **~/.pymolrc-settings.py** which is automatically recognised and loaded by PyMOL on startup. You have to manually delete the file if you no more want those settings to be loaded. 

This answers a [feature request on sourceforge](https://sourceforge.net/tracker/?func=detail&aid=1009951&group_id=4546&atid=354546). 

## Usage
    
    
    save_settings [ filename ]
    

## See Also

  * [pymolrc](/index.php/Pymolrc "Pymolrc")



Retrieved from "[https://pymolwiki.org/index.php?title=Save_settings&oldid=13926](https://pymolwiki.org/index.php?title=Save_settings&oldid=13926)"


---

## Save2traj

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Note:** See [save_traj](/index.php/Save_traj "Save traj") for an improved version of this script. 

## Contents

  * 1 Description
  * 2 USAGE
  * 3 EXAMPLES
  * 4 Script



### Description

WARNING: This script is still under development so please use at your own risk! 

This script can be used to generate a DCD trajectory file after you have loaded multiple states into a single object. 

Currently, the only supported trajectory format is DCD but there will be support for other formats in the future. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    save2traj (selection, name)
    

"name" corresponds to the output name and ".dcd" will be appended to the end of the name automatically. 

The trajectory file that is generated could be numbered differently than your PDB file since PyMOL modifies the RESID values. Thus, it may be necessary to open your original PDB file and save the molecule (with modified RESIDs) to a new file. Then, load the trajectory into that new file. Alternatively, one could retain the RESIDs as well before importing a PDB structure. 

### EXAMPLES

save2traj 1bna & c. A+B, output 

This will select chains A and B from the object 1bna and save the selection from each state into the trajectory file called "output.dcd". 

### Script

**WARNING: This script is still under development so please use at your own risk!**
    
    
    from pymol import cmd
    from pymol import stored
    import struct
    
    def save2traj (selection,name,format="dcd"):
    
        #Author: Sean Law
        #Michigan State University
    
        #Get NATOMS, NSTATES
        NATOMS=cmd.count_atoms(selection)
        NSTATES=cmd.count_states()
    
        #Determine Trajectory Format
        format=format.lower()
        if (format == "charmm"):
            name=name+".dcd"
        elif (format == "amber"):
            print "The amber format has not been implemented yet"
            return
            name=name+".trj"
        elif (format == "trj"):
            print "The amber format has not been implemented yet"
            return
            name=name+".trj"
            format="amber"
        else:
            name=name+".dcd"
            format="charmm"
    
        f=open(name,'wb')
    
        #Write Trajectory Header Information
        if (format == "charmm"):
            writeCHARMMheader(f,NSTATES,NATOMS)
        elif (format == "amber"):
            print "The amber format has not been implemented yet"
            return
        else:
            print "Unknown format"
            return
    
        #Write Trajectory Coordinates
        if (format == "charmm"):
            fmt=str(NATOMS)+'f'
            for state in range(cmd.count_states()):
                stored.xyz=[]
                cmd.iterate_state(state,selection,"stored.xyz.append([x,y,z])")
                for xyz in range (3):
                    coor=[]
                    for atom in range (NATOMS):
                        coor.append(stored.xyz[atom][xyz])
                    writeFortran(f,coor,fmt,length=NATOMS*4)
        elif (format == "amber"):
            print "The amber format has not been implemented yet"
            return
        else:
            print "Unknown format"
            return
    
        f.close()
        return
    cmd.extend("save2traj",save2traj)
    
    def writeCHARMMheader (f,NSTATES,NATOMS):
        header=['CORD']
        fmt='4s9i1f10i'
        icontrol=[]
        for i in range(20):
            #Initialize icontrol
            icontrol.insert(i,0)
        icontrol[0]=NSTATES
        icontrol[1]=1 
        icontrol[2]=1
        icontrol[3]=NSTATES
        icontrol[7]=NATOMS*3-6
        icontrol[9]=2.045473
        icontrol[19]=27
    
        for i in range(20):
            header.append(icontrol[i])
    
        writeFortran(f,header,fmt)
    
        #Title
        fmt='i80s80s'
        title=[2]
        title.append('* TITLE')
        while (len(title[1])<80):
            title[1]=title[1]+' '
        title.append('* Generated by savetraj.py (Author: Sean Law)')
        while (len(title[2])< 80):
            title[2]=title[2]+' '
        
        writeFortran(f,title,fmt,length=160+4)
    
        #NATOM
        fmt='i'
        writeFortran(f,[NATOMS],fmt)
    
        return
    cmd.extend("writeCHARMMheader",writeCHARMMheader)
    
    def readtraj (name):
        f=open(name,'rb')
        #Read Header
        fmt='4s9i1f10i'
        header=readFortran(f,fmt)
        frames=header[1]
    
        #Read Title
        readCHARMMtitle(f)
    
        #Read NATOMS
        fmt='i'    
        [NATOMS]=readFortran(f,fmt)
        
        #Read COORDINATES
        for frame in range(frames):
            fmt=str(NATOMS)+'f'
            x=readFortran(f,fmt)
            y=readFortran(f,fmt)
            z=readFortran(f,fmt)
    
        f.close()
        return
    cmd.extend("readtraj",readtraj)
    
    def writeFortran (f,buffer,fmt,length=0):
        
        if (length == 0):
            length=len(buffer)*4
        f.write(struct.pack('i',length)) #Fortran unformatted
        f.write(struct.pack(fmt,*(buffer)))
        f.write(struct.pack('i',length)) #Fortran unformatted
    
        return
    cmd.extend("writeFortran",writeFortran)
    
    def readFortran (f,fmt):
        [bytes]=struct.unpack('i',f.read(4))
        buffer=struct.unpack(fmt,f.read(bytes))
        [bytes]=struct.unpack('i',f.read(4))
    
        return buffer
    cmd.extend("readFortran",readFortran)
    
    def readCHARMMtitle(f):
        [bytes]=struct.unpack('i',f.read(4))
        [lines]=struct.unpack('i',f.read(4))
        fmt=''
        for line in range(lines):
            fmt=fmt+'80s'
        buffer=(struct.unpack(fmt,f.read(80*lines)))
        [bytes]=struct.unpack('i',f.read(4))
        #print buffer
        return
    cmd.extend("readCHARMMtitle",readCHARMMtitle)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Save2traj&oldid=12891](https://pymolwiki.org/index.php?title=Save2traj&oldid=12891)"


---

## SaveGroup

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script creates the command "saveGroup". saveGroup will save the specified group to disk as one file with all group objects in it, or as many files one pdb in each file. To save to one file just call 
    
    
    saveGroup groupName
    

to save to many files call 
    
    
    saveGroup groupName, oneFile=False
    

. 

# The Code
    
    
    import pymol
    from pymol import cmd, stored
    
    
    def saveGroup(g, oneFile=None):
        """
        Save all files inside group 'g' to either
        one file or a list of files
        
        PARAMS
        g
            name of the group to save
    
        oneFile
            if not specified or None, saves each protein in the group
            to its own file, if oneFile=True, then saves all the files
            in the group to just one file.
            
        RETURNS
            None
        """
    
        oneFile = (oneFile!=None)
    
        if cmd.get_type(g) != "object:group":
            print "Error: please provide a group name to save."
            return
    
        stored.models = set()
        cmd.iterate(g, 'stored.models.add(model)')
        
        if oneFile:
            cmd.save( g + ".pdb", g)
        else:
            for x in stored.models:
                cmd.save( x + ".pdb", x)
    
    
    cmd.extend("saveGroup", saveGroup)
    

# See Also

[group](/index.php/Group "Group") [toGroup](/index.php/ToGroup "ToGroup")

Retrieved from "[https://pymolwiki.org/index.php?title=SaveGroup&oldid=8682](https://pymolwiki.org/index.php?title=SaveGroup&oldid=8682)"


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

## Select sites

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/select_sites.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/select_sites.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3") & [Troels Linnet](/index.php/User:Tlinnet "User:Tlinnet")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
select_sites make named selections from SITE, LINK, SSBOND, HET records.  
A super fast way to annotate the protein, according to the authors input i pdb file. 

See:   
<http://www.wwpdb.org/procedure.html#toc_10>

For instance, trypsin (<http://www.rcsb.org/pdb/files/1SGT.pdb?headerOnly=YES>): 
    
    
    SITE     1 CAT  3 HIS A  57  ASP A 102  SER A 195
    SITE     1 AC1  6 ASP A 165  ALA A 177A GLU A 180  GLU A 230
    SITE     2 AC1  6 HOH A 259  HOH A 261
    

This requires that the authors have performed annotation in the pdb file. 

## Usage
    
    
    import select_sites
    select_sites [ selection [, filename [, prefix [, nice [, quiet ]]]]]
    

nice = 0 or 1: make colored sticks representation for sites {default :1} 

## Example
    
    
    fetch 1sgt, async=0
    import select_sites
    select_sites
    

Or fetch and select_sites in one go 
    
    
    import select_sites
    sites 1sgt
    

## See Also

  * [uniprot_features](/index.php/Uniprot_features "Uniprot features")



Retrieved from "[https://pymolwiki.org/index.php?title=Select_sites&oldid=13927](https://pymolwiki.org/index.php?title=Select_sites&oldid=13927)"


---

## SelectClipped

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script selects those atoms that are between the near and far [clipping plane](/index.php/Clip "Clip"). 

## Usage
    
    
    select_clipped [name [, selection [, state]]]
    

## The Code
    
    
    from pymol import cmd
    
    def clipped_by(at,v):
       x,y,z = at.coord
       nz = v[2]*(x-v[12])+v[5]*(y-v[13])+v[8]*(z-v[14])-v[11]
       return nz > v[15] and nz < v[16]
    
    def clipped(selection="all",state=1):
       v = cmd.get_view()
       return [ i.id for i in cmd.get_model(selection,state).atom if clipped_by(i,v) ]
    
    def select_clipped(name='clipped', selection='all', state=1):
        '''
    DESCRIPTION
    
        Select atoms between clipping planes
    
    USAGE
    
        select_clipped [name [, selection [, state]]]
        '''
        state = int(state)
        cmd.select(name, 'none')
        for model in cmd.get_object_list('(' + selection + ')'):
            cmd.select_list('__tmp', model, clipped('(%s) and (%s)' % (model, selection), state))
            cmd.select(name, '__tmp', merge=1)
        cmd.delete('__tmp')
    
    cmd.extend('select_clipped', select_clipped)
    

To use this script save it to "clipped.py". Then load the script from PyMOL ("run clipped.py"). 

## Example
    
    
    # for example, fetch 1cll
    
    fetch 1cll, async=0
    
    # orient it to clip out some data
    
    set_view (\
        -0.555155039,    0.671998382,    0.490123481,\
         0.050991829,    0.615659416,   -0.786360741,\
        -0.830182374,   -0.411559790,   -0.376052886,\
         0.000000000,    0.000000000, -188.930786133,\
        13.649702072,   25.646595001,   12.279129982,\
       172.813293457,  205.048278809,  -20.000000000 )
    
    # run the script
    
    select_clipped 1cllClipped, 1cll
    
    # setup the view
    
    orient
    
    # remove the surface flag on these atoms
    
    flag exfoliate, 1cllClipped
    
    # show everything as sticks
    
    show sticks
    
    # show those atoms not exfoliated as surface
    
    show surface
    

Retrieved from "[https://pymolwiki.org/index.php?title=SelectClipped&oldid=9414](https://pymolwiki.org/index.php?title=SelectClipped&oldid=9414)"


---

## Selection Exists

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    # The function below will return true if the selection is defined.
    #
    from pymol import cmd
    from types import *
    
    def sele_exists(sele):
        sess = cmd.get_session()
        for i in sess["names"]:
            if type(i) is ListType:
                if sele==i[0]:
                    return 1
        return 0
    
    # Note from Warren:
    #
    # simpler, faster:
    
    def sele_exists(sele):
       return sele in cmd.get_names("selections")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Selection_Exists&oldid=6347](https://pymolwiki.org/index.php?title=Selection_Exists&oldid=6347)"


---

## SelInside

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

selInside selects all atoms of a certain type that are located within the bounding box defined by some selection. See examples. 

# Usage
    
    
    # find all atoms inside the box spanned by (0,0,0) and (10,10,10)
    # for the example PDB, 1hpv
    fetch 1hpv
    pseudoatom b1, pos=[0,0,0]
    pseudoatom b2, pos=[10,10,10]
    selInside b1 or b2
    

  


# The Code
    
    
    import pymol
    from pymol import cmd, stored
    from random import randint
    
    def selInside(bounding_object, target="*", radius=None):
        """
        Selects all atoms of a given target type inside some selection that acts like a bounding box.
    
        PARAMS:
            bounding_object -- [PyMOL object/selection] the object/selection inside of which we want to choose atoms
            target -- [PyMOL selection] that defines the atoms to select.  Use "solvent" to select
                      only solvent atoms inside the selection
            radius -- [float] expand the box by this many angstroms on all sides
    
        RETURNS:
            (string) name of the selection created
    
        NOTES:
           assumes a single-state setup; for multi-state you must modify this script
        """
    
        # get the extents of the bounding object into two arrays
        (bMin, bMax) = cmd.get_extent(bounding_object)
    
        # if the user supplied a size by which we need to expand the box, we do that here
        if radius!=None:
            try:
                # try to convert a string to a float
                radius = float(str(radius))
            except ValueError:
                print "ERROR: Invalid radius declared.  Please make sure the radius is a number or left blank"
                return
            # add/subtract the radius to the bounds
            bMin = map( lambda f: f-radius, bMin)
            bMax = map( lambda f: f+radius, bMax)
     
        # initialize a selection; selections starting with underscores are hidden by PyMOL
        cmd.select("__toRemove", "None")
    
        # for each atom in the target selection, if it's inside the box, add it to the selection
        for at in cmd.get_model(target).atom:
            if at.coord[0] > bMin[0] and at.coord[0] < bMax[0] and at.coord[1] > bMin[1] and at.coord[1] < bMax[1] and at.coord[2] > bMin[2] and at.coord[2] < bMax[2] :
               cmd.select("__toRemove", "__toRemove or index %d" % at.index)
    
        # prepare the return values
        rName = "insideBoundingObj"
        cmd.set_name("__toRemove", rName)
    
        return rName
    
    cmd.extend("selInside", selInside)
    

Retrieved from "[https://pymolwiki.org/index.php?title=SelInside&oldid=7671](https://pymolwiki.org/index.php?title=SelInside&oldid=7671)"


---

## Set toggle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Syntax
  * 3 Examples
  * 4 Script
    * 4.1 See Also



## Overview

Set_toggle binds a key to switch on and off a setting. 

## Syntax
    
    
    set_toggle key, setting
    

See [Set_Key](/index.php/Set_Key "Set Key") to find a list of key that can be redefined and [Settings](/index.php/Settings "Settings") to find a list of settings. You can also use auto-completion. 

## Examples
    
    
    # Associates "F1" to grid_mode setting
    set_toggle F1, grid_mode
    

## Script
    
    
    # @AUTHOR: Joseph ANDRE
    # Copyright (c) 2011, Joseph ANDRE
    # All rights reserved.
    #
    # Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
    # conditions are met:
    #
    #     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following
    #     * disclaimer.
    #     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    #     * disclaimer in the documentation and/or other materials provided with the distribution.
    #     * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived
    #     * from this software without specific prior written permission.
    #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
    # NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
    # THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
    # OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    #
    # DATE  : 2011-07-14
    # REV   : 1
    
    from pymol import cmd
    
    def toggle(setting=[]):
    	setting_value=cmd.get(setting)
    	if cmd.toggle_dict.get(setting_value, 0):
    		cmd.unset(setting)
    	else:
    		cmd.set(setting)
    
    def set_toggle(key, setting):
    	cmd.set_key(key, toggle, [setting])
    
    cmd.auto_arg[1].update({
       'set_toggle'         : cmd.auto_arg[0]['set'],
    })
    	
    cmd.extend('set_toggle', set_toggle)
    

### See Also

[Set_Key ](/index.php/Set_Key "Set Key")

Retrieved from "[https://pymolwiki.org/index.php?title=Set_toggle&oldid=9187](https://pymolwiki.org/index.php?title=Set_toggle&oldid=9187)"


---

## Show aromatics

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    #
    # Show side chain sticks for aromatic residues
    #
    # (script for PyMol)
    # Tom Stout, 08/05/2004
    #
    
    mstop
    dss
    hide all
    #zoom all
    #orient
    show cartoon,all
    color gray,all
    
    select aromatics,(resn phe+tyr+trp+his)
    show sticks, (aromatics and (!name c+n+o))
    color green,aromatics
    disable aromatics
    set cartoon_smooth_loops,0
    

Retrieved from "[https://pymolwiki.org/index.php?title=Show_aromatics&oldid=6381](https://pymolwiki.org/index.php?title=Show_aromatics&oldid=6381)"


---

## Show bumps

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/show_bumps.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/show_bumps.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The red discs in the image below show vdW overlaps or steric clashing. This script will create those bumps for your chosen molecule. 

[![Sb.png](/images/a/a3/Sb.png)](/index.php/File:Sb.png)

This is a script to visualize VDW clashes, just like the [mutagenesis](/index.php/Mutagenesis "Mutagenesis") wizard does. 

Posted on [pymol-users mailing list](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg09190.html). 

## Usage
    
    
    show_bumps [ selection [, name ]]
    

## Related Settings

  * sculpt_vdw_vis_max (default: 0.3)
  * sculpt_vdw_vis_mid (default: 0.1)
  * sculpt_vdw_vis_min (default: -0.1)



Retrieved from "[https://pymolwiki.org/index.php?title=Show_bumps&oldid=13928](https://pymolwiki.org/index.php?title=Show_bumps&oldid=13928)"


---

## Show charged

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    # Show side chain sticks for charged residues
    #
    # (script for PyMol)
    # Tom Stout, 08/05/2004
    #
    
    mstop
    dss
    hide all
    #zoom all
    #orient
    show cartoon,all
    color gray,all
    
    select pos,(resn arg+lys+his)
    show sticks, (pos and !name c+n+o)
    color marine,pos
    disable pos
    select neg,(resn glu+asp)
    show sticks, (neg and !name c+n+o)
    color red,neg
    disable neg
    set cartoon_smooth_loops,0
    

Retrieved from "[https://pymolwiki.org/index.php?title=Show_charged&oldid=6383](https://pymolwiki.org/index.php?title=Show_charged&oldid=6383)"


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

## Show hydrophilic

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    # Show side chain sticks for hydrophilic residues
    #
    # (script for PyMol)
    # Tom Stout, 08/05/2004
    #
    
    mstop
    dss
    hide all
    #zoom all
    #orient
    show cartoon,all
    color gray,all
    
    select hydrophilics,(resn arg+lys+his+glu+asp+asn+gln+thr+ser+cys)
    show sticks, (hydrophilics and !name c+n+o)
    color green,hydrophilics
    disable hydrophilics
    set cartoon_smooth_loops,0
    

Retrieved from "[https://pymolwiki.org/index.php?title=Show_hydrophilic&oldid=6384](https://pymolwiki.org/index.php?title=Show_hydrophilic&oldid=6384)"


---

## Show hydrophobics

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    #
    #
    #
    # Show side chain sticks for hydrophobic residues
    #
    # (script for PyMol)
    # Tom Stout, 08/05/2004
    #
    
    mstop
    dss
    hide all
    #zoom all
    #orient
    show cartoon,all
    color gray,all
    
    select hydrophobes,(resn ala+gly+val+ile+leu+phe+met)
    show sticks, (hydrophobes and (!name c+n+o))
    color orange,hydrophobes
    disable hydrophobes
    set cartoon_smooth_loops,0
    

Retrieved from "[https://pymolwiki.org/index.php?title=Show_hydrophobics&oldid=6382](https://pymolwiki.org/index.php?title=Show_hydrophobics&oldid=6382)"


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

## Slerpy

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Manual for Slerpy.py
    * 1.1 General Use
    * 1.2 Important concepts
    * 1.3 Quick Start Tutorial
    * 1.4 Tips and Tricks
      * 1.4.1 Converting scenes to movies
      * 1.4.2 Starting off right
      * 1.4.3 Live pymol presentations during a talk
      * 1.4.4 Pausing on a view
      * 1.4.5 Morphing and multi-state models
    * 1.5 Command Reference
    * 1.6 Why is it called slerpy?
  * 2 Script Files
    * 2.1 Installation
    * 2.2 Code
      * 2.2.1 slerpy.py
      * 2.2.2 movs.py



# Manual for Slerpy.py

An extension to pymol that creates a moderately easy to use environment for doing keyframe animation. 

  


## General Use

At the pymol command line type: 
    
    
    import slerpy
    

This will load the extended commands. All commands in slerpy begin with the letter s. Pymol's tab autocomplete feature will work on the additional commands. 

## Important concepts

The main function of slerpy is to record a series of pymol views. A movie can then be created by interpolating between these views. A pymol view consists mostly of the camera orientation (that is, the orientation of the viewers eye with respect to the molecule). It also includes information about the clipping planes. 

It is important to realize that most slerpy commands act on the "current" view. You can navigate among views by using the sn, sp, and sgo commands. If you've changed the view around with the mouse or you just want to know the number of the current view you can get back to the current view with the sc command. 

Pymol views do not contain information about how pymol objects and selections are displayed. If you just want to create a movie which moves around a single representation of a molecule then all you need to do is record the set of views that define the tour of the molecule. 

If, on the other hand, you want to change the representation or change which items are displayed you will need to add actions to some of your views. An action is any set of pymol commands. Actions can be associated with any view in the series recorded by slerpy. The _sscene_ command inserts a new view and simultaneously creates a pymol scene and the action to display it. The scene will include all of the objects and representations visible at the time the command was issued. 

In order to control the rate of motion between the defined views in a slerpy movie, you can control the number of frames used in each interpolation. When a view is saved in slerpy it is associated by default with a transition of 50 frames to the next view. The number of frames in the transition can be altered with the _ssetf_ command in slerpy. 

The views and actions stored by slerpy can (and should) be saved to a key file with the _swrite_ command. They can then be retrieved with the sread command. Note that swrite saves the current pymol state in a .pse file but sread does not read in the .pse file. If you're starting a new pymol session to continue work on an existing movie you should load the pse file before doing an sread. 

## Quick Start Tutorial

If you haven't installed slerpy yet see Installation

Step 1
    You probably want to start off in a nice clean working directory that just has the coordinate files you want to work with.

Step 2
    Read in your molecule(s) and create the various selections and representations that you want to include in the movie.

Step 3
    At the pymol prompt, type:
    
    
    import slerpy
    

Step 4
    Get your molecule in exactly the orientation and representation that you want to use for the beginning of your movie.

Step 5
    Type:
    
    
    sinsert
    

Step 6
    Using the mouse, move your molecule to the next orientation that you want to use. When you record the movie, the camera orientation will be interpolated between each consecutive pair of views. This can include changes in rotation, zooming, clipping etc.

    Loop back to Step 5. Continue this until you've got all your orientations stored.

    You can check how any set of transitions will look at any time by using the _sshow_ command (see Command Reference for details).

    You can adjust the rate of any transition using the _ssetf_ command

Step 7A
    Add any actions to your views using the saction command. Any number of pymol commands can be strung together separated by semicolons. If, for example, you want to change your protein from a cartoon to a surface and add a ligand when you get to view 5 you would do the following (assuming you've defined the pymol selections prot and lig):
    
    
    sgo 5
    saction "hide cartoon, prot; show surface, prot; show sticks, lig"
    

Step 7B (Alternative using scenes)
    
    
    
    sgo 5
    

    Now use the gui to create the representation you want and then:
    
    
    sscene
    

Step 8
    Save everything. Save your slerpy views and actions as well as a pse file of your current pymol session:
    
    
    swrite mymovie
    

This will create mymovie.key, which has all of the views, frame counts, actions etc. and mymovie.pse, the associated pymol session file. 

Step 9
    Record the movie! Type:
    
    
    srecord
    

    You can then play the movie by typing the standard pymol command _mplay_ or by clicking the play button in pymol.

Step 10
    If you save the pymol session again, the pse file will contain the movie which can then be shown immediately after startup without running slerpy.py. Note that pymol will warn you when you load a pse file that contains a movie.

Step 11
    If you want to, the movie can be exported using the mpng command (see the pymol [documentation](http://pymol.sourceforge.net/newman/ref/S1000comref.html#2_105)). Also, see the useful article [Making Movies](/index.php/Making_Movies "Making Movies").

## Tips and Tricks

### Converting scenes to movies

You can just step through the scenes and type sscene for each one. This will create a duplicate slerpy scene for each of the scenes you'd already saved but that's not such a disaster. Be sure to swrite when you're done. 

Note that there is some overhead associated with recalling scenes. To avoid pauses at view transitions, I prefer to actually issue the set of show and hide commands that will generate the scene rather than using the above method. 

### Starting off right

It's a bit of a pain, but I like to associate with the first frame of the movie an action list that hides everything and then turns on all the objects that I want to have visible at the beginning. This ensures that when your movie loops back to the beginning it will look the same as it did the first time through. For example: 
    
    
    sgo 0
    saction "hide everything; show lines, prot; show surface, activesite; show sticks, ligand"
    

Alternatively, start your slerpy work with an _sscene_. 

### Live pymol presentations during a talk

Be sure to run your movie once it's been opened and before your presentation if you're presenting in pymol. This will ensure that any objects that don't appear until the middle of the movie are available in memory and won't need to be rebuilt while your audience waits. 

Of course showing your movie from within pymol allows you to show movies in stereo if you've got a presentation system that allows this. If you pass out stereo glasses first it's guaranteed that everyone will remember your talk... 

### Pausing on a view

Just sgo to the view you want to stay on for a while and do an sinsert. This will insert a new view with the same orientation etc as the one you were just on. You can adjust the length of the pause by changing the number of frames for the transistion between these two identical views using the ssetf command. 

### Morphing and multi-state models

Morphs and multi-conformation models are represented in pymol as single objects with multiple states. Cycling through the states as an animation is very straightforward in pymol. Slerpy allows you to include this animation over the states of an object as a transition between slerpy views within the context of larger movie. This is done using the _smorph_ command (see Command Reference). _smorph_ allows you to specify the starting and ending pymol state numbers to use during the transition from the current to the next view. It will often be appropriate to have the movie continue after the morph using the final state of the morphing model, that is, the conformation to which you morphed. Since the default for a typical slerpy view is to show the first state of an object, you'll probably need to have available, in addition to your multi-state model, a single-state object containing the final conformation. You can then hide the multistate object and show the single-state final conformation as an action associated with the final view of your morphing sequence. 

To generate morphs you can use the [rigimol](http://www.delanoscientific.com/rigimol.html) program provided by Warren as part of the incentive pymol package for subscribers or use [lsqman](http://alpha2.bmc.uu.se/usf/mol_morph.html) from Gerard Kleywegt. 

## Command Reference

Note that it is essential to understand that slerpy uses the concept of a current or active view. This is the element in the list to which most commands are applied. It is not necessarily the view that is currently visible on the screen. It is advisable to use the sc command frequently to make sure you really know which view your command is being applied to. 

**saction _string_ :** Assoiciates the pymol commands in _string_ with the current view. _string_ must be enclosed in quotes and must contain only valid pymol commands separated by semicolons. 

**sappend** : Add the currently displayed view at the end of the list of views in slerpy 

**sappendaction _string_ :** Adds the action(s) in _string_ to the list of actions associated with the current view 

**sc** : Go to the slerpy active view 

**scrossfade _startobject, endobject_ :** During the transition from the current view to the next view _startobject_ will fade out and _endobject_ will fade in. The two objects must be shown as sticks. They must be objects, not merely selections, as of pymol 0.99. 

**sdelaction** : Deletes the action associated with the current view. Be sure you're sure which view you're on before you use this. This will also delete any actions etc. associated with the view so be careful. 

**sdelete** : Remove the slerpy active view from the list of views to be interpolated. 

**sdeletesetting:** Remove the setting interpolation for the current view. 

**sdumpactions:** List all actions by frame. 

**sgo _n_ :** Change the slerpy active view to view _n_. 

**sinsert** : Add the currently displayed view after the slerpy active view. 

**sinterpsetting _setting, selection, startval, endval_ :** The pymol setting _setting_ will be interpolated linearly from _startval_ to _endval_ during the transition from the current view to the next view. You can only interpolate one setting per transition. This is the hard way to, for example, fade out an object: 
    
    
    sinterpsetting stick_transparency, lig, 0.0, 1.0
    

**slist** : List all views stored for interpolation. Also lists the number of frames in each transition. 

**smorph _startmodel,endmodel_ :** The transition from the current frame to the next frame will consist of one frame per pymol state starting with state _startmodel_ and ending with state _endmodel_. Subsequent frames (i.e. from subsequent views) will revert to state 1. The state numbers apply to currently visible objects so you will most likely want to have an object with your starting conformation, an object with your multi-state morphing model, and an object with your final conformation. You would then sgo to the frame where you want the morph to start and add an action to hide the starting conformation object and show the multi-model morphing object, do an smorph 1,30 or whatever the number of states in your morph is, append another frame and give it an action where the multi-state model is hidden and the final conformation is shown. 

**sn** : Go to the next view 

**sp** : Go to the previous view. 

**sread _filename_ :** Restore all the information written with swrite. Does not read in the pse file (to avoid inadvertantly writing over some new selections, scenes or whatever). 

**srecord** : Records the movie 

**sreplace** : Replace the slerpy current view with the currently displayed view. 

**sscene** : Add a the currently displayed view after the slerpy active view and create a scene to go with it. The current pymol state, including which objects are displayed and how they are shown will be captured. 

**ssetf** : Set the number of frames to use in the transition from the slerpy active view to the next view 

**ssetupview _n_ :** This attempts to make sure that all objects are displayed (or not) as they would be when view _n_ is arrived at in the movie. It goes through and executes all of the sactions from all preceeding views. For some reason this doesn't always work in complex cases. 

**sshow _n,m_ :** This command records and shows a segment of the movie showing the transitions starting with view n and ending with view m. If the arguments m and n are omitted the transition from the current view to the next view will be shown. 

**swrite _filename_ :** Writes all the information used by slerpy to a file _filename_.key. It also writes a pymol session file _filename_.pse containing all of your current objects, selections etc. The format of the .key file is supposed to be human readable and it can be convenient to edit this file rather than using saction over and over. After editing it by hand be sure to do an sread. 

## Why is it called slerpy?

slerp is an acronym of sorts for spherical linear interpolation. The transition between two views approximates spherical linear interpolation of the their camera positions using quaternions. See the Wikipedia [article](http://en.wikipedia.org/wiki/Slerp) on slerps. 

# Script Files

## Installation

To install, copy (by pasting into a text editor) the following two python code segments to files named slerpy.py and movs.py, respectively. Place these files in the pymol/modules/pymol/ directory created by the pymol installation. 

If you don't have write priviledges in the pymol installation directories you can just copy these files to the working directory from which pymol will be run. 

## Code

### slerpy.py

    The command definitions for slerpy
    
    
    ################################################################################
    #slerpy.py - Command definition routines for slerpy
    #Copyright (C) 2006 Joel Bard
    #
    #This program is free software; you can redistribute it and/or
    #modify it under the terms of the GNU General Public License
    #as published by the Free Software Foundation; either version 2
    #of the License, or (at your option) any later version.
    #
    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.
    #
    #You should have received a copy of the GNU General Public License
    #along with this program; if not, write to the Free Software
    #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
    #################################################################################
    
    from pymol import cmd
    import movs
    
    def readViews( filename ):
    #obsolete
        vfile = open( filename, 'r')
        views = [] #a list of views each of which is a list of 18 numbers
        vstrings = vfile.readlines()
        for vstring in vstrings:
            vals = vstring.split()
            view = [float(v) for v in vals]
            views.append( view )
        vfile.close()
        return views
    
    def readFrames( filename ):
    #obsolete
        ffile = open( filename, 'r' )
        frames = []
        fstrings = ffile.readlines()
        for fstring in fstrings:
            frames.append( int(fstring) )
        ffile.close()
        return frames
    
    def readActions( filename ):
    #obsolete
        #actions are stored in file where
        #for each line, the first 4 chars are the view index
        #associated with the action and the rest of the
        #line is the pymol command to be executed
        #upon reading, a dictionary is returned
        #with view indices as keys and actions as values
        actions = {}
        try:
            afile = open( filename, 'r' )
        except:
            print "No actions for this project"
            return actions
        astrings = afile.readlines()
        for astring in astrings:
            try:
                aindex = int(astring[:4])
                action = astring[4:]
                actions[ aindex ] = action[:-1]
            except:
                print "empty line"
        afile.close()
        return actions
    
    def readModels( filename ):
    #obsolete
        models = {}
        try:
            mfile = open( filename, 'r' )
        except:
            print "No models for this project"
            return models
        mstrings = mfile.readlines()
        for mstring in mstrings:
            try:
                mindex = int(mstring[:4])
                model = mstring[4:]
                models[ mindex ] = model[:-1]
            except:
                print "empty line"
        mfile.close()
        return models
    
    def readSettings( filename ):
    #obsolete
        settings = {}
        try:
            sfile = open( filename, 'r' )
        except:
            print "No settings for this project"
            return settings
        sstrings = sfile.readlines()
        for sstring in sstrings:
            try:
                sindex = int(sstring[:4])
                scommas = sstring[4:]
                settingName,selection,startVal,endVal = scommas.split(',')
                setting = [settingName,selection,float(startVal),float(endVal)]
                settings[sindex] = setting
            except:
                print "unable to parse setting"
        sfile.close()
        return settings
                
    def readScenes( filename ):
    #obsolete
        global scene_counter
    
        scene_counter = 0
        scenes = {}
        try:
            sfile = open( filename, 'r' )
        except:
            print "No scenes file for this project"
            return scenes
        sstrings = sfile.readlines()
        for sstring in sstrings:
            try:
                sindex = int(sstring[:4])
                scene = sstring[4:]
                scenes[ sindex ] = scene[:-1]
                scene_counter += 1
                #print "reading scene", sstring, sindex, scene
            except:
                print "empty line"
        sfile.close()
        return scenes
    
    def read_all( fileroot ):
    #obsolete in favor of readKeyViewFile
        global views
        global frames
        global actions
        global cur_view
        global cur_index
        global scenes
        global models
        global settings
        
        views = readViews( fileroot+".txt" )
        frames = readFrames( fileroot+".frm")
        actions = readActions( fileroot+".act")
        scenes = readScenes( fileroot+".scn")
        models = readModels( fileroot+".mod")
        settings = readSettings( fileroot+".set")
        cur_view = views[0]
        cur_index = 0
        show_cur()
        
    def print_view( view ):
        for i in range(0,6):
            for j in range(0,3):
                print "%12.6f"% view[ 3*i+j ] ,
            print
    
    def show_next():
        global cur_index
        global cur_view
        global views
        if cur_index == len( views )-1:
            print "No more views."
            return
        cur_index += 1
        cur_view = views[ cur_index ]
        cmd.set_view(cur_view)
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ",cur_index
    
    def show_prev():
        global cur_index
        global cur_view
        global views
        if cur_index == 0:
            print "No more views."
            return
        cur_index -= 1
        cur_view = views[ cur_index ]
        cmd.set_view(cur_view)
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ",cur_index
    
    def show_cur():
        global cur_index
        global cur_view
        global views
        cur_view = views[ cur_index ]
        cmd.set_view(cur_view)
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ",cur_index
    
    def go_to_view( arg=0 ):
        global cur_index
        global cur_view
        global views
        n = int( arg )
        if n < 0 or n >= len(views):
            print "Index out of range."
            return
        cur_index = n
        cur_view = views[n]
        cmd.set_view( cur_view )
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ", cur_index
    
    def insert_current():
        #insert the current view into the list after the view
        #in views[cur_index]
        #set frames to default
        global cur_index
        global cur_view
        global views
        global frames
        global actions
        global scenes
        global settings
        global models
        global fades
        
        cur_index += 1
        cur_view = cmd.get_view()
        views.insert( cur_index, [cv for cv in cur_view] )
        frames.insert( cur_index, 50 )
        
        #deal with actions dictionary
        actions = incKeyAbove( actions, cur_index )
        scenes = incKeyAbove( scenes, cur_index )
        settings = incKeyAbove( settings, cur_index )
        models = incKeyAbove( models, cur_index )
        fades = incKeyAbove( fades, cur_index )
    
        print "New view:"
        print_view( cur_view )
        print "Inserted view with index", cur_index, "and a 50 frame transition"
    
    def append_view():
        global views
        global frames
        global cur_index
        global cur_view
    
        cur_index = len(views)
        cur_view = cmd.get_view()
        views.append( [cv for cv in cur_view] )
        frames.append( 50 )
    
        print "New view: "
        print_view( cur_view )
        print "Appended view with index", cur_index, "and a 50 frame transition"
        print "The current view is", cur_index
        
    def incKeyAbove( dict, index ):
        tempDict = {}
        for key, val in dict.iteritems():
            if key >= index:
                newkey = key + 1
                tempDict[newkey] = val
            else:
                tempDict[key] = val
        return tempDict
    
    def decKeyAbove( dict, index ):
        tempDict = {}
        for key, val in dict.iteritems():
            if key > index:
                newkey = key - 1
                tempDict[newkey] = val
            else:
                tempDict[key] = val
        return tempDict
        
    def delete_current():
        #remove the current view from the list
        #show the previous view
        global cur_index
        global cur_view
        global views
        global actions
        global scenes
        global settings
        global models
        global frames
        global fades
        
        del views[cur_index]
        del frames[cur_index]
        if cur_index in actions:
            del actions[cur_index]
        if cur_index in scenes:
            del scenes[cur_index]
        if cur_index in settings:
            del settings[cur_index]
            
        #deal with dictionaries
        actions = decKeyAbove( actions, cur_index )
        scenes = decKeyAbove( scenes, cur_index )
        settings = decKeyAbove( settings, cur_index )
        models = decKeyAbove( models, cur_index )                               
        fades = decKeyAbove( fades, cur_index )
    
        print "View number",cur_index,"deleted."
        if cur_index > 0:
            cur_index -= 1
        cur_view = views[cur_index]
        cmd.set_view( cur_view )
        print "Current view is number",cur_index
    
    def delete_settings():
        global settings
        global cur_index
        del settings[cur_index]
    
    def replace_current():
        global cur_index
        global cur_view
        global views
        cur_view = cmd.get_view()
        views[cur_index] = [cv for cv in cur_view]
    
    def insert_scene():
        global views
        global actions
        global settings
        global frames
        global cur_index
        global cur_view
        global scenes
        global scene_counter
        global models
        global fades
     
        cur_index += 1
        cur_view = cmd.get_view()
        views.insert( cur_index, [cv for cv in cur_view] )
        frames.insert( cur_index, 50 )
        
        #deal with dictionaries
        actions = incKeyAbove( actions, cur_index )
        scenes = incKeyAbove( scenes, cur_index )
        settings = incKeyAbove( settings, cur_index )
        models = incKeyAbove( models, cur_index )
        fades = incKeyAbove( fades, cur_index )
    
        #this stuff has to be done after the above
        #find a free scene name
        i = 1
        while True:
            for sname in scenes.values():
                print "|"+sname+"|"
                if sname == "slerpy_"+str(i):
                    break
            else:
                break
            i += 1
        newname = "slerpy_"+str(i)
            
        scene_counter += 1
        cmd.scene( newname, "store" )
        scenes[ cur_index ] = newname
        actions[ cur_index ] = "scene "+newname
        
        print "New view:"
        print_view( cur_view )
        print "Inserted view with index", cur_index, "and a 50 frame transition"
        print "Added scene",newname
        
    def write_views( filename ):
    #deprecated in favor of key files
        global views
        global frames
        global scenes
        global actions
        global settings
        global models
        global fades
    
        viewfile = open( filename+".txt", 'w')
        for view in views:
            for v in view:
                viewfile.write( str(v) + " " )
            viewfile.write('\n')
        viewfile.close()
        
        framefile = open( filename+".frm", 'w' )
        for frame in frames:
            framefile.write( str( frame )+'\n' )
        framefile.close()
    
        actionfile = open( filename+".act", 'w' )
        for key,action in actions.iteritems():
            keystring = str( key )
            actionfile.write( keystring.rjust(4)+action + '\n' )
        actionfile.close()
    
        scenefile = open( filename+".scn", 'w' )
        for key,scene in scenes.iteritems():
            keystring = str( key )
            scenefile.write( keystring.rjust(4)+scene + '\n' )
        scenefile.close()
    
        modelfile = open( filename+".mod", 'w' )
        for key,model in models.iteritems():
            keystring = str( key )
            modelfile.write( keystring.rjust(4)+model + '\n' )
        modelfile.close()
    
        settingsfile = open( filename+".set", 'w' )
        for key,setting in settings.iteritems():
            keystring = str( key )
            settingName, selection, startVal, endVal = setting
            settingsfile.write( "%s%s,%s,%f,%f\n" % (keystring.rjust(4), settingName, selection, startVal, endVal))
        settingsfile.close()
        cmd.save( filename+".pse" )
            
        print "Wrote files", filename+".txt,",filename+".frm,",filename+".pse, and",filename+".act"
    
    def writeKeyViewFile( filename ):
        global views
        global frames
        global scenes
        global actions
        global settings
        global models
        global fades
     
        keyviewfile = open( filename + ".key", 'w' )
        for i,view in enumerate(views):
            keyviewfile.write( "VIEW: %4d " % i )
            for v in view:
                keyviewfile.write( str(v) + " " )
            keyviewfile.write('\n')
            keyviewfile.write( "FRAMES: %d\n" % frames[i] )
            if i in actions:
                keyviewfile.write( "ACTIONS: %s\n" % actions[i] )
            if i in scenes:
                keyviewfile.write( "SCENES: %s\n" % scenes[i] )
            if i in models:
                keyviewfile.write( "MODELS: %s\n" % models[i] )
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                keyviewfile.write( "SETTINGS: %s, %s, %f, %f\n" % (settingName, selection, startVal, endVal))
            if i in fades:
                startVisSelection, endVisSelection, sticksOnly = fades[i]
                keyviewfile.write( "FADES: %s, %s, %d\n" % (startVisSelection, endVisSelection, sticksOnly))            
            keyviewfile.write("\n")
        keyviewfile.close()
        cmd.save( filename + ".pse" )
        print "Wrote files " , filename + ".key", filename + ".pse"
    
    def readKeyViewFile( filename ):
        global views
        global frames
        global scenes
        global actions
        global settings
        global models
        global fades
        global scene_counter
    
        views = []
        frames = []
        actions = {}
        scenes = {}
        models = {}
        settings = {}
        fades = {}
        scene_counter = 0
    
        if filename.endswith(".key"): filename = filename[:-4]
        keyviewfile = open( filename + ".key", 'r' )
        viewstrings = keyviewfile.readlines()
        keyviewfile.close()
        viewindex = 0
        for line in viewstrings:
            if line.startswith("VIEW: "):
                viewindex = int( line[6:10] )
                vals = line[10:].split()
                view = [float(v) for v in vals]
                views.append( view )
            if line.startswith("FRAMES: "):
                frames.append( int( line[8:] ) )
            if line.startswith("ACTIONS: "):
                actions[ viewindex ] = line[9:-1]
            if line.startswith("SCENES: "):
                scenes[ viewindex ] = line[8:-1]
                scene_counter += 1
            if line.startswith("MODELS: "):
                models[ viewindex ] = line[8:-1]
            if line.startswith("SETTINGS: "):
                settingName,selection,startVal,endVal = line[10:-1].split(',')
                settings[ viewindex ] = [settingName,selection,float(startVal),float(endVal)]
            if line.startswith( "FADES: " ):
                startVisSelection, endVisSelection, sticksOnly = line[7:-1].split(',')
                fades[ viewindex ] = [startVisSelection, endVisSelection, int(sticksOnly) ]
        cur_view = views[0]
        cur_index = 0
        show_cur()
    
    def set_frames_current( arg=50 ):
        global frames
        global cur_index
        frames[cur_index] = int(arg)
    
    def list_frames():
        global frames
        global views
        global actions
        global models
        global settings
        
        f=0
        for i,view in enumerate(views[:-1]):
            if i in actions:
                a = actions[i]
            else:
                a = ""
            if i in models:
                m = "States: " + models[i]
            else:
                m = ""
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                s = "Settings: %s %s %f %f" % (settingName, selection, startVal, endVal)
            else:
                s = ""
            print "View",i,"to",i+1,"Frames ",f,"to",f+frames[i],a,m,s
            f += frames[i]
    
    def add_action_current( cmd ):
        global cur_index
        global actions
        actions[cur_index] = cmd[1:-1] #strip off quotes
    
    def append_action_current( cmd ):
        global cur_index
        global actions
        actions[cur_index] += ";" + cmd[1:-1]
    
    def clear_action_current():
        global actions
        global cur_index
        del actions[cur_index]
    
    def list_actions():
        global actions
        for i,a in actions.iteritems():
            print i,a
    
    def morph_models( start_model, end_model ):
        global cur_index
        global frames
        global models
        models[cur_index] = "%s -%s" % (start_model, end_model)
        frames[cur_index] = abs(int(end_model) - int(start_model)) + 1
    
    def interpolate_settings( setting, selection, startval, endval ):
        global cur_index
        global settings
        settingEntry = [setting, selection, float(startval), float(endval)]
        settings[cur_index] = settingEntry 
    
    def crossfade( startVisSelection, endVisSelection, sticksOnly = 1 ):
    #cross fade the specified objects, vary stick transparency only if stickOnly=1
        global cur_index
        global fades
        fades[cur_index] = [startVisSelection, endVisSelection, int(sticksOnly) ]
    
    def setup_view( index ):
        for i in range( int(index)+1 ):
            if i in actions:
                print "Executing %s from actions %d" % (actions[i],i)
                cmd.do( actions[i] )
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                action = "set %s, %f, %s;" % (settingName, endVal, selection)
                print "Executing %s from settings %d" % (action,i)
                cmd.do( action )
            if i in fades:
                startVisSelection, endVisSelection, sticksOnly = fades[i]
                action = "set stick_transparency, 0, %s; set stick_transparency, 1, %s;" % (endVisSelection, startVisSelection)
                print "Executing %s from fades %d" % (action, i)
                cmd.do( action )
    
    def show_transition(start_index=0, end_index=0):
        #show the transition from the current view to the next view
        global frames
        global views
        global cur_index
        global actions
        global models
        if start_index == 0 and end_index == 0:
            if cur_index >= len(views)-1:
                print "Current view is last in sequence."
                return
            start_index=cur_index
            end_index=cur_index+1
        else:
            start_index = int(start_index)
            end_index = int(end_index)
        ftot = 0
        setcommand = ""
        i = start_index
        for nframes in frames[start_index:end_index]:
            #ftot += nframes
            if i in models:
                setcommand += " " + models[i] + " "
            else:
                setcommand += " 1 x%i" % nframes
            i += 1
            
    #    cmd.mset("1 x%i" % ftot)
        cmd.mset( setcommand )
        start_frame = 1
        #first do all actions that happen up to this point to make sure
        #things look the way they should
        first_action = ""
        for i in range( start_index ):
            if i in actions:
                first_action += actions[i] + ';'
                #print "Executing %s from actions %d" % (actions[i],i)
                #cmd.do( actions[i] )
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                action = "set %s, %f, %s;" % (settingName, endVal, selection)
                first_action += action
                #print "Executing %s from settings %d" % (action,i)
                #cmd.do( action )
            if i in fades:
                startVisSelection, endVisSelection, sticksOnly = fades[i]
                action = "set stick_transparency, 0, %s; set stick_transparency, 1, %s;" % (endVisSelection, startVisSelection)
                first_action += action
                #print "Executing %s from fades %d" % (action, i)
                #cmd.do( action )
        for i in range( start_index, end_index ):
            if i in settings:
                movs.animate_transition( views[i], views[i+1], frames[i], start_frame, settings[i] )
            elif i in fades:
                movs.animate_transition( views[i], views[i+1], frames[i], start_frame, fades[i] )
            else:
                movs.animate_transition( views[i], views[i+1], frames[i], start_frame )
            #add an action
            if start_frame == 1:
                mdo_cmd = first_action
                if i in actions:
                    mdo_cmd += actions[i]+";"
                #mdo_cmd += "set_view("+str(views[i])+")"
                print mdo_cmd
                cmd.mdo(start_frame, mdo_cmd)
            elif i in actions:
                mdo_cmd = actions[i]+";set_view("+str(views[i])+")"
                cmd.mdo(start_frame, mdo_cmd)
                #print mdo_cmd
            start_frame += frames[i]
        cmd.frame(1)
        cmd.mplay()
    
    def make_all():
        #make the whole movie
        global views
        global frames
        global models
        
        #first get total number of frames
        ftot = 0
        setcommand = ""
        for i,nframes in enumerate(frames[:-1]):
            ftot += nframes
            if i in models:
                setcommand += " " + models[i] + " "
            else:
                setcommand += " 1 x%i" % nframes
    
        #initialize movie
        #cmd.mset("1 x%i" % ftot)
        #cmd.mset("1 x50 1 -30 30 x20")
        cmd.mset( setcommand )
    
        #loop through views
        start_view = views[0][:]
        first_frame = 1
        for i,view in enumerate(views[1:]):
            end_view = view[:]
            if i in settings:
                movs.animate_transition( start_view, end_view, frames[i], first_frame, settings[i] )
            elif i in fades:
                movs.animate_transition( start_view, end_view, frames[i], first_frame, fades[i] )
            else:
                movs.animate_transition( start_view, end_view, frames[i], first_frame )
            #add an action
            if i in actions:
                mdo_cmd = actions[i]#+";set_view ("+str( views[i] )+")"
                print mdo_cmd
                cmd.mdo(first_frame, mdo_cmd)
            first_frame += frames[i]
            start_view = end_view[:]
        cmd.frame(1)
    
    ## views = readViews( "viewfile.txt" )
    ## frames = readFrames( "viewfile.frm" )
    ## actions = readActions( "viewfile.act" )
    ##print "Length ",len(views)
    #for v in views:
     #   print v
    #cur_view = views[0]
    views = []
    frames = []
    models = {}
    actions = {}
    scenes = {}
    settings = {}
    fades = {}
    scene_counter = 0
    cur_index = -1
    cmd.set( "scene_animation_duration","0" )
    #cmd.set_view( cur_view )
    
    cmd.extend("sn", show_next )
    cmd.extend("sp", show_prev )
    cmd.extend("sc", show_cur )
    cmd.extend("sinsert", insert_current )
    cmd.extend("sdelete", delete_current )
    cmd.extend("sreplace", replace_current )
    cmd.extend("sappend", append_view )
    cmd.extend("sscene", insert_scene )
    cmd.extend("sgo", go_to_view )
    cmd.extend("sreadold", read_all )
    cmd.extend("swriteold", write_views )
    cmd.extend("slist", list_frames )
    cmd.extend("ssetf", set_frames_current )
    cmd.extend("sshow", show_transition )
    cmd.extend("srecord", make_all )
    cmd.extend("saction", add_action_current )
    cmd.extend("sdelaction", clear_action_current )
    cmd.extend("sdumpactions", list_actions )
    cmd.extend("sappendaction", append_action_current )
    cmd.extend("smorph", morph_models )
    cmd.extend("sinterpsetting", interpolate_settings )
    cmd.extend("sdeletesetting", delete_settings )
    cmd.extend("scrossfade", crossfade )
    cmd.extend("swrite", writeKeyViewFile )
    cmd.extend("sread", readKeyViewFile )
    cmd.extend("ssetupview", setup_view )
    

### movs.py

    Math and animation routines for slerpy
    
    
    ##################################################################################
    #movs.py - Math and animation routines for slerpy
    #Copyright (C) 2006 Joel Bard
    #
    #This program is free software; you can redistribute it and/or
    #modify it under the terms of the GNU General Public License
    #as published by the Free Software Foundation; either version 2
    #of the License, or (at your option) any later version.
    #
    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.
    #
    #You should have received a copy of the GNU General Public License
    #along with this program; if not, write to the Free Software
    #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
    #################################################################################
    
    from pymol import cmd,stored
    from math import *
    
    def rmat2quat( M ):
        #M is a list of 9 values being the elements of the rotation matrix in row order
        T = M[0] + M[4] + M[8] + 1
        print "Trace ",T
        if T>0:
            S = 0.5 / sqrt(T)
            W = 0.25/S
            X = (M[7] - M[5])*S
            Y = (M[2] - M[6])*S
            Z = (M[3] - M[1])*S
        elif M[0] > M[4] and M[0] > M[8]:
            S = sqrt( 1.0 + M[0] - M[4] - M[8]) * 2
            X = 0.25 * S
            Y = (M[1] + M[3])/S
            Z = (M[2] + M[6])/S
            W = (M[5] - M[7])/S
        elif M[4] > M[8]:
            S = sqrt( 1.0 + M[4] - M[0] - M[8] ) * 2
            X = (M[1] + M[3])/S
            Y = 0.25 * S
            Z = (M[5] + M[7])/S
            W = (M[2] - M[6])/S
        else:
            S = sqrt( 1.0 + M[8] - M[0] - M[4]) * 2
            X = (M[2] + M[6])/S
            Y = (M[5] + M[7])/S
            Z = 0.25 * S
            W = (M[1] - M[3])/S
        return [X,Y,Z,W]
    
    def quat2rmat( Q ):
        #Q is a list of 4 values being the quaternion X Y Z W
        X=Q[0]
        Y=Q[1]
        Z=Q[2]
        W=Q[3]
        xx = X*X
        xy = X*Y
        xz = X*Z
        xw = X*W
        yy = Y*Y
        yz = Y*Z
        yw = Y*W
        zz = Z*Z
        zw = Z*W
    
        M= [1.0]*9
        M[0] = 1 - 2 * ( yy + zz )
        M[1] = 2 * ( xy - zw )
        M[2] = 2 * ( xz + yw )
        M[3] = 2 * ( xy + zw )
        M[4] = 1 - 2 * ( xx + zz )
        M[5] = 2 * ( yz - xw )
        M[6] = 2 * ( xz - yw )
        M[7] = 2 * ( yz + xw )
        M[8] = 1 - 2 * ( xx + yy )
        return M
    
    def quatconj( Q ):
        return [-Q[0],-Q[1],-Q[2],Q[3]]
    
    def quatmag( Q ):
        s = 0.0
        QC = quatconj(Q)
        for x in range(4):
            s += Q[x]*Q[x]
        print s
        return sqrt(s)
    
    def quatnorm( Q ):
        m = quatmag( Q )
        return [q/m for q in Q]
    
    def quatdotprod( q1, q2 ):
        dp = 0
        for i in range(4):
            dp += q1[i]*q2[i]
        return dp
    
    def vectnorm( V ):
        mag = 0.0
        for x in V:
            mag += x*x
        mag = sqrt(mag)
        return [x/mag for x in V]
    
    def quat2axisangle( Q ):
        #returns list where 0..2 are rot axis and 3 is angle
        qn = quatnorm( Q )
        cos_a = Q[3]
        angle = acos( cos_a ) * 2
        sin_a = sqrt( 1.0 - cos_a * cos_a )
        if fabs( sin_a ) < 0.000005:
            sin_a = 1
        ax_an = [ q/sin_a for q in Q[0:3] ]
        ax_an.append( angle )
        return ax_an
    
    def axisangle2quat( ax_an ):
        #ax_an is a list with axis coordinates followed by rotation angle
        axn = vectnorm( ax_an[:3] )
        angle = ax_an[3]
        sin_a = sin( angle / 2 )
        cos_a = cos( angle / 2 )
        Q = [ x * sin_a for x in axn ]
        Q.append( cos_a )
        return Q
        
    def rmat2axisangle( M ):
        q = rmat2quat( M )
        return quat2axisangle( q )
    
    def axisangle2rmat( a ):
        q = axisangle2quat( a )
        return quat2rmat( q )
    
    def animate_transition( start_view, end_view, nframes, first_frame, settings = [] ):
        #print "Views"
        #print start_view,'\n',end_view
    
        cview = start_view[:]
        cmd.set_view( start_view )
    
        #get initial and final quaternions for interpolation
        #print start_view[0:9]
        #get quaternions
        qstart = rmat2quat( start_view[0:9] )
        qend = rmat2quat( end_view[0:9] )
        
        #test for long way vs. short way
        if quatdotprod( qstart,qend ) < 0:
            qstart = [-q for q in qstart]
    
        axan_start = quat2axisangle( qstart )
        axan_end = quat2axisangle( qend )
            
        axan_cur = axan_start[:]
        frame_start = first_frame
        frame_end = frame_start + nframes
        doFade = 0
        doSetting = 0
        if len( settings ) == 4:
            settingName, selection, startVal, endVal = settings
            settingStep = (endVal-startVal)/float(nframes)
            print "Setting step ", settingStep
            doSetting = 1
        elif len( settings ) == 3:
            startVisSelection, endVisSelection, sticksOnly = settings
            settingStep = 1.0/float(nframes)
            doFade = 1
        for f in range( frame_start , frame_end):
            #get rotmat
            #using angle axis
    
            for i in range(4):
                axan_cur[i] = axan_cur[i] + (axan_end[i]-axan_start[i])/nframes
            newmat = axisangle2rmat( axan_cur )
            #print cview
            for i in range(9):
                cview[i] = newmat[i]
    
            mdo_cmd = "set_view (["
            for i in range(18):
                if i>8:
                    cview[i] = cview[i]+(end_view[i]-start_view[i])/nframes
                mdo_cmd += "%12.7f,"% cview[i]
            mdo_cmd = mdo_cmd[:-1]+"])"
            if doSetting:       
                val = float(f-frame_start)*settingStep + startVal 
                print val;
                mdo_cmd += "; set %s, %f, %s" % (settingName, val, selection)
                print mdo_cmd;
            #print "mdo ", mdo_cmd
            if doFade:
                val = float(f-frame_start)*settingStep
                otherVal = 1.0 - val
                mdo_cmd += "; set stick_transparency, %f, %s; set stick_transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
                if not sticksOnly:
                    #comment out surface transparency interpolation due to problem with transparent sticks in front of 
                    #transparent surfaces (get black holes)
                   # mdo_cmd += "; set transparency, %f, %s; set transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
                    mdo_cmd += "; set cartoon_transparency, %f, %s; set cartoon_transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
            cmd.mdo(f,mdo_cmd)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Slerpy&oldid=8080](https://pymolwiki.org/index.php?title=Slerpy&oldid=8080)"


---

## Spectrum states

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/spectrum_states.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrum_states.py)  
Author(s)  | [Takanori Nakane](/index.php?title=User:TakanoriNakane&action=edit&redlink=1 "User:TakanoriNakane \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.viewing](http://pymol.org/psicoredirect.php?psico.viewing)  
  
**spectrum_states** colors each state in a multi-state object different. Colors are along a user defined palette, just like with the [spectrumany](/index.php/Spectrumany "Spectrumany") command. 

## Usage
    
    
    spectrum_states [ selection [, representations [, color_list ]]]
    

## Example
    
    
    import spectrum_states
    
    # fetch a morph from molmovdb.org
    load http://molmovdb.org/uploads/284066-6299/movie.pdb.gz
    
    # show ribbon in all states
    as ribbon
    set all_states
    
    # color states from dark gray to red
    spectrum_states *, ribbon, gray20 orange red
    

## See Also

  * [spectrumany](/index.php/Spectrumany "Spectrumany")
  * [split_states](/index.php/Split_states "Split states")
  * [color_obj](/index.php/Color_Objects "Color Objects")



Retrieved from "[https://pymolwiki.org/index.php?title=Spectrum_states&oldid=13929](https://pymolwiki.org/index.php?title=Spectrum_states&oldid=13929)"


---

## Spectrumany

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/spectrumany.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrumany.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.viewing](http://pymol.org/psicoredirect.php?psico.viewing)  
  
[![](/images/8/8b/SpectrumanyExample.png)](/index.php/File:SpectrumanyExample.png)

[](/index.php/File:SpectrumanyExample.png "Enlarge")

Coloring with a gradient of different green shades

This script works similar to the [spectrum](/index.php/Spectrum "Spectrum") command, but instead of predefined palettes, any color sequence can be used. 

The color sequence is given by a space separated list of colors, so palette "red_white_blue" is the same as color sequence "red white blue". 

# Example
    
    
    fetch 2x19
    
    # these two produce the same result
    spectrum count, red_white_blue, chain B
    spectrumany count, red white blue, chain B
    
    # gradient of different green colors
    spectrumany count, smudge palegreen limegreen limon green forest, chain B
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spectrumany&oldid=13930](https://pymolwiki.org/index.php?title=Spectrumany&oldid=13930)"


---

## Spectrumbar

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/spectrumbar.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/spectrumbar.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
    * 1.1 USAGE
    * 1.2 EXAMPLES
  * 2 Script



## Description

Often times one may color their structure based upon some predefined color spectrum. However, this alone is not helpful when dealing with publication images. The purpose of this script is to allow the user to draw their own custom spectrum bar to accompany their structural images. 

Further work will be done to improve the functionality of this script. Please feel free to contact the author for function requests. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "spectrumbar" (which can be changed through the options) which is a cylinder colored according to the user-specified colors. 

Note: It is strongly recommended to turn off the specular reflections before ray tracing 

### EXAMPLES

[![](/images/b/ba/Sb_red.png)](/index.php/File:Sb_red.png) [](/index.php/File:Sb_red.png "Enlarge")Fig.1 - spectrumbar red | [![](/images/b/ba/Sb_red.png)](/index.php/File:Sb_red.png) [](/index.php/File:Sb_red.png "Enlarge")Fig.2 - spectrumbar 1,0,0 | [![](/images/f/f5/Sb_red_orange_yellow.png)](/index.php/File:Sb_red_orange_yellow.png) [](/index.php/File:Sb_red_orange_yellow.png "Enlarge")Fig.3 - spectrumbar red, orange, yellow  
---|---|---  
[![](/images/0/00/Sb_default.png)](/index.php/File:Sb_default.png) [](/index.php/File:Sb_default.png "Enlarge")Fig.4 - spectrumbar red,orange,yellow,green,blue,purple | [![](/images/0/00/Sb_default.png)](/index.php/File:Sb_default.png) [](/index.php/File:Sb_default.png "Enlarge")Fig.5 - spectrumbar 1,0,0, orange, yellow, 0,1,0, 0,0,1, purple | [![](/images/5/58/Sb_radius.png)](/index.php/File:Sb_radius.png) [](/index.php/File:Sb_radius.png "Enlarge")Fig.6 - spectrum red,radius=2  
[![](/images/6/63/Sb_long_red_orange_long_yellow.png)](/index.php/File:Sb_long_red_orange_long_yellow.png) [](/index.php/File:Sb_long_red_orange_long_yellow.png "Enlarge")Fig 7 - spectrumbar red, red, red, orange, yellow, yellow | [![](/images/b/ba/Sb_rounded.png)](/index.php/File:Sb_rounded.png) [](/index.php/File:Sb_rounded.png "Enlarge")Fig 8 - spectrumbar red, radius=2, ends=rounded  
  
## Script

load the script using the [run](/index.php/Run "Run") command 
    
    
    from pymol.cgo import *
    from math import *
    from pymol import cmd
    from re import *
    
    def spectrumbar (*args, **kwargs):
    
        """
        Author Sean M. Law
        University of Michigan
        seanlaw_(at)_umich_dot_edu
    
        USAGE
    
        While in PyMOL
    
        run spectrumbar.py
    
        spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)
    
        Parameter     Preset         Type     Description
        RGB_Colors    [1.0,1.0,1.0]  N/A      RGB colors can be specified as a
                                              triplet RGB value or as PyMOL
                                              internal color name (i.e. red)
        radius        1.0            float    Radius of cylindrical spectrum bar
        name          spectrumbar    string   CGO object name for spectrum bar
        head          (0.0,0.0,0.0)  float    Starting coordinate for spectrum bar
        tail          (10.0,0.0,0.0) float    Ending coordinate for spectrum bar
        length        10.0           float    Length of spectrum bar
        ends          square         string   For rounded ends use ends=rounded
    
        Examples:
    
        spectrumbar red, green, blue
        spectrumbar 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0
    
        The above two examples produce the same spectrumbar!
    
        spectrumbar radius=5.0
        spectrumbar length=20.0
    
        """
    
        rgb=[1.0, 1.0, 1.0]
        name="spectrumbar"
        radius=1.0
        ends="square"
        x1=0
        y1=0
        z1=0
        x2=10
        y2=0
        z2=0
        num=re.compile('[0-9]')
        abc=re.compile('[a-z]')
    
        for key in kwargs:
            if (key == "radius"):
                radius = float(kwargs["radius"])
            elif (key == "name"):
                name=kwargs["name"]
            elif (key == "head"):
                head=kwargs["head"]
                head=head.strip('" []()')
                x1,y1,z1=map(float,head.split(','))
            elif (key == "tail"):
                tail=kwargs["tail"]
                tail=tail.strip('" []()')
                x2,y2,z2=map(float,tail.split(','))
            elif (key == "length"):
                if (abc.search(kwargs["length"])):
                    print "Error: The length must be a value"
                    return
                else:
                    x2=float(kwargs["length"]);
            elif (key == "ends"):
                ends=kwargs["ends"]
            elif (key != "_self"):
                print "Ignoring unknown option \""+key+"\""
            else:
                continue
    
        args=list(args)
        if (len(args)>=1):
            rgb=[]
        while (len(args)>=1):
            if (num.search(args[0]) and abc.search(args[0])):
                if (str(cmd.get_color_tuple(args[0])) != "None"):
                    rgb.extend(cmd.get_color_tuple(args.pop(0)))
                else:
                    return
            elif (num.search(args[0])):
                rgb.extend([float(args.pop(0))])
            elif (abc.search(args[0])):
                if (str(cmd.get_color_tuple(args[0])) != "None"):
                    rgb.extend(cmd.get_color_tuple(args.pop(0)))
                else:
                    return
            else:
                print "Error: Unrecognized color format \""+args[0]+"\""
                return
        
        if (len(rgb) % 3):
            print "Error: Missing RGB value"
            print "Please double check RGB values"
            return
    
        dx=x2-x1
        dy=y2-y1
        dz=z2-z1
        if (len(rgb) == 3):
            rgb.extend([rgb[0]])
            rgb.extend([rgb[1]])
            rgb.extend([rgb[2]])
        t=1.0/(len(rgb)/3.0-1)
        c=len(rgb)/3-1
        s=0
        bar=[]
        
        while (s < c):
            if (len(rgb) >0):
                r=rgb.pop(0)
                g=rgb.pop(0)
                b=rgb.pop(0)
            if (s == 0 and ends == "rounded"):
                bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz, radius])
            bar.extend([CYLINDER])
            bar.extend([x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz])
            bar.extend([x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz])
            bar.extend([radius, float(r), float(g), float(b)])
            if (len(rgb) >= 3):
                bar.extend([float(rgb[0]), float(rgb[1]), float(rgb[2])])
                r=rgb[0]
                g=rgb[1]
                b=rgb[2]
            else:
                bar.extend([float(r), float(g), float(b)])
            if (s == c-1 and ends == "rounded"):
                bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz, radius])
            s=s+1
    
        cmd.delete(name)
        cmd.load_cgo(bar,name)
        
    
        return
    cmd.extend("spectrumbar",spectrumbar)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Spectrumbar&oldid=13949](https://pymolwiki.org/index.php?title=Spectrumbar&oldid=13949)"


---

## Split chains

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.editing](http://pymol.org/psicoredirect.php?psico.editing)  
  
**split_chains** is a script that creates for each chain in selection its own object. 

_This is a core command since PyMOL 1.6.0.0_

## Contents

  * 1 Usage
  * 2 Example
  * 3 The Script
  * 4 See Also



## Usage
    
    
    split_chains [ selection [, prefix ]]
    

## Example
    
    
    fetch 2yko 2ykp 2ykq, async=0
    split_chains
    delete 2yko 2ykp 2ykq
    alignto 2yko_A
    

## The Script
    
    
    from pymol import cmd, CmdException
    
    def split_chains(selection='(all)', prefix=None):
        '''
    DESCRIPTION
    
        Create a single object for each chain in selection
    
    SEE ALSO
    
        split_states, http://pymolwiki.org/index.php/Split_object
        '''
        count = 0
        models = cmd.get_object_list('(' + selection + ')')
        for model in models:
            for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
                if chain == '':
                    chain = "''"
                count += 1
                if not prefix:
                    name = '%s_%s' % (model, chain)
                else:
                    name = '%s%04d' % (prefix, count)
                cmd.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
            cmd.disable(model)
    
    cmd.extend('split_chains', split_chains)
    

## See Also

  * [split_states](/index.php/Split_States "Split States")
  * [split_object](/index.php/Split_object "Split object")



Retrieved from "[https://pymolwiki.org/index.php?title=Split_chains&oldid=11391](https://pymolwiki.org/index.php?title=Split_chains&oldid=11391)"


---

## Split Movement

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Moves two parts of one object into different directions. 

## Source Code
    
    
    fetch 1FJ1, async=0
     
    # split object
    
    
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
    

Retrieved from "[https://pymolwiki.org/index.php?title=Split_Movement&oldid=8398](https://pymolwiki.org/index.php?title=Split_Movement&oldid=8398)"


---

## Split Object Along Axis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

I have a number of small molecules that I am looking at. Many of these small molecules have distinct sections which are connected by a single bond. In a separate program I am writing, I would like to rotate each section around the given bond. I wanted a script that would allow me to select a bond, and then generate 2 selections: one for the selection of all atoms that are on one side of this bond, and the other selection for the atoms on the other side of the bond. If there is a 2nd bond which connects the 2 sections of atoms, the script will still run and terminate, but the two selections won't make much sense! 

To use, just call splitObjAlongAxis(axisSelection, outFileName) where: axisSelection is a selection that includes 2 atoms. These atoms define the axis which we are splitting across. outFileName is an optional file name where the script will spit three lines. The first line will be the two atom id's of the axis. The second line will be a list of atom id's in the selection generated, and the third line will be a list of atom id's in the second selection generated. 

## Source Code
    
    
    def findAttachedAtoms(objModel, atomIndex, excludeAtomIndex, atoms):
        atoms += [objModel.atom[atomIndex]]
        
        for b in objModel.bond:
            if((atomIndex in b.index) and (excludeAtomIndex not in b.index)):
                newAtomIndex = -1
                if(b.index[0] == atomIndex):
                    newAtomIndex = b.index[1]
                else:
                    newAtomIndex = b.index[0]
    
                newAtom = objModel.atom[newAtomIndex]
                if(newAtom not in atoms):
                    atoms += [newAtom]
                    findAttachedAtoms(objModel, newAtomIndex, excludeAtomIndex, atoms)
    
    
    def getNewSelection(objModel, atomId, excludeAtomId):
        atomIndex = -1
        excludeAtomIndex = -1
    
        for i in range(len(objModel.atom)):
            if(objModel.atom[i].id == atomId):
                atomIndex = i
            if(objModel.atom[i].id == excludeAtomId):
                excludeAtomIndex = i
    
        atoms = []
        findAttachedAtoms(objModel, atomIndex, excludeAtomIndex, atoms)
    
        return atoms
    
    def outputSelections(leftSel, rightSel, axisModel, outFileName, append):
        mode = 'w'
        if(append):
            mode = 'a'
    
        of = open(outFileName, mode)
        line = str(axisModel.atom[0].id) + ',' + str(axisModel.atom[1].id) + '\n'
        of.write(line)
    
        line = ''
        for a in leftSel:
            line += str(a.id) + ','
        of.write(line[0:-1] + '\n')
    
        line = ''
        for a in rightSel:
            line += str(a.id) + ','
        of.write(line[0:-1] + '\n')
    
        of.close()
    
    
    def splitObjAlongAxis(objName, axisSelection, outFileName = '', append = False):
        objModel = cmd.get_model(objName)
        axisModel = cmd.get_model(axisSelection)
    
        leftSel = []
        rightSel = []
    
        try:
            leftSel = getNewSelection(objModel, axisModel.atom[0].id, axisModel.atom[1].id)
            rightSel = getNewSelection(objModel, axisModel.atom[1].id, axisModel.atom[0].id)
        except:
            print 'Error in splitObjAlongAxis!'
    
        selection = ''
        for a in leftSel:
            selection += 'id ' + str(a.id) + ' or '
        selection = selection[0:-3]
        cmd.select('left',selection)
    
        selection = ''
        for a in rightSel:
            selection += 'id ' + str(a.id) + ' or '
        selection = selection[0:-3]
        cmd.select('right',selection)
    
        if(outFileName != ''):
            outputSelections(leftSel, rightSel, axisModel, outFileName, append)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Split_Object_Along_Axis&oldid=8186](https://pymolwiki.org/index.php?title=Split_Object_Along_Axis&oldid=8186)"


---

## Split object

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

split_object takes a multi-molecular object and converts it into a multi-state object. Similar to [split_states](/index.php/Split_States "Split States") but for molecules instead of states and the result of this is a multi-state object. 
    
    
    import pymol
    
    def split_object(target_obj=None, source_obj=None, max_iter=500, quiet=1, _self=cmd):
        """
    DESCRIPTION
    
        Splits a multi-molecular object into one multi-state object
    
    ARGUMENTS
    
        target_obj
    
            (string) name of target object
            
        source_obj
    
            (string) name of source object
        
        max_iter
    
            (int) maximum number of object to process; set to 0 to unlimit
        
        """
        if source_object==None:
            print "Error: Please provide a source object."
            return
    
        # ensure the user gave us one object; save for prefix
    
        obj_list = _self.get_object_list(target_obj)
    
        if len(obj_list)>1:
            print " Error: Please provide only one object at a time."
            return
    
        if target_object==None:
            target_object = _self.get_unused_name(source_obj, alwaysnumber=0)
    
        # grab unused selection name
            
        s = cmd.get_unused_name("_store")
    
        # make original selection which we'll pare down
    
        cmd.select(s, source_obj)
    
        count = 0
    
        while cmd.count_atoms(s) and count<max_iter:
            count+=1
    
            # create the object from the first molecular
            # object inside pfx
            cmd.create(pfx, "bm. first " + s, 1, count)
    
            # remove the first molecular object from
            # the source selection and re-iterate
            cmd.select(s, "%s and not bm. first %s" % (s,s))
    
        if not quiet:
            print " Created new object %s." % target_obj
    
    cmd.extend("split_object", split_object)
    

## See Also

  * [split_chains](/index.php/Split_chains "Split chains")
  * [split_states](/index.php/Split_States "Split States")



Retrieved from "[https://pymolwiki.org/index.php?title=Split_object&oldid=9424](https://pymolwiki.org/index.php?title=Split_object&oldid=9424)"


---

## Split selection

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Given an initial selection, Split_selection will create two new selections. One, called 'lo,' will have all the atoms with IDs lower than your input atom (actually input residue's alpha carbon); and the second selection is 'hi,' which contains the atoms with IDs higher than the initial (residue's alpha carbon) atom. 

  * [![The original selection](/images/a/a7/Sel_orig.png)](/index.php/File:Sel_orig.png "The original selection")

The original selection 

  * [![The lower selection](/images/f/f9/Sel_lo.png)](/index.php/File:Sel_lo.png "The lower selection")

The lower selection 

  * [![The higher selection](/images/5/56/Sel_hi.png)](/index.php/File:Sel_hi.png "The higher selection")

The higher selection 




# The Code
    
    
    import pymol
    from pymol import cmd
    
    def get_index_list(s):
    	"""
    	Given an atom selection, return the list of atom IDs in this selection
    	"""
    	return map(lambda x: x.index, cmd.get_model(s).atom)
    
    def get_single_index(s):
    	"""
    	Get the ID of the first alpha carbon in the selection, S
    	"""
    	# assume CA
    	return get_index_list( "n. CA and (br. %s)" % s)[0]
    
    def split_selection(s):
    	"""
    	PARAMS
    	    s
    	        An atom selection.
    
    	RETURNS
    	    None
    
    	SIDE EFFECTS
    	    Creates two new selections, called lo and hi.  Lo will have all atoms in the same molecule
    	    with atom IDs less than the alpha carbon ID in S.  Hi will have all the atoms in the same
    	    molecule with IDs greater than the atom ID of the alpha carbon in S.
    
    	EXAMPLE
    	    run /path/to/this/script/split_selection.py
    	    fetch 1cll
    	    select i. 79
    	    split_selection (sele)
    	    # now look at the 'hi' and 'lo' selections.
    
    	AUTHOR: Jason Vertrees, 2010.
    	"""
    	l = get_index_list("bm. " + s)
    	m = min(l)
    	M = max(l)
    	# assume using alpha carbons
    	selected_index = get_single_index( "n. CA and (br. sele)" )
    
    	low_sel_name = cmd.get_unused_name("lo")
    	hi_sel_name  = cmd.get_unused_name("hi")
    	
    	cmd.select(low_sel_name, "ID %d-%d" % (m,selected_index-1))
    	cmd.select(hi_sel_name, "ID %d-%d" % (selected_index+1,M))
    	
    cmd.extend("split_selection", split_selection)
    

# See Also

[Script_Library](/index.php?title=Script_Library&action=edit&redlink=1 "Script Library \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Split_selection&oldid=8147](https://pymolwiki.org/index.php?title=Split_selection&oldid=8147)"


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

## Stereo ray

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/stereo_ray.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/stereo_ray.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This scripts provides the stereo_ray command, which saves two ray traced images in PNG format. The left image is rotated +3 degrees and the right image is rotated -3 degrees. 

_Note that the end result is almost identical to using "walleye"[stereo](/index.php/Stereo "Stereo") mode, which might be more convenient and is readily available in PyMOL._

## Contents

  * 1 Usage
  * 2 Example
  * 3 Manually
  * 4 See Also



## Usage
    
    
    stereo_ray filename [, width [, height ]]
    

## Example
    
    
    import stereo_ray 
    stereo_ray output, 1000, 600
    

This will create two images, "output_l.png" and "output_r.png". Just paste the two images next to each other (in some image editing program) and you're done. 

[![](/images/a/ad/Stereo_ex_l.png)](/index.php/File:Stereo_ex_l.png)

[](/index.php/File:Stereo_ex_l.png "Enlarge")

Left Image

[![](/images/4/49/Stereo_ex_r.png)](/index.php/File:Stereo_ex_r.png)

[](/index.php/File:Stereo_ex_r.png "Enlarge")

Right Image

[![](/images/7/7b/Stereo_complete.png)](/index.php/File:Stereo_complete.png)

[](/index.php/File:Stereo_complete.png "Enlarge")

Combined Images

  


## Manually

To obtain the same result without the script, use the [ray](/index.php/Ray "Ray") and the [png](/index.php/Png "Png") command like this: 
    
    
    ray angle=+3
    png left-image.png
    ray angle=-3
    png right-image.png
    

Obtain almost the same result using the [stereo](/index.php/Stereo "Stereo") command: 
    
    
    stereo walleye
    png combined.png, ray=1
    

## See Also

  * [stereo](/index.php/Stereo "Stereo")
  * [stereo_mode](/index.php/Stereo_Mode "Stereo Mode")



Retrieved from "[https://pymolwiki.org/index.php?title=Stereo_ray&oldid=13931](https://pymolwiki.org/index.php?title=Stereo_ray&oldid=13931)"


---

## Stereo Figures

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The stereo command sets the current stereo mode. Stereo mode is a convenient way to "see" 3D with two images from slightly different angles. 

_There are corresponding**stereo** and [stereo_mode](/index.php/Stereo_mode "Stereo mode") settings which are controlled by the **stereo** command, so you don't need to set them directly._

## Usage
    
    
    stereo [ toggle ]
    

Valid values for the **toggle** argument are: on, swap, off, quadbuffer, crosseye, walleye, geowall, sidebyside, byrow, bycolumn, checkerboard, custom, anaglyph, dynamic, clonedynamic (see also [stereo_mode](/index.php/Stereo_mode "Stereo mode")) 

[![](/images/b/b8/Stereo_on.png)](/index.php/File:Stereo_on.png)

[](/index.php/File:Stereo_on.png "Enlarge")

Example of 1ESR shown in cross-eyed stereo

## Example
    
    
    fetch 1ESR, async=0
    as cartoon
    set cartoon_smooth_loops
    spectrum
    bg white
    stereo crosseye
    

## See Also

  * [stereo_mode](/index.php/Stereo_mode "Stereo mode")
  * [stereo_angle](/index.php/Stereo_angle "Stereo angle")
  * [stereo_shift](/index.php?title=Stereo_shift&action=edit&redlink=1 "Stereo shift \(page does not exist\)")
  * [stereo_ray](/index.php/Stereo_ray "Stereo ray")



Retrieved from "[https://pymolwiki.org/index.php?title=Stereo&oldid=10665](https://pymolwiki.org/index.php?title=Stereo&oldid=10665)"


---

## Supercell

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.xtal](http://pymol.org/psicoredirect.php?psico.xtal)  
  
[![](/images/4/48/SupercellExample.png)](/index.php/File:SupercellExample.png)

[](/index.php/File:SupercellExample.png "Enlarge")

Example with 2 unit cells in c-direction, created with: supercell 1,1,2,2x19

**supercell** can display multiple copies of the unit cell. Can also fill the unit cell (and its copies) with symmetry mates. 

See [thread on pymol-users mailing list](http://sourceforge.net/mailarchive/forum.php?thread_name=l2vdcf611bd1004140816zeca28714mf76b9f72008099ab%40mail.gmail.com&forum_name=pymol-users). 

Requires [numpy](http://numpy.scipy.org). 

## Example
    
    
    # import psico.xtal, or run the standalone script:
    run https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/xtal.py
    
    fetch 2x19, async=0
    supercell 2,1,1, 2x19, green
    supercell 1,1,2, 2x19, orange, name=super2
    

## See Also

  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [symexp](/index.php/Symexp "Symexp")



Retrieved from "[https://pymolwiki.org/index.php?title=Supercell&oldid=13313](https://pymolwiki.org/index.php?title=Supercell&oldid=13313)"


---

## SuperSym

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/SuperSymPlugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/SuperSymPlugin.py)  
Author(s)  | [Stuart Ballard](/index.php?title=User:Srballard&action=edit&redlink=1 "User:Srballard \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![](/images/3/34/SuperSymExample.png)](/index.php/File:SuperSymExample.png)

[](/index.php/File:SuperSymExample.png "Enlarge")

Symmetry partners for 1hpv showing 6-1 screw axis

[![](/images/d/da/SuperSymExample2.png)](/index.php/File:SuperSymExample2.png)

[](/index.php/File:SuperSymExample2.png "Enlarge")

Full cell of symmetry partners with symmetry axes displayed

SuperSym is a PyMOL plugin providing a large number of tools for visualization of space groups; unit cells; and symmetry axes, operators, and partners. 

The original source code is available from <https://sourceforge.net/projects/supersym/>

## Contents

  * 1 Dependencies
  * 2 Bugs
  * 3 Acknowledgments
  * 4 Installing SuperSym
  * 5 Feedback
  * 6 The Menu



## Dependencies

  * [cctbx](/index.php/CCTBX "CCTBX")
  * numpy



**PyMOL, cctbx and numpy must all be compiled with the same Python distribution!** See [CCTBX](/index.php/CCTBX "CCTBX"). 

This plugin was developed in 2010 for PyMOL version 1.2r1 and has not been updated since. 

## Bugs

Symmetry axes are not defined for all space groups, and do not display properly for some. 

## Acknowledgments

Primary coding and development was done by Stuart Ballard. All comments, questions, and issues should be directed to him at srballard@wisc.edu. 

Code for unit cell and symmetry axis building is borrowed from scripts created by Robert Campbell and Ralf W. Grosse-Kunstleve, available at <http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/>. Some of this code has been modified for use in SuperSym. 

[FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") is utilized for some of SuperSym's graphics generation, with some modifications. 

## Installing SuperSym

To install SuperSym v1.2, download SuperSymPlugin12.py from <https://sourceforge.net/projects/supersym/>. In PyMOL, go to: 

  * Plugin > Manage Plugins > Install...



A file selector dialog will appear. Select SuperSymPlugin12.py. PyMOL will direct you to restart, and upon doing so SuperSym will be accessible through the Plugin menu. 

To use functions of SuperSym directly, without creating a drop-down menu, use the run command in PyMOL on SuperSymPlugin12.py. 

Note: previous errors resulting from incorrect naming of the plugin file have been resolved in v1.2. 

## Feedback

Please post any comments, complaints, bug fix requests, useful tricks, or cool adaptations of SuperSym here. 

## The Menu

  * **Default Symmetry Partner Set**
    * See **Build Symmetry Partners > Cell [0,0,0] (default)**
  * **Draw Unit Cell**
    * Creates a cgo object with unit cell axes as cylinders. This functions similarly to _show cell_ , but the cell axes are cylinders instead of lines, allowing for printing.
  * **Build Symmetry Partners >**
    * All options in this submenu generate sets of symmetry partners
    * **Cell [0,0,0] (default)**
      * Generates a suite of symmetry partners for a given object for the default unit cell, which is lattice position [0,0,0]
    * **Cell [x,y,z] (custom)**
      * Generates a suite of symmetry partners for a given object for a lattice position which you specify
    * **2x2x2 Block**
      * Generates 8 sets of symmetry partners for a given object, filling lattice positions [0,0,0] through [1,1,1]
    * **3x3x3 Block**
      * Generates 27 sets of symmetry partners for a given object, filling lattice positions [-1,-1,-1] through [1,1,1]. This option may take a long time to execute
    * **By Partner**
      * Generates only those symmetry partners which the user specifies by their defining symmetry operators
  * **Coloring >**
    * **Default Rainbow**
      * Colors all symmetry objects with a specified by their symmetry operations automatically
    * **Select color for each operation**
      * Select symmetry partners to color by their defining symmetry operation and select the color for each
    * **Select one color for custom set of operations**
      * Select a set of symmetry partners defined by symmetry operations and select one color for all of them
  * **Graphics >**
    * **Lines**
      * Convenience function to display symmetry partners as lines
    * **Ribbon**
    * Convenience function to display symmetry partners as ribbons
    * **Cartoon**
      * Convenience function to display symmetry partners as cartoons
    * **Sphere Surface (best for printing)**
      * Uses the findSurfaceResidues function and shows surface residues as spheres. If printing, this option saves at least 60% of materials relative to regular surfaces, with minimal loss in resolution
    * **Surface (high load render)**
      * Displays symmetry partners as surfaces. This option may take a very long time to execute
  * **Symmetry Axes >**
    * **Build Axes**
      * Builds all symmetry axes for the given object. This functionality will be customizable and extended in future versions
  * **Move symmetry partners**
    * Merely displays instructions for using built in hotkeys to move symmetry partners
  * **About**
    * Developer info
  * **Help**
    * Reference to this page



Retrieved from "[https://pymolwiki.org/index.php?title=SuperSym&oldid=13315](https://pymolwiki.org/index.php?title=SuperSym&oldid=13315)"


---

## Symmetry Axis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    ###########################################################
    #
    #  Pymol script copyright Matthew O'Meara and Xavier Ambroggio 2007
    #
    #  Last updated Nov 29, 2007
    #
    #  Draw an axis given a point and a direction.  Optionally give color,
    #  length and width.
    #
    #  Usage: from the pymol menu click file->run...  then find this file.
    #  Then at the prompt, type
    #
    #           draw_axis x,y,z, i,k,j
    #
    #  where (x,y,z) is the point and (i,k,j) is the direction
    #
    #
    #  Also one can run the script as follows
    #
    #
    #           draw_axis x,y,z, i,k,j, length r,g,b, width,
    #
    #  where (r,g,b) is the color for the line in red, green, colors,
    #  width is the thickness of the line and length is the length of the
    #  line.
    #
    #
    #  For a fun example of the script run from the command line after the
    #  script is loaded
    #
    #           draw_axis_example
    #
    #
    
    from pymol.cgo import *    # get constants
    from pymol import cmd
    
    import math
    
    class Counter:
       def __init__(self):
           self.state = 1
    counter = Counter()
    
    def draw_axis(x=None, y=None, z=None, i=None, j=None, k=None, length=20.0, r=1.0, g=1.0, b=1.0, width=1.0 ):
       if x == None or y == None or z == None or i == None or j == None or k== None :
           print 'Usage: draw_axis x,y,z, i,k,j, length, r,g,b, width'
           print 'draw a line centered at (x,y,z) with the direction vector (i,j,k)'
           print 'length, color (r,g,b), and width arguments are optional'
    #        print 'For a fun example of the command, run draw_axis_example'
       else :
           x,y,z = float(x), float(y), float(z)
           i,j,k = float(i), float(j), float(k)
           r,g,b = float(r), float(g), float(b)
           width = float(width)
           length = float(length) / 2.0
    
           x1,y1,z1 = (x+i*length,y+j*length,z+k*length)
           x2,y2,z2 = (x-i*length,y-j*length,z-k*length)
    
           obj = [
               LINEWIDTH, width,
               BEGIN, LINES,
    
               COLOR,   r,  g,  b,
               VERTEX, x1, y1, z1,
               VERTEX, x2, y2, z2,
    
               END
               ]
    
           cmd.load_cgo(obj,'axis'+str(counter.state))
           counter.state += 1
    
    cmd.extend("draw_axis", draw_axis)
    
    # a simple example
    #draw_line(x=18.232,  y=17.150,  z=9.488,
    #          i=-.226639,j=0.708772,k=-.668039,
    #          r=1,       b=1,       g=1,
    #          width=1,   length=1)
    
    
    
    
    # a more complex example
    
    #import random
    #def example1(n, f):
    #    """draw a gradient field with n segments with the function f(x,y,z)=(i,j,k)"""
    #    for i in range(n):
    #        scale = 4
    #        x,y,z = [random.random()*scale for i in range(3)]
    #        i,j,k = f(x,y,z)
    
    #        draw_axis(x,y,z,i,j,k,abs(i),abs(j),abs(k))
    
    
    #def f(x,y,z):
    #    return (2*x,pow(z,2)+x,y-z)
    
    #cmd.extend("draw_axis_example", lambda :example1(1000,f))
    

Retrieved from "[https://pymolwiki.org/index.php?title=Symmetry_Axis&oldid=6392](https://pymolwiki.org/index.php?title=Symmetry_Axis&oldid=6392)"


---

## Theseus

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.fitting](http://pymol.org/psicoredirect.php?psico.fitting)  
  
theseus is a wrapper for the [theseus](http://www.theseus3d.org/) program for maximum likelihood superpositioning of macromolecular structures. It produces results very similar to [xfit](/index.php/Xfit "Xfit"). 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 See Also



## Installation

For Linux and macOS, all dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c schrodinger theseus
    

## Usage

Superpose two objects: 
    
    
    theseus mobile, target [, match [, cov [, cycles [, mobile_state [, target_state [, exe [, preserve ]]]]]]]
    

Ensemble-fit states of a multi-state object: 
    
    
    intra_theseus selection [, state [, cov [, cycles [, exe [, preserve ]]]]]
    

## Arguments

  * **mobile** = string: atom selection for mobile atoms
  * **target** = string: atom selection for target atoms
  * **match** = string: in, like, align, none or the name of an alignment object (see [local_rms](/index.php?title=Local_rms&action=edit&redlink=1 "Local rms \(page does not exist\)") help for details) {default: align}
  * **cov** = 0/1: 0 is variance weighting, 1 is covariance weighting (slower) {default: 0}
  * **cycles** = int: number of weights refinement cycles {default: 200}



For _intra_theseus_ : 

  * **selection** = string: atoms to fit
  * **state** = integer: keep transformation of this state unchanged {default: 1}



## Examples
    
    
    import psico.fitting
    fetch 1adz, async=0
    set all_states
    
    intra_theseus 1adz
    

## See Also

  * [xfit](/index.php/Xfit "Xfit")
  * [intra_xfit](/index.php/Intra_xfit "Intra xfit")
  * [intra_fit](/index.php/Intra_fit "Intra fit")
  * [intra_rms_cur](/index.php/Intra_rms_cur "Intra rms cur")
  * [align](/index.php/Align "Align")



Retrieved from "[https://pymolwiki.org/index.php?title=Theseus&oldid=12762](https://pymolwiki.org/index.php?title=Theseus&oldid=12762)"


---

## Tiff2ccp4

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will convert tiff stacks to CCP4 maps for reading in PyMOL. 

_See also the[load_img_stack](/index.php/Load_img_stack "Load img stack") script for loading image stacks directly into PyMOL_. 

If someone knows how to determine the file's dimensions from the file itself (not the image size, but the physical dimension size) please let me know. 

Note: requires Python 2.7 or later (for argparse). 

# Usage
    
    
    #
    # convert input.tiff to output.ccp4 which has cell dimensions
    # x=40 A, y=80.5 A, Z=125.2 A
    #
    python tiff2ccp4.py input.tiff output.ccp4 -x 40 -y 80.5 -z 125.2
    
    #
    # print help
    #
    python tiff2ccp4.py -h
    

# The Code
    
    
    #
    # tiff2ccp4.py -- convert a TIFF stack to a CCP4 map
    #
    # Author: Jason Vertrees
    # Date  : 2011-09-16
    #
    # Notes : To get help, type "python tiff2ccp4.py -h"
    #
    import struct
    import sys
    import string
    
    # key   = field_type
    # value = # bytes for that field
    field_types = {                                                                                                                
        1 : 1,                                                                                                                     
        2 : 1,                                                                                                                     
        3 : 2,                                                                                                                     
        4 : 4,                                                                                                                     
        5 : 8,                                                                                                                     
        6 : 1,                                                                                                                     
        7 : 1,                                                                                                                     
        8 : 2,                                                                                                                     
        9 : 4,                                                                                                                     
        10: 8,                                                                                                                     
        11: 4,                                                                                                                     
        12: 8 }                                                                                                                    
                                                                                                                                   
    field_fmt = {                                                                                                                  
        1 : "B",                                                                                                                   
        2 : "c",                                                                                                                   
        3 : "H",                                                                                                                   
        4 : "I",                                                                                                                   
        5 : "II",                                                                                                                  
        6 : "i",                                                                                                                   
        7 : "B",                                                                                                                   
        8 : "h",                                                                                                                   
        9 : "i",                                                                                                                   
        10: "ii",                                                                                                                  
        11: "f",                                                                                                                   
        12: "d" }                                                                                                                  
                                                                                                                                   
    tag_types = {                                                                                                                  
        "NewSubfileType" : 254,                                                                                                    
        "ImageWidth"  : 256,  # # cols (number of pixels per scanline)                                                             
        "ImageLength" : 257,  # # rows (number of scanlines)                                                                       
        "BitsPerSample" : 258,                                                                                                     
        "PhotometricInterpretation" : 262,                                                                                         
        "ImageDescription" : 270,                                                                                                  
        "StripOffsets" : 273,                                                                                                      
        "SamplesPerPixel" : 277,                                                                                                   
        "RowsPerStrip" : 278,                                                                                                      
        "StripByteCounts" : 279,                                                                                                   
        "XResolution" : 282,                                                                                                       
        "YResolution" : 283,                                                                                                       
        "ZResolution" : 284,                                                                                                       
        "ResolutionUnit" : 296,                                                                                                    
        }
    
    
    class TIFFStack:
    
        def __init__(self,filename):
    
            self._filename = filename
            self._f = open(self._filename, 'rb')
            
            self._file_size = self.get_file_size()
    
            # read and store the file header
    
            self._header = self.get_header()
    
            # read all IFDs and store
    
            self._IFDs = self.get_IFDs()
    
    
        def __str__(self):
            s  = ""
            s += "[TIFFStack]\n"
            s += "Filename   : %s\n" % self._filename
            s += "File size  : %s\n" % self._file_size
            s += "Header     : %s\n" % self._header
            s += "IFDs       : %s\n" % self._IFDs
            return s
    
    
        def get_file_size(self):
            """
            file size in bytes
            """
            if not self._f:
                print "Invalid: Must have open file handle."
                return -1
    
            # less typing
    
            f = self._f
    
            # get file size
    
            curP = f.tell()
            f.seek(0,2)
            sz = f.tell()
    
            # return to previous location
    
            f.seek(curP)
    
            return sz
    
    
        def get_header(self):
            return TIFFHeader(self._f)
    
    
        def get_IFDs(self):
            return IFD_Store(self._f, self._header)
    
        def get_data(self):
    
            f = self._f
    
            # bits per sample
    
            bps = 258
    
            # StripOffsets
    
            offset = 273
    
            # byte counts per strip
    
            byte_count = 279
    
            sample_format = 339
    
            data = []
    
            ifds = self._IFDs._store
    
            for curIFD in ifds:
    
                curOffset = curIFD[offset]._val
                curLength = curIFD[byte_count]._val
                curBPS = curIFD[bps]._val
                bytesPerSample = curBPS / 8.0
                
                fld = self.get_field_fmt(curIFD)
    
                # use "B" for now; support 16-bit later using tag 339 (SampleFormat)
                unpackStr = self._header._e_flg + fld * int(curLength/bytesPerSample)
    
                f.seek(curOffset)
    
                data.extend(struct.unpack(unpackStr, f.read(curLength)))
    
            return data
            
        def get_field_fmt(self,ifd):
            """
            Determines the Python struct code given
            BitsPerSample and SampleFormat from the IFD
            """
            bits_per_sample, sample_fmt = 258, 339
    
            bps = ifd[bits_per_sample]._val
    
            fmt = None
    
            if sample_fmt in ifd.keys():
                fmt = ifd[sample_fmt]._val
            else:
                fmt = 1
    
            if bps==8:
                # 1-byte unsigned int
                if fmt==1:
                    return "B"
                # 1-byte signed
                elif fmt==2:
                    return "b"
            elif bps==16:
                # 2-byte unsigned
                if fmt==1:
                    return "H"
                elif fmt==2:
                    return "h"
    
        def get_size(self):
            fX, fY = 256, 257
            x = self._IFDs._store[0][fX]._val
            y = self._IFDs._store[0][fY]._val
            z = len(self._IFDs._store)
            
            return x, y, z
    
        def get_axes_scale(self):
            # for now
            return 1,1,1
    
    
        def asCCP4(self,filename,dimX=-1,dimY=-1,dimZ=-1):
            data = self.get_data()
            
            # pixels size x,y,z
    
            nX, nY, nZ = self.get_size()
    
            # dimension scaling
    
            if -1 in (dimX,dimY,dimZ):
                dimX, dimY, dimZ = self.get_axes_scale()
                dimX *= nX
                dimY *= nY
                dimZ *= nZ
    
            m, avg, M = min(data), sum(data)/len(data), max(data)
            
            outFile = open(filename, 'wb')
            
            outFile.write(struct.pack('i',nX)) # col
            outFile.write(struct.pack('i',nY)) # row
            outFile.write(struct.pack('i',nZ)) # section
            outFile.write(struct.pack('i',2))  # mode = 2
            
            outFile.write(struct.pack('i',1))  # number of first col 
            outFile.write(struct.pack('i',1))  # '' row
            outFile.write(struct.pack('i',1))  # '' section
            
            outFile.write(struct.pack('i',nX)) # nIntervals
            outFile.write(struct.pack('i',nY)) #
            outFile.write(struct.pack('i',nZ)) #
            
            outFile.write(struct.pack('f',dimX)) # length in X
            outFile.write(struct.pack('f',dimY)) # length in Y
            outFile.write(struct.pack('f',dimZ)) # length in Z
            
            outFile.write(struct.pack('f',90.)) # alpha 
            outFile.write(struct.pack('f',90.)) # beta
            outFile.write(struct.pack('f',90.)) # gamma
            
            outFile.write(struct.pack('i',1)) # axis cols
            outFile.write(struct.pack('i',2)) # axis rows
            outFile.write(struct.pack('i',3)) # axis section
            
            outFile.write(struct.pack('f', m)) # min density
            outFile.write(struct.pack('f', avg))  # max density
            outFile.write(struct.pack('f', M))  # mean density
            
            outFile.write(struct.pack('i',0)) # space gp ?
    
            # header info; blank for us
            for x in range(24,257): outFile.write(struct.pack('i',0))
    
            # assume constant data in file
            norm = 255.
            bps = tag_types["BitsPerSample"]
            max_bits = self._IFDs._store[0][bps]._val
            norm = float(2**max_bits-1.)
    
            # read/write data
            for x in data:
                outFile.write(struct.pack('f', x/norm))
    
            outFile.close()  
    
    
    class TIFFHeader:
        def __init__(self,fileHandle):
            self._endian, self._e_flg = self.get_endian(fileHandle)
            self._magic_number = self.get_magic_number(fileHandle)
            self._first_IFD = self.get_first_IFDOffset(fileHandle)
            
        def __str__(self):
            s  = "\n"
            s += "  [TIFF Header]\n"
            s += "  Endian        : %s\n" % self._endian
            s += "  Endian Flag   : %s\n" % self._e_flg
            s += "  Magic Number  : %s\n" % self._magic_number
            s += "  IFD[0] Offset : %s" % self._first_IFD
            return s
    
        # for struct.unpackx
        def _1byte(self,n=1):
            return self._e_flg + "C"*n
        def _2byte(self,n=1):
            return self._e_flg + "H"*n
        def _4byte(self,n=1):
            return self._e_flg + "I"*n
        def _IFDbyte(self,n=1):
            return self._e_flg + "HHII"*n
    
        def get_endian(self,fileHandle):
            f = fileHandle
            f.seek(0)
            code = struct.unpack("cc", f.read(2))
            code = "".join(code)
            flg = ""
            if code=="II":
                flg = "<"
            elif code=="MM":
                flg = ">"
            else:
                print "This file is not a valid TIFF (bad endian tag)."
                flg = "?"
            return code,flg
    
        def get_magic_number(self,fileHandle):
            f = fileHandle
            f.seek(2)
            # read the magic number
            idx = 0
            if self._endian == "II": idx = 1
            _42 = struct.unpack(self._2byte(), f.read(2))[0]
            if _42!=42:
                print "Error: Improperly formatted TIFF file (bad magic number)."
                return None
            return _42
    
        def get_first_IFDOffset(self,fileHandle):
            f=fileHandle
            f.seek(4)
            off = struct.unpack(self._4byte(), f.read(4))[0]
            return off
    
    
    class IFD_Store:
        def __init__(self,fileHandle,fileHeader):
            self._f = fileHandle
            self._header = fileHeader
            self._first_offset = self._header._first_IFD
            # array of IFDs
            self._store = self.read_ifds()
    
        def __str__(self):
            s  = "\n"
            s += "  [IFD_Store]\n"
            s += "  First Offset : %d\n" % self._first_offset
            s += "  Store :\n"
            for st in self._store:
                s += "\n\n  IFD Store =>\n"
                for k in st:
                    s += "    %s => %s" % (k,st[k])
            return s
    
        def read_ifds(self):
    
            f = self._f
    
            pos = self._first_offset
    
            ifds = []
    
            while pos!=0:
    
                ifds.append({})
                
                # get number of IFD_Entries
    
                f.seek(pos)
            
                # read number of IFDs
                num_entries = struct.unpack(self._header._2byte(), f.read(2))
                num_entries = num_entries[0]
    
                # read all into IFD[x]
                for x in range(num_entries):
                    # pull the current record from file
                    curParams = struct.unpack(self._header._IFDbyte(), f.read(12))
    
                    # format the data if necessary
                    if (self._header._e_flg==">" and sys.byteorder=="little") and \
                            curParams[0] not in (270,50838,50839):
                        scale = 32 - (8 * field_types[curParams[1]])
                        scaledData = curParams[3] >> scale
                    else:
                        scaledData = curParams[3]
    
    
                    ifds[-1][curParams[0]] = IFDEntry(curParams[0], curParams[1], curParams[2], scaledData)
    
                # read next offset
                pos = struct.unpack(self._header._4byte(), f.read(4))[0]
    
            return ifds
    
    class IFDEntry:
        def __init__(self,theTag=None,theType=None,theCount=None,theValue=None):
            self._tag = theTag
            self._type = theType
            self._count = theCount
            self._val = theValue
    
        def __str__(self):
            s  = "\n"
            s += "  [IFD_Entry]\n"
            s += " Tag    : %s\n" % self._tag
            s += " Type   : %s\n" % self._type
            s += " Count  : %s\n" % self._count
            s += " Val/Off: %s\n" % self._val
            return s
    
    if __name__=="__main__":
    
        """
        Running,
    
        python tiff2ccp4.py
    
        will convert all TIFF stacks in the current directory to 
        CCP4 maps.
        """
        import argparse
        from string import split
    
        parser = argparse.ArgumentParser(description="Convert a TIFF Stack to a CCP4 Map")
        parser.add_argument("input", type=str, help="input file name (usually .tif, .tiff)")
        parser.add_argument("output", type=str, help="output file name (usually .ccp4)")
        parser.add_argument('-x',"--x", help="length of x-dimension",default=-1.,type=float)
        parser.add_argument('-y',"--y", help="length of y-dimension",default=-1.,type=float)
        parser.add_argument('-z',"--z", help="length of z-dimension",default=-1.,type=float)
        
        args = parser.parse_args()
    
        s = TIFFStack(args.input)
        s.asCCP4(args.output, args.x, args.y, args.z)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Tiff2ccp4&oldid=12428](https://pymolwiki.org/index.php?title=Tiff2ccp4&oldid=12428)"


---

## TMalign

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/tmalign.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/tmalign.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.fitting](http://pymol.org/psicoredirect.php?psico.fitting)  
  
**tmalign** is a python module (or script) that provides wrappers to TMalign, TMscore and MMalign. The executables can be downloaded from <http://zhanglab.ccmb.med.umich.edu/TM-align/> and should be saved to any directory in [PATH](http://en.wikipedia.org/wiki/PATH_\(variable\)). 

The module also provides the command **alignwithanymethod** which is useful to quickly test different alignment methods (with their respective default values). 

## Contents

  * 1 Installation
  * 2 Usage
  * 3 Examples
  * 4 See Also



## Installation

All dependencies are available from [Anaconda Cloud](https://anaconda.org): 
    
    
    conda install -c schrodinger pymol
    conda install -c schrodinger pymol-psico
    conda install -c speleo3 tmalign
    

## Usage
    
    
    tmalign mobile, target [, args [, exe [, ter [, transform [, object ]]]]]
    

**tmscore** and **mmalign** usage is equivalent. 
    
    
    alignwithanymethod mobile, target [, methods ]
    

## Examples
    
    
    fetch 2xwu 2x19 3gjx, async=0
    
    # TMscore example
    tmscore 2x19 and chain B, 2xwu and chain B
    
    # TMalign example with alignment object
    tmalign 3gjx and chain A, 2xwu and chain B, object=aln
    
    # full path to executable
    mmalign 3gjx, 2x19, exe=/usr/local/src/TM/MMalign
    

Example for **alignwithanymethod** (tmalign and [cealign](/index.php/Cealign "Cealign") will nicely align the beta sheet, but [align](/index.php/Align "Align") and [super](/index.php/Super "Super") will fail to find a nice superposition): 
    
    
    fetch 1C0M  1BCO, async=0
    remove 1C0M and not chain B
    alignwithanymethod 1C0M, 1BCO
    

  * [![1C0M_align01](/images/5/59/1C0M_align01.png)](/index.php/File:1C0M_align01.png "1C0M_align01")

1C0M_align01 

  * [![1C0M_super01](/images/4/44/1C0M_super01.png)](/index.php/File:1C0M_super01.png "1C0M_super01")

1C0M_super01 

  * [![1C0M_cealign01](/images/0/08/1C0M_cealign01.png)](/index.php/File:1C0M_cealign01.png "1C0M_cealign01")

1C0M_cealign01 

  * [![1C0M_tmalign01](/images/2/26/1C0M_tmalign01.png)](/index.php/File:1C0M_tmalign01.png "1C0M_tmalign01")

1C0M_tmalign01 




## See Also

  * [align](/index.php/Align "Align")
  * [super](/index.php/Super "Super")
  * [cealign](/index.php/Cealign "Cealign")



Retrieved from "[https://pymolwiki.org/index.php?title=TMalign&oldid=13932](https://pymolwiki.org/index.php?title=TMalign&oldid=13932)"


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

## Transform odb

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**transform_odb** transforms the coordinates of a selection and creates a new object with the transformed coordinates. The transformation matrix is read from a specified "O"-style tranformation matrix file (.odb) written by "O" or by any of the Uppsala Software Factory programs (from Gerard Klegweit) such as LSQMAN. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 USER COMMENTS
  * 4 SOURCE
  * 5 SEE ALSO



### USAGE
    
    
    transform_odb name, (selection), matrix_file [, transpose]
    

  * name = new or modified object that will contain transformed coordinates
  * selection = selection of atoms to transform
  * matrix_file = file name or path of .odb file containing a transformation matrix data block
  * transpose (default 0]



### EXAMPLES
    
    
    transform_odb moved_helix, ( mol1 and resi 200:220 ),  move_helix.odb
    

### USER COMMENTS

Please send questions or bug reports to Mark Saper, [mailto:saper@umich.edu](mailto:saper@umich.edu)

### SOURCE
    
    
    from pymol import cmd
    import pymol
    import os
    import re
    
    def __init__(self):
    	cmd.extend('transform_odb', transform_odb)
    
    # Creates a new object name from selection after transforming it with O-style matrix
    # found in matrix_file
    # Author: Mark Saper <saper@umich.edu>
    
    def transform_odb( name, selection, matrix_file='matrix.odb',  transpose=0):
    
    	# open the file for reading
    	matrix_file = os.path.expanduser(matrix_file)
    	matrix_file = os.path.expandvars(matrix_file)
    	theInFile = open ( matrix_file,"r")
    	
    	# what we'll store the results in
    	theMatrix = []
    	 
    	# read every line in the file and ...
    	for theCurrLine in theInFile.readlines():
    	   if (theCurrLine) and (theCurrLine[0] != '!') and (theCurrLine[0] != '.'):
    		  # if the line isn't blank, make a list of items seperated by tabs
    		  theNewRow = theCurrLine.split ()
    		  # add it in the matrix
    		  theMatrix.extend ( theNewRow )
    	
    	# change matrix to pymol unsupported format
    	
    	theMatrix = [ theMatrix[0], theMatrix[3], theMatrix[6], theMatrix[9],
    					  theMatrix[1], theMatrix[4], theMatrix[7], theMatrix[10],
    					  theMatrix [2], theMatrix [5], theMatrix[8], theMatrix[11], 
    					  0.0, 0.0, 0.0, 0.0 ]
    	theMatrix = [ float(x) for x in theMatrix]	
    	
    	# close the file
    	theInFile.close ()
    	
    	r = cmd.create ( name, selection)
    	r = cmd.transform_object( name, theMatrix, transpose=transpose)
    
    	return r
    
    cmd.extend('transform_odb', transform_odb)
    

### SEE ALSO

[transform_selection](/index.php/Transform_selection "Transform selection"), [transform_object](/index.php/Transform_object "Transform object")

Retrieved from "[https://pymolwiki.org/index.php?title=Transform_odb&oldid=7664](https://pymolwiki.org/index.php?title=Transform_odb&oldid=7664)"


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

## TransformSelectionByCameraView

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This little script was posted to the PyMol list. It will orient the molucule like cmd.orient() does, but does so by the camera view. 

<https://sourceforge.net/p/pymol/mailman/message/10097639/>
    
    
    # transform selection coordinates by the camera view
    #
    # The script answers this:
    #   Thanks!
    #   But translate[x,y,z] only translate the molecule.
    #   What I want  is to put longest length of molecule in the X axes, the 
    #   second Y axes, the third z axes.
    #   Just like what orient command does which change the view of camera but 
    #   not the coordinates.
    #   Now I want the coordinates also change after orient it.
    #
    cv=list(cmd.get_view())
    
    cmd.transform_selection("all", \
      cv[0:3]+[0.0]+ \
      cv[3:6]+[0.0]+ \
      cv[6:9]+[0.0]+ \
      cv[12:15]+[1.0], transpose=1)
    
    cmd.reset()
    

Retrieved from "[https://pymolwiki.org/index.php?title=TransformSelectionByCameraView&oldid=12248](https://pymolwiki.org/index.php?title=TransformSelectionByCameraView&oldid=12248)"


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

## Uniprot features

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/uniprot_features.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/uniprot_features.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
uniprot_features fetches the sequence annotation (features) from [uniprot.org](http://www.uniprot.org) and create named selections for each feature. 

Requires residue numbering (resi) to match uniprot sequence! 

## Usage
    
    
    uniprot_features uniprot_id [, selection [, withss [, prefix ]]]
    

## Example

Show the active site of [human pepsin A](http://www.uniprot.org/uniprot/PEPA_HUMAN): 
    
    
    fetch 1flh, async=0
    
    # fix residue numbering to match uniprot sequence
    alter all, resi=resv+62
    
    # create selections
    uniprot_features PEPA4_HUMAN
    
    # show sticks for active site
    as cartoon
    set cartoon_side_chain_helper
    show sticks, feature_active_site
    zoom feature_active_site
    

### See also

[select_sites](/index.php/Select_sites "Select sites")

Retrieved from "[https://pymolwiki.org/index.php?title=Uniprot_features&oldid=13934](https://pymolwiki.org/index.php?title=Uniprot_features&oldid=13934)"


---

## VisLoad

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  |   
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD](http://www.opensource.org/licenses/BSD-2-Clause)  
  
visLoad will load an object and show it in your desired representation. PyMOL by default loads things as lines (or spheres, a setting you may change), but not others. 

# Usage

  * Save the code to "visLoad.py".
  * Run the code from PyMOL or put it in your .pymolrc
  * Use visLoad whenever you would use load
  * To change representations, update "cartoon" to something else, or add more intermediate commands.



# The Code
    
    
    import os
    from os import path
    from pymol import cmd
    
    def visLoad(filename, object=None, *args, **kwargs):
            if object==None:
                    object = os.path.basename(filename).split(".")[0]
            cmd.set("suspend_updates")
            try:
                    cmd.load(filename, object, *args, **kwargs)
                    cmd.show_as("cartoon", object)
            finally:
                    cmd.set("suspend_updates", "off")
    
    cmd.extend("visLoad", visLoad)
    

# See Also

[Load](/index.php/Load "Load"), [Save](/index.php/Save "Save"), [Show_as](/index.php/Show_as "Show as")

Settings: [auto_show_lines](/index.php?title=Auto_show_lines&action=edit&redlink=1 "Auto show lines \(page does not exist\)"), [auto_show_nonbonded](/index.php?title=Auto_show_nonbonded&action=edit&redlink=1 "Auto show nonbonded \(page does not exist\)"), [auto_show_spheres](/index.php?title=Auto_show_spheres&action=edit&redlink=1 "Auto show spheres \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=VisLoad&oldid=10114](https://pymolwiki.org/index.php?title=VisLoad&oldid=10114)"


---

## Wfmesh

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/wfmesh.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/wfmesh.py)  
Author(s)  | [Dan Kulp](/index.php/User:Tmwsiy "User:Tmwsiy")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 DESCRIPTION
  * 2 IMAGES
  * 3 SETUP
  * 4 NOTES / STATUS
  * 5 EXAMPLES
  * 6 REFERENCES
  * 7 ADDITIONAL RESOURCES



### DESCRIPTION

This script will create an object for any Wavefront(.OBJ) mesh file. This is a way to extend the number of objects you can use. Also, you have more control over the coloring, transformations, etc than the CGOs. Although there are a number of these obj files on the web, you can also easily created them with open source tools ([OpenFX](http://www.openfx.org), [Crossroads3D](http://synapses.mcg.edu/tools/xroads/xroads.stm)). It takes literally, 2 min to get an object created and then loaded into pymol. Simply open OpenFX Designer, click File->Insert->Model, then choose any of the models (or create your own of course!), then export it as .3ds file. Then open the .3ds file from Crossroads3D and export as Wavefront OBJ. 

  * createWFMesh - create a mesh object from Wavefront (*.obj) formated file



### IMAGES

[![](/images/3/38/Starwars_pymol.png)](/index.php/File:Starwars_pymol.png)

[](/index.php/File:Starwars_pymol.png "Enlarge")

Star Wars Anyone?

[![](/images/8/8e/Torus_pymol.png)](/index.php/File:Torus_pymol.png)

[](/index.php/File:Torus_pymol.png "Enlarge")

A Torus, as an example shape you could do. Notice polygon normals are being used...need smoothing!

### SETUP

Simply "import wfmesh.py" 

### NOTES / STATUS

  * Tested on Pymolv0.97, Windows platform, should work on linux as well.
  * Coloring is fixed for grey and sections of mesh are stored, but not used.
  * Simple opengl calls; not optimized (display lists, etc) or anything.
  * Vertex Normal code is broken, so normals are per polygon right now.
  * Post problems in the discussion page, on 'my talk' page or just email me : dwkulp@mail.med.upenn.edu



### EXAMPLES
    
    
    import wfmesh
    cd /home/tlinnet/Software/pymol/Pymol-script-repo/files_for_examples
    wfmesh.createWFObj('torus.obj','Torus',flip=1)    # Flip = 1, if OBJ created by openFX, crossroads3D combination
    wfmesh.createWFObj("torus.obj","Torus2",translate=[5,5,0], flip=1)
    wfmesh.createWFObj("ship.obj","Ship")
    

### REFERENCES

[OpenFX](http://www.openfx.org) [Crossroads3D](http://synapses.mcg.edu/tools/xroads/xroads.stm)

### ADDITIONAL RESOURCES

Torus.obj [Torus.zip](http://pymolwiki.org/images/9/98/Torus.zip)

Retrieved from "[https://pymolwiki.org/index.php?title=Wfmesh&oldid=13935](https://pymolwiki.org/index.php?title=Wfmesh&oldid=13935)"


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

## XML-RPC server

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes a custom user implementation of an XML-RPC server. 

For PyMOL's built-in server, see [RPC](/index.php/RPC "RPC"). 

This server will enable any method in the 'cmd' module to be called remotely, as well as custom methods of the "pymol_interface" class (several examples of which are shown below). Other modules (such as util) could be wrapped using the same trick. 
    
    
    from pymol import cmd, util
    from threading import Thread
    from SimpleXMLRPCServer import SimpleXMLRPCServer
    import cPickle
    import os
    import re
    import string
    
    def pymol_startup () :
      print "Loading PyMOL XML-RPC extensions"
      server = pymol_interface()
    
    class pymol_xmlrpc_server (SimpleXMLRPCServer) :
      def __init__ (self, addr, interface) :
        self.interface = interface
        SimpleXMLRPCServer.__init__(self, addr, logRequests=0)
    
      def _dispatch (self, method, params) :
        if not self.interface.enable_xmlrpc :
          return -1
        result = -1
        func = None
        if hasattr(self.interface, method) :
          func = getattr(self.interface, method)
        elif hasattr(cmd, method) :
          func = getattr(cmd, method)
        if not callable(func) :
          print "%s is not a callable object" % method
        else :
          result = func(*params)
          if result is None :
            result = -1
        return result
    
    class pymol_interface (object) :
      def __init__ (self) :
        self.enable_xmlrpc = True
        # the port can be set via an environment variable - although it could just as easily be passed to __init__
        port = string.atoi(os.environ.get("PYMOL_XMLRPC_PORT", "9123"))
        self._server = pymol_xmlrpc_server(("localhost", port), self)
        t = threading.Thread(target=self._server.serve_forever)
        t.setDaemon(1)
        t.start()
        print "Started XML-RPC server on port %d" % port
    
      def close_all (self) :
        cmd.delete("*")
    
      def disable_all (self) :
        cmd.disable("*")
    
      def load_current_model_and_maps (self,
          pdb_file,
          fwt_map,
          delfwt_map) :
        model_name = os.path.basename(os.path.splitext(pdb_file)[0])
        cmd.load(pdb_file, model_name, state=1)
        cmd.load(fwt_map, "2fofc", state=1, format="ccp4")
        cmd.load(delfwt_map, "fofc", state=1, format="ccp4")
        cmd.isomesh("m1", "2fofc", 1.0, model_name, 5.0)
        cmd.isomesh("m2", "fofc", 3.0, model_name, 5.0)
        cmd.isomesh("m3", "fofc", -3.0, model_name, 5.0)
        cmd.color("marine", "m1")
        cmd.color("green", "m2")
        cmd.color("red", "m3")
    
      def setup_roving (self, f_map="2fofc", diff_map="fofc") :
        cmd.set("roving_detail", 20)
        cmd.set("roving_isomesh", 20)
        cmd.set("roving_origin", 1)
        cmd.set("roving_sticks", 0)
        cmd.set("roving_ribbon", 0)
        cmd.set("roving_lines", 0)
        cmd.set("roving_spheres", 0)
        cmd.set("roving_nb_spheres", 0)
        cmd.set("roving_polar_contacts", 0)
        cmd.set("roving_polar_cutoff", 0)
        cmd.set("stick_radius", 0.12)
        cmd.set("roving_map1_name", f_map)
        cmd.set("roving_map1_level", 1)
        cmd.set("roving_map2_name", diff_map)
        cmd.set("roving_map3_name", diff_map)
        cmd.set("roving_map2_level", 3)
        cmd.set("roving_map3_level", -3)
        cmd.refresh()
    
      def refresh_maps (self, f_map="2fofc", diff_map="fofc") :
        cmd.isomesh("rov_m1", f_map, 1.0, "center", 20)
        cmd.isomesh("rov_m2", diff_map, 3.0, "center", 20)
        cmd.isomesh("rov_m3", diff_map, -3.0, "center", 20)
        cmd.color("skyblue", "rov_m1")
        cmd.color("green", "rov_m2")
        cmd.color("red", "rov_m3")
    
      def show_selection (self, selection_string, zoom=True): #"True") :
        cmd.hide("sticks")
        util.cbay()
        try :
          cmd.show("sticks", selection_string)
        except Exception, e :
          return "Invalid selection."
        cmd.color("white", selection_string)
        if zoom :
          cmd.zoom(selection_string, buffer=5)
    
      def recenter (self, x, y, z) :
        view = list(cmd.get_view())
        view[12] = x
        view[13] = y
        view[14] = z
        cmd.set_view(view)
    
    if __name__ == "pymol" : # optional, if the module will be a command line argument to pymol
      pymol_startup()
    

The corresponding client code is shown below; not surprisingly, it is very simple to the example above. Because of the use of getattr() to retrieve methods in the 'cmd' API, we do not need to explicitly wrap cmd.load(), cmd.show(), or cmd.hide() for this code to work, and any other method in 'cmd' can be similarly accessed, as well as our custom methods. 
    
    
    from xmlrpclib import ServerProxy
    pymol = ServerProxy(uri="http://localhost:9123/RPC2")
    pymol.load("model.pdb")
    pymol.show("cartoon")
    pymol.hide("lines")
    pymol.load_current_model_and_maps("refined.pdb", "2mFo-DFc.ccp4", "mFo-DFc.ccp4")
    

Retrieved from "[https://pymolwiki.org/index.php?title=XML-RPC_server&oldid=12798](https://pymolwiki.org/index.php?title=XML-RPC_server&oldid=12798)"


---

## Zero residues

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 OVERVIEW
    * 1.1 Ordering from Sequence End
  * 2 EXAMPLE
  * 3 INSTALL
  * 4 CODE



## OVERVIEW

This script will renumber all the residues such that the first one is numbered 0. This is often helpful when dealing with alignments. 

### Ordering from Sequence End

If you want to change numbering based off the last residue's number in the sequence, just replace **first** in the code with **last**. 

## EXAMPLE
    
    
    zero_residues 1AFK
    zero_residues *
    
    # make the first residue's number, 30.
    zero_residues 1AFK, 30
    

## INSTALL

Copy the source code below, to "zero_residues.py" and then simply run the file. The command, "zero_residues" will now be defined and can be used as in the examples above. 

## CODE
    
    
    from pymol import cmd, stored
    
    def zero_residues(sel1,offset=0,chains=0):
            """
    DESCRIPTION
    
        Renumbers the residues so that the first one is zero, or offset
    
    USAGE
    
        zero_residues selection [, offset [, chains ]]
    
    EXAMPLES
    
        zero_residues protName            # first residue is 0
        zero_residues protName, 5         # first residue is 5
        zero_residues protName, chains=1  # each chain starts at 0
        zero_residues *
            """
            offset = int(offset)
    
            # variable to store the offset
            stored.first = None
            # get the names of the proteins in the selection
    
            names = ['(model %s and (%s))' % (p, sel1)
                            for p in cmd.get_object_list('(' + sel1 + ')')]
    
            if int (chains):
                    names = ['(%s and chain %s)' % (p, chain)
                                    for p in names
                                    for chain in cmd.get_chains(p)]
    
            # for each name shown
            for p in names:
                    # get this offset
                    ok = cmd.iterate("first %s and polymer and n. CA" % p,"stored.first=resv")
                    # don't waste time if we don't have to
                    if not ok or stored.first == offset:
                            continue;
                    # reassign the residue numbers
                    cmd.alter("%s" % p, "resi=str(int(resi)-%s)" % str(int(stored.first)-offset))
                    # update pymol
    
            cmd.rebuild()
    
    # let pymol know about the function
    cmd.extend("zero_residues", zero_residues)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Zero_residues&oldid=10939](https://pymolwiki.org/index.php?title=Zero_residues&oldid=10939)"


---

## Category:Biochemical Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Pages in category "Biochemical Scripts"

The following 26 pages are in this category, out of 26 total. 

### A

  * [AAindex](/index.php/AAindex "AAindex")
  * [Average b](/index.php/Average_b "Average b")



### C

  * [Color cbcobj](/index.php/Color_cbcobj "Color cbcobj")
  * [Color h](/index.php/Color_h "Color h")
  * [CreateAtom](/index.php/CreateAtom "CreateAtom")
  * [Cyspka](/index.php/Cyspka "Cyspka")



### D

  * [DistancesRH](/index.php/DistancesRH "DistancesRH")
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")



### F

  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [Forster distance calculator](/index.php/Forster_distance_calculator "Forster distance calculator")



### L

  * [ListSelection2](/index.php/ListSelection2 "ListSelection2")
  * [Load new B-factors](/index.php/Load_new_B-factors "Load new B-factors")



### P

  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")
  * [Polarpairs](/index.php/Polarpairs "Polarpairs")
  * [Propka](/index.php/Propka "Propka")
  * [Pucker](/index.php/Pucker "Pucker")
  * [Pytms](/index.php/Pytms "Pytms")



### Q

  * [Quickdisplays](/index.php/Quickdisplays "Quickdisplays")



### R

  * [Resicolor](/index.php/Resicolor "Resicolor")



### S

  * [Show aromatics](/index.php/Show_aromatics "Show aromatics")
  * [Show charged](/index.php/Show_charged "Show charged")
  * [Show contacts](/index.php/Show_contacts "Show contacts")
  * [Show hydrophilic](/index.php/Show_hydrophilic "Show hydrophilic")
  * [Show hydrophobics](/index.php/Show_hydrophobics "Show hydrophobics")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Biochemical_Scripts&oldid=6380](https://pymolwiki.org/index.php?title=Category:Biochemical_Scripts&oldid=6380)"


---

## Category:Math Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  * [Perp Maker](/index.php/Perp_Maker "Perp Maker") \-- Creates a perpendicular plane through the center of your protein with respect to the camera's current position. (If you translate the protein towards the camera a bit, you get a nice surface, sometimes.) A stupid little script I wrote in response to a request a few months ago (and it doesn't even conform to the request!) Load a protein, run the script (read the documentation in the script). (Jason Vertrees/[Tree](/index.php/User:Tree "User:Tree"))
  * [Axes](/index.php/Axes "Axes") \-- Creates a 3D-CGO object that shows the three coordinate axes.
  * [Symmetry Axis](/index.php/Symmetry_Axis "Symmetry Axis") \-- Draw a 3D-CGO line given a point and a direction.
  * [CGO Text](/index.php/CGO_Text "CGO Text") \-- Creates a 3D-CGO text object.
  * [Slerpy](/index.php/Slerpy "Slerpy") \-- Pymol command extensions for key frame animation movie making.
  * [Plane Wizard](/index.php/Plane_Wizard "Plane Wizard") \-- Wizard to draw planes between three picked points.
  * [Bounding Box](/index.php/Bounding_Box "Bounding Box") \-- Create a bounding box around a selection (Python script; requires numarray and Scientific; gilleain)
  * [Ellipsoid](/index.php/Ellipsoid "Ellipsoid") \-- Create callback object (opengl) ellipsoids. (Python script; gilleain)
  * [TransformSelectionByCameraView](/index.php/TransformSelectionByCameraView "TransformSelectionByCameraView") \-- Transforms the selection by the camera view.
  * [modevectors](/index.php/Modevectors "Modevectors") \-- A wonderful script that allows you to draw highly customizable vectors between two states (NMA, TMD, etc)
  * [Center Of Mass](/index.php/Center_Of_Mass "Center Of Mass") \-- Given a selection of atoms (of equal weight) - Calculates the center of mass and represents it with a CGO sphere



## Pages in category "Math Scripts"

The following 29 pages are in this category, out of 29 total. 

### A

  * [Axes](/index.php/Axes "Axes")



### B

  * [BbPlane](/index.php/BbPlane "BbPlane")
  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat")
  * [Bounding Box](/index.php/Bounding_Box "Bounding Box")



### C

  * [Cart to frac](/index.php/Cart_to_frac "Cart to frac")
  * [Center of mass](/index.php/Center_of_mass "Center of mass")
  * [Cgo arrow](/index.php/Cgo_arrow "Cgo arrow")
  * [Cgo grid](/index.php/Cgo_grid "Cgo grid")
  * [CGO Text](/index.php/CGO_Text "CGO Text")
  * [CgoCircle](/index.php/CgoCircle "CgoCircle")
  * [Contact Surface](/index.php/Contact_Surface "Contact Surface")
  * [Cubes](/index.php/Cubes "Cubes")



### D

  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")
  * [Draw Protein Dimensions](/index.php/Draw_Protein_Dimensions "Draw Protein Dimensions")
  * [DrawBoundingBox](/index.php/DrawBoundingBox "DrawBoundingBox")
  * [Dump2CGO](/index.php/Dump2CGO "Dump2CGO")



### E

  * [Ellipsoid](/index.php/Ellipsoid "Ellipsoid")



### M

  * [Mark center](/index.php/Mark_center "Mark center")
  * [Modevectors](/index.php/Modevectors "Modevectors")



### P

  * [Perp maker](/index.php/Perp_maker "Perp maker")
  * [Plane Wizard](/index.php/Plane_Wizard "Plane Wizard")



### R

  * [Radius of gyration](/index.php/Radius_of_gyration "Radius of gyration")
  * [RotationAxis](/index.php/RotationAxis "RotationAxis")



### S

  * [Slerpy](/index.php/Slerpy "Slerpy")
  * [Supercell](/index.php/Supercell "Supercell")
  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [Symmetry Axis](/index.php/Symmetry_Axis "Symmetry Axis")



### T

  * [TransformSelectionByCameraView](/index.php/TransformSelectionByCameraView "TransformSelectionByCameraView")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Math_Scripts&oldid=6390](https://pymolwiki.org/index.php?title=Category:Math_Scripts&oldid=6390)"


---

## Category:ObjSel Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  * [ConnectedCloud](/index.php/ConnectedCloud "ConnectedCloud")



    

    Find connected clouds of objects in PyMOL

  * [Color Objects](/index.php/Color_Objects "Color Objects")



    

    Colors all objects differently (Python script).

  * [FindSeq](/index.php/FindSeq "FindSeq")



    

    Find sequence or regular expression of amino acids in a protein.

  * [Get Coordinates I](/index.php/Get_Coordinates_I "Get Coordinates I")



    

    Retrieves atom coordinates as Python objects.

  * [Get Coordinates II](/index.php/Get_Coordinates_II "Get Coordinates II")



    

    Retrieves atom coordinates as Python array (list object).

  * [grepsel](/index.php/Grepsel "Grepsel")



    

    Make named selections using regular expressions (protein sequence).

  * [List Selection](/index.php/List_Selection "List Selection")



    

    Prints a list of all residues in a selection (both Python and .pml).

  * [List Colors](/index.php/List_Colors "List Colors")



    

    Lists the color of all residues in a selection (both Python and .pml).

  * [List Secondary Structures](/index.php/List_Secondary_Structures "List Secondary Structures")



    

    Secondary structures (both predefined and those calculated with the 'dss' command) can be exported as a long string ('HHHHLLLLSSS').

  * [Selection Exists](/index.php/Selection_Exists "Selection Exists")



    

    Python method that returns true if a selection of a given name exists.

  * [Split Movement](/index.php/Split_Movement "Split Movement")



    

    Moves two parts of one object into different directions.

  * [Split Object Along Axis](/index.php/Split_Object_Along_Axis "Split Object Along Axis")



    

    Splits an object into 2 pieces along a given axis.

  * [toGroup](/index.php/ToGroup "ToGroup")



    

    Convert a multistate object into a group of single state objects.

  * [Zero_residues](/index.php/Zero_residues "Zero residues")



    

    Renumber residues such that the first residue is 0. Useful for alignments.

## Pages in category "ObjSel Scripts"

The following 39 pages are in this category, out of 39 total. 

### A

  * [AlphaToAll](/index.php/AlphaToAll "AlphaToAll")
  * [Annotate v](/index.php/Annotate_v "Annotate v")



### C

  * [Cluster Count](/index.php/Cluster_Count "Cluster Count")
  * [CollapseSel](/index.php/CollapseSel "CollapseSel")
  * [Color cbcobj](/index.php/Color_cbcobj "Color cbcobj")
  * [Color Objects](/index.php/Color_Objects "Color Objects")
  * [ConnectedCloud](/index.php/ConnectedCloud "ConnectedCloud")
  * [Count molecules in selection](/index.php/Count_molecules_in_selection "Count molecules in selection")



### D

  * [DistancesRH](/index.php/DistancesRH "DistancesRH")



### E

  * [Expand To Surface](/index.php/Expand_To_Surface "Expand To Surface")



### F

  * [Find buried waters](/index.php/Find_buried_waters "Find buried waters")
  * [FindObjectsNearby](/index.php/FindObjectsNearby "FindObjectsNearby")
  * [Findseq](/index.php/Findseq "Findseq")
  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [Flatten obj](/index.php/Flatten_obj "Flatten obj")



### G

  * [Get Coordinates I](/index.php/Get_Coordinates_I "Get Coordinates I")
  * [Get Coordinates II](/index.php/Get_Coordinates_II "Get Coordinates II")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")
  * [GetNamesInSel](/index.php/GetNamesInSel "GetNamesInSel")
  * [Grepsel](/index.php/Grepsel "Grepsel")



### L

  * [List Colors](/index.php/List_Colors "List Colors")
  * [List Secondary Structures](/index.php/List_Secondary_Structures "List Secondary Structures")
  * [List Selection](/index.php/List_Selection "List Selection")
  * [ListSelection2](/index.php/ListSelection2 "ListSelection2")



### P

  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")



### R

  * [Removealt](/index.php/Removealt "Removealt")
  * [Rotkit](/index.php/Rotkit "Rotkit")



### S

  * [Save sep](/index.php/Save_sep "Save sep")
  * [SaveGroup](/index.php/SaveGroup "SaveGroup")
  * [SelectClipped](/index.php/SelectClipped "SelectClipped")
  * [Selection Exists](/index.php/Selection_Exists "Selection Exists")
  * [SelInside](/index.php/SelInside "SelInside")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")
  * [Split Movement](/index.php/Split_Movement "Split Movement")
  * [Split Object Along Axis](/index.php/Split_Object_Along_Axis "Split Object Along Axis")
  * [Split selection](/index.php/Split_selection "Split selection")



### T

  * [ToGroup](/index.php/ToGroup "ToGroup")



### Z

  * [Zero residues](/index.php/Zero_residues "Zero residues")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:ObjSel_Scripts&oldid=8152](https://pymolwiki.org/index.php?title=Category:ObjSel_Scripts&oldid=8152)"


---

## Category:Structural Biology Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page is a sub-category of scripts for Structural Biology Applications. 

  * [Cealign](/index.php/Cealign "Cealign") \-- Implementation of the CE Structure Alignment algorithm as a PyMOL plugin.
  * [Kabsch](/index.php/Kabsch "Kabsch") \-- Kabsch alignment of two sets of vectors. (Part 2 of a protein alignment.)
  * [LigAlign](/index.php/LigAlign "LigAlign") \-- Ligand-based active site alignment and comparison.
  * [WriteSS](/index.php/WriteSS "WriteSS") \-- Writes secondary structural elements, for each residues, to a file.
  * [ss](/index.php/Ss "Ss") \-- Simple command to summarise the Secondary Structure as a list of "start-end type" like sses.
  * [iterate_sses](/index.php/Iterate_sses "Iterate sses") \-- Slightly more complex version of "ss" that allows the user to pass in a function to act on the sse list.
  * [Helicity_check](/index.php/Helicity_check "Helicity check") \-- helicity_check show the evolution of O - N distances over an amino acid sequence
  * [Measure Distance](/index.php/Measure_Distance "Measure Distance") \-- Measures the distance between two atoms (Python script).
  * [Translate_And_Measure](/index.php/Translate_And_Measure "Translate And Measure") \-- prints **overlap** if any of the atoms in molA or molB were within 4 Angstrom after translating by 1 along X
  * [motif](/index.php/Motif "Motif") \-- Designed for easy display of backbone motifs (nests, catgrips, etc).
  * [Show NMR constrains](/index.php/Show_NMR_constrains "Show NMR constrains") \-- This script will display the NMR constrains used for a structure calculation atop a structure. Usage: Save this as "NMRcnstr.py" load your protein in PyMOL, and run the script. type upl('fname') or cns('fname') where fname is the filename with the NMR constrains you want to display.
  * [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure") \-- Grow a helix,strand or loop from ends of proteins
  * [DynoPlot](/index.php/DynoPlot "DynoPlot") \-- Generic dynamic plotting utility; Interactive Ramachandran Plots.
  * [Rotamer Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle") \-- Toggle through the most common side chain rotamers and/or color by rotamer probability (Dunbrack's Backbone-dependent library)
  * [DisplacementMap](/index.php/DisplacementMap "DisplacementMap") \-- Finds the best residue pair for FRET and EPR measurements. Given and open and closed form of a protein (and after pymol alignment), it returns



a data-matrix for displacement between residues and a gnuplot plot file. It returns best positive and negative delta displacement suggestions. The suggestions can be mutated into cysteines and labelled for FRET/EPR measurements. 

## Pages in category "Structural Biology Scripts"

The following 61 pages are in this category, out of 61 total. 

### A

  * [AAindex](/index.php/AAindex "AAindex")
  * [Angle between domains](/index.php/Angle_between_domains "Angle between domains")
  * [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices")
  * [AutoMultiFit](/index.php/AutoMultiFit "AutoMultiFit")
  * [Average b](/index.php/Average_b "Average b")



### B

  * [BbPlane](/index.php/BbPlane "BbPlane")
  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat")
  * [Bondpack](/index.php/Bondpack "Bondpack")



### C

  * [Cart to frac](/index.php/Cart_to_frac "Cart to frac")
  * [Ccp4 contact](/index.php/Ccp4_contact "Ccp4 contact")
  * [Ccp4 ncont](/index.php/Ccp4_ncont "Ccp4 ncont")
  * [Ccp4 pisa](/index.php/Ccp4_pisa "Ccp4 pisa")
  * [Centroid](/index.php/Centroid "Centroid")
  * [Color by conservation](/index.php/Color_by_conservation "Color by conservation")
  * [Color By Mutations](/index.php/Color_By_Mutations "Color By Mutations")
  * [Colorblindfriendly](/index.php/Colorblindfriendly "Colorblindfriendly")
  * [Colorbydisplacement](/index.php/Colorbydisplacement "Colorbydisplacement")
  * [ColorByRMSD](/index.php/ColorByRMSD "ColorByRMSD")
  * [Consistent View/ Map Inspect](/index.php/Consistent_View/_Map_Inspect "Consistent View/ Map Inspect")
  * [Contact Surface](/index.php/Contact_Surface "Contact Surface")
  * [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure")
  * [Cyspka](/index.php/Cyspka "Cyspka")



### D

  * [Displacementmap](/index.php/Displacementmap "Displacementmap")
  * [DistancesRH](/index.php/DistancesRH "DistancesRH")
  * [Draw Protein Dimensions](/index.php/Draw_Protein_Dimensions "Draw Protein Dimensions")
  * [Dssp](/index.php/Dssp "Dssp")
  * [Dssr block](/index.php/Dssr_block "Dssr block")
  * [DynoPlot](/index.php/DynoPlot "DynoPlot")



### E

  * [Elbow angle](/index.php/Elbow_angle "Elbow angle")



### F

  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")



### H

  * [Helicity check](/index.php/Helicity_check "Helicity check")
  * [HighlightAlignedSS](/index.php/HighlightAlignedSS "HighlightAlignedSS")



### I

  * [InterfaceResidues](/index.php/InterfaceResidues "InterfaceResidues")
  * [Intra xfit](/index.php/Intra_xfit "Intra xfit")
  * [Iterate sses](/index.php/Iterate_sses "Iterate sses")



### L

  * [LigAlign](/index.php/LigAlign "LigAlign")
  * [Load aln](/index.php/Load_aln "Load aln")
  * [Load new B-factors](/index.php/Load_new_B-factors "Load new B-factors")



### M

  * [Mcsalign](/index.php/Mcsalign "Mcsalign")
  * [Measure Distance](/index.php/Measure_Distance "Measure Distance")
  * [Motif](/index.php/Motif "Motif")
  * [Msms surface](/index.php/Msms_surface "Msms surface")



### P

  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")
  * [PDB plugin](/index.php/PDB_plugin "PDB plugin")
  * [Propka](/index.php/Propka "Propka")



### R

  * [Resicolor plugin](/index.php/Resicolor_plugin "Resicolor plugin")
  * [RmsdByResidue](/index.php/RmsdByResidue "RmsdByResidue")
  * [Rotamer Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle")
  * [RotationAxis](/index.php/RotationAxis "RotationAxis")



### S

  * [Show contacts](/index.php/Show_contacts "Show contacts")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")
  * [Sidechaincenters](/index.php/Sidechaincenters "Sidechaincenters")
  * [Ss](/index.php/Ss "Ss")
  * [Sst](/index.php/Sst "Sst")



### N

  * [Nmr cnstr](/index.php/Nmr_cnstr "Nmr cnstr")



### T

  * [Transformations](/index.php/Transformations "Transformations")
  * [Translate And Measure](/index.php/Translate_And_Measure "Translate And Measure")



### W

  * [WriteSS](/index.php/WriteSS "WriteSS")



### X

  * [Xfit](/index.php/Xfit "Xfit")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Structural_Biology_Scripts&oldid=8930](https://pymolwiki.org/index.php?title=Category:Structural_Biology_Scripts&oldid=8930)"


---

## Category:System Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  * [FilterByMol](/index.php/FilterByMol "FilterByMol") \-- Filter a directory of files, and save their ligands to disk (by molecule).
  * [LoadDir](/index.php/LoadDir "LoadDir") \-- loads all the files of a type you specify from the path you specify and puts them into the group you specify (or none).
  * [Process_All_Files_In_Directory](/index.php/Process_All_Files_In_Directory "Process All Files In Directory") \-- Do something to all files in a directory. The examples show how to print the disulfide bond lengths, then in general all sulfur distances (not necessarily bound).
  * [PythonTerminal](/index.php/PythonTerminal "PythonTerminal") \-- Allows execution of python commands from the PyMOL command line.
  * [pdbsurvey](/index.php/Pdbsurvey "Pdbsurvey") \-- Surveys the pdb for recently added structures that are relevant to a user-specified keywords list (in a text file)
  * [Read PDB-String](/index.php/Read_PDB-String "Read PDB-String") \-- Parses a string in PDB format to a PyMOL object.
  * [Monitor file continuously](/index.php/Monitor_file_continuously "Monitor file continuously") \-- Watch a file on a separate thread and re-load when it is modified.
  * [XML-RPC server](/index.php/XML-RPC_server "XML-RPC server") \-- An API for controlling PyMOL remotely (from the same computer or on the network). Adapted from the server built in to PyMOL.



## Pages in category "System Scripts"

The following 12 pages are in this category, out of 12 total. 

### F

  * [FetchLocal](/index.php/FetchLocal "FetchLocal")
  * [FilterByMol](/index.php/FilterByMol "FilterByMol")
  * [Find symbol](/index.php/Find_symbol "Find symbol")



### L

  * [LoadDir](/index.php/LoadDir "LoadDir")



### M

  * [Monitor file continuously](/index.php/Monitor_file_continuously "Monitor file continuously")



### P

  * [Pdbsurvey](/index.php/Pdbsurvey "Pdbsurvey")
  * [Pml2py](/index.php/Pml2py "Pml2py")
  * [Process All Files In Directory](/index.php/Process_All_Files_In_Directory "Process All Files In Directory")
  * [PythonTerminal](/index.php/PythonTerminal "PythonTerminal")



### R

  * [Read PDB-String](/index.php/Read_PDB-String "Read PDB-String")



### T

  * [Tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")



### X

  * [XML-RPC server](/index.php/XML-RPC_server "XML-RPC server")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:System_Scripts&oldid=10467](https://pymolwiki.org/index.php?title=Category:System_Scripts&oldid=10467)"


---

## Category:ThirdParty Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Pages in category "ThirdParty Scripts"

The following 12 pages are in this category, out of 12 total. 

### C

  * [Ccp4 contact](/index.php/Ccp4_contact "Ccp4 contact")
  * [Ccp4 ncont](/index.php/Ccp4_ncont "Ccp4 ncont")
  * [Ccp4 pisa](/index.php/Ccp4_pisa "Ccp4 pisa")



### P

  * [PoseView](/index.php/PoseView "PoseView")
  * [PovRay](/index.php/PovRay "PovRay")
  * [Pymol2glmol](/index.php/Pymol2glmol "Pymol2glmol")



### R

  * [Rasmolify](/index.php/Rasmolify "Rasmolify")
  * [RUCAP UM-5](/index.php/RUCAP_UM-5 "RUCAP UM-5")



### S

  * [Save Mopac](/index.php/Save_Mopac "Save Mopac")



### T

  * [Tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")
  * [Transform odb](/index.php/Transform_odb "Transform odb")



### W

  * [Wfmesh](/index.php/Wfmesh "Wfmesh")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:ThirdParty_Scripts&oldid=6376](https://pymolwiki.org/index.php?title=Category:ThirdParty_Scripts&oldid=6376)"


---

## Category:UI Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  * [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz") \-- A script used in conjunction with head tracking software to provide an immersive virtual experience.
  * [ObjectByArrows](/index.php/ObjectByArrows "ObjectByArrows") \-- Navigate between objects by using the arrow keys
  * [PowerMate Dial OS X](/index.php/PowerMate_Dial_OS_X "PowerMate Dial OS X") \-- Script and instructions to use the PowerMate dial on Mac OS X.
  * [mouse_modes](/index.php/Mouse_modes "Mouse modes") \-- customize the default mouse bindings for Viewing or Editing modes. - _by EHP_
  * [Key Wait](/index.php/Key_Wait "Key Wait") \-- Process key events in a Python script.
  * [grepset](/index.php/Grepset "Grepset") \-- List all settings matching a given keyword. - _by EHP_
  * [apropos](/index.php/Apropos "Apropos") \-- List all commands matching a given keyword or whose docs contain the keyword. - _by EHP_
  * [Stereo_Ray](/index.php/Stereo_Ray "Stereo Ray") \-- This script will create two resolution specific ray traced images rotated appropriately for inclusion into a single file to represent a stereo view of the desired macromolecule.
  * [Save2traj](/index.php/Save2traj "Save2traj") \-- Save object and object states into a trajectory file



## Pages in category "UI Scripts"

The following 30 pages are in this category, out of 30 total. 

### A

  * [Aa codes](/index.php/Aa_codes "Aa codes")
  * [Apropos](/index.php/Apropos "Apropos")



### D

  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")
  * [Dynamic mesh](/index.php/Dynamic_mesh "Dynamic mesh")



### F

  * [Format bonds](/index.php/Format_bonds "Format bonds")
  * [Frame slider](/index.php/Frame_slider "Frame slider")



### G

  * [Get colors](/index.php/Get_colors "Get colors")
  * [Grepset](/index.php/Grepset "Grepset")



### I

  * [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz")
  * [Inertia tensor](/index.php/Inertia_tensor "Inertia tensor")



### K

  * [Key Wait](/index.php/Key_Wait "Key Wait")



### M

  * [Make Figures](/index.php/Make_Figures "Make Figures")
  * [Mouse modes](/index.php/Mouse_modes "Mouse modes")
  * [Movie color fade](/index.php/Movie_color_fade "Movie color fade")
  * [Movie fade](/index.php/Movie_fade "Movie fade")
  * [Movit](/index.php/Movit "Movit")



### O

  * [ObjectByArrows](/index.php/ObjectByArrows "ObjectByArrows")
  * [ObjectFocus](/index.php/ObjectFocus "ObjectFocus")



### P

  * [PDB Web Services Script](/index.php/PDB_Web_Services_Script "PDB Web Services Script")
  * [PowerMate Dial OS X](/index.php/PowerMate_Dial_OS_X "PowerMate Dial OS X")
  * [Pytms](/index.php/Pytms "Pytms")



### Q

  * [Quickdisplays](/index.php/Quickdisplays "Quickdisplays")



### R

  * [RUCAP UM-5](/index.php/RUCAP_UM-5 "RUCAP UM-5")



### S

  * [Save settings](/index.php/Save_settings "Save settings")
  * [Save2traj](/index.php/Save2traj "Save2traj")
  * [Set toggle](/index.php/Set_toggle "Set toggle")
  * [Spectrumbar](/index.php/Spectrumbar "Spectrumbar")
  * [Stereo Figures](/index.php/Stereo_Figures "Stereo Figures")
  * [Stereo ray](/index.php/Stereo_ray "Stereo ray")



### V

  * [VisLoad](/index.php/VisLoad "VisLoad")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:UI_Scripts&oldid=7362](https://pymolwiki.org/index.php?title=Category:UI_Scripts&oldid=7362)"


---

