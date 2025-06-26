# Category: ObjSel Scripts

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

