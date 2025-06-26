# Category: Biochemical Scripts

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

## Polarpairs

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Find polar pairs with the [cmd.find_pairs](/index.php/Find_pairs "Find pairs") function and return them as a list of atom pairs. It finds (more or less) the same contacts like [cmd.distance(mode=2)](/index.php/Distance "Distance"). 

## Example
    
    
    pairs = polarpairs('chain A', 'chain B')
    for p in pairs:
        dist = cmd.get_distance('(%s`%s)' % p[0], '(%s`%s)' % p[1])
        print(p, 'Distance: %.2f' % (dist))
    

## The Script
    
    
    '''
    (c) 2011 Thomas Holder, MPI for Developmental Biology
    '''
    
    from pymol import cmd
    
    def polarpairs(sel1, sel2, cutoff=4.0, angle=63.0, name='', state=1, quiet=1):
        '''
    ARGUMENTS
    
        sel1, sel2 = string: atom selections
    
        cutoff = float: distance cutoff
    
        angle = float: h-bond angle cutoff in degrees. If angle="default", take
        "h_bond_max_angle" setting. If angle=0, do not detect h-bonding.
    
        name = string: If given, also create a distance object for visual representation
    
    SEE ALSO
    
        cmd.find_pairs, cmd.distance
        '''
        cutoff = float(cutoff)
        quiet = int(quiet)
        state = int(state)
        if angle == 'default':
            angle = cmd.get('h_bond_max_angle', cmd.get_object_list(sel1)[0])
        angle = float(angle)
        mode = 1 if angle > 0 else 0
        x = cmd.find_pairs('(%s) and donors' % sel1, '(%s) and acceptors' % sel2,
                state, state,
                cutoff=cutoff, mode=mode, angle=angle) + \
            cmd.find_pairs('(%s) and acceptors' % sel1, '(%s) and donors' % sel2,
                state, state,
                cutoff=cutoff, mode=mode, angle=angle)
        x = sorted(set(x))
        if not quiet:
            print('Settings: cutoff=%.1fangstrom angle=%.1fdegree' % (cutoff, angle))
            print('Found %d polar contacts' % (len(x)))
        if len(name) > 0:
            for p in x:
                cmd.distance(name, '(%s`%s)' % p[0], '(%s`%s)' % p[1])
        return x
    
    cmd.extend('polarpairs', polarpairs)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Polarpairs&oldid=12897](https://pymolwiki.org/index.php?title=Polarpairs&oldid=12897)"


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

