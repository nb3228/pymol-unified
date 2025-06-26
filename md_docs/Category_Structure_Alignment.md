# Category: Structure Alignment

## Align

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/6/6e/After_alignment.png)](/index.php/File:After_alignment.png)

[](/index.php/File:After_alignment.png "Enlarge")

Two proteins after structure alignment

**align** performs a sequence alignment followed by a structural superposition, and then carries out zero or more cycles of refinement in order to reject structural outliers found during the fit. align does a good job on proteins with decent sequence similarity (identity >30%). For comparing proteins with lower sequence identity, the [super](/index.php/Super "Super") and [cealign](/index.php/Cealign "Cealign") commands perform better. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Alignment Objects
  * 4 RMSD
  * 5 Examples
  * 6 PyMOL API
  * 7 Notes
  * 8 See Also



## Usage
    
    
    align mobile, target [, cutoff [, cycles
        [, gap [, extend [, max_gap [, object
        [, matrix [, mobile_state [, target_state
        [, quiet [, max_skip [, transform [, reset ]]]]]]]]]]]]]
    

## Arguments

  * **mobile** = string: atom selection of mobile object
  * **target** = string: atom selection of target object
  * **cutoff** = float: outlier rejection cutoff in RMS {default: 2.0}
  * **cycles** = int: maximum number of outlier rejection cycles {default: 5}
  * **gap, extend, max_gap** : sequence alignment parameters
  * **object** = string: name of alignment object to create {default: (no alignment object)}
  * **matrix** = string: file name of substitution matrix for sequence alignment {default: BLOSUM62}
  * **mobile_state** = int: object state of mobile selection {default: 0 = all states}
  * **target_state** = int: object state of target selection {default: 0 = all states}
  * **quiet** = 0/1: suppress output {default: 0 in command mode, 1 in API}
  * **max_skip** = ?
  * **transform** = 0/1: do superposition {default: 1}
  * **reset** = ?



## Alignment Objects

An alignment object can be created with the **object=**_somename_ argument. An alignment object provides: 

  * aligned [sequence viewer](/index.php/Seq_view "Seq view")
  * graphical representation of aligned atom pairs as lines in the 3D viewer
  * can be [saved](/index.php/Save "Save") to a clustalw sequence alignment file



## RMSD

The RMSD of the aligned atoms (after outlier rejection!) is reported in the text output. The **all-atom RMSD** can be obtained by setting **cycles=0** and thus not doing any outlier rejection. The RMSD can also be captured with a python script, see the API paragraph below. Note that the output prints "RMS" but it is in fact "RMSD" and the units are Angstroms. 

## Examples
    
    
    fetch 1oky 1t46, async=0
    
    # 1) default with outlier rejection
    align 1oky, 1t46
    
    # 2) with alignment object, save to clustalw file
    align 1oky, 1t46, object=alnobj
    save alignment.aln, alnobj
    
    # 3) all-atom RMSD (no outlier rejection) and without superposition
    align 1oky, 1t46, cycles=0, transform=0
    

## PyMOL API
    
    
    cmd.align( string mobile, string target, float cutoff=2.0,
               int cycles=5, float gap=-10.0, float extend=-0.5,
               int max_gap=50, string object=None, string matrix='BLOSUM62',
               int mobile_state=0, int target_state=0, int quiet=1,
               int max_skip=0, int transform=1, int reset=0 )
    

This returns a list with 7 items: 

  1. RMSD after refinement
  2. Number of aligned atoms after refinement
  3. Number of refinement cycles
  4. RMSD before refinement
  5. Number of aligned atoms before refinement
  6. Raw alignment score
  7. Number of residues aligned



## Notes

  * The molecules you want to align need to be in **two different objects**. Else, PyMOL will answer with: _ExecutiveAlign: invalid selections for alignment._ You can skirt this problem by making a temporary object and aligning your original to the copy.
  * By defaults, **all states** (like in NMR structures or trajectories) are considered, this might yield a bad or suboptimal alignment for a single state. Use the **mobile_state** and **target_state** argument to be explicit in such cases.



## See Also

  * [super](/index.php/Super "Super"), [cealign](/index.php/Cealign "Cealign"), [fit](/index.php/Fit "Fit"), [pair_fit](/index.php/Pair_fit "Pair fit")
  * [rms](/index.php/Rms "Rms"), [rms_cur](/index.php/Rms_cur "Rms cur"),
  * [intra_fit](/index.php/Intra_fit "Intra fit"), [intra_rms](/index.php/Intra_rms "Intra rms"), [intra_rms_cur](/index.php/Intra_rms_cur "Intra rms cur")
  * [extra_fit](/index.php/Extra_fit "Extra fit")
  * [align_all.py and super_all.py](http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/)
  * [tmalign](/index.php/Tmalign "Tmalign")
  * [Color_by_conservation](/index.php/Color_by_conservation "Color by conservation")
  * [Get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")
  * [mcsalign](/index.php/Mcsalign "Mcsalign") (psico)



Retrieved from "[https://pymolwiki.org/index.php?title=Align&oldid=12703](https://pymolwiki.org/index.php?title=Align&oldid=12703)"


---

## Cealign

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/d/d8/Cealign_ex1.png)](/index.php/File:Cealign_ex1.png)

[](/index.php/File:Cealign_ex1.png "Enlarge")

cealign superposition of 1c0mB and 1bco

cealign aligns two proteins using the CE algorithm. It is very robust for proteins with little to no sequence similarity (twilight zone). For proteins with decent structural similarity, the [super](/index.php/Super "Super") command is preferred and with decent sequence similarity, the [align](/index.php/Align "Align") command is preferred, because these commands are much faster than cealign. 

_This command is new in PyMOL 1.3, see the[cealign plugin](/index.php/Cealign_plugin "Cealign plugin") for manual installation._

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 See Also



## Usage
    
    
    cealign target, mobile [, target_state [, mobile_state
        [, quiet [, guide [, d0 [, d1 [, window [, gap_max
        [, transform [, object ]]]]]]]]]]
    

## Arguments

**Note** : The **mobile** and **target** arguments are swapped with respect to the [align](/index.php/Align "Align") and [super](/index.php/Super "Super") commands. 

  * **target** = string: atom selection of target object (CE uses only CA atoms)
  * **mobile** = string: atom selection of mobile object (CE uses only CA atoms)
  * **target_state** = int: object state of target selection {default: 1}
  * **mobile_state** = int: object state of mobile selection {default: 1}
  * **quiet** = 0/1: suppress output {default: 0 in command mode, 1 in API}
  * **guide** = 0/1: only use "guide" atoms (CA, C4') {default: 1}
  * **d0, d1, window, gap_max** : CE algorithm parameters
  * **transform** = 0/1: do superposition {default: 1}
  * **object** = string: name of alignment object to create {default: (no alignment object)}



## Example
    
    
    fetch 1c0mB 1bco, async=0
    as ribbon
    cealign 1bco, 1c0mB, object=aln
    

## See Also

  * [super](/index.php/Super "Super")
  * [align](/index.php/Align "Align")
  * [cealign plugin](/index.php/Cealign_plugin "Cealign plugin")



Retrieved from "[https://pymolwiki.org/index.php?title=Cealign&oldid=12745](https://pymolwiki.org/index.php?title=Cealign&oldid=12745)"


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

## Fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Fit superimposes the model in the first selection on to the model in the second selection. Only _matching atoms_ in both selections will be used for the fit. 

## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    fit mobile, target [, mobile_state [, target_state [, quiet [, matchmaker [, cutoff [, cycles [, object ]]]]]]]
    

### ARGUMENTS

  * **mobile** = string: atom selection
  * **target** = string: atom selection
  * **mobile_state** = integer: object state {default=0, all states)
  * **target_state** = integer: object state {default=0, all states)
  * **matchmaker** = integer: how to match atom pairs {default: 0} 
    * -1: assume that atoms are stored in the identical order
    * 0/1: match based on all atom identifiers (segi,chain,resn,resi,name,alt)
    * 2: match based on ID
    * 3: match based on rank
    * 4: match based on index (same as -1 ?)
  * **cutoff** = float: outlier rejection cutoff (only if cycles>0) {default: 2.0}
  * **cycles** = integer: number of cycles in outlier rejection refinement {default: 0}
  * **object** = string: name of alignment object to create {default: None}



### EXAMPLES
    
    
    fit ( mutant and name ca ), ( wildtype and name ca )
    

If atom identifiers (like segi, chain, ...) in mobile and target do not match, you need to alter them (or use [Pair_Fit](/index.php/Pair_Fit "Pair Fit") instead): 
    
    
    fetch 1a00, async=0
    extract hbaA, chain A
    extract hbaC, chain C
    # hbaA and hbaC are the same protein, but have different chain identifiers
    alter hbaC, chain='A'
    alter hbaC, segi='A'
    # now both have identical atom identifiers
    fit hbaC, hbaA
    

### SEE ALSO

[Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

Retrieved from "[https://pymolwiki.org/index.php?title=Fit&oldid=13194](https://pymolwiki.org/index.php?title=Fit&oldid=13194)"


---

## Get raw alignment

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_raw_alignment** is an API only function that returns a list of lists of (object,index) tuples containing the raw per-atom alignment relationships. Alignment objects can be created by passing the "object" argument to [align](/index.php/Align "Align") or [super](/index.php/Super "Super"). 

_Please note:_

  * _The order of the atom tuples are not necessarily in the order in which the two (or more) selections were passed to[cmd.align](/index.php/Align "Align")._
  * _Will not return atom tuples of hidden objects (see also[hide_underscore_names](/index.php/Hide_underscore_names "Hide underscore names"))_
  * _Reimplemented in PyMOL 2.3, order of returned atom tuples can differ to previous versions_



## PYMOL API
    
    
    cmd.get_raw_alignment(string name)
    

## EXAMPLE
    
    
    # start a python block
    python
    
    # get two structures
    cmd.fetch('2xwu 2x19', async=0)
    
    # align and get raw alignment
    cmd.align('/2xwu//B//CA', '/2x19//B//CA', cycles=0, transform=0, object='aln')
    raw_aln = cmd.get_raw_alignment('aln')
    
    # print residue pairs (atom index)
    for idx1, idx2 in raw_aln:
        print('%s`%d -> %s`%d' % tuple(idx1 + idx2))
    
    #end the python block
    python end
    

To print residue numbers instead of atom indices: 
    
    
    # continued from previous example
    python
    
    idx2resi = {}
    cmd.iterate('aln', 'idx2resi[model, index] = resi', space={'idx2resi': idx2resi})
    
    # print residue pairs (residue number)
    for idx1, idx2 in raw_aln:
        print('%s -> %s' % (idx2resi[idx1], idx2resi[idx2]))
    
    python end
    

## SEE ALSO

  * [align](/index.php/Align "Align")
  * [find_pairs](/index.php/Find_pairs "Find pairs") returns a list of lists of (object,index) tuples as well
  * [set_raw_alignment](/index.php/Set_raw_alignment "Set raw alignment")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_raw_alignment&oldid=12848](https://pymolwiki.org/index.php?title=Get_raw_alignment&oldid=12848)"


---

## Intra fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**intra_fit** fits all states of an object to an atom selection in the specified state. It returns the rms values to python as an array. 

## Contents

  * 1 USAGE
  * 2 PYMOL API
  * 3 EXAMPLES
    * 3.1 Simple Selection
    * 3.2 Fitting NMR Ensembles
  * 4 PYTHON EXAMPLE
  * 5 USER EXAMPLES
  * 6 USER COMMENTS
  * 7 SEE ALSO



### USAGE
    
    
    intra_fit (selection),state
    

### PYMOL API
    
    
    cmd.intra_fit( string selection, int state )
    

### EXAMPLES

#### Simple Selection
    
    
    intra_fit ( name ca )
    

#### Fitting NMR Ensembles

Warren provided a great example on the PyMOL list. If the NMR ensemble has all the structures loaded as multiple states (which is the default behavoir (see [split_states](/index.php/Split_states "Split states")) for multimodel PDBs) then we can simply call: 
    
    
    print cmd.intra_fit(selection, state)
    

which will fit all of the other states to the indicated states and return the RMS values as a list. 

A simple example is: 
    
    
    fetch 1i8e, async=0
    print cmd.intra_fit("1i8e////CA", 1)
    
    # results:
    [-1.0, 1.1934459209442139, 1.2950557470321655, 0.71329855918884277,
    0.76704370975494385, 0.78973227739334106, 0.99323123693466187,
    1.0165935754776001, 0.6535714864730835, 0.95591926574707031,
    1.1299723386764526, 0.28637325763702393, 0.69836461544036865,
    0.40816938877105713, 1.1637680530548096]
    

Note that a RMS value of -1.0 is returned for the target state. 

### PYTHON EXAMPLE
    
    
    from pymol import cmd
    rms = cmd.intra_fit("(name ca)",1)
    

### USER EXAMPLES

### USER COMMENTS

See [Rms](/index.php/Rms "Rms") for selection caveats for this group of commands. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit"), [Extra_fit](/index.php/Extra_fit "Extra fit"), [Align](/index.php/Align "Align")

Retrieved from "[https://pymolwiki.org/index.php?title=Intra_fit&oldid=12573](https://pymolwiki.org/index.php?title=Intra_fit&oldid=12573)"


---

## Intra rms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**intra_rms** calculates rms fit values for all states of an object over an atom selection relative to the indicated state. Coordinates are left unchanged. The rms values are returned as a python array. 

## Contents

  * 1 PYMOL API
  * 2 PYTHON EXAMPLE
  * 3 USER COMMENTS
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.intra_rms( string selection, int state)
    

### PYTHON EXAMPLE
    
    
    from pymol import cmd
    rms = cmd.intra_rms("(name ca)",1)
    

### USER COMMENTS

See [Rms](/index.php/Rms "Rms") for selection caveats for this group of commands. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

Retrieved from "[https://pymolwiki.org/index.php?title=Intra_rms&oldid=11476](https://pymolwiki.org/index.php?title=Intra_rms&oldid=11476)"


---

## Intra rms cur

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**intra_rms_cur** calculates rms values for all states of an object over an atom selection relative to the indicated state without performing any fitting. The rms values are returned as a python array. 

## Contents

  * 1 PYMOL API
  * 2 PYTHON EXAMPLE
  * 3 USER EXAMPLES/COMMENTS
  * 4 SEE ALSO



### PYMOL API
    
    
    cmd.intra_rms_cur( string selection, int state)
    

### PYTHON EXAMPLE
    
    
    from pymol import cmd
    rms = cmd.intra_rms_cur("(name ca)",1)
    

### USER EXAMPLES/COMMENTS

See [Rms](/index.php/Rms "Rms") for selection caveats for this group of commands. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

Retrieved from "[https://pymolwiki.org/index.php?title=Intra_rms_cur&oldid=11477](https://pymolwiki.org/index.php?title=Intra_rms_cur&oldid=11477)"


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

## Matrix Copy

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Matrix_copy](/index.php/Matrix_copy "Matrix copy") copies the [object matrix](/index.php/Object_Matrix "Object Matrix") from one object to another. 

This command is often used after a protein structure alignment to bring other related objects into the same frame of reference. 

## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Example
  * 3 Example
  * 4 See Also



# Usage
    
    
    matrix_copy source_name, target_name
    

## Arguments

  * source_name = string: name of object to take matrix from
  * target_name = string: name(s) of object(s) to copy matrix to
  * source_mode = integer: 0: raw coordinates, 1: object TTT matrix, 2: state matrix, 3: camera matrix transformation {default: -1: [matrix_mode](/index.php/Matrix_mode "Matrix mode") setting}
  * target_mode = integer: (see source_mode)
  * source_state = integer: object state {default: 1}
  * target_state = integer: object state {default: 1}
  * target_undo = 1/0: ??? {default: 1}



# Example

Here's a practical example. We grab two proteins and their density maps. We then align one to the other and then use matrix_copy to move over the density map: 
    
    
    # fetch two proteins and their maps
    
    fetch 1rx1 3dau, async=0
    fetch 1rx1 3dau, type=2fofc, async=0
    
    # align them proteins
    
    align 1rx1, 3dau
    
    # copy 1rx1's matrix to 1rx1_2fofc
    # it's density map
    
    matrix_copy 1rx1, 1rx1_2fofc
    
    # show the result
    
    enable *
    show extent, *2fofc
    

  


# Example
    
    
    # load two molecules
    load mol1.pdb
    load mol2.pdb
    
    # load their maps
    load map1.ccp4
    load map2.ccp4
    
    # align the two molecules
    align mol2////CA, mol1////CA
    
    # copy mol2's map to mol2
    matrix_copy mol2, map2
    
    # show the isomesh
    isomesh mesh1, map1
    isomesh mesh2, map2
    

# See Also

[ Object Matrix](/index.php/Object_Matrix "Object Matrix") [Matrix_reset](/index.php/Matrix_reset "Matrix reset"), [align](/index.php/Align "Align"), [fit](/index.php/Fit "Fit"), and [pair_fit](/index.php/Pair_fit "Pair fit"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Matrix_Copy&oldid=9441](https://pymolwiki.org/index.php?title=Matrix_Copy&oldid=9441)"


---

## Pair fit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[Pair_Fit](/index.php/Pair_Fit "Pair Fit") fits a set of atom pairs between two models. Each atom in each pair must be specified individually, which can be tedious to enter manually. Script files are recommended when using this command. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 NOTES
  * 4 USER EXAMPLES/COMMENTS
  * 5 SEE ALSO



### USAGE
    
    
    pair_fit (selection), (selection), [ (selection), (selection) [ ...] ]
    

### EXAMPLES
    
    
    # superimpose protA residues 10-25 and 33-46 to protB residues 22-37 and 41-54:
    pair_fit protA///10-25+33-46/CA, protB///22-37+41-54/CA
    # superimpose ligA atoms C1, C2, and C4 to ligB atoms C8, C4, and C10, respectively:
    pair_fit ligA////C1, ligB////C8, ligA////C2, ligB////C4, ligA////C3, ligB////C10
    

### NOTES

So long as the atoms are stored in PyMOL with the same order internally, you can provide just two selections. Otherwise, you may need to specify each pair of atoms separately, two by two, as additional arguments to pair_fit. 

Script files are usually recommended when using this command. 

### USER EXAMPLES/COMMENTS

An description of selection caveats for these commands may be found at [Rms](/index.php/Rms "Rms"). 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms](/index.php/Rms "Rms"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur")

Retrieved from "[https://pymolwiki.org/index.php?title=Pair_fit&oldid=11478](https://pymolwiki.org/index.php?title=Pair_fit&oldid=11478)"


---

## Rms

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Rms computes a RMS fit between two atom selections, but does not tranform the models after performing the fit. 

## Contents

  * 1 USAGE
  * 2 EXAMPLES
  * 3 USER COMMENTS
  * 4 SEE ALSO



### USAGE
    
    
    rms (selection), (target-selection)
    

### EXAMPLES
    
    
    fit ( mutant and name ca ), ( wildtype and name ca )
    

### USER COMMENTS

To determine the RMS without any fitting, see [Rms_Cur](/index.php/Rms_Cur "Rms Cur")

[Fit](/index.php/Fit "Fit"), Rms, [Rms_Cur](/index.php/Rms_Cur "Rms Cur") are finicky and only work when all atom identifiers match: segi, chain, resn, name, alt. If they don't then you'll need to use the alter command to change the identifiers to that they do -- typically that means clearing out the SEGI field, renaming chains, and sometimes renumbering. 

I tried made two selections A, and D as 
    
    
    PyMOL>sel A, 1gh2 and n. CA and i. 65-99
    Selector: selection "A" defined with 35 atoms.
    PyMOL>sel D, 1kao and n. CA and i. 64-98
    Selector: selection "D" defined with 35 atoms
    

which as you can see both yield 35 atoms. Now, 
    
    
    rms_cur A, D
    

won't work, due to the aforementioned reason. To fix this, one needs to do, 
    
    
    alter all,segi=""
    alter all,chain=""
    alter D, resi=str(int(resi)+1)  # I don't actually use this line
    

and now 
    
    
    rms_cur A, D
    

should work. 

### SEE ALSO

[Fit](/index.php/Fit "Fit"), [Rms_Cur](/index.php/Rms_Cur "Rms Cur"), [Intra_Fit](/index.php/Intra_Fit "Intra Fit"), [Intra_Rms](/index.php/Intra_Rms "Intra Rms"), [Intra_Rms_Cur](/index.php/Intra_Rms_Cur "Intra Rms Cur"), [Pair_Fit](/index.php/Pair_Fit "Pair Fit")

[Warren DeLano's comment on rms_* and commands.](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg00001.html)

Retrieved from "[https://pymolwiki.org/index.php?title=Rms&oldid=11473](https://pymolwiki.org/index.php?title=Rms&oldid=11473)"


---

## Rms cur

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

rms_cur computes the RMS difference between two atom selections without performing any fitting. 

By default, only matching atoms in both selections will be used for the fit (same chain, residue number, atoms names etc.). Alternate location mess up the match! 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
    * 3.1 Example 1: Identical Identifiers
    * 3.2 Example 2: Homologues
  * 4 See Also



## Usage
    
    
    rms_cur mobile, target [, mobile_state [, target_state [, quiet [, matchmaker [, cutoff [, cycles [, object ]]]]]]]
    

## Arguments

See [fit](/index.php/Fit "Fit"). 

## Examples

### Example 1: Identical Identifiers

Alternate location mess up the match, so remove them first. 
    
    
    fetch 1p36 1kw7, async=0
    remove not alt ""+"A"
    rms_cur 1p36, 1kw7
    

### Example 2: Homologues

Use [align](/index.php/Align "Align") or [super](/index.php/Super "Super") to create an alignment object (without fitting) and then use the alignment object in the atom selection and turn off identifier matching with **matchmaker=-1**. 
    
    
    fetch 1oky 1t46, async=0
    
    # create alignment object
    align 1oky, 1t46, cycles=0, transform=0, object=aln
    
    # RMSD of alignment object
    rms_cur 1oky & aln, 1t46 & aln, matchmaker=-1
    

## See Also

  * [align](/index.php/Align "Align")
  * [fit](/index.php/Fit "Fit")
  * [rms](/index.php/Rms "Rms")
  * [intra_fit](/index.php/Intra_fit "Intra fit")
  * [intra_rms](/index.php/Intra_rms "Intra rms")
  * [intra_rms_cur](/index.php/Intra_rms_cur "Intra rms cur")
  * [pair_fit](/index.php/Pair_fit "Pair fit")



Retrieved from "[https://pymolwiki.org/index.php?title=Rms_cur&oldid=12632](https://pymolwiki.org/index.php?title=Rms_cur&oldid=12632)"


---

## Set raw alignment

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

set_raw_alignment is an API-only function to create an alignment object from lists of atom indices. 

_New in PyMOL 2.3_

## Arguments

  * **name** = str: alignment object name
  * **raw** = list: list of lists (columns) with `(model, index)` tuples
  * **guide** = str: name of guide object {default: first in list}
  * **state** = int: object state {default: 1}



## Example
    
    
    cmd.align('1t46', '1oky', object='aln')
    raw = cmd.get_raw_alignment('aln')
    cmd.delete('aln')
    cmd.set_raw_alignment('alnnew', raw)
    

## See Also

  * [get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")
  * [align](/index.php/Align "Align")



Retrieved from "[https://pymolwiki.org/index.php?title=Set_raw_alignment&oldid=12847](https://pymolwiki.org/index.php?title=Set_raw_alignment&oldid=12847)"


---

## Super

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**super** aligns two selections. It does a **sequence-independent** (unlike [align](/index.php/Align "Align")) structure-based dynamic programming alignment followed by a series of refinement cycles intended to improve the fit by eliminating pairing with high relative variability (just like [align](/index.php/Align "Align")). **super** is more robust than **align** for proteins with low sequence similarity. 

## Contents

  * 1 Usage
  * 2 Caveats
  * 3 User Scripts
    * 3.1 Write rmsd to file
    * 3.2 Write rmsd to file and looping
  * 4 See Also



## Usage

See [align](/index.php/Align "Align") command. 

## Caveats

  * **Alternative Conformations:** If super ever tells you no matched atoms, then instead of 
        
        super p1, p2
        

try 
        
        super p1 & alt A+'', p2 & alt B+''
        




## User Scripts

### Write rmsd to file

**pymol_rmsd_test.pml**
    
    
    reinitialize
    
    fetch 1F9J, async=0
    fetch 1YX5, async=0
    
    extract 1F9J_A, 1F9J and chain A
    extract 1YX5_B, 1YX5 and chain B
    
    test=cmd.super("1F9J_A","1YX5_B")
    
    python
    writefile=open("rmsd_file.txt","a")
    writefile.write(' '.join('%s' % x for x in test))
    writefile.write('\n')
    writefile.close()
    python end
    

In terminal 
    
    
    pymol -c pymol_rmsd_test.pml ; tail -n 1 rmsd_file.txt
    

### Write rmsd to file and looping

**pymol_pymol_loop.sh**
    
    
    #!/bin/csh -f
    set x = 1
    
    while ( $x <= 20 )
    	set prot="energy_${x}.pdb"
    	pymol -c pymol_super.pml $prot
    	@ x = $x + 1
    end
    

  
**pymol_super.pml**
    
    
    reinitialize
    import sys
    
    python 
    prot1="XXXX"
    prot2="YYYY_trimmed"
    
    prot3=sys.argv[3]
    #prot3="energy_54.pdb"
    prot3name=prot3.split(".pdb")[0]
    print prot3, prot3name
    python end
    
    cd /XXXXX
    
    cmd.load("%s.pdb"%prot1)
    cmd.load("%s.pdb"%prot2)
    cmd.load("./ensemblesize2_numstruct/%s"%prot3)
    #show_as cartoon, all
    
    align1=cmd.super("%s"%prot3name,"%s"%prot1)
    print align1
    
    python
    writefile=open("pymol_rmsd_file.txt","a")
    writefile.write('%s %s '%(prot3name, prot1))
    writefile.write(' '.join('%s' % x for x in align1))
    writefile.write(' ')
    python end
    
    align2=cmd.super("%s"%prot3name,"%s"%prot2)
    print align2
    
    python
    writefile=open("pymol_rmsd_file.txt","a")
    writefile.write('%s %s '%(prot3name, prot2))
    writefile.write(' '.join('%s' % x for x in align2))
    writefile.write('\n')
    writefile.close()
    python end
    

In terminal 
    
    
    chmod +x pymol_loop.sh
    ./pymol_loop.sh
    

## See Also

  * [align](/index.php/Align "Align")
  * [Cealign](/index.php/Cealign "Cealign")
  * [Get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")



Retrieved from "[https://pymolwiki.org/index.php?title=Super&oldid=12180](https://pymolwiki.org/index.php?title=Super&oldid=12180)"


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

