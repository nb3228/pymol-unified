# Category: ThirdParty Scripts

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

## PovRay

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL can export input files for [POV-Ray](http://www.povray.org/) with the ".pov" file extension: 
    
    
    set stick_ball
    save input.pov
    

It can also use POV-Ray directly for rendering with the [ray](/index.php/Ray#Renderer "Ray") command: 
    
    
    ray renderer=1
    

Since PyMOL 1.7.4, round stick caps are only exported correctly with [stick_ball](/index.php/Stick_ball "Stick ball")=on. 

## Nice PovRay settings

I typically use the make_pov.py script and "run" it from pymol once to load the function, and then I do _make_pov('povray.inp')_ to create the povray.inp file. Then I edit that file to insert some lines like: 
    
    
    fog {
      distance 10
      fog_type 2
      fog_alt 10.
      fog_offset -160.
      up <0.,1.,.4>
    colour rgbt<1.0, 1.0, 1.0, 0.1>
    turbulence 0.8
    }
    

In this case I'm not really doing depth-cueing but adding fog at the lower background edge (there were two planes defining the background and a surface below the molecule) rising up towards the front upper edge of the scene. 

"fog_type 2" means a "rising fog" along the "up" vector. fog_type 1 is a constant fog. To get pure depth cueing, you would want "up" to be along the <0., 0., 1.> vector (I think!). You'll need to play around with the distance and fog_offset parameters. You wouldn't necessarily want the "turbulence" parameter in there either. Check out "Atmospheric Effects" in the povray documentation for many more details: <http://www.povray.org/documentation/view/201/>

## make_pov.py v1
    
    
    # make_pov.py
    # Do "run make_pov.py" from within pymol and then execute the script
    # with "make_pov('povray.inp')" to create the povray.inp file.
    #
    # written by Robert Campbell 2003
    #
    from pymol import cmd
    
    def make_pov(file):
    	(header,data) = cmd.get_povray()
    	povfile=open(file,'w')
    	povfile.write(header)
    	povfile.write(data)
    	povfile.close()
    

## make_pov.py v2

This is a more extended version of an earlier extension of the version by Robert Campbell. The scene is written in two parts, a .pov file containing all meta data, such as the lights, camera and #defaults, and an include file (.inc) which contains the structure. In this way you have maximum control over your scene without having to edit a huge povray file. You may even want to consider splitting your scene up in separate parts (taken from the same perspective), which you combine in a global .pov file using #include statements. This will give even more control with regards to modifications to the scene. If 'clip' is set to near|far|both, then the corresponding clipping plane(s) is/are included in a CSG difference object. Note that the result may increase the render time significantly unless the scene is simple. 

Once you run **run make_pov.py** , run **make_pov** to execute the script. 

NB. the .pov file contains a commented statement with regards to a povray macro file, which allows transforming scenes and objects from model space to camera space and vice versa. The macro file is given below. 
    
    
    # make_pov.py
    # Do "run make_pov.py" from within pymol and then execute the script
    # with "make_pov('povray.inp')" to create the povray.inp file.
    #                                                                                                   
    # Original script written by Robert Campbell
    # Modified by Tsjerk A. Wassenaar
    #
    
    from pymol import cmd
    
    def make_pov(file, name="PymolObject", meta=True, clip=False ):
            f1, f2 = file, file[:-4] + '.inc'
    
            (header,data) = cmd.get_povray()
            povfile = open(f1,'w')
            if meta: povfile.write(header)
            povview = cmd.get_view()
    
            if clip:
                    objtype = "difference"
                    objclip = ""
                    if clip in ["near","both"]:
                            objclip = objclip + "plane { z, -%f }" % povview[15] 
                    if clip in ["far","both"]:
                            objclip = objclip + "plane { z, -%f }" % povview[16] 
            else:
                    objtype = "object"
                    objclip = ""
    
            povfile.write("""\n
    // Uncomment the following lines if you have the pymolmacro.inc include file and want to use it.
    /*
    #include \"pymolmacro.inc\"
    PYMOL_VIEW( %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f )
    */
    
    """ % povview)
    
            povfile.write("""
    #declare %s = union { #include "%s" }
    %s { %s %s }
    """ % (name, f2, objtype, name, objclip ) )
    
            povfile.close()
            povfile = open(f2,'w')
            povfile.write(data)
            povfile.close()
    
    cmd.extend('make_pov',make_pov)
    
    
    
    //
    //  PYMOLMACRO.INC v0.2 
    //
    //  (c)2005 Tsjerk Wassenaar, University of Groningen
    //
    //  This include file for Povray contains
    //  just a few macros which together allow
    //  the conversion between the model space
    //  (cartesian coordinates) and the Pymol
    //  camera space.
    //
    //  With these macros one can easily combine
    //  a Pymol scene with objects defined in the
    //  coordinate space of the original
    //  structure file.
    //
    //  The input consists of the output of the
    //  get_view() command in Pymol. This output
    //  consists of 18 floating point numbers
    //  defining a rotation matrix and shift
    //  vectors for the origin of rotation and
    //  for the camera position.
    //
    //  The macro PYMOL_VIEW loads a
    //  view obtained from Pymol.
    //
    //  It #declares two transformation statements:
    //
    //  FROM_PYMOL_VIEW
    //  TO_PYMOL_VIEW
    //
    //  The first can be used to transform the Pymol
    //  scene back to model (normal) space, the latter
    //  is used to transform other objects to appear in
    //  the scene on the correct position.
    //
    //  Additionally four macros are defined to transform
    //  vectors (points) from one space to another:
    //
    //  VEC2PYMOLSPACE( <x, y, z> )
    //  VEC2CARTSPACE( <x, y, z> )
    //  VEC2PYMOLVEC( <x, y, z> )
    //  VEC2CARTVEC( <x, y, z> ) 
    //
    //  *NEW*
    //
    //  If the view from pymol is stored as an array:
    //
    //  #declare M = array[18] {...}
    //
    //  then the macros
    //
    //  SET_PYMOL_VIEW 
    //    and 
    //  UNSET_PYMOL_VIEW
    //
    //  can be used directly to transform objects to and from that view:
    //  object { ... SET_PYMOL_VIEW( M ) }
    //
    //  This is especially useful if multiple views are defined 
    //  and the scene was set in one:
    //
    //  #declare VIEW1 = M;
    //  #declare VIEW2 = N;
    //  union { #include "file.inc" UNSET_PYMOL_VIEW( M ) SET_PYMOL_VIEW( N ) }
    //
    //  NOTE: transform statements are combined by POV-Ray prior to 
    //  transformations of objects, so there's little need to implement a macro
    //  SWITCH_PYMOL_VIEW( M, N )
    //  although that would appear simpler in the scenes  
    
    //  Tsjerk A. Wassenaar
    //  February 16, 2005
    //  April 5, 2005
    //  September 2, 2009
    //
    
    // Determinant of a matrix
    //------------------------
    #macro DET( M )
    
      #local a = M[0] * ( M[4]*M[8] - M[5]*M[7] ); 
      #local b = M[1] * ( M[3]*M[8] - M[5]*M[6] ); 
      #local c = M[2] * ( M[3]*M[7] - M[4]*M[6] );
    
      (a - b + c)
    
    #end // of DET()
    
    
    // The inverse of a matrix
    //------------------------
    #macro INV( m11, m12, m13, m21, m22, m23, m31, m32, m33 )
    
      #local M = array[9] { m11, m12, m13, m21, m22, m23, m31, m32, m33 };
      #local invdet = 1/DET( M );
    	
      #local t11 = invdet * ( m22*m33 - m23*m32 ); 
      #local t12 = invdet * ( m13*m32 - m12*m33 ); 
      #local t13 = invdet * ( m12*m23 - m13*m22 ); 
      #local t21 = invdet * ( m23*m31 - m21*m33 ); 
      #local t22 = invdet * ( m11*m33 - m13*m31 ); 
      #local t23 = invdet * ( m13*m21 - m11*m23 );
      #local t31 = invdet * ( m21*m32 - m22*m31 );
      #local t32 = invdet * ( m12*m31 - m11*m32 );
      #local t33 = invdet * ( m11*m22 - m12*m21 );
    
      t11, t12, t13, t21, t22, t23, t31, t32, t33, 0, 0, 0
    
    #end // of INV()
    
    #macro M_INV( M )
      #local invdet = 1/DET( M );
    	
      #local t11 = invdet * ( M[4]*M[8] - M[5]*M[7] ); 
      #local t21 = invdet * ( M[2]*M[7] - M[1]*M[8] ); 
      #local t31 = invdet * ( M[1]*M[5] - M[2]*M[4] ); 
    
      #local t12 = invdet * ( M[5]*M[6] - M[3]*M[8] ); 
      #local t22 = invdet * ( M[0]*M[8] - M[2]*M[6] ); 
      #local t32 = invdet * ( M[2]*M[3] - M[0]*M[5] );
    
      #local t13 = invdet * ( M[3]*M[7] - M[4]*M[6] );
      #local t23 = invdet * ( M[1]*M[6] - M[0]*M[7] );
      #local t33 = invdet * ( M[0]*M[4] - M[1]*M[3] );
    
      array[9] {t11, t12, t13, t21, t22, t23, t31, t32, t33}
    #end
    
    #macro MV_MUL( M, V )
        < M[0]*V.x + M[1]*V.y + M[2]*V.z,
          M[3]*V.x + M[4]*V.y + M[5]*V.z,
          M[6]*V.x + M[7]*V.y + M[8]*V.z >
    #end
    
    #macro SET_PYMOL_VIEW( M )
      transform {
        translate -< M[12], M[13], M[14] >
        matrix < M[0], M[1],  M[2],
    	     M[3], M[4],  M[5], 
    	     M[6], M[7],  M[8], 
    	     M[9], M[10], M[11] >
      } 
    #end // of SET_PYMOL_VIEW
    
    #macro UNSET_PYMOL_VIEW( M )
      transform {
        translate -< M[9], M[10], M[11] >
        matrix < INV( M[0], M[1], M[2], M[3], M[4], M[5], M[6], M[7], M[8] ) > 
        translate < M[12], M[13], M[14] >
      } 
    #end // of UNSET_PYMOL_VIEW
    
    #macro C2P_VEC( M, vec)
      #local nvec = vec - <M[12],M[13],M[14]>;
      #local nvec =
        < M[0]*nvec.x + M[1]*nvec.y + M[2]*nvec.z,
          M[3]*nvec.x + M[4]*nvec.y + M[5]*nvec.z,
          M[6]*nvec.x + M[7]*nvec.y + M[8]*nvec.z >; 
      nvec + <M[9],M[10],M[11]>
    #end
    
    #macro P2C_VEC( M, vec)
      MV_MUL( M_INV(M), vec - <M[9],M[10],M[11]> ) + <M[12],M[13],M[14]>
      //#local nvec = vec - <M[9],M[10],M[11]>;
      //#local N = M_INV( M ) ;
      //#local nvec =
      //  < N[0]*nvec.x + N[1]*nvec.y + N[2]*nvec.z,
      //    N[3]*nvec.x + N[4]*nvec.y + N[5]*nvec.z,
      //    N[6]*nvec.x + N[7]*nvec.y + N[8]*nvec.z >; 
      //nvec
    #end
    
    
    #macro PYMOL_VIEW( r11, r12, r13,     // 3x3 Rotation matrix ( Model space to Camera space )
    		   r21, r22, r23, 
    		   r31, r32, r33,
    		    c1,  c2,  c3,     // Camera position ( Model space )
    		    o1,  o2,  o3,     // Origin of rotation ( Model space )
    		    s1,  s2,  or)     // Slab near and far, orthoscopic flag ( discarded )
    
      #declare PYMOLVIEW_RMATRIX = array[9] { r11, r12, r13, 
    					  r21, r22, r23, 
    					  r31, r32, r33 }
      #declare PYMOLVIEW_CAMPOS  = < c1, c2, c3 >;
      #declare PYMOLVIEW_ORGPOS  = < o1, o2, o3 >;
    
      #declare TO_PYMOL_VIEW = transform {
        translate -< o1, o2, o3 >
        matrix < r11, r12, r13,
    	     r21, r22, r23, 
    	     r31, r32, r33, 
    	      c1,  c2,  c3 >
      }
    
      #declare FROM_PYMOL_VIEW = transform {
        translate -< c1, c2, c3>
        matrix < INV( r11, r12, r13, r21, r22, r23, r31, r32, r33 ) >
        translate  < o1, o2, o3>
      }
    
      #macro VEC2PYMOLSPACE(vec)
        #local nvec = vec - PYMOLVIEW_ORGPOS;
        #local nvec =
          < PYMOLVIEW_RMATRIX[0]*nvec.x + PYMOLVIEW_RMATRIX[3]*nvec.y + PYMOLVIEW_RMATRIX[6]*nvec.z,
            PYMOLVIEW_RMATRIX[1]*nvec.x + PYMOLVIEW_RMATRIX[4]*nvec.y + PYMOLVIEW_RMATRIX[7]*nvec.z,
            PYMOLVIEW_RMATRIX[2]*nvec.x + PYMOLVIEW_RMATRIX[5]*nvec.y + PYMOLVIEW_RMATRIX[8]*nvec.z >; 
        nvec + PYMOLVIEW_CAMPOS
      #end
    
      #macro VEC2CARTSPACE(vec)
    
        #local nvec = vec - PYMOLVIEW_CAMPOS;
    
        #local R = PYMOLVIEW_RMATRIX;
        #local invdet = 1/DET( R );
    
        #local T = array[9];
    
        #local T[0] = invdet * ( R[4]*R[8] - R[5]*R[7] ); 
        #local T[1] = invdet * ( R[2]*R[7] - R[1]*R[8] ); 
        #local T[2] = invdet * ( R[1]*R[5] - R[2]*R[4] ); 
        #local T[3] = invdet * ( R[5]*R[6] - R[3]*R[8] ); 
        #local T[4] = invdet * ( R[0]*R[8] - R[2]*R[6] ); 
        #local T[5] = invdet * ( R[2]*R[3] - R[0]*R[5] );
        #local T[6] = invdet * ( R[3]*R[7] - R[4]*R[6] );
        #local T[7] = invdet * ( R[1]*R[6] - R[0]*R[7] );
        #local T[8] = invdet * ( R[0]*R[4] - R[1]*R[3] );
    
        < T[0]*nvec.x + T[3]*nvec.y + T[6]*nvec.z + PYMOLVIEW_ORGPOS.x,
          T[1]*nvec.x + T[4]*nvec.y + T[7]*nvec.z + PYMOLVIEW_ORGPOS.y,
          T[2]*nvec.x + T[5]*nvec.y + T[8]*nvec.z + PYMOLVIEW_ORGPOS.z >
      #end
    
      #macro VEC2PYMOLVEC(vec)
        < PYMOLVIEW_RMATRIX[0]*vec.x + PYMOLVIEW_RMATRIX[3]*vec.y + PYMOLVIEW_RMATRIX[6]*vec.z,
          PYMOLVIEW_RMATRIX[1]*vec.x + PYMOLVIEW_RMATRIX[4]*vec.y + PYMOLVIEW_RMATRIX[7]*vec.z,
          PYMOLVIEW_RMATRIX[2]*vec.x + PYMOLVIEW_RMATRIX[5]*vec.y + PYMOLVIEW_RMATRIX[8]*vec.z >
      #end
    
      #macro VEC2CARTVEC(vec)
        #local nvec = vec - PYMOLVIEW_CAMPOS;
      #end
    
      #macro CAM2PYMOLCAM()
    
      #end
    
      #macro CAM2CARTCAM()
    
      #end
    #end
    

Retrieved from "[https://pymolwiki.org/index.php?title=PovRay&oldid=12676](https://pymolwiki.org/index.php?title=PovRay&oldid=12676)"


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

