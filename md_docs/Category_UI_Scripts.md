# Category: UI Scripts

## Aa codes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Just a quick little script to allow you to convert from 1-to-3 letter codes and 3-to-1 letter codes in PyMOL. Copy the code below and drop it into your .pymolrc file. Then, each time you load PyMOL, "one_letter" and "three_letter" will be defined. 

## Contents

  * 1 The Code
    * 1.1 Simple
    * 1.2 Simple and Clever
    * 1.3 Using BioPython
    * 1.4 Using PyMOL
  * 2 Example Usage



# The Code

## Simple
    
    
    # one_letter["SER"] will now return "S"
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}
    
    # three_letter["S"] will now return "SER"
    three_letter = dict([[v,k] for k,v in one_letter.items()])
    
    three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
    'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
    'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
    'G':'GLY', 'P':'PRO', 'C':'CYS'}
    

## Simple and Clever

Here's another way to accomplish this 
    
    
    # The real convenience in there is that you can easily construct any
    # kind of hash by just adding a matching list, and zipping.
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    aa123 = dict(zip(aa1,aa3))
    aa321 = dict(zip(aa3,aa1))
    
    # Then to extract a sequence, I tend to go for a construction like:
    sequence = [ aa321[i.resn] for i in cmd.get_model(selection + " and n. ca").atom ]
    

## Using BioPython

If you have BioPython you can use the following, which includes also many three-letter codes of modified amino acids: 
    
    
    from Bio.PDB import to_one_letter_code as one_letter
    
    
    
    from Bio.Data.SCOPData import protein_letters_3to1 as one_letter
    

## Using PyMOL
    
    
    from pymol.exporting import _resn_to_aa as one_letter
    

# Example Usage
    
    
    # we used to have to do the following to get the amino acid name
    from pymol import stored
    stored.aa = ""
    cmd.iterate("myselection", "stored.aa=resn")
    
    # now we can just call
    three_letter[string.split(cmd.get_fastastr("myselection"),'\n')[1]]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Aa_codes&oldid=12634](https://pymolwiki.org/index.php?title=Aa_codes&oldid=12634)"


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

