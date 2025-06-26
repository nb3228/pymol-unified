# Category: Mac

## MAC Install

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes how to install PyMOL on Mac OS X. 

## Contents

  * 1 Incentive PyMOL
    * 1.1 Launching from Command Line
    * 1.2 X11 Hybrid
    * 1.3 Stereo on Second Monitor
  * 2 Open-Source PyMOL
    * 2.1 Package managers
    * 2.2 Install from Source
    * 2.3 Install APBS with Fink
    * 2.4 Stereo issues
  * 3 See Also



## Incentive PyMOL

[Schrödinger](http://www.schrodinger.com) provides pre-compiled PyMOL to paying sponsors. The bundle also includes ready-to-use [APBS](/index.php/APBS "APBS"), [RigiMOL](/index.php/Morph "Morph"), an MPEG encoder for movie export, and a small molecule energy minimization engine. 

Download: <https://pymol.org/>

Installation: Drag **PyMOL.app** on the **/Applications** shortcut. (In principle, you could drag it into any Finder window and run it from there, it doesn’t have to live in /Applications). 

Uninstallation: Move **/Applications/PyMOL.app** to Trash 

### Launching from Command Line

The unix executable resides at **/Applications/PyMOL.app/Contents/MacOS/PyMOL**

### X11 Hybrid

_Applies to PyMOL 1.x, not to PyMOL 2.x_

MacPyMOL can optionally run with the same two-window GUI which PyMOL uses on Windows and Linux. This GUI has some additional features, like the [Plugin](/index.php/Plugins "Plugins") menu and the [Builder](/index.php/Builder "Builder"). 

Requires [XQuartz](http://xquartz.macosforge.org/). 

There are two ways to launch the X11 interface: 

  1. Rename or copy/duplicate **/Applications/MacPyMOL.app** to **/Applications/MacPyMOLX11Hybrid.app** or to **/Applications/PyMOLX11Hybrid.app**
  2. Launch the unix executable with the **-m** flag: **/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL -m**



### Stereo on Second Monitor

The trick to getting MacPyMOL to work in stereo on the second monitor is to force it to initially open on that display by providing an appropriate "-X #" (and perhaps -Y #) option on launch. That way the OpenGL context will be created with stereo support. 
    
    
    ./MacPyMOL.app/Contents/MacOS/MacPyMOL -X -1000
    ./MacPyMOL.app/Contents/MacOS/MacPyMOL -X -1000 -Y 100
    

**Source:** [Warren DeLano; PyMOL Users Archive](https://sourceforge.net/p/pymol/mailman/message/11671952/)

## Open-Source PyMOL

### Package managers

Open-Source PyMOL is available [free of charge](https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE) and may be readily installed via the [Homebrew](http://brew.sh/) (recommended), [MacPorts](https://www.macports.org/), or [Fink](http://www.finkproject.org/) package managers. 
    
    
    # Homebrew (recommended)
    brew install brewsci/bio/pymol
    
    # Fink
    fink install pymol-py27
    
    # MacPorts
    sudo port install pymol
    

You may need to make sure that the dependencies are installed with the required flags, e.g. for MacPorts: 
    
    
    # MacPorts
    sudo port install tcl -corefoundation
    sudo port install tk -quartz
    

If PyMOL complains that it wasn't able to find X11, try starting xquartz first, then run pymol from the console. 

### Install from Source

If you want the latest PyMOL code (warning: might include experimental changes), then follow the [Linux installation instructions](/index.php/Linux_Install#Install_from_source "Linux Install"). You will need an environment like Fink, MacPorts or Homebrew to install the dependencies. Make sure you use the appropriate python interpreter (e.g. **/sw/bin/python2.7** when using Fink). 

To run PyMOL with a native PyQt library (linked against macOS OpenGL framework, not against XQuartz), it needs to be built with the `--osx-frameworks` option: 
    
    
    python setup.py --osx-frameworks install
    

### Install APBS with Fink

To use the electrostatics plugin, you will need [APBS](http://apbs.sourceforge.net/) and its dependencies. These are also available as Fink packages, and include [APBS](http://pdb.finkproject.org/pdb/package.php/apbs), [maloc](http://pdb.finkproject.org/pdb/package.php/maloc) and [pdb2pqr](http://pdb.finkproject.org/pdb/package.php/pdb2pqr). If you have multiple processors available, you might wish to install the [MPI version of APBS](http://pdb.finkproject.org/pdb/package.php/apbs-mpi-openmpi). 

Issuing the command 
    
    
    fink install apbs
    

will install apbs and its required dependencies for you. The fink pymol package is already preconfigured to do the right thing to use apbs as a plugin. 

### Stereo issues

Some older Macs seem to crash with stereo graphics. If this happens to you, a workaround is to launch PyMOL explicitly in Mono mode with `pymol -M`. You can also set up an alias in your ~/.profile: 
    
    
    alias pymol='pymol -M'
    

## See Also

  * [pymolrc](/index.php/Pymolrc "Pymolrc")
  * [Linux Install](/index.php/Linux_Install "Linux Install")
  * [Windows Install](/index.php/Windows_Install "Windows Install")
  * [FreeMOL installation for MPEG movie export](/index.php/MovieSchool_6#Exporting_your_Movie "MovieSchool 6")
  * [Bill Scott’s](/index.php/User:Wgscott "User:Wgscott") [MacOSX-specific .pymolrc file](/index.php/MacOSX-specific_.pymolrc_file "MacOSX-specific .pymolrc file") and his crystallographic software [wiki](http://xanana.ucsc.edu/~wgscott/xtal/wiki/index.php/Main_Page) and [website](http://chemistry.ucsc.edu/~wgscott/xtal/), including instructions on [how to install precompiled binary packages using fink](http://xanana.ucsc.edu/~wgscott/xtal/wiki/index.php/Getting_your_fink_installation_to_use_packages_that_I_have_pre-compiled).
  * [Launching_PyMOL#MacOS_X](/index.php/Launching_PyMOL#MacOS_X "Launching PyMOL")



Retrieved from "[https://pymolwiki.org/index.php?title=MAC_Install&oldid=12800](https://pymolwiki.org/index.php?title=MAC_Install&oldid=12800)"


---

## MacOSX-specific .pymolrc file

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Here is an example of a ~/.pymolrc file that [I](/index.php/User:Wgscott "User:Wgscott") use on Mac OS X. There are two main things going on here: 

  * It defines a bunch of aliases, including one that ribbonizes protein/nucleic acid comples (this could benefit from updating)


  * It runs a series of pythonized pymol commands to permit interaction with the [ PowerMate Dial](/index.php/MAC_Install#PowerMate_Dial "MAC Install")



  

    
    
     
    _ feedback push
    _ feedback disable,all,everything
    
    alias clear, mstop; mclear; hide all
    alias nogui, set internal_gui=0
    alias gui, set internal_gui=1
    alias shiny, set spec_power=250; set spec_refl=1.5; set antialias=1; ray
    alias grab, os.system("open -a Grab")
    alias stop, quit
    alias exit, quit
    alias white, bg_color white; set depth_cue=0; set ray_trace_fog=0
    
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
    
    
    alias ribbonize,\
        hide; show cartoon; set cartoon_fancy_helices, 1; \
        hide cartoon, ( resname A or resname G or resname C or resname U or resname T ); \
        show sticks, ( resname A or resname G or resname C or resname U or resname T ); \
        util.rainbow; color hotpink, ( resname A or resname G or resname C or resname U or resname T )
    
    
    # END COMMANDS
    _ feedback pop
    
    print "PowerMate Dial interface has been enabled"
    
    print "Finished reading ~/.pymolrc"
    

Retrieved from "[https://pymolwiki.org/index.php?title=MacOSX-specific_.pymolrc_file&oldid=5375](https://pymolwiki.org/index.php?title=MacOSX-specific_.pymolrc_file&oldid=5375)"


---

## Mouse Controls

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## MacPyMol And PyMolX11 Hybrid

We strongly recommend (require?) use of a three-button mouse with MacPyMOL and PyMOLX11Hybrid. However, it is possible to use MacPyMOL in a limited way on Macs that have a single button mouse thanks to some built-in mouse remapping in the Mac OS X GLUT implementation. Her is how that works... 

If the Mac is hooked up to a 1-button mouse (only) when MacPyMOL or HybridX11PyMOL are launched, then Mac OS X itself will furnish translations as follows: 
    
    
    Rotate:       Click & Drag
    XY-Translate: Option-Click & Dreg
    Zoom:         Control-Click & Drag
    Select:       Click & Release
    Box-Select:   Shift-Click & Drag
    Box-Deselect: Shift-Option-Click & Drag
    Clipping:    Control-Shift-Click & Drag with...
     * near plane controlled by vertical motion
     * far plane controlled by horizontal motion
    

Note that not all mouse actions are possible with a one-button mouse. For example, I don't think molecular editing is possible with a one-button mouse. 

# See Also

[Configuring the Mouse](/index.php/Config_Mouse "Config Mouse")

Retrieved from "[https://pymolwiki.org/index.php?title=Mouse_Controls&oldid=8663](https://pymolwiki.org/index.php?title=Mouse_Controls&oldid=8663)"


---

