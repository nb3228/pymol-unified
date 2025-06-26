# Category: Technical Issues

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

## Windows Install

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes how to install PyMOL on Microsoft Windows. 

## Contents

  * 1 Incentive PyMOL
  * 2 Open-Source PyMOL
    * 2.1 One-click Installer
    * 2.2 Pre-compiled Open-Source PyMOL
    * 2.3 Building from source
    * 2.4 Extend PyMOL with additional scripts
  * 3 See Also



## Incentive PyMOL

[Schrödinger](http://www.schrodinger.com) provides an installer to paying sponsors (EXE for PyMOL 2.0, MSI for previous version). The bundle also includes ready-to-use [APBS](/index.php/APBS "APBS"), [RigiMOL](/index.php/Morph "Morph"), an MPEG encoder for movie export, and a small molecule energy minimization engine. 

Download: <https://pymol.org/#download>

## Open-Source PyMOL

Open-Source PyMOL is available [free of charge](https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE). It also allows sponsors to create highly customized PyMOL installations which might not be possible with the MSI installer. 

### One-click Installer

A convenient and easy-to-use one-click installer is available that runs out-of-the-box without the need for an existing Python installation. You can download it [here](https://github.com/kullik01/pymol-open-source-windows-setup/releases/download/v3.1.0/PyMOL_Open_source_v3.1.0a0_WINx64_setup.exe) or go to the [GitHub repository](https://github.com/kullik01/pymol-open-source-windows-setup) and then to [Releases](https://github.com/kullik01/pymol-open-source-windows-setup/releases/tag/v3.1.0). It features an alternate splash screen and a slightly updated menu stylesheet. 

### Pre-compiled Open-Source PyMOL

Pre-compiled Open-Source PyMOL is available free from [Christoph Gohlke of the Laboratory for Fluorescence Dynamics, University of California, Irvine](https://www.cgohlke.com/). 

  1. Install the latest version of Python 3 for Windows (e.g., by going to <http://www.python.org/downloads/> and choosing the x64 EXE installer). Use the standard options, which should mean that the installation directory is most likely C:\Users\<Your Username>\AppData\Local\Programs\Python\Python38). Make sure the option to add environment variables is selected or add the folder of python.exe to system PATH.
  2. Install [the current Microsoft Visual C++ Redistributable for Visual Studio 2015, 2017 and 2019](https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads). Otherwise the installed PyMOL binary may fail to run (without any error message!). 
     1. Have a look in "Add/Remove" programs if it's installed already
  3. Download the [appropriate wheel files](https://github.com/cgohlke/pymol-open-source-wheels/), along with all requirement wheel files ([Numpy+MKL](https://github.com/cgohlke/numpy-mkl-wheels/)) into a single file directory, e.g., `C:\Users\<Your Username>\Downloads`
  4. Example of filenames 2023-01-12 
     1. pymol_launcher-2.5-cp311-cp311-win_amd64.whl
     2. pymol-2.6.0a0-cp311-cp311-win_amd64.whl
     3. numpy-1.22.4+mkl-cp311-cp311-win_amd64.whl



Navigate to the installation directory in a CMD window (Not PowerShell!) 
    
    
    cd C:\Users\<Your Username>\Downloads
    

(or where ever you put the files) and begin the installation using the command: ("%CD%" only works in CMD) 
    
    
    python -m pip pmw
    python -m pip install --no-index --find-links="%CD%" pymol_launcher-2.5-cp311-cp311-win_amd64.whl
    where.exe pymol
    pymol
    

PyMOL.exe should now be in C:\Users\<Your Username>\AppData\Local\Programs\Python 

To use the newer single-window Qt interface, also install the optional PyQt5 dependency for your Python installation: 
    
    
    python -m pip install pyqt5
    

To update PyMOL update the files in the PyMOL install directory and run: 
    
    
    pip install --upgrade --no-deps pymol.whl
    

where `pymol.whl` is replaced by the PyMOL wheel file name (not the launcher, the launcher should not require updating). 

### Building from source

This [GitHub repository](https://github.com/urban233/pymol-open-source-windows-build) provides the tools to build PyMOL for Windows from source. It automatically downloads all necessary dependencies via the vcpkg package manager. You can either build the complete wheel file or just the _cmd C/C++ extension module using the provided automation script. A pre-built wheel file is available for download directly from the [Releases](https://github.com/urban233/pymol-open-source-windows-build/releases/tag/v3.1.0a0) tab. 

### Extend PyMOL with additional scripts

If you now want to extend the capabilities of PyMOL, and take advantage of all the available plugins+scripts "out there", then do the following.   


  1. First install "numpy" as an available module to Python. [Select appropriate installer from here](https://github.com/cgohlke/numpy-mkl-wheels/)
  2. Download the script/plugin collection [ Pymol-script-repo](/index.php/Git "Git") from [a .zip file from here](https://github.com/Pymol-Scripts/Pymol-script-repo/zipball/master)


    
    
    git clone <https://github.com/Pymol-Scripts/Pymol-script-repo>
    

  1. Unpack it to here: **C:\Python27\Lib\site-packages\pymol\pymol_path\Pymol-script-repo** Double check that the folder name is correct and the same.



Open "Notepad" and write. 
    
    
    # Add paths to sys.path so PyMOL can find modules and scripts
    import sys, os
    pymol_git = os.path.abspath(os.path.join(os.environ['PYMOL_PATH'], 'Pymol-script-repo'))
    os.environ['PYMOL_GIT_MOD'] = os.path.join(pymol_git,'modules')
    sys.path.append(pymol_git)
    sys.path.append(os.environ['PYMOL_GIT_MOD'])
    
    # Make setting changes to Plugin Manager
    import pymol.plugins
    pymol.plugins.preferences = {'instantsave': False, 'verbose': False}
    pymol.plugins.autoload = {'apbs_tools': False}
    pymol.plugins.set_startup_path([os.path.join(pymol_git, 'plugins'), os.path.join(sys.prefix, 'Lib', 'site-packages', 'pmg_tk', 'startup')])
    pymol.plugins.preferences = {'instantsave': True, 'verbose': False}
    

**Then "File- >Save as->All files-> C:\Python27\Lib\site-packages\pymol\pymol_path\run_on_startup.py**

Now start pymol, and enjoy all the plugins available from the menu. 

**PyMOL** shortcut  
Make a **pymol** directory in your homepath. **mkdir %HOMEPATH%\pymol** Then make sure, PyMOL starts here, when you open the shortcut.  
Make a shortcut to the .cmd file, and modify it.   
Target: C:\python27\PyMOL\pymol.cmd   
Start in: %HOMEPATH%\pymol 

## See Also

  * [pymolrc](/index.php/Pymolrc "Pymolrc")
  * [Linux Install](/index.php/Linux_Install "Linux Install")
  * [MAC Install](/index.php/MAC_Install "MAC Install")



Retrieved from "[https://pymolwiki.org/index.php?title=Windows_Install&oldid=13838](https://pymolwiki.org/index.php?title=Windows_Install&oldid=13838)"


---

## XFree86 Configuration

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Iiyama HM204DT + Nvidia Quadro FX 1100 on Fedora Core 2 (updated as of 2005-02-17)

` `

`
    
    
    Section "Module"
           Load  "dbe"
           Load  "extmod"
           Load  "fbdevhw"
           Load  "glx"
           Load  "record"
           Load  "freetype"
           Load  "type1"
           #Load  "dri"    # THIS MUST BE COMMENTED OUT
    EndSection
    

``
    
    
    Section "Monitor"
      Identifier   "hm204dt"
      VendorName   "Iiyama"
      ModelName    "HM204DT"
      HorizSync    30 - 142
      VertRefresh  50 - 200
      Option      "dpms"
    
      # 1400x1100 @ 120.00 Hz (GTF) hsync: 141.48 kHz; pclk: 275.04 MHz
      Modeline "1400x1100_120"  275.04  1400 1520 1672 1944  1100 1101 1104 1179  -HSync +Vsync
    EndSection
    

``
    
    
    Section "Device"
           Identifier  "quadrofx"
           Driver      "nvidia"
           VendorName  "PNY"
           BoardName   "NVIDIA Quadro FX 1100"
           Option "AllowDFPStereo" "true"
           Option "Stereo" "3"
    EndSection
    

``
    
    
    Section "Screen"
           Identifier "Screen0"
           Device     "quadrofx"
           Monitor    "hm204dt"
           DefaultDepth     24
           SubSection "Display"
                   Viewport   0 0
                   Depth     24
                   Modes   "1400x1100_120" "1280x1024"
           EndSubSection
    EndSection
    

```

``

Retrieved from "[https://pymolwiki.org/index.php?title=XFree86_Configuration&oldid=4398](https://pymolwiki.org/index.php?title=XFree86_Configuration&oldid=4398)"


---

