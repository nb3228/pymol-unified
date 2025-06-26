# installation

installation

### Table of Contents

  * Installation of PyMOL

    * Microsoft Windows

      * Required Hardware

      * Required Software

      * Installation

    * Mac OS X

      * Required Hardware

      * Required Software

      * Installation

    * Linux

      * Required Hardware

      * Required Software

      * Installation

      * Alternate Approach: Compilation From Source




# Installation of PyMOL

The PyMOL installation process is slightly different for each platform. 

First, you must [download the latest build](http://pymol.org/download "http://pymol.org/download") for your platform. 

## Microsoft Windows

### Required Hardware

  * A “wheel” or three-button mouse

  * A graphics card that supports OpenGL (e.g. GeForce, Radeon, Quadro, or FireGL).

  * 512 MB of RAM (1 GB recommended for routine use).




### Required Software

  * PyMOL should run on Windows Vista, Windows 7, 8 and 10.

  * Current installers such as pymol-X_XrX-bin-win32.msi can be downloaded from [here](http://pymol.org/download "http://pymol.org/download"). Note that the X_XX in the preceeding file name corresponds to the PyMOL version number, which changes over time.




### Installation

  1. Download the installer to your hard disk and double-click to launch.

  2. Answer questions asked or click “Next” as appropriate to complete the installation.

  3. Launch the program by selecting “PyMOL” from “All Programs” in the “Start” menu.




## Mac OS X

There are several different ways to run PyMOL on the Macintosh. The following information pertains to [MacPyMOL](/dokuwiki/doku.php?id=macpymol "macpymol"), the native-like port of PyMOL for [Cocoa](http://en.wikipedia.org/wiki/Cocoa_\(API\) "http://en.wikipedia.org/wiki/Cocoa_\(API\)"). 

### Required Hardware

  * A “wheel” mouse or “mighty” mouse configured for three-button usage.

  * 512 MB of RAM (1 GB recommended for routine use).




### Required Software

  * MacPyMOL requires MacOS X 10.9 or later.

  * MacPyMOL-vX_XX.dmg can be downloaded [here](http://pymol.org/dsc/ip "http://pymol.org/dsc/ip"). Note that the X_XX in the preceeding file name corresponds to the PyMOL version number (which changes over time).




### Installation

  1. Double-click the MacPyMOL installer to open (mount) the DMG (disk iamge).

  2. Drag the MacPyMOL icon to the Applications icon. MacPyMOL is now installed.

  3. You can now start MacPyMOL by clicking on the MacPyMOL icon in the Applications folder.




## Linux

### Required Hardware

  * A “wheel” mouse.

  * A graphics card that supports OpenGL from a vendor that provides Linux drivers (e.g. nVidia, ATI, or 3DLabs).

  * An x86-compatible CPU (if not, then you must compile from source).

  * 512 MB of RAM (1 GB recommended for routine use).




### Required Software

  * PyMOL should run with any Linux based on the GNU C library, version 2.3 or greater.

  * PyMOL requires installation of the accelrated OpenGL graphics drivers for your specific graphics hardware.

  * pymol-X_XX-bin-linux-x86-glibc23.tar.bz2 can be downloaded [here](http://pymol.org/dsc/ip "http://pymol.org/dsc/ip"). Note that the X_XX in the preceeding file name corresponds to the PyMOL version number (which changes over time).




### Installation

  * Extract the archive



    
    
    tar -jxf pymol-X_XX-bin-linux-x86-glibc23.tar.bz2

to create a pymol directory. 

  * Run the setup script from within the new directory



    
    
    cd pymol
    ./setup.sh

to create a pymol launch script. 

  * Launch by running the newly created script.



    
    
    ./pymol

### Alternate Approach: Compilation From Source

Information about how to compile the PyMOL from the open source code can be found on the [PyMOL Wiki](http://pymolwiki.org/index.php/Linux_Install "http://pymolwiki.org/index.php/Linux_Install"). 

**Regarding Compilation under Microsoft Windows**

Microsoft Windows, as a proprietary closed-source operating system, is not a supported compilation environment for Open-Source PyMOL. While Windows compilation is of course possible and allowed, it is defined as “beyond scope” for the open-source project. Thus, to use current versions of PyMOL on Windows, you must either sponsor the project in order to access precompiled executables or create, maintain, and support your own port for Windows, preferably via [Cygwin ](http://www.cygwin.com "http://www.cygwin.com"). 

**Regarding Compilation under Mac OS X**

Although Mac OS X is a proprietary closed-source operating system, compilation of PyMOL is supported on Mac OS X under the X11/Fink environment, since that setup is directly compatible with what you'd find on Linux or FreeBSD. Do note however that the resulting executable will not have the integrated single-window user interface available in [MacPyMOL](/dokuwiki/doku.php?id=macpymol "macpymol"). 

installation.txt · Last modified: 2017/09/18 08:50 by holder
  *[MB]: Megabyte
  *[GB]: Gigabyte
  *[ OS]: Operating System
  *[OS]: Operating System
