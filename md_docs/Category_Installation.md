# Category: Installation

## Linux Install

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes how to install PyMOL on Linux. 

## Contents

  * 1 Incentive PyMOL
  * 2 Open-Source PyMOL in Linux Distros
  * 3 Install from source
    * 3.1 Requirements
    * 3.2 Get latest source from Git
    * 3.3 Libraries in non-standard places
    * 3.4 Compile and install
    * 3.5 Troubleshooting
  * 4 Customized Installations



## Incentive PyMOL

[Schrödinger](http://www.schrodinger.com) provides pre-compiled (64 bit) PyMOL to paying sponsors. The bundle also includes ready-to-use [APBS](/index.php/APBS "APBS"), [RigiMOL](/index.php/Morph "Morph"), an MPEG encoder for movie export, and a small molecule energy minimization engine. 

Download: <https://pymol.org/>

## Open-Source PyMOL in Linux Distros

Many Linux distributions provide binary packages for open-source PyMOL. They often do not provide the latest version, but if the provided package fits your needs this is the most convenient way to install PyMOL. 

Command line install examples for some popular distributions (note that all of these commands must be run as root or superuser): 
    
    
    # Arch/Manjaro
    pacman -S pymol
    
    # CentOS with EPEL
    rpm -i http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-5.noarch.rpm
    yum --enablerepo=epel install pymol
    
    # Debian/Ubuntu/Mint
    apt-get install pymol
    
    # Fedora
    dnf install pymol
    
    # Gentoo
    emerge -av pymol
    
    # openSUSE (12.1 and later)
    zypper install pymol
    
    # Sabayon
    equo i -av pymol
    

## Install from source

Installation from source gives you the latest version and is the generic way to install PyMOL. 

Please also consult the [INSTALL](https://github.com/schrodinger/pymol-open-source/blob/master/INSTALL) file. 

### Requirements

Libraries as well as development files (headers) of the following software is required: 

  * C++11 compiler (e.g. GCC 4.7 or higher)
  * [Python](http://www.python.org/) 3.6+
  * OpenGL driver (I use [NVidia](http://www.nvidia.com/object/unix.html))
  * GLEW
  * libpng
  * freetype
  * libxml2 (optional, for COLLADA export, disable with `--no-libxml`)
  * [msgpack-c](https://github.com/msgpack/msgpack-c) 2.1.5+ (optional, for fast [MMTF](http://mmtf.rcsb.org/) loading, new in SVN r4167, disable with `--use-msgpackc=no`)
  * PyQt5, PyQt4, or PySide (optional, will fall back to Tk interface if compiled with `--glut`)
  * [glm](https://glm.g-truc.net/)
  * [mmtf-cpp](https://github.com/rcsb/mmtf-cpp) (optional, for [MMTF](http://mmtf.rcsb.org/) export, disable with `--use-msgpackc=no`)
  * libnetcdf (optional, disable with `--no-vmd-plugins`)



Optional/deprecated: 

  * [Pmw](https://github.com/schrodinger/pmw-patched) (Python Megawidgets, for legacy GUI/plugins)
  * GLUT/freeglut (enable with `--glut`)



On many Linux systems, one of the following commands installs all requirements (and must be run as root): 
    
    
    # Debian/Ubuntu/Mint
    apt-get install git build-essential python3-dev libglew-dev \
      libpng-dev libfreetype6-dev libxml2-dev \
      libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev
    
    # CentOS
    yum install gcc gcc-c++ kernel-devel python-devel tkinter python-pmw glew-devel \
      freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel \
      msgpack-devel netcdf-devel
    
    # Fedora
    dnf install gcc gcc-c++ kernel-devel python3-devel glew-devel PyQt5 msgpack-devel \
      freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel
    
    # Gentoo
    emerge -av dev-lang/python dev-python/pmw media-libs/glew \
      media-libs/freeglut media-libs/libpng media-libs/freetype media-libs/glm
    
    # openSUSE
    zypper install python-devel freeglut-devel gcc-c++ glew-devel libpng-devel python-pmw glm
    
    # Sabayon
    equo i -av dev-lang/python dev-python/pmw media-libs/glew \
      media-libs/freeglut media-libs/libpng media-libs/freetype
    

### Get latest source from Git
    
    
    git clone https://github.com/schrodinger/pymol-open-source.git
    git clone https://github.com/rcsb/mmtf-cpp.git
    mv mmtf-cpp/include/mmtf* pymol-open-source/include/
    cd pymol-open-source
    

_The master branch requires Python 3.6+. Use the legacy[py2](https://github.com/schrodinger/pymol-open-source/tree/py2) branch for Python 2.7 compatibility._

### Libraries in non-standard places

Optional: You may use the colon-delimited `$PREFIX_PATH` variable to point `setup.py` to non-standard locations of libraries and headers (those locations should have `include` and `lib` directories). 

### Compile and install

This will install PyMOL as normal user into `$HOME/pymol-open-source-build`. 
    
    
    #!/bin/bash -e
    
    prefix=$HOME/pymol-open-source-build
    
    # Example for dependencies in non-standard places
    # export PREFIX_PATH="$HOME/extra/glew-2.0.0:$HOME/extra/libpng-1.6.5:/opt/local"
    
    python3 setup.py build install \
        --home=$prefix
    

Now launch PyMOL like this: 
    
    
    $HOME/pymol-open-source-build/bin/pymol
    

### Troubleshooting

  * Do do a "clean" build, remove the "build" directory (sometimes necessary if "git pull" changed header files)


  * If you are using Ubuntu with a NVIDIA graphic card and generic drivers you may experience bad rendering, black pixelation and other graphical oddities. A guide to installing NVIDIA proprietary drivers can be found under [Ubuntu community Nvidia Drivers How To](https://help.ubuntu.com/community/BinaryDriverHowto/Nvidia)



## Customized Installations

  * [Troels Linnet's installations scripts](/index.php/User:Tlinnet/Linux_Install "User:Tlinnet/Linux Install"): Detailed installation scripts for Ubuntu, Mint 12 and RHEL 6, including MPEG support from FREEMOL



Retrieved from "[https://pymolwiki.org/index.php?title=Linux_Install&oldid=13303](https://pymolwiki.org/index.php?title=Linux_Install&oldid=13303)"


---

## User:Tlinnet/Linux Install

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

These are customized installation scripts. For more general instructions, see [Linux Install](/index.php/Linux_Install "Linux Install"). 

## Contents

  * 1 Ubuntu/Mint 12 compile and install with MPEG support
  * 2 Red Hat Enterprise Linux RHEL 6 compile and install with MPEG support for x86_64 bit
  * 3 Install script
    * 3.1 Change MPEG settings



## Ubuntu/Mint 12 compile and install with MPEG support

First install dependencies. 
    
    
    sudo apt-get install subversion git
    sudo apt-get install python-imaging python-pygame apbs
    sudo apt-get install build-essential python-dev python-pmw libglew-dev freeglut3-dev libpng-dev libfreetype6-dev
    bash installpymol.sh
    

## Red Hat Enterprise Linux RHEL 6 compile and install with MPEG support for x86_64 bit

Install the EPEL repository: The .rpm will be download from here   
<http://fedoraproject.org/wiki/EPEL>   
And the CentOS repository for additional .rpm are located here   
<http://mirror.centos.org/centos/6/os/x86_64/Packages/>   

    
    
    bash
    cd && wget http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-5.noarch.rpm
    sudo rpm -i epel-release-6-5.noarch.rpm
    rm epel-release-6-5.noarch.rpm</pre>
    

Then add Centos repository 
    
    
    rpm --import http://mirror.centos.org/centos/6/os/x86_64/RPM-GPG-KEY-CentOS-6
    sudo gedit /etc/yum.repos.d/centos.repo
    

Write these lines. 
    
    
    [centos]
    name=Centos for RHEL/ CentOS $releasever - $basearch
    baseurl=http://mirror.centos.org/centos/6/os/x86_64
    enabled=1
    

Then check configuration. Do NOT UPGRADE PACKAGES from CentOS. It will replace RHEL packages.   
If you do upgrade packages, you will now see Centos login and such, but not much is problematic. 
    
    
    sudo yum repolist all
    sudo yum install subversion git
    sudo yum install python-imaging pygame apbs
    sudo yum install gcc-c++ python-devel python-pmw glew-devel freeglut-devel libpng-devel freetype-devel
    

Then disable centos repository. Set **enabled=0**
    
    
    sudo gedit /etc/yum.repos.d/centos.repo
    

Now install pymol 
    
    
    bash installpymol.sh
    

Make sure, that you will have **~/bin** as part of your path 
    
    
    gedit ~/.cshrc
    

Make sure this lines is in the file. Only if you use tcsh shell. 
    
    
    setenv PATH ${PATH}:{$HOME}/bin
    

## Install script

Make a text file "installpymol.sh" and make it executable 
    
    
    chmod u+x installpymol.sh
    

Put this in the file, modify the first line 
    
    
    #!/bin/bash -e
    prefix=$PWD/pymol
    modules=$prefix/modules
    svnpymol=svnpymol
    svnfreemol=svnfreemol
    pymolscriptrepo=Pymol-script-repo
    update=$prefix/updatepymol.sh
    
    ###################################################
    mkdir -p $prefix
    mkdir -p $HOME/bin
     
    ###### Checkout pymol svn
    svn co svn://svn.code.sf.net/p/pymol/code/trunk/pymol $prefix/$svnpymol
    ###### Build and install pymol
    cd $prefix/$svnpymol
    python setup.py build install --home=$prefix --install-lib=$modules --install-scripts=$prefix
     
    ########## Setup freemol - for MPEG support ############
    svn co svn://bioinformatics.org/svnroot/freemol/trunk $prefix/$svnfreemol
    cd $prefix/$svnfreemol/src/mpeg_encode
    export FREEMOL=$prefix/$svnfreemol/freemol
    ./configure
    make
    make install
     
    ########## Install Pymol-script-repo ############
    git clone git://github.com/Pymol-Scripts/Pymol-script-repo.git $prefix/$pymolscriptrepo
     
    ## Make a shortcut to an extended pymol execution
    echo "#!/bin/bash" > $prefix/pymolMPEG.sh
    echo "export FREEMOL=$prefix/$svnfreemol/freemol" >> $prefix/pymolMPEG.sh
    echo "export PYMOL_GIT_MOD=$prefix/$pymolscriptrepo/modules" >> $prefix/pymolMPEG.sh
    echo "export PYTHONPATH=$prefix/$pymolscriptrepo/modules"':$PYTHONPATH' >> $prefix/pymolMPEG.sh
    echo "export PYTHONPATH=$prefix/$pymolscriptrepo"':$PYTHONPATH' >> $prefix/pymolMPEG.sh
    echo '#export PYTHONPATH=$PYTHONPATH:/sbinlab2/software/x64/lib64/python2.6/site-packages/PIL' >> $prefix/pymolMPEG.sh
    echo '#export PYTHONPATH=$PYTHONPATH:/sbinlab2/software/x64/lib64/python2.6/site-packages/lib-dynload' >> $prefix/pymolMPEG.sh
    echo '#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sbinlab2/software/x64/lib/pymollib' >> $prefix/pymolMPEG.sh
    echo '#export LIBGL_ALWAYS_INDIRECT=no' >> $prefix/pymolMPEG.sh
    tail -n +2 $prefix/pymol >> $prefix/pymolMPEG.sh
    chmod ugo+x $prefix/pymolMPEG.sh 
     
    ## Make a startup files, which is always executed on startup.
    t="'"
    echo "import sys,os" > $modules/pymol/pymol_path/run_on_startup.py
    echo "import pymol.plugins" >> $modules/pymol/pymol_path/run_on_startup.py
    echo "pymol.plugins.preferences = {'instantsave': False, 'verbose': False}" >> $modules/pymol/pymol_path/run_on_startup.py
    echo "pymol.plugins.autoload = {'apbs_tools': False}" >> $modules/pymol/pymol_path/run_on_startup.py
    echo "pymol.plugins.set_startup_path( [$t$prefix/$pymolscriptrepo/plugins$t,$t$modules/pmg_tk/startup$t] )" >> $modules/pymol/pymol_path/run_on_startup.py
    echo "pymol.plugins.preferences = {'instantsave': True, 'verbose': False}" >> $modules/pymol/pymol_path/run_on_startup.py
    
    ## Make a update script
    cat > "$update" <<EOF
    #!/bin/bash -e
    prefix=$prefix
    modules=$modules
    svnpymol=svnpymol
    svnfreemol=svnfreemol
    pymolscriptrepo=Pymol-script-repo
    
    ###### Checkout pymol svn
    svn co svn://svn.code.sf.net/p/pymol/code/trunk/pymol \$prefix/\$svnpymol
    ###### Build and install pymol
    cd \$prefix/\$svnpymol
    python setup.py build install --home=\$prefix --install-lib=\$modules --install-scripts=\$prefix
    
    ########## Setup freemol - for MPEG support ############
    svn co svn://bioinformatics.org/svnroot/freemol/trunk \$prefix/\$svnfreemol
    cd \$prefix/$svnfreemol/src/mpeg_encode
    export FREEMOL=\$prefix/\$svnfreemol/freemol
    ./configure
    make
    make install
    
    ########## Update Pymol-script-repo ############
    cd \$prefix/\$pymolscriptrepo
    git pull origin master
    cd \$prefix
    EOF
    chmod +x $update
    

### Change MPEG settings

Change settings in 
    
    
    $HOME/Software/pymol/svnfreemol/freemol/libpy/freemol/mpeg_encode.py
    

For example, change in line 205:  
FRAME_RATE 24  
(Note, only legal values is allowed: 23.976, 24, 25, 29.97, 30, 50 ,59.94, 60) 

Then restart PyMOL. 

Retrieved from "[https://pymolwiki.org/index.php?title=User:Tlinnet/Linux_Install&oldid=11163](https://pymolwiki.org/index.php?title=User:Tlinnet/Linux_Install&oldid=11163)"


---

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

