# Category: Windows

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

