# Category: Launching

## Command Line Options

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL has lots of options for controlling it immediately from startup. Run '**help launching'** from the command line for detailed listings of options for your current version of PyMol. Here is a sample output which should be pretty consistent for recent versions. 

# Command Line Options
    
    
    shell> pymol --help
    
    
    
    PyMOL 1.8.6.0 Incentive Product
    Copyright (C) Schrodinger, LLC
    
    Usage: pymol [OPTIONS]... [FILES]... [-- CUSTOM SCRIPT ARGUMENTS]
    
    Options
    
      --help    display this help and exit
      --version display PyMOL version and exit
      --retina  use retina resolution (MacPyMOL) and set display_scale_factor=2
    
      -1        config_mouse one_button
      -2        config_mouse two_button
      -a N      alias for -A
      -A N      application configuration:
        -A1     simple viewer window          (-qxiF -X 68 -Y 100)
        -A3     internal GUI only, no splash  (-qx -X 68 -Y 100)
        -A4     used by PYMOLVIEWER           (-X 68 -Y 100)
        -A5     helper application            (-QxiICUF -X 68 -Y 100)
        -A6     full screen presentation      (-qxieICUPF)
      -b[N]     benchmark wizard
      -B        (DEPRECATED)
      -c        launch in command-line only mode for batch processing
      -C        don't terminate on Ctrl-C
      -d cmd    execute PyMOL command
      -D N      defer_builds_mode=N
      -e        full screen
      -E N      multisampling (GL_MULTISAMPLE_ARB)
      -f N      internal_feedback=N
      -F        internal_feedback=0
      -g file   save image (png) or movie (mpg)
      -G        game mode (DEPRECATED)
      -h        generic helper application (no controls, no feedback)
      -H N      window height in pixels
      -i        internal_gui=0
      -I        auto_reinitialize=1 (Mac only)
      -j        side-by-side stereo (stereo_mode=4)
      -J        cd to user's home directory
      -k        don't load pymolrc or plugins
      -K        keep alive: when running without a GUI, don't quit after the input
                is exhausted
      -l file   run python script in thread (spawn)
      -L file   load file after everything else (only if something was loaded before)
      -m        INTERNAL - do not use (mac external GUI)
      -M        force mono
      -n        INTERNAL - do not use (incentive_product=1)
      -N name   UNSUPPORTED - external gui type (default: pmg_tk) (same as -w)
      -o        disable security protections
      -O N      sphere_mode=N
      -p        read commands from STDIN
      -P        handle scenes as if the session were opened in presentation mode
      -q        supress startup message
      -Q        quiet, suppress all text output
      -r file   run python script
      -R        launch RPC Server
      -s file   log to file
      -S        force stereo
      -t N      stereo_mode=N
      -T name   UNSUPPORTED - Tcl/Tk GUI skin
      -u file   resume log file (execute existing content and append new log output)
      -U        UNSUPPORTED reuse the helper application
      -v        UNUSED
      -V N      external GUI window height in pixels
      -w name   UNSUPPORTED - external gui type (default: pmg_tk) (same as -N)
      -W N      window width in pixels
      -x        no external gui
      -X N      window x position on screen
      -y        exit on error
      -Y N      window y position on screen
      -z N      window_visible=N
      -Z N      zoom_mode=N
    

# Supported File Formats

<file> can have one of the following extensions, and all files provided will be loaded or run after PyMOL starts. 

PyMOL Supported File Formats  Extension | Format   
---|---  
**.pml** | PyMOL command script to be run on startup   
**.py, .pym, .pyc** | Python program to be run on startup   
**.pdb** | Protein Data Bank format file to be loaded on startup   
**.mmod** | Macromodel format to be loaded on startup   
**.mol** | MDL MOL file to be loaded on startup   
**.sdf** | MDL SD file to be parsed and loaded on startup   
**.xplor** | X-PLOR Map file (ASCII) to be loaded on startup   
**.ccp4** | CCP4 map file (BINARY) to be loaded on startup   
**.cc1, .cc2** | ChemDraw 3D cartesian coordinate file   
**.pkl** | Pickled ChemPy Model (class "chempy.model.Indexed")   
**.r3d** | Raster3D file   
**.cex** | CEX file (Metaphorics)   
**.top** | AMBER topology file   
**.crd** | AMBER coordinate file   
**.rst** | AMBER restart file   
**.trj** | AMBER trajectory   
**.pse** | PyMOL session file   
**.phi** | Delphi/Grasp Electrostatic Potential Map   
  
# Examples
    
    
    # load the 3 pdb files and execte the following commands
    pymol mol1.pdb mol2.pdb mol3.pdb -d 'as ribbon;spectrum count;set seq_view'
    
    # load the 3 pdbs and run the script
    pymol mol1.pdb mol2.pdb mol3.pdb my_script.pml
    
    # Run my_program using PyMOL as the Python interpreter
    # passing the 3 pdb files to my_program.py
    pymol my_program.py -- mol1.pdb mol2.pdb mol3.pdb
    

Retrieved from "[https://pymolwiki.org/index.php?title=Command_Line_Options&oldid=12570](https://pymolwiki.org/index.php?title=Command_Line_Options&oldid=12570)"


---

## Launching From a Script

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The recommended way to run PyMOL-Python scripts is by using PyMOL as the interpreter. This is supported by all versions of PyMOL, including the pre-compiled bundles provided by Schrödinger. 

Example from a shell: 
    
    
    shell> pymol -cq script.py
    

With arguments (_sys.argv_ becomes _["script.py", "foo", "bar"]_): 
    
    
    shell> pymol -cq script.py -- foo bar
    

Example from a running PyMOL instance: 
    
    
    PyMOL> run script.py
    

For advanced users, the following PyMOL versions also allow to run PyMOL from an existing Python process: 

  * [PyMOL 2.0](https://pymol.org/2/#download) based on Anaconda (using Anaconda's python, which is included in bundles provided by Schrödinger)
  * Open-Source PyMOL
  * Schrödinger-provided "[Mac alternative X11-only build](http://pymol.org/download)" of the 1.8.x series



After importing the **pymol** module, PyMOL's event loop has to be started with a call to **pymol.finish_launching()** _(not supported on macOS)_. 

## Contents

  * 1 Usage
  * 2 Example 1
  * 3 Example 2
  * 4 STDOUT
  * 5 Independent PyMOL Instances
  * 6 Independent PyMOL Instances (Context Manager)
  * 7 See Also



## Usage

With PyMOL 2.1, calling any pymol.cmd function will automatically start a backend process without the GUI in the main thread. "finish_launching" should not be necessary, and will launch PyMOL in a new thread with an event loop, which will cause 100% CPU usage (at least with "-c"). 
    
    
    from pymol import cmd
    cmd.fragment('ala')
    cmd.zoom()
    cmd.png('/tmp/test.png', 300, 200)
    

Since PyMOL 1.7.4, the following form is sufficient: 

_This is not supported on macOS (see[bug report](https://github.com/schrodinger/pymol-open-source/issues/28))_
    
    
    import pymol
    pymol.finish_launching(['pymol', '-q'])
    

Before 1.7.4, either "-K" was needed as an additional argument, or arguments had to be assigned to __main__.pymol_argv or pymol.pymol_argv. 

## Example 1

Here is an example script that launches PyMol for stereo viewing on a [VisBox](http://www.visbox.com/boxMain.html). It runs PyMol fullscreen stereo, and disables the internal gui. The environment (PYTHON_PATH and PYMOL_PATH) must already be set up for this example to work (see Example 2 below for how to setup within the script). 
    
    
    #!/usr/bin/env python
     
    # Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
    import __main__
    __main__.pymol_argv = [ 'pymol', '-qei' ]
     
    import pymol
     
    # Call the function below before using any PyMOL modules.
    pymol.finish_launching()  # not supported on macOS
     
    from pymol import cmd
    cmd.stereo('walleye')
    cmd.set('stereo_shift', 0.23)
    cmd.set('stereo_angle', 1.0)
    

## Example 2

This script launches PyMOL without any GUI for scripting only. It enables tab-completion on the python command line and does the PyMOL environment setup (you need to adjust the **moddir** variable!). _Hint: You may save this as "pymol-cli" executable._
    
    
    #!/usr/bin/python2.6 -i
    
    import sys, os
    
    # autocompletion
    import readline
    import rlcompleter
    readline.parse_and_bind('tab: complete')
    
    # pymol environment
    moddir='/opt/pymol-svn/modules'
    sys.path.insert(0, moddir)
    os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')
    
    # pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
    import pymol
    pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
    pymol.finish_launching()
    cmd = pymol.cmd
    

## STDOUT

PyMOL captures **sys.stdout** and **sys.stderr** , to control it with it's own [feedback](/index.php/Feedback "Feedback") mechanism. To prevent that, save and restore both streams, e.g.: 
    
    
    import sys
    stdout = sys.stdout
    stderr = sys.stderr
    pymol.finish_launching(['pymol', '-xiq'])  # not supported on macOS
    sys.stdout = stdout
    sys.stderr = stderr
    

## Independent PyMOL Instances

It's possible to have multiple independent instances. 
    
    
    import pymol2
    p1 = pymol2.PyMOL()
    p1.start()
    
    p2 = pymol2.PyMOL()
    p2.start()
    
    p1.cmd.fragment('ala')
    p1.cmd.zoom()
    p1.cmd.png('/tmp/ala.png', 1000, 800, dpi=150, ray=1)
    
    p2.cmd.fragment('ser')
    p2.cmd.zoom()
    p2.cmd.png('/tmp/ser.png', 1000, 800, dpi=150, ray=1)
    
    p1.stop()
    p2.stop()
    

## Independent PyMOL Instances (Context Manager)

PyMOL 2.2 adds context manager support. 
    
    
    import pymol2
    with pymol2.PyMOL() as p1:
        p1.cmd.fragment('ala')
        p1.cmd.zoom()
        p1.cmd.png('/tmp/ala.png', 1000, 800, dpi=150, ray=1)
    

## See Also

  * [Command Line Options](/index.php/Command_Line_Options "Command Line Options")
  * [Jupyter](/index.php/Jupyter "Jupyter")



Retrieved from "[https://pymolwiki.org/index.php?title=Launching_From_a_Script&oldid=13294](https://pymolwiki.org/index.php?title=Launching_From_a_Script&oldid=13294)"


---

## Launching PyMOL

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL can be started from the command line or from a Python script. Also, PyMOL has lots of [command line options](/index.php/Command_Line_Options "Command Line Options") that one can pass in to affect PyMOL's behavior. For example, one can load a few proteins and even execute a script all from the commandline. Also, there's a [ .pymolrc](/index.php/Pymolrc "Pymolrc") file that initializes PyMOL each time PyMOL is executed. 

  * More on launching PyMOL 
    * [Command Line Options](/index.php/Command_Line_Options "Command Line Options")
    * [Launching From a Script](/index.php/Launching_From_a_Script "Launching From a Script")



## Contents

  * 1 Invoking PyMOL and reading startup commands from a file
    * 1.1 Linux
    * 1.2 Windows
    * 1.3 MacOS X
      * 1.3.1 Launching
      * 1.3.2 IDLE
      * 1.3.3 Reading the pymolrc file
    * 1.4 Launching PyMOL from an external application
    * 1.5 Running PyMOL in batch mode
    * 1.6 Suppressing PyMOL output
    * 1.7 Launching Python scripts
    * 1.8 Overwriting Default Settings



# Invoking PyMOL and reading startup commands from a file

## Linux

Assuming the executable is in your **$PATH** , simply issue 
    
    
    pymol
    

together with any [Command_Line_Options](/index.php/Command_Line_Options "Command Line Options") and arguments (pdb files, pse files, map files and so forth) you require. 

Whenever PyMol starts, a user-created '~/.pymolrc' file containing commands is run. All you need to do is create ".pymolrc" and place it in your home directory. Alternatively, you can instead create ".pymolrc.py" which contains actual Python code instead of just PyMOL commands. See an [example .pymolrc](/index.php/Inchoates_pymolrc "Inchoates pymolrc")). 

## Windows

On Windows, use 'pymolrc', 'pymolrc.py' or 'pymolrc.pym'. For global defaults (all users), you can place a .pymolrc file in C:\Program Files\DeLano Scientific\PyMOL (or C:\Program Files\PyMOL\PyMOL). You can launch PyMOL from the applications menu or from the icon on your desktop (if you placed one there). 

Alternatively, if you want to change the default directory for data, you can just put the path to the new default directory in the "Start In" field of the "shortcut link" that is used to launch Pymol. To do this, right click on the "shortcut link" that is used to launch Pymol, click on "Properties" and then, in the "Start In" field, type in the path to the new default data directory. Then click "OK" and the changes are complete. 

## MacOS X

### Launching

  * Double-click the application's icon
  * Issue the unix command 
        
        open -a MacPyMOL
        

  * Directly invoking the unix executable 
        
        /Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL
        

  * If MacPyMOL is set as the default application to open pdb, pse and other such files, double-clicking any of those files (or issuing the unix open command) will start also pymol and load the file you clicked on.



  
The third one assumes the application has been placed in **/Applications** , so adjust the absolute path if you have it elsewhere. If you directly invoke the unix executable to launch pymol, it has the advantage that you can pass [Command_Line_Options](/index.php/Command_Line_Options "Command Line Options") and arguments to it in the usual way. You might wish to **make an alias** : 
    
    
    alias pymol=/Applications/Xtal/MacPyMOL.app/Contents/MacOS/MacPyMOL
    

(leave out the equal sign for tcsh) **or[symbolic link](http://en.wikipedia.org/wiki/Symbolic_link)**
    
    
    sudo ln -s /Applications/Xtal/MacPyMOL.app/Contents/MacOS/MacPyMOL /usr/local/bin/pymol
    

or **use this shell script** [I](/index.php/User:Wgscott "User:Wgscott") wrote: [pymol](http://xanana.ucsc.edu/Library/init/zsh/local-functions/xtal/pymol) shell script (and zsh function) to invoke pymol on the command line. It uses OS X 10.4's mdfind to locate the executable. 

### IDLE

MacPyMOL can also be run from within [IDLE](http://www.python.org/idle/doc/idle2.html). Make sure you use Fink's IDLE, not Mac IDLE: 
    
    
    which idle2.5     # should return the following
    /sw/bin/idle2.5   # Fink's IDLE.
    

### Reading the pymolrc file

In each case, PyMol will read the contents of the user's ~/.pymolrc file and/or ~/.pymolrc.py file (as with Linux). Here is one example of a [MacOSX-specific .pymolrc file](/index.php/MacOSX-specific_.pymolrc_file "MacOSX-specific .pymolrc file") with a script enabling interaction via the [PowerMate Dial](/index.php/MAC_Install#PowerMate_Dial "MAC Install"). A couple of simple lines to put in your .pymolrc might be to respectively change to your favorite directory and increase the window size: 
    
    
    cd ~/Documents/structures/
    viewport 750,750
    

## Launching PyMOL from an external application

If **PYMOL_PATH** , **LD_LIBRARY_PATH** , and **TCL_LIBRARY** are correctly defined, then you can launch PyMOL from an external Python program as shown in _examples/devel/start_pymol.py_. **nb.** : This approach is not recommended, since the PyMOL launching process is subject to change without warning. The recommended approach is to just use PyMOL as your python interpreter: 
    
    
    pymol -r <script.py>
    pymol -qcr <script.py>
    

## Running PyMOL in batch mode

To perform PyMOL commands from stdin (file, pipe) without opening an OpenGL window use the following [Command_Line_Options](/index.php/Command_Line_Options "Command Line Options"): 
    
    
    pymol -cq
    

## Suppressing PyMOL output

To suppress most of PyMOL's normal chatter, just type on the PyMOL command line: 
    
    
    feedback disable,all,actions
    feedback disable,all,results
    

or, from Python (API): 
    
    
    cmd.feedback("disable","all","actions")
    cmd.feedback("disable","all","results")
    

## Launching Python scripts

Running a Python script from PyMOL, usually the command: 
    
    
    run script.py
    

Is enough. Of course, the file script.py needs to be in the working directory. For more detailed examples, see the commands to launch Python scripts when starting PyMOL. Asynchronous means, that a new Python thread is started: 
    
    
    pymol example.py     # synchronous, in PyMOL module
    pymol -r example.py  # synchronous in __main__ module
    pymol -l example.py  # asychronous in a new module
    

You can also launch python programs from within PyMOL with the commands: 
    
    
    run example.py        # synchronous in pymol module
    run example.py,main   # synchronous in __main__ module
    
    spawn example.py        # asychronous in a new module
    spawn example.py,global # asychronous in the PyMOL module
    spawn example.py,main   # asychronous in the __main__ module
    

## Overwriting Default Settings

If you don't like the cartoon default color and want to change default settings for it, from green to slate, you add the command line to do this in **$HOME/.pymolrc**. Here, 'set cartoon_color, slate' 

Windows users can do this with a pymolrc file. 

Retrieved from "[https://pymolwiki.org/index.php?title=Launching_PyMOL&oldid=9017](https://pymolwiki.org/index.php?title=Launching_PyMOL&oldid=9017)"


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

## Plugin Manager

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Since version 1.5.0.5 PyMOL comes with a **Plugin Manager** , which can be used to load plugins such as those in the [PyMOL Script Repo](https://github.com/Pymol-Scripts/Pymol-script-repo#).  


## Contents

  * 1 Features
    * 1.1 Appending additional paths to the Plugin Manager
  * 2 Screenshots
  * 3 See Also



## Features

  * Install/Uninstall plugins
  * Disable plugins and load them on demand to optimize PyMOL's startup time
  * Configure the plugin search path



### Appending additional paths to the Plugin Manager

Should your scripts be located in several alternative locations, it is possible to append additional directories to the **Plugin Manager**.  


    Plugin > Plugin Manager > Settings > Add new directory...

## Screenshots

[![Plugin manager installed.png](/images/3/3d/Plugin_manager_installed.png)](/index.php/File:Plugin_manager_installed.png) [![Plugin manager install new.png](/images/1/1b/Plugin_manager_install_new.png)](/index.php/File:Plugin_manager_install_new.png) [![Plugin manager settings.png](/images/7/7b/Plugin_manager_settings.png)](/index.php/File:Plugin_manager_settings.png)

## See Also

  * [Pymol-script-repo](/index.php/Pymol-script-repo "Pymol-script-repo")



Retrieved from "[https://pymolwiki.org/index.php?title=Plugin_Manager&oldid=12823](https://pymolwiki.org/index.php?title=Plugin_Manager&oldid=12823)"


---

## Pymolrc

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

When Pymol [launches](/index.php/Launching_PyMOL "Launching PyMOL"), it will read custom settings and functions from a pymolrc file, if it exists. This is simply a script in PyMOL command syntax or in python syntax, depending on the suffix.  
Creating pymolrc files is a convenient way of loading individualized settings. 

## Contents

  * 1 Open a pymolrc File for Editing
  * 2 Technical Details
    * 2.1 Multiple pymolrc files
  * 3 Example
  * 4 Reload a pymolrc File
  * 5 See Also



## Open a pymolrc File for Editing

PyMOL 2.0 has a built-in editor: **File > Edit pymolrc**

On Windows: **Start > Run** and then paste 
    
    
    notepad "%HOMEDRIVE%%HOMEPATH%\pymolrc.pml"
    

On Unix/Linux-type system (including Mac OS X): Open a terminal and type 
    
    
    nano ~/.pymolrc
    

## Technical Details

  * The leading dot is optional, so everything that starts with **.pymolrc** or **pymolrc** will be found by PyMOL. Files with dot take precedence over files without dot.
  * Files ending on **.py** (or .pym) will be parsed as python scripts, files ending on **.pml** or without suffix will be parsed as PyMOL command syntax.
  * Several directories are searched, in order: 
    * current working directory
    * $HOME
    * $HOMEDRIVE + $HOMEPATH (on Windows)
    * $PYMOL_PATH
  * $PYMOL_PATH/run_on_startup* is loaded before any **pymolrc** files



### Multiple pymolrc files

PyMOL will even load **multiple pymolrc files** , however only either with dot or without dot, and only from the same directory.  
You can have multiple scripts: e.g. **pymolrc-settings.pml** and **pymolrc-misc.pml** in your home directory. 

  * **pymolrc-settings.pml** can e.g. be used to define 'permanent' custom [Settings](/index.php/Settings "Settings") that you rarely change
  * **pymolrc-misc.pml** can e.g. be used to define more transient custom [Settings](/index.php/Settings "Settings"), such as [Working Directory](/index.php/Cd "Cd") or [Fetch Path](/index.php/Fetch_Path "Fetch Path")



You can query which pymolrc files have been loaded: 
    
    
    PyMOL> print invocation.options.deferred
    

## Example
    
    
    # simple test: change background color of PyMOL window
    bg blue
    
    # this will run the script in the specified location
    run /path/to/home/pymol/load_sep.py
    
    # your favorite settings
    set movie_loop, 0
    set two_sided_lighting, 1
    
    set label_size, 60
    set label_outline_color, 1
    set label_color, 0
    set label_position, [0, 0, 10]
    
    # for images:
    #   antialias =1 smooths jagged edges, 0 turns it off
    set antialias, 1
    
    #   stick_radius -adjust thickness of atomic bonds
    set stick_radius, 0.3
    
    # save fetched PDB files here
    set fetch_path, /your/fetch/path
    
    # Personal short-cut to color_obj function
    import color_obj
    cmd.extend("co",color_obj.color_obj)
    

## Reload a pymolrc File

To reload a pymolrc file (e.g. after editing .pymolrc, or after running [reinitialize](/index.php/Reinitialize "Reinitialize")), [run it](/index.php/Running_Scripts "Running Scripts") like any other script: 
    
    
    @~/.pymolrc
    

or 
    
    
    run ~/.pymolrc.py
    

## See Also

  * [save_settings](/index.php/Save_settings "Save settings")
  * [Launching PyMOL](/index.php/Launching_PyMOL "Launching PyMOL")



Retrieved from "[https://pymolwiki.org/index.php?title=Pymolrc&oldid=12769](https://pymolwiki.org/index.php?title=Pymolrc&oldid=12769)"


---

## Remote Desktop

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes solutions to run PyMOL on a remote computer. 

## VirtualGL

[VirtualGL](https://www.virtualgl.org/) must be installed on the client (vglconnect) and the server (vglrun). It renders on the GPU of the remote machine. 
    
    
    localcomputer $ vglconnect remotecomputer
    remotecomputer $ vglrun pymol
    

## X Forwarding

SSH X forwarding requires an [X server](https://en.wikipedia.org/wiki/X_Window_System) on the client. It uses indirect rendering (GLX) and renders on the GPU of the client machine. This is limited to [use_shaders=off](/index.php?title=Use_shaders&action=edit&redlink=1 "Use shaders \(page does not exist\)") and can be quite slow. 
    
    
    localcomputer $ ssh -Y remotecomputer
    remotecomputer $ pymol
    

On macOS, since [XQuartz 2.7.10](https://www.xquartz.org/releases/XQuartz-2.7.10.html), IGLX has to be enabled first: 
    
    
    defaults write org.macosforge.xquartz.X11 enable_iglx -bool true
    killall Xquartz
    

## See Also

  * [RPC](/index.php/RPC "RPC")



Retrieved from "[https://pymolwiki.org/index.php?title=Remote_Desktop&oldid=12835](https://pymolwiki.org/index.php?title=Remote_Desktop&oldid=12835)"


---

## RPC

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Remote control of PyMOL is possible with [XML-RPC](http://en.wikipedia.org/wiki/XML-RPC). 

Launch PyMOL with the "-R" option to start the RPC server. 
    
    
    pymol -R
    

Connect from any client, using python: 
    
    
    try:
        # Python 2
        import xmlrpclib
    except ImportError:
        # Python 3:
        import xmlrpc.client as xmlrpclib
    
    srv = xmlrpclib.ServerProxy('http://localhost:9123')
    
    srv.do('fetch 2xwu')
    srv.do('as cartoon')
    srv.do('spectrum')
    
    # or
    srv.fetch('2xwu')
    srv.show_as('cartoon')
    srv.spectrum()
    

## See Also

  * [Launching From a Script](/index.php/Launching_From_a_Script "Launching From a Script")
  * [Launching PyMOL](/index.php/Launching_PyMOL "Launching PyMOL")



Retrieved from "[https://pymolwiki.org/index.php?title=RPC&oldid=12797](https://pymolwiki.org/index.php?title=RPC&oldid=12797)"


---

