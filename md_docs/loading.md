# loading

loading

### Table of Contents

  * Loading Molecules

    * Formats

    * Methods

  * See Also




# Loading Molecules

## Formats

At the present time, Open-Source PyMOL has support for the following molecular formats: 

  * PDB format with support for various extensions for specifying connectivity, bond valences, multiple models, ensembles, and trajectories. 



  * MOL format for loading of individual small molecule structures.



  * SDF format for loading of multiple small molecule structures.



  * MOL2 format for loading systems containing both small molecules and proteins.



  * XYZ variants, including multimodel XYZ.




PyMOL also has limited support for 

  * ChemDraw3D CC1 and CC2.



  * The old Macromodel MMD/MMOD format.



  * Amber and Gromacs trajectory files, through the VMD plugin architecture.




Incentive PyMOL Builds also support the following proprietary formats: 

  * MOE Molecular Operating Environment (Chemical Computing Group)



  * MAE Maestro (Schrödinger, LLC)




## Methods

Molecular structures can be loaded into PyMOL in the following ways: 

  * via Menus

    * File Menu: Open   





[![:file_open_win.png](/dokuwiki/lib/exe/fetch.php?w=300&tok=01e858&media=loading_file_open.png)](/dokuwiki/lib/exe/detail.php?id=loading&media=loading_file_open.png "loading_file_open.png")

  * Plugin Menu: PDB Loader Service   





[![:plugin_pdb_loader.png](/dokuwiki/lib/exe/fetch.php?w=100&tok=3dd0a8&media=loading_pdb_loader.png)](/dokuwiki/lib/exe/detail.php?id=loading&media=loading_pdb_loader.png "loading_pdb_loader.png")

  * via PyMOL Commands

    * [load](/dokuwiki/doku.php?id=command:load "command:load") imports structures from disk. 
          
          load 1vgc.pdb
          load 1nb6.pdb
          load 1fjf.pdb

    * [fetch](/dokuwiki/doku.php?id=command:fetch "command:fetch") imports structures directly from the [Protein Data Bank](http://www.pdb.org "http://www.pdb.org"). 
          
          fetch 1vgc 1nb6 1fjf

  * via the command shell during program launch:
        
        Unix: pymol 1vgc.pdb 1nb6.pdb 1fjf.pdb
        Windows: "c:\program files\pymol\pymol\pymolwin.exe" +2 1vgc.pdb 1nb6.pdb 1fjf.pdb

  * via Python using [cmd.load](/dokuwiki/doku.php?id=api:cmd:load "api:cmd:load").




# See Also

[fetch](/dokuwiki/doku.php?id=command:fetch "command:fetch") | [load](/dokuwiki/doku.php?id=command:load "command:load")

loading.txt · Last modified: 2013/08/19 21:00 (external edit)
