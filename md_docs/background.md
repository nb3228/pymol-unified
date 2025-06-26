# background

background

### Table of Contents

  * Project Background Information

    * Motivation

    * Words of Caution

    * Strengths

    * Weaknesses




# Project Background Information

## Motivation

Why Open-Source PyMOL? 

PyMOL started as one scientist's answer to frustration encountered with then-existing visualization and modeling software while practicing computational chemistry at a small pharmaceutical company. 

All who have studied the remarkable complexity of a protein structure will agree that visualization is essential to understanding structural biology. Nevertheless, most researchers who use visualization packages ultimately run up against limitations inherent in them which make it difficult or impossible to get exactly what you need. Such limitations in a closed-source commercial software package cannot be easily surmounted, and the same is still true for free programs which aren't available in source-code form. 

Only open-source software allows you to surmount problems by directly changing and enhancing the way software operates, and it places virtually no restrictions on your power and opportunity to innovate. For these reasons, we believe that open-source software is an intrinsically superior research product and will provide greatest benefit to computer-assisted scientific research over the long term. 

Launched in December 1999, PyMOL was originally designed to: 

  1. visualize multiple conformations of a single structure [trajectories or docked ligand ensembles] 

  2. interface with external programs

  3. provide professional strength graphics under both Windows and Unix,

  4. prepare publication quality images, and 

  5. fit into a tight budget. 




By 2006, all of these goals have long since been realized, and PyMOL now includes many other useful capabilities. However, PyMOL also has many weaknesses, and it will undoubtedly take years of continued effort to surmount them. Nevertheless, we hope that you will find PyMOL to be a valuable tool for your work today, and we encourage you to provide any suggestions you have for making it even better. 

## Words of Caution

PyMOL was created in an efficient but highly pragmatic manner, with heavy emphasis on delivering specific features to end users. Expediency has almost always taken precedence over elegance, and adherence to established software development practices is inconsistent at best. To date, PyMOL has been about getting specific jobs done as fast as possible and by whatever means are available. Though PyMOL can some important needs today, it is just the first generation of a multi-generational software development effort, where increased attention is paid to good software engineering practices in subsequent iterations, as resources permit. 

Though PyMOL has succeeded in becoming a community-based open-source tool, it has so far failed to become a true multi-developer project due to complex intricacy of the source code. As a result, open-source developers should view PyMOL a semi-opaque tool, best extended through indepenent Python modules, rather than as a C application that could be improved upon or extended directly. 

## Strengths

* **Cross-Platform.** A single code base supports both Unix, Macintosh, and Windows, using OpenGL and Python and a small set of Open-source external dependencies. 

* **Command-Line and GUI Control.** Real world applications require both. 

* **Atom Selections.** Arbitrary logical expressions enable focused visualization and editing. 

* **Molecular Splits/Joins.** Structures can be sliced, diced, and reassembled on the fly and written out to standard files (i.e. PDB). 

* **Movies.** Creating movies can be as simple as loading multiple PDB files and hitting the play button. 

* **Surfaces.** As good if not better than Grasp, and mesh surfaces are supported too. 

* **Cartoon Ribbons.** PyMOL's cartoons are almost as nice as Molscript but are much easier to create and render. 

* **Scripting.** The best way to control PyMOL is through reusable scripts, which can be written in the PyMOL command language or in Python. 

* **Rendering.** A built-in ray tracer gives you shadows and depth on any scene. You also render externally using PovRay or VRML. 

* **Output.** PNG files output from PyMOL can be directly imported into PowerPoint. 

* **Conformational Editing.** Click and drag interface allows you to edit conformations naturally. 

* **Sculpting.** Allows the molecule to adapt to your changes. 

* **Expandability.** The PyMOL Python API provides a solid way to extend and interface. 

## Weaknesses

* **User Interface.** Development has been focused on capabilities, not on easy-of-use for new users or non-specialists who are averse to using scripts. 

* **Documentation.** Only recently has any documentation become available. 

* **Object-Orientation.** There is a single monolithic, functional API. 

* **Electrostatics.** PyMOL is not yet a replacement for Delphi/Grasp. 

* **No Mechanics Engine** Although PyMOL sports potent molecular editing features, you cannot yet perform a molecular “clean-up” without use of external tools. 

background.txt · Last modified: 2013/08/19 21:00 (external edit)
  *[ GUI]: Graphical User Interface
  *[API]: Application Programming Interface
