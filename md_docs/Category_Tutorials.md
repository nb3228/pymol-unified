# Category: Tutorials

## Advanced Scripting

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

On this page, we discuss more complex scripting. Python is great, but it is much slower at mathematics than C/C++/Java/FORTRAN. For that reason, you may find it more useful to export your data to another language, operate on it there and then import the results back into PyMOL. We discuss the Python API and the general operating procedure for successfully writing your own scripts. 

## Contents

  * 1 Advanced Scripting
    * 1.1 Python, PyMOL and C
      * 1.1.1 C++
      * 1.1.2 Calling the External Function
      * 1.1.3 Unpacking the Data
      * 1.1.4 Reference Counting
      * 1.1.5 More Complex Unpacking
    * 1.2 Sending the Results back to Python/PyMOL
    * 1.3 Initialization
    * 1.4 Installing Your Module
      * 1.4.1 Overview
      * 1.4.2 Setup.py
    * 1.5 Notes
    * 1.6 Conclusion
    * 1.7 Example
    * 1.8 See Also



# Advanced Scripting

Python while incredibly useful, is much slower at math than some other strictly typed languages and sometimes we have libraries built in other languages. It's faster, for complicated problems to package your data, send it to C, do some math, and pass the results back to Python than to just do everything in Python. The beauty of the Python API, is that we can do just that. 

This is more advanced scripting, and requires some knowledge of the [Python API](https://docs.python.org/2/extending/extending.html), and some outside language. The example shown here is in C. The C++ extensions are very similar. 

  


### Python, PyMOL and C

Here, I will show you how to write a C-module that plugs into Python and talks nicely with PyMOL. The example actually shows how to make a generic C-function and use it in Python. 

First, let's assume that we want to call a function, let's call it **funName**. Let's assume **funName** will take a Python list of lists and return a list---for example passing the C++ program the XYZ coordinates of each atom, and returning a list of certain atoms with some property. I will also assume we have **funName.h** and **funName.c** for C code files. I have provided this, a more complex example, to show a real-world problem. If you were just sending an integer or float instead of packaged lists, the code is simpler; if you understand unpacking the lists then you'll certainly understand unpacking a simple scalar. 

##### C++

If you tell Python that you're using C++ code (see the setup below) then it'll automatically call the C++ compiler instead of the C compiler. There are [warnings](http://docs.python.org/ext/cplusplus.html) you may want to be aware of though. 

My experience with this has been pretty easy. I simple renamed my ".c" files to ".cpp", caught the few errors (darn it, I didn't typecast a few pointers from malloc) and the code compiled fine. My experience with this is also quite limited, YMMV. 

#### Calling the External Function

So, to start, let's look at the Python code that will call the C-function: 
    
    
    #
    # -- in someCode.py
    #
    # Call funName.  Pass it a list () of lists.  (sel1 and sel2 are lists.)
    # Get the return value into rValFromC.
    #
    rValFromC = funName( (sel1, sel2) );
    

where **sel1** and **sel2** could be any list of atom coordinates, say, from PyMOL. (See above.) 

Ok, this isn't hard. Now, we need to see what the code that receives this function call in C, looks like. Well, first we need to let C know we're integrating with Python. So, in your [header file](http://docs.python.org/api/includes.html) of **funName.h** we put: 
    
    
    // in funName.h
    #include <Python.h>
    

Next, by default your C-function's name is **funName_funName** (and that needs to be setup, I'll show how, later). So, let's define funName: 
    
    
    static PyObject*
    funName_funName(PyObject* self, PyObject* args)
    {
    ...more code...
    

This is the generic call. **funName** is taking two pointers to [PyObjects](http://docs.python.org/api/common-structs.html). It also returns a PyObject. This is how you get the Python data into and out of C. It shows up in "[args](http://docs.python.org/api/arg-parsing.html)" array of packaged Python objects and we then unpack it into C, using some [helper methods](http://docs.python.org/api/arg-parsing.html). Upon completion of unpacking, we perform our C/C++ procedure with the data, package up the results using the [Python API](http://docs.python.org/api/api.html), and send the results back to Python/PyMOL. 

#### Unpacking the Data

Let's unpack the data in **args**. Remember, **args** has a Python [list of lists](http://docs.python.org/lib/typesseq.html). So, to unpack that we do the following inside of funName: 
    
    
    static PyObject*
    funName_funName(PyObject* self, PyObject* args)
    {
           PyObject *listA, *listB;
    
           if ( ! PyArg_ParseTuple(args, "(OO)", &listA, &listB) ) {
                    printf("Could not unparse objects\n");
                    return NULL;
            }
    
            // let Python know we made two lists
            Py_INCREF(listA);
            Py_INCREF(listB);
     ... more code ...
    

Line 4 creates the two C objects that we will unpack the lists into. They are pointers to PyObjects. Line 6 is where the magic happens. We call, **[PyArg_ParseTuple](http://docs.python.org/api/arg-parsing.html)** passing it the args we got from Python. The **(OO)** is Python's code for _I'm expecting two _O_ bjects inside a list _()__. Were it three objects, then **(OOO)**. The first object will be put into **& listA** and the second into **& listB**. The exact [argument building specifications](http://docs.python.org/api/arg-parsing.html) are very useful. 

#### Reference Counting

Next, we check for success. Unpacking could fail. If it does, complain and quit. Else, **listA** and **listB** now have data in them. To avoid memory leaks we need to [manually keep track of PyObjects](http://docs.python.org/api/countingRefs.html) we're tooling around with. That is, I can create PyObjects in C (being sneaky and not telling Python) and then when Python quits later on, it'll not know it was supposed to clean up after those objects (making a leak). To, we let Python know about each list with **Py_INCREF(listA)** and **Py_INCREF(listB)**. This is [reference counting](http://docs.python.org/api/countingRefs.html). 

Now, just for safety, let's check the lists to make sure they actually were passed something. A tricky user could have given us empty lists, looking to hose the program. So, we do: 
    
    
         // handle empty selections (should probably do this in Python, it's easier)
         const int lenA = PyList_Size(listA);
         if ( lenA < 1 ) {
                 printf("ERROR: First selection didn't have any atoms.  Please check your selection.\n");
                 // let Python remove the lists
                 Py_DECREF(listA);
                 Py_DECREF(listB);
                 return NULL;
          }
    

We check the list size with, **[PyList_Size](http://docs.python.org/api/listObjects.html)** and if it's 0 -- we quit. But, before quitting we give control of the lists back to Python so it can clean up after itself. We do that with **Py_DECREF**. 

#### More Complex Unpacking

If you're dealing with simple scalars, then you might be able to skip this portion. 

Now, we should have access to the data the user sent us, in **listA** and **listB,** and it should be there and be clean. But, not forgetting that **listA** and **listB** are list of 3D coordinates, let's unpack them further into sets of coordinates. Because we know the length of the lists, we can do something like the following: 
    
    
           // make space for the current coords; pcePoint is just a float[3]
           pcePoint coords = (pcePoint) malloc(sizeof(cePoint)*length);
    
           // loop through the arguments, pulling out the
           // XYZ coordinates.
           int i;
           for ( i = 0; i < length; i++ ) {
                   PyObject* curCoord = PyList_GetItem(listA,i);
                   Py_INCREF(curCoord);
           
                   PyObject* curVal = PyList_GetItem(curCoord,0);
                   Py_INCREF(curVal);
                   coords[i].x = PyFloat_AsDouble(curVal);
                   Py_DECREF(curVal);
    
                   curVal = PyList_GetItem(curCoord,1);
                   Py_INCREF(curVal);
                   coords[i].y = PyFloat_AsDouble(curVal);
                   Py_DECREF(curVal);
    
                   curVal = PyList_GetItem(curCoord,2);
                   Py_INCREF(curVal);
                   coords[i].z = PyFloat_AsDouble(curVal);
                   Py_DECREF(curVal);
    
                   Py_DECREF(curCoord);
            }
    
     ... more code ...
    

Where, **pcePoint** is just a float[3]. Line 2 just gets some memory ready for the 3xlenght list of coordinates. Then, for each item for 1..length, we unpack the list using **[PyList_GetItem](http://docs.python.org/api/listObjects.html)** , into **curCoord**. This then gets further unpacked into the float[3], **coords**. 

We now have the data in C++/C data structures that the user passed from PyMOL. Now, perform your task in C/C++ and then return the data to PyMOL. 

### Sending the Results back to Python/PyMOL

Once you're done with your calculations and want to send your data back to PyMOL, you need to package it up into a Python object, using the Python API, and then return it. You should be aware of the expected return value and how you're packaging the results. If you user calls, 
    
    
    (results1,results2) = someCFunction(parameters1,parameters2)
    

then you need to package a list with two values. To build values for returning to PyMOL, use **[Py_BuildValue](http://www.python.org/doc/1.5.2p2/ext/buildValue.html)**. Py_BuildValue takes a string indicating the type, and then a list of values. [Building values](http://docs.python.org/ext/buildValue.html) for return has been documented very well. Consider an example: if I want to package an array of integers, the type specifier for two ints for Py_BuildValue is, "[i,i]", so my call could be: 
    
    
    # Package the two ints into a Python pair of ints.
    PyObject* thePair = Py_BuildValue( "[i,i]", int1, in2 );
    
    # Don't forget to tell Python about the object.
    Py_INCREF(thePair);
    

If you need to make a list of things to return, you iterate through a list and make a bunch of **thePairs** and add them to a Python list as follows: 
    
    
    # Make the python list
    PyObject* theList = PyList_New(0);
    # Tell Python about it
    Py_INCREF(theList);
    
    for ( int i = 0; i < someLim; i++ ) {
      PyObject* thePair = Py_BuildValue( "[i,i]", int1, in2 );
      Py_INCREF(thePair);
      PyList_Append(theList,thePair);
    

To add a list of lists, just make an outer list, 
    
    
    PyObject* outerList = PyList_New(0);
    

and iteratively add to it your inner lists: 
    
    
    PyObject* outerList = PyList_New(0);
    Py_INCREF(outerList);
    
    for ( int i = 0; i < someLim; i++ ) {
      // make the inner list, called curList;
      curList = PyObject* curList = PyList_New(0);
      Py_INCREF(curList);
    
      // fill the inner list, using PyList_Append with some data, shown above
      ...
    
      PyList_Append(outerList,curList);
    

Great, now we can extract data from Python, use it in C/C++, and package it back up for returning to Python. Now, we need to learn about the minimal baggage needed for C to operate with Python. Keep reading; almost done. 

### Initialization

We need to discuss how our functions will be called from Python. First, we need to create a [method table](http://docs.python.org/ext/methodTable.html). 
    
    
    static PyMethodDef CEMethods[] = {
            {"ccealign", ccealign_ccealign, METH_VARARGS, "Align two proteins using the CE Algorithm."},
            {NULL, NULL, 0, NULL}     /* Always use this as the last line in your table. */
    };
    

**[METH_VARARGS](http://docs.python.org/ext/methodTable.html)** can also be **METH_KEYWORDS** , where the former tells C that it should expect a simple tuple or list which we will unpack with **PyArg_ParseTuple** , and the latter tells C that it should expect to unpack the variables by name with the **PyArg_ParseTupleAndKeywords**. When using **METH_KEYWORDS** your function needs to accept a third parameter, a **Py_Object*** that is the dictionary of names for unpacking. For more information check out the [Python method table docs](http://docs.python.org/ext/methodTable.html). 

Each module undergoes initialization. By default the modules initialization function is: **initNAME()**. So, in our example above, **initccealign()". During this initialization step, we need to call[Py_InitModule](http://docs.python.org/ext/methodTable.html). For or above example, we'd have,**
    
    
    PyMODINIT_FUNC
    initccealign(void)
    {
        (void) Py_InitModule("ccealign", CEMethods);
    }
    

Finally, the main function that starts the whole shebang should look something like: 
    
    
    int
    main(int argc, char* argv[])
    {
            Py_SetProgramName(argv[0]);
            Py_Initialize();
            initccealign();
            return(EXIT_SUCCESS);
    }
    

At this point, you should have a fully functioning program in C/C++ intergrated with PyMOL/Python. 

### Installing Your Module

#### Overview

The [Python distutils pacakge](http://www.python.org/doc/2.2.3/ext/distributing.html) is a great method for distributing your modules over various platforms. It handles platform specific issues as well as simplifying the overall install process. For us, those module-builders, we need to create the distuils' setup.py script, and given the above -- that's the last step. 

More detailed information can be found one the Python documentation page for [installing C/C++ modules](http://docs.python.org/ext/building.html). There is also information on [how to build source and binary distribution packages](http://www.python.org/doc/2.2.3/ext/distributing.html). 

For example of how powerful disutils is, I have [cealign] setup to install as simply as: 
    
    
    python setup.py build cealign
    python setup.py install cealign
    

PyMOL also uses distutils for it's source-install. If more people understood distutils, I think they would install PyMOL from source since you get all the latest features. 

#### Setup.py

The setup file needs to know the following (at the very least): what source files comprise the project, what include directories to scan, the project name. You can also add more metadata such as version number, author, author_email, url, etc. For this example, let's assume we have the following directory structure, 
    
    
    .
    |-- build
    |-- dist
    |-- doc
    |   `-- funName
    |-- src
    |   |-- etc
    |   |   `-- tnt
    |   |       |-- doxygen
    |   |       |   `-- html
    |   |       `-- html
    |   `-- tnt
    

and we want to include all the _.cpp_ files from the **src** directory, and all the include files in **tnt**. We start setup.py as follows, 
    
    
    #
    # -- setup.py -- your module's install file
    #
    
    # import distutils
    from distutils.core import setup, Extension
    # for pasting together file lists
    from glob import glob
    # for handling path names in a os independent way
    from os.path import join;
    
    # grab all of the .h and .cpp files in src/ and src/tnt
    srcList = [ x for x in glob(join("src", "*.cpp")) ]
    # set the include directories
    incDirs = [ join( "src", "tnt") ]
    

Ok, now Python knows which files to include. Now we need to create a new [Extension](http://docs.python.org/dist/module-distutils.extension.html). We can simply call, 
    
    
    # create the extension given the function name, ''funName,'' the souce list and include directories.
    ccealignMods = Extension( 'funName', sources=srcList, include_dirs=incDirs  )
    

Lastly, all we have to do is call the final setup function, with the extension we just created and some metadata (if we want): 
    
    
    setup( name="funName",
            version="0.1-alpha",
            description="funName: A simple example to show users how to make C/C++ modules for PyMOL",
            author="Your Name Here",
            author_email="Your Email Goes Here",
            url="The URL of your work",
            ext_modules=[ccealignMods]
             )
    

And voila -- we're done. The users should now be able to execute, 
    
    
    python setup.py build
    # remove the brackets if you need to be root to install, see [Linux_Install#Installing_a_Script_Without_Superuser_Access Installing PyMOL w/o Superuser access] for an example.
    [sudo] python setup.py install
    

## Notes

  * discuss the pains of debugging



## Conclusion

I hope you found this helpful and will spur you to actually write some PyMOL modules or help you overcome the speed limitations inherent in Python's math (in comparison to other strictly-typed languages). 

I'm happy to hear any comments or questions you may have. [Tree](/index.php/User:Inchoate "User:Inchoate") 09:14, 19 May 2008 (CDT) 

## Example

See the source code for [cealign](/index.php/Cealign "Cealign"). 

  


## See Also

[stored](/index.php?title=Stored&action=edit&redlink=1 "Stored \(page does not exist\)"), [iterate_state](/index.php/Iterate_state "Iterate state"), [identify](/index.php/Identify "Identify"). 

Retrieved from "[https://pymolwiki.org/index.php?title=Advanced_Scripting&oldid=11745](https://pymolwiki.org/index.php?title=Advanced_Scripting&oldid=11745)"


---

## Biochemistry student intro

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Biochemistry course
  * 2 Install PyMOL to your computer
    * 2.1 Windows
    * 2.2 mac
    * 2.3 Extending pymol
  * 3 Find a suitable protein data file
    * 3.1 Read more about your protein
  * 4 Background
  * 5 Start your PyMOL
  * 6 The "Faster" way
    * 6.1 Movie of 1LEV
  * 7 Movie of Epidermal Growth Factor
    * 7.1 GUI - Scene loop
    * 7.2 GUI - Camera loop
    * 7.3 By a movie file
  * 8 Movie of Potassium Channel KcsA-Fab complex in high concentration of K+
  * 9 Export movies
    * 9.1 Collect images and create movie in "Windows Live Movie Maker"
  * 10 See Also



## Biochemistry course

This tutorial was written directly for biochemistry students at Copenhagen University, 2012/2013.  


This is the very first introduction to the powerful molecular visualizer PyMOL.  
We will only cover the very basic steps to get a image of your enzyme and put in your written student article. 

If you want to read about PyMOL, then try this introduction [Practical_Pymol_for_Beginners](/index.php/Practical_Pymol_for_Beginners "Practical Pymol for Beginners")

## Install PyMOL to your computer

You first want to install PyMOL.   
Do this, by following this guide.  


### Windows

[Windows Pre-compiled_PyMOL](http://www.pymolwiki.org/index.php/Windows_Install#Pre-compiled_PyMOL). Consider this little nice texteditor: [Notepad++](http://notepad-plus-plus.org/download/)   


### mac

[pymol for mac](http://erlendsson.dk/PSA/). Download and install MacPyMOL-v1.3r1-edu.tar.bz2 

### Extending pymol

You don't have to follow the steps of extending PyMOL, but if you are a little technical, and want to become friends with PyMOL over time, then consider it.  


## Find a suitable protein data file

We need to find a Protein Databank File (PDB), which describe the x,y,z coordinates of your enzyme.  
These are stored at the homepage: <http://www.rcsb.org>   
Find a suitable file by searching for: **porcine kidney fructose 1,6-bisphosphatase**   
The PDB file **1LEV** , seems suitable. 

### Read more about your protein

There exist homepages, to get more info about your protein.  
These homepages collect material from several sources, and present them in easy format. 

<http://www.proteopedia.org/wiki/index.php/1lev>   
<http://pdbwiki.org/wiki/1lev>

## Background

This tutorial is designed to give you a basic working knowledge of making pretty and informative pictures of protein structures using PyMOL.   
This tutorial does not cover all the functions of PyMol, but tries to focus on the most important ones.   
A great resource for more advanced use is this wiki and the large number of good tutorials found online, which can be accessed via google.  


When you open PyMOL you will see two windows opening.  
The upper window entitled **PyMOL Tcl/Tk GUI** controls the general settings and functionalities of the program while the lower one entitled **PyMOL viewer** , contains the settings that are related to the current display of the molecule.   
Each of these two windows contain a command line, where we can enter commands into the program.   
One can, however, get really far without ever having to worry about writing any commands. 

First, we will need to load a structure file into the program. Protein structures are deposited at [www.pdb.org](/index.php?title=Www.pdb.org&action=edit&redlink=1 "Www.pdb.org \(page does not exist\)")  
There are several ways to open this file:  
Either we can download the file from PDB and save it to our computer.  
We can then open the PDB file (if you have it on the harddisk) using the menu **File > Open**.  
or we can write: **fetch** followed by the PDB id in either of the two command lines. Ex: **fetch 1LEV**   
The Viewer window should now contain the PDB file displayed with lines. 

In the right hand side of the viewer window there is a selection menu that currently contains two lines:  
**all** and the four character name of the PDB file, ex. **1lev**.  
Later on we will make new selections that will appear in this menu.  
Each of the two lines have 5 buttons labeled **A** (Actions), **S** (Show), **H** (Hide), **L** (Label), **C** (Color).  
Press these buttons to see the options available in each menu. 

**Hide the line representation**. Press **H > everything** in the **all** line. This will remove the line representation of the molecule.  
**Visualize the molecule in the cartoon representation** by selecting **S > cartoon**.  
Try also to show the molecule as ribbon, sticks or surface to get a feeling for the different representations.  
The cartoon option gives us the best overview of the overall structure, so hide the other representations and show the protein only in cartoon mode.  
Color everything white selecting **C > grey > white**. 

Ok, now we are going to try to move the molecule around.   
Try to move the mouse around while pressing the left mouse button. This will rotate the molecule.  
Try to move it around holding the middle button (moves the view) or the right button (zooms in and out).  
If you have a one or two button mouse you can change the mode in the **Mouse** menu in the Tcl/Tk GUI window.   
The box in the lower right corner will show how to move and rotate the molecule in this case. 

## Start your PyMOL

Start your shortcut to PyMOL "C:\Python27\PyMOL\PyMOL.exe"  


Now click and do the following: 

  1. Write in command line: **fetch 1lev**
  2. Right Menu: -> 1lev -> "A" -> preset -> Publication
  3. Top Menu: Display -> Sequence
  4. Top Menu: Display -> Sequence Mode -> Chains
  5. In Sequence, select so all "F" is marked.
  6. Right Menu: -> (sele) -> "A" -> remove atoms. (1LEV is crystallized in dimer formation. So we only need to view 1 chain)
  7. Top Menu: Display -> Sequence Mode -> Residue Codes
  8. In Sequence, select so only substrate (F6P) is marked (A/338)
  9. Right Menu: -> (sele) -> "A" -> zoom
  10. Right Menu: -> (sele) -> "A" -> rename selection -> f6p
  11. Right Menu: -> (f6p) -> "C" -> by element -> Select to Carbon is not green
  12. Right Menu: -> (f6p) -> "A" -> find polar contacts -> to others excluding solvent
  13. In Sequence, select so only MN is marked (A/340)
  14. Right Menu: -> (sele) -> "A" -> rename selection -> mn
  15. Right Menu: -> (mn) -> "S" -> spheres
  16. Write in console: **select act_site, byres f6p around 3.5**
  17. Right Menu: -> (act_site) -> "S" -> sticks
  18. Right Menu: -> (act_site) -> "L" -> residues
  19. Top Menu: Display -> Background -> White
  20. Find a good view, and push "Ray" in the top right of the grey command console.
  21. Top Menu: File -> Save Image As -> PNG
  22. Put it into your student article



## The "Faster" way

The real power of PyMOL, comes into power, when your write a PyMOL command file.  
Here you write which commands pymol should execute, and so it only take 1 second to get the same.  
The commands are stored in a ".pml" file. 

Open Notepad, and then: **File- >Save as->All files-> C:\Users\YOU\pymol\1lev.pml**
    
    
    # Best to restart PyMOL every time from fresh
    reinitialize
    cd C:\Users\DIG\pymol
    
    fetch 1lev, async=0
    preset.publication(selection='all')
    remove chain F
    select substrates, organic
    select f6p, resn F6P
    zoom f6p
    util.cbac('f6p')
    
    select act_site, byres f6p around 3.5
    show sticks, act_site
    
    distance pol_cont, f6p, act_site, mode=2
    
    select cli, /1lev//A/CLI # OR: select cli, resn CLI
    select mn, name MN # OR: select mn, symbol Mn # OR: select mn, inorganic
    show spheres, mn
    
    label act_site and name CB, resn+resi
    
    zoom pol_cont
    viewport 1024,768
    bg_color white
    ray
    png 1lev.png
    

  1. "#" Line starting with hashes is not read by PyMOL. Use a comment field.
  2. Best to restart PyMOL every time from fresh
  3. Go to your working directory
  4. Get the pdb file from the RCSB server. async=0 makes sure it waits for completion of download before continuing.



Then you just open the .pml file with PyMOL.   
Or start PyMOL, and write in command: **@1lev.pml**   
Or start PyMOL, top menu -> File -> Run... -> C:\Users\YOU\pymol\1lev.pml 

### Movie of 1LEV

File: **1lev_movie.pml**  
Note, you need the **movie.pml** file in same directory, see [Biochemistry_student_intro#By_a_movie_file](/index.php/Biochemistry_student_intro#By_a_movie_file "Biochemistry student intro"). 
    
    
    fetch 1lev, async=0
    preset.publication(selection='all')
    remove chain F
    select substrates, organic
    select f6p, resn F6P
    zoom f6p
    util.cbac('f6p')
     
    select act_site, byres f6p around 3.5
    show sticks, act_site
     
    distance pol_cont, f6p, act_site, mode=2
     
    select cli, /1lev//A/CLI # OR: select cli, resn CLI
    select mn, name MN 
    # OR: select mn, symbol Mn # OR: select mn, inorganic
    show spheres, mn
     
    label act_site and name CB, resn+resi
     
    zoom pol_cont
    viewport 1024,768
    bg_color white
    
    ###### Movie
    # Reset
    hide sticks, act_site
    disable pol_cont
    disable cli
    hide spheres, mn
    disable mn
    
    zoom 1lev
    scene F1, store, Publication
    
    zoom substrates
    scene F2, store, substrates
     
    zoom f6p
    scene F3, store, f6p
     
    show sticks, act_site
    enable pol_cont
    show spheres, mn
    
    zoom pol_cont
    scene F4, store, polar contacts
    
    @movie.pml
    #ray
    #png 1lev.png
    

## Movie of Epidermal Growth Factor

Lets make a movie of the Molecule of the month on RCSB.org.  
Let's take a membrane protein, [1nql@rcsb](http://www.rcsb.org/pdb/explore/explore.do?structureId=1nql), [pdbwiki](http://pdbwiki.org/wiki/1nql), [proteopedia](http://www.proteopedia.org/wiki/index.php/1nql)  
[June 2010 Molecule of the Month by David Goodsell](http://www.rcsb.org/pdb/101/motm.do?momID=126)

Open Notepad, and then: **File- >Save as->All files-> C:\Users\YOU\pymol\1nql.pml**
    
    
    reinitialize
    
    fetch 1nql, type=pdb1, multiplex=1,async=0
    
    # So we get buttons for scenes
    set scene_buttons, 1
    viewport 1280,800
    bg_color white
    set fog_start, 0.60
    
    #### Scene 1, publication ####
    show_as cartoon, all
    preset.publication(selection='all')
    
    # See http://pymolwiki.org/index.php/Single-word_Selectors
    extract substrates, organic
    util.cbac('substrates')
    select others, inorganic
    select water, solvent
    # See http://pymolwiki.org/index.php/Selection_Macros
    # Select organic molecules
    select nag, ////NAG
    select bma, ////BMA
    disable bma
    # Hide organic
    hide everything, substrates
    
    # Save scene 1
    
    zoom 1nql
    scene F1, store, Publication
    
    #### Scene 2, show cysteines ####
    # Select sulfurs, since they play a role
    select sulf_cys, resn cys and not (name O or name N or name C)
    disable sulf_cys
    show sticks, sulf_cys
    color sulfur, sulf_cys and elem S
    set_view (\
         0.351475894,    0.052040517,    0.934747040,\
        -0.784766018,   -0.528064251,    0.324480295,\
         0.510492086,   -0.847605526,   -0.144762829,\
         0.000000000,    0.000000000, -258.137542725,\
        47.065788269,  -10.656063080,    0.561561584,\
       217.718627930,  298.556427002,  -20.000000000 )
    scene F2, store, Cysteines
    
    #### Scene 3, show organis ####
    show sticks, substrates
    set_view (\
        -0.350075662,   -0.229279995,   -0.908223689,\
         0.627038240,    0.662946343,   -0.409051389,\
         0.695891201,   -0.712693632,   -0.088314489,\
        -0.000047103,    0.000066929, -180.624176025,\
        29.827386856,   19.096229553,   14.554395676,\
       144.383483887,  216.888458252,  -20.000000000 )
    scene F3, store, substrates
    
    #### Scene 4, polar contacts ####
    select prot_cont, byres substrates around 3.5
    distance subs_bond, substrates, prot_cont, mode=2
    show sticks, prot_cont
    set_view (\
        -0.564916074,    0.749067366,    0.346062332,\
         0.049069319,   -0.388159811,    0.920283496,\
         0.823683023,    0.536867142,    0.182521686,\
         0.000190482,   -0.000248071,  -91.693801880,\
        25.880264282,   29.162174225,   10.565655708,\
        55.450801849,  127.955612183,  -20.000000000 )
    scene F4, store, polar contacts
    
    #@movie.pml
    

### GUI - Scene loop

  1. Click **Movie** ->**Program** ->**Scene Loop** ->**Steady** ->**Program** ->**4 seconds each**
  2. Click the **Play button** at the lover right corner. Or write **mplay**.
  3. Click **Movie** ->**Reset**
  4. Try also the other **Scene Loop** method.



### GUI - Camera loop

  1. Click **Mouse** ->**# Button Motions**. In the selection menu, under **all** , you now have a **M** button!.
  2. **Movie** ->**Reset**
  3. Click **F1**
     1. sel. menu **all-[M]** ->**Store with scene F1**
     2. **Movie** ->**Program** ->**Camera Loop** ->**Y-Roll** ->**4 seconds**
     3. **Movie** ->**Append** ->**2 seconds**
     4. Click lower right **Full forward button** **- >**
  4. Click **F2**
     1. sel. menu **all-[M]** ->**Store with scene F2**
     2. **Movie** ->**Program** ->**Camera Loop** ->**X-Rock** ->**60\. deg. over 4 sec.**
     3. **Movie** ->**Append** ->**2 seconds**
     4. Click lower right **Full forward button** **- >**
  5. Click **F3**
     1. sel. menu **all-[M]** ->**Store with scene F3**
     2. **Movie** ->**Program** ->**Camera Loop** ->**X-Roll** ->**4 seconds**
     3. **Movie** ->**Append** ->**2 seconds**
     4. Click lower right **Full forward button**
  6. Click **F4**
     1. sel. menu **all-[M]** ->**Store with scene F4**
     2. **Movie** ->**Append** ->**1 seconds**
     3. Click lower right **Full forward button**
     4. sel. menu **all-[M]** ->**Store with scene F4**
     5. **Movie** ->**Program** ->**Camera Loop** ->**Nutate** ->**30\. deg. over 8 sec.**
     6. **Movie** ->**Append** ->**4 seconds**
     7. Click lower right **Full backward button**



Play 

### By a movie file

Open Notepad, and then: **File- >Save as->All files-> C:\Users\YOU\pymol\movie.pml**
    
    
    set movie_panel, 1
    mset 1 x1000
    
    scene F1
    python
    
    f=1
    cmd.frame(f)
    cmd.mview("store",scene="F1")
    f=f+99; print f
    cmd.frame(f)
    cmd.mview("store",scene="F1")
    
    f=f+49; print f
    cmd.frame(f)
    cmd.mview("store",scene="F2")
    f=f+99; print f
    cmd.frame(f)
    cmd.mview("store",scene="F2")
    
    f=f+49; print f
    cmd.frame(f)
    cmd.mview("store",scene="F3")
    f=f+99; print f
    cmd.frame(f)
    cmd.mview("store",scene="F3")
    
    f=f+49; print f
    cmd.frame(f)
    cmd.mview("store",scene="F4")
    f=f+99; print f
    cmd.frame(f)
    cmd.mview("store",scene="F4")
    
    f=f+49; print f
    cmd.frame(f)
    cmd.turn('y',50)
    cmd.mview("store")
    
    f=f+99; print f
    cmd.frame(f)
    cmd.turn('y',-100)
    cmd.mview("store")
    
    f=f+49; print f
    cmd.frame(f)
    cmd.turn('y',50)
    cmd.mview("store")
    
    f=f+49; print f
    cmd.frame(f)
    cmd.mview("store",scene="F4")
    
    python end
    
    frame 1
    mplay
    

Write in pymol: **@movie.pml**

## Movie of Potassium Channel KcsA-Fab complex in high concentration of K+

Let's take another membrane protein, [1k4c@rcsb](http://www.rcsb.org/pdb/explore/explore.do?structureId=1k4c), [pdbwiki](http://pdbwiki.org/wiki/1k4c), [proteopedia](http://www.proteopedia.org/wiki/index.php/1k4c)  


Open Notepad, and then: **File- >Save as->All files-> C:\Users\YOU\pymol\1k4c.pml**
    
    
    reinitialize
    
    fetch 1k4c, type=pdb1, multiplex=1,async=0
    
    # So we get buttons for scenes
    set scene_buttons, 1
    #viewport 320,200
    bg_color white
    set fog_start, 0.60
    
    #### Scene 1, publication ####
    group 1k4c, 1k4c_* 
    preset.publication(selection='1k4c')
    
    # See http://pymolwiki.org/index.php/Single-word_Selectors
    extract substrates, organic
    util.cbac('substrates')
    extract ions, inorganic
    show spheres, ions
    disable ions
    
    extract water, solvent
    show nonbonded, water
    color grey, water
    disable water
    # See http://pymolwiki.org/index.php/Selection_Macros
    # Select organic molecules
    # Hide organic
    hide everything, substrates
    
    # Save scene 1
    zoom 1k4c
    scene F1, store, Publication
    
    #### Scene 2, show cysteines ####
    # Select sulfurs, since they play a role
    util.cbag('1k4c_0001')
    util.cbao('1k4c_0002')
    util.cbas('1k4c_0003')
    util.cbaw('1k4c_0004')
    
    select sulf_cys, resn cys and not (name O or name N or name C)
    disable sulf_cys
    show sticks, sulf_cys
    color sulfur, sulf_cys and elem S
    color red, sulf_cys and not elem S
    set_view (\
         0.443872631,   -0.125126541,    0.887310863,\
         0.892753899,   -0.023616746,   -0.449926227,\
         0.077252157,    0.991857469,    0.101223812,\
         0.000020482,   -0.000164529, -253.516693115,\
       164.839996338,  190.414398193,   -7.539838314,\
       168.313964844,  338.728668213,  -20.000000000 )
    scene F2, store, Cysteines
    
    #### Scene 3, show organis ####
    show sticks, substrates
    set_view (\
        -0.975404024,    0.029795930,   -0.218399704,\
        -0.214692041,    0.096028417,    0.971945703,\
         0.049933217,    0.994929075,   -0.087270327,\
         0.000279146,   -0.000071049, -226.327316284,\
       160.358612061,  152.090927124,  -33.061843872,\
       141.134780884,  311.549621582,  -20.000000000 )
    scene F3, store, substrates
    
    #### Scene 4, polar contacts ####
    select prot_cont, byres 1k4c_* within 3.5 of (substrates or ions)
    util.cbay('prot_cont') 
    distance subs_bond, substrates, prot_cont, mode=2
    show sticks, prot_cont
    enable ions
    set_view (\
        -0.961359501,   -0.272781283,   -0.037127689,\
        -0.273679018,    0.961574137,    0.021663748,\
         0.029792555,    0.030986462,   -0.999071598,\
         0.000279146,   -0.000071049, -226.327316284,\
       160.358612061,  152.090927124,  -33.061843872,\
       141.134780884,  311.549621582,  -20.000000000 )
    scene F4, store, polar contacts
    
    #@movie.pml
    

## Export movies

Write in command 
    
    
    frame 1
    

Then go to File->Save Movie As -> MPEG   
It takes a little time, 5 min. Be patient. 

* * *

You can also export it like PNG images, and then collect them.   
Then go to File->Save Movie As -> PNG Images 

If you wan't to have a really nice movie, you want to have all images ray-traced. Then do 
    
    
    set ray_trace_frames, 1
    

But note, this takes a **looong** time, so wait until you are SURE that you have the final movie. Try first to make and collect a video without. 

#### Collect images and create movie in "Windows Live Movie Maker"

Start or install [Windows Live Movie Maker](http://windows.microsoft.com/da-dk/windows7/products/features/movie-maker)

  1. Add all images from folder
  2. Ctrl+a to select all images
  3. Click in tab "Edit", set Duration to: 0,03
  4. Click in tab "Home", "Save movie", "Windows Phone (large)"



You are done 

## See Also

  * [ 3d protein image in pdf](/index.php/3d_pdf "3d pdf")
  * [ Movies in pdf](/index.php/Movie_pdf "Movie pdf")
  * [Example pymol movie](http://youtu.be/r6o35ZRsdqw)



Retrieved from "[https://pymolwiki.org/index.php?title=Biochemistry_student_intro&oldid=11555](https://pymolwiki.org/index.php?title=Biochemistry_student_intro&oldid=11555)"


---

## MovieSchool

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**PyMOLWiki's Movie School**

## Welcome

Here you will learn all the latest and greatest techniques for making movies in PyMOL. These make you look like a **master of PyMOL** when you can show them off in a presentation or even in MPEG/AVI format. For your ease of learning I have commented all of the movie scripts and tried to ensure that they work in PyMOL if you just copy/paste the code. 

If you have any comments on this tutorial, please feel free to add them to the discussion page. If you have corrections, improvements or additions, as always, please feel free to improve the pages. 

  * NB: While this is pretty expansive (hence the 6 lessons) it is still incomplete and surely rough in some places.
  * NB: Make sure you have the newest [PyMOL v1.2r0](http://www.pymol.org/).



## Movie Making Roadmap

  * [Introduction & Movie Making for the Impatient](/index.php/MovieSchool_1 "MovieSchool 1")



    

    The basics; **movies in 60 seconds**.

  * [Terminology & Movie Related Commands](/index.php/MovieSchool_2 "MovieSchool 2")



    

    Definitions; commands; settings;

  * [PyMOL GUI for Movie Making](/index.php/MovieSchool_3 "MovieSchool 3")



    

    Using the GUI to make movies

  * [Simple Movie Examples](/index.php/MovieSchool_4 "MovieSchool 4")



    

    Easy copy/paste examples

  * [Putting it All Together](/index.php/MovieSchool_5 "MovieSchool 5")



    

    More advanced examples; using states, scenes, and object motions

  * [Exporting & Using your Movie](/index.php/MovieSchool_6 "MovieSchool 6")



    

    Using your movie, once it's made

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool&oldid=7196](https://pymolwiki.org/index.php?title=MovieSchool&oldid=7196)"


---

## MovieSchool 1

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Movie Making
  * 2 Your First Movie
  * 3 Movie Making for the Impatient
    * 3.1 Simple Camera Motions
      * 3.1.1 Simple 360 Scene Rotation
      * 3.1.2 Simple 30, 60, 90, 120, 180, Scene Rocking
      * 3.1.3 Nutate
      * 3.1.4 Zooming Around an Object
      * 3.1.5 Real-world Example



## Movie Making

While PyMOL's capability to produce static images is quite powerful, there are some stories that are better told through movies, than static images alone. This little page will provide the necessary ideas, links, code and examples for making movies in PyMOL. 

## Your First Movie

Movies can be very simple, for example, animating an NMR ensemble: 
    
    
    # Your first movie.
    fetch 1nmr
    mplay
    
    # to stop the movie when you're ready
    # type 'mstop'.
    

What PyMOL did here was to [fetch](/index.php/Fetch "Fetch") the file from the PDB and load it into an object with 20 states. Somewhere between then and issuing [mplay](/index.php/Mplay "Mplay") PyMOL created 20 frames for your object and assigned one state to each frame. This created the animated effect as we scroll through the frames. 

## Movie Making for the Impatient

If you don't have time to or care to make more complex movies and want a movie _now_ then read this section. It mostly involves the GUI for making movies, so get your mouse ready. 

### Simple Camera Motions

#### Simple 360 Scene Rotation

To make a movie that simply rotates around your scene 360 degrees, in the menu click on: 

    

    **Movie- >Program->Camera->X-Roll->N Seconds**

for an N-second movie rolling over the X-axis. Chose 

    

    **Movie- >Program->Camera->Y-Roll->N Seconds**

for the Y-axis roll. 

Done! Press **play**. 

#### Simple 30, 60, 90, 120, 180, Scene Rocking

This will show up to a 30, 60, 90, 120 or 180 rocking 'wedge' of the scene. If you don't want to rotate all the way around, use this. 

    

    **Movie- >Program->Camera->X-Rock->X-Degrees over N-Seconds**

Done! Press **play**. 

#### Nutate

Nutating is like a wiggle-rock; try it and see. 

    

    **Movie- >Program->Camera->X-Rock->X-Degrees over N-Seconds**

Done! Press **play**. 

#### Zooming Around an Object

This is also known as camera movement. Let's make a simple program that just zooms in and out on a some atom. 

    

    **Build- >Residue->Tryptophan**
    **Scene- >Store->F1**
    _Click any atom and zoom on it_
    **Scene-Store- >F2**
    **Movie- >Program->Scene Loop->** and if you want to nutate while at the scene choose **Nutate** otherwise choose **Y-Rock**.

Done! Press **play**. 

  * _Hint_ : Make multiple scenes, making sure to store each one. Then use the Scene Loop to automatically program a movie for you!



#### Real-world Example

First load the tutorial PDB: 
    
    
    load $TUT/1hpv.pdb
    

    

    **Action- >Preset->Technical** (on the object in the viewer gui)
    **Scene- >Store->F1**
    Type _zoom i. 200_ to zoom on the ligand
    **Scene- >Store->F2**
    **Movie- >Program->Scene Loop->Y-Rock->4 Seconds Each**

Done! Press **play**. 

[ ← MovieSchool Home](/index.php/MovieSchool "MovieSchool") [ Next Lesson →](/index.php/MovieSchool_2 "MovieSchool 2")

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool_1&oldid=10903](https://pymolwiki.org/index.php?title=MovieSchool_1&oldid=10903)"


---

## MovieSchool 2

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Terminology
    * 1.1 Basic Movie Terminology
  * 2 What is a Movie?
  * 3 Movie Making Commands
    * 3.1 frame
    * 3.2 set state
      * 3.2.1 States & Frames (optional reading)
    * 3.3 mset
    * 3.4 mview
      * 3.4.1 Very Basic Mview Syntax



## Terminology

Imagine a complex movie for a moment, a movie that has camera motions, objects moving independently or in concert with other objects, changing colors and representations. To handle camera motions PyMOL must know at all times where the camera is located and what it's pointed toward (as well as clipping planes). For objects to move around or be rotated without regard to the camera (the objects themselves rotate/translate, not just the camera), PyMOL has to store the coordinates and matrices for these objects, too. Changing colors and representations for each object must somehow also be stored. So, as you can see this is a multidimensional problem: at each time point in your movie, PyMOL must remember positions and representations, as well as make it easy for you to transition between them (interpolation). 

Despite these complexities, PyMOL tries to enable movie making for even novice users. Let's start by defining a few PyMOL concepts—states, frames and scenes. 

### Basic Movie Terminology

[![](/images/b/b1/Movies.png)](/index.php/File:Movies.png)

[](/index.php/File:Movies.png "Enlarge")

Schematic of how frames, scenes and states link together.

**object**

    

    An object is any PyMOL-object loaded into PyMOL, like atoms, molecules, complexes, etc. When you load an PDB from disk/net it is loaded into PyMOL as an object.
    [All pages regarding objects](/index.php/Category:Objects_and_Selections "Category:Objects and Selections")

**selection**

    

    A selection is a specifically chosen set of atoms, molecules, complexes etc. in PyMOL. A selection is not an object, it's a subset of stuff from a (collection of) object(s). Selections can be named and when named have are distinguished from objects by having parentheses around their names. For example, _foo_ would be an object and _(foo)_ would be some selection. When you pick an atom (and get the default **(sele)** selection) or issue the ever-popular [Selection](/index.php/Select "Select") command, you get a selection.
    [All pages regarding selections](/index.php/Category:Selections "Category:Selections")

**states**

    

    A state is a particular conformation (set of coordinates) for a given object. For example an NMR ensemble could contain the same molecule, but in 20 different states. PyMOL can make movies from states. States **do not store representations** in PyMOL (eg. cartoons vs. sticks).
    See also [All pages regarding states](/index.php/Category:States "Category:States")

**scenes**

    

    A **scene** as I understand stores the orientation/position of the camera, all object activity information, all atom-wise visibility, color, representations, and the global frame index.
    See [Category:Scenes](/index.php/Category:Scenes "Category:Scenes") and [Scene](/index.php/Scene "Scene") for more information. The [Scene](/index.php/Scene "Scene") page has information on how to quickly add/remove and change scenes--_suggested reading!_

**interpolation**

    

    A scene is the staged representations of objects and the orientation of the camera.
    See also [All pages regarding scenes](/index.php/Category:Scenes "Category:Scenes")

**frames**

    

    A frame can be thought of as a single frame in a movie reel. A frame stores state information and scene information.
    See also [All pages regarding frames](/index.php?title=Category:Frames&action=edit&redlink=1 "Category:Frames \(page does not exist\)")

**Movie Panel**

    

    The movie panel is a frame indicator strip at the bottom of the screen. It shows a little icon for which frame you're currently on, and whether or not the camera has been set for that frame.
    See [movie_panel](/index.php/Movie_panel "Movie panel") for more information.

## What is a Movie?

Now the we have the appropriate terminology to talk about movies in PyMOL, we can discuss what a movie really is. A movie in PyMOL is a series of frames stitched together in some way so as to create the desired animation. Once a movie is made in PyMOL, we have a few options for exporting it to other formats like png arrays or MPEG moves. 

## Movie Making Commands

This tutorial assumes you have some basic knowledge about how to use PyMOL (eg. mousing, choosing and setting your representations, etc). If you're not yet at this level, please check out [the Tutorial Category](/index.php/Category:Tutorials "Category:Tutorials") of pages (most notably the beginner tutorials). 

I think it's help to think of the movie as a set of frames, like in a movie reel, so let's start there. (Each command below links to the command's PyMOL wiki page, so feel free to click through for more info.) 

### [frame](/index.php/Frame "Frame")

This command tells PyMOL to set the current frame to whichever you desire. To use it, just issue the command, 
    
    
    frame X
    

where **X** is some integer number indicating the frame you want to go to. If you issue a frame number greater than the number of frames, PyMOL sets the frame to the highest-numbered frame you have (similarly for negative numbers or numbers smaller than the minimum numbered frame). 

Let's try a quick example with [frame](/index.php/Frame "Frame"), 
    
    
    # create an empty 90 frame movie
    mset 1 x90
    # turn on the movie panel (bottom of screen)
    set movie_panel, on
    # goto frame one
    frame 1         
    # try some intermediate frames; notice the blue indicator in the movie panel         
    frame 10
    frame 50
    frame 90
    # try going beyond the end and see what PyMOL does
    frame -1
    frame 100
    # play through the empty movie
    mplay
    # stop the movie
    mstop
    

  


### [set state](/index.php/States "States")

Again, states are particular stored conformations of objects. Here we use PyMOL to set and get the states, and see how PyMOL mapped them to our earlier movie example. 

This command has a similar idea of [frame](/index.php/Frame "Frame"), but works a little differently. Instead of typing, 
    
    
    # invalid command
    state 40
    

in PyMOL we [set](/index.php/Set "Set") the [state](/index.php/States "States"): 
    
    
    # how to set a state in PyMOL
    set state, stateNo, objectName
    
    # for example
    # set state to 40
    set state, 40
    
    # also, get state
    get state
    

#### States & Frames (optional reading)

As an example, look at the code from the "first movie": 
    
    
    fetch 1nmr
    mplay
    # issue mstop, to stop the movie
    

We can do a couple things now, let's try counting the number of states and frames PyMOL now knows about: 
    
    
    # how many frames does PyMOL know about?
    count_frames
    
    # what about states?
    count_states
    

and now let's see how PyMOL mapped frames to states. Using the above commands and a little Python, let's see how PyMOL mapped the frames to states: 
    
    
    python
    for x in range(1,cmd.count_frames()+1):
      cmd.frame(x)
      print "Frame => %s; and State => %s" % ( str(x), str(cmd.get('state')))
    python end
    

which should show a 1-1 mapping of states to frames. 

### [mset](/index.php/Mset "Mset")

[mset](/index.php/Mset "Mset") is a very powerful command. This command tells PyMOL how to assign states to frames. So, now you see why it's necessary to clearly distinguish between and use the two. Let's learn how to use [mset](/index.php/Mset "Mset"). 

The syntax for [mset](/index.php/Mset "Mset") can be a little tricky at first. I would write the syntax as: 
    
    
    mset stateSpec, frameSpec
    

which assigns the states in **stateSpec** to the frames in **frameSpec**., where **stateSpec** is any mset-valid state specification. PyMOl supports to patterns for **stateSpec**. You can do simply supply a number eg 
    
    
    mset 1
    

or you can specify a range of states—like 1 through 55 as 
    
    
    # setting states 1 through 55
    # caution: notice the space: 1 -55, not 1-55 (this is a PyMOL parser caveat)
    mset 1 -55
    

Simple enough. Now for **frameSpec** you can specify a single frame number like so or you can specify _how many frames PyMOL should use to map to your specified states_ with the **xNumber** command. This will make sense with an example 
    
    
    # Recall: mset stateSpec frameSpec
    # so we are setting STATE 1 across a span of 90 frames
    mset 1 x90
    
    # Recall: mset stateSpec frameSpec
    # so we are setting states 1..120 to the next 120 frames
    mset 1 -120 x120
    

NB: Actually the syntax is a little more complicated than this as PyMOL's mset command has the ability to remember in which frame the prior specification left the movie. So, you can sort of chain the specifications. Type _help mset_ in PyMOL for more info or see [mset](/index.php/Mset "Mset"). 

### [mview](/index.php/Mview "Mview")

The [mview](/index.php/Mview "Mview") command can be intimidating, but all we need to know about it at present is that it can store the (1) camera position or (2) a given object position. The idea is to essentially make 'way points' in your movie and have PyMOL interpolate the in-between positions/coordinates, etc. For example, if I wanted to make a 100-frame movie of a zoom into an object, I could store and manually set 100 camera positions, or I could do the starting position and the final position and ask PyMOL to just interpolate that over 100 frames. The latter idea is obviously much simpler. So simple in fact, let's make a super-quick movie that does exactly what I just mentioned—100 frames of a slow zoom into some object. Start with a fresh PyMOL session ([reinitialize](/index.php/Reinitialize "Reinitialize")) and then copy/paste the following: 
    
    
    # let's initialize the movie to 100 frames, all of state 1
    # it's ONLY state 1, because we're only moving the camera around, not
    # changing structure coordinates of the leucine:
    mset 1 x100
    
    # show a leucine
    frag leu
    
    # position the residue
    orient
    
    # let's store the current camera position
    mview store
    
    # now set our way point to be frame 100
    frame 100
    
    # now let's zoom into some atom on the fragment
    zoom ID 10
    
    # now save this view at frame 100
    mview store
    
    # last thing is to tell PyMOL to interpolate the 100 frame zoom
    # so we don't have to do those 100 snapshots:
    mview reinterpolate
    
    # voila, you have a movie.  To watch it go back to frame 1 and play it
    frame 1
    mplay
    
    # mstop when you're ready
    

'But, hold on!' you might say. Why is it so herky-jerky? We have smooth zooming but then a snap and back to frame one! Well, we never gave PyMOL any number of frames to interpolate the change from frame 100 back to frame 1 (since it wraps). If we wanted an a **zoom in** that was equally as fast as the **zoom out** we would simply **replace line #16** with 
    
    
    # now set our way point to be frame 100
    frame 100
    

but, if we wanted a slow, zoom in and a _fast_ zoom out we could do 
    
    
    # now set our way point to be frame 100
    frame 80
    

instead which would only give PyMOL 20 frames with which to zoom us out. Try it! 

#### Very Basic Mview Syntax

_This is a simple overview, see[mview](/index.php/Mview "Mview") for complete details._

Using [mview](/index.php/Mview "Mview") as you can see from above is pretty simple. The very basic syntax is: 
    
    
    # store the camera OR some object given by objName (if it's supplied)
    mview store[, object=objName]
    
    # reinterpolate (link together the positions) for the saved camera or object
    mview reinterpolate[, object=objName]
    

If you leave off the **object=objName** then you're storing the **camera information only** —and so none of your objects will be moving anywhere—just the camera. If you include an object name, then it stores that object's position information in the current frame. The **mview store** tells PyMOL to store the camera or objects coordinates, while the **mview reinterpolate** command tells PyMOL to link together the saved positions for the camera or the selected object in a smooth, cool way. 

[ ← Previous Lesson](/index.php/MovieSchool_1 "MovieSchool 1") [ Next Lesson →](/index.php/MovieSchool_3 "MovieSchool 3")

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool_2&oldid=11154](https://pymolwiki.org/index.php?title=MovieSchool_2&oldid=11154)"


---

## MovieSchool 3

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Movie Making & the PyMOL GUI
    * 1.1 Basic TK GUI
      * 1.1.1 Movie Menu
      * 1.1.2 Scene Menu
        * 1.1.2.1 Scene's Recall Submenu
      * 1.1.3 Mouse Menu
    * 1.2 openGL GUI
    * 1.3 Mouse Modes
      * 1.3.1 Motions Mode



## Movie Making & the PyMOL GUI

For those people who prefer to use the mouse alot in PyMOL, this page is for you. The GUI is also quite useful for doing some basic utility functions like appending the rocking/rolling of a molecule in your movie. 

  * These screenshots are all current as of PyMOL 1.2b8.



### Basic TK GUI

[![](/images/d/dd/Tkgui1.png)](/index.php/File:Tkgui1.png)

The basic TK GUI. Menu options are File, Edit, Build, etc.

#### Movie Menu

This image shows the movie menu in PyMOL. 

[![](/images/d/d1/Tkgui2.png)](/index.php/File:Tkgui2.png)

The movie menu in PyMOL

  * **Append** just adds X seconds to the end of your movie. This is the equivalent of **[madd](/index.php/Madd "Madd") 1 x(FPS*time)**. So, you can very easily add anywhere from 1/4s to 60 seconds to your video. You can use this multiple times.



[![](/images/3/3e/Tkgui3.png)](/index.php/File:Tkgui3.png)

The Movie->Append submenu

. 

  * **Frame Rate** submenu that controls the frame rate of your movie in frames per second. It also controls whether or not the frame openGL screen frame rate is shown.



[![](/images/2/2f/Tkgui4.png)](/index.php/File:Tkgui4.png)

"Frame Rate" submenu

#### Scene Menu

For storing our scenes, we can issue the _[scene](/index.php/Scene "Scene") store_ command or **[mview](/index.php/Mview "Mview") store scene=sceneName** command. You can also use the mouse and store/recall/clear scenes. 

[![](/images/2/2e/Tkgui5.png)](/index.php/File:Tkgui5.png)

Scene Menu

##### Scene's Recall Submenu

This menu shows the first 12 scenes stored by F-key that you can recall (if they were stored). Please take note of the **[Buttons](/index.php/Scene_buttons "Scene buttons")** menu. This make a small menu of scenes in the lower left hand corner of your openGL screen. Very handy for roving through, moving around and switching to scenes. 

[![](/images/4/4b/Tkgui6.png)](/index.php/File:Tkgui6.png)

Scene->Recall submenu

#### Mouse Menu

There is a new mouse mode called [Three_Button_Motions](/index.php/Three_Button_Motions "Three Button Motions"). So we now have View->Edit-Motions as mouse modes. See [Config_Mouse](/index.php/Config_Mouse "Config Mouse"). 

[![](/images/7/79/Tkgui7.png)](/index.php/File:Tkgui7.png)

Mouse Menu

### openGL GUI

### Mouse Modes

If you've chosen **Mouse- >3 Button All** Modes from the Tk GUI you can now access the new mouse modes by clicking on the mouse info panel in the lower righ of the openGL window. Click until the modes cycle through Edit->View->Motions->Edit->... 

  * [![Editing mode pay careful attention to the button on top.](/images/1/1a/Glgui1.png)](/index.php/File:Glgui1.png "Editing mode pay careful attention to the button on top.")

Editing mode pay careful attention to the button on top. 

  * [![Viewing mode pay careful attention to the button on top.](/images/e/e2/Glgui2.png)](/index.php/File:Glgui2.png "Viewing mode pay careful attention to the button on top.")

Viewing mode pay careful attention to the button on top. 

  * [![Motions mode pay careful attention to the button on top.](/images/7/72/Glgui3.png)](/index.php/File:Glgui3.png "Motions mode pay careful attention to the button on top.")

Motions mode pay careful attention to the button on top. 




#### Motions Mode

In motions mode, if you click on **all- >M** you get options for storing, clearing, smoothing, interpolating, etc _the camera_ motions and positions. 

  * [![All->Motions shows a menu for Camera Motion/Position options](/images/e/e4/Glgui4.png)](/index.php/File:Glgui4.png "All->Motions shows a menu for Camera Motion/Position options")

All->Motions shows a menu for Camera Motion/Position options 

  * [![All->Motions->Store with Scene allows you to assign the chosen scene with the current frame.](/images/e/e1/Glgui5.png)](/index.php/File:Glgui5.png "All->Motions->Store with Scene allows you to assign the chosen scene with the current frame.")

All->Motions->Store with Scene allows you to assign the chosen scene with the current frame. 

  * [![All->Motions->Store with State allows you to assign a chosen state to this scene.](/images/f/f6/Glgui6.png)](/index.php/File:Glgui6.png "All->Motions->Store with State allows you to assign a chosen state to this scene.")

All->Motions->Store with State allows you to assign a chosen state to this scene. 

  * [![All->Motions->Smooth helps with smoothing interpolation between scenes](/images/c/ce/Glgui7.png)](/index.php/File:Glgui7.png "All->Motions->Smooth helps with smoothing interpolation between scenes")

All->Motions->Smooth helps with smoothing interpolation between scenes 

  * [![ALA->Motions menu shows our options, similar to the camera, but this time for the ALA fragment.](/images/4/49/Glgui8.png)](/index.php/File:Glgui8.png "ALA->Motions menu shows our options, similar to the camera, but this time for the ALA fragment.")

_ALA_ ->Motions menu shows our options, similar to the camera, but this time for the _ALA_ fragment. 




[ ← Previous Lesson](/index.php/MovieSchool_2 "MovieSchool 2") [ Next Lesson →](/index.php/MovieSchool_4 "MovieSchool 4")

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool_3&oldid=6993](https://pymolwiki.org/index.php?title=MovieSchool_3&oldid=6993)"


---

## MovieSchool 4

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Simple Movie Examples
    * 1.1 Multiple Zooming
    * 1.2 Animating an Alignment
    * 1.3 BB Inspector



## Simple Movie Examples

We now the ability to make some pretty simple, but cool movies. So, let's try a few. 

### Multiple Zooming

Let's try making a movie where we zoom into each ligand that's not water. In order to make this movie, I had to find a protein with suitable ligands, so you can do the same for your own protein. Just replace the hard-coded residue numbers. 

**Goal:** Make a movie that zoom into the three ligands, stays on that ligand for 2 seconds, then moves to the next. I also want smooth zoom out at the end. Don't let the length of this movie script throw you off, you've seen all of the movie commands and the initial commands are just loading the and making it look good. 
    
    
    # setup PyMOL for the movie
    reinitialize
    set matrix_mode, 1
    set movie_panel, 1
    
    # load the PDB, make selections for the ligands and
    # make the protein look snazzy.
    #load /spc/pdb/2jep.pdb
    fetch 2jep, async=0
    remove resn HOH
    orient
    select l1, c. A and i. 1397
    select l2, c. A and i. 1396
    select l3, c. B and i. 1396
    as cartoon
    color grey
    show sticks, het
    color magnesium, het
    
    # At 30 FPS this is then a 16 second movie.
    # We look at the structure for 2 seconds, zoom in to each ligand
    # and look at it for another 2 seconds, then, zoom out and look again
    # at everything for another 2 seconds.
    
    # initialize the 480 frame movie
    mset 1 x480
    
    # zoom all ('scene #1')
    frame 1
    mview store
    # stay here for 2 seconds
    frame 60
    mview store
    
    # zoom on ligand 1  ('scene #2')
    frame 120
    zoom l1
    mview store
    # stay here for 2 seconds
    frame 180
    mview store
    
    # zoom on ligand 2  ('scene #3')
    frame 240
    zoom l2
    mview store
    # stay for 2 seconds
    frame 300
    mview store
    
    # zoom to ligand 3  ('scene #4')
    frame 360
    zoom l3
    mview store
    # stay for 2 seconds
    frame 420
    mview store
    
    # zoom out  ('back to scene #1')
    frame 480
    zoom
    mview store
    
    # interpolate the frames
    mview reinterpolate
    
    # play the awesome movie!
    mplay
    
    # stop when you want
    # mstop
    

  


### Animating an Alignment
    
    
    # setup PyMOL for the movie
    reinitialize
    set matrix_mode, 1
    set movie_panel, 1
    
    # load the PDBs
    fetch 1cll 1ggz, async=0
    
    # orient the scene
    as cartoon
    orient
    
    # make 100-frame movie
    mset 1 x100
    # goto frame 1
    frame 1
    
    # store the camera position and object
    # positions in frame 1
    mview store
    mview store, object=1cll
    mview store, object=1ggz
    
    # goto frame 90
    frame 90
    # align the two proteins
    super 1cll, 1ggz
    # we rezoom to center the camera on the 
    # two aligned proteins
    zoom
    # store the camera positions
    mview store
    # store the new object position(s)
    mview store, object=1cll
    mview store, object=1ggz
    
    # have PyMOL stitch together the scenes.
    mview reinterpolate
    mview reinterpolate, object=1cll
    mview reinterpolate, object=1ggz
    
    # rewind
    frame 1
    # get some popcorn!  :-)
    mplay
    

### BB Inspector

This movie will walk down the alpha carbons inspecting each one for 1/3 of a second. :-) What would be cool would be to calculate the difference vector from i. n+1 to i. n and then walk that path. Anyhow, here you go: 
    
    
    # usual setup
    reinitialize
    set matrix_mode, 1
    set movie_panel, 1
    
    # fetch 1CLL to work on; this will only work on 1cll
    # or any other protein with 144 AAs starting at resi 4.
    fetch 1cll, async=0
    as lines, n. C+O+N+CA
    color marine
    zoom i. 4+5
    
    # 10 frames per AA
    mset 1 x1440
    mview store
    
    # this code maps the zooming of
    # one AA and it's neighbor to 10 frames
    python
    for x in range(0,144):
      cmd.frame((10*x)+1)
      cmd.zoom( "n. CA and i. " + str(x) + "+" + str(x+1))
      cmd.mview("store")
    python end
    
    # goto the end and interpolate all the frames
    frame 288
    mview store
    mview reinterpolate
    
    # I know, it's not smooth & cool like the other ones
    mplay
    

[ ← Previous Lesson](/index.php/MovieSchool_3 "MovieSchool 3") [ Next Lesson →](/index.php/MovieSchool_5 "MovieSchool 5")

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool_4&oldid=8025](https://pymolwiki.org/index.php?title=MovieSchool_4&oldid=8025)"


---

## MovieSchool 5

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Putting It All Together
    * 1.1 Camera Motions
      * 1.1.1 Camera Motions Movie Example
    * 1.2 Object Motions
    * 1.3 Camera & Object Motions
    * 1.4 Representations
    * 1.5 Motions & Representations
    * 1.6 Extras



## Putting It All Together

Now that we understand the basic PyMOL commands for movie making, we build a few ideas--which have already been hinted at--that lead to the final goal: allowing you to generate PyMOL movies to tell the stories you want to tell. We start simple and build up the complexity. Here's an outline of the ideas: 

  * **camera motions** —just moving the camera around your scene
  * **object motions** —keeping the camera still, but moving the objects around the scene
  * **camera & object motions **—moving both the camera and objects around
  * **representations** —changing representations in a movie (eg. sticks to cartoons to surface; hiding/showing, etc.)
  * **motions & representations **—adding all the motions and representations together
  * **extras** —pseudoatom labels; scene messages, etc.
  * **final example movies** —some examples combining the above knowledge



If you've read through most the above article on movie making, then these sections should be more of a review. There are some tricks in here that might be worth reading, however. 

Many of the movie scripts below assume that you have readied PyMOL for movie generation. To do that use the following code: 
    
    
    reinitialize
    set matrix_mode, 1
    set movie_panel, 1
    set scene_buttons, 1
    set cache_frames, 1
    config_mouse three_button_motions, 1
    
    # initialize a 100 frame movie
    mset 1 x100
    

### Camera Motions

One of our first movies above was a very simple zoom on an atom in an amino acid. The first 'scene' was the entire amino acid and the 2nd 'scene' was the zoomed in atom. We just connected the two scenes and asked PyMOL to make the transition between the two smooth. This is the idea of camera motions in PyMOL. (You may not know that when you click on a protein and rotate it or drag it in PyMOL you're actually moving the camera, not the protein.) 

Assuming you had readied PyMOL to make your movie to get a camera motion you do the following: 

  * **save the camera's first position information** : Once you have the protein/object aligned and shown in the representation of your choice, set the **camera** position information


    
    
    # goto the first frame
    frame 1
    # store the CAMERA positions ONLY
    mview store
    

  * **save the cameras second position information** : Now using the mouse (or scripted commands) move the camera to its new position--say zooming in on an important ligand or catalytic residue. Once that's done, tell PyMOL to store the new camera position in this frame:


    
    
    # goto the first frame
    frame 88
    # store the CAMERA positions ONLY
    mview store
    
    # now link the two stored camera positions together:
    mview reinterpolate
    
    # play your movie
    mplay
    

Hints: 

  * **mview reinterpolate, power=1" will turn off the smoothed starting and stopping of camera motions between the scenes. The smoothing gives a nicer feel to transitions. Try both, see which you prefer.**
  * check out mview's other options--like 'wrap'
  * using the three_button_motions option, you can pretty much make the entire movie w/your mouse: All->M->Store is the same as, _mview store_ , and All->M->Reinterpolate is the same as _mview reinterpolate_.
  * **[mview](/index.php/Mview "Mview") reset** unstores the camera position information for this frame.



#### Camera Motions Movie Example

This idea should be pretty sound at this point, but examples rock, so here's another. The pattern is: 

  * frame XYZ
  * script the view
  * mview store


    
    
    # setup PyMOL for movies
    reinitialize
    set matrix_mode, 1
    set movie_panel, 1
    set scene_buttons, 1
    set cache_frames, 1
    config_mouse three_button_motions, 1
    
    fetch 1te1, async=0
    extract AA, c. A
    extract BB, c. B
    color marine, AA
    color grey, BB
    as surface, BB
    as cartoon, AA
    
    mset 1 x620
    
    orient
    
    wizard message, "Can you see the blue protein inhibiting the gray protein?"
    
    frame 1
    mview store
    frame 30
    mview store
    
    ### cut below here and paste into script ###
    set_view (\
         0.307660401,    0.011366921,    0.951428533,\
         0.930296898,   -0.213488042,   -0.298277378,\
         0.199727684,    0.976880252,   -0.076255992,\
         0.000000000,    0.000000000, -196.781448364,\
        27.129878998,   68.309677124,   51.827075958,\
       155.143981934,  238.418914795,  -20.000000000 )
    ### cut above here and paste into script ###
    
    # slowly show the inhibition
    frame 120
    mview store
    
    # wait 3 seconds
    frame 180
    mview store
    
    # define the inhib as the binding loop
    select inhib, AA and i. 148-155
    select (none)
    
    # slowly zoom in
    frame 300
    zoom inhib
    mview store
    
    # stop a second
    frame 330
    mview store
    
    # look around the binding pocket
    frame 390
    turn y, 150
    mview store
    
    # wrap back more quickly...
    frame 420
    turn y, -150
    mview store
    
    # one more gratuitous view
    frame 500
    ### cut below here and paste into script ###
    set_view (\
         0.943371952,    0.309539229,   -0.119302809,\
        -0.044248745,   -0.239008784,   -0.970008850,\
        -0.328769624,    0.920357347,   -0.211777285,\
         0.000000000,    0.000000000,  -30.773454666,\
        35.418403625,   72.805625916,   52.437019348,\
        20.233829498,   41.313076019,  -20.000000000 )
    ### cut above here and paste into script ###
    mview store
    
    frame 560
    mview store
    
    mview reinterpolate
    mplay
    

### Object Motions

Now that we're experts at moving the PyMOL camera around, let's start moving objects around while keeping the camera steady. To do this **you must have matrix_mode set to 1** , otherwise PyMOL won't save your object's repositioning. 

Let's use the same proteins as from the above inhibitor example. This time, let's try to get a simple movie that shows one of the proteins and then have the other one fly in to do the inhibiting. 
    
    
    # setup PyMOL for movies
    reinitialize
    set matrix_mode, 1
    set movie_panel, 1
    set scene_buttons, 1
    set cache_frames, 1
    config_mouse three_button_motions, 1
    
    # download the complex and setup it up
    fetch 1te1, async=0
    extract AA, c. A
    extract BB, c. B
    color marine, AA
    color grey, BB
    as surface, BB
    as cartoon, AA
    
    # intialize the movie 
    mset 1 x410
    
    # orient the scene
    set_view (\
         0.423117876,    0.061672822,    0.903973043,\
         0.789699256,   -0.514252067,   -0.334546506,\
         0.444237947,    0.855418444,   -0.266292989,\
         0.000107866,   -0.000027858, -196.784057617,\
        28.171787262,   70.919288635,   52.095287323,\
       155.143981934,  238.418914795,  -20.000000000 )
    
    # move the inhibitor off the screeen
    translate [0,0,100], object=AA
    
    # first movie scene
    frame 1
    wizard message, "Let's watch the binder float it, while the camera doesn't move."
    mview store, object=AA
    mview store, object=BB
    
    # 2 second pause for the user to catch up
    frame 60
    mview store, object=AA
    mview store, object=BB
    
    frame 300
    # slide the inhibitor in from over the camera.  :-)
    translate [0,0,-100], object=AA
    mview store, object=AA
    mview interpolate, object=AA
    
    # store & wait 2 seconds...
    frame 360
    mview store, object=AA
    mview store, object=BB
    mview reinterpolate, object=AA
    mview reinterpolate, object=BB
    
    # 'explode' apart
    frame 380
    translate [-70, 70, 70], object=AA
    translate [70, -70, -70], object=BB
    mview store, object=AA
    mview store, object=BB
    mview reinterpolate, object=AA
    mview reinterpolate, object=BB
    
    mplay
    

Hints: 

  * Use the mouse to get the 'right orientation & zoom'. Then use [get_view](/index.php/Get_view "Get view") to get the view matrix. Finally, store that camera-position view matrix in your script. Works every time.
  * For object motions, the command **translate [-70,70,70], object=AA** would be the same as using the mouse and moving the object AA -70 units on the X-axis, +70 units on the Y and 70 on the Z. If you don't use the **object=** you will not get the desired effect.
  * For the above 'explosion' you can get quick motions by interpolating a large change over just a few frames.



### Camera & Object Motions

Now let's combine the above two sections into one movie that has both object and camera motions. This should be cool... 
    
    
    # setup PyMOL for movies
    reinitialize
    set movie_auto_interpolate, off
    set matrix_mode, 1
    set movie_panel, 1
    set scene_buttons, 1
    set cache_frames, 1
    config_mouse three_button_motions, 1
     
    # download the complex and set it up
    fetch 1te1, async=0
    extract AA, c. A
    extract BB, c. B
    color marine, AA
    color grey, BB
    as surface, BB
    as cartoon, AA
    
    mset 1 x120
    # overview of the scene
    frame 1
    mview store
    mview store, object=AA
    mview store, object=BB
    
    # zoom into the binding pocket- setting the view means
    # that this will be a camera motion from frames 1 to 120.
    frame 120
    set_view (\
         0.993863702,    0.110482253,   -0.005255031,\
         0.054543663,   -0.530888498,   -0.845684826,\
        -0.096224494,    0.840209842,   -0.533656776,\
         0.000000000,    0.000000000,  -50.366786957,\
        34.781314850,   71.208221436,   52.535022736,\
        39.709556580,   61.024017334,  -20.000000000 )
    mview store
    mview store, object=AA
    mview store, object=BB
    
    # wiggle the inhibitor, like it's trying to escape
    python
    for x in range(20):
      cmd.madd("1 x3"); cmd.frame(1000);
      cmd.rotate("x", 2.0, object="AA")
      cmd.mview("store", object="AA")
      cmd.mview("store")
      cmd.mview("interpolate", object="AA")
      cmd.mview("reinterpolate")
    
      cmd.madd("1 x3"); cmd.frame(1000);
      cmd.rotate("x", -2.0, object="AA")
      cmd.mview("store", object="AA")
      cmd.mview("store")
      cmd.mview("interpolate", object="AA")
      cmd.mview("reinterpolate")
    python end
    
    mview store
    
    madd 1 x60
    frame 240
    mview store
    mview store, object=AA
    mview store, object=BB
    
    mview interpolate, object=AA
    mview interpolate, object=BB
    mview reinterpolate
    
    mplay
    

### Representations

Scenes are the only way to change representations (eg sticks to cartoon). In PyMOL to show representation changes we need to have a list of scenes, that we then assign to a given frame. Once this is done we can reinterpolate through scenes to have beautifully smooth transitions. 

PyMOL makes it very easy to setup your scenes--for that _look_ \--and save them in a stack (and, now, even move them around). We covered scenes already in the tutorial, so please check that if you need more help on scenes, or see [Category:Scenes](/index.php/Category:Scenes "Category:Scenes") for more commands and hints. 

In PyMOL, to attach a scene to a frame (or, technically a frame to a scene) we simply do the following: 
    
    
    mview store, scene=sceneName
    

Let's take a look at a simple movie that changes the representation of some object. This will show a tryptophan going from lines to sticks and back. 
    
    
    set scene_buttons, 1
    set movie_panel, 1
    
    # make a 90 frame movie, all STATE 1.
    mset 1 x90
    
    # load a trypotphan fragment
    frag trp
    
    # Tell PyMOL to call this current view '001'.
    scene 001, store
    # goto frame 1 and store this scene & camera
    frame 1
    mview store, scene=001
    
    # setup the next 'view'
    as sticks
    scene 002, store
    
    # goto frame 60 and show sticks
    frame 45
    mview store, scene=002
    
    mview reinterpolate
    mplay
    

To show you how easy adding camera motions + representations is, simply insert after line 17 (_as sticks_) a new line 18 that only has 
    
    
    orient
    

Here's what's happening in the movie, above. Lines 1—2 set PyMOL up to show you [scene buttons](/index.php/Scene_buttons "Scene buttons") and the [movie panel](/index.php/Movie_panel "Movie panel"). Line 5 makes a 90 **frame** movie that only spans the first state. Line 8 makes the TRP fragment from PyMOL's stored knowledge of residues. Line 11 asks PyMOL to store this current 'view' as a [scene](/index.php/Scene "Scene") and call that scene _001_. Now that a scene is made, we need to associate with a frame. So, we goto frame 1 in line 13. In line 14 we actually make the scene-to-frame assignment with the [mview](/index.php/Mview "Mview") store command, specifying the mapping between scene 001 and frame 1. Lines 17 and 18 change the representation from lines to sticks and then orients the TRP residue (which causes the camera to move). In line 19, we save this new camera position and the sticks representation of TRP in the scene called _002_. Next, in line 22, we go somewhere ahead in the movie, here 1/2 way to the end. Line 23 assigns scene _002_ to frame 45. Lastly, because we changed the camera position, we need to reniterpolate it's movement, and that's done in line 24. As you know [mplay](/index.php/Mplay "Mplay") just plays the movie. 

  * **Hint** : Clear frames with **mview clear**
  * **Hint** : Clear scenes with **scene sceneName, delete** ; you can delete all scenes with **scene *, delete**.



  
Now that we have that, let's learn more complex tricks with motion and representations! 

### Motions & Representations

Let's look at the inhibitor complex again, (pdb 1TE1). This time, we would like to make a more compex movie that shows moving objects (not the camera) and changing representations. I do not want wrapping b/c I slide the molecule off the sceen. To ensure that, I allow 0 transition frames, by setting the last frame in the movie by **mview store**. 
    
    
    # b/c there will be motions, we need matrix_mode
    reinitialize
    set matrix_mode, 1
    set scene_buttons, 1
    set movie_panel, 1
    
    # setup the movie for 410 frames
    mset 1 x240
    
    # load the complex & set it up
    fetch 1te1, async=0
    as cartoon
    remove resn HOH
    extract AA, c. A
    extract BB, c. B
    
    # zoom in on the binding pocket
    ### cut below here and paste into script ###
    set_view (\
         0.478582859,   -0.269358903,    0.835705996,\
         0.805654645,   -0.243718535,   -0.539927721,\
         0.349110991,    0.931690335,    0.100370787,\
         0.000000000,    0.000000000,  -49.810981750,\
        35.418403625,   72.805625916,   52.437019348,\
        39.271358490,   60.350605011,  -20.000000000 )
    ### cut above here and paste into script ###
    
    # store object AA's position
    frame 1
    scene 001, store, message="This doesn't quite give us a good feel for what's going on."
    mview store, object=AA
    mview store, scene=001
    
    # change up the scene a bit
    frame 60
    color grey, AA
    color marine, BB
    scene 002, store, message="Recoloring helps a little."
    mview store, scene=002
    
    # show chain B as surface
    frame 120
    as surface, BB
    mview store, object=AA
    scene 003, store, message="Surface helps alot."
    mview store, 120, scene=003
    
    # show more
    frame 180
    select inhib, AA and i. 148-155; 
    show sticks, inhib
    color magenta, inhib
    select nn, BB within 7 of inhib; deselect;
    set transparency, 0.65, nn
    show sticks, nn
    set_bond stick_color, chartreuse, nn
    scene 004, store, message="Coloring, transparency, and other objects help more..."
    mview store, object=AA
    mview store, 180
    
    # move AA away
    frame 240
    translate [20, 10, 80], object=AA
    mview store, object=AA
    mview reinterpolate, object=AA
    mview store, 240
    mview reinterpolate
    

  


  * **Hint** : Use **mview store, frameNo, scene=sceneName** to save complex scenes. If you do:


    
    
    frame X
    scene Y
    mview store X, scene=Y
    

then the PyMOL GUI may not have caught up to the correct place before saving the scene information. 

  * **Hint** : Things can get confusing with frames, states and whatnot. A good idea is to ALWAYS set your frame when setting a new scene. If you don't **mview store** might save a different frame number than you think you're on. Or, you can save the frame number too with, **mview store, frameNo, scene=XXX, object=YYY**.
  * **Hint** : Use [deselect](/index.php/Deselect "Deselect") in a script to hide the dots from making a new selection.
  * **Hint** : Use [set_bond](/index.php/Set_bond "Set bond") to set bond properties, like colors, on [sticks](/index.php/Sticks "Sticks") representations.



### Extras

[ ← Previous Lesson](/index.php/MovieSchool_4 "MovieSchool 4") [ Next Lesson →](/index.php/MovieSchool_6 "MovieSchool 6")

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool_5&oldid=9503](https://pymolwiki.org/index.php?title=MovieSchool_5&oldid=9503)"


---

## MovieSchool 6

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Exporting your Movie

This wiki already has lots of information on how to convert your PyMOL movie to another format. Check out [those](/index.php/Making_Movies "Making Movies") [pages](/index.php/Software_Codecs "Software Codecs"). 

Once you've setup your movie as in any of the previous examples, you have a couple options for making a movie. 

**Export from PyMOL** : Newer PyMOLs support 

  * File→Save Movie As→MPEG
  * File→Save Movie As→PNG Images



from the menu. 

**mpeg** : This needs a bit of work. The following instructions work with trunk (1.2r2pre). Compile pymol trunk. 

Checkout [freemol](http://freemol.org) source from 
    
    
    svn co svn://bioinformatics.org/svnroot/freemol/trunk freemol-trunk
    

This creates a directory called **freemol-trunk**. Inside **freemol-trunk** is a directory called **freemol**. Now you need to setenv (or export) a variable called FREEMOL pointing to the above folder - **freemol**). In TCSH and related C-shells: 
    
    
    setenv FREEMOL /my/long/path/to/freemol-trunk/freemol
    cd freemol-trunk/src/mpeg_encode
    

or in bash: 
    
    
    export FREEMOL=/my/long/path/to/freemol-trunk/freemol
    cd freemol-trunk/src/mpeg_encode
    

And then 
    
    
    cd freemol-trunk/src/mpeg_encode
    

(Ignore all other directories inside directory **src**. These were unnecessary, at least for me.) 
    
    
    ./configure
    make
    make install
    

Add these two lines 
    
    
    FREEMOL=/my/long/path/to/freemol-trunk/freemol
    export FREEMOL
    

to your pymol script (_i.e._ /sw/bin/pymol). Launch pymol and **Save Movie as MPEG** should work now. No need of any complicated codecs. 

**mpng** : You can still use the good old [mpng](/index.php/Mpng "Mpng") option to save all your frames to disk. You can then compile them into a MPEG (see below). 

**Old Style** : One of the older scripting styles was to make minor changes and dump PNGs. This is essentially obviated with PyMOL's new movie-making functionality. The **old style** was to simply call cmd.png every time you made a scene change. 

Hints: 

  * Movie not ray traced? Make sure you set ray_trace_frames to 1.



### Codecs

See [Software_Codecs](/index.php/Software_Codecs "Software Codecs") for information on how to stitch together movies from PNGs and optimize them for great crisp-looking movies. 

  
[ ← Previous Lesson](/index.php/MovieSchool_5 "MovieSchool 5")

[ Back to Start](/index.php/MovieSchool "MovieSchool")

Retrieved from "[https://pymolwiki.org/index.php?title=MovieSchool_6&oldid=11133](https://pymolwiki.org/index.php?title=MovieSchool_6&oldid=11133)"


---

## Mutagenesis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Mutagenesis

PyMol has a **Mutagenesis Wizard** to make mutagenesis very easy for the end user. 

In [rotkit](/index.php/Rotkit "Rotkit"), a function has been made to call a mutagenesis. 

As of PyMOL version [2.2](https://pymol.org/dokuwiki/?id=media:new22), users may now perform base mutations in nucleotide chains. 

## Walk-through

To mutate a residue follow these easy steps: 

  1. Load a PDB file  
  


[![Mutag1.png](/images/7/78/Mutag1.png)](/index.php/File:Mutag1.png)

  2. Under the **Wizard** menu select **Mutagenesis**  
  


[![Mutag2.png](/images/1/19/Mutag2.png)](/index.php/File:Mutag2.png)

  3. In the PyMol viewer window select a residue  
  


[![Mutag3.png](/images/1/19/Mutag3.png)](/index.php/File:Mutag3.png)

  4. Select **No Mutation** and select resultant residue  
  


[![Mutag4.png](/images/b/bb/Mutag4.png)](/index.php/File:Mutag4.png)

  5. Selecting the rotamer you think better fits your structure.   
  
Several side chain orientations (rotamers) are possible. You can use the back-and-forth movie controls (lower right corner) to display (in white) each of the rotamers available for this residue in PyMOL, whose current and total numbers are shown in the (green) Frame info. The rotamers are ordered according to their frequencies of occurrence in proteins, shown as a red percentage at the mutation object, which exists while mutagenesis is being performed.   
  

  6. Select Apply 
  7. Select Done 



## Explanation of colour codes

From [[1]](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg06196.html): 

The visible disks & colors in the Mutagenesis Wizard indicate pairwise overlap of atomic van der Waals radii. The intent is to provide a qualitative feedback regarding contacts and bumps. 

Short green lines or small green disks are shown when atoms are almost in contact or slightly overlapping. Large red disks indicate signficant van der Waals overlap. Everything else lies between those extremes. 

Retrieved from "[https://pymolwiki.org/index.php?title=Mutagenesis&oldid=12784](https://pymolwiki.org/index.php?title=Mutagenesis&oldid=12784)"


---

## Plugins Tutorial

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This tutorial is about writing your own plugin for PyMOL 2.x. 

See the [plugins](/index.php/Plugins "Plugins") page for how to install and use exisiting plugins. 

## Contents

  * 1 Writing Plugins: Learn By Example
    * 1.1 Plugin Files
    * 1.2 Registering your Plugin
    * 1.3 Creating a GUI
    * 1.4 Filling the GUI with Widgets
    * 1.5 Make Buttons do something
    * 1.6 Deploy the final plugin
    * 1.7 Full Source
  * 2 Extending Plugins to the Command Line
  * 3 See Also



## Writing Plugins: Learn By Example

This tutorial shows how to write a PyMOL plugin with PyQt. **The full source of the demo plugin is[available on github](https://github.com/Pymol-Scripts/pymol2-demo-plugin)**. 

The demo plugin adds a dialog to render images at a custom resolution. 

### Plugin Files

We will create a plugin which consists of [multiple files inside a directory](/index.php/PluginArchitecture#More_than_one_file_per_plugin "PluginArchitecture"): 
    
    
    pymol2-demo-plugin/
    ├── demowidget.ui
    └── __init__.py (defines "__init_plugin__(app=None)" function)
    

### Registering your Plugin

First you must add your plugin to the _Plugins_ menu. This is done in the `__init_plugin__` function of your plugin. A callback (here: `run_plugin_gui`) is added. 
    
    
    # file pymol2-demo-plugin/__init__.py
    def __init_plugin__(app=None):
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Demo "Render" Plugin', run_plugin_gui)
    

_Legacy note_ : The `__init_plugin__` function takes one argument, a reference to the main Tk app to [support legacy plugins witten with Tkinter](/index.php/PluginArchitecture#init_plugin "PluginArchitecture") (unused with PyQt plugins). 

### Creating a GUI

The callback may do arbitrary stuff. Here we're going to create a dialog and show it to the user. 
    
    
    # global reference to avoid garbage collection of our dialog
    dialog = None
    
    def run_plugin_gui():
        # pymol.Qt provides the PyQt5 interface, but may support PyQt4 and/or PySide as well
        from pymol.Qt import QtWidgets
    
        global dialog
    
        if dialog is None:
            # create a new (empty) Window
            dialog = QtWidgets.QDialog()
    
            # TODO: FILL DIALOG WITH WIDGETS HERE
    
        dialog.show()
    

### Filling the GUI with Widgets

[![](/images/8/81/Designer-demo-plugin.png)](/index.php/File:Designer-demo-plugin.png)

[](/index.php/File:Designer-demo-plugin.png "Enlarge")

Screenshot of Qt Designer

The most convenient and maintainable way to create user interfaces with Qt is by using the [Qt Designer](http://doc.qt.io/qt-5/qtdesigner-manual.html). Follow these steps to get started: 

  1. Create a new "Widget" form (New Form > templates/forms > Widget > Create)
  2. Drag widgets like "Push Button" or "Line Edit" into your form
  3. Name your widgets ("objectName" in the "Property Editor"), this name will become a Python variable in your code
  4. Save as a UI file (`*.ui`) inside the `pymol2-demo-plugin` directory



PyMOL provides a utility function to load a UI file into a parent widget: `pymol.Qt.utils.loadUi(filename, widget)`
    
    
        # filename of our UI file
        uifile = os.path.join(os.path.dirname(__file__), 'demowidget.ui')
    
        # load the UI file into our dialog
        from pymol.Qt.utils import loadUi
        form = loadUi(uifile, dialog)
    

### Make Buttons do something

We need to connect the [signals](http://doc.qt.io/qt-5/signalsandslots.html) of those widgets we created with the Qt Designer, like our buttons' `clicked` signal. Our example has a "Ray" button (_objectName=button_ray_) that we'll connect to the `run()` function. 
    
    
        def run():
            # get form data
            height = form.input_height.value()
    
            # some debugging feedback
            print('User entered height', height)
    
            # TODO: DO SOMETHING WITH FORM DATA
    
        form.button_ray.clicked.connect(run)
    

### Deploy the final plugin

The `pymol2-demo-plugin` directory can be zipped for deployment (see [PluginArchitecture#More than one file per plugin](/index.php/PluginArchitecture#More_than_one_file_per_plugin "PluginArchitecture")). 

### Full Source

The full source of the demo plugin is [available on github](https://github.com/Pymol-Scripts/pymol2-demo-plugin)**.**

## Extending Plugins to the Command Line

While not really applicable to our simple "Render" plugin (it only wraps existing `ray` and `png` commands), most plugins should also provide their functionality as new PyMOL commands using [cmd.extend()](/index.php/Extend "Extend"). 

## See Also

  * [PluginArchitecture](/index.php/PluginArchitecture "PluginArchitecture")
  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins](/index.php/Plugins "Plugins")
  * [cmd.extend()](/index.php/Extend "Extend")



Retrieved from "[https://pymolwiki.org/index.php?title=Plugins_Tutorial&oldid=13173](https://pymolwiki.org/index.php?title=Plugins_Tutorial&oldid=13173)"


---

## Practical Pymol for Beginners

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Although PyMOL has a powerful and flexible interface, it is complex, and can appear daunting to new users. This guide is intended to introduce the PyMOL interface and basic tasks without leaving the mouse behind. 

## Contents

  * 1 The PyMOL Interface
    * 1.1 About the command line
  * 2 Getting started: explore a protein
    * 2.1 Related Commands
  * 3 What else can this thing do?
    * 3.1 Saving an image
    * 3.2 Selecting parts of an object
    * 3.3 Whoops: Getting unstuck
    * 3.4 Related Commands
    * 3.5 Saving work
      * 3.5.1 Sessions
      * 3.5.2 Molecules
      * 3.5.3 Images
      * 3.5.4 Movies
    * 3.6 Scripting
      * 3.6.1 Introduction and Very Simple Scripting
      * 3.6.2 The Python MiniShell
      * 3.6.3 Learning More...
  * 4 Coming Soon
    * 4.1 Logging



## The PyMOL Interface

When PyMOL is opened, two windows appear. The smaller window (called the "External GUI" in PyMOL documentation) contains the menu bar (**File** , **Edit** , **Help** , **Display** , etc), shortcut buttons for common commands, and the command line. 

[![](/images/2/20/Viewer_guide.png)](/index.php/File:Viewer_guide.png)

[](/index.php/File:Viewer_guide.png "Enlarge")

The PyMol Viewer Window

The second window is the PyMOL Viewer, which is where all the magic happens. In the Viewer, 3D models are displayed, and the user interacts (eg rotates) and manipulates the model. 

The objects that PyMOL renders in 3D are loaded from coordinate files that describe (in great detail) locations of individual atoms in the molecule. PyMOL can display more than one object at a time, and provides an Object Control Panel to adjust viewing modes, colors, labels, hiding, and just about anything else relating to objects. After each object name is a set of command buttons which control the object. Here are the buttons and some of their options: 

  * **A** \- _Actions_ : Rename, duplicate, remove, apply presets (like "ball-and-stick" or "publication"), perform computations
  * **S** \- _Show_ : Change the way things appear, eg change to stick or cartoon view.
  * **H** \- _Hide_ : Things that are shown using **S** accumulate, and don't automatically replace the last view. **H** is the opposite of **S** and hides unwanted representations.
  * **L** \- _Label_ : Label atoms, residues, etc.
  * **C** \- _Color_ : Change the color of atoms and groups.



The lower-right corner of the Viewer contains a guide to using the mouse, as well as a powerful selection tool. There is also another command line at the bottom of the Viewer (**PyMOL >**). 

### About the command line

The PyMOL command line is a great tool that lets the experienced user change all sorts of options that simply don't appear in the point-and-click graphical interface. It can also be a lot faster. Combined with scripting, it is a powerful option for automating tasks and making intricate sets of changes. But, it's complex, and page upon page of PyMOL documentation cover these commands, so we're going to ignore them as much as possible. 

Although this guide may include some text commands and links to more advanced documentation, they're purely optional and meant to be informative. 

To run any text command, type it in at a **PyMOL >** command line and hit _[Enter]_. 

## Getting started: explore a protein

PyMOL is great for casual visualization of biological molecules. In this example, a PDB file describing a protein is loaded and its style and color are tweaked. 

[![](/images/7/75/Kchannel-rainbow.png)](/index.php/File:Kchannel-rainbow.png)

[](/index.php/File:Kchannel-rainbow.png "Enlarge")

The end result will look something like this

[![](/images/5/54/Mouse-3view.png)](/index.php/File:Mouse-3view.png)

[](/index.php/File:Mouse-3view.png "Enlarge")

Default buttons for viewing with a 3-button mouse

  1. Obtain a PDB coordinates file for your favorite protein. (The [RCSB Protein Data Bank](http://www.pdb.org/) is a public structure repository containing over 40,000 protein structures in PDB format available for download, not a bad place to look.) For this example, we're using the potassium channel from _Streptomyces Lividans_ ([1BL8](http://www.pdb.org/pdb/files/1bl8.pdb)). 
     * or just type, 
           
           fetch 1bl8
           

and you can skip the next step (see [Fetch](/index.php/Fetch "Fetch") command), as PyMOL will open the file for you.
  2. Open the PDB file using **File** => **Open...** from the menu bar. The protein's structure will appear, probably rendered as simple bonding lines.
  3. The right side of the Viewer shows the loaded PDB as an object, as well as its command buttons. Each button contains a submenu with more options. Click **S** , then **cartoon** to show the protein's secondary structure in popular cartoon form. 
     * Notice that the lines view is still visible on top of the cartoon view. To hide the lines, click **H** then **lines**.
  4. To change the color of each protein chain (as defined in the coordinate file), click **C** then select **chainbows** from the **by chain** menu. "Chainbows" colors residues in each protein chain as a rainbow that begins with blue and ends with red. 
     * Another common coloring method assigns a single color to each chain. Click **C** then select **by chain** from the **by chain** menu.
  5. Click and drag the protein to change the view. A list of mouse buttons is below the object control panel. 
     * Rota: Rotate
     * Move: Translate object along an axis
     * MoveZ: aka Zoom
     * Sele: Select
     * Slab:
     * Cent:
     * PkAt:



### Related Commands

[Load](/index.php/Load "Load"), [Fetch](/index.php/Fetch "Fetch"), [Color](/index.php/Color "Color"), [Show](/index.php/Show "Show"), [Show_as](/index.php/Show_as "Show as"), [Cartoon](/index.php/Cartoon "Cartoon"), [Lines](/index.php/Lines "Lines"), [Rotate](/index.php/Rotate "Rotate"), [Select](/index.php/Select "Select"), [Center](/index.php/Center "Center")

## What else can this thing do?

So, now what? Good question. PyMOL is a powerful program, and everyone uses it for something different. The remainder of this guide is devoted to common tasks that come in handy. 

### Saving an image

You've found the perfect view, and you'd like to [Save](/index.php/Save "Save") it? 

  1. Change the size of the viewer and zoom in/out to make the image the right size. Images are saved exactly as they appear in the viewer, and the size of the viewer determines the size of the image. 
     * For example, images for printing and presenting should be larger than images for a posting on a website.
  2. Because PyMOL was designed to run on older computers, the maximum quality is not enabled by default. To change this, select **Display** , **Quality** , **Maximum Quality**. Notice that round things are rounder, curves are curvier, and color shading is smoother.
  3. Once you've found the appropriate view, save an image using **File** , **Save Image...** An picture of the current view will be saved in PNG format.



**Tip:** Using the _[ray](/index.php/Ray "Ray")_ command before saving an image will create a higher quality version with shadows, etc. This can take time, depending on the size of the image and speed of the computer, but the images created after ray tracing are usually spectacular. However, the ray tracing disappears if the view is changed in any way. 

### Selecting parts of an object

Sometimes it might be useful to select part of an object and modify it directly; for example, selecting active-site residues in a protein to highlight them in another color. 

In the lower-right corner of the Viewer, next to the animation controls, is an **S** button (not to be confused with the **S** how button in the object control panel) that activates the selection tool. The selection tool can be changed to select atoms or residues by clicking _Selecting Residues(or whatever)_ until the right mode appears. 

Once selecting is activated, a list of parts to select appears at the top of the Viewer. Select things clicking or dragging across a range. 

Selections can be controlled individually in the object control panel, just like any other object. To save a selection, select **rename** from the **A** menu. 

### Whoops: Getting unstuck

PyMOL is a program meant to be explored and played with, and sometimes funny things happen in the process. A few common problems: 

  * _The model disappeared:_ Sometimes while rotating and moving a model, it can get lost. Right-click the background of the viewer, and select **reset** from the _Main Pop-Up_. The model should return to view; if it doesn't, make sure the object is being drawn using the **S** menu.
  * _The model has funny colors, labels, etc and they won't go away:_ The **H** menu in the object control panel will remove unwanted details; however, sometimes it's difficult to know exactly what to remove. Select **H** , then **everything** to hide all details and start fresh.
  * _Things are really messed up:_ Use **File** , **Reinitalize** to reset PyMOL to its initial state, but all work will be lost.



### Related Commands

[Save](/index.php/Save "Save"), [Viewport](/index.php/Viewport "Viewport"), [Zoom](/index.php/Zoom "Zoom"), [Save](/index.php/Save "Save"), [Ray](/index.php/Ray "Ray"), [Select](/index.php/Select "Select")

### Saving work

PyMOL supports saving your work in various formats. You can save, images, molecules, sessions, movies, etc. 

#### Sessions

A PyMOL sessions retains the state of your PyMOL instance. You can save the session to disk and reload it later. You can setup a complicated scene, with transitions and more, and simply save it as a PyMOL Session (.pse) file. Use, **File** =>**Save Session** or **Save Session As...**. 

Loading sessions is just as easy: **File** =>**Load** and choose your session file. 

#### Molecules

You can save a molecule by selecting **File** =>**Save Molecule**. A dialog box will appear from which you can select your molecule to save. You can also save an object or selection using the [Save](/index.php/Save "Save") command. It's very easy: 
    
    
    save  fileName, objSel
    

where fileName is something like "1abc.pdb", and objSel can be any object or selection. For example, to save just the waters to disk do: 
    
    
    save wat.pdb, resn HOH
    

#### Images

You can save images that you've rendered (with [Ray](/index.php/Ray "Ray")) or drawn (with [Draw](/index.php/Draw "Draw")) again using the [Save](/index.php/Save "Save") command or by **File** =>**Save Image**. You can save in [Png](/index.php/Png "Png"), VRML-2 and the POVRay formats. 

You can save images to disk, through the command line, by using the [Png](/index.php/Png "Png") command. 

#### Movies

PyMOL allows you to save movies you've created, too. You can automatically save the MPEG or save a series of PNG files--for stitching together later. This is a new command, and I don't know too much about it. Use **File** =>**Save Movie**. 

### Scripting

#### Introduction and Very Simple Scripting

Scripting in PyMOL ranges from very simple to very intricate. You can make a simple loop (say rotating a molecule one-degree at a time) or execute full featured scripts. Once you've got the basics of PyMOL down, learning scripting will greatly enhance your productivity and capabilities. 

Because of the way PyMOL is built on top of the Python interpreter, any command that PyMOL doesn't recognize it passes on to Python to execute. This is a **very** handy feature--you essentially have a live Python interpreter you can interact with, which makes your life much easier. Take the following for example: 
    
    
    f = 10.
    for x in range(0,100,10):
      cmd.set("spec_direct_power", float(float(x) / f))
      cmd.png("spec_dir_power" + str(x) + ".png", ray=1)
    

This simple script of 4 lines will create 10 images each one with the [Spec_direct_power](/index.php/Spec_direct_power "Spec direct power") changed (see the [Spec_direct_power](/index.php/Spec_direct_power "Spec direct power") for the output of this script; the animated GIF). Did you notice that **f** on line 1 and **for** and **x** on line 2 are not PyMOL commands or symbols? You can simply write Python code that interacts with PyMOL. Brilliant! 

#### The Python MiniShell

Taking this one level higher, you can write little code snippets, like 10-20+ lines of code that perform some specific task and wrap these in the `python` and `python end` commands. (If your script ever makes a mistake, you can abort the endeavor with `end` instead of `python end`. The power here is that none of your command are executed until you type `python end`. Confused? Here's the above example in using the wrapper: 
    
    
    python
    f = 10.
    for x in range(0,100,10):
      cmd.set("spec_direct_power", float(float(x) / f))
      cmd.png("spec_dir_power" + str(x) + ".png", ray=1)
    python end
    

The `[python](/index.php/Python "Python")` command gives you complete access to the Python shell and `python end` brings you back into PyMOL's shell. Also, note that PyMOL saves information across instantiations of the `python` command. For example, 
    
    
    # enter mini python shell
    python
    ff = 10.
    python end
    
    # now we're back in the normal PyMOL shell, where PyMOL knows about the value
    print(ff)
    
    # in the mini shell, Python still knows about ff.
    python 
    print(ff)
    python end
    

#### Learning More...

To learn more about scripting check out: 

  * [Biochemistry_student_intro](/index.php/Biochemistry_student_intro "Biochemistry student intro") Basic use of GUI and script
  * [Simple_Scripting](/index.php/Simple_Scripting "Simple Scripting") introduction
  * [Advanced_Scripting](/index.php/Advanced_Scripting "Advanced Scripting") pages
  * [Popular Script Library](/index.php/Category:Script_Library "Category:Script Library").



## Coming Soon

### Logging

Retrieved from "[https://pymolwiki.org/index.php?title=Practical_Pymol_for_Beginners&oldid=12744](https://pymolwiki.org/index.php?title=Practical_Pymol_for_Beginners&oldid=12744)"


---

## Running Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page is a short description for beginners on how to run scripts in PyMOL, e.g. from the [Script Library](/index.php/Category:Script_Library "Category:Script Library"). For help on writing scripts, look here: [Simple_Scripting](/index.php/Simple_Scripting "Simple Scripting"), [Advanced_Scripting](/index.php/Advanced_Scripting "Advanced Scripting"). 

## Contents

  * 1 Saving the source code
  * 2 Python or Pymol?
  * 3 Python Modules
  * 4 Example: Color_Objects



## Saving the source code

First, you have to save the source code. Copy the text in the box to a text editor and save it in text format under an arbitrary file name, let's say we save it as script.txt in the PyMOL working directory. 

## Python or Pymol?

Then, you have to find out wheter the script you want to run is a python script or a PyMOL script. Look at the source code: if the lines look just as you type them in the command line in PyMOL, it is a [PyMOL script](/index.php/PML "PML"). Run it with @, in our example, type 
    
    
    @script.txt
    

in the PyMOL command line. You can find examples in the script library: [Split_Movement](/index.php/Split_Movement "Split Movement") (loads a structure from the internet and moves parts in different directions), [Show_charged](/index.php/Show_charged "Show charged") (Selects charged residues and highlights them). Any PyMOL log file would also be an example for a pymol script. 

If, in contrast, you find words as "import" in one of the first lines, or "def" or "cmd" in the text body, you have a python script. Run it with run. In the example, type 
    
    
    run script.txt
    

in the PyMOL command line. 

Most python scripts in the script library don't start action immediately but define a new command instead. You first have to run the script and then you can use the command. Many script pages provide hints for usage. If not, look at the bottom of the script for a line like this: 
    
    
    cmd.extend("new_command",new_command)
    

The text in quotation marks is the new command. Enter it in the PyMOL command line. 

You can find many examples for python scripts in the script library, e.g.: [Color_Objects](/index.php/Color_Objects "Color Objects") (colors all the objects differently) [Resicolor](/index.php/Resicolor "Resicolor") (Colors residues according to their property) 

## Python Modules

A python **module** is a python script that runs in it's own namespace, by using the **import** syntax instead of **run**. The file must be located in any directory of the [sys.path](http://docs.python.org/library/sys.html#sys.path) variable. This is the recommended way for scripts from the [Pymol-script-repo](/index.php/Git_intro "Git intro"). 
    
    
    import color_obj   # skip the .py extension!
    

## Example: Color_Objects
    
    
    PyMOL>run color_obj.py
    PyMOL>color_obj 
     
    Colouring objects using PyMOL defined colours
     
       obj_0001 red
       obj_0002 green
       obj_0003 blue
       obj_0004 yellow
       obj_0005 violet
       obj_0006 cyan
       obj_0007 salmon
       obj_0008 lime
       obj_0009 pink
       obj_0010 slate
       obj_0011 magenta
       obj_0012 orange
       obj_0013 marine
       obj_0014 olive
       ...
    

Retrieved from "[https://pymolwiki.org/index.php?title=Running_Scripts&oldid=13250](https://pymolwiki.org/index.php?title=Running_Scripts&oldid=13250)"


---

## Simple Scripting

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

_This page discusses writing your own simple scripts for PyMOL. If you're looking to download scripts try the[Script Library](/index.php/Category:Script_Library "Category:Script Library"). If you need help with running a script try [Running_Scripts](/index.php/Running_Scripts "Running Scripts") _

One of the more powerful features of PyMOL is that it supports Python scripting. This gives you the power to use most of the [Python libraries](http://docs.python.org/api/api.html) to write programs and then send the results back into PyMOL. Some useful extensions to PyMOL can be found in our [Script Library](/index.php/Category:Script_Library "Category:Script Library"). 

## Contents

  * 1 General Scripts
    * 1.1 Getting PyMOL Data into your Script
      * 1.1.1 Getting Data From your Script into PyMOL
    * 1.2 Example
    * 1.3 Basic Script Body
  * 2 Python Version Supported by PyMOL



# General Scripts

General PyMOL scripting is done in Python. It's really quite simple, just write your function (following a couple simple rules) and then let PyMOL know about it by using the **[cmd.extend](/index.php/Extend "Extend")** command. Here's the simple recipe for writing your own simple scripts for PyMOL: 

**To write them** : 

  1. Write the function, let's call it **doSimpleThing** , in a Python file, let's call the file **pyProgram.py**.
  2. Add the following command to the end of the **pyProgram.py** file 
         
         cmd.extend("doSimpleThing",doSimpleThing)
         




**To use them** : 

  1. simply import the script into PyMOL: 
         
         run /home/userName/path/toscript/pyProgram.py
         

  2. Then, just type the name of the command: _doSimpleThing_ and pass any needed arguments.



That's it. Your script can, through Python, import any modules you need and also edit modify objects in PyMOL. 

## Getting PyMOL Data into your Script

To get PyMOL data into your script you will need to somehow get access to the PyMOL objects and pull out the data. For example, if you want the atomic coordinates of a selection of alpha carbon atoms your Python function may do something like this (see also [iterate_state](/index.php/Iterate_state "Iterate state")): 
    
    
    # Import PyMOL's stored module.  This will allow us with a 
    # way to pull out the PyMOL data and modify it in our script.
    # See below.
    from pymol import stored
    
    def functionName( userSelection ):
        # this array will be used to hold the coordinates.  It
        # has access to PyMOL objects and, we have access to it.
        stored.alphaCarbons = []
    
        # let's just get the alpha carbons, so make the
        # selection just for them
        userSelection = userSelection + " and n. CA"
    
        # iterate over state 1, or the userSelection -- this just means
        # for each item in the selection do what the next parameter says.
        # And, that is to append the (x,y,z) coordinates to the stored.alphaCarbon
        # array.
        cmd.iterate_state(1, selector.process(userSelection), "stored.alphaCarbons.append([x,y,z])")
    
        # stored.alphaCarbons now has the data you want.
    
        ... do something to your coordinates ...
    

### Getting Data From your Script into PyMOL

Usually this step is easier. To get your data into PyMOL, it's usually through modifying some object, rotating a molecule, for example. To do that, you can use the [alter](/index.php/Alter "Alter") or [alter_state](/index.php?title=Alter_state&action=edit&redlink=1 "Alter state \(page does not exist\)") commands. Let's say for example, that we have translated the molecular coordinates from the last example by some vector (we moved the alpha carbons). Now, we want to make the change and see it in PyMOL. To write the coordinates back we do: 
    
    
    # we need to know which PyMOL object to modify.  There could be many molecules and objects
    # in the session, and we don't want to ruin them.  The following line, gets the object
    # name from PyMOL
    objName = cmd.identify(sel2,1)[0][0]
    
    # Now, we alter each (x,y,z) array for the object, by popping out the values
    # in stored.alphaCarbons.  PyMOL should now reflect the changed coordinates.
    cmd.alter_state(1,objName,"(x,y,z)=stored.alphaCarbons.pop(0)")
    

## Example

Here's a script I wrote for [cealign](/index.php/Cealign "Cealign"). It takes two selections **of equal length** and computes the optimal overlap, and aligns them. See [Kabsch](/index.php/Kabsch "Kabsch") for the original code. Because this tutorial is for scripting and not optimal superposition, the original comments have been removed. 
    
    
    def optAlign( sel1, sel2 ):
            """
            @param sel1: First PyMol selection with N-atoms
            @param sel2: Second PyMol selection with N-atoms
            """
    
            # make the lists for holding coordinates
            # partial lists
            stored.sel1 = []
            stored.sel2 = []
            # full lists
            stored.mol1 = []
            stored.mol2 = []
    
            # -- CUT HERE
            sel1 = sel1 + " and N. CA"
            sel2 = sel2 + " and N. CA"
            # -- CUT HERE
    
            # This gets the coordinates from the PyMOL objects
            cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
            cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
    
            # ...begin math that does stuff to the coordinates...
            mol1 = cmd.identify(sel1,1)[0][0]
            mol2 = cmd.identify(sel2,1)[0][0]
            cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
            cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
            assert( len(stored.sel1) == len(stored.sel2))
            L = len(stored.sel1)
            assert( L > 0 )
            COM1 = numpy.sum(stored.sel1,axis=0) / float(L)
            COM2 = numpy.sum(stored.sel2,axis=0) / float(L)
            stored.sel1 = stored.sel1 - COM1
            stored.sel2 = stored.sel2 - COM2
            E0 = numpy.sum( numpy.sum(stored.sel1 * stored.sel1,axis=0),axis=0) + numpy.sum( numpy.sum(stored.sel2 * stored.sel2,axis=0)
    ,axis=0)
            reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
            if reflect == -1.0:
                    S[-1] = -S[-1]
                    V[:,-1] = -V[:,-1]
            RMSD = E0 - (2.0 * sum(S))
            RMSD = numpy.sqrt(abs(RMSD / L))
            U = numpy.dot(V, Wt)
            # ...end math that does stuff to the coordinates...
    
            # update the _array_ of coordinates; not PyMOL the coords in the PyMOL object
            stored.sel2 = numpy.dot((stored.mol2 - COM2), U) + COM1
            stored.sel2 = stored.sel2.tolist()
    
            # This updates PyMOL.  It is removing the elements in 
            # stored.sel2 and putting them into the (x,y,z) coordinates
            # of mol2.
            cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
    
            print "RMSD=%f" % RMSD
    
            cmd.orient(sel1 + " and " + sel2)
    
    # The extend command makes this runnable as a command, from PyMOL.
    cmd.extend("optAlign", optAlign)
    

## Basic Script Body

Want an easy block of working code to start your function from? Just copy/paste the following into your Python editor and get going! 
    
    
    #
    # -- basicCodeBlock.py
    #
    from pymol import cmd, stored
    
    def yourFunction( arg1, arg2 ):
        '''
    DESCRIPTION
    
        Brief description what this function does goes here
        '''
        #
        # Your code goes here
        #
        print "Hello, PyMOLers"
        print "You passed in %s and %s" % (arg1, arg2)
        print "I will return them to you in a list.  Here you go."
        return (arg1, arg2)
    
    cmd.extend( "yourFunction", yourFunction );
    

# Python Version Supported by PyMOL

Scripts used within PyMOL can only be written using the current version of Python that is supported by your version of PyMOL. To determine which version of Python you can use, type the following command into PyMOL: 
    
    
    print sys.version
    

Note that this version of Python is not necessarily related to the version that you may have installed on your system. 

This command can also be used to ensure that code you are distributing can be supported by the user's system. 

Retrieved from "[https://pymolwiki.org/index.php?title=Simple_Scripting&oldid=9103](https://pymolwiki.org/index.php?title=Simple_Scripting&oldid=9103)"


---

## Visualizing a computed structure - a commented example

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

  1. Obtain the [File.xyz.tar](/images/7/71/File.xyz.tar "File.xyz.tar") and untar it. It contains the coordinates, which were originally published in [J. Am. Chem. Soc.](http://pubs.acs.org/journals/jacsat/index.html) in the Supporting Information of this [article](http://pubs.acs.org/cgi-bin/abstract.cgi/jacsat/2000/122/i37/abs/ja991878x.html).
  2. Open the file with PyMOL and save it as a pdb file.
  3. Save the following pymolscript [Script.pml.tar](/images/2/26/Script.pml.tar "Script.pml.tar") and untar it to script.pml. Open it with an editor and adjust the Path_To_The_PDB. Open pymol and run the script with "@PATH_Of_The_Script/script.pml".
  4. You will get the following image, you can save it with "png filename.png" [![File.png](/images/b/b0/File.png)](/index.php/File:File.png)
  5. Use the script to modify your own pdb-file.



Retrieved from "[https://pymolwiki.org/index.php?title=Visualizing_a_computed_structure_-_a_commented_example&oldid=12178](https://pymolwiki.org/index.php?title=Visualizing_a_computed_structure_-_a_commented_example&oldid=12178)"


---

