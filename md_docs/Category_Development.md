# Category: Development

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

## Ideas

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

> * * *
> 
> * * *
> 
> **DISCLAIMER:** This page was originally developed for the 2008 Google Summer of Code competition and is kept only for historical reasons. If you want to discuss an idea or feature implementation just send and e-mail to the [PyMOL_mailing_list](/index.php/PyMOL_mailing_list "PyMOL mailing list"). If you want to request a new feature or report a bug, please fill a ticket at [sourceforge](http://sourceforge.net/projects/pymol/%7C)
> 
> * * *
> 
> * * *

## Contents

  * 1 Ideas for PyMOL Development
  * 2 Where to Start
  * 3 Choosing a Topic
  * 4 Integration Ideas (Linking Out to Useful Open-Source Tools)
  * 5 High-Level Enhancement Ideas (Mostly Python-oriented)
  * 6 Low-Level Enhancement Ideas (Mostly C-oriented)
  * 7 Difficult C-level Code Refactoring Ideas
  * 8 Ideas Involving Proprietary APIs
  * 9 Ideas for Plugin Developers
  * 10 More Ideas (Please add your own!)



### Ideas for PyMOL Development

## Where to Start

Always start with Python and only delve down into the C code when absolutely necessary. Although PyMOL is mostly a C-based application, much of the that code is opaque, fragile, and unforgiving. Although C code refactoring is an important project goal, such work may not be ideal since once mistake could potentially to destabilize the entire platform. 

Fortunately, the Python interpreter and the PyMOL command and selection languages make it possible to extend PyMOL safely and quickly. Even when performance is critical, Python should be the interface between external C, C++, and Java code and PyMOL's internal C data structures. 

## Choosing a Topic

The best open-source code is usually written by an end-users attempting to meet their own pressing needs. So if you have already have a specific need which relates to PyMOL, then we strongly encourage you to follow up on that first! 

If you are looking for ideas, then try to seek out enhancements and/or integrations that will impact the largest potential user base. For example, imagine what new things might be useful to virtually all medicinal chemists, all structural biologists, all movie-makers, all paper-writers, and so forth. 

The ideas below are organized by category. Right now, integration with other open-source projects seems like the approach most likely to yield significant benefit, so those ideas are first. 

## Integration Ideas (Linking Out to Useful Open-Source Tools)

In most cases, depending on the need, integration can be accomplished through standalone Python scripts, through new PyMOL commands, through PyMOL Wizards, or via Tkinter plugins. 

  * APBS (electrostatics calculations): Improve the existing plugin. Michael Lerner is currently leading this effort. See [APBS](/index.php/APBS "APBS")



    

    yea ([Tree](/index.php/User:Inchoate "User:Inchoate") [Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)") [Vcpmartin](/index.php?title=User:Vcpmartin&action=edit&redlink=1 "User:Vcpmartin \(page does not exist\)") [Siderator](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)") [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)
    Feel free to mention specific plugin enhancements you'd like! [Michael Lerner](/index.php/User:Mglerner "User:Mglerner")

  * RDKit (cheminformatics, depiction, UFF cleanup, etc.): Lots of potential here, however C++ coding may be necessary for more advanced integration tasks. [RDKit home](http://www.rdkit.org)



    

    yea ([Markvanraaij](/index.php?title=User:Markvanraaij&action=edit&redlink=1 "User:Markvanraaij \(page does not exist\)")) / nay (0)

  * GIMP (image manipulation): Streamline & document the process of exporting images from PyMOL into GIMP and preparing them for submission to scientific Journals.



    

    yea (0) / nay (0)

  * Jmol (publishing visualizations inside of web pages): Liason between PyMOL & Jmol projects to develop a shared molecular visualization data model compatible with both applications (meaning, being able to import/export sessions from one to the other).



    

    yea ([tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)"))

  * Firefox (plugin): Develop an PyMOL plugin compatible with Firefox.



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)") [assaff](/index.php?title=User:Assaff&action=edit&redlink=1 "User:Assaff \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay ([tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)"))

  * MMTK (molecular mechanics -- Python/flexible): Develop the ability to round-trip molecular systems from PyMOL, into MMTK, and back.



    

    yea (0) / nay (0)

  * GROMACS (molecular mechanics -- C/fast) - Maybe some ideas can be shared with this guy. [Gromacs GUI](http://www.kde-apps.org/content/show.php/Gromacs+GUI+?content=47665)



    

    yea ([Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)") [Michael Lerner](/index.php/User:Mglerner "User:Mglerner")) / nay (0)

  * OpenOffice (escape Microsoft hegemony): Develop an PyMOL plugin.



    

    yea (0) / nay (0)

  * IPython integration (interactive shell): a robust alternative to the PyMOL command line?



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)") [assaff](/index.php?title=User:Assaff&action=edit&redlink=1 "User:Assaff \(page does not exist\)")) / nay (0)

  * R (statistics): PyMOL a 3D viewer environment for visualizating & manipulating large statistical data sets?



    

    yea (0) / nay (0)

Are there other key open-source packages we might specifically target for integration with PyMOL, either through GSoC or beyond? 

## High-Level Enhancement Ideas (Mostly Python-oriented)

  * Work on  MolViz



    

    yea (0) / nay (0)

  * Develop new plugins which automate routine tasks.



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)")) / nay (0)

  * Improve the Python API documentation.



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay (0)

  * Flesh out the new "from pymol2 import PyMOL" instance-based PyMOL API.



    

    yea (0) / nay (0)

  * Develop alternate Tkinter "skins" (for custom OEM-like applications).



    

    yea (0) / nay (0)

  * Develop a Tkinter/TOGL widget which holds a PyMOL viewer instance.



    

    yea (0) / nay (0)

  * Develop a PyQt widget which holds a PyMOL viewer instance.



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)")) / nay (0)

  * Create a plugin-manager GUI in the style of Firefox, Rythmbox, Gedit, Eclipse. A GUI where it is easy to turn off/on plugins, configure them and see help-contents for them. Maybe also some way to paste a url to install a new Plugin.



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)

  * Add a plugin for a GUI window with the same functionality as the "Control Panel" window in SwissPDB Viewer.



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)")) / nay (0)

  * Extend and modify the PyMOL command language so as to be compatible with existing RasMol and/or Jmol scripts.



    

    yea ([tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw"))

  * Enhance the Mutagenesis Wizard in order to support Nucleic acids and/or Sugars.



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay (0)

  * Better tab completion for commands



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)") [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) /nay(0)

## Low-Level Enhancement Ideas (Mostly C-oriented)

  * Enable editing of displayed sequence alignments.



    

    yea ([Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)"), [Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [slaw](/index.php/User:Slaw "User:Slaw"), [Speleo3](/index.php/User:Speleo3 "User:Speleo3")) / nay (0) 

    Would this then feed back to the structural alignment? [Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)")
    A search function in the sequence would also be nice. [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)")

  * Add multi-line textual annotations



    

    yea (0) / nay (0)

  * Support additional annotation object including: arrow, lines, and blobs.



    

    yea ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)")) / nay (0)

  * Add display of secondary structure into the sequence viewer.



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)

  * Enable per-object Z clipping (especially in the ray tracer)



    

    yea ([Gregori](/index.php?title=User:Gregori&action=edit&redlink=1 "User:Gregori \(page does not exist\)"), [Xevi](/index.php/User:Xevi "User:Xevi"), [Johnm](/index.php?title=User:Johnm&action=edit&redlink=1 "User:Johnm \(page does not exist\)"), [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0) 

    I would go a step further with fully customizable selection-based clipping planes (XYZ, color and transparency) ([Gregori](/index.php?title=User:Gregori&action=edit&redlink=1 "User:Gregori \(page does not exist\)"))

  * Highlight H-bonds, salt bridges, Pi-stacking, Pi-cations, etc.



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)"), [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")) / nay (0)

  * Build in a simple forcefield and energy minimizer suitable for use with Mutagenesis. [optimize](/index.php/Optimize "Optimize")



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")) / nay (0)

  * Incorporate a suite of standard NMR visualizations (restraint violations, per-residue RMS, etc.)



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)"), [Vcpmartin](/index.php?title=User:Vcpmartin&action=edit&redlink=1 "User:Vcpmartin \(page does not exist\)"), [rpetrenko](/index.php?title=User:Rpetrenko&action=edit&redlink=1 "User:Rpetrenko \(page does not exist\)"), [SteffenG](/index.php?title=User:SteffenG&action=edit&redlink=1 "User:SteffenG \(page does not exist\)")) / nay (0)

  * Enumeration and display of low-energy conformers. see [optimize](/index.php/Optimize "Optimize")



    

    yea ([Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)")) / nay (0) 

    This could be done by integrating RDKit, I think. [Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)")

  * Automated structure grafting (poor-man's homology modeling).



    

    yea ([tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay ([Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)"), [Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)")]) 

    Perhaps a plugin to Modeller instead? ([Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"))

  * Import of alignment files.



    

    yea ([Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")) / nay (0)

  * Implement IMD (Interactive Molecular Dynamics) Interface, see <http://www.ks.uiuc.edu/Research/vmd/imd/>



    

    yea ( [slaw](/index.php/User:Slaw "User:Slaw")) / nay (0)

  * Add buttons for **Set ChainID** and **Renumber Residues From...** to Edit menu or Actions (wrapper around Alter command)



    

    yea ([Sheehanj](/index.php?title=User:Sheehanj&action=edit&redlink=1 "User:Sheehanj \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)

## Difficult C-level Code Refactoring Ideas

  * Assemble a test suite which thoroughly exercises the existing code (a prerequisite to major refactoring).



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)"),[Karo](/index.php/User:Karo "User:Karo")) / nay (0)

  * Catch & handle memory-allocation failures gracefully (instead of crashing).



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)")) / nay (0)

  * Replace PyMOL's memory management & custom containers with a simple runtime object model.



    

    yea (0) / nay (0)

  * Separate the View and the Controllers from the Model so that they can all run asynchronously (on multiple cores).



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)")) / nay (0)

  * Enable generalized undo of changes made to the Model.



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)")"[Vcpmartin](/index.php?title=User:Vcpmartin&action=edit&redlink=1 "User:Vcpmartin \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")", [Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)")) / nay (0)

  * Clean up the internal matrix handling code.



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)"),[Karo](/index.php/User:Karo "User:Karo")) / nay (0)

## Ideas Involving Proprietary APIs

Since these involve closed-source APIs and infrastructure, they aren't suitable for open-source development efforts. However, such requests are noted here for the sake of complete coverage. 

  * Create a Windows port with "native" look & feel. <\- Could this be done in PyQT or PyGTK?. Then it would look "native", but be cross-platform and not proprietary.



    

    yea (0) / nay ( [slaw](/index.php/User:Slaw "User:Slaw"))

  * Integrate directly via Mathematica via MathLink.



    

    yea (0) / nay ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"))

  * Further enhance JyMOL (Java-JNI/wrapped PyMOL)



    

    yea ([Tree](/index.php/User:Inchoate "User:Inchoate"), [Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay (0)

  * Integrate with Matlab.



    

    yea (0) / nay ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"))

  * Quicklook Plugin on the Mac



    

    yea ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)")) / nay (0)

## Ideas for Plugin Developers

The range of things that PyMOL plugins can do has grown by leaps and bounds over the last several years. It may be time for some new infrastructure to facilitate plugin development. 

  * Plugins that span multiple files. This means 1. the ability to have them 2. the ability for users to easily install them.



    

    yea (1) / nay (0) ([Michael Lerner](/index.php/User:Mglerner "User:Mglerner"))

  * Config files that persist between sessions. This would be useful for things like remembering the locations of external programs, etc., without requiring the user to modify their .pymolrc.



    

    yea (1) / nay (0) ([Michael Lerner](/index.php/User:Mglerner "User:Mglerner"))

  * Restart PyMOL. After installing a plugin, PyMOL needs to be restarted to use the plugin. This is tedious when testing plugins under development.



    

    yea (1) /nay (0) ([Hongbo Zhu](/index.php/User:Hongbo_zhu "User:Hongbo zhu"))

## More Ideas (Please add your own!)

  * I often need to "replicate" previously PyMOL-made figures using newer coordinates/structures. It would be a big help (to me), if PyMOL could produce an exported text file from the old PyMOL .pse file of graphical settings that I could then modify and apply as a pymol script to the new structure -- I'm thinking of something like an old-style MOLSCRIPT.INP file...



    

    yea ([tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)

  * [MolViz](http://molviz.cs.toronto.edu/molviz) is a project to incorporate head tracking input into [PyMol](http://pymol.sourceforge.net/). This is accomplished through a [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz") script written in Python to control the molecule's position using the existing [PyMol API](http://www.pymolwiki.org/index.php/Category:Commands). Related projects would include: 
    * Improving the existing [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz") PyMol plugin for more precise control of the environment.
    * Developing new input drivers for the Wiimote form of control. This would require some bluetooth hacking.
    * Implementing some other forms of input for head tracking, such as fisheye head tracking, IR webcam tracking, etc (refer to the end of this [[video](http://www.youtube.com/watch?v=ncShaY4VSac)] for a better description).



    

    yea (0) / nay (0)

  * Provide a 2D chemical depiction of the current 3D view.



    

    yea ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)")) / nay ([tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) 

    RDKit?

  * Spreadsheet view with additional information (e.g. IC50's).



    

    yea ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)")) / nay ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)"))

  * Create additional documentation, screen casts, & tutorials.



    

    yea ([Markvanraaij](/index.php?title=User:Markvanraaij&action=edit&redlink=1 "User:Markvanraaij \(page does not exist\)") [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay (0)

  * Export 3D PDF images.



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw")) / nay ([Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)"))

  * Add extra "Single Word Selectors" like "nucleic", "protein", "water", "ions", "backbone" (for nucleic acids or proteins), "mainchain", "sidechain"



    

    yea ( [slaw](/index.php/User:Slaw "User:Slaw")) / nay (0)

  * Add functionality that allows you to select atoms based on their location (i.e. select (x_coordinate < 10) and (z_coordinate > 0))



    

    yea (0) / nay (0)

  * Implemented: ~~set pdb_mirror option to use PDB mirrors other than RCSB for fetching structures (PDBe, PDBj); the EBI mirror is much faster from Europe for example.~~



    

    yea ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)
    There is already a [fetch_host](/index.php/Fetch_Host "Fetch Host") setting --[Speleo3](/index.php/User:Speleo3 "User:Speleo3") 17:05, 29 October 2010 (UTC)

  * have the ability to link the TK console to the viewer so that users don't have to constantly alt+tab between what they want



    

    yea (0) / nay (0)

  * have the ability to disable typing in the viewer and automatically type in the TK console (I like being able to cut/paste/home/end)



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")) / nay (0)

  * I see a lot of "can pymol do this" threads - any ideas of a good UI for a page of "things PyMOL can do?"



    

    yea (0) / nay (0)

  * make an option where I can turn on a coordinate grid - perhaps an object that is a cuboid grid around any object in the view so I can still alter how it's rendered?



    

    yea (0) / nay (0)

  * iPhone / Nexus One app(s)



    

    yea (0) / nay ( [slaw](/index.php/User:Slaw "User:Slaw"))

  * single-color bonds between nonidentical or any spherical atoms colored specifically



    

    yea (0) / nay (0)

  * double bonds as two parallel cylinders



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)")) / nay (0)

  * export scenes as [WebGL](http://www.khronos.org/webgl/) / Could make mobile apps, presentation plugins obsolete



    

    yea ([Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [Cowsandmilk](/index.php?title=User:Cowsandmilk&action=edit&redlink=1 "User:Cowsandmilk \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)

  * export images in vector format



    

    yea ([Vvostri](/index.php?title=User:Vvostri&action=edit&redlink=1 "User:Vvostri \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)") [slaw](/index.php/User:Slaw "User:Slaw") ) / nay (0)

  * Keynote plugin



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")) / nay (0)

  * Electron density contour sliders (plugin: [isoslider](/index.php/Isoslider "Isoslider"))



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)"), [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)")) / nay (0) 

    May be even map the contouring level onto the wheel [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)")

  * Implemented: ~~Internal support for FFT, i.e. reading of map coefficients~~



    

    yea ([Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")) / nay (0)

  * Implemented: ~~Automatic electron density map generation from PDB mmcif files or user-supplied MTZ files~~



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [Aschreyer](/index.php?title=User:Aschreyer&action=edit&redlink=1 "User:Aschreyer \(page does not exist\)"), [Johnm](/index.php?title=User:Johnm&action=edit&redlink=1 "User:Johnm \(page does not exist\)"),[Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)")) / nay (0)

  * Display of crystallographic symmetry and NCS axes, with possibility of showing symbols indicating what kind of axes they are



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)"), [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)")) / nay (0) 

    Support for labeling cell axes and origin, etc... [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)"), [Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)")

  * Automatic symmetry expansion to show overall crystal packing; automatic generation of scenes showing crystal packing interface details



    

    yea ([Luca Jovine](/index.php?title=User:Lucajovine&action=edit&redlink=1 "User:Lucajovine \(page does not exist\)"), [tstout](/index.php?title=User:Tstout&action=edit&redlink=1 "User:Tstout \(page does not exist\)"), [Carsten Schubert](/index.php?title=User:Siderator&action=edit&redlink=1 "User:Siderator \(page does not exist\)")) / nay (0)

  * Update selections automatically: add a new option to selection action drop-down menu called "update". It will re-evaluate the selection expression and update the content of the selection automatically.



Retrieved from "[https://pymolwiki.org/index.php?title=Ideas&oldid=12156](https://pymolwiki.org/index.php?title=Ideas&oldid=12156)"


---

## PluginArchitecture

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes PyMOL's plugin architechture since version 1.5.0.6. 

## Contents

  * 1 API Entry Points
  * 2 __init_plugin__
    * 2.1 Example for a PyQt5 GUI (PyMOL 2.x)
    * 2.2 Example for a legacy Tkinter GUI (PyMOL 1.x)
  * 3 More than one file per plugin
  * 4 Config files
  * 5 PyMOL OS Fellowship Project 2011/2012
  * 6 See Also



## API Entry Points

A plugin is a Python module which uses PyMOL's API. The following entrypoints add functionality to PyMOL: 

  1. [`pymol.cmd.extend(func)`](/index.php/Extend "Extend") registers a function as a PyMOL command
  2. `pymol.plugins.addmenuitemqt(label, callback)` adds a plugin menu item, should be called inside `__init_plugin__`
  3. `MyPlugin.__init_plugin__(app)` is called during plugin initialization (if defined by the plugin)



## __init_plugin__

A plugin can implement an `__init_plugin__(app)` function (_previously`__init__`, deprecated_) which is called during plugin initialization. It is passed an object which is compatible with the legacy `PMGApp` instance and provides access to the Tkinter parent (`app.root`). 

### Example for a PyQt5 GUI (PyMOL 2.x)

_See also the[Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")_
    
    
    def __init_plugin__(app=None):
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('My Qt Plugin', myqtdialog)
    
    def myqtdialog():
        from pymol.Qt import QtWidgets
        QtWidgets.QMessageBox.information(None, 'My Plugin Title', 'Hello World')
    

### Example for a legacy Tkinter GUI (PyMOL 1.x)

It is important that all Tkinter objects are created with the `app.root` parent, otherwise legacy plugins won't work in PyMOL 2.0. 
    
    
    def __init_plugin__(app):
        app.menuBar.addmenuitem('Plugin', 'command',
            label='My Tk Plugin',
            command=lambda: mytkdialog(app.root))
    
    def mytkdialog(parent):
        try:
            import tkMessageBox  # Python 2
        except ImportError:
            import tkinter.messagebox as tkMessageBox  # Python 3
        tkMessageBox.showinfo(parent=parent, title='My Plugin Title', message='Hello World')
    

## More than one file per plugin

Plugins can be either a single Python file, or a directory with an `__init__.py` file. 

**Single file layout:**
    
    
    MyPlugin.py (defines "__init_plugin__(app=None)" function)
    

**Directory layout:**
    
    
    MyPlugin/
    ├── data
    │   ├── datasheet.txt
    │   └── image.png
    ├── submodule1.py
    ├── submodule2.py
    └── __init__.py (defines "__init_plugin__(app=None)" function)
    

They can be zipped (`.zip`) for distribution, the [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager") handles the unzipping during installation. 
    
    
    zip -r MyPlugin-1.0.zip MyPlugin
    

## Config files

PyMOL offers a basic API to store custom settings in `~/.pymolpluginsrc.py`. Only basic types (`str`, `int`, `float`, `list`, `tuple`, `dict`, etc.) are suppored (those who produce executable code with `repr()`). 

**Store:**
    
    
    pymol.plugins.pref_set(key, value)
    
    # needed if pymol.plugins.pref_get("instantsave") == False
    pymol.plugins.pref_save()
    

**Load:**
    
    
    value = pymol.plugins.pref_get(key, default=None)
    

**Example:**
    
    
    apbsbinary = pymol.plugins.pref_get("APBS_BINARY_LOCATION", "/usr/bin/apbs")
    

## PyMOL OS Fellowship Project 2011/2012

[Thomas Holder](/index.php/User:Speleo3 "User:Speleo3") was working on the new plugin system as part of his [2011-2012 Fellowship](http://pymol.org/fellowship). The improvements were incorporated into PyMOL 1.5.0.5. 

  * graphical [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * multiple files per plugin/module (see #More than one file per plugin)
  * user plugin directory
  * rarely used plugins can be disabled for better performance and enabled on demand
  * metadata support (version, author, description, citation, tags, ...)
  * install from online repositories
  * settings API for plugins (see #Config files)



## See Also

  * [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")
  * [Plugins Tutorial](/index.php/Plugins_Tutorial "Plugins Tutorial")
  * [Plugins](/index.php/Plugins "Plugins")



Retrieved from "[https://pymolwiki.org/index.php?title=PluginArchitecture&oldid=12829](https://pymolwiki.org/index.php?title=PluginArchitecture&oldid=12829)"


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

