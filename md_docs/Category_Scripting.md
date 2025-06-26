# Category: Scripting

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

## Auto arg

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

cmd.auto_arg controls auto-completion of command arguments in PyMOL's command line. It is an array of dictionaries. When Pymol's console attempts to auto-complete the _n_ -th argument of command _abc_ , it will look at `cmd.auto_arg[n]['abc']` (_n_ starts from 0). This is a list with three elements. For arguments which do not have an entry in cmd.auto_arg, the default is to auto-complete file names. 

  1. The first element is a lambda function which returns an Shortcut object. Shortcut object wraps the list of candidate strings, given to its constructor.
  2. The second element is a string, which describs the argument.
  3. The third element is a string, which will be added after autocompletion (postfix).



## Contents

  * 1 PYTHON EXAMPLE
  * 2 Pre-defined lambdas
  * 3 More Examples
  * 4 SEE ALSO



## PYTHON EXAMPLE
    
    
    cmd.auto_arg[0]['test']=[lambda: cmd.Shortcut(['abc','bcd']), '1st argument', ', ']
    

This code defines the auto-completion list for the first argument of command _test_. When you type 'test ' and press TAB, PyMOL will show you two candidates as: 
    
    
    parser: matching 1st argument:
         abc  bcd 
    

If you type 'a' and press TAB again, PyMOL will complete it to "test abc, ". Note that ', ' is added. 

## Pre-defined lambdas

In most cases, you just want to complete from object names or setting lists. Then you don't have to write your own lambda function. Just use **cmd.object_sc** for choosing from objects. 
    
    
    cmd.auto_arg[0]['test'] = [ cmd.object_sc, 'object', '']
    

Or even simpler, borrow from one of the core commands: 
    
    
    # "extract" also uses cmd.object_sc in its first argument
    cmd.auto_arg[0]['test'] = cmd.auto_arg[0]['extract']
    

Use **cmd.selection_sc** for setting names (like 'line_width'), **cmd.map_sc** for maps and **cmd.repres_sc** for representations ('sticks', 'lines', etc). 

## More Examples

The first argument of the [save](/index.php/Save "Save") command auto-completes to filenames by default. If you want for example filenames **and** object names to auto-complete, use this (from [pymol-users mailing list](http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10950.html)): 
    
    
    import glob
    names_filenames_sc = lambda: cmd.Shortcut(cmd.get_names() + glob.glob('*'))
    cmd.auto_arg[0]['save'] = [names_filenames_sc, 'filename or object name', '']
    

## SEE ALSO

  * [extend](/index.php/Extend "Extend")
  * parser.py, cmd.py, shortcut.py, completing.py



Retrieved from "[https://pymolwiki.org/index.php?title=Auto_arg&oldid=10872](https://pymolwiki.org/index.php?title=Auto_arg&oldid=10872)"


---

## Example Scripts

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Multiple Object Manipulation Scripts
    * 1.1 Dis/Enable Objects
    * 1.2 Comments in Scripts
      * 1.2.1 Sources



## Multiple Object Manipulation Scripts

### Dis/Enable Objects

If someone has many (possibly hundreds) of objects -- say distance objects -- one can turn them all on or off by using Python scripts within PyMol. The following script will extend the commands "sd" and "hd" for "show distance" and "hide distance," respectively. The first script, 
    
    
    from pymol import cmd 
    num_dist = 100 
            
    def show_dist(): 
        """ show all of my distance objects """ 
        for i in range(num_dist): 
            cmd.enable('_dist%s'%i) 
            
    def hide_dist(): 
        """ hide all of my distance objects """ 
        for i in range(num_dist): 
            cmd.disable('_dist%s'%i) 
            
    cmd.extend('sd',show_dist) 
    cmd.extend('hd',hide_dist)
    

works on 100 objects. We can extend the idea with more elegant scripting to work w/o forcing us to keep track of the number of objects: 
    
    
    def show_dist(): 
        dists = [name for name in cmd.get_names() if cmd.get_type(name) == 'object:distance'] 
        for name in dists: cmd.enable(name)
    

Now, just running "sd" would show all the distance objects. 

  


### Comments in Scripts

The hash-mark, "#" is the Python comment. However, when scripting in Python for PyMol note that the hash character (#) works to include comments on separate lines, but not at the end of a line that contains commands. So you can do 
    
    
    # Create separate dimer
    create dimer,(chain A,B)
    

but not 
    
    
    create dimer,(chain A,B)  # Create separate dimer
    

instead use: 
    
    
    create dimer,(chain A,B);  # Create separate dimer
    

#### Sources

Taken from the PyMol Users list. Python source by Michael Lerner. 

Retrieved from "[https://pymolwiki.org/index.php?title=Example_Scripts&oldid=13570](https://pymolwiki.org/index.php?title=Example_Scripts&oldid=13570)"


---

## Extend

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page is about an API function. For the selection operator with the same name, see [extend (selection operator)](/index.php/Extend_\(selection_operator\) "Extend \(selection operator\)"). 

Extend is an API-only function which binds a user defined function as a command into the PyMOL scripting language. 

## Contents

  * 1 Details
  * 2 PyMOL API
  * 3 Simple Example
  * 4 Advanced Example
  * 5 Notes and recommendations
  * 6 See Also



## Details

  * All command arguments are passed as strings to the Python function. This may require type conversion before those arguments can be used by the function, for example for numbers (int, float).
  * If the function has a **quiet** argument, then PyMOL will pass **quiet=0** to the command. Most PyMOL core commands have a default **quiet=1** argument and have no (or little) output when used in a Python script, but would print a more verbose feedback when invoked as a command.
  * If the function has a **_self** argument, then PyMOL will assign the global **cmd** module to **_self** , or if using the **pymol2** module, an instance of **pymol2.cmd2.Cmd**. This wrapper allows running multiple PyMOL instances in the same process.



## PyMOL API
    
    
    pymol.cmd.extend(string name, function function)
    

Simplified (since PyMOL 1.6), takes the function name as command name and can be used as a function [decorator](https://en.wikipedia.org/wiki/Python_syntax_and_semantics#Decorators): 
    
    
    pymol.cmd.extend(function function)
    

## Simple Example

Put the following in a Python script (e.g. file.py) and "run" the script (e.g. with "File > Run..."). 
    
    
    from pymol import cmd
    
    def foo(moo=2):
        print(moo)
    
    cmd.extend('foo', foo)
    

Or with decorator syntax (since PyMOL 1.6): 
    
    
    from pymol import cmd
    
    @cmd.extend
    def foo(moo=2):
        print(moo)
    

The following would now work within PyMOL's command line: 
    
    
    PyMOL>foo
    2
    PyMOL>foo 3
    3
    PyMOL>foo moo=5
    5
    PyMOL>foo ?
    Usage: foo [ moo ]
    

## Advanced Example

This example provides a command help, does proper type conversion of arguments and would work with **pymol2**. 
    
    
    from pymol import cmd
    
    @cmd.extend
    def adjust_vdw(selection="all", factor=0.5, quiet=1, _self=cmd):
        """
    DESCRIPTION
    
        The "adjust_vdw" command alters all vdw radii in the given selection
        by the given scale factor.
        """
        factor = float(factor)
        quiet = int(quiet)
    
        n = _self.alter(selection, "vdw = vdw * %f" % (factor))
    
        if n > 0:
          _self.rebuild()
    
        if not quiet:
            if n == 0:
                print(" No atoms in selection")
            else:
                print(" Updated VDW radii of %d atoms by factor %.2f" % (n, factor))
    

## Notes and recommendations

  * So far reasonable, provide default values for arguments (like **selection="all"**).
  * Provide a [docstring](https://en.wikipedia.org/wiki/Docstring#Python) in the same format like all other PyMOL core commands (with DESCRIPTION, USAGE, etc.), that will make calling "[help](/index.php?title=Help&action=edit&redlink=1 "Help \(page does not exist\)") yourcommand" a familiar experience for a user.
  * Take advantage of PyMOL's selection language, a command which takes a **selection** argument will feel more familiar to a user than a command which takes for example **object, chain** and **resid** arguments.
  * Pass boolean values as **0/1** instead of named values like **true, on, yes / false, ...** , that will be consistent with PyMOL's core commands.
  * You may want to set [auto_arg](/index.php/Auto_arg "Auto arg") to enable auto-completion for your command's arguments.
  * For security reasons, new PyMOL commands created using "extend" are not saved or restored in sessions (PSE files).



## See Also

  * [alias](/index.php/Alias "Alias")
  * [set_key](/index.php/Set_Key "Set Key")
  * [api](/index.php?title=Api&action=edit&redlink=1 "Api \(page does not exist\)")
  * [auto_arg](/index.php/Auto_arg "Auto arg")



Retrieved from "[https://pymolwiki.org/index.php?title=Extend&oldid=12776](https://pymolwiki.org/index.php?title=Extend&oldid=12776)"


---

## Find pairs

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**find_pairs** is an API only function which returns a list of atom pairs. Atoms are represented as (model,index) tuples. 

Can be restricted to hydrogen-bonding-like contacts. WARNING: Only checks atom orientation, not atom type (so would hydrogen bond between carbons for example), so make sure to provide appropriate atom selections. 

## Contents

  * 1 PyMOL API
  * 2 ARGUMENTS
  * 3 NOTE
  * 4 SEE ALSO



## PyMOL API
    
    
    list cmd.find_pairs(string selection1, string selection2, int state1=1, int state2=1,
        float cutoff=3.5, int mode=0, float angle=45)
    

## ARGUMENTS

  * selection1, selection2 = string: atom selections


  * state1, state2 = integer: state-index (only positive values allowed) {default: 1}


  * cutoff = float: distance cutoff {default: 3.5}


  * mode = integer: if mode=1, do coarse hydrogen bonding assessment {default: 0}


  * angle = float: hydrogen bond angle cutoff, only if mode=1 {default: 45.0}



## NOTE

Although this does a similar job like "[distance](/index.php/Distance "Distance")", it uses a completely different routine and the "mode" argument has different meanings! 

## SEE ALSO

  * [distance](/index.php/Distance "Distance")
  * Selection operators [within, around and expand](/index.php/Selection_Algebra "Selection Algebra")



Retrieved from "[https://pymolwiki.org/index.php?title=Find_pairs&oldid=9178](https://pymolwiki.org/index.php?title=Find_pairs&oldid=9178)"


---

## Get coords

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

get_coords is an API only function that returns the coordinates of a selection as a numpy array. 

_New in PyMOL 1.7.4_

The order of coordinates is that of the internal atom ordering, like in [iterate](/index.php/Iterate "Iterate"), and unlike in [get_coordset](/index.php/Get_coordset "Get coordset"). The function considers the object state and TTT matrices, like in [get_model](/index.php/Get_Model "Get Model"). 

## PyMOL API
    
    
    cmd.get_coords(str selection='all', int state=1)
    

**Returns:** ndarray, shape (N, 3), where N is the number of atoms if **state > 0**, or (natoms * nstates) if **state=0**. 

## Arguments

  * **selection** = str: atom selection {default: all}
  * **state** = int: state index or all states if state=0 {default: 1}



## See Also

  * [load_coords](/index.php/Load_coords "Load coords")
  * [get_coordset](/index.php/Get_coordset "Get coordset")
  * [Get Coordinates I](/index.php/Get_Coordinates_I "Get Coordinates I")



For state > 0, the function is equivalent to: 
    
    
    numpy.array(cmd.get_model(selection, state).get_coord_list())
    

Retrieved from "[https://pymolwiki.org/index.php?title=Get_coords&oldid=12260](https://pymolwiki.org/index.php?title=Get_coords&oldid=12260)"


---

## Get coordset

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

get_coordset is an API only function that returns the coordinates of one object-state (one "coordinate set") as a numpy array. 

_New in PyMOL 1.7.4_

The order of coordinates might differ from functions like [iterate](/index.php/Iterate "Iterate") or [get_coords](/index.php/Get_coords "Get coords"). The function does **not** take the state or TTT matrices into account. 

## Contents

  * 1 PyMOL API
  * 2 Arguments
  * 3 Examples
  * 4 See Also



## PyMOL API
    
    
    cmd.get_coordset(str name, int state=1, int copy=1)
    

**Returns:** ndarray, shape (N, 3), where N is the number of atoms. 

## Arguments

  * **name** = str: object name
  * **state** = int: state index {default: 1}
  * **copy** = 0/1: {default: 1} **Warning: only use copy=0 if you know what you're doing.** copy=0 will return a numpy array which is a wrapper of the internal coordinate set memory and provides read-write access. If the internal memory gets freed or reallocated, this wrapper will become invalid.



## Examples

Move an object to its geometric center: 
    
    
    cs = cmd.get_coordset('1ubq', copy=0)
    cs -= cs.mean(0)
    cmd.rebuild()
    

## See Also

  * [load_coordset](/index.php?title=Load_coordset&action=edit&redlink=1 "Load coordset \(page does not exist\)")
  * [get_coords](/index.php/Get_coords "Get coords")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_coordset&oldid=12633](https://pymolwiki.org/index.php?title=Get_coordset&oldid=12633)"


---

## Get object list

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**cmd.get_object_list** is an API only command which returns the list of all molecular objects in selection. It will not return selection names or objects which are not of type molecule. 

## PyMOL API
    
    
    list cmd.get_object_list(string selection='(all)')
    

_Important: Always put the selection in parenthesis to force PyMOL to interpret it as an atom selection. Otherwise this might fail on name wildcards, groups, etc._

For proper atom selections this should be identical to: 
    
    
    list cmd.get_names('objects', 0, string selection='(all)')
    

## Examples
    
    
    PyMOL> fetch 2x19 2xwu 1d7q, async=0
    PyMOL> print(cmd.get_object_list('solvent'))
    ['2x19', '2xwu']
    PyMOL> print(cmd.get_object_list('hydro'))
    ['1d7q']
    PyMOL> print(cmd.get_object_list('2x*'))    # FAIL!!! (fixed in PyMOL 1.7.6)
    ['2x19']
    PyMOL> print(cmd.get_object_list('(2x*)'))  # correct, with parenthesis
    ['2x19', '2xwu']
    

## See Also

  * [get_names](/index.php/Get_names "Get names")
  * [get_names_of_type](/index.php/Get_names_of_type "Get names of type")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_object_list&oldid=12545](https://pymolwiki.org/index.php?title=Get_object_list&oldid=12545)"


---

## Get raw alignment

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_raw_alignment** is an API only function that returns a list of lists of (object,index) tuples containing the raw per-atom alignment relationships. Alignment objects can be created by passing the "object" argument to [align](/index.php/Align "Align") or [super](/index.php/Super "Super"). 

_Please note:_

  * _The order of the atom tuples are not necessarily in the order in which the two (or more) selections were passed to[cmd.align](/index.php/Align "Align")._
  * _Will not return atom tuples of hidden objects (see also[hide_underscore_names](/index.php/Hide_underscore_names "Hide underscore names"))_
  * _Reimplemented in PyMOL 2.3, order of returned atom tuples can differ to previous versions_



## PYMOL API
    
    
    cmd.get_raw_alignment(string name)
    

## EXAMPLE
    
    
    # start a python block
    python
    
    # get two structures
    cmd.fetch('2xwu 2x19', async=0)
    
    # align and get raw alignment
    cmd.align('/2xwu//B//CA', '/2x19//B//CA', cycles=0, transform=0, object='aln')
    raw_aln = cmd.get_raw_alignment('aln')
    
    # print residue pairs (atom index)
    for idx1, idx2 in raw_aln:
        print('%s`%d -> %s`%d' % tuple(idx1 + idx2))
    
    #end the python block
    python end
    

To print residue numbers instead of atom indices: 
    
    
    # continued from previous example
    python
    
    idx2resi = {}
    cmd.iterate('aln', 'idx2resi[model, index] = resi', space={'idx2resi': idx2resi})
    
    # print residue pairs (residue number)
    for idx1, idx2 in raw_aln:
        print('%s -> %s' % (idx2resi[idx1], idx2resi[idx2]))
    
    python end
    

## SEE ALSO

  * [align](/index.php/Align "Align")
  * [find_pairs](/index.php/Find_pairs "Find pairs") returns a list of lists of (object,index) tuples as well
  * [set_raw_alignment](/index.php/Set_raw_alignment "Set raw alignment")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_raw_alignment&oldid=12848](https://pymolwiki.org/index.php?title=Get_raw_alignment&oldid=12848)"


---

## Load brick

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

load_brick is an API-only function to load an array of data as a [map](/index.php?title=Map&action=edit&redlink=1 "Map \(page does not exist\)") into PyMOL. 

[![](/images/1/16/Load_brick_example.png)](/index.php/File:Load_brick_example.png)

[](/index.php/File:Load_brick_example.png "Enlarge")

Volume visualization of a map (artificial data)

## Example
    
    
    from pymol import cmd
    from chempy.brick import Brick
    import numpy
    
    # grid spacing in Angstrom
    yourgridspacing = (0.2, 0.2, 0.2)
    
    # numpy array with some artificial data
    yourdata = numpy.zeros((100, 100, 100))
    yourdata[20:80, 20:80, 20:80] = 1.0
    yourdata[40:60, 40:60, 40:60] = 2.0
    
    # load into PyMOL
    brick = Brick.from_numpy(yourdata, yourgridspacing)
    cmd.load_brick(brick, 'yourmap')
    
    # visualize
    cmd.volume('yourvolume', 'yourmap', [\
          0.5, 'blue', 0.4, \
          0.6, 'blue', 0.0, \
          1.4, 'red', 0.0, \
          1.5, 'red', 0.9, \
        ])
    

## See Also

  * [volume](/index.php/Volume "Volume")
  * [isomesh](/index.php/Isomesh "Isomesh")



Retrieved from "[https://pymolwiki.org/index.php?title=Load_brick&oldid=12806](https://pymolwiki.org/index.php?title=Load_brick&oldid=12806)"


---

## Load coords

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

load_coords is an API only function to load selection coordinates. 

_CHANGED IN VERSION 1.7.3: This used to be the[load_coordset](/index.php?title=Load_coordset&action=edit&redlink=1 "Load coordset \(page does not exist\)") function. [load_coordset](/index.php?title=Load_coordset&action=edit&redlink=1 "Load coordset \(page does not exist\)") may load coordinates in different order (original order from PDB file) than load_coords (atom sorted order)._

## Contents

  * 1 PyMOL API
  * 2 Arguments
  * 3 Example
  * 4 See Also



## PyMOL API
    
    
    cmd.load_coords(sequence coords, str selection, int state=1)
    

## Arguments

  * **coords** = list: Nx3 float array
  * **selection** = str: atom selection
  * **state** = int: object state {default: 1}



## Example
    
    
    cmd.fragment('ala')
    coords = cmd.get_coords('ala')
    coords += 5.0
    cmd.load_coords(coords, 'ala')
    

## See Also

  * [get_coords](/index.php/Get_coords "Get coords")
  * [load_coordset](/index.php?title=Load_coordset&action=edit&redlink=1 "Load coordset \(page does not exist\)")
  * [get_coordset](/index.php/Get_coordset "Get coordset")



Retrieved from "[https://pymolwiki.org/index.php?title=Load_coords&oldid=11929](https://pymolwiki.org/index.php?title=Load_coords&oldid=11929)"


---

## Pymol.stored

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The **pymol.stored** helper variable serves as a namespace for user defined globals which are accessible from [iterate](/index.php/Iterate "Iterate")-like commands. The [iterate](/index.php/Iterate "Iterate") commands by default expose the **pymol** module namespace as the [globals](https://docs.python.org/2/library/functions.html#eval) dictionary, so **pymol.stored** is accessible as **stored** , and (user defined) members as **stored._membername_**. 

## Example

Count number of atoms with a PyMOL script: 
    
    
    stored.count = 0
    iterate all, stored.count += 1
    print("number of atoms:", stored.count)
    

Count number of atoms with a Python script: 
    
    
    from pymol import cmd
    from pymol import stored
    stored.count = 0
    cmd.iterate("all", "stored.count += 1")
    print("number of atoms:", stored.count)
    

## Problems

There is no guarantee that other scripts or plugins don't use the same **stored** member names. It is recommended that properly written plugins use the ["space" argument](/index.php/Iterate#"space"_argument "Iterate") with **cmd.iterate** and **cmd.alter** to define their own global namespace. 

## See Also

  * [iterate](/index.php/Iterate "Iterate")
  * [label](/index.php/Label "Label")



Retrieved from "[https://pymolwiki.org/index.php?title=Pymol.stored&oldid=12561](https://pymolwiki.org/index.php?title=Pymol.stored&oldid=12561)"


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

## Script Highlighting

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Syntax Highlighting for writing Pymol Scripts
    * 1.1 Sublime Text
    * 1.2 Atom
    * 1.3 Vim



## Syntax Highlighting for writing Pymol Scripts

When editing pymol scripts, it can be useful to highlight syntax features. This can make it more clear when both reading and writing scripts. Highlighting of python scripts is available by default on most text editors, but support for highlighting of pymol scripts requires third-party extensions.. There are pymol script highlighters available for a number of text editors. 

### Sublime Text

[![Highlighting in Sublime Text](/images/9/9a/Sublime_text_highlighting.png)](/index.php/File:Sublime_text_highlighting.png "Highlighting in Sublime Text")

Syntax highlighting for Sublime Text can be installed via [Package Control](https://sublime.wbond.net/packages/Pymol%20Language) and will update automatically. Alternatively, the most recent edits (pre-release) and the source file can be downloaded from [Github](https://github.com/bbarad/pymol_syntax). This package is compatible with both Sublime Text 2 and Sublime Text 3. 

The syntax file should work with other editors, such as TextMate, with minimal changes. 

### Atom

Sublime text syntax highlighting has been ported to the [Atom editor](https://atom.io) from the sublime text version and should behave identically. The package, named "language-pymol," can be found on the [Atom site](https://atom.io/packages/language-pymol) and installed via the atom package manager. 

### Vim

There is also a different syntax highlighting package available for vim. Behavior differs from the Sublime Text and Atom implementations. 

  * <https://github.com/speleo3/PyMol-syntax>
  * <http://www.vim.org/scripts/script.php?script_id=2814> (original version, outdated)



Retrieved from "[https://pymolwiki.org/index.php?title=Script_Highlighting&oldid=11931](https://pymolwiki.org/index.php?title=Script_Highlighting&oldid=11931)"


---

## Scripting FAQs

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Script within a Structure

Q: Is there a way to embed a script within a PDB structure? 

A: Yes, there are two ways 

  * Use the [Read_Pdbstr](/index.php/Read_Pdbstr "Read Pdbstr") command like,


    
    
    delete all
    cmd.read_pdbstr("""HETATM 1985  O00 MOH   132      18.797   6.477 -12.112 0.00  0.00           O\
    HETATM 1988  H03 MOH   132      18.437   7.229 -11.665  0.00  0.00    H\
    HETATM 1989  C04 MOH   132      17.737   5.662 -12.563  0.00  0.00    C\
    HETATM 1990  H05 MOH   132      18.129   4.785 -13.080  0.00  0.00    H\
    HETATM 1991  H06 MOH   132      17.096   6.211 -13.253  0.00  0.00    H\
    HETATM 1992  H07 MOH   132      17.130   5.322 -11.722  0.00  0.00    H""","mymolecule")
    show sticks
    

Don't forget the backslashes there -- it's one long command. [Example](http://chips.csb.ki.se/pymol/msg02104.html)

  * Or, as Warren noted, you can do with with a "p1m" file, which is like a "pml" file, but data can be included as well. Also note that p1m files are intended for web publishing, so embedded Python is disallowed. In p1m files, there is an "embed" command that enables this for PDB, MOL, MOL2, SDF, and XPLOR data. It works like this:


    
    
    embed tag, format
    ****REPLACE THIS LINE WITH YOUR DATA FILE****
    embed end
    
    
    
    load_embedded tag
    

where tag is some identifier of your choice. 

## Scripting and Command Line Options

PyMol will allow you to use command line options in a script. That means, the following 
    
    
    pymol -qrc script.py arg1 arg2
    

will work if you add a double-hyphen before the first argument in order to signal Pymol to stop interpreting arguments: 
    
    
    pymol -qrc script.py -- arg1 arg2
    

Then your script can get these arguments as follows: 
    
    
    from sys import argv
    my_argv = argv[1:]
    print my_argv[0], my_argv[1]
    

Retrieved from "[https://pymolwiki.org/index.php?title=Scripting_FAQs&oldid=11505](https://pymolwiki.org/index.php?title=Scripting_FAQs&oldid=11505)"


---

## Set raw alignment

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

set_raw_alignment is an API-only function to create an alignment object from lists of atom indices. 

_New in PyMOL 2.3_

## Arguments

  * **name** = str: alignment object name
  * **raw** = list: list of lists (columns) with `(model, index)` tuples
  * **guide** = str: name of guide object {default: first in list}
  * **state** = int: object state {default: 1}



## Example
    
    
    cmd.align('1t46', '1oky', object='aln')
    raw = cmd.get_raw_alignment('aln')
    cmd.delete('aln')
    cmd.set_raw_alignment('alnnew', raw)
    

## See Also

  * [get_raw_alignment](/index.php/Get_raw_alignment "Get raw alignment")
  * [align](/index.php/Align "Align")



Retrieved from "[https://pymolwiki.org/index.php?title=Set_raw_alignment&oldid=12847](https://pymolwiki.org/index.php?title=Set_raw_alignment&oldid=12847)"


---

## Set state order

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

set_state_order is an API only function to set the order of states for an object. 

_New in PyMOL 1.7.4_

## PyMOL API
    
    
    cmd.set_state_order(str name, list order)
    

## Arguments

  * **name** = str: object name
  * **order** = list of int: index array (1-based state indices)



## Example
    
    
    # reverse the order of a 20 model object
    cmd.set_state_order('1nmr', range(20, 0, -1))
    

Retrieved from "[https://pymolwiki.org/index.php?title=Set_state_order&oldid=11928](https://pymolwiki.org/index.php?title=Set_state_order&oldid=11928)"


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

## Category:Script Library

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

| PyMOL Script Library  |   
---|---|---  
[Running Scripts](/index.php/Running_Scripts "Running Scripts") | [Script Requests](/index.php/Category_talk:Script_Library#Requests "Category talk:Script Library") | [Policy](/index.php/Category_talk:Script_Library#Policy "Category talk:Script Library")  
[Structural Biology](/index.php/Category:Structural_Biology_Scripts "Category:Structural Biology Scripts") | [Objects and Selections](/index.php/Category:ObjSel_Scripts "Category:ObjSel Scripts") | [Math/Geometry/CGO](/index.php/Category:Math_Scripts "Category:Math Scripts")  
  
  * [Rotamer Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle")
  * [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure")
  * [WriteSS](/index.php/WriteSS "WriteSS")
  * [Iterate sses](/index.php/Iterate_sses "Iterate sses")
  * [Helicity check](/index.php/Helicity_check "Helicity check")
  * [Measure Distance](/index.php/Measure_Distance "Measure Distance")
  * [Translate And Measure](/index.php/Translate_And_Measure "Translate And Measure")
  * [Motif](/index.php/Motif "Motif")
  * [Nmr cnstr](/index.php/Nmr_cnstr "Nmr cnstr")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [InterfaceResidues](/index.php/InterfaceResidues "InterfaceResidues")
  * [ColorByRMSD](/index.php/ColorByRMSD "ColorByRMSD")
  * [Centroid](/index.php/Centroid "Centroid")
  * [Average b](/index.php/Average_b "Average b")
  * [Ss](/index.php/Ss "Ss")
  * [Ccp4 ncont](/index.php/Ccp4_ncont "Ccp4 ncont")
  * [AutoMultiFit](/index.php/AutoMultiFit "AutoMultiFit")
  * [Consistent View/ Map Inspect](/index.php/Consistent_View/_Map_Inspect "Consistent View/ Map Inspect")
  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [BbPlane](/index.php/BbPlane "BbPlane")
  * [AngleBetweenHelices](/index.php/AngleBetweenHelices "AngleBetweenHelices")
  * [AAindex](/index.php/AAindex "AAindex")
  * [Sidechaincenters](/index.php/Sidechaincenters "Sidechaincenters")
  * [Displacementmap](/index.php/Displacementmap "Displacementmap")
  * [HighlightAlignedSS](/index.php/HighlightAlignedSS "HighlightAlignedSS")
  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat")
  * [Color By Mutations](/index.php/Color_By_Mutations "Color By Mutations")
  * [Ccp4 contact](/index.php/Ccp4_contact "Ccp4 contact")
  * [Ccp4 pisa](/index.php/Ccp4_pisa "Ccp4 pisa")
  * [Propka](/index.php/Propka "Propka")
  * [DynoPlot](/index.php/DynoPlot "DynoPlot")
  * [Cart to frac](/index.php/Cart_to_frac "Cart to frac")
  * [Colorbydisplacement](/index.php/Colorbydisplacement "Colorbydisplacement")
  * [Resicolor plugin](/index.php/Resicolor_plugin "Resicolor plugin")
  * [Color by conservation](/index.php/Color_by_conservation "Color by conservation")
  * [Cyspka](/index.php/Cyspka "Cyspka")
  * [Elbow angle](/index.php/Elbow_angle "Elbow angle")
  * [Transformations](/index.php/Transformations "Transformations")
  * [DistancesRH](/index.php/DistancesRH "DistancesRH")
  * [Contact Surface](/index.php/Contact_Surface "Contact Surface")
  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")
  * [Bondpack](/index.php/Bondpack "Bondpack")
  * [Colorblindfriendly](/index.php/Colorblindfriendly "Colorblindfriendly")
  * [Load new B-factors](/index.php/Load_new_B-factors "Load new B-factors")
  * [RotationAxis](/index.php/RotationAxis "RotationAxis")
  * [Draw Protein Dimensions](/index.php/Draw_Protein_Dimensions "Draw Protein Dimensions")
  * [Dssr block](/index.php/Dssr_block "Dssr block")
  * [PDB plugin](/index.php/PDB_plugin "PDB plugin")
  * [RmsdByResidue](/index.php/RmsdByResidue "RmsdByResidue")
  * [Angle between domains](/index.php/Angle_between_domains "Angle between domains")
  * [Mcsalign](/index.php/Mcsalign "Mcsalign")
  * [Msms surface](/index.php/Msms_surface "Msms surface")
  * [Dssp](/index.php/Dssp "Dssp")
  * [LigAlign](/index.php/LigAlign "LigAlign")
  * [Xfit](/index.php/Xfit "Xfit")
  * [Intra xfit](/index.php/Intra_xfit "Intra xfit")
  * [Sst](/index.php/Sst "Sst")
  * [Load aln](/index.php/Load_aln "Load aln")
  * [Show contacts](/index.php/Show_contacts "Show contacts")
  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")

| 

  * [ConnectedCloud](/index.php/ConnectedCloud "ConnectedCloud")
  * [Zero residues](/index.php/Zero_residues "Zero residues")
  * [List Selection](/index.php/List_Selection "List Selection")
  * [List Colors](/index.php/List_Colors "List Colors")
  * [List Secondary Structures](/index.php/List_Secondary_Structures "List Secondary Structures")
  * [Selection Exists](/index.php/Selection_Exists "Selection Exists")
  * [Color Objects](/index.php/Color_Objects "Color Objects")
  * [Get Coordinates I](/index.php/Get_Coordinates_I "Get Coordinates I")
  * [Get Coordinates II](/index.php/Get_Coordinates_II "Get Coordinates II")
  * [Grepsel](/index.php/Grepsel "Grepsel")
  * [Findseq](/index.php/Findseq "Findseq")
  * [ToGroup](/index.php/ToGroup "ToGroup")
  * [Save sep](/index.php/Save_sep "Save sep")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [AlphaToAll](/index.php/AlphaToAll "AlphaToAll")
  * [Expand To Surface](/index.php/Expand_To_Surface "Expand To Surface")
  * [Removealt](/index.php/Removealt "Removealt")
  * [GetNamesInSel](/index.php/GetNamesInSel "GetNamesInSel")
  * [FindObjectsNearby](/index.php/FindObjectsNearby "FindObjectsNearby")
  * [CollapseSel](/index.php/CollapseSel "CollapseSel")
  * [SelInside](/index.php/SelInside "SelInside")
  * [Split selection](/index.php/Split_selection "Split selection")
  * [Split Object Along Axis](/index.php/Split_Object_Along_Axis "Split Object Along Axis")
  * [Split Movement](/index.php/Split_Movement "Split Movement")
  * [SaveGroup](/index.php/SaveGroup "SaveGroup")
  * [SelectClipped](/index.php/SelectClipped "SelectClipped")
  * [Rotkit](/index.php/Rotkit "Rotkit")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")
  * [Find buried waters](/index.php/Find_buried_waters "Find buried waters")
  * [Count molecules in selection](/index.php/Count_molecules_in_selection "Count molecules in selection")
  * [DistancesRH](/index.php/DistancesRH "DistancesRH")
  * [Cluster Count](/index.php/Cluster_Count "Cluster Count")
  * [ListSelection2](/index.php/ListSelection2 "ListSelection2")
  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")
  * [Flatten obj](/index.php/Flatten_obj "Flatten obj")
  * [Annotate v](/index.php/Annotate_v "Annotate v")
  * [Color cbcobj](/index.php/Color_cbcobj "Color cbcobj")
  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")

| 

  * [Perp maker](/index.php/Perp_maker "Perp maker")
  * [Symmetry Axis](/index.php/Symmetry_Axis "Symmetry Axis")
  * [CGO Text](/index.php/CGO_Text "CGO Text")
  * [Plane Wizard](/index.php/Plane_Wizard "Plane Wizard")
  * [Bounding Box](/index.php/Bounding_Box "Bounding Box")
  * [Ellipsoid](/index.php/Ellipsoid "Ellipsoid")
  * [TransformSelectionByCameraView](/index.php/TransformSelectionByCameraView "TransformSelectionByCameraView")
  * [Modevectors](/index.php/Modevectors "Modevectors")
  * [Center of mass](/index.php/Center_of_mass "Center of mass")
  * [DrawBoundingBox](/index.php/DrawBoundingBox "DrawBoundingBox")
  * [Slerpy](/index.php/Slerpy "Slerpy")
  * [CgoCircle](/index.php/CgoCircle "CgoCircle")
  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [Supercell](/index.php/Supercell "Supercell")
  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [BbPlane](/index.php/BbPlane "BbPlane")
  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat")
  * [Radius of gyration](/index.php/Radius_of_gyration "Radius of gyration")
  * [Cart to frac](/index.php/Cart_to_frac "Cart to frac")
  * [Mark center](/index.php/Mark_center "Mark center")
  * [Dump2CGO](/index.php/Dump2CGO "Dump2CGO")
  * [Axes](/index.php/Axes "Axes")
  * [Contact Surface](/index.php/Contact_Surface "Contact Surface")
  * [Cgo arrow](/index.php/Cgo_arrow "Cgo arrow")
  * [Cubes](/index.php/Cubes "Cubes")
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")
  * [Cgo grid](/index.php/Cgo_grid "Cgo grid")
  * [RotationAxis](/index.php/RotationAxis "RotationAxis")
  * [Draw Protein Dimensions](/index.php/Draw_Protein_Dimensions "Draw Protein Dimensions")

  
[User Interface](/index.php/Category:UI_Scripts "Category:UI Scripts") | [Biochemical Scripts](/index.php/Category:Biochemical_Scripts "Category:Biochemical Scripts") | [System Scripts](/index.php/Category:System_Scripts "Category:System Scripts")  
  
  * [ObjectByArrows](/index.php/ObjectByArrows "ObjectByArrows")
  * [PowerMate Dial OS X](/index.php/PowerMate_Dial_OS_X "PowerMate Dial OS X")
  * [Mouse modes](/index.php/Mouse_modes "Mouse modes")
  * [Key Wait](/index.php/Key_Wait "Key Wait")
  * [Grepset](/index.php/Grepset "Grepset")
  * [Apropos](/index.php/Apropos "Apropos")
  * [Stereo ray](/index.php/Stereo_ray "Stereo ray")
  * [PDB Web Services Script](/index.php/PDB_Web_Services_Script "PDB Web Services Script")
  * [Stereo Figures](/index.php/Stereo_Figures "Stereo Figures")
  * [Movit](/index.php/Movit "Movit")
  * [Spectrumbar](/index.php/Spectrumbar "Spectrumbar")
  * [Save2traj](/index.php/Save2traj "Save2traj")
  * [VisLoad](/index.php/VisLoad "VisLoad")
  * [Aa codes](/index.php/Aa_codes "Aa codes")
  * [Set toggle](/index.php/Set_toggle "Set toggle")
  * [Movie fade](/index.php/Movie_fade "Movie fade")
  * [Save settings](/index.php/Save_settings "Save settings")
  * [Dynamic mesh](/index.php/Dynamic_mesh "Dynamic mesh")
  * [ObjectFocus](/index.php/ObjectFocus "ObjectFocus")
  * [Make Figures](/index.php/Make_Figures "Make Figures")
  * [Inertia tensor](/index.php/Inertia_tensor "Inertia tensor")
  * [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz")
  * [RUCAP UM-5](/index.php/RUCAP_UM-5 "RUCAP UM-5")
  * [Format bonds](/index.php/Format_bonds "Format bonds")
  * [Get colors](/index.php/Get_colors "Get colors")
  * [Movie color fade](/index.php/Movie_color_fade "Movie color fade")
  * [Quickdisplays](/index.php/Quickdisplays "Quickdisplays")
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")
  * [Pytms](/index.php/Pytms "Pytms")
  * [Frame slider](/index.php/Frame_slider "Frame slider")

| 

  * [Show aromatics](/index.php/Show_aromatics "Show aromatics")
  * [Show hydrophobics](/index.php/Show_hydrophobics "Show hydrophobics")
  * [Show charged](/index.php/Show_charged "Show charged")
  * [Show hydrophilic](/index.php/Show_hydrophilic "Show hydrophilic")
  * [Pucker](/index.php/Pucker "Pucker")
  * [Resicolor](/index.php/Resicolor "Resicolor")
  * [CreateAtom](/index.php/CreateAtom "CreateAtom")
  * [Color h](/index.php/Color_h "Color h")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [Average b](/index.php/Average_b "Average b")
  * [AAindex](/index.php/AAindex "AAindex")
  * [Polarpairs](/index.php/Polarpairs "Polarpairs")
  * [Propka](/index.php/Propka "Propka")
  * [Forster distance calculator](/index.php/Forster_distance_calculator "Forster distance calculator")
  * [Cyspka](/index.php/Cyspka "Cyspka")
  * [DistancesRH](/index.php/DistancesRH "DistancesRH")
  * [ListSelection2](/index.php/ListSelection2 "ListSelection2")
  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")
  * [Quickdisplays](/index.php/Quickdisplays "Quickdisplays")
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")
  * [Load new B-factors](/index.php/Load_new_B-factors "Load new B-factors")
  * [Pytms](/index.php/Pytms "Pytms")
  * [Show contacts](/index.php/Show_contacts "Show contacts")
  * [Color cbcobj](/index.php/Color_cbcobj "Color cbcobj")
  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")

| 

  * [FilterByMol](/index.php/FilterByMol "FilterByMol")
  * [LoadDir](/index.php/LoadDir "LoadDir")
  * [Process All Files In Directory](/index.php/Process_All_Files_In_Directory "Process All Files In Directory")
  * [PythonTerminal](/index.php/PythonTerminal "PythonTerminal")
  * [Pdbsurvey](/index.php/Pdbsurvey "Pdbsurvey")
  * [Read PDB-String](/index.php/Read_PDB-String "Read PDB-String")
  * [Pml2py](/index.php/Pml2py "Pml2py")
  * [Tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")
  * [Monitor file continuously](/index.php/Monitor_file_continuously "Monitor file continuously")
  * [XML-RPC server](/index.php/XML-RPC_server "XML-RPC server")
  * [FetchLocal](/index.php/FetchLocal "FetchLocal")
  * [Find symbol](/index.php/Find_symbol "Find symbol")

  
[Third Party](/index.php/Category:ThirdParty_Scripts "Category:ThirdParty Scripts") |  |   
  
  * [Rasmolify](/index.php/Rasmolify "Rasmolify")
  * [Wfmesh](/index.php/Wfmesh "Wfmesh")
  * [Transform odb](/index.php/Transform_odb "Transform odb")
  * [PovRay](/index.php/PovRay "PovRay")
  * [Ccp4 ncont](/index.php/Ccp4_ncont "Ccp4 ncont")
  * [Ccp4 contact](/index.php/Ccp4_contact "Ccp4 contact")
  * [Ccp4 pisa](/index.php/Ccp4_pisa "Ccp4 pisa")
  * [Tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")
  * [Pymol2glmol](/index.php/Pymol2glmol "Pymol2glmol")
  * [Save Mopac](/index.php/Save_Mopac "Save Mopac")
  * [PoseView](/index.php/PoseView "PoseView")
  * [RUCAP UM-5](/index.php/RUCAP_UM-5 "RUCAP UM-5")

  
  
## Subcategories

This category has the following 7 subcategories, out of 7 total. 

### B

  * [Biochemical Scripts](/index.php/Category:Biochemical_Scripts "Category:Biochemical Scripts")



### M

  * [Math Scripts](/index.php/Category:Math_Scripts "Category:Math Scripts")



### O

  * [ObjSel Scripts](/index.php/Category:ObjSel_Scripts "Category:ObjSel Scripts")



### S

  * [Structural Biology Scripts](/index.php/Category:Structural_Biology_Scripts "Category:Structural Biology Scripts")
  * [System Scripts](/index.php/Category:System_Scripts "Category:System Scripts")



### T

  * [ThirdParty Scripts](/index.php/Category:ThirdParty_Scripts "Category:ThirdParty Scripts")



### U

  * [UI Scripts](/index.php/Category:UI_Scripts "Category:UI Scripts")



## Pages in category "Script Library"

The following 197 pages are in this category, out of 197 total. 

### A

  * [AAindex](/index.php/AAindex "AAindex")
  * [AlphaToAll](/index.php/AlphaToAll "AlphaToAll")
  * [Annotate v](/index.php/Annotate_v "Annotate v")
  * [Apropos](/index.php/Apropos "Apropos")
  * [AutoMultiFit](/index.php/AutoMultiFit "AutoMultiFit")
  * [Average b](/index.php/Average_b "Average b")
  * [Axes](/index.php/Axes "Axes")



### B

  * [B2transparency](/index.php/B2transparency "B2transparency")
  * [BbPlane](/index.php/BbPlane "BbPlane")
  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat")
  * [Bondpack](/index.php/Bondpack "Bondpack")
  * [Bounding Box](/index.php/Bounding_Box "Bounding Box")



### C

  * [CalcArea](/index.php/CalcArea "CalcArea")
  * [Cart to frac](/index.php/Cart_to_frac "Cart to frac")
  * [Ccp4 contact](/index.php/Ccp4_contact "Ccp4 contact")
  * [Ccp4 ncont](/index.php/Ccp4_ncont "Ccp4 ncont")
  * [Ccp4 pisa](/index.php/Ccp4_pisa "Ccp4 pisa")
  * [Cealign plugin](/index.php/Cealign_plugin "Cealign plugin")
  * [Center of mass](/index.php/Center_of_mass "Center of mass")
  * [Centroid](/index.php/Centroid "Centroid")
  * [Cgo arrow](/index.php/Cgo_arrow "Cgo arrow")
  * [Cgo grid](/index.php/Cgo_grid "Cgo grid")
  * [CGO Text](/index.php/CGO_Text "CGO Text")
  * [CgoCircle](/index.php/CgoCircle "CgoCircle")
  * [Check Key](/index.php/Check_Key "Check Key")
  * [Cluster Count](/index.php/Cluster_Count "Cluster Count")
  * [CMPyMOL](/index.php/CMPyMOL "CMPyMOL")
  * [CollapseSel](/index.php/CollapseSel "CollapseSel")
  * [Color by conservation](/index.php/Color_by_conservation "Color by conservation")
  * [Color By Mutations](/index.php/Color_By_Mutations "Color By Mutations")
  * [Color cbcobj](/index.php/Color_cbcobj "Color cbcobj")
  * [Color h](/index.php/Color_h "Color h")
  * [Color Objects](/index.php/Color_Objects "Color Objects")
  * [Colorblindfriendly](/index.php/Colorblindfriendly "Colorblindfriendly")
  * [Colorbydisplacement](/index.php/Colorbydisplacement "Colorbydisplacement")
  * [ColorByRMSD](/index.php/ColorByRMSD "ColorByRMSD")
  * [ConnectedCloud](/index.php/ConnectedCloud "ConnectedCloud")
  * [Consistent View/ Map Inspect](/index.php/Consistent_View/_Map_Inspect "Consistent View/ Map Inspect")
  * [Contact Surface](/index.php/Contact_Surface "Contact Surface")
  * [Count molecules in selection](/index.php/Count_molecules_in_selection "Count molecules in selection")
  * [CreateAtom](/index.php/CreateAtom "CreateAtom")
  * [CreateSecondaryStructure](/index.php/CreateSecondaryStructure "CreateSecondaryStructure")
  * [Cubes](/index.php/Cubes "Cubes")
  * [Cyspka](/index.php/Cyspka "Cyspka")



### D

  * [Displacementmap](/index.php/Displacementmap "Displacementmap")
  * [DistancesRH](/index.php/DistancesRH "DistancesRH")
  * [Distancetoatom](/index.php/Distancetoatom "Distancetoatom")
  * [Draw Protein Dimensions](/index.php/Draw_Protein_Dimensions "Draw Protein Dimensions")
  * [DrawGridBox](/index.php/DrawGridBox "DrawGridBox")
  * [DrawBoundingBox](/index.php/DrawBoundingBox "DrawBoundingBox")
  * [Dssp](/index.php/Dssp "Dssp")
  * [Dump2CGO](/index.php/Dump2CGO "Dump2CGO")
  * [Dynamic mesh](/index.php/Dynamic_mesh "Dynamic mesh")
  * [DynoPlot](/index.php/DynoPlot "DynoPlot")



### E

  * [Elbow angle](/index.php/Elbow_angle "Elbow angle")
  * [Ellipsoid](/index.php/Ellipsoid "Ellipsoid")
  * [Expand To Surface](/index.php/Expand_To_Surface "Expand To Surface")
  * [Extra fit](/index.php/Extra_fit "Extra fit")



### F

  * [FetchLocal](/index.php/FetchLocal "FetchLocal")
  * [FilterByMol](/index.php/FilterByMol "FilterByMol")
  * [Find buried waters](/index.php/Find_buried_waters "Find buried waters")
  * [Find symbol](/index.php/Find_symbol "Find symbol")
  * [FindObjectsNearby](/index.php/FindObjectsNearby "FindObjectsNearby")
  * [Findseq](/index.php/Findseq "Findseq")
  * [FindSurfaceCharge](/index.php/FindSurfaceCharge "FindSurfaceCharge")
  * [FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues")
  * [Flatten obj](/index.php/Flatten_obj "Flatten obj")
  * [FocalBlur](/index.php/FocalBlur "FocalBlur")
  * [Format bonds](/index.php/Format_bonds "Format bonds")
  * [Forster distance calculator](/index.php/Forster_distance_calculator "Forster distance calculator")
  * [Frame slider](/index.php/Frame_slider "Frame slider")



### G

  * [Get colors](/index.php/Get_colors "Get colors")
  * [Get Coordinates I](/index.php/Get_Coordinates_I "Get Coordinates I")
  * [Get Coordinates II](/index.php/Get_Coordinates_II "Get Coordinates II")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")
  * [GetNamesInSel](/index.php/GetNamesInSel "GetNamesInSel")
  * [Grepsel](/index.php/Grepsel "Grepsel")
  * [Grepset](/index.php/Grepset "Grepset")



### H

  * [Hbplus](/index.php/Hbplus "Hbplus")
  * [Colorama](/index.php/Colorama "Colorama")
  * [Helicity check](/index.php/Helicity_check "Helicity check")
  * [HighlightAlignedSS](/index.php/HighlightAlignedSS "HighlightAlignedSS")



### I

  * [ImmersiveViz](/index.php/ImmersiveViz "ImmersiveViz")
  * [Inertia tensor](/index.php/Inertia_tensor "Inertia tensor")
  * [InterfaceResidues](/index.php/InterfaceResidues "InterfaceResidues")
  * [Intra xfit](/index.php/Intra_xfit "Intra xfit")
  * [Iterate sses](/index.php/Iterate_sses "Iterate sses")



### J

  * [Jump](/index.php/Jump "Jump")



### K

  * [Kabsch](/index.php/Kabsch "Kabsch")
  * [Key Wait](/index.php/Key_Wait "Key Wait")



### L

  * [LatticeGenerator](/index.php/LatticeGenerator "LatticeGenerator")
  * [Launching From a Script](/index.php/Launching_From_a_Script "Launching From a Script")
  * [LigAlign](/index.php/LigAlign "LigAlign")
  * [List Colors](/index.php/List_Colors "List Colors")
  * [List Secondary Structures](/index.php/List_Secondary_Structures "List Secondary Structures")
  * [List Selection](/index.php/List_Selection "List Selection")
  * [ListSelection2](/index.php/ListSelection2 "ListSelection2")
  * [Load aln](/index.php/Load_aln "Load aln")
  * [Load img stack](/index.php/Load_img_stack "Load img stack")
  * [Load new B-factors](/index.php/Load_new_B-factors "Load new B-factors")
  * [LoadDir](/index.php/LoadDir "LoadDir")



### M

  * [Make Figures](/index.php/Make_Figures "Make Figures")
  * [MakeVinaCommand](/index.php/MakeVinaCommand "MakeVinaCommand")
  * [Mark center](/index.php/Mark_center "Mark center")
  * [Mcsalign](/index.php/Mcsalign "Mcsalign")
  * [Measure Distance](/index.php/Measure_Distance "Measure Distance")
  * [Modevectors](/index.php/Modevectors "Modevectors")
  * [Monitor file continuously](/index.php/Monitor_file_continuously "Monitor file continuously")
  * [Motif](/index.php/Motif "Motif")
  * [Mouse modes](/index.php/Mouse_modes "Mouse modes")
  * [Movie color fade](/index.php/Movie_color_fade "Movie color fade")
  * [Movie fade](/index.php/Movie_fade "Movie fade")
  * [Movit](/index.php/Movit "Movit")
  * [Msms surface](/index.php/Msms_surface "Msms surface")



### O

  * [ObjectByArrows](/index.php/ObjectByArrows "ObjectByArrows")
  * [ObjectFocus](/index.php/ObjectFocus "ObjectFocus")



### P

  * [Pairwise distances](/index.php/Pairwise_distances "Pairwise distances")
  * [PDB Web Services Script](/index.php/PDB_Web_Services_Script "PDB Web Services Script")
  * [Pdbsurvey](/index.php/Pdbsurvey "Pdbsurvey")
  * [Perp maker](/index.php/Perp_maker "Perp maker")
  * [Plane Wizard](/index.php/Plane_Wizard "Plane Wizard")
  * [Plot noe](/index.php/Plot_noe "Plot noe")
  * [Pml2py](/index.php/Pml2py "Pml2py")
  * [PoseView](/index.php/PoseView "PoseView")
  * [PowerMate Dial OS X](/index.php/PowerMate_Dial_OS_X "PowerMate Dial OS X")
  * [Process All Files In Directory](/index.php/Process_All_Files_In_Directory "Process All Files In Directory")
  * [Propka](/index.php/Propka "Propka")
  * [Psico](/index.php/Psico "Psico")
  * [Pucker](/index.php/Pucker "Pucker")
  * [Pymol2glmol](/index.php/Pymol2glmol "Pymol2glmol")
  * [Pytms](/index.php/Pytms "Pytms")



### Q

  * [Quick dist](/index.php/Quick_dist "Quick dist")
  * [Quickdisplays](/index.php/Quickdisplays "Quickdisplays")



### R

  * [Radius of gyration](/index.php/Radius_of_gyration "Radius of gyration")
  * [Rasmolify](/index.php/Rasmolify "Rasmolify")
  * [Read PDB-String](/index.php/Read_PDB-String "Read PDB-String")
  * [Removealt](/index.php/Removealt "Removealt")
  * [Renumber](/index.php/Renumber "Renumber")
  * [Resicolor](/index.php/Resicolor "Resicolor")
  * [Resicolor plugin](/index.php/Resicolor_plugin "Resicolor plugin")
  * [RmsdByResidue](/index.php/RmsdByResidue "RmsdByResidue")
  * [Rotamer Toggle](/index.php/Rotamer_Toggle "Rotamer Toggle")
  * [RotationAxis](/index.php/RotationAxis "RotationAxis")
  * [Rotkit](/index.php/Rotkit "Rotkit")
  * [RUCAP UM-5](/index.php/RUCAP_UM-5 "RUCAP UM-5")



### S

  * [Save Mopac](/index.php/Save_Mopac "Save Mopac")
  * [Save pdb with anisou](/index.php/Save_pdb_with_anisou "Save pdb with anisou")
  * [Save sep](/index.php/Save_sep "Save sep")
  * [Save settings](/index.php/Save_settings "Save settings")
  * [Save2traj](/index.php/Save2traj "Save2traj")
  * [SaveGroup](/index.php/SaveGroup "SaveGroup")
  * [Ss](/index.php/Ss "Ss")
  * [Select sites](/index.php/Select_sites "Select sites")
  * [SelectClipped](/index.php/SelectClipped "SelectClipped")
  * [Selection Exists](/index.php/Selection_Exists "Selection Exists")
  * [SelInside](/index.php/SelInside "SelInside")
  * [Set toggle](/index.php/Set_toggle "Set toggle")
  * [Show aromatics](/index.php/Show_aromatics "Show aromatics")
  * [Show bumps](/index.php/Show_bumps "Show bumps")
  * [Show charged](/index.php/Show_charged "Show charged")
  * [Show contacts](/index.php/Show_contacts "Show contacts")
  * [Show hydrophilic](/index.php/Show_hydrophilic "Show hydrophilic")
  * [Show hydrophobics](/index.php/Show_hydrophobics "Show hydrophobics")
  * [Nmr cnstr](/index.php/Nmr_cnstr "Nmr cnstr")
  * [ShowLigandWaters](/index.php/ShowLigandWaters "ShowLigandWaters")
  * [Sidechaincenters](/index.php/Sidechaincenters "Sidechaincenters")
  * [Slerpy](/index.php/Slerpy "Slerpy")
  * [Spectrum states](/index.php/Spectrum_states "Spectrum states")
  * [Spectrumany](/index.php/Spectrumany "Spectrumany")
  * [Spectrumbar](/index.php/Spectrumbar "Spectrumbar")
  * [Split chains](/index.php/Split_chains "Split chains")
  * [Split Movement](/index.php/Split_Movement "Split Movement")
  * [Split Object Along Axis](/index.php/Split_Object_Along_Axis "Split Object Along Axis")
  * [Split object](/index.php/Split_object "Split object")
  * [Split selection](/index.php/Split_selection "Split selection")
  * [Sst](/index.php/Sst "Sst")
  * [Stereo ray](/index.php/Stereo_ray "Stereo ray")
  * [Stereo Figures](/index.php/Stereo_Figures "Stereo Figures")
  * [Supercell](/index.php/Supercell "Supercell")
  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [Symmetry Axis](/index.php/Symmetry_Axis "Symmetry Axis")



### T

  * [Theseus](/index.php/Theseus "Theseus")
  * [Tiff2ccp4](/index.php/Tiff2ccp4 "Tiff2ccp4")
  * [TMalign](/index.php/TMalign "TMalign")
  * [ToGroup](/index.php/ToGroup "ToGroup")
  * [Transform odb](/index.php/Transform_odb "Transform odb")
  * [Transformations](/index.php/Transformations "Transformations")
  * [TransformSelectionByCameraView](/index.php/TransformSelectionByCameraView "TransformSelectionByCameraView")
  * [Translate And Measure](/index.php/Translate_And_Measure "Translate And Measure")



### U

  * [Uniprot features](/index.php/Uniprot_features "Uniprot features")



### V

  * [VisLoad](/index.php/VisLoad "VisLoad")



### W

  * [Wfmesh](/index.php/Wfmesh "Wfmesh")
  * [WriteSS](/index.php/WriteSS "WriteSS")



### X

  * [Xfit](/index.php/Xfit "Xfit")
  * [XML-RPC server](/index.php/XML-RPC_server "XML-RPC server")



### Z

  * [Zero residues](/index.php/Zero_residues "Zero residues")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Script_Library&oldid=10250](https://pymolwiki.org/index.php?title=Category:Script_Library&oldid=10250)"


---

