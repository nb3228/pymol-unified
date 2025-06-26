# Category: System Scripts

## FetchLocal

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Note
  * 3 Syntax
  * 4 The Code
  * 5 See Also



## Overview

Try fetching local copy of PDB files before going on internet to fetch a remote copy. 

Sometimes you have a local copy of the PDB structure files, or [biological unit](http://www.rcsb.org/pdb/static.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html) files. They are, however, usually organized into sub-directories (for a good reason). Setting [Fetch_Path](/index.php/Fetch_Path "Fetch Path") to the top directory will not make PyMOL to use local files before go on internet to fetch a remote copy. This script extends the functionality of the PyMOL command [Fetch](/index.php/Fetch "Fetch") by searching local sub-directories for target files first. 

The search order of the script: 1. local copy of PDB; 2. [Fetch_Path](/index.php/Fetch_Path "Fetch Path"); 3. remote repository ([Fetch_Host](/index.php/Fetch_Host "Fetch Host")). 

## Note

The script only copes with `.pdb` and `.pdb1` (or `.pdb2`, `.pdb3`, ...) file types. 

## Syntax

See [Fetch](/index.php/Fetch "Fetch")

## The Code
    
    
    """ 2012_09_18, Hongbo Zhu <hongbo.zhu.cn gmail>
        DESCRIPTION: Look for pdb* on local disk before fetching remotely.
    """
    import os
    import cmd
    
    #####################
    # user configuration
    localpdbdir = '/your/dir/to/PDB/divided'
    localbudir  = '/your/dir/to/PDB_BU/divided' # set to '' (empty string) if N/A
    
    def fetchlocal(code, name='', state=0, finish=1, discrete=-1,
                   multiplex=-2, zoom=-1, type='pdb', async=-1, path=None,
                   file=None, quiet=1):  
        """ Default parameters are the same as in function definition for fetch() 
            in file ./modules/pymol/importing.py
        """
        if type.startswith('pdb'): # pdb files, or biological unit files .pdb1 .pdb2 etc.
            if type == 'pdb': # pdb files
                localdir = localpdbdir
            else:
                localdir = localbudir
            for c in string.split(str(code)):
                subdir = c.lower()[1:3]
                if type == 'pdb': 
                    if os.path.isfile('%s/%s/pdb%s.ent.gz' % (localdir,subdir,c.lower())):
                        localfn = '%s/%s/pdb%s.ent.gz' % (localdir,subdir,c.lower())
                    elif os.path.isfile('%s/%s/pdb%s.ent' % (localdir,subdir,c.lower())):
                        localfn = '%s/%s/pdb%s.ent.gz' % (localdir,subdir,c.lower())
                    else:
                        localfn = None
                else:             
                    if os.path.isfile('%s/%s/%s.%s.gz' % (localdir,subdir,c.lower(),type)):
                        localfn =  '%s/%s/%s.%s.gz' % (localdir,subdir,c.lower(),type)
                    elif os.path.isfile('%s/%s/%s.%s' % (localdir,subdir,c.lower(),type)):
                        localfn = '%s/%s/%s.%s' % (localdir,subdir,c.lower(),type)
                    else:
                        localfn = None
    
                if localfn is not None:
                    print 'Load local file instead of fetching from internet: ', localfn 
                    cmd.load(filename=localfn, object=name, format='pdb', state=state, 
                             finish=finish, discrete=discrete, multiplex=multiplex, 
                             zoom=zoom, quiet=quiet)
                else:  # otherwise hand it over to pymol fetch function
                    cmd.fetch(c,name,state,finish,discrete,multiplex,zoom,type,async,path,file,quiet)
        else:
            cmd.fetch(code,name,state,finish,discrete,multiplex,zoom,type,async,path,file,quiet)
    
    
    cmd.extend('fetchlocal', fetchlocal)
    

# See Also

[Fetch](/index.php/Fetch "Fetch"), [Fetch_Path](/index.php/Fetch_Path "Fetch Path"), [Fetch_Host](/index.php/Fetch_Host "Fetch Host"), [Psico](/index.php/Psico "Psico")

Retrieved from "[https://pymolwiki.org/index.php?title=FetchLocal&oldid=12456](https://pymolwiki.org/index.php?title=FetchLocal&oldid=12456)"


---

## FilterByMol

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This script filters through all the PDBs in the parent dir (you can easily the the directory it scans). For each molecule, it saves **just** the ligands/heteroatoms (excluding the waters). This gives you a simple way to filter through a database of proteins looking only at their ligands. 

This script, as noted below, works on the objects at the level of a **molecule**. While we can [iterate](/index.php/Iterate "Iterate") over atom number (ID), residue number (resi), etc we do not have any such "MOLID". So, we provide this simple workaround. You might need this file because if you have a residue (like #111 from 3BEP) that consists of a molecule and an atom then there's no other way to save the separate pieces (of molecule/atom) into two (or more files). As you can see in the following listing, if we iterate over the hetero atoms (and not waters) in 3BEP we get, 
    
    
    PyMOL>iterate bymol het, print resi, resn, ID, chain, segi, alt
    111 5CY 6473 C  
    111 5CY 6474 C  
    111 5CY 6476 C  
    111 5CY 6477 C  
    111 5CY 6478 C  
    111 5CY 6479 C  
    111 5CY 6480 C  
    111 5CY 6481 C  
    111 5CY 6482 C  
    111 5CY 6483 C  
    111 5CY 6484 C  
    111 5CY 6485 C  
    111 5CY 6486 C  
    111 5CY 6487 C  
    111 5CY 6488 C  
    111 5CY 6489 C  
    111 5CY 6490 C
    

which does not allow us to separate the two pieces. 

## The Code
    
    
    python
    
    #
    # This simple script will filter through all PDBs in a directory, and for each one
    # save all the ligands/heterotoms (that aren't waters) to their own file.  This
    # script operates at the level of molecules, not residues, atoms, etc.  Thus, if
    # you have a ligand that PyMOL is treating as ONE residue, but is actually two
    # separate molecules, or a molecule and an atom, then you will get multiple files.
    #
    
    from glob import glob
    from os import path
    from pymol import stored
    
    theFiles = glob("../*.pdb");
    
    for f in theFiles:
        # load the file
        cmd.load(f);
        # remove the protein and waters
        cmd.remove("polymer or resn HOH");
    
        cmd.select("input", "all")
        cmd.select("processed", "none")
        mol_cnt = 0
    
        while cmd.count_atoms("input"):
            # filter through the selections, updating the lists
            cmd.select("current","bymolecule first input")
            cmd.select("processed","processed or current")
            cmd.select("input","input and not current")
    
            # prepare the output parameters
            curOut = path.basename(f).split(".")[0] + "_" + str(mol_cnt).zfill(5) + "_het.pdb"
            curSel = "current"
            
            # save the file
            cmd.save( curOut, curSel );
            print "Saved " + curSel + " to " + curOut
            
            mol_cnt = mol_cnt + 1;
    
        # remove all to move to next molecule
        cmd.delete("*");        
    
    python end
    

Retrieved from "[https://pymolwiki.org/index.php?title=FilterByMol&oldid=7665](https://pymolwiki.org/index.php?title=FilterByMol&oldid=7665)"


---

## Find symbol

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script looks for the given symbol in a given module(s). This can be handy if you forgot which module something is in or just want to query for a substring in a set of modules. 

# Usage
    
    
    # search pymol and pymol.cmd for any symbol starting with 'mov'
    
    fs("mov")
    ['cmd.get_movie_length',
     'cmd.get_movie_locked',
     'cmd.get_movie_playing',
     'cmd.mmove',
     'cmd.move',
     'cmd.movie',
     'cmd.moving',
     'cmd.remove',
     'cmd.remove_picked']
    
    # Search PyMOL's CMD module for something called align
    
    fs("align", "cmd")
    ['cmd.align', 'cmd.alignto', 'cmd.cealign', 'cmd.get_raw_alignment']
    

# The Code
    
    
    import pymol
    import inspect
    import pprint
    
    def fs(needle,haystack=["pymol","cmd"]):
        """
        This script will find the 'needle' in the 'haystack,' where the former is
        a symbol to be found in the latter, which is a module.
        """
        
        if type("") == type(haystack):
    
            haystack = [haystack,]
    
        for mod in haystack:
    
            found_list = map(lambda x: "%s.%s" % (mod,x), [name for name,obj in inspect.getmembers(eval(mod)) if needle in name])
    
        pprint.pprint(found_list)
    
        return found_list
    

Retrieved from "[https://pymolwiki.org/index.php?title=Find_symbol&oldid=11434](https://pymolwiki.org/index.php?title=Find_symbol&oldid=11434)"


---

## LoadDir

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Load all files of the suffix **suff** from the directory **dirName** , where **suff** and **dirName** are function parameters. 

_This script is superseded by the[loadall](/index.php/Loadall "Loadall") command, which was added in PyMOL 1.7.2._

## Install

  1. copy the source below to the a file called "loadDir.pml" somewhere on your computer
  2. load the file with "run /your/path/toLoadDir/loadDir.pml"
  3. run loadDir. See examples below.



## Examples
    
    
    # load the script
    run ~/loadDir.pml
    
    # load all SD files from /tmp
    loadDir /tmp, sdf
    loadDir /tmp, .sdf
    loadDir /tmp, *.sdf
    # even stupid stuff works; hopefully as one would want.
    # load all PDBs from /tmp
    loadDir /tmp, foo.woo.pdb
    
    # load all the PDBs in all the directories under ./binders/ERK
    loadDir ./binders/ERK/*, .pdb
    
    # load the PDBs into groups: now we can load all the files in the tree under
    # ./ERK into the group "ERK" and the files from ./SYK into the group "SYK"
    loadDir ./binders/ERK/*, .pdb, group=ERKb
    loadDir ./binders/SYK/*, .pdb, group=SYKb
    

## The Code
    
    
    from glob import glob
    from os.path import sep, basename
    
    def loadDir(dirName=".", suff="pdb", group=None):
            """
            Loads all files with the suffix suff (the input parameter) from the directory dirName).
    
            dirName:        directory path
            suff:           file suffix.  Should be simply "pdb" or "sdf" or similar.  Will accept the
                            wildcard and dot in case the user doesn't read this.  So, "*.pdb", ".pdb",
                            and "pdb" should work.  The suffix can be anything valid that PyMOL knows
                            how to natively load.
            group:          groupName to add the files to.
    
            example:
                    # load all the PDBs in the current directory
                    loadDir
    
                    # load all SD files from /tmp
                    loadDir /tmp, "sdf"
    
            notes:
                    make sure you call this script w/o quotes around your parameters:
                            loadDir ., .pdb
                    as opposed to
                            loadDir ".", "*.pdb"
                    Use the former.
            """
    
            g = dirName + sep + "*." + suff.split(".")[-1]
    
            for c in glob( g ):
                    cmd.load(c)
    
                    if ( group != None ):
                            cmd.group( group, basename(c).split(".")[0], "add" )
    
    cmd.extend("loadDir", loadDir)
    

Retrieved from "[https://pymolwiki.org/index.php?title=LoadDir&oldid=12540](https://pymolwiki.org/index.php?title=LoadDir&oldid=12540)"


---

## Monitor file continuously

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script can be used to continuously check the modification timestamp on a file (any format, although this example assumes it's a PDB file), and re-loads it whenever the timestamp changes. As written it is intended to be started from the command line, but this is not a requirement. 

# The Code
    
    
    from pymol import cmd
    import threading
    import time
    import os
    import sys
    
    class pymol_file_monitor (object) :
      def __init__ (self,
          file_name, 
          time_wait=1) : # time in seconds between mtime check
        self.file_name = file_name
        self.time_wait = time_wait
        self.watch = True # this can be toggled elsewhere to stop updating
        self.mtime = 0
        t = threading.Thread(target=self.check_file)
        t.setDaemon(1)
        t.start()
        print "Watching file %s" % file_name
    
      def check_file (self) :
        while (self.watch) :
          if (os.path.exists(self.file_name)) :
            print "checking..."
            mtime = os.path.getmtime(self.file_name)
            if (mtime > self.mtime) :
              self.mtime = mtime
              print "Re-loading %s" % self.file_name
              cmd.load(self.file_name, state=1)
          time.sleep(self.time_wait)
    
    if (__name__ == "pymol") :
      monitor = pymol_file_monitor("status.pdb")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Monitor_file_continuously&oldid=10469](https://pymolwiki.org/index.php?title=Monitor_file_continuously&oldid=10469)"


---

## Pdbsurvey

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Hi everyone ! 

I have created a script because I was tired of browsing the pdb for endless searches of structures relevant to my field. This scripts takes as input a text file in which you copied your favourite keywords, and the number of days you want to search back (by default, it will look at the structures added the last 50 days). It generates a report text file that contains the pdb id and name of the relevant structures that have been recently added. All you need to do is add 'run pdbsurvey.py' to your startup file, and create a text file called 'keywords.txt' with your keywords separated by an end-of-line character. Then you're ready to go. Just hit 'pdbsurvey' from within your PyMol instance, and the program returns the report file. The pdb ftp server is updated weekly. 
    
    
    from pymol import cmd
    
    def pdbsurvey(days=50):
    
        """USAGE : pdbsurvey (<days>)
        Surveys the updates added to the PDB (ftp.rcsb.org) in the last
        50 days (or otherwise specified when calling this function) for
        entries that contain the words specified in the file
        keywords.txt.
        """
        print days
    
        import ftplib
        import time
        import os
    
    
    
        def todaymerge():
            """Puts today's date in a pdb format string.
            """
            date=time.localtime()
            fyear="%i" %(date[0])
            fmonth="%i" %(date[1])
            if date[1]<10:
                fmonth="0"+"%i" %(date[1])
            fday="%i" %(date[2])
            if date[2]<10:
                fday="0"+"%i" %(date[2])
            dateS=fyear+fmonth+fday
            return dateS
    
        def file2list(filename):
            """Low-level routine to brainlessly implement
            file.read().
            """
            fq=open(filename,'rb')
            linesS=fq.read()
            fq.close()
            LIST=linesS.splitlines()
            return LIST
    
        def connect2pdb():
            """Opens an anonymous socket to ftp://ftp.rcsb.org
            """
            f=ftplib.FTP()
            f.connect ('ftp.rcsb.org')
            f.login ()
            print "Remote connection established","\n"
            return f
    
        def decrementdate(dateS):
            """given a string date (pdb format yyyymmdd)
            this routine returns a string of the day before
            (sadly assuming that every month has 31 days, but
            no big deal here...).
            """
            #decompose dateS into components
            yearS=dateS[0]+dateS[1]+dateS[2]+dateS[3]
            monthS=dateS[4]+dateS[5]
            dayS=dateS[6]+dateS[7]
    
            #convert each into integers
            yearI=int(yearS)
            monthI=int(monthS)
            dayI=int(dayS)
    
            #actual stuff
            dayI=dayI-1
            if dayI==0:
                dayI=31
                monthI-=1
                if monthI==0:
                    monthI=12
                    yearI-=1
            dayS="%i" %(dayI)
            monthS="%i" %(monthI)
            yearS="%i" %(yearI)
            if dayI<10:
                dayS="0"+dayS
            if monthI<10:
                monthS="0"+monthS
            #and finally...
            dateS=yearS+monthS+dayS
            return dateS
    
        def findlastdir(dateS,f,days):
            """Puts the names of the "recent" directories in the
            list named "directoriesL".
            """
            directoriesL=[]
            for p in range(days):
                dateS=decrementdate(dateS)
                attempt="/pub/pdb/data/status/"+dateS
                try:
                    f.cwd(attempt)
                    directoriesL.append(attempt)
                except:
                    pass
            return directoriesL
    
        def compilinfile(directoriesL,f):
            """lists all structures in the added.pdb files
            contained in the directories specified in directoriesL
            """
            command="RETR added.pdb"
            handle=open("donotedit.dat","wrb")
            for k in directoriesL:
                f.cwd(k)
                print "Currently in directory ",f.pwd()
                f.retrbinary(command,handle.write)
            handle.close()
            return len(directoriesL)
    
        def listparser():
            """Extracts the pdbids from donotedit.dat file,
            and stacks them into the list pdbidsL
            """
            linesL=file2list("donotedit.dat")
            pdbidsL=[]
            for iter in linesL:
                pdbidsL.append(iter[57:61])
            pdbidsL.sort()
            return pdbidsL
    
        def currentrelease(f):
            """Stores the content of cmpd_res.idx file
            This file contains the equivalencies pdbid<->title
            for all current entries of the PDB.
            """
            command="RETR cmpd_res.idx"
            f.cwd("/pub/pdb/derived_data/index/")
            print "Currently in directory ",f.pwd()
            fq=open("dictionnary.dat",'wrb')
            f.retrbinary(command,fq.write)
            fq.close()
            dictL=file2list("dictionnary.dat")
            return dictL
    
        def extract(pdbidsL,dictL):
            """Populates dictionnaryD with pdb entries found in the
            latest releases.
            """
            dictionnaryD={}
            problemL=[]
            extractL=[dictionnaryD,problemL]
            for i in dictL:
                tempS=i[0:4].lower()
                for ii in pdbidsL:
                    if ii == tempS:
                        title=i[14:216]
                        extractL[0][ii]=title
            if len(extractL[0]) != len(pdbidsL):
                print "Dimension mismatch, seeking troublemaker..."
                for i in pdbidsL:
                    if not any(i == ii for ii in extractL[0]):
                        extractL[1].append(i)
            return extractL
    
        def disconnectpdb(f):
            """Diconnects the current ftp session
            """
            f.quit()
            print "Remote connection terminated","\n"
            return f
    
        def releventries(dictionnaryD):
            """Generates a cleaned dictionary with only entries
            that have one or more keywords specified in the local
            user-defined keywords.txt file
            """
            keywL=file2list("keywords.txt")
            relevdicD={}
            for i in keywL:
                for elem,temp in dictionnaryD.items():
                    if i in temp:
                        relevdicD[elem]=temp
            return relevdicD
    
        def diskcleanup(filelist=["donotedit.dat","dictionnary.dat"]):
            """Lo-level disk cleanup to free up memory without the user
            """
            for filename in filelist:
                command='DEL '+filename
                os.system(command)
            return "clean"
    
    
    
    
        print "Welcome in the auto-PDB updater !"
    
        print "Survey of updates made since",days,"days ago."
        
        print "Acquisition of local time..."
        dateS=todaymerge()                                                 #Initializes dateS
        print "today is ",dateS
        print "Connecting to remote ftp server..."
        f=connect2pdb()                                                    #Connect anonymously to ftp.rcsb.org
    
        print "Acquisition of latest added remote directories..."
        directoriesL=findlastdir(dateS,f,days)                             #Lists recent directories in directoriesL
        if len(directoriesL)==0:
            print "No updates have been found since",days,"ago. Starting over with 50 days ago."
            directoriesL=findlastdir(dateS,f,50)
    
        print "Acquisition of latest addedremote files..."
        updatesnumberI=compilinfile(directoriesL,f)                        #Concatenates the corresponding added.pdb into donotedit.dat
    
        print "Parsing of latest entries..."
        pdbidsL=listparser()                                               #Recent names now present in the pdbidsL list (one name per element)
    
        print "Acquisition of the current pdb distribution..."
        dictL=currentrelease(f)                                            #Populates dictL with the current entries of the PDB
    
        print "Parsing of the current pdb distribution into [code,title] tuples..."
        extractL=extract(pdbidsL,dictL)                                    #generates the dictionnary of latest releases key:PDBid ; definition:pdbtitle
    
        print "Disconnection from the remote ftp server..."
        f=disconnectpdb(f)                                                 #Closes the ftp instance
    
        print "Extraction of the relevant entries..."
        relevdicD=releventries(extractL[0])                               #Generates a subset of dictionnary D with criterion being "has keywords contained in keywords.txt in its title"
    
        print "Cleaning program-generated temporary files..."
        clean=diskcleanup()                                                #Cleans the mess generated by the program
    
        reportL=[]
        reportL.append("\n")
        reportL.append("###############REPORT########################################\n")
        reportL.append("\n")
        lendictS="%i" %(len(dictL))
        chmilblik = 'The current pdb version (as of '+dateS+") has "+lendictS+" entries.\n"
        reportL.append(chmilblik)
        line="The most recent directory is : "+directoriesL[0]+".\n"
        reportL.append(line)
        updatesnumberS="%i" %(updatesnumberI)
        entriesnumber="%i" %(len(extractL[0]))
        line="The "+updatesnumberS+" last updates ("+entriesnumber+" entries) have been examined.\n"
        reportL.append(line)
        diclengthS="%i" %(len(relevdicD))
        line=diclengthS+" are relevant to you :\n"
        reportL.append(line)
        for i,x in relevdicD.items():
            entry=i+" : "+x+"\n"
            reportL.append(entry)
        problemS=""
        for i in extractL[1]:
            problemS=i+";"+problemS
        problemS="["+problemS
        problemS=problemS.strip(";")
        problemS=problemS+"]"
        lineS="The entries "+problemS+" raised problems,"
        reportL.append(lineS)
        reportL.append("they should be examined manually.")
        reportL.append("\n")
        reportL.append("###############END OF REPORT#################################\n")
        report=open("report.aut","w")
        for elem in reportL:
            print elem
            report.writelines(elem+'\n')
        report.close()
        command2='start keywords.txt'
        command3='start report.aut'
        os.system(command2)
        os.system(command3)
    
    cmd.extend("pdbsurvey",pdbsurvey)
    

Thank you for any feedback, ways of improving it,... 

Retrieved from "[https://pymolwiki.org/index.php?title=Pdbsurvey&oldid=8010](https://pymolwiki.org/index.php?title=Pdbsurvey&oldid=8010)"


---

## Pml2py

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script converts a pml script to a python script. 

See [pymol-users mailing list](http://sourceforge.net/mailarchive/message.php?msg_id=27331786) (Subject: Convert pml script to Pymol Python script, Fri, 8 Apr 2011). 
    
    
    from __future__ import print_function
    
    import sys
    from pymol import cmd, parsing
    
    def pml2py(filename, out=sys.stdout):
        '''
    DESCRIPTION
    
        Convert a pml script to python syntax.
    
    USAGE
    
        pml2py infile [, outfile]
        '''
        def quote(args):
            args = iter(args)
            for arg in args:
                if '=' not in arg:
                    prefix = ''
                else:
                    prefix, arg = arg.split('=', 1)
                    prefix += '='
                    arg = arg.lstrip()
                yield prefix + repr(arg)
    
        class stackiter:
            def __init__(self, collection):
                self.iterator = iter(collection)
                self.stack = []
            def __iter__(self):
                return self
            def next(self):
                if len(self.stack):
                    return self.stack.pop()
                return next(self.iterator)
            __next__ = next
            def push(self, v):
                self.stack.append(v)
    
        basestring = (bytes, unicode) if (bytes is str) else (bytes, str)
        if isinstance(out, basestring):
            out = open(out, 'w')
    
        print('''
    # automatically converted from "%s"
    import pymol
    from pymol import *
    ''' % (filename), file=out)
    
        handle = stackiter(open(filename, 'rU'))
        for line in handle:
            while line.endswith('\\\n'):
                line = line[:-2] + next(handle)
    
            a = line.split(None, 1)
            if len(a) > 1 and a[0] == '_':
                line = a[1]
                a = line.split(None, 1)
    
            try:
                name = a[0]
                if name.startswith('/'):
                    line = line.lstrip()[1:]
                    raise
                name = cmd.kwhash.shortcut.get(name, name)
                kw = cmd.keyword[name]
                assert kw[4] != parsing.PYTHON
                func = kw[0]
            except:
                out.write(line)
                continue
            
            # PyMOL stuff without named python function
            if func.__name__ == '<lambda>' or name.startswith('@'):
                print('cmd.do(%s)' % (repr(line.strip())), file=out)
                continue
    
            # code blocks
            if name == 'python':
                for line in handle:
                    if line.split() == ['python', 'end']:
                        break
                    out.write(line)
                continue
            if name == 'skip':
                for line in handle:
                    if line.split() == ['skip', 'end']:
                        break
                continue
    
            # split args
            tok = ','
            if kw[4] == parsing.MOVIE:
                tok = kw[3]
                split_mx = 1
            else:
                split_mx = kw[4] - parsing.LITERAL
            if len(a) > 1:
                a[1] = parsing.split(a[1], '#', 1)[0] # strip of comment
                if split_mx < 0:
                    # strip of additional commands
                    a[1:] = parsing.split(a[1], ';', 1)
                    if len(a) > 2:
                        handle.push(a[2])
                if split_mx == 0:
                    args = [a[1]]
                else:
                    args = [i.strip() for i in parsing.split(a[1], tok, max(0, split_mx))]
            else:
                args = []
    
            # old syntax: set property=value
            if kw[4] == parsing.LEGACY and '=' in args[0]:
                args = [i.strip() for i in args[0].split('=', 1)] + args[1:]
    
            # register alias
            if name == 'alias':
                cmd.alias(*args)
    
            # use 'cmd' module if possible
            try:
                test = getattr(cmd, func.__name__)
                assert func == test
                module = 'cmd'
            except:
                module = func.__module__
    
            print('%s.%s(%s)' % (module, func.__name__, ', '.join(quote(args))), file=out)
    
    cmd.extend('pml2py', pml2py)
    
    # vi:expandtab:smarttab
    

Retrieved from "[https://pymolwiki.org/index.php?title=Pml2py&oldid=13247](https://pymolwiki.org/index.php?title=Pml2py&oldid=13247)"


---

## Process All Files In Directory

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Explanation

For a given directory with PDB files in it, the following code will output, for each PDB, the bound disulfide bond lengths like this: 
    
    
    1alk.pdb
      A CYS 168 SG
      A CYS 178 SG
        1.975
      A CYS 286 SG
      A CYS 336 SG
        1.995
      B CYS 168 SG
      B CYS 178 SG
        1.996
      B CYS 286 SG
      B CYS 336 SG
        2.032
    1btu.pdb
       CYS 42 SG
       CYS 58 SG
        2.039
       CYS 136 SG
       CYS 201 SG
        2.031
       CYS 168 SG
       CYS 182 SG
        2.001
       CYS 191 SG
       CYS 220 SG
        2.019
    ...
    

# Bound Disulfides
    
    
    from pymol import cmd
    from glob import glob
    
    for file in glob("*.pdb"):
        print file
        cmd.load(file,'prot')
        for a in cmd.index("elem s and bound_to elem s"):
            if cmd.select("s1","%s`%d"%a) and \
               cmd.select("s2","elem s and bound_to %s`%d"%a):
                if cmd.select("(s1|s2) and not ?skip"):
                    cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
                    print '   ',round(cmd.dist("tmp","s1","s2"),3)
                    cmd.select("skip","s1|s2|?skip")
        cmd.delete("all")
    

# All Sulfur Distances

Note that the above is for bonded sulfurs in disulfides. For all intra-cysteine gamma sulfur distances, you'd want to do something more like: 
    
    
    1alk.pdb
      A CYS 168 SG
      A CYS 178 SG
        1.975
      A CYS 168 SG
      A CYS 286 SG
        35.845
      A CYS 168 SG
      A CYS 336 SG
        35.029
      A CYS 168 SG
      B CYS 168 SG
        63.64
      A CYS 168 SG
      B CYS 178 SG
        63.775
      A CYS 168 SG
      B CYS 286 SG
        39.02
      A CYS 168 SG
      B CYS 336 SG
        39.314
    1btu.pdb
       CYS 42 SG
       CYS 58 SG
        2.039
       CYS 42 SG
       CYS 136 SG
    
    
    
    from pymol import cmd
    from glob import glob
    
    for file in glob("*.pdb"):
        print file
        cmd.load(file,'prot')
        for a in cmd.index("CYS/SG"):
            for b in cmd.index("CYS/SG"):
                if a[1]<b[1]:
                    cmd.select("s1","%s`%d"%a)
                    cmd.select("s2","%s`%d"%b)
                    if cmd.select("(s1|s2) and not ?skip"):
                        cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
                        print '   ',round(cmd.dist("tmp","s1","s2"),3)
                        cmd.select("skip","s1|s2|?skip")
        cmd.delete("all")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Process_All_Files_In_Directory&oldid=8031](https://pymolwiki.org/index.php?title=Process_All_Files_In_Directory&oldid=8031)"


---

## PythonTerminal

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL allows the execution of python commands from the PyMOL command line. It is very useful for both debugging, and for discovering new functions. 

  * Any expression that is not recognized as a PyMOL command is passed to the underlying python interpreter
  * To force a one-line expression to be executed as python, begin the line with a slash (/)
  * Use the [python](/index.php/Python "Python") command to input multi-line python code



## Examples
    
    
    # there is no "print" command in PyMOL, so this will go to the python interpreter
    print "Hello World (1)"
    
    # same, but force it to be python
    /print "Hello World (2)"
    
    # no lets trick this system by introducing a PyMOL command named "print"
    cmd.extend('print', lambda msg: sys.stdout.write("You gave me `%s`\n" % (msg)))
    
    # see what happens
    print "Hello World (3)"
    
    # this will still go to the python interpreter
    /print "Hello World (4)"
    

Retrieved from "[https://pymolwiki.org/index.php?title=PythonTerminal&oldid=9299](https://pymolwiki.org/index.php?title=PythonTerminal&oldid=9299)"


---

## Read PDB-String

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    delete all
    cmd.read_pdbstr("""HEADER    CREATED BY CONVERTPROSPECT              27-JAN-02   2tnf  \
    REMARK   1                                                                       \
    ATOM      1  N   PRO A   9       1.895  67.213 -38.182  1.00  0.00           N   \
    ATOM      2  CA  PRO A   9       1.703  68.680 -38.402  1.00  0.00           C   \
    ....
    ATOM   1153  C   GLY A 157       6.927  59.108 -38.901  1.00  6.00           C   \
    ATOM   1154  O   GLY A 157       6.700  59.292 -37.676  1.00  6.00           O   \
    TER    1155      GLY A 157                                                       \
    MASTER                                                                           \
    END                                                                              \
    ""","2tnfa")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Read_PDB-String&oldid=6362](https://pymolwiki.org/index.php?title=Read_PDB-String&oldid=6362)"


---

## Tiff2ccp4

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will convert tiff stacks to CCP4 maps for reading in PyMOL. 

_See also the[load_img_stack](/index.php/Load_img_stack "Load img stack") script for loading image stacks directly into PyMOL_. 

If someone knows how to determine the file's dimensions from the file itself (not the image size, but the physical dimension size) please let me know. 

Note: requires Python 2.7 or later (for argparse). 

# Usage
    
    
    #
    # convert input.tiff to output.ccp4 which has cell dimensions
    # x=40 A, y=80.5 A, Z=125.2 A
    #
    python tiff2ccp4.py input.tiff output.ccp4 -x 40 -y 80.5 -z 125.2
    
    #
    # print help
    #
    python tiff2ccp4.py -h
    

# The Code
    
    
    #
    # tiff2ccp4.py -- convert a TIFF stack to a CCP4 map
    #
    # Author: Jason Vertrees
    # Date  : 2011-09-16
    #
    # Notes : To get help, type "python tiff2ccp4.py -h"
    #
    import struct
    import sys
    import string
    
    # key   = field_type
    # value = # bytes for that field
    field_types = {                                                                                                                
        1 : 1,                                                                                                                     
        2 : 1,                                                                                                                     
        3 : 2,                                                                                                                     
        4 : 4,                                                                                                                     
        5 : 8,                                                                                                                     
        6 : 1,                                                                                                                     
        7 : 1,                                                                                                                     
        8 : 2,                                                                                                                     
        9 : 4,                                                                                                                     
        10: 8,                                                                                                                     
        11: 4,                                                                                                                     
        12: 8 }                                                                                                                    
                                                                                                                                   
    field_fmt = {                                                                                                                  
        1 : "B",                                                                                                                   
        2 : "c",                                                                                                                   
        3 : "H",                                                                                                                   
        4 : "I",                                                                                                                   
        5 : "II",                                                                                                                  
        6 : "i",                                                                                                                   
        7 : "B",                                                                                                                   
        8 : "h",                                                                                                                   
        9 : "i",                                                                                                                   
        10: "ii",                                                                                                                  
        11: "f",                                                                                                                   
        12: "d" }                                                                                                                  
                                                                                                                                   
    tag_types = {                                                                                                                  
        "NewSubfileType" : 254,                                                                                                    
        "ImageWidth"  : 256,  # # cols (number of pixels per scanline)                                                             
        "ImageLength" : 257,  # # rows (number of scanlines)                                                                       
        "BitsPerSample" : 258,                                                                                                     
        "PhotometricInterpretation" : 262,                                                                                         
        "ImageDescription" : 270,                                                                                                  
        "StripOffsets" : 273,                                                                                                      
        "SamplesPerPixel" : 277,                                                                                                   
        "RowsPerStrip" : 278,                                                                                                      
        "StripByteCounts" : 279,                                                                                                   
        "XResolution" : 282,                                                                                                       
        "YResolution" : 283,                                                                                                       
        "ZResolution" : 284,                                                                                                       
        "ResolutionUnit" : 296,                                                                                                    
        }
    
    
    class TIFFStack:
    
        def __init__(self,filename):
    
            self._filename = filename
            self._f = open(self._filename, 'rb')
            
            self._file_size = self.get_file_size()
    
            # read and store the file header
    
            self._header = self.get_header()
    
            # read all IFDs and store
    
            self._IFDs = self.get_IFDs()
    
    
        def __str__(self):
            s  = ""
            s += "[TIFFStack]\n"
            s += "Filename   : %s\n" % self._filename
            s += "File size  : %s\n" % self._file_size
            s += "Header     : %s\n" % self._header
            s += "IFDs       : %s\n" % self._IFDs
            return s
    
    
        def get_file_size(self):
            """
            file size in bytes
            """
            if not self._f:
                print "Invalid: Must have open file handle."
                return -1
    
            # less typing
    
            f = self._f
    
            # get file size
    
            curP = f.tell()
            f.seek(0,2)
            sz = f.tell()
    
            # return to previous location
    
            f.seek(curP)
    
            return sz
    
    
        def get_header(self):
            return TIFFHeader(self._f)
    
    
        def get_IFDs(self):
            return IFD_Store(self._f, self._header)
    
        def get_data(self):
    
            f = self._f
    
            # bits per sample
    
            bps = 258
    
            # StripOffsets
    
            offset = 273
    
            # byte counts per strip
    
            byte_count = 279
    
            sample_format = 339
    
            data = []
    
            ifds = self._IFDs._store
    
            for curIFD in ifds:
    
                curOffset = curIFD[offset]._val
                curLength = curIFD[byte_count]._val
                curBPS = curIFD[bps]._val
                bytesPerSample = curBPS / 8.0
                
                fld = self.get_field_fmt(curIFD)
    
                # use "B" for now; support 16-bit later using tag 339 (SampleFormat)
                unpackStr = self._header._e_flg + fld * int(curLength/bytesPerSample)
    
                f.seek(curOffset)
    
                data.extend(struct.unpack(unpackStr, f.read(curLength)))
    
            return data
            
        def get_field_fmt(self,ifd):
            """
            Determines the Python struct code given
            BitsPerSample and SampleFormat from the IFD
            """
            bits_per_sample, sample_fmt = 258, 339
    
            bps = ifd[bits_per_sample]._val
    
            fmt = None
    
            if sample_fmt in ifd.keys():
                fmt = ifd[sample_fmt]._val
            else:
                fmt = 1
    
            if bps==8:
                # 1-byte unsigned int
                if fmt==1:
                    return "B"
                # 1-byte signed
                elif fmt==2:
                    return "b"
            elif bps==16:
                # 2-byte unsigned
                if fmt==1:
                    return "H"
                elif fmt==2:
                    return "h"
    
        def get_size(self):
            fX, fY = 256, 257
            x = self._IFDs._store[0][fX]._val
            y = self._IFDs._store[0][fY]._val
            z = len(self._IFDs._store)
            
            return x, y, z
    
        def get_axes_scale(self):
            # for now
            return 1,1,1
    
    
        def asCCP4(self,filename,dimX=-1,dimY=-1,dimZ=-1):
            data = self.get_data()
            
            # pixels size x,y,z
    
            nX, nY, nZ = self.get_size()
    
            # dimension scaling
    
            if -1 in (dimX,dimY,dimZ):
                dimX, dimY, dimZ = self.get_axes_scale()
                dimX *= nX
                dimY *= nY
                dimZ *= nZ
    
            m, avg, M = min(data), sum(data)/len(data), max(data)
            
            outFile = open(filename, 'wb')
            
            outFile.write(struct.pack('i',nX)) # col
            outFile.write(struct.pack('i',nY)) # row
            outFile.write(struct.pack('i',nZ)) # section
            outFile.write(struct.pack('i',2))  # mode = 2
            
            outFile.write(struct.pack('i',1))  # number of first col 
            outFile.write(struct.pack('i',1))  # '' row
            outFile.write(struct.pack('i',1))  # '' section
            
            outFile.write(struct.pack('i',nX)) # nIntervals
            outFile.write(struct.pack('i',nY)) #
            outFile.write(struct.pack('i',nZ)) #
            
            outFile.write(struct.pack('f',dimX)) # length in X
            outFile.write(struct.pack('f',dimY)) # length in Y
            outFile.write(struct.pack('f',dimZ)) # length in Z
            
            outFile.write(struct.pack('f',90.)) # alpha 
            outFile.write(struct.pack('f',90.)) # beta
            outFile.write(struct.pack('f',90.)) # gamma
            
            outFile.write(struct.pack('i',1)) # axis cols
            outFile.write(struct.pack('i',2)) # axis rows
            outFile.write(struct.pack('i',3)) # axis section
            
            outFile.write(struct.pack('f', m)) # min density
            outFile.write(struct.pack('f', avg))  # max density
            outFile.write(struct.pack('f', M))  # mean density
            
            outFile.write(struct.pack('i',0)) # space gp ?
    
            # header info; blank for us
            for x in range(24,257): outFile.write(struct.pack('i',0))
    
            # assume constant data in file
            norm = 255.
            bps = tag_types["BitsPerSample"]
            max_bits = self._IFDs._store[0][bps]._val
            norm = float(2**max_bits-1.)
    
            # read/write data
            for x in data:
                outFile.write(struct.pack('f', x/norm))
    
            outFile.close()  
    
    
    class TIFFHeader:
        def __init__(self,fileHandle):
            self._endian, self._e_flg = self.get_endian(fileHandle)
            self._magic_number = self.get_magic_number(fileHandle)
            self._first_IFD = self.get_first_IFDOffset(fileHandle)
            
        def __str__(self):
            s  = "\n"
            s += "  [TIFF Header]\n"
            s += "  Endian        : %s\n" % self._endian
            s += "  Endian Flag   : %s\n" % self._e_flg
            s += "  Magic Number  : %s\n" % self._magic_number
            s += "  IFD[0] Offset : %s" % self._first_IFD
            return s
    
        # for struct.unpackx
        def _1byte(self,n=1):
            return self._e_flg + "C"*n
        def _2byte(self,n=1):
            return self._e_flg + "H"*n
        def _4byte(self,n=1):
            return self._e_flg + "I"*n
        def _IFDbyte(self,n=1):
            return self._e_flg + "HHII"*n
    
        def get_endian(self,fileHandle):
            f = fileHandle
            f.seek(0)
            code = struct.unpack("cc", f.read(2))
            code = "".join(code)
            flg = ""
            if code=="II":
                flg = "<"
            elif code=="MM":
                flg = ">"
            else:
                print "This file is not a valid TIFF (bad endian tag)."
                flg = "?"
            return code,flg
    
        def get_magic_number(self,fileHandle):
            f = fileHandle
            f.seek(2)
            # read the magic number
            idx = 0
            if self._endian == "II": idx = 1
            _42 = struct.unpack(self._2byte(), f.read(2))[0]
            if _42!=42:
                print "Error: Improperly formatted TIFF file (bad magic number)."
                return None
            return _42
    
        def get_first_IFDOffset(self,fileHandle):
            f=fileHandle
            f.seek(4)
            off = struct.unpack(self._4byte(), f.read(4))[0]
            return off
    
    
    class IFD_Store:
        def __init__(self,fileHandle,fileHeader):
            self._f = fileHandle
            self._header = fileHeader
            self._first_offset = self._header._first_IFD
            # array of IFDs
            self._store = self.read_ifds()
    
        def __str__(self):
            s  = "\n"
            s += "  [IFD_Store]\n"
            s += "  First Offset : %d\n" % self._first_offset
            s += "  Store :\n"
            for st in self._store:
                s += "\n\n  IFD Store =>\n"
                for k in st:
                    s += "    %s => %s" % (k,st[k])
            return s
    
        def read_ifds(self):
    
            f = self._f
    
            pos = self._first_offset
    
            ifds = []
    
            while pos!=0:
    
                ifds.append({})
                
                # get number of IFD_Entries
    
                f.seek(pos)
            
                # read number of IFDs
                num_entries = struct.unpack(self._header._2byte(), f.read(2))
                num_entries = num_entries[0]
    
                # read all into IFD[x]
                for x in range(num_entries):
                    # pull the current record from file
                    curParams = struct.unpack(self._header._IFDbyte(), f.read(12))
    
                    # format the data if necessary
                    if (self._header._e_flg==">" and sys.byteorder=="little") and \
                            curParams[0] not in (270,50838,50839):
                        scale = 32 - (8 * field_types[curParams[1]])
                        scaledData = curParams[3] >> scale
                    else:
                        scaledData = curParams[3]
    
    
                    ifds[-1][curParams[0]] = IFDEntry(curParams[0], curParams[1], curParams[2], scaledData)
    
                # read next offset
                pos = struct.unpack(self._header._4byte(), f.read(4))[0]
    
            return ifds
    
    class IFDEntry:
        def __init__(self,theTag=None,theType=None,theCount=None,theValue=None):
            self._tag = theTag
            self._type = theType
            self._count = theCount
            self._val = theValue
    
        def __str__(self):
            s  = "\n"
            s += "  [IFD_Entry]\n"
            s += " Tag    : %s\n" % self._tag
            s += " Type   : %s\n" % self._type
            s += " Count  : %s\n" % self._count
            s += " Val/Off: %s\n" % self._val
            return s
    
    if __name__=="__main__":
    
        """
        Running,
    
        python tiff2ccp4.py
    
        will convert all TIFF stacks in the current directory to 
        CCP4 maps.
        """
        import argparse
        from string import split
    
        parser = argparse.ArgumentParser(description="Convert a TIFF Stack to a CCP4 Map")
        parser.add_argument("input", type=str, help="input file name (usually .tif, .tiff)")
        parser.add_argument("output", type=str, help="output file name (usually .ccp4)")
        parser.add_argument('-x',"--x", help="length of x-dimension",default=-1.,type=float)
        parser.add_argument('-y',"--y", help="length of y-dimension",default=-1.,type=float)
        parser.add_argument('-z',"--z", help="length of z-dimension",default=-1.,type=float)
        
        args = parser.parse_args()
    
        s = TIFFStack(args.input)
        s.asCCP4(args.output, args.x, args.y, args.z)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Tiff2ccp4&oldid=12428](https://pymolwiki.org/index.php?title=Tiff2ccp4&oldid=12428)"


---

## XML-RPC server

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page describes a custom user implementation of an XML-RPC server. 

For PyMOL's built-in server, see [RPC](/index.php/RPC "RPC"). 

This server will enable any method in the 'cmd' module to be called remotely, as well as custom methods of the "pymol_interface" class (several examples of which are shown below). Other modules (such as util) could be wrapped using the same trick. 
    
    
    from pymol import cmd, util
    from threading import Thread
    from SimpleXMLRPCServer import SimpleXMLRPCServer
    import cPickle
    import os
    import re
    import string
    
    def pymol_startup () :
      print "Loading PyMOL XML-RPC extensions"
      server = pymol_interface()
    
    class pymol_xmlrpc_server (SimpleXMLRPCServer) :
      def __init__ (self, addr, interface) :
        self.interface = interface
        SimpleXMLRPCServer.__init__(self, addr, logRequests=0)
    
      def _dispatch (self, method, params) :
        if not self.interface.enable_xmlrpc :
          return -1
        result = -1
        func = None
        if hasattr(self.interface, method) :
          func = getattr(self.interface, method)
        elif hasattr(cmd, method) :
          func = getattr(cmd, method)
        if not callable(func) :
          print "%s is not a callable object" % method
        else :
          result = func(*params)
          if result is None :
            result = -1
        return result
    
    class pymol_interface (object) :
      def __init__ (self) :
        self.enable_xmlrpc = True
        # the port can be set via an environment variable - although it could just as easily be passed to __init__
        port = string.atoi(os.environ.get("PYMOL_XMLRPC_PORT", "9123"))
        self._server = pymol_xmlrpc_server(("localhost", port), self)
        t = threading.Thread(target=self._server.serve_forever)
        t.setDaemon(1)
        t.start()
        print "Started XML-RPC server on port %d" % port
    
      def close_all (self) :
        cmd.delete("*")
    
      def disable_all (self) :
        cmd.disable("*")
    
      def load_current_model_and_maps (self,
          pdb_file,
          fwt_map,
          delfwt_map) :
        model_name = os.path.basename(os.path.splitext(pdb_file)[0])
        cmd.load(pdb_file, model_name, state=1)
        cmd.load(fwt_map, "2fofc", state=1, format="ccp4")
        cmd.load(delfwt_map, "fofc", state=1, format="ccp4")
        cmd.isomesh("m1", "2fofc", 1.0, model_name, 5.0)
        cmd.isomesh("m2", "fofc", 3.0, model_name, 5.0)
        cmd.isomesh("m3", "fofc", -3.0, model_name, 5.0)
        cmd.color("marine", "m1")
        cmd.color("green", "m2")
        cmd.color("red", "m3")
    
      def setup_roving (self, f_map="2fofc", diff_map="fofc") :
        cmd.set("roving_detail", 20)
        cmd.set("roving_isomesh", 20)
        cmd.set("roving_origin", 1)
        cmd.set("roving_sticks", 0)
        cmd.set("roving_ribbon", 0)
        cmd.set("roving_lines", 0)
        cmd.set("roving_spheres", 0)
        cmd.set("roving_nb_spheres", 0)
        cmd.set("roving_polar_contacts", 0)
        cmd.set("roving_polar_cutoff", 0)
        cmd.set("stick_radius", 0.12)
        cmd.set("roving_map1_name", f_map)
        cmd.set("roving_map1_level", 1)
        cmd.set("roving_map2_name", diff_map)
        cmd.set("roving_map3_name", diff_map)
        cmd.set("roving_map2_level", 3)
        cmd.set("roving_map3_level", -3)
        cmd.refresh()
    
      def refresh_maps (self, f_map="2fofc", diff_map="fofc") :
        cmd.isomesh("rov_m1", f_map, 1.0, "center", 20)
        cmd.isomesh("rov_m2", diff_map, 3.0, "center", 20)
        cmd.isomesh("rov_m3", diff_map, -3.0, "center", 20)
        cmd.color("skyblue", "rov_m1")
        cmd.color("green", "rov_m2")
        cmd.color("red", "rov_m3")
    
      def show_selection (self, selection_string, zoom=True): #"True") :
        cmd.hide("sticks")
        util.cbay()
        try :
          cmd.show("sticks", selection_string)
        except Exception, e :
          return "Invalid selection."
        cmd.color("white", selection_string)
        if zoom :
          cmd.zoom(selection_string, buffer=5)
    
      def recenter (self, x, y, z) :
        view = list(cmd.get_view())
        view[12] = x
        view[13] = y
        view[14] = z
        cmd.set_view(view)
    
    if __name__ == "pymol" : # optional, if the module will be a command line argument to pymol
      pymol_startup()
    

The corresponding client code is shown below; not surprisingly, it is very simple to the example above. Because of the use of getattr() to retrieve methods in the 'cmd' API, we do not need to explicitly wrap cmd.load(), cmd.show(), or cmd.hide() for this code to work, and any other method in 'cmd' can be similarly accessed, as well as our custom methods. 
    
    
    from xmlrpclib import ServerProxy
    pymol = ServerProxy(uri="http://localhost:9123/RPC2")
    pymol.load("model.pdb")
    pymol.show("cartoon")
    pymol.hide("lines")
    pymol.load_current_model_and_maps("refined.pdb", "2mFo-DFc.ccp4", "mFo-DFc.ccp4")
    

Retrieved from "[https://pymolwiki.org/index.php?title=XML-RPC_server&oldid=12798](https://pymolwiki.org/index.php?title=XML-RPC_server&oldid=12798)"


---

