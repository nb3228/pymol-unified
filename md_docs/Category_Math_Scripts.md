# Category: Math Scripts

## Axes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Axes with text labels
    
    
    # axes.py
    from pymol.cgo import *
    from pymol import cmd
    from pymol.vfont import plain
    
    # create the axes object, draw axes with cylinders coloured red, green,
    #blue for X, Y and Z
    
    obj = [
       CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
       CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
       CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
       ]
    
    # add labels to axes object (requires pymol version 0.8 or greater, I
    # believe
    
    cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    
    # then we load it into PyMOL
    cmd.load_cgo(obj,'axes')
    

## Axes with nice cones

[![Axes demo.png](/images/e/e9/Axes_demo.png)](/index.php/File:Axes_demo.png)

This script draws a simple cartesian coordinate system. 
    
    
    from pymol.cgo import *
    from pymol import cmd
    
    w = 0.06 # cylinder width 
    l = 0.75 # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter
    
    obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
           CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
           CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    
    cmd.load_cgo(obj, 'axes')
    

## Axes which always stay in the lower left corner
    
    
    from pymol import cmd
    from chempy import cpv
    
    class PutCenterCallback(object):
        prev_v = None
    
        def __init__(self, name, corner=0):
            self.name = name
            self.corner = corner
            self.cb_name = cmd.get_unused_name('_cb')
    
        def load(self):
            cmd.load_callback(self, self.cb_name)
    
        def __call__(self):
            if self.name not in cmd.get_names('objects'):
                import threading
                threading.Thread(None, cmd.delete, args=(self.cb_name,)).start()
                return
    
            v = cmd.get_view()
            if v == self.prev_v:
                return
            self.prev_v = v
    
            t = v[12:15]
    
            if self.corner:
                vp = cmd.get_viewport()
                R_mc = [v[0:3], v[3:6], v[6:9]]
                off_c = [0.15 * v[11] * vp[0] / vp[1], 0.15 * v[11], 0.0]
                if self.corner in [2,3]:
                    off_c[0] *= -1
                if self.corner in [3,4]:
                    off_c[1] *= -1
                off_m = cpv.transform(R_mc, off_c)
                t = cpv.add(t, off_m)
    
            z = -v[11] / 30.0
            m = [z, 0, 0, 0, 0, z, 0, 0, 0, 0, z, 0, t[0] / z, t[1] / z, t[2] / z, 1]
            cmd.set_object_ttt(self.name, m)
    
    def axes(name='axes'):
        '''
    DESCRIPTION
    
        Puts coordinate axes to the lower left corner of the viewport.
        '''
        from pymol import cgo
    
        cmd.set('auto_zoom', 0)
    
        w = 0.06 # cylinder width
        l = 0.75 # cylinder length
        h = 0.25 # cone hight
        d = w * 1.618 # cone base diameter
    
        obj = [cgo.CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
               cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
               cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
               cgo.CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
               cgo.CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
               cgo.CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    
        PutCenterCallback(name, 1).load()
        cmd.load_cgo(obj, name)
    
    cmd.extend('axes', axes)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Axes&oldid=12542](https://pymolwiki.org/index.php?title=Axes&oldid=12542)"


---

## BbPlane

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/bbPlane.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/bbPlane.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate"), [Blaine Bell](/index.php?title=User:Bell&action=edit&redlink=1 "User:Bell \(page does not exist\)") and [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO plane between the backbone atoms of two neighboring residues. This is to show the planarity of the atoms. The image style this is meant to represent can be found many places, like "Introduction to Protein Structure" by Branden and Tooze (2nd ed. pp. 8). 

  * [![Close up of planar atoms](/images/8/8d/BbPlane3.png)](/index.php/File:BbPlane3.png "Close up of planar atoms")

Close up of planar atoms 

  * [![A few more](/images/8/89/BbPlane1.png)](/index.php/File:BbPlane1.png "A few more")

A few more 

  * [![Global view](/images/b/b8/BbPlane2.png)](/index.php/File:BbPlane2.png "Global view")

Global view 




# Examples
    
    
    # download the source and save as bbPlane.py
    run bbPlane.py
    fetch 1cll
    # make planes for residues 4-9
    bbPlane i. 4-10
    

Retrieved from "[https://pymolwiki.org/index.php?title=BbPlane&oldid=13888](https://pymolwiki.org/index.php?title=BbPlane&oldid=13888)"


---

## BiologicalUnit

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Note:** PyMOL 1.8 can load biological units from mmCIF files with the **[assembly](/index.php/Assembly "Assembly")** setting. 

This script can be used to re-create biological units for proteins. (This was created as a workaround of PyMOL's semi-functioning [Symexp](/index.php/Symexp "Symexp") command.) It's also a fun script to play with for learning about symmetry. 

  * [![single unit, before running this program](/images/a/a9/Before.png)](/index.php/File:Before.png "single unit, before running this program")

single unit, before running this program 

  * [![after the expansion of the 60 units](/images/c/c6/After.png)](/index.php/File:After.png "after the expansion of the 60 units")

after the expansion of the 60 units 

  * [![Example creating the biological unit for PDB 1RMV created from fiber diffraction.](/images/7/75/1rmv_fiber.png)](/index.php/File:1rmv_fiber.png "Example creating the biological unit for PDB 1RMV created from fiber diffraction.")

Example creating the biological unit for PDB 1RMV created from fiber diffraction. 




## Contents

  * 1 Usage
  * 2 The Code
  * 3 Notes
  * 4 See Also



# Usage
    
    
    load /path/to/some/pdbFile.pdb
    symMat = readSymmetry("/path/to/some/pdbFile.pdb","pdbFile")
    biologicalUnit("mates", "pdbFile", symMat)
    

# The Code
    
    
    #
    # Jason Vertrees <Jason-dot-Vertrees-at-schrodinger_dot_com>, 2010.
    #
    import pymol
    from pymol import cmd
    
    def readSymmetry(inFile, verbose=None):
      """
      This function will read "inFile" and glean the
      symmetry operations, if any, from it.
      
      PARAMS
        inFile
          (string) path to PDB file
          
        verbose
          (boolean) if verbose is not None, print more
          
      RETURNS
        matrix
          Array of lists.  One 16-element list per symmetry operation.  Feed this matrix
          into manualSymExp in order to make the other symmetry mates in the biological unit
      """
      # a remark-350 lines has:
      # REMARK 350 BIOMTn TAGn X Y Z Tx
      REM, TAG, BIOMT, OPNO, X, Y, Z, TX = range(8)
      
      thePDB = open(inFile, 'rb').readlines()
      
      matrices = []
      curTrans = -1
      
      # The transformation is,
      # output = U*input + Tx
      
      for l in thePDB:
        tokens = l.split()
        if len(tokens)!=8:
          continue
        if tokens[REM]=="REMARK" and tokens[TAG]=="350" and tokens[BIOMT].startswith("BIOMT"):
          if tokens[OPNO]!=curTrans:
            # new transformation matrix
            matrices.append([])
          
          matrices[-1].append( map( lambda s: float(s), tokens[X:]))
          curTrans = tokens[OPNO]
    
      if verbose!=None:
        print "Found %s symmetry operators in %s." % (len(matrices), inFile)
      return matrices
    
    
    def biologicalUnit(prefix, objSel, matrices ):
      """
      Manually expands the object in "objSel" by the symmetry operations provided in "matrices" and
      prefixes the new objects with "prefix".
      
      PARAMS
        prefix
          (string) prefix name for new objects
        
        objSel
          (string) name of object to expand
          
        matrices
          (list of 16-element lists) array of matrices from readSymmetry
          
        RETUNRS
          None
      
        SIDE EFFECTS
          Creates N new obects each rotated and translated according to the symmetry operators, where N
          equals len(matrices).
      """
      for m in matrices:
        n = cmd.get_unused_name(prefix)
        cmd.create(n, objSel)
        s1 = "%s + (x*%s + y*%s + z*%s)" % (m[0][3], m[0][0], m[0][1], m[0][2])
        s2 = "%s + (x*%s + y*%s + z*%s)" % (m[1][3], m[1][0], m[1][1], m[1][2])
        s3 = "%s + (x*%s + y*%s + z*%s)" % (m[2][3], m[2][0], m[2][1], m[2][2])
        cmd.alter_state(1, n, "(x,y,z) = (%s, %s, %s)" % (s1, s2, s3) )
    

# Notes

This is slow compared to [Symexp](/index.php/Symexp "Symexp"); use the above for learning, playing and when [Symexp](/index.php/Symexp "Symexp") doesn't work as advertised. 

# See Also

  * [BiologicalUnit/Quat](/index.php/BiologicalUnit/Quat "BiologicalUnit/Quat") (alternative implementation)
  * [Symexp](/index.php/Symexp "Symexp")
  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [PDB Tutorial Biol. Units](http://www.rcsb.org/pdb/static.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html)
  * [Wikipedia article](http://en.wikipedia.org/wiki/Fiber_diffraction)
  * [assembly](/index.php/Assembly "Assembly")



Retrieved from "[https://pymolwiki.org/index.php?title=BiologicalUnit&oldid=12433](https://pymolwiki.org/index.php?title=BiologicalUnit&oldid=12433)"


---

## BiologicalUnit/Quat

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**Note:** PyMOL 1.8 can load biological units from mmCIF files with the **[assembly](/index.php/Assembly "Assembly")** setting. 

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [(biomolecule) psico.xtal (biomolecule)](http://pymol.org/psicoredirect.php?psico.xtal)  
  
This script does more or less the same as [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit"). It creates the **quat** ernary structure (BIOMOLECULE 1 assembly) from the REMARK 350 header. 

This script is convenient to use because it searches automatically for the PDB file in the current directory, in [fetch_path](/index.php/Fetch_Path "Fetch Path") and (if available) in your local PDB mirror. 

_Also available from[psico](/index.php/Psico "Psico"), but the command is called **biomolecule** instead of **quat**. _

# Example Usage
    
    
    fetch 3bw1, type=pdb
    quat 3bw1
    as cartoon
    

# The Code
    
    
    '''
    (c) 2010-2011 Thomas Holder, MPI for Developmental Biology
    
    Module for reading REMARK records from PDB files and in particular
    generate quaterny structure from REMARK 350.
    '''
    
    import sys, os
    from pymol import cmd, stored
    
    local_mirror_divided = '/mnt/bio/db/pdb.divided'
    
    def pdbremarks(filename):
        '''
        Read REMARK lines from PDB file. Return dictionary with remarkNum as key
        and list of lines as value.
        '''
        remarks = dict()
        if not isinstance(filename, basestring):
            f = filename
        elif filename[-3:] == '.gz':
            import gzip
            f = gzip.open(filename)
        else:
            f = open(filename)
        for line in f:
            recname = line[0:6]
            if recname == 'REMARK':
                num = int(line[7:10])
                lstring = line[11:]
                remarks.setdefault(num, []).append(lstring)
        return remarks
    
    def quat350(rem350):
        '''
        Get transformation matrices for biomolecule 1 from REMARK 350.
        '''
        biomt = dict()
        chains = tuple()
        seenbiomolecule = False
        for line in rem350:
            if line.startswith('BIOMOLECULE:'):
                if seenbiomolecule:
                    break
                seenbiomolecule = True
            elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
                chains = tuple(chain.strip() for chain in line[30:].split(','))
            elif line.startswith('                   AND CHAINS:'):
                chains += tuple(chain.strip() for chain in line[30:].split(','))
            elif line.startswith('  BIOMT'):
                row = int(line[7])
                num = int(line[8:12])
                vec = line[12:].split()
                vec = map(float, vec)
                biomt.setdefault(chains, dict()).setdefault(num, []).extend(vec)
        return biomt
    
    def quat(name=None, filename=None, prefix=None, quiet=0):
        '''
    DESCRIPTION
    
        Read REMARK 350 from `filename` and create biological unit
        (quaternary structure)
    
    USAGE
    
        quat [name [, filename [, prefix]]]
    
    ARGUMENTS
    
        name = string: name of object and basename of PDB file, if
        filename is not given {default: first loaded object}
    
        filename = string: file path {default: <name>.pdb}
    
        prefix = string: prefix for new objects {default: <name>}
    
    EXAMPLE
    
        fetch 1rmv, type=pdb
        quat 1rmv
        '''
        quiet = int(quiet)
        if name is None:
            name = cmd.get_object_list()[0]
        if prefix is None:
            prefix = name
        if filename is None:
            candidates = [
                '%s.pdb' % (name),
                '%s/%s.pdb' % (cmd.get('fetch_path'), name),
                '%s/%s/pdb%s.ent.gz' % (local_mirror_divided, name[1:3], name),
            ]
            for filename in candidates:
                if os.path.exists(filename):
                    break
            else:
                print 'please provide filename'
                return
            if not quiet:
                print 'loading from %s' % (filename)
        remarks = pdbremarks(filename)
        if 350 not in remarks:
            print 'There is no REMARK 350 in', filename
            return
        quat = quat350(remarks[350])
        for chains in quat:
            matrices = quat[chains]
            for num in matrices:
                mat = matrices[num][0:12]
                mat.extend([0,0,0,1])
                copy = '%s_%d' % (prefix, num)
                if not quiet:
                    print 'creating %s' % (copy)
                cmd.create(copy, '/%s//%s' % (name, '+'.join(chains)))
                cmd.alter(copy, 'segi="%d"' % (num))
                cmd.transform_object(copy, mat)
        cmd.disable(name)
        cmd.group('%s_quat' % (prefix), '%s_*' % (prefix))
    
    cmd.extend('quat', quat)
    
    # vi:expandtab:smarttab
    

# See Also

  * [BiologicalUnit](/index.php/BiologicalUnit "BiologicalUnit")
  * [Supercell](/index.php/Supercell "Supercell")
  * [assembly](/index.php/Assembly "Assembly")



Retrieved from "[https://pymolwiki.org/index.php?title=BiologicalUnit/Quat&oldid=12743](https://pymolwiki.org/index.php?title=BiologicalUnit/Quat&oldid=12743)"


---

## Bounding Box

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    import math
    from pymol import querying
    from pymol.cgo import *
    from pymol import cmd
    
    #NOTE!! : These packages (numarray, Scientific) must be present in pymol's 
    #'$PYMOLDIR/py23/lib/python2.3/site-packages/' directory!!
    # OR..if you are using Mac PyMOL 0.99beta19 or later then the system installs will be used
    from numarray import *
    from numarray.ma import average
    from numarray import linear_algebra as la
    
    from Scientific.Geometry import Vector, Tensor, Transformation
    
    def boundingBox(selectionName, colourRGB=[1,1,1]):
            """
            The main function to call : 
    
                    run "box.py"
                    boundingBox("peptide")
    
            Should make a box around "peptide" (assuming it exists!). For a different colour use:
    
                    boundingBox("peptide", colourRGB=[1, 0, 1])
    
            Or whatever. The box should be a cgo called 'peptide-box' (for this example).
            """
            model = querying.get_model(selectionName)
            coords = model.get_coord_list()
    
            #find the least square plane through the coordinates
            eigenvalues, eigenvectors, centroid = findPlaneWithEigenVectors(coords)
            normal = eigenvectors[eigenvalues.argmin()]
            eigenvectors.remove(normal)
    
            #orient the axes correctly
            x, y, normal = orientAxes(eigenvectors[0], eigenvectors[1], normal)
    
            #determine the dimensions and the structure's orientation wrt the origin
            minDimensions, rotation = findMinDimensionsAndRotation(coords, centroid, x, y, normal)
    
            #'create' the box(IE: make the corners) and 'draw' it (IE: make a cgo)
            box = makeBox(minDimensions, rotation, centroid)
            drawBox(box, selectionName, colourRGB)
    
    def findPlaneWithEigenVectors(coords):
            centroid = average(coords)
            coords -= centroid
            B = matrixmultiply(transpose(coords), coords)
            eigenvalues, eigenvectors = la.eigenvectors(B)
            #return eigenvalues, [Vector(e) for e in eigenvectors], Vector(centroid)
            return eigenvalues, [Vector([i for i in e]) for e in eigenvectors], Vector(centroid) #not sure why I had to make this change!
    
    def orientAxes(x, y, z):
            XcrossY = x.cross(y)
            #ZXY = around(math.degrees(z.angle(XcrossY)))
            ZXY = int(around(math.degrees(z.angle(XcrossY)))) #again, a bit of a hack!
            if (ZXY == 180): x, y = y, x
            return x, y, z
    
    def makeBox(dimensions, rotation, centroid):
            x, y, z = dimensions
            v = [[]] * 8
    
            #make a cuboid with the lower corner on the origin
            v[0] = [0, 0, 0]                # [0, 0, 0]
            v[1] = [x, 0, 0]                # [1, 0, 0]
            v[2] = [x, y, 0]                # [1, 1, 0]
            v[3] = [x, 0, z]                # [1, 0, 1]
            v[4] = [0, y, 0]                # [0, 1, 0]
            v[5] = [0, 0, z]                # [0, 0, 1]
            v[6] = [0, y, z]                # [0, 1, 1]
            v[7] = [x, y, z]                # [1, 1, 1]
    
            #move to the origin, THEN move to the centroid of the points, then rotate
            translationToOrigin = Transformation.Translation(-Vector(x/2, y/2, z/2))
            translationToCentroid = Transformation.Translation(centroid)
            transform = translationToCentroid * rotation * translationToOrigin
    
            #use the Transformation to transform the corners of the box
            v = [transform(Vector(i)) for i in v]
    
            bot =  [v[0], v[1], v[2], v[4]] # O, x, xy, y
            top =  [v[7], v[3], v[5], v[6]] # xyz, xz, z, yz
            minL = [v[0], v[4], v[6], v[5]] # O, y, yz, z
            minR = [v[0], v[5], v[3], v[1]] # O, z, xz, x
            maxL = [v[4], v[2], v[7], v[6]] # y, xy, xyz, yz
            maxR = [v[3], v[1], v[2], v[7]] # xz, x, xy, xyz
            box = [bot, minR, minL, maxR, maxL, top]
    
            return box
    
    def drawBox(box, name, colourRGB):
            boxObj = []
            for side in box:
                    boxObj.append(BEGIN)
                    boxObj.append(LINE_STRIP)
                    boxObj.append(COLOR)
                    boxObj.extend(colourRGB)
                    for point in side:
                            boxObj.append(VERTEX)
                            boxObj.extend(point)
                    boxObj.append(END)
    
            cmd.set('auto_zoom', 0)
            cmd.load_cgo(boxObj, "%s-box" % name)
            cmd.set('auto_zoom', 1)
    
    def findMinDimensionsAndRotation(coords, centroid, x, y, z):
            O = Vector(0, 0, 0)
            X = Vector(1, 0, 0)
            Y = Vector(0, 1, 0)
            Z = Vector(0, 0, 1)
    
            #Create a Transformation t = |x, y, z| . |X, Y, Z| ^ -1
            mfinal = array([array(X), array(Y), array(Z)])
            morig = array([array(x), array(y), array(z)])
            rotmatrix = matrixmultiply(morig, transpose(mfinal))
            tTotal = Transformation.Rotation(Tensor(rotmatrix))
    
            #Transform the coordinates and find the min, max dimensions
            coordArray = array([array(tTotal(Vector(c))) for c in coords])
            minDimensions = [max(coordArray[:,i]) - min(coordArray[:,i]) for i in range(3)]
    
            #Now, compose the inverse rotation used to move the bounding box to the right place
            tForward = Transformation.Rotation(Tensor(matrixmultiply(mfinal, transpose(morig))))
    
            return minDimensions, tForward
    

Retrieved from "[https://pymolwiki.org/index.php?title=Bounding_Box&oldid=6395](https://pymolwiki.org/index.php?title=Bounding_Box&oldid=6395)"


---

## Cart to frac

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will convert the real space orthonormal cartesian coordinates of a selection to the fractional coordinates given the loaded cell symmetry. 

  * Thanks to [Wikipedia](http://en.wikipedia.org/wiki/Fractional_coordinates).



# Example
    
    
    # load the script
    
    run cart_to_frac.py
    
    # fetch a protein and test
    
    fetch 1rx1, async=0
    
    # get the coordinates for the organic small
    # molecule as fractional coordinates
    
    print cart_to_frac("org")
    

# The Code
    
    
    from pymol import cmd
    
    def cart_to_frac(objSel,quiet=0,_self=cmd):
        """
        Returns an array of fractional atomic coordinates
        for a given object or selection.
    
        PARAMS
          objSel -- any object or selection
    
          quiet -- suppress output (default, no)
    
          _self -- core CMD object; or none
    
        RETURNS
          Python list of fractional coordinates
    
        NOTES/EXAMPLES
          cart_to_frac org
    
          x = cart_to_frac("solvent", quiet=1)
        """
        import numpy
        from numpy import cos, sin, sqrt
    
        a2r = numpy.pi / 180.0
    
        # get the model and coordinates
    
        m = _self.get_model(objSel)
        cart_coord = numpy.matrix(m.get_coord_list())
    
        # get the symmetry information
        try:
            a,b,c,alpha,beta,gamma,gp = _self.get_symmetry(objSel)
        except:
            print "Error-Failed to get symmetry. Please ensure you have a"
            print "valid object with proper crystal symmetry loaded."
            return None
    
        # convert to radians
    
        alpha = a2r * alpha
        beta  = a2r * beta
        gamma = a2r * gamma
        
        # (scaled) volume of the cell
        
        v = sqrt(1 -cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma) + 2*cos(alpha)*cos(beta)*cos(gamma))
    
        tmat = numpy.matrix( [
          [ 1.0 / a, -cos(gamma)/(a*sin(gamma)), (cos(alpha)*cos(gamma)-cos(beta)) / (a*v*sin(gamma))  ],
          [ 0.0,     1.0 / (b*sin(gamma)),         (cos(beta) *cos(gamma)-cos(alpha))/ (b*v*sin(gamma))  ],
          [ 0.0,     0.0,                        sin(gamma) / (c*v)                                    ] ]
          )
    
        r = cart_coord * tmat.T
    
        if not quiet:
            for (x,y,z) in r.tolist():
                print '%8.5f %8.5f %8.5f' % (x,y,z)
    
        # return the Nx3 results
        return r
        
    cmd.extend("cart_to_frac", cart_to_frac)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cart_to_frac&oldid=12442](https://pymolwiki.org/index.php?title=Cart_to_frac&oldid=12442)"


---

## Center of mass

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/center_of_mass.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/center_of_mass.py)  
Author(s)  | [Henschel](/index.php?title=User:Henschel&action=edit&redlink=1 "User:Henschel \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Examples
  * 4 See Also



### Description

This script calculates the true center-of-mass (COM) or the center-of-geometry (COG) for a given selection and returns the x, y, z values in the form of a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") (rather than a CGO sphere). The benefit of using a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") is that it can be selected and used in calculations. In addition, this script also iteratively goes through all states of a selection if more than one state exists and appends the corresponding COM/COG values as states into the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). The script itself is quite simple yet robust enough to be applied in different settings. As well, the calculation of the COM/COG is handled independently from the formation of the [Pseudoatom](/index.php/Pseudoatom "Pseudoatom") and can be called as an independent function where applicable. 

### Usage
    
    
    com selection [,state=None [,mass=None [,object=None]]]
    

  


### Examples
    
    
    import center_of_mass
    fetch 1c3y, finish=1, multiplex=0
    
    com 1c3y, state=1
    #Create a pseudoatom representing the 1c3y COG and store it as "1c3y_COM"
    #The "1c3y_COM" object will contain 1 state only
    
    com 1c3y, state=1, object=COG
    #Create a pseudoatom representing the 1c3y COG and store it as "COG"
    #The "COG" object will contain 1 state only
    
    com 1c3y, state=1, object=COM, mass=1
    #Create a pseudoatom representing the 1c3y COM and store it as "COM"
    #The "COM" object will contain 1 state only
    
    com 1c3y, object=COM, mass=1
    #Create a single pseudoatom containing the COM for each state found in 1c3y and store it as "COM"
    #The "COM" object will contain MULTIPLE states!
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass") (new command in PyMOL 1.7.2)
  * [COM](/index.php/COM "COM")
  * [get_extent](/index.php/Get_extent "Get extent")
  * [get_position](/index.php/Get_position "Get position")



Retrieved from "[https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892](https://pymolwiki.org/index.php?title=Center_of_mass&oldid=13892)"


---

## Cgo arrow

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/cgo_arrow.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_arrow.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO arrow between two picked atoms. 

## Usage
    
    
    cgo_arrow [ atom1 [, atom2 [, radius [, gap [, hlength [, hradius [, color [, name ]]]]]]]]
    

## Example
    
    
    run cgo_arrow.py
    
    fetch 1rx1, async=0
    preset.pretty("*")
    
    cgo_arrow [34.9, 68.4, 19.1], A/164/O3X, gap=1.0
    

[![Cgo arrow example.png](/images/1/17/Cgo_arrow_example.png)](/index.php/File:Cgo_arrow_example.png)

Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_arrow&oldid=13894](https://pymolwiki.org/index.php?title=Cgo_arrow&oldid=13894)"


---

## Cgo grid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/cgo_grid.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cgo_grid.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![cgo_grid creates flowing mesh objects](/images/2/27/Cgo_grid.gif)](/index.php/File:Cgo_grid.gif "cgo_grid creates flowing mesh objects")

## Contents

  * 1 About cgo_grid
  * 2 Usage
  * 3 Arguments
  * 4 Examples
  * 5 SEE ALSO



## About cgo_grid

**cgo_grid** will generate a flowing mesh object using the points provided or the current view. By default is will generate a flowing membrane. The shape is affected substantially by the arguments!  


  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage

**cgo_grid** has many arguments, but not all need to be set necessarily (see arguments or examples). 
    
    
    cgo_grid [ pos1 [, pos2 [, pos3 [, length_x [, length_z [, npoints_x [, npoints_z
    [, nwaves_x [, nwaves_z [, offset_x [, offset_z [, gain_x [, gain_z
    [, thickness [, color [, nstates [, startframe [, endframe
    [, mode [, view [, name [, quiet ]]]]]]]]]]]]]]]]]]]]]]
    

|   
---|---  
  
  


## Arguments
    
    
        pos1 = single atom selection (='pk1') or list of 3 floats {default: [0,0,0]}
    
        pos2 = single atom selection (='pk2') or list of 3 floats {default: [1,0,0]}
    
        pos3 = single atom selection (='pk3') or list of 3 floats {default: [0,0,1]}
    
        --> the plane is defined by pos1 (origin) and vectors to pos2 and pos3, respectively
    
        length_x = <float>: length of membrane {default: 30}
        length_z = <float>: length of membrane {default: ''} # same as length_x
    
        npoints_x = <int>: number of points(lines) along x-direction
                    {default: ''} #will be set to give a ~1 unit grid
        npoints_z = <int>: number of points(lines) along z-direction
                    {default: ''} #will be set to give a ~1 unit grid
                    {minimum: 1 # automatic}
    
        nwaves_x =   <float>: number of complete sin waves along object x-axis
                     {default: 2}
        nwaves_z =  <float>: number of complete sin waves along object z-axis
                    {default: ''} # same as nwaves_x
                    define separately to adjust number of waves in each direction
    
        offset_x = <float> phase delay (in degrees) of sin wave in x-axis
                 can be set to affect shape and starting amplitude {default: 0}
        offset_z = <float> phase delay (in degrees) of sin wave in z-axis
                 can be set to affect shape and starting amplitude
                 {default: ''} # same as  offset_x
        offset_x and offset_z can be used together to phase
        otherwise identical objects
    
        gain_x = <float>: multiplication factor for y-amplitude for x-direction
                 {default: 1}
        gain_z = <float>: multiplication factor for y-amplitude for z-direction
                 {default: ''} #=gain_x
    
        thickness = <float>: line thickness {default: 2}
    
        color = color name <string> (e.g. 'skyblue') OR
                rgb-value list of 3 floats (e.g. [1.0,1.0,1.0]) OR
                {default: ''} // opposite of background
                input illegal values for random coloring
    
        nstates =  <int>: number of states; {default: 60}
                   this setting will define how many states
                   the object will have (per wave) and how fluent and fast the
                   animation will be.
                   Higher values will promote 'fluent' transitions,
                   but decrease flow speed.
                       Note: Frame animation cycles thought the states one at a time
                       and needs to be set accordingly. Can also be used to phase
                       otherwise identical objects.
                   Set to 1 for static object {automatic minimum}
    
        startframe: specify starting frame <int> or set (='') to use current frame
                    set to 'append' to extend movie from the last frame {default: 1}
          endframe: specify end frame <int> or set (='') to use last frame
                    if 'append' is used for startframe,
                    endframe becomes the number of frames to be appended instead
                    {default: 1}
                    Note: if start- and endframe are the same, movie animation will
                    be skipped, the object will be loaded and can be used afterwards
    
        mode: defines positioning {default: 0}:
        0: pos1 is center
        1: pos1 is corner
    
        view {default: 0}:
        '0': off/ uses provided points to create CGO
        '1': overrides atom selections and uses current orienatation for positioning
             - pos1 = origin/center
             - pos2 = origin +1 in camera y
             - pos3 = origin +1 in camera z
    
        name: <string> name of cgo object {default: ''} / automatic
    
        quiet: <boolean> toggles output
    

| Concept sketch:  
[![cgo_grid concept sketch](/images/d/d1/Cgo_grid.png)](/index.php/File:Cgo_grid.png "cgo_grid concept sketch")  
---|---  
  
  


## Examples

The behaviour or shape of the cgo_grid object are substantially influenced by the arguments 
    
    
    delete all
    set_view (\
         0.263772517,   -0.113038681,    0.957937598,\
        -0.040853567,    0.990910411,    0.128179103,\
        -0.963716805,   -0.072944686,    0.256756991,\
         0.000000000,    0.000000000, -131.816467285,\
         0.000000000,    0.000000000,    0.000000000,\
       -50.008331299,  353.641235352,  -20.000000000 )
    
    #membrane
    cgo_grid color=blue
    
    #swimming worm, random color
    cgo_grid \
    pos1=[0,-5,0], pos2=[1,-5,1], pos3=[0,-5,1],\
    length_x=15,\
    npoints_z=1,\
    gain_x=2,\
    gain_z=0,\
    thickness=20,\
    color=3,\
    mode=1,\
    name="worm"
    
    #Moving Ladder
    cgo_grid \
    length_x=15,\
    pos1=[0,10,0], pos2=[0,10,1], pos3=[0,9,0],\
    npoints_x=2, npoints_z=30,\
    name="ladder"
    
    #Roof
    cgo_grid \
    nstates=1,\
    npoints_x=15,\
    npoints_z=15,\
    gain_x=20,\
    gain_z=20,\
    nwaves_x=0.5,\
    thickness=5,\
    color=cyan,\
    name="roof"
    
    #Boxes
    cgo_grid \
    pos1=[0,-10,0], pos2=[1,-10,0], pos3=[0,-10,1],\
    nstates=1,\
    npoints_x=50,\
    npoints_z=50,\
    nwaves_x=0,\
    color=[0.00 , 0.53 , 0.22],\
    thickness=5,\
    name="bottom"
    
    cgo_grid \
    nstates=1,\
    npoints_x=2,\
    npoints_z=2,\
    nwaves_x=0,\
    color=gray60,\
    thickness=10,\
    name="top"
    
    cgo_grid \
    pos1=[-15,-10,15], pos2=[-14,-10,15], pos3=[-15,-9,15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=-2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="front"
    
    cgo_grid \
    pos1=[-15,-10,-15], pos2=[-14,-10,-15], pos3=[-15,-9,-15],\
    nstates=1,\
    npoints_x=5,\
    npoints_z=5,\
    gain_x=0,\
    gain_z=2,\
    length_z=10,\
    nwaves_x=0.5,\
    color=gray60,\
    thickness=5,\
    mode=1,\
    name="back"
    
    set ray_trace_frames, 0
    set movie_loop,1
    mplay
    
    # play around with the ARGUMENTS! :-)
    

|   
---|---  
  
  


## SEE ALSO

  * [CgoCircle](/index.php/CgoCircle "CgoCircle")



Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938](https://pymolwiki.org/index.php?title=Cgo_grid&oldid=13938)"


---

## CGO Text

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/5/5f/Cgo_text.png)](/index.php/File:Cgo_text.png)

[](/index.php/File:Cgo_text.png "Enlarge")

CGO Text Example from following example
    
    
    # draw text using cgo
    from pymol import cmd
    from pymol.cgo import *
    from pymol.vfont import plain
    
    cgo = []
    
    axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]
    
    pos = [0.0,0.0,0.0]
    wire_text(cgo,plain,pos,'Hello World',axes)
    
    pos = [0.0,-3.0,0.0]
    cyl_text(cgo,plain,pos,'Hello Universe',0.10,axes=axes)
    
    cmd.set("cgo_line_radius",0.03)
    cmd.load_cgo(cgo,'txt')
    cmd.zoom("all",2.0)
    

Retrieved from "[https://pymolwiki.org/index.php?title=CGO_Text&oldid=10878](https://pymolwiki.org/index.php?title=CGO_Text&oldid=10878)"


---

## CgoCircle

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script will create a CGO circle with the origin at the specified X,Y,Z coordinates. Also, you can specify the radius and the colors. See the examples. 

If you want to draw a circle around an object or selection, use **circleSelection**. If you want pure flexibility over your circle then use **cgoCircle**. 

  
**There are two functions here:**

**cgoCircle**

    

    — creates a CGO circle at some user-specified location

**circleSelection**

    

    —creates a circle around the named object or selection.

  


  * [![Drawn circle.](/images/5/52/Circle1.png)](/index.php/File:Circle1.png "Drawn circle.")

Drawn circle. 

  * [![CGO circle.](/images/8/87/Circle2.png)](/index.php/File:Circle2.png "CGO circle.")

CGO circle. 

  * [![Circles of specified radius.](/images/4/47/Circle3.png)](/index.php/File:Circle3.png "Circles of specified radius.")

Circles of specified radius. 

  * [![Circle with specified width.](/images/d/da/CircleR.png)](/index.php/File:CircleR.png "Circle with specified width.")

Circle with specified width. 

  * [![Circle with a line width of 150. Pores anyone?](/images/d/d3/CircleR2.png)](/index.php/File:CircleR2.png "Circle with a line width of 150. Pores anyone?")

Circle with a line width of 150. Pores anyone? 




# Usage
    
    
    import math
    import pymol
    from pymol.cgo import *
    import random
    
    def cgoCircle(x, y, z, r=8.0, cr=1.0, cg=0.4, cb=0.8, w=2.0):
      """
      Create a CGO circle
    
      PARAMS
            x, y, z
              X, Y and Z coordinates of the origin
    
            r
              Radius of the circle
    
            cr, cg, cb
              Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].
    
            w
              Line width of the circle
    
      RETURNS
            the CGO object (it also loads it into PyMOL, too).
    
      """
      x = float(x)
      y = float(y)
      z = float(z)
      r = abs(float(r))
      cr = abs(float(cr))
      cg = abs(float(cg))
      cb = abs(float(cb))
      w = float(w)
    
      obj = [ BEGIN, LINES, COLOR, cr, cg, cb ]
      for i in range(180):
            obj.append( VERTEX )
            obj.append(r*math.cos(i) + x )
            obj.append(r*math.sin(i) + y )
            obj.append(z)
            obj.append( VERTEX )
            obj.append(r*math.cos(i+0.1) + x )
            obj.append(r*math.sin(i+0.1) + y )
            obj.append(z)
      obj.append(END)
     
      cName = cmd.get_unused_name("circle_")
      cmd.load_cgo( obj, cName )
      cmd.set("cgo_line_width", w, cName )
      return obj
    
    
    def circleSelection( selName, r=None, cr=1.0, cg=0.4, cb=0.8, w=2.0 ):
      """
      circleSelection -- draws a cgo circle around a given selection or object
    
      PARAMS
            selName
              Name of the thing to encircle.
    
            r
              Radius of circle.
              DEFAULT: This cript automatically defines the radius for you.  If
              you select one atom and the resultant circle is too small, then
              you can override the script's calculation of r and specify your own.
    
            cr, cg, cb
              red, green and blue coloring, each a value in the range [0.0, 1.0]
    
      RETURNS
            The circle object.
    
      """
      ((minX, minY, minZ), (maxX, maxY, maxZ)) = cmd.get_extent(selName)
    
      if r==None:
            r = max( [maxX-minX, maxY-minY, maxZ-minZ] )
    
      stored.coords = []
      cmd.iterate_state(1, selName, "stored.coords.append([x,y,z])")
      l = len(stored.coords)
    
      centerX = sum(map(lambda x: x[0], stored.coords)) / l
      centerY = sum(map(lambda x: x[1], stored.coords)) / l
      centerZ = sum(map(lambda x: x[2], stored.coords)) / l
    
      return cgoCircle( centerX, centerY, centerZ, r, cr, cg, cb, w )
    
    
    cmd.extend( "cgoCircle", cgoCircle )
    cmd.extend( "circleSelection", circleSelection )
    

# Updates

  * Line width option
  * better circle naming



Retrieved from "[https://pymolwiki.org/index.php?title=CgoCircle&oldid=10881](https://pymolwiki.org/index.php?title=CgoCircle&oldid=10881)"


---

## Contact Surface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This script calculates individual or global contact areas between a receptor molecule and a (multimodel) bundle of docked ligand structures. The exact contact surface area values (in Angstrom^2) are printed to the screen and also appended to a file called contactareas.txt. If only a single global contact surface is calculated, a selection named "contact" is created that includes all receptor atoms within 3.9A of any ligand atom to illustrate the _approximate_ contact surface. 

  
The parameters are: 

**receptor** _(string)_

    

    The name of the selection/object representing the receptor protein

**ligand** _(string)_

    

    The name of the selection/object representing the ligand
    Note that this may be another protein!

**states** _(integer)_ , default:0 

    

    Calculate contact surface between the receptor and the first n states of the ligand.
    If states = 0, the script calculates a global contact surface which takes all possible ligand states into account.

  


# Usage
    
    
    contact_surface receptor, ligand, [states=0]
    

  


# The Code
    
    
    #contact_surface v.3.0
    #Copyleft Martin Christen, 2013
    
    from pymol import cmd,stored
    def contact_surface(receptor,ligand,states=0):
    
    	"""
    	AUTHOR
    	Martin Christen
    	
    	DESCRIPTION
    	This script calculates individual or global contact surfaces between a
    	receptor molecule and a bundle of docked ligand structures (which have
    	to be loaded into PyMOL as a multimodel object).
    	
    	The exact contact surface area values (in Angstrom^2) are printed to
    	the screen and also appended to a file called contactareas.txt
    	
    	If only a single global contact surface is calculated, a selection
    	named "contact" is created that includes all receptor atoms within
    	3.9A of any ligand atom.
    	
    	USAGE
    	contact_surface receptor, ligand, [states=0]
    	
    	PARAMETERS
    	
    	receptor (string)
    	The name of the selection/object representing the receptor protein
    	
    	ligand (string)
    	The name of the selection/object representing the ligand.
    	Note that this may be another protein!
    	
    	states (integer)
    	Calculate contact surface between the receptor and the first n states
    	of the ligand. If states = 0 (default), the script calculates a global
    	contact surface which takes  all possible ligand states into account.
    	"""
    	# sanity check the number of states
    	states = abs(int(states))
    	
    	# make sure all atoms within an object occlude one another
    	cmd.flag('ignore','none')
    	
    	# use solvent-accessible surface with high sampling density
    	cmd.set('dot_solvent','1')
    	cmd.set('dot_density','3')
    	
    	#if the 'states' parameter = 0 create a superposition of all ligand states
    	if states == 0:
    		cmd.split_states(ligand)
    		cmd.group('ligandtemp',ligand+"_*")
    		cmd.create(ligand+"_all",'ligandtemp')
    		cmd.delete('ligandtemp')
    		
    		#create complex
    		cmd.create('complextemp',ligand+"_all "+receptor)
    	
    		#measure area
    		ligand_area=cmd.get_area(ligand+"_all")
    		receptor_area=cmd.get_area(receptor)
    		complex_area=cmd.get_area('complextemp')
    		#normalize since the area is counted TWICE (once on receptor and once on ligand)
    		contact_area=((ligand_area + receptor_area) - complex_area) / 2
    		#delete complex
    		cmd.delete('complextemp')
    		
    		#create the contact surface
    		cmd.select('contact',"("+receptor+" and ("+ligand+"_all around 3.9))")
    		
    		#print contact surface area
    		f=open('contactareas.txt','a')
    		print "%s - %s : " % (receptor,ligand),
    		print >>f, "%-s\t%-s\t" % (receptor,ligand),
    		print >>f, "%-s" % (contact_area)
    		print contact_area
    		f.close()
    		print "The GLOBAL contact area between "+receptor+ " and "+ligand+" is (A^2):"
    		print ((ligand_area + receptor_area) - complex_area) / 2
    	
    	#If 'states' <> 0 calculate the contact areas to the first 'states' ligand states.
    	#No individual contact surface objects are created to avoid overloading PyMOL.
    	else:
    		#create an object for each ligand state
    		cmd.split_states(ligand)
    		
    		#sanity check: do not exceed that maximum number of states
    		if states > cmd.count_states(ligand):
    			states = cmd.count_states(ligand)
    		
    		#calculate contact surface area
    		print "The contact areas between "+receptor+" and "+ligand+" [states 1 - "+str(states)+"] are (A^2):"
    		#start looping
    		for s in range(1,states+1):
    			#create complex
    			cmd.create("tmp",ligand,s,1)
    			cmd.create('complextemp',"tmp "+receptor)
    			#measure areas
    			ligand_area=cmd.get_area('tmp')
    			receptor_area=cmd.get_area(receptor)
    			complex_area=cmd.get_area('complextemp')
    			#normalize since the area is counted TWICE (once on receptor and once on ligand)
    			contact_area=((ligand_area + receptor_area) - complex_area)/2
    			#delete temporary files
    			cmd.delete('tmp')
    			cmd.delete(ligand+"_*")
    			cmd.delete('complextemp')
    			#print contact surface area
    			f=open('contactareas.txt','a')
    			print "%s - %s_%-5s: " % (receptor,ligand,s),
    			print >>f, "%-s\t%-s_%-5s\t" % (receptor,ligand,s),
    			print >>f, "%-s" % (contact_area)
    			print contact_area
    			f.close()
    
    cmd.extend("contact_surface",contact_surface)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Contact_Surface&oldid=10991](https://pymolwiki.org/index.php?title=Contact_Surface&oldid=10991)"


---

## Cubes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [scripts/cubes.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/cubes.py)  
Author(s)  | [Thomas Holder](/index.php/User:Speleo3 "User:Speleo3")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the **cubes** and **tetrahedra** commands which create CGOs with cube or tetrahedron representations for all atoms in selection. 

## Usage
    
    
    cubes [ selection [, name [, state [, scale [, atomcolors ]]]]]
    
    
    
    tetrahedra [ selection [, name [, state [, scale [, atomcolors ]]]]]
    

## Example
    
    
    fetch 1rx1, async=0
    cubes solvent & b < 50
    tetrahedra solvent & b > 50
    

Retrieved from "[https://pymolwiki.org/index.php?title=Cubes&oldid=13899](https://pymolwiki.org/index.php?title=Cubes&oldid=13899)"


---

## Distancetoatom

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | Python Module   
---|---  
Download  | [scripts/distancetoatom.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/distancetoatom.py)  
Author(s)  | [Andreas Warnecke](/index.php/User:Andwar "User:Andwar") and [Jared Sampson](/index.php/User:Jaredsampson "User:Jaredsampson")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![get distance to atoms within cutoff](/images/c/c2/Distancetoatom_0.png)](/index.php/File:Distancetoatom_0.png "get distance to atoms within cutoff")

## Contents

  * 1 About distancetoatom
  * 2 Usage
    * 2.1 Arguments
  * 3 Examples
  * 4 SEE ALSO



## About distancetoatom

**distancetoatom** prints all distances between a specified atom, coordinate or group selection center and all atoms within cutoff distance that are part of the selection.  
All coordinates and distances can be saved in a csv-style text file report and/or can be stored in a (custom) atom property, if defined. 

  * For instruction on setting up plugin import see [Git intro](/index.php/Git_intro "Git intro") or [Plugin Manager](/index.php/Plugin_Manager "Plugin Manager")



## Usage
    
    
        distancetoatom [ origin [, cutoff [, filename [, selection [, state [, property_name [, coordinates [, decimals [, sort [, quiet ]]]]]]]]]]
    

|   
---|---  
  
### Arguments

Arguments for distancetoatom   
---  
Keyword  | Default  | What it does   
origin  | pk1  | "origin" defines the coordinates for the origin and can be: 1\. a list with coordinates [x,y,z]  
2\. a single atom selection string  
3\. a multi-atom selection string (center will be used)   
cutoff  | 10  | "cutoff" defines the maximum distance, atoms further away are not considered   
filename  | None  | "filename" is the name of a report file that will be created (optional). set to e.g. 'report.txt' to create a report. This file is CSV style and can be imported into EXCEL.  
(omit or set to "", None, 0 or False to disable)   
selection  | all  | only atoms within "selection" (and cutoff distance) will be used for calculation   
state  | 0  | "state" is the state that will be used for calculations use 0 (omit), to automatically use the current state   
property_name  | p.dist  | "property_name" defines a (custom) property in which the calculated distance will be stored. This can be used e.g. for labeling or spectrum coloring etc. (cf. examples)  
NB! older PyMOL versions only support b or q and the distance will overwrite this property if set.   
coordinates  | 0  | "coordinates" toggles whether, besides distance, the atom coordinates will be included in the report.   
decimals  | 3  | "decimals" is the number of decimal places calculated. Note that PDB files do not support a higher resolution than 3.   
sort  | 1  | "sort" defines how the output will be sorted: 1: ascending (default)  
0: no sorting (by selection)  
-1: descending   
quiet  | 1  | toggle verbosity   
  
  


## Examples
    
    
    # make function available to PyMOL
    import distancetoatom
    
    ##############################
    # Example 1: print all distances
    frag LYS
    edit name CA
    distancetoatom
    
    ##############################
    # Example 2: print all distances to file
    fetch 1cll, async=0
    distancetoatom origin=/1cll//A/CA`150/CA, filename=distances.txt, selection=elem O, property_name=b
    # the file is CSV-format and can be imported, e.g. to EXCEL
    
    # format
    hide everything
    select byres (resi 150 expand 5 and not resn hoh) 
    show_as sticks, sele
    util.cbaw sele
    show_as spheres, resi 150
    zoom sele
    
    ##########
    # Label by stored distance (property_name)
    label sele and elem O, "%.2f"%b
    
    
    ##############################
    # Example 3: color an object by distance to a coordinate
    
    fetch 1cll, async=0
    distancetoatom [0,0,0], cutoff=1000, property_name=b
    # the distance is stored in the b-factor in this case
    # newer PyMOL versions support custom properties, e.g. "p.dist"
    spectrum b, rainbow_rev
    

|  [![Example 2](/images/5/5c/Distancetoatom_2.png)](/index.php/File:Distancetoatom_2.png "Example 2")  
  
[![Example 3](/images/3/3a/Distancetoatom_1.png)](/index.php/File:Distancetoatom_1.png "Example 3")  
---|---  
  
  


## SEE ALSO

  * [Distance](/index.php/Distance "Distance")
  * [Get Distance](/index.php/Get_Distance "Get Distance")
  * [Get raw distances](/index.php/Get_raw_distances "Get raw distances")



Retrieved from "[https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939](https://pymolwiki.org/index.php?title=Distancetoatom&oldid=13939)"


---

## Draw Protein Dimensions

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/Draw_Protein_Dimensions.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/Draw_Protein_Dimensions.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw the dimensions of a protein based on an Inertia Axis Aligned Bounding Box (IABB). 

  
The idea behind this script is to calculate an approximate minimal bounding box (MBB) to extract the cell dimensions of a protein. To calculate the MBB is not trivial and usually the Axis Aligned Bounding Box (AABB) does not show up the real dimensions of the protein. This script calculates the inertia tensor of the object, extract the eigenvectors and use them to rotate the molecule (using as rotation matrix the transpose of the eigenvector matrix). The result is a molecule oriented with the inertia axis aligned with the cartesian axis. A new Bounding Box is calculated, which is called Inertia Axis Aligned Bounding Box (IABB), with a volume always lower than the AABB volume, and in many cases may correspond with the MBB. 

As always with these type of things, you have to use at your own risk. I did not try all the possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado (at) gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run path-to-the-script/Draw_Protein_Dimensions.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_Protein_Dimensions _selection_**

This will draw the cell dimensions of your selection based on a IABB box. It also generates the IABB box and the inertia axis, you just need to do "show cgo" to display them. 

You could also try: 

**draw_BB _selection_**

This will draw the AABB and IABB boxes with their cell dimensions and show up their volumes, you can compare them. 

  


# Examples
    
    
    # download the source and save as Draw_Protein_Dimensions.py
    run Draw_Protein_Dimensions.py
    fetch 2vak
    # calculate the dimensions of the full ASU
    draw_Protein_Dimensions 2vak
    

  * [![Dimensions of 2vak ASU](/images/7/77/2vak_ASU.png)](/index.php/File:2vak_ASU.png "Dimensions of 2vak ASU")

Dimensions of 2vak ASU 



    
    
    # you can extract only one chain and calculates the dimensions
    draw_Protein_Dimensions obj01
    

  * [![Dimensions of one protomer \(chain A - 2vak\)](/images/e/ed/2vak_A.png)](/index.php/File:2vak_A.png "Dimensions of one protomer \(chain A - 2vak\)")

Dimensions of one protomer (chain A - 2vak) 



    
    
    # You can also draw the Bounding boxes (AABB and IABB) to compare them.
    fetch 2vak
    draw_BB 2vak
    

  


  * [![Axis-aligned bounding box](/images/4/46/2vak_AABB.png)](/index.php/File:2vak_AABB.png "Axis-aligned bounding box")

Axis-aligned bounding box 

  * [![Inertia-axis-aligned bounding box](/images/4/41/2vak_IABB.png)](/index.php/File:2vak_IABB.png "Inertia-axis-aligned bounding box")

Inertia-axis-aligned bounding box 




Retrieved from "[https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902](https://pymolwiki.org/index.php?title=Draw_Protein_Dimensions&oldid=13902)"


---

## DrawBoundingBox

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Example
  * 3 Installation
  * 4 See Also



# Overview

Draws a bounding box around a given selection. 

  * [![Example of a bounding box](/images/5/51/DrawMinBB.png)](/index.php/File:DrawMinBB.png "Example of a bounding box")

Example of a bounding box 

  * [![Two bounding boxes, one with 0 padding, the other with 10 Ang padding.](/images/6/63/Bb_with_padding.png)](/index.php/File:Bb_with_padding.png "Two bounding boxes, one with 0 padding, the other with 10 Ang padding.")

Two bounding boxes, one with 0 padding, the other with 10 Ang padding. 




# Example
    
    
    run ~/drawBoundingBox.py
    fetch 1jsd
    drawBoundingBox 1jsd, r=0.33, g=0.80
    
    # example from above w/padding, draw it light blue
    drawBoundingBox padding=10, r=0.5, g=0.8, b=1.0
    

# Installation

Copy the source code to your computer, and execute it by issuing the command "run /path/to/drawBoundingBox.py" in PyMOL. 
    
    
    # -*- coding: utf-8 -*-                                                                                     
    from pymol.cgo import *                                                                                     
    from pymol import cmd                                                                                       
    from random import randint                                                                                  
    
    #############################################################################
    #                                                                            
    # drawBoundingBox.py -- Draws a box surrounding a selection 
    #
    #                                                                            
    # AUTHOR: Jason Vertrees                                                   
    # DATE  : 2/20/2009                                                          
    # NOTES : See comments below.                                                
    #                                                                            
    #############################################################################
    def drawBoundingBox(selection="(all)", padding=0.0, linewidth=2.0, r=1.0, g=1.0, b=1.0):     
            """                                                                  
            DESCRIPTION                                                          
                    Given selection, draw the bounding box around it.          
    
            USAGE:
                    drawBoundingBox [selection, [padding, [linewidth, [r, [g, b]]]]]
    
            PARAMETERS:
                    selection,              the selection to enboxen.  :-)
                                            defaults to (all)
       
                    padding,                defaults to 0
    
                    linewidth,              width of box lines
                                            defaults to 2.0
    
                    r,                      red color component, valid range is [0.0, 1.0]
                                            defaults to 1.0                               
    
                    g,                      green color component, valid range is [0.0, 1.0]
                                            defaults to 1.0                                 
    
                    b,                      blue color component, valid range is [0.0, 1.0]
                                            defaults to 1.0                                
    
            RETURNS
                    string, the name of the CGO box
    
            NOTES
                    * This function creates a randomly named CGO box that minimally spans the protein. The
                    user can specify the width of the lines, the padding and also the color.                            
            """                                                                                                    
    
            ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)
    
            print "Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)
    
            minX = minX - float(padding)
            minY = minY - float(padding)
            minZ = minZ - float(padding)
            maxX = maxX + float(padding)
            maxY = maxY + float(padding)
            maxZ = maxZ + float(padding)
    
            if padding != 0:
                     print "Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)
    
            boundingBox = [
                    LINEWIDTH, float(linewidth),
    
                    BEGIN, LINES,
                    COLOR, float(r), float(g), float(b),
    
                    VERTEX, minX, minY, minZ,       #1
                    VERTEX, minX, minY, maxZ,       #2
    
                    VERTEX, minX, maxY, minZ,       #3
                    VERTEX, minX, maxY, maxZ,       #4
    
                    VERTEX, maxX, minY, minZ,       #5
                    VERTEX, maxX, minY, maxZ,       #6
    
                    VERTEX, maxX, maxY, minZ,       #7
                    VERTEX, maxX, maxY, maxZ,       #8
    
    
                    VERTEX, minX, minY, minZ,       #1
                    VERTEX, maxX, minY, minZ,       #5
    
                    VERTEX, minX, maxY, minZ,       #3
                    VERTEX, maxX, maxY, minZ,       #7
    
                    VERTEX, minX, maxY, maxZ,       #4
                    VERTEX, maxX, maxY, maxZ,       #8
    
                    VERTEX, minX, minY, maxZ,       #2
                    VERTEX, maxX, minY, maxZ,       #6
    
    
                    VERTEX, minX, minY, minZ,       #1
                    VERTEX, minX, maxY, minZ,       #3
    
                    VERTEX, maxX, minY, minZ,       #5
                    VERTEX, maxX, maxY, minZ,       #7
    
                    VERTEX, minX, minY, maxZ,       #2
                    VERTEX, minX, maxY, maxZ,       #4
    
                    VERTEX, maxX, minY, maxZ,       #6
                    VERTEX, maxX, maxY, maxZ,       #8
    
                    END
            ]
    
            boxName = "box_" + str(randint(0,10000))
            while boxName in cmd.get_names():
                    boxName = "box_" + str(randint(0,10000))
    
            cmd.load_cgo(boundingBox,boxName)
            return boxName
    
    cmd.extend ("drawBoundingBox", drawBoundingBox)
    

# See Also

[Bounding_Box](/index.php/Bounding_Box "Bounding Box")

Retrieved from "[https://pymolwiki.org/index.php?title=DrawBoundingBox&oldid=10935](https://pymolwiki.org/index.php?title=DrawBoundingBox&oldid=10935)"


---

## Dump2CGO

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

Dumps a PyMOL object to a CGO object. 

# The Code
    
    
    from pymol import cmd
    from pymol.cgo import *
    
    def dump2surfaceCGO():
        CGOobj = []
        dumpedFile = open("dump.tmp").read()
        for block in dumpedFile.split('\n\n'):
            CGOobj.append(BEGIN)
            CGOobj.append(TRIANGLES)
    
            for line in block.split('\n'):
                if line == '':
                    continue
    
                vals = line.split()
                CGOobj.append(NORMAL)
                CGOobj.append(float(vals[3]))
                CGOobj.append(float(vals[4]))
                CGOobj.append(float(vals[5]))
                CGOobj.append(VERTEX)
                CGOobj.append(float(vals[0]))
                CGOobj.append(float(vals[1]))
                CGOobj.append(float(vals[2]))
    
            CGOobj.append(END)
        return CGOobj
    
    def dump2meshCGO():
        CGOobj = []
        dumpedFile = open("dump.tmp").read()
        for block in dumpedFile.split('\n\n'):
            CGOobj.append(BEGIN)
            CGOobj.append(LINE_STRIP)
    
            for line in block.split('\n'):
                if line == '':
                    continue
    
                CGOobj.append(VERTEX)
                vals = line.split()
    
                CGOobj.append(float(vals[0]))
                CGOobj.append(float(vals[1]))
                CGOobj.append(float(vals[2]))
    
            CGOobj.append(END)
        return CGOobj
    
    def getType(objname):
        session = cmd.get_session()['names']
        for obj in session:
            if obj == None:
                continue
            if obj[0] != objname:
                continue
            return obj[4]
        return -1
    
    
    def dump2CGO(obj):
        cmd.dump("dump.tmp", obj)
        type = getType(obj)
        cgo = []
        if (type == 3): # Mesh
            cgo = dump2meshCGO()
        elif (type == 7): #Surface
            cgo = dump2surfaceCGO()
        else:
            print("Unknown type")
            return
    
        cmd.load_cgo(cgo, "CGO " + obj)
    
    cmd.extend('dump2CGO', dump2CGO)
    cmd.auto_arg[0]['dump2CGO'] = [cmd.object_sc, 'object', '']
    

Retrieved from "[https://pymolwiki.org/index.php?title=Dump2CGO&oldid=12733](https://pymolwiki.org/index.php?title=Dump2CGO&oldid=12733)"


---

## Ellipsoid

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This script provides methods that create [cgos](/index.php?title=Cgo&action=edit&redlink=1 "Cgo \(page does not exist\)") as triangles. It uses code that is ported from [this c++ code](http://www.gamedev.net/reference/articles/article1172.asp) and seems to be correct! 

Here is the script. The last four lines show this in use, by making an ellipse and a toroid and loading them into pymol. This is done most easily by something like "cmd.load_cgo(makeEllipsoid(1, 1, 1, 2, 3, 4), 'ellipsoid')" which makes an ellipsoid at x, y, z = 1, 1, 1 and dimensions 2, 3, 4 and called 'ellipsoid'. 
    
    
    from pymol.cgo import BEGIN, COLOR, TRIANGLES, VERTEX, NORMAL, END
    from pymol import cmd
    
    def signOfFloat(f):
            if f < 0: return -1
            if f > 0: return 1
            return 0
    
    def sqC(v, n):
            return signOfFloat(math.cos(v)) *  math.pow(math.fabs(math.cos(v)), n)
    
    def sqCT(v, n, alpha):
            return alpha + sqC(v, n)
    
    def sqS(v, n):
            return signOfFloat(math.sin(v)) * math.pow(math.fabs(math.sin(v)), n)
    
    def sqEllipsoid(x, y, z, a1, a2, a3, u, v, n, e):
            x = a1 * sqC(u, n) * sqC(v, e) + x
            y = a2 * sqC(u, n) * sqS(v, e) + y
            z = a3 * sqS(u, n) + z
            nx = sqC(u, 2 - n) * sqC(v, 2 - e) / a1
            ny = sqC(u, 2 - n) * sqS(v, 2 - e) / a2
            nz = sqS(u, 2 - n) / a3
            return x, y, z, nx, ny, nz
    
    def sqToroid(x, y, z, a1, a2, a3, u, v, n, e, alpha):
            a1prime = 1.0 / (a1 + alpha)
            a2prime = 1.0 / (a2 + alpha)
            a3prime = 1.0 / (a3 + alpha)
            x = a1prime * sqCT(u, e, alpha) * sqC(v, n)
            y = a2prime * sqCT(u, e, alpha) * sqS(v, n)
            z = a3prime * sqS(u, e)
            nx = sqC(u, 2 - e) * sqC(v, 2 - n) / a1prime
            ny = sqC(u, 2 - e) * sqS(v, 2 - n) / a2prime
            nz = sqS(u, 2 - e) / a3prime
            return x, y, z, nx, ny, nz
    
    def makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, n, e, u1, u2, v1, v2, u_segs, v_segs, color=[0.5, 0.5, 0.5]):
    
            r, g, b = color
    
            # Calculate delta variables */
            dU = (u2 - u1) / u_segs
            dV = (v2 - v1) / v_segs
    
            o = [ BEGIN, TRIANGLES ]
    
            U = u1
            for Y in range(0, u_segs):
                    # Initialize variables for loop */
                    V = v1
                    for X in range(0, v_segs):
                            # VERTEX #1 */
                            x1, y1, z1, n1x, n1y, n1z = sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, e)
                            x2, y2, z2, n2x, n2y, n2z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, e)
                            x3, y3, z3, n3x, n3y, n3z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e)
                            x4, y4, z4, n4x, n4y, n4z = sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, e)
    
                            o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
    
                            # Update variables for next loop */
                            V += dV
                    # Update variables for next loop */
                    U += dU
            o.append(END)
            return o
    
    def makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, n, e, u1, u2, v1, v2, u_segs, v_segs, color=[0.5, 0.5, 0.5]):
    
            r, g, b = color
    
            # Calculate delta variables */
            dU = (u2 - u1) / u_segs
            dV = (v2 - v1) / v_segs
    
            o = [ BEGIN, TRIANGLES ]
    
            U = u1
            for Y in range(0, u_segs):
                    # Initialize variables for loop */
                    V = v1
                    for X in range(0, v_segs):
                            # VERTEX #1 */
                            x1, y1, z1, n1x, n1y, n1z = sqToroid(x, y, z, a1, a2, a3, U, V, n, e, alpha)
                            x2, y2, z2, n2x, n2y, n2z = sqToroid(x, y, z, a1, a2, a3, U + dU, V, n, e, alpha)
                            x3, y3, z3, n3x, n3y, n3z = sqToroid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e, alpha)
                            x4, y4, z4, n4x, n4y, n4z = sqToroid(x, y, z, a1, a2, a3, U, V + dV, n, e, alpha)
    
                            o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                            o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                            o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                            o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
    
                            # Update variables for next loop */
                            V += dV
                    # Update variables for next loop */
                    U += dU
            o.append(END)
            return o
    
    def makeEllipsoid(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeCylinder(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeSpindle(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeDoublePyramid(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 2.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makePillow(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 0.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)
    
    def makeRoundCube(x, y, z, a1, a2, a3):
                    return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.2, 0.2, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
    
    def makeToroid(x, y, z, a1, a2, a3, alpha):
                    return makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, 1.0, 1.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)
    
    x, y, z, rx, ry, rz = 1, 1, 1, 1, 2, 3
    cmd.load_cgo(makeEllipsoid(x, y, z, rx, ry, rz), 'ellipsoid-cgo')
    x, y, z, rx, ry, rz = 1, 1, 1, 8, 2, 2
    cmd.load_cgo(makeToroid(x, y, z, rx, ry, rz, 3), 'toroid-cgo')
    

Retrieved from "[https://pymolwiki.org/index.php?title=Ellipsoid&oldid=6396](https://pymolwiki.org/index.php?title=Ellipsoid&oldid=6396)"


---

## Mark center

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Module](/index.php/Running_Scripts#Python_Modules "Running Scripts")  
---|---  
Download  | [mark_center.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/mark_center.py)  
Author(s)  | [User:Inchoate](/index.php/User:Inchoate "User:Inchoate")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
Mark_center is a small callback that puts a cross at the center of space, where the camera is pointing. 

You could also use this to setup arbitrary CGOs (points, lines, etc) in 3D space. 
    
    
    from pymol import cmd
    
    def crosshair_put_center():
       t = cmd.get_position()
       m = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0] + t + [1]
       cmd.set_object_ttt('crosshair', m)
    
    cmd.load_callback(crosshair_put_center, '_crosshair_cb')
    
    cmd.pseudoatom('crosshair', pos=(0,0,0))
    cmd.show_as('nonbonded', 'crosshair')
    

Retrieved from "[https://pymolwiki.org/index.php?title=Mark_center&oldid=12541](https://pymolwiki.org/index.php?title=Mark_center&oldid=12541)"


---

## Modevectors

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/modevectors.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/modevectors.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
### DESCRIPTION

[![](/images/3/34/Modevectors.png)](/index.php/File:Modevectors.png)

[](/index.php/File:Modevectors.png "Enlarge")

Example showing modevectors in action. (See Examples below).

**modevectors.py** is a PyMol script that was originally written to visualize results obtained from Normal Mode Analysis (NMA) by drawing arrows or vectors from a starting structure to a final structure. However, this script can be used to visualize the direction of motion between two specified states (e.g. targeted MD, comparing open and closed structures, etc). The strength of this script is that the arrows are highly customizable with numerous options to choose from (see script below). It is important to note that the two states MUST be identical except for the X, Y, and Z coordinates. That is, you cannot draw vectors between homologous structures. The default settings sets the background to white and displays the arrows along with the first object frame (in cartoon representation). 

  
Note: Ray tracing these CGO arrows appears to be a memory intensive process so please save a session before ray tracing and use "pymol -qc" where necessary. If anybody can come up with a method to improve this please e-mail me (see address in script). 

Update: A new way of drawing the arrows has been implemented for PyMOL v.1.1 and has been updated in this code (which automatically detects the version that you are using). The new method uses the built in cone primitive rather than concentric circles decreasing in size. Thus, the new method is much less memory intensive and is a significant improvement over the old method. Additionally, you may want to turn off [ray shadow](/index.php/Ray_shadow "Ray shadow") depending on the [scene](/index.php/Scene "Scene") that you are viewing as shadows can be confusing to others when seen in publications. 

### USAGE

load the script using the [run](/index.php/Run "Run") command 
    
    
    modevectors first_obj_frame, last_obj_frame [,first_state=1 [,last_state=1 [,outname=modevectors \
      [,head=1.0 [,tail=0.3 \[,head_length=1.5 [,headrgb=(1.0,1.0,1.0) [,tailrgb=(1.0,1.0,1.0) [,cutoff=4.0 
      [,skip=0 [,cut=0.5 [,atom=CA [,stat=show [,factor=1.0 [,notail=0]]]]]]]]]]]]]]
    

Please see the script comments for further custom options. Once the script completes, it will generate a new object called "modevectors" (which can be changed through the options). 

  


### EXAMPLES
    
    
    modevectors 1E3M, 1W7A
    modevectors 1E3M, 1W7A, outname="arrows"
    modevectors 1E3M, 1W7A, atom="P"
    

Copy/paste the following code to see an example of modevectors. This uses a multistate protein and the arrows are connected between the first and last states. 
    
    
    import modevectors
    # fetch the PDBs from pdb.org
    fetch 1c3y, finish=1, multiplex=0, async=0
    # separate the first and last states of the NMR ensemble to individual objects
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
    # run the modevectors code
    modevectors 1c3y_0001, 1c3y_0023
    # just setup a nice representation
    as cartoon, 1c3y_0001 or 1c3y_0023
    show cgo, modevectors
    color marine, 1c3y_0001
    color purpleblue, 1c3y_0023
    

The following set of examples will illustrate the power of each optional argument. Each example should be compared to the default figure in the table below. 

[![](/images/3/35/Mv_default.png)](/index.php/File:Mv_default.png) [](/index.php/File:Mv_default.png "Enlarge")Fig.1 - Default Settings | [![](/images/b/b0/Mv_fig2.png)](/index.php/File:Mv_fig2.png) [](/index.php/File:Mv_fig2.png "Enlarge")Fig.2 - Arrow direction drawn from CA to CA | [![](/images/4/4f/Mv_fig3.png)](/index.php/File:Mv_fig3.png) [](/index.php/File:Mv_fig3.png "Enlarge")Fig.3 - Arrow head radius  
---|---|---  
[![](/images/4/4f/Mv_fig4.png)](/index.php/File:Mv_fig4.png) [](/index.php/File:Mv_fig4.png "Enlarge")Fig.4 - Arrow tail radius | [![](/images/f/f1/Mv_fig5.png)](/index.php/File:Mv_fig5.png) [](/index.php/File:Mv_fig5.png "Enlarge")Fig.5 - Arrow head length | [![](/images/e/e4/Mv_fig6.png)](/index.php/File:Mv_fig6.png) [](/index.php/File:Mv_fig6.png "Enlarge")Fig.6 - Arrow head RGB  
[![](/images/f/fb/Mv_fig7.png)](/index.php/File:Mv_fig7.png) [](/index.php/File:Mv_fig7.png "Enlarge")Fig.7 - Arrow tail RGB | [![](/images/b/bc/Mv_fig8.png)](/index.php/File:Mv_fig8.png) [](/index.php/File:Mv_fig8.png "Enlarge")Fig.8 - Arrow length cutoff | [![](/images/7/7f/Mv_fig9.png)](/index.php/File:Mv_fig9.png) [](/index.php/File:Mv_fig9.png "Enlarge")Fig.9 - Skip every other arrow  
[![](/images/d/d7/Mv_fig10.png)](/index.php/File:Mv_fig10.png) [](/index.php/File:Mv_fig10.png "Enlarge")Fig.10 - Subtract value from vector length | [![](/images/4/42/Mv_fig11.png)](/index.php/File:Mv_fig11.png) [](/index.php/File:Mv_fig11.png "Enlarge")Fig.11 - Draw arrows from CB atoms | [![](/images/b/b1/Mv_fig12.png)](/index.php/File:Mv_fig12.png) [](/index.php/File:Mv_fig12.png "Enlarge")Fig.12 - Shorten by 50%  
| [![](/images/a/a3/Mv_fig13.png)](/index.php/File:Mv_fig13.png) [](/index.php/File:Mv_fig13.png "Enlarge")Fig.13 - Multiple options being used (see final example below)  
      
    
    reinitialize
    import modevectors
    
    fetch 1c3y, async=0
    
    split_states 1c3y, 1, 1
    split_states 1c3y, 23, 23
    hide
     
    #This is the default setting (Fig.1)
    modevectors 1c3y_0001, 1c3y_0023
     
    #This is shows that the direction of the arrows is drawn from the 1c3y_0001 towards 1c3y_0023 (Fig.2)
    show cartoon, 1c3y_0023
    color red, 1c3y_0023
     
    #This controls the base radius of the cone/arrow head in angstroms (Fig.3)
    modevectors 1c3y_0001, 1c3y_0023, head=2.5
     
    #This controls the radius of the cylinders/arrow tails in angstroms (Fig.4)
    modevectors 1c3y_0001, 1c3y_0023, tail=0.75
     
    #This controls the length of the cone/arrow head in angstroms (Fig.5)
    #Note that this option does NOT increase the vector length but simply changes the tail length
    modevectors 1c3y_0001, 1c3y_0023, head_length=3.0
     
    #This controls the colour of the cone/arrow head using RGB values (Fig.6)
    modevectors 1c3y_0001, 1c3y_0023, headrgb=(1.0,0.0,0.0)
     
    #This controls the colour of the cylinder/arrow tails using RGB values (Fig.7)
    modevectors 1c3y_0001, 1c3y_0023, tailrgb=(1.0,0.0,0.0)
     
    #This controls the which vectors to show based on a specific cutoff value in angstroms.  Vector lengths that are less 
    #than the cutoff value will not be displayed (Fig.8)
    modevectors 1c3y_0001, 1c3y_0023, cutoff=30.0
     
    #This controls how many vectors to skip (integer value) and is useful when there are too many vectors showing.  
    #Skip=1 will show every other arrow (Fig.9)
    modevectors 1c3y_0001, 1c3y_0023, skip=1
     
    #This controls how much to cut from each vector (in angstroms).  Note that arrows will point in the opposite directions 
    #when too much is cutoff (resulting in a negative vector length) (Fig.10) and should be used with caution!
    modevectors 1c3y_0001, 1c3y_0023, cut=15.0
     
    #This controls which atom to draw a vector from (Fig.11).  Note that this is case-sensitive and is really only useful 
    #when atom=CA or when atom=P (for DNA backbones)
    modevectors 1c3y_0001, 1c3y_0023, atom=CB
     
    #This controls how much to multiple the length of each vector by (percentage increase/decrease) (Fig.12)
    #This example halves the length of each vector (50%)
    modevectors 1c3y_0001, 1c3y_0023, factor=0.5
     
    #This hides the statistics which count atoms skipped, atoms counted (number of arrows showing), and number of atoms 
    #that did not meet the cutoff and are not shown
    modevectors 1c3y_0001, 1c3y_0023, stat=hide
     
    #This example shows multiple options being used (Fig.13)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, tail=1.0, head_length=2.0, headrgb=(1.0,0.0,0.0), tailrgb=(0.0,0.0,1.0),
    cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8
     
    #Finally, this example hides all arrow tails and only uses arrow heads via the notail option(No Figure)
    modevectors 1c3y_0001, 1c3y_0023, head=2.0, cutoff=0.0,skip=0,cut=0,atom=CA,factor=0.8, notail=1
    

Retrieved from "[https://pymolwiki.org/index.php?title=Modevectors&oldid=13914](https://pymolwiki.org/index.php?title=Modevectors&oldid=13914)"


---

## Perp maker

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/perp_maker.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/perp_maker.py)  
Author(s)  | [Jason Vertrees](/index.php/User:Inchoate "User:Inchoate")  
License  | [MIT](http://opensource.org/licenses/MIT)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
The **perp_maker.py** script creates perpendicular planes. 

Nothing to do with cops. Given a simple PyMol scene, attempts to create a CGO background triangle perpendicular to the vector created - which is parallel to the line segment drawn through the camera point and current center of mass - as obtained by "get_position," or "get_view." 

## Usage

  * Load your scene
  * Orient the scene as you wish
  * [Run](/index.php/Running_Scripts "Running Scripts") the script



Could it be any simpler? 

You can rotate and move the plane using the editing features (A > drag, Shift+Drag with right, middle or left mouse button). 

Retrieved from "[https://pymolwiki.org/index.php?title=Perp_maker&oldid=13916](https://pymolwiki.org/index.php?title=Perp_maker&oldid=13916)"


---

## Plane Wizard

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Introduction
  * 2 Original
  * 3 Modification
    * 3.1 Gallery
    * 3.2 plane.py
    * 3.3 Examples



## Introduction

This wizard has a simple purpose - to draw a cgo plane that passes through three points picked by the user. Most of the wizard itself was copied from the measure wizard. 

To use, put it in the same directory as the other wizards. This is not quality code, and there may be bugs, but it seems to work okay. 

## Original
    
    
    import pymol
    from pymol import cmd
    from pymol.wizard import Wizard
    from chempy import cpv
    from pymol.cgo import *
    
    def makePrimitive(cgo, name):
        az = cmd.get('auto_zoom', quiet=1)
        cmd.set('auto_zoom', 0, quiet=1)
        cmd.load_cgo(cgo, name)
        cmd.set('auto_zoom', az, quiet=1)
    
    def point(p):
        x, y, z = p
        return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]
    
    def line(p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]
    
    def plane(corner1, corner2, corner3, corner4, normal):
        planeObj = []
        planeObj.extend(point(corner1))
        planeObj.extend(point(corner2))
        planeObj.extend(point(corner3))
        planeObj.extend(point(corner4))
        planeObj.extend(line(corner1, corner2))
        planeObj.extend(line(corner2, corner3))
        planeObj.extend(line(corner3, corner4))
        planeObj.extend(line(corner4, corner1))
    
        planeObj.extend([COLOR, 0.8, 0.8, 0.8])
        planeObj.extend([BEGIN, TRIANGLE_STRIP])
        planeObj.append(NORMAL)
        planeObj.extend(normal)
        for corner in [corner1, corner2, corner3, corner4, corner1]:
            planeObj.append(VERTEX)
            planeObj.extend(corner)
        planeObj.append(END)
        return planeObj
    
    def planeFromPoints(point1, point2, point3, facetSize):
        v1 = cpv.normalize(cpv.sub(point2, point1))
        v2 = cpv.normalize(cpv.sub(point3, point1))
        normal = cpv.cross_product(v1, v2)
        v2 = cpv.cross_product(normal, v1)
        x = cpv.scale(v1, facetSize)
        y = cpv.scale(v2, facetSize)
        center = point2
        corner1 = cpv.add(cpv.add(center, x), y)
        corner2 = cpv.sub(cpv.add(center, x), y)
        corner3 = cpv.sub(cpv.sub(center, x), y)
        corner4 = cpv.add(cpv.sub(center, x), y)
        return plane(corner1, corner2, corner3, corner4, normal)
    
    
    class PlaneWizard(Wizard):
    
        def __init__(self):
            Wizard.__init__(self)
    
            # some attributes to do with picking
            self.pick_count = 0
            self.object_count = 0
            self.object_prefix = "pw"
    
            # the plane facet size (the 'radius' of the section of plane we show)
            self.facetSize = 5
    
            self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
            cmd.set("mouse_selection_mode",0) # set selection mode to atomic
            cmd.deselect()
    
        def reset(self):
            cmd.delete(self.object_prefix + "*")
            cmd.delete("sele*")
            cmd.delete("_indicate*")
            cmd.unpick()
            cmd.refresh_wizard()
    
        def delete_all(self):
            cmd.delete("plane*")
    
        def cleanup(self):
            cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
            self.reset()
            self.delete_all()
    
        def get_prompt(self):
            self.prompt = None
            if self.pick_count == 0:
                self.prompt = [ 'Please click on the first atom...']
            elif self.pick_count == 1:
                self.prompt = [ 'Please click on the second atom...' ]
            elif self.pick_count == 2:
                self.prompt = [ 'Please click on the third atom...' ]
            return self.prompt
    
        def do_select(self, name):
            # "edit" only this atom, and not others with the object prefix
            try:
                cmd.edit("%s and not %s*" % (name, self.object_prefix))
                self.do_pick(0)
            except pymol.CmdException, pmce:
                print pmce
    
        def pickNextAtom(self, atom_name):
            # transfer the click selection to a named selection
            cmd.select(atom_name, "(pk1)")
    
            # delete the click selection
            cmd.unpick()
    
            # using the magic of indicate, highlight stuff
            indicate_selection = "_indicate" + self.object_prefix
            cmd.select(indicate_selection, atom_name)
            cmd.enable(indicate_selection)
    
            self.pick_count += 1
            self.error = None
    
            # necessary to force update of the prompt
            cmd.refresh_wizard()
    
        def do_pick(self, picked_bond):
    
            # this shouldn't actually happen if going through the "do_select"
            if picked_bond:
                self.error = "Error: please select bonds, not atoms"
                print self.error
                return
    
            atom_name = self.object_prefix + str(self.pick_count)
            if self.pick_count < 2:
                self.pickNextAtom(atom_name)
            else:
                self.pickNextAtom(atom_name)
    
                point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
                point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
                point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
                plane = planeFromPoints(point1, point2, point3, self.facetSize)
    
                planeName = "plane-%02d" % self.object_count
                self.object_count += 1
                makePrimitive(plane, planeName)
                cmd.show("cgo", "plane*")
    
                self.pick_count = 0
                self.reset()
    
        def get_panel(self):
            return [
                [ 1, 'Plane Wizard',''],
                [ 2, 'Reset','cmd.get_wizard().reset()'],
                [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
                [ 2, 'Done','cmd.set_wizard()'],
            ]
    
    # create an instance
    
    wiz = PlaneWizard()
    
    # make this the active wizard
    
    cmd.set_wizard(wiz)
    

## Modification

Small modifications to the same code as above. 

### Gallery

  * [![Drawing a plane in a protein](/images/2/28/Plane_img1.png)](/index.php/File:Plane_img1.png "Drawing a plane in a protein")

Drawing a plane in a protein 

  * [![Drawing 6 planes in a box](/images/d/d3/Plane_img2.png)](/index.php/File:Plane_img2.png "Drawing 6 planes in a box")

Drawing 6 planes in a box 

  * [![Drawing 6 planes in a box](/images/2/25/Plane_img3.png)](/index.php/File:Plane_img3.png "Drawing 6 planes in a box")

Drawing 6 planes in a box 




### plane.py

Make a **plane.py** file in the same directory where you are working 
    
    
    '''
    Described at PyMOL wiki:
    https://pymolwiki.org/index.php/Plane_Wizard
    
    Authors : Troels Schwarz-Linnet
    Date    : Dec 2016
    Modified: From previous contributors. 
    '''
    
    import pymol
    from pymol import cmd
    from pymol.wizard import Wizard
    from chempy import cpv
    from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, TRIANGLE_STRIP, NORMAL, VERTEX, END, ALPHA
    
    def makePrimitive(cgo, name):
        az = cmd.get('auto_zoom', quiet=1)
        cmd.set('auto_zoom', 0, quiet=1)
        cmd.load_cgo(cgo, name)
        cmd.set('auto_zoom', az, quiet=1)
    
    def point(p):
        x, y, z = p
        return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]
    
    
    def line(p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]
    
    
    def plane(corner1, corner2, corner3, corner4, normal, settings):
        planeObj = []
        planeObj.extend(point(corner1))
        planeObj.extend(point(corner2))
        planeObj.extend(point(corner3))
        planeObj.extend(point(corner4))
        planeObj.extend(line(corner1, corner2))
        planeObj.extend(line(corner2, corner3))
        planeObj.extend(line(corner3, corner4))
        planeObj.extend(line(corner4, corner1))
    
        # Make settings
        if 'ALPHA' in settings:
            planeObj.extend([ALPHA, settings['ALPHA']])
    
        if 'COLOR' in settings:
            planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
        else:
            planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish
    
        planeObj.extend([BEGIN, TRIANGLE_STRIP])
        planeObj.append(NORMAL)
    
        if 'INVERT' in settings:
            if settings['INVERT']==True:
                planeObj.extend(cpv.negate(normal))
            else:
                planeObj.extend(normal)
        else:
            planeObj.extend(normal)
    
    
        for corner in [corner1, corner2, corner3, corner4, corner1]:
            planeObj.append(VERTEX)
            planeObj.extend(corner)
        planeObj.append(END)
        return planeObj
    
    
    def planeFromPoints(p1, p2, p3, vm1=None, vm2=None, center=True, settings={}):
        v1 = cpv.sub(p1, p2)
        v2 = cpv.sub(p3, p2)
        normal = cpv.cross_product(v1, v2)
    
        if 'translate' in settings:
            vtran = cpv.scale(cpv.normalize(normal), settings['translate'])
            p1_t = cpv.sub(p1, vtran)
            p2_t = cpv.sub(p2, vtran)
            p3_t = cpv.sub(p3, vtran)
            print("New coordinates are:")
            print_info("New", p1_t, p2_t, p3_t)
            print("New coordinates are for normalized plane:")
            v1_t = cpv.normalize(cpv.sub(p1_t, p2_t))
            v2_t = cpv.normalize(cpv.sub(p3_t, p2_t))
            normal_t = cpv.normalize(cpv.cross_product(v1_t, v2_t))
            v2_t = cpv.normalize(cpv.cross_product(normal_t, v1_t))
            p1_t2 = cpv.add(v1_t, p2_t)
            p3_t2 = cpv.add(v2_t, p2_t)
            print_info("Newnormal", p1_t2, p2_t, p3_t2)
    
        if vm1!=None:
            v1 = cpv.scale(cpv.normalize(v1), vm1)
        if vm2!=None:
            v2 = cpv.scale(cpv.normalize(v2), vm2)
    
        centrum = p2
        if center:
            corner1 = cpv.add(cpv.add(centrum, v1), v2)
            corner2 = cpv.sub(cpv.add(centrum, v1), v2)
            corner3 = cpv.sub(cpv.sub(centrum, v1), v2)
            corner4 = cpv.add(cpv.sub(centrum, v1), v2)
        else:
            corner1 = cpv.add(cpv.add(centrum, v1), v2)
            corner2 = cpv.add(centrum, v1)
            corner3 = centrum
            corner4 = cpv.add(centrum, v2)
    
        return plane(corner1, corner2, corner3, corner4, normal, settings)
    
    
    def print_info(name, coor1, coor2, coor3):
        cs1 = (map(float, [ '%.2f' % elem for elem in coor1 ]) )
        cs2 = (map(float, [ '%.2f' % elem for elem in coor2 ]) )
        cs3 = (map(float, [ '%.2f' % elem for elem in coor3 ]) )
        print("You can also use the function calls with these coordinates")
        print("plane.make_plane_points(name='%s', l1=%s, l2=%s, l3=%s)"%(name, cs1, cs2, cs3))
    
    
    def make_plane(name,a1='(pk1)',a2='(pk2)',a3='(pk3)', vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
        """
        DESCRIPTION
        Create a CGO plane from three atomic coordinates
    
        USAGE
        make_plane name, a1, a2, a3
    
        where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
        """
        # get coordinates from atom selections
        coor1 = cmd.get_model(a1).get_coord_list()[0]
        coor2 = cmd.get_model(a2).get_coord_list()[0]
        coor3 = cmd.get_model(a3).get_coord_list()[0]
    
        # Help with alternative
        print_info(name, coor1, coor2, coor3)
    
        # Get the plane
        plane = planeFromPoints(p1=coor1, p2=coor2, p3=coor3, vm1=vm1, vm2=vm2, center=center, settings=settings)
        makePrimitive(plane, name)
        #cmd.show("cgo", "plane*")
    
        if makepseudo:
            cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=coor1)
            cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=coor2)
            cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=coor3)
    
    # Extend function to be called inside pymol
    cmd.extend("make_plane", make_plane)
    
    def make_plane_points(name,l1=None,l2=None,l3=None, vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
        """
        DESCRIPTION
        Create a CGO plane from three atomic coordinates
    
        USAGE
        make_plane name, l1, l2, l3
    
        where each xys is a list with floats of x,y,z coordinates
        """
        if l1==None or l2==None or l3==None:
            print("Please provide a list of xyz floats for each 3 positions")
            return
        if type(l1) is not list or type(l2) is not list or type(l3) is not list:
            print(type(l1),type(l2),type(l3))
            print("Please provide 3 list of xyz floats for each 3 positions")
            return
    
        plane = planeFromPoints(p1=l1, p2=l2, p3=l3, vm1=vm1, vm2=vm2, center=center, settings=settings)
        makePrimitive(plane, name)
    
        if makepseudo:
            cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=l1)
            cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=l2)
            cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=l3)
    
    # Extend function to be called inside pymol
    cmd.extend("make_plane_points", make_plane_points)
    
    class PlaneWizard(Wizard):
        def __init__(self):
            Wizard.__init__(self)
    
            # some attributes to do with picking
            self.pick_count = 0
            self.object_count = 0
            self.object_prefix = "pw"
    
            self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
            cmd.set("mouse_selection_mode",0) # set selection mode to atomic
            cmd.deselect()
    
        def reset(self):
            cmd.delete(self.object_prefix + "*")
            cmd.delete("sele*")
            cmd.delete("_indicate*")
            cmd.unpick()
            cmd.refresh_wizard()
    
        def delete_all(self):
            cmd.delete("plane*")
    
        def cleanup(self):
            cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
            self.reset()
            self.delete_all()
    
        def get_prompt(self):
            self.prompt = None
            if self.pick_count == 0:
                self.prompt = [ 'Please click on the first atom...']
            elif self.pick_count == 1:
                self.prompt = [ 'Please click on the second atom...' ]
            elif self.pick_count == 2:
                self.prompt = [ 'Please click on the third atom...' ]
            return self.prompt
    
        def do_select(self, name):
            # "edit" only this atom, and not others with the object prefix
            try:
                cmd.edit("%s and not %s*" % (name, self.object_prefix))
                self.do_pick(0)
            except pymol.CmdException, pmce:
                print pmce
    
        def pickNextAtom(self, atom_name):
            # transfer the click selection to a named selection
            cmd.select(atom_name, "(pk1)")
    
            # delete the click selection
            cmd.unpick()
    
            # using the magic of indicate, highlight stuff
            indicate_selection = "_indicate" + self.object_prefix
            cmd.select(indicate_selection, atom_name)
            cmd.enable(indicate_selection)
    
            self.pick_count += 1
            self.error = None
    
            # necessary to force update of the prompt
            cmd.refresh_wizard()
    
        def do_pick(self, picked_bond):
    
            # this shouldn't actually happen if going through the "do_select"
            if picked_bond:
                self.error = "Error: please select bonds, not atoms"
                print self.error
                return
    
            atom_name = self.object_prefix + str(self.pick_count)
            if self.pick_count < 2:
                self.pickNextAtom(atom_name)
            else:
                self.pickNextAtom(atom_name)
    
                point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
                point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
                point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
                plane = planeFromPoints(point1, point2, point3)
    
                planeName = "plane-%02d" % self.object_count
    
                print_info(planeName, point1, point2, point3)
    
                self.object_count += 1
                makePrimitive(plane, planeName)
                cmd.show("cgo", "plane*")
    
                self.pick_count = 0
                self.reset()
    
        def get_panel(self):
            return [
                [ 1, 'Plane Wizard',''],
                [ 2, 'Reset','cmd.get_wizard().reset()'],
                [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
                [ 2, 'Done','cmd.set_wizard()'],
            ]
    

### Examples

Plane in a protein 
    
    
    reinitialize
    
    fetch 1ubq, async=0
    show_as cartoon, all 
    
    import plane
    # Start the wizard, and do manual picking
    cmd.set_wizard(plane.PlaneWizard())
    

Or by selection 
    
    
    reinitialize
    
    fetch 1ubq, async=0
    show_as cartoon, all 
    
    import plane
    
    # Set alpha level and color
    #violetpurple: 0.55, 0.25, 0.60 #yellow: 1.0, 1.0, 0.0, #blue: 0.0, 0.0, 1.0 #orange: 1.0, 0.5, 0.0 #forest: 0.2, 0.6, 0.2 #red: 1.0, 0.0, 0.0
    dict = {'ALPHA':0.6, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':True}
    plane.make_plane(name='test', a1='/1ubq//A/24/CA', a2='/1ubq//A/29/CA', a3='/1ubq//A/40/CA', center=False, settings=dict)
    

Or by atom coordinates 
    
    
    reinitialize
    
    fetch 1ubq, async=0
    show_as cartoon, all 
    
    import plane
    
    # Or  make from atom coordinates
    #plane.make_plane_points(name='test', l1=[35.03, 21.72, 17.07], l2=[37.47, 27.39, 10.67], l3=[37.74, 31.64, 23.71])
    
    # Define plane, 10 angstrom in length
    #plane.make_plane_points(name='p1', l1=[0.0, 10.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 0.0, 10.0], center=False, makepseudo=False)
    

Or make a color cube 
    
    
    reinitialize
    #violetpurple: 0.55, 0.25, 0.60 #yellow: 1.0, 1.0, 0.0, #blue: 0.0, 0.0, 1.0 #orange: 1.0, 0.5, 0.0 #forest: 0.2, 0.6, 0.2 #red: 1.0, 0.0, 0.0
    
    # YZ Plane, #purple
    dict = {'ALPHA':0.4, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':True}
    plane.make_plane_points(name='p1', l1=[0.0, 10.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 0.0, 10.0], center=False, makepseudo=False, settings=dict)
    # YZ Plane, shifted in X, #yellow
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 1.0, 0.0]}
    plane.make_plane_points(name='p6', l1=[10.0, 10.0, 0.0], l2=[10.0, 0.0, 0.0], l3=[10.0, 0.0, 10.0], center=False, makepseudo=False, settings=dict)
    
    # XZ Plane, blue
    dict = {'ALPHA':0.4, 'COLOR':[0.0, 0.0, 1.0]}
    plane.make_plane_points(name='p2', l1=[10.0, 0.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 0.0, 10.0], center=False, makepseudo=False, settings=dict)
    # XZ Plane, shifted in Y, #orange
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.5, 0.0], 'INVERT':True}
    plane.make_plane_points(name='p5', l1=[10.0, 10.0, 0.0], l2=[0.0, 10.0, 0.0], l3=[0.0, 10.0, 10.0], center=False, makepseudo=False, settings=dict)
    
    # XY Plane, forest
    dict = {'ALPHA':0.4, 'COLOR':[0.2, 0.6, 0.2], 'INVERT':True}
    plane.make_plane_points(name='p4', l1=[10.0, 0.0, 0.0], l2=[0.0, 0.0, 0.0], l3=[0.0, 10.0, 0.0], center=False, makepseudo=False, settings=dict)
    # XY Plane, shifted in Z, red
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.0, 0.0]}
    plane.make_plane_points(name='p3', l1=[10.0, 0.0, 10.0], l2=[0.0, 0.0, 10.0], l3=[0.0, 10.0, 10.0], center=False, makepseudo=False, settings=dict)
    
    zoom all
    

Or make a color cube, by initial fixed plane 
    
    
    import plane
    
    python
    from chempy import cpv
    p1 = [23.76, -47.69, 45.23]
    p2 = [34.96, -18.57, -1.25]
    p3 = [90.76, -4.31, 21.69]
    v1 = cpv.sub(p1, p2)
    v2 = cpv.sub(p3, p2)
    normal = cpv.cross_product(v1, v2)
    normal_norm = cpv.normalize(normal)
    v3 = cpv.scale(normal_norm, 40)
    p4 = cpv.add(p2, v3)
    cmd.pseudoatom("mine", color="tv_blue", pos=p4)
    python end
    
    # XY Plane, forest
    dict = {'ALPHA':0.4, 'COLOR':[0.2, 0.6, 0.2], 'INVERT':True}
    plane.make_plane_points(name='p4', l1=p1, l2=p2, l3=p3, center=False, makepseudo=False, settings=dict)
    # XY Plane, shifted in Z, red
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.0, 0.0]}
    plane.make_plane_points(name='p3', l1=cpv.add(p1, v3), l2=cpv.add(p2, v3), l3=cpv.add(p3, v3), center=False, makepseudo=False, settings=dict)
    
    # XZ Plane, blue
    dict = {'ALPHA':0.4, 'COLOR':[0.0, 0.0, 1.0]}
    plane.make_plane_points(name='p2', l1=p1, l2=p2, l3=cpv.add(cpv.sub(p3, v2), v3), center=False, makepseudo=False, settings=dict)
    # XZ Plane, shifted in Y, #orange
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 0.5, 0.0], 'INVERT':True}
    plane.make_plane_points(name='p5', l1=cpv.add(p1, v2), l2=cpv.add(p2, v2), l3=cpv.add(p3, v3), center=False, makepseudo=False, settings=dict)
    
    # YZ Plane, #purple
    dict = {'ALPHA':0.4, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':True}
    plane.make_plane_points(name='p1', l1=cpv.add(cpv.sub(p1, v1), v3), l2=p2, l3=p3, center=False, makepseudo=False, settings=dict)
    ## YZ Plane, shifted in X, #yellow
    dict = {'ALPHA':0.4, 'COLOR':[1.0, 1.0, 0.0]}
    plane.make_plane_points(name='p6', l1=cpv.add(p1, v3), l2=cpv.add(p2, v1), l3=cpv.add(p3, v1), center=False, makepseudo=False, settings=dict)
    
    zoom all
    

Retrieved from "[https://pymolwiki.org/index.php?title=Plane_Wizard&oldid=12500](https://pymolwiki.org/index.php?title=Plane_Wizard&oldid=12500)"


---

## Radius of gyration

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.querying#gyradius](http://pymol.org/psicoredirect.php?psico.querying#gyradius)  
  
This script does calcualate the [radius of gyration](http://en.wikipedia.org/wiki/Radius_of_gyration) of a molecule. 

Thanks to Tsjerk Wassenaar for posting this script on the [pymol-users mailing list](http://sourceforge.net/mailarchive/message.php?msg_id=27288491)! 

_The function was added to[psico](/index.php/Psico "Psico") as "gyradius"._
    
    
    from pymol import cmd
    import math
    
    def rgyrate(selection='(all)', quiet=1):
        '''
    DESCRIPTION
    
        Radius of gyration
    
    USAGE
    
        rgyrate [ selection ]
        '''
        try:
            from itertools import izip
        except ImportError:
            izip = zip
        quiet = int(quiet)
        model = cmd.get_model(selection).atom
        x = [i.coord for i in model]
        mass = [i.get_mass() for i in model]
        xm = [(m*i,m*j,m*k) for (i,j,k),m in izip(x,mass)]
        tmass = sum(mass)
        rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in izip(x,xm))
        mm = sum((sum(i)/tmass)**2 for i in izip(*xm))
        rg = math.sqrt(rr/tmass - mm)
        if not quiet:
            print("Radius of gyration: %.2f" % (rg))
        return rg
    
    cmd.extend("rgyrate", rgyrate)
    
    # vi:expandtab
    

### Example

[![Radius of gyration Example.png](/images/2/2b/Radius_of_gyration_Example.png)](/index.php/File:Radius_of_gyration_Example.png)

[](/index.php/File:Radius_of_gyration_Example.png "Enlarge")

This example uses the radius of gyration and the [center of mass](/index.php/Center_Of_Mass "Center Of Mass") to display a semitransparent sphere. 
    
    
    from pymol import cmd, stored, util
    import centerOfMass, radiusOfGyration
    
    cmd.set('sphere_transparency', 0.5)
    
    cmd.fetch('2xwu')
    cmd.hide('everything')
    cmd.show('cartoon', 'chain A')
    
    r = radiusOfGyration.rgyrate('chain A and polymer')
    centerOfMass.com('chain A and polymer', object='com', vdw=r)
    
    util.cbc()
    

### See Also

  * [centerofmass](/index.php/Centerofmass "Centerofmass")
  * [Center Of Mass](/index.php/Center_Of_Mass "Center Of Mass")
  * [COM](/index.php/COM "COM")



Retrieved from "[https://pymolwiki.org/index.php?title=Radius_of_gyration&oldid=13172](https://pymolwiki.org/index.php?title=Radius_of_gyration&oldid=13172)"


---

## RotationAxis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/draw_rotation_axis.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/draw_rotation_axis.py)  
Author(s)  | [Pablo Guardado Calvo](/index.php?title=User:PabloGuardado&action=edit&redlink=1 "User:PabloGuardado \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script will draw a CGO cylinder representing a rotation axis for a given transformation. It's very useful for drawing the axes of rotational symmetry in an oligomeric assembly. 

The idea is to align two molecules/domains/chains/selections (using cmd.super) and extract the trasformation (TTT) matrix (T). The direction of the rotation axis, and a point are obtained from T and are used to create a cgo object representing the axis. The script generates two measures: one in the graphical screen (the rotation axis value, and the norm of the translation vector along the rotation axis) and some basic information in the command-line (the transformation matrix, the rotation angle, distance between centroids, and some pml lines that you can use to reproduce the axis...) 

As always with these type of things, you have to use it at your own risk. I did not try all possible combinations, but if you find a bug, do not hesitate to contact me (pablo.guardado@gmail.com) or try to modify the code for yourself to correct it. 

To load the script just type: 

**run _path-to-the-script_ /draw_rotation_axis.py**

or if you want something more permanent add the previous line to your .pymolrc file 

The script works just typing: 

**draw_axis('selection1', 'selection2')**

Please, pay attention to the apostrophes around the selections, you MUST use them. Also works with chains: 

**draw_axis('chain A', 'chain B')**

or objects: 

**draw_axis('obj01', 'obj02')**

Also, you can play a bit with the length, width and color of the axis. 

**draw_axis('selection1', 'selection2', scale_factor, width, r1, g1, b1, r2, g2, b2)**

**scale_factor** = to control the length of the axis, the default is 20 

**width** = to control the width of the axis. Default is 0.6 

**r1, g1, b1** = first RGB color code. Default is 1, 1, 1 

**r2, g2, b2** = second RGB color code to create the gradient. Default is 1, 0, 0. 

To create a single color axis, just made r1,g1,b1=r2,g2,b2 

# Examples
    
    
    # download the source and save as draw_rotation_axis.py
    run draw_rotation_axis.py
    fetch 2vak
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/d/df/2vak4.png)](/index.php/File:2vak4.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/3/37/2vak2.png)](/index.php/File:2vak2.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Another example, calculate the rotation axis of an homodimer
    run draw_rotation_axis.py
    fetch 3rfg
    # calculate rotation axis between chains A and B
    draw_axis('chain A', 'chain B')
    

  * [![Rotation axis between chains A and B](/images/5/55/3rfg1.png)](/index.php/File:3rfg1.png "Rotation axis between chains A and B")

Rotation axis between chains A and B 

  * [![Some basic information is printed in the screen](/images/d/d6/3rfg4.png)](/index.php/File:3rfg4.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 



    
    
    # Clearly, each of the domains are made up with motifs with internal symmetry
    draw_axis('resi 200-236 and chain A', 'resi 238-274 and chain A', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # Also, you can create first the selections and use them to calculate the axis
    sele selection1, resi 200-236 and chain A
    sele selection2, resi 238-274 and chain A
    draw_axis('selection1', 'selection2', 20, 0.6, 1, 0, 0, 1, 0, 0)
    # will produce the same result.
    

  * [![Internal rotation axis](/images/6/60/3rfg2.png)](/index.php/File:3rfg2.png "Internal rotation axis")

Internal rotation axis 

  * [![Some basic information is printed in the screen](/images/0/03/3rfg3.png)](/index.php/File:3rfg3.png "Some basic information is printed in the screen")

Some basic information is printed in the screen 




Retrieved from "[https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923](https://pymolwiki.org/index.php?title=RotationAxis&oldid=13923)"


---

## Slerpy

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Manual for Slerpy.py
    * 1.1 General Use
    * 1.2 Important concepts
    * 1.3 Quick Start Tutorial
    * 1.4 Tips and Tricks
      * 1.4.1 Converting scenes to movies
      * 1.4.2 Starting off right
      * 1.4.3 Live pymol presentations during a talk
      * 1.4.4 Pausing on a view
      * 1.4.5 Morphing and multi-state models
    * 1.5 Command Reference
    * 1.6 Why is it called slerpy?
  * 2 Script Files
    * 2.1 Installation
    * 2.2 Code
      * 2.2.1 slerpy.py
      * 2.2.2 movs.py



# Manual for Slerpy.py

An extension to pymol that creates a moderately easy to use environment for doing keyframe animation. 

  


## General Use

At the pymol command line type: 
    
    
    import slerpy
    

This will load the extended commands. All commands in slerpy begin with the letter s. Pymol's tab autocomplete feature will work on the additional commands. 

## Important concepts

The main function of slerpy is to record a series of pymol views. A movie can then be created by interpolating between these views. A pymol view consists mostly of the camera orientation (that is, the orientation of the viewers eye with respect to the molecule). It also includes information about the clipping planes. 

It is important to realize that most slerpy commands act on the "current" view. You can navigate among views by using the sn, sp, and sgo commands. If you've changed the view around with the mouse or you just want to know the number of the current view you can get back to the current view with the sc command. 

Pymol views do not contain information about how pymol objects and selections are displayed. If you just want to create a movie which moves around a single representation of a molecule then all you need to do is record the set of views that define the tour of the molecule. 

If, on the other hand, you want to change the representation or change which items are displayed you will need to add actions to some of your views. An action is any set of pymol commands. Actions can be associated with any view in the series recorded by slerpy. The _sscene_ command inserts a new view and simultaneously creates a pymol scene and the action to display it. The scene will include all of the objects and representations visible at the time the command was issued. 

In order to control the rate of motion between the defined views in a slerpy movie, you can control the number of frames used in each interpolation. When a view is saved in slerpy it is associated by default with a transition of 50 frames to the next view. The number of frames in the transition can be altered with the _ssetf_ command in slerpy. 

The views and actions stored by slerpy can (and should) be saved to a key file with the _swrite_ command. They can then be retrieved with the sread command. Note that swrite saves the current pymol state in a .pse file but sread does not read in the .pse file. If you're starting a new pymol session to continue work on an existing movie you should load the pse file before doing an sread. 

## Quick Start Tutorial

If you haven't installed slerpy yet see Installation

Step 1
    You probably want to start off in a nice clean working directory that just has the coordinate files you want to work with.

Step 2
    Read in your molecule(s) and create the various selections and representations that you want to include in the movie.

Step 3
    At the pymol prompt, type:
    
    
    import slerpy
    

Step 4
    Get your molecule in exactly the orientation and representation that you want to use for the beginning of your movie.

Step 5
    Type:
    
    
    sinsert
    

Step 6
    Using the mouse, move your molecule to the next orientation that you want to use. When you record the movie, the camera orientation will be interpolated between each consecutive pair of views. This can include changes in rotation, zooming, clipping etc.

    Loop back to Step 5. Continue this until you've got all your orientations stored.

    You can check how any set of transitions will look at any time by using the _sshow_ command (see Command Reference for details).

    You can adjust the rate of any transition using the _ssetf_ command

Step 7A
    Add any actions to your views using the saction command. Any number of pymol commands can be strung together separated by semicolons. If, for example, you want to change your protein from a cartoon to a surface and add a ligand when you get to view 5 you would do the following (assuming you've defined the pymol selections prot and lig):
    
    
    sgo 5
    saction "hide cartoon, prot; show surface, prot; show sticks, lig"
    

Step 7B (Alternative using scenes)
    
    
    
    sgo 5
    

    Now use the gui to create the representation you want and then:
    
    
    sscene
    

Step 8
    Save everything. Save your slerpy views and actions as well as a pse file of your current pymol session:
    
    
    swrite mymovie
    

This will create mymovie.key, which has all of the views, frame counts, actions etc. and mymovie.pse, the associated pymol session file. 

Step 9
    Record the movie! Type:
    
    
    srecord
    

    You can then play the movie by typing the standard pymol command _mplay_ or by clicking the play button in pymol.

Step 10
    If you save the pymol session again, the pse file will contain the movie which can then be shown immediately after startup without running slerpy.py. Note that pymol will warn you when you load a pse file that contains a movie.

Step 11
    If you want to, the movie can be exported using the mpng command (see the pymol [documentation](http://pymol.sourceforge.net/newman/ref/S1000comref.html#2_105)). Also, see the useful article [Making Movies](/index.php/Making_Movies "Making Movies").

## Tips and Tricks

### Converting scenes to movies

You can just step through the scenes and type sscene for each one. This will create a duplicate slerpy scene for each of the scenes you'd already saved but that's not such a disaster. Be sure to swrite when you're done. 

Note that there is some overhead associated with recalling scenes. To avoid pauses at view transitions, I prefer to actually issue the set of show and hide commands that will generate the scene rather than using the above method. 

### Starting off right

It's a bit of a pain, but I like to associate with the first frame of the movie an action list that hides everything and then turns on all the objects that I want to have visible at the beginning. This ensures that when your movie loops back to the beginning it will look the same as it did the first time through. For example: 
    
    
    sgo 0
    saction "hide everything; show lines, prot; show surface, activesite; show sticks, ligand"
    

Alternatively, start your slerpy work with an _sscene_. 

### Live pymol presentations during a talk

Be sure to run your movie once it's been opened and before your presentation if you're presenting in pymol. This will ensure that any objects that don't appear until the middle of the movie are available in memory and won't need to be rebuilt while your audience waits. 

Of course showing your movie from within pymol allows you to show movies in stereo if you've got a presentation system that allows this. If you pass out stereo glasses first it's guaranteed that everyone will remember your talk... 

### Pausing on a view

Just sgo to the view you want to stay on for a while and do an sinsert. This will insert a new view with the same orientation etc as the one you were just on. You can adjust the length of the pause by changing the number of frames for the transistion between these two identical views using the ssetf command. 

### Morphing and multi-state models

Morphs and multi-conformation models are represented in pymol as single objects with multiple states. Cycling through the states as an animation is very straightforward in pymol. Slerpy allows you to include this animation over the states of an object as a transition between slerpy views within the context of larger movie. This is done using the _smorph_ command (see Command Reference). _smorph_ allows you to specify the starting and ending pymol state numbers to use during the transition from the current to the next view. It will often be appropriate to have the movie continue after the morph using the final state of the morphing model, that is, the conformation to which you morphed. Since the default for a typical slerpy view is to show the first state of an object, you'll probably need to have available, in addition to your multi-state model, a single-state object containing the final conformation. You can then hide the multistate object and show the single-state final conformation as an action associated with the final view of your morphing sequence. 

To generate morphs you can use the [rigimol](http://www.delanoscientific.com/rigimol.html) program provided by Warren as part of the incentive pymol package for subscribers or use [lsqman](http://alpha2.bmc.uu.se/usf/mol_morph.html) from Gerard Kleywegt. 

## Command Reference

Note that it is essential to understand that slerpy uses the concept of a current or active view. This is the element in the list to which most commands are applied. It is not necessarily the view that is currently visible on the screen. It is advisable to use the sc command frequently to make sure you really know which view your command is being applied to. 

**saction _string_ :** Assoiciates the pymol commands in _string_ with the current view. _string_ must be enclosed in quotes and must contain only valid pymol commands separated by semicolons. 

**sappend** : Add the currently displayed view at the end of the list of views in slerpy 

**sappendaction _string_ :** Adds the action(s) in _string_ to the list of actions associated with the current view 

**sc** : Go to the slerpy active view 

**scrossfade _startobject, endobject_ :** During the transition from the current view to the next view _startobject_ will fade out and _endobject_ will fade in. The two objects must be shown as sticks. They must be objects, not merely selections, as of pymol 0.99. 

**sdelaction** : Deletes the action associated with the current view. Be sure you're sure which view you're on before you use this. This will also delete any actions etc. associated with the view so be careful. 

**sdelete** : Remove the slerpy active view from the list of views to be interpolated. 

**sdeletesetting:** Remove the setting interpolation for the current view. 

**sdumpactions:** List all actions by frame. 

**sgo _n_ :** Change the slerpy active view to view _n_. 

**sinsert** : Add the currently displayed view after the slerpy active view. 

**sinterpsetting _setting, selection, startval, endval_ :** The pymol setting _setting_ will be interpolated linearly from _startval_ to _endval_ during the transition from the current view to the next view. You can only interpolate one setting per transition. This is the hard way to, for example, fade out an object: 
    
    
    sinterpsetting stick_transparency, lig, 0.0, 1.0
    

**slist** : List all views stored for interpolation. Also lists the number of frames in each transition. 

**smorph _startmodel,endmodel_ :** The transition from the current frame to the next frame will consist of one frame per pymol state starting with state _startmodel_ and ending with state _endmodel_. Subsequent frames (i.e. from subsequent views) will revert to state 1. The state numbers apply to currently visible objects so you will most likely want to have an object with your starting conformation, an object with your multi-state morphing model, and an object with your final conformation. You would then sgo to the frame where you want the morph to start and add an action to hide the starting conformation object and show the multi-model morphing object, do an smorph 1,30 or whatever the number of states in your morph is, append another frame and give it an action where the multi-state model is hidden and the final conformation is shown. 

**sn** : Go to the next view 

**sp** : Go to the previous view. 

**sread _filename_ :** Restore all the information written with swrite. Does not read in the pse file (to avoid inadvertantly writing over some new selections, scenes or whatever). 

**srecord** : Records the movie 

**sreplace** : Replace the slerpy current view with the currently displayed view. 

**sscene** : Add a the currently displayed view after the slerpy active view and create a scene to go with it. The current pymol state, including which objects are displayed and how they are shown will be captured. 

**ssetf** : Set the number of frames to use in the transition from the slerpy active view to the next view 

**ssetupview _n_ :** This attempts to make sure that all objects are displayed (or not) as they would be when view _n_ is arrived at in the movie. It goes through and executes all of the sactions from all preceeding views. For some reason this doesn't always work in complex cases. 

**sshow _n,m_ :** This command records and shows a segment of the movie showing the transitions starting with view n and ending with view m. If the arguments m and n are omitted the transition from the current view to the next view will be shown. 

**swrite _filename_ :** Writes all the information used by slerpy to a file _filename_.key. It also writes a pymol session file _filename_.pse containing all of your current objects, selections etc. The format of the .key file is supposed to be human readable and it can be convenient to edit this file rather than using saction over and over. After editing it by hand be sure to do an sread. 

## Why is it called slerpy?

slerp is an acronym of sorts for spherical linear interpolation. The transition between two views approximates spherical linear interpolation of the their camera positions using quaternions. See the Wikipedia [article](http://en.wikipedia.org/wiki/Slerp) on slerps. 

# Script Files

## Installation

To install, copy (by pasting into a text editor) the following two python code segments to files named slerpy.py and movs.py, respectively. Place these files in the pymol/modules/pymol/ directory created by the pymol installation. 

If you don't have write priviledges in the pymol installation directories you can just copy these files to the working directory from which pymol will be run. 

## Code

### slerpy.py

    The command definitions for slerpy
    
    
    ################################################################################
    #slerpy.py - Command definition routines for slerpy
    #Copyright (C) 2006 Joel Bard
    #
    #This program is free software; you can redistribute it and/or
    #modify it under the terms of the GNU General Public License
    #as published by the Free Software Foundation; either version 2
    #of the License, or (at your option) any later version.
    #
    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.
    #
    #You should have received a copy of the GNU General Public License
    #along with this program; if not, write to the Free Software
    #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
    #################################################################################
    
    from pymol import cmd
    import movs
    
    def readViews( filename ):
    #obsolete
        vfile = open( filename, 'r')
        views = [] #a list of views each of which is a list of 18 numbers
        vstrings = vfile.readlines()
        for vstring in vstrings:
            vals = vstring.split()
            view = [float(v) for v in vals]
            views.append( view )
        vfile.close()
        return views
    
    def readFrames( filename ):
    #obsolete
        ffile = open( filename, 'r' )
        frames = []
        fstrings = ffile.readlines()
        for fstring in fstrings:
            frames.append( int(fstring) )
        ffile.close()
        return frames
    
    def readActions( filename ):
    #obsolete
        #actions are stored in file where
        #for each line, the first 4 chars are the view index
        #associated with the action and the rest of the
        #line is the pymol command to be executed
        #upon reading, a dictionary is returned
        #with view indices as keys and actions as values
        actions = {}
        try:
            afile = open( filename, 'r' )
        except:
            print "No actions for this project"
            return actions
        astrings = afile.readlines()
        for astring in astrings:
            try:
                aindex = int(astring[:4])
                action = astring[4:]
                actions[ aindex ] = action[:-1]
            except:
                print "empty line"
        afile.close()
        return actions
    
    def readModels( filename ):
    #obsolete
        models = {}
        try:
            mfile = open( filename, 'r' )
        except:
            print "No models for this project"
            return models
        mstrings = mfile.readlines()
        for mstring in mstrings:
            try:
                mindex = int(mstring[:4])
                model = mstring[4:]
                models[ mindex ] = model[:-1]
            except:
                print "empty line"
        mfile.close()
        return models
    
    def readSettings( filename ):
    #obsolete
        settings = {}
        try:
            sfile = open( filename, 'r' )
        except:
            print "No settings for this project"
            return settings
        sstrings = sfile.readlines()
        for sstring in sstrings:
            try:
                sindex = int(sstring[:4])
                scommas = sstring[4:]
                settingName,selection,startVal,endVal = scommas.split(',')
                setting = [settingName,selection,float(startVal),float(endVal)]
                settings[sindex] = setting
            except:
                print "unable to parse setting"
        sfile.close()
        return settings
                
    def readScenes( filename ):
    #obsolete
        global scene_counter
    
        scene_counter = 0
        scenes = {}
        try:
            sfile = open( filename, 'r' )
        except:
            print "No scenes file for this project"
            return scenes
        sstrings = sfile.readlines()
        for sstring in sstrings:
            try:
                sindex = int(sstring[:4])
                scene = sstring[4:]
                scenes[ sindex ] = scene[:-1]
                scene_counter += 1
                #print "reading scene", sstring, sindex, scene
            except:
                print "empty line"
        sfile.close()
        return scenes
    
    def read_all( fileroot ):
    #obsolete in favor of readKeyViewFile
        global views
        global frames
        global actions
        global cur_view
        global cur_index
        global scenes
        global models
        global settings
        
        views = readViews( fileroot+".txt" )
        frames = readFrames( fileroot+".frm")
        actions = readActions( fileroot+".act")
        scenes = readScenes( fileroot+".scn")
        models = readModels( fileroot+".mod")
        settings = readSettings( fileroot+".set")
        cur_view = views[0]
        cur_index = 0
        show_cur()
        
    def print_view( view ):
        for i in range(0,6):
            for j in range(0,3):
                print "%12.6f"% view[ 3*i+j ] ,
            print
    
    def show_next():
        global cur_index
        global cur_view
        global views
        if cur_index == len( views )-1:
            print "No more views."
            return
        cur_index += 1
        cur_view = views[ cur_index ]
        cmd.set_view(cur_view)
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ",cur_index
    
    def show_prev():
        global cur_index
        global cur_view
        global views
        if cur_index == 0:
            print "No more views."
            return
        cur_index -= 1
        cur_view = views[ cur_index ]
        cmd.set_view(cur_view)
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ",cur_index
    
    def show_cur():
        global cur_index
        global cur_view
        global views
        cur_view = views[ cur_index ]
        cmd.set_view(cur_view)
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ",cur_index
    
    def go_to_view( arg=0 ):
        global cur_index
        global cur_view
        global views
        n = int( arg )
        if n < 0 or n >= len(views):
            print "Index out of range."
            return
        cur_index = n
        cur_view = views[n]
        cmd.set_view( cur_view )
        print "Matrix: "
        print_view( cur_view )
        print "Showing view number ", cur_index
    
    def insert_current():
        #insert the current view into the list after the view
        #in views[cur_index]
        #set frames to default
        global cur_index
        global cur_view
        global views
        global frames
        global actions
        global scenes
        global settings
        global models
        global fades
        
        cur_index += 1
        cur_view = cmd.get_view()
        views.insert( cur_index, [cv for cv in cur_view] )
        frames.insert( cur_index, 50 )
        
        #deal with actions dictionary
        actions = incKeyAbove( actions, cur_index )
        scenes = incKeyAbove( scenes, cur_index )
        settings = incKeyAbove( settings, cur_index )
        models = incKeyAbove( models, cur_index )
        fades = incKeyAbove( fades, cur_index )
    
        print "New view:"
        print_view( cur_view )
        print "Inserted view with index", cur_index, "and a 50 frame transition"
    
    def append_view():
        global views
        global frames
        global cur_index
        global cur_view
    
        cur_index = len(views)
        cur_view = cmd.get_view()
        views.append( [cv for cv in cur_view] )
        frames.append( 50 )
    
        print "New view: "
        print_view( cur_view )
        print "Appended view with index", cur_index, "and a 50 frame transition"
        print "The current view is", cur_index
        
    def incKeyAbove( dict, index ):
        tempDict = {}
        for key, val in dict.iteritems():
            if key >= index:
                newkey = key + 1
                tempDict[newkey] = val
            else:
                tempDict[key] = val
        return tempDict
    
    def decKeyAbove( dict, index ):
        tempDict = {}
        for key, val in dict.iteritems():
            if key > index:
                newkey = key - 1
                tempDict[newkey] = val
            else:
                tempDict[key] = val
        return tempDict
        
    def delete_current():
        #remove the current view from the list
        #show the previous view
        global cur_index
        global cur_view
        global views
        global actions
        global scenes
        global settings
        global models
        global frames
        global fades
        
        del views[cur_index]
        del frames[cur_index]
        if cur_index in actions:
            del actions[cur_index]
        if cur_index in scenes:
            del scenes[cur_index]
        if cur_index in settings:
            del settings[cur_index]
            
        #deal with dictionaries
        actions = decKeyAbove( actions, cur_index )
        scenes = decKeyAbove( scenes, cur_index )
        settings = decKeyAbove( settings, cur_index )
        models = decKeyAbove( models, cur_index )                               
        fades = decKeyAbove( fades, cur_index )
    
        print "View number",cur_index,"deleted."
        if cur_index > 0:
            cur_index -= 1
        cur_view = views[cur_index]
        cmd.set_view( cur_view )
        print "Current view is number",cur_index
    
    def delete_settings():
        global settings
        global cur_index
        del settings[cur_index]
    
    def replace_current():
        global cur_index
        global cur_view
        global views
        cur_view = cmd.get_view()
        views[cur_index] = [cv for cv in cur_view]
    
    def insert_scene():
        global views
        global actions
        global settings
        global frames
        global cur_index
        global cur_view
        global scenes
        global scene_counter
        global models
        global fades
     
        cur_index += 1
        cur_view = cmd.get_view()
        views.insert( cur_index, [cv for cv in cur_view] )
        frames.insert( cur_index, 50 )
        
        #deal with dictionaries
        actions = incKeyAbove( actions, cur_index )
        scenes = incKeyAbove( scenes, cur_index )
        settings = incKeyAbove( settings, cur_index )
        models = incKeyAbove( models, cur_index )
        fades = incKeyAbove( fades, cur_index )
    
        #this stuff has to be done after the above
        #find a free scene name
        i = 1
        while True:
            for sname in scenes.values():
                print "|"+sname+"|"
                if sname == "slerpy_"+str(i):
                    break
            else:
                break
            i += 1
        newname = "slerpy_"+str(i)
            
        scene_counter += 1
        cmd.scene( newname, "store" )
        scenes[ cur_index ] = newname
        actions[ cur_index ] = "scene "+newname
        
        print "New view:"
        print_view( cur_view )
        print "Inserted view with index", cur_index, "and a 50 frame transition"
        print "Added scene",newname
        
    def write_views( filename ):
    #deprecated in favor of key files
        global views
        global frames
        global scenes
        global actions
        global settings
        global models
        global fades
    
        viewfile = open( filename+".txt", 'w')
        for view in views:
            for v in view:
                viewfile.write( str(v) + " " )
            viewfile.write('\n')
        viewfile.close()
        
        framefile = open( filename+".frm", 'w' )
        for frame in frames:
            framefile.write( str( frame )+'\n' )
        framefile.close()
    
        actionfile = open( filename+".act", 'w' )
        for key,action in actions.iteritems():
            keystring = str( key )
            actionfile.write( keystring.rjust(4)+action + '\n' )
        actionfile.close()
    
        scenefile = open( filename+".scn", 'w' )
        for key,scene in scenes.iteritems():
            keystring = str( key )
            scenefile.write( keystring.rjust(4)+scene + '\n' )
        scenefile.close()
    
        modelfile = open( filename+".mod", 'w' )
        for key,model in models.iteritems():
            keystring = str( key )
            modelfile.write( keystring.rjust(4)+model + '\n' )
        modelfile.close()
    
        settingsfile = open( filename+".set", 'w' )
        for key,setting in settings.iteritems():
            keystring = str( key )
            settingName, selection, startVal, endVal = setting
            settingsfile.write( "%s%s,%s,%f,%f\n" % (keystring.rjust(4), settingName, selection, startVal, endVal))
        settingsfile.close()
        cmd.save( filename+".pse" )
            
        print "Wrote files", filename+".txt,",filename+".frm,",filename+".pse, and",filename+".act"
    
    def writeKeyViewFile( filename ):
        global views
        global frames
        global scenes
        global actions
        global settings
        global models
        global fades
     
        keyviewfile = open( filename + ".key", 'w' )
        for i,view in enumerate(views):
            keyviewfile.write( "VIEW: %4d " % i )
            for v in view:
                keyviewfile.write( str(v) + " " )
            keyviewfile.write('\n')
            keyviewfile.write( "FRAMES: %d\n" % frames[i] )
            if i in actions:
                keyviewfile.write( "ACTIONS: %s\n" % actions[i] )
            if i in scenes:
                keyviewfile.write( "SCENES: %s\n" % scenes[i] )
            if i in models:
                keyviewfile.write( "MODELS: %s\n" % models[i] )
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                keyviewfile.write( "SETTINGS: %s, %s, %f, %f\n" % (settingName, selection, startVal, endVal))
            if i in fades:
                startVisSelection, endVisSelection, sticksOnly = fades[i]
                keyviewfile.write( "FADES: %s, %s, %d\n" % (startVisSelection, endVisSelection, sticksOnly))            
            keyviewfile.write("\n")
        keyviewfile.close()
        cmd.save( filename + ".pse" )
        print "Wrote files " , filename + ".key", filename + ".pse"
    
    def readKeyViewFile( filename ):
        global views
        global frames
        global scenes
        global actions
        global settings
        global models
        global fades
        global scene_counter
    
        views = []
        frames = []
        actions = {}
        scenes = {}
        models = {}
        settings = {}
        fades = {}
        scene_counter = 0
    
        if filename.endswith(".key"): filename = filename[:-4]
        keyviewfile = open( filename + ".key", 'r' )
        viewstrings = keyviewfile.readlines()
        keyviewfile.close()
        viewindex = 0
        for line in viewstrings:
            if line.startswith("VIEW: "):
                viewindex = int( line[6:10] )
                vals = line[10:].split()
                view = [float(v) for v in vals]
                views.append( view )
            if line.startswith("FRAMES: "):
                frames.append( int( line[8:] ) )
            if line.startswith("ACTIONS: "):
                actions[ viewindex ] = line[9:-1]
            if line.startswith("SCENES: "):
                scenes[ viewindex ] = line[8:-1]
                scene_counter += 1
            if line.startswith("MODELS: "):
                models[ viewindex ] = line[8:-1]
            if line.startswith("SETTINGS: "):
                settingName,selection,startVal,endVal = line[10:-1].split(',')
                settings[ viewindex ] = [settingName,selection,float(startVal),float(endVal)]
            if line.startswith( "FADES: " ):
                startVisSelection, endVisSelection, sticksOnly = line[7:-1].split(',')
                fades[ viewindex ] = [startVisSelection, endVisSelection, int(sticksOnly) ]
        cur_view = views[0]
        cur_index = 0
        show_cur()
    
    def set_frames_current( arg=50 ):
        global frames
        global cur_index
        frames[cur_index] = int(arg)
    
    def list_frames():
        global frames
        global views
        global actions
        global models
        global settings
        
        f=0
        for i,view in enumerate(views[:-1]):
            if i in actions:
                a = actions[i]
            else:
                a = ""
            if i in models:
                m = "States: " + models[i]
            else:
                m = ""
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                s = "Settings: %s %s %f %f" % (settingName, selection, startVal, endVal)
            else:
                s = ""
            print "View",i,"to",i+1,"Frames ",f,"to",f+frames[i],a,m,s
            f += frames[i]
    
    def add_action_current( cmd ):
        global cur_index
        global actions
        actions[cur_index] = cmd[1:-1] #strip off quotes
    
    def append_action_current( cmd ):
        global cur_index
        global actions
        actions[cur_index] += ";" + cmd[1:-1]
    
    def clear_action_current():
        global actions
        global cur_index
        del actions[cur_index]
    
    def list_actions():
        global actions
        for i,a in actions.iteritems():
            print i,a
    
    def morph_models( start_model, end_model ):
        global cur_index
        global frames
        global models
        models[cur_index] = "%s -%s" % (start_model, end_model)
        frames[cur_index] = abs(int(end_model) - int(start_model)) + 1
    
    def interpolate_settings( setting, selection, startval, endval ):
        global cur_index
        global settings
        settingEntry = [setting, selection, float(startval), float(endval)]
        settings[cur_index] = settingEntry 
    
    def crossfade( startVisSelection, endVisSelection, sticksOnly = 1 ):
    #cross fade the specified objects, vary stick transparency only if stickOnly=1
        global cur_index
        global fades
        fades[cur_index] = [startVisSelection, endVisSelection, int(sticksOnly) ]
    
    def setup_view( index ):
        for i in range( int(index)+1 ):
            if i in actions:
                print "Executing %s from actions %d" % (actions[i],i)
                cmd.do( actions[i] )
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                action = "set %s, %f, %s;" % (settingName, endVal, selection)
                print "Executing %s from settings %d" % (action,i)
                cmd.do( action )
            if i in fades:
                startVisSelection, endVisSelection, sticksOnly = fades[i]
                action = "set stick_transparency, 0, %s; set stick_transparency, 1, %s;" % (endVisSelection, startVisSelection)
                print "Executing %s from fades %d" % (action, i)
                cmd.do( action )
    
    def show_transition(start_index=0, end_index=0):
        #show the transition from the current view to the next view
        global frames
        global views
        global cur_index
        global actions
        global models
        if start_index == 0 and end_index == 0:
            if cur_index >= len(views)-1:
                print "Current view is last in sequence."
                return
            start_index=cur_index
            end_index=cur_index+1
        else:
            start_index = int(start_index)
            end_index = int(end_index)
        ftot = 0
        setcommand = ""
        i = start_index
        for nframes in frames[start_index:end_index]:
            #ftot += nframes
            if i in models:
                setcommand += " " + models[i] + " "
            else:
                setcommand += " 1 x%i" % nframes
            i += 1
            
    #    cmd.mset("1 x%i" % ftot)
        cmd.mset( setcommand )
        start_frame = 1
        #first do all actions that happen up to this point to make sure
        #things look the way they should
        first_action = ""
        for i in range( start_index ):
            if i in actions:
                first_action += actions[i] + ';'
                #print "Executing %s from actions %d" % (actions[i],i)
                #cmd.do( actions[i] )
            if i in settings:
                settingName, selection, startVal, endVal = settings[i]
                action = "set %s, %f, %s;" % (settingName, endVal, selection)
                first_action += action
                #print "Executing %s from settings %d" % (action,i)
                #cmd.do( action )
            if i in fades:
                startVisSelection, endVisSelection, sticksOnly = fades[i]
                action = "set stick_transparency, 0, %s; set stick_transparency, 1, %s;" % (endVisSelection, startVisSelection)
                first_action += action
                #print "Executing %s from fades %d" % (action, i)
                #cmd.do( action )
        for i in range( start_index, end_index ):
            if i in settings:
                movs.animate_transition( views[i], views[i+1], frames[i], start_frame, settings[i] )
            elif i in fades:
                movs.animate_transition( views[i], views[i+1], frames[i], start_frame, fades[i] )
            else:
                movs.animate_transition( views[i], views[i+1], frames[i], start_frame )
            #add an action
            if start_frame == 1:
                mdo_cmd = first_action
                if i in actions:
                    mdo_cmd += actions[i]+";"
                #mdo_cmd += "set_view("+str(views[i])+")"
                print mdo_cmd
                cmd.mdo(start_frame, mdo_cmd)
            elif i in actions:
                mdo_cmd = actions[i]+";set_view("+str(views[i])+")"
                cmd.mdo(start_frame, mdo_cmd)
                #print mdo_cmd
            start_frame += frames[i]
        cmd.frame(1)
        cmd.mplay()
    
    def make_all():
        #make the whole movie
        global views
        global frames
        global models
        
        #first get total number of frames
        ftot = 0
        setcommand = ""
        for i,nframes in enumerate(frames[:-1]):
            ftot += nframes
            if i in models:
                setcommand += " " + models[i] + " "
            else:
                setcommand += " 1 x%i" % nframes
    
        #initialize movie
        #cmd.mset("1 x%i" % ftot)
        #cmd.mset("1 x50 1 -30 30 x20")
        cmd.mset( setcommand )
    
        #loop through views
        start_view = views[0][:]
        first_frame = 1
        for i,view in enumerate(views[1:]):
            end_view = view[:]
            if i in settings:
                movs.animate_transition( start_view, end_view, frames[i], first_frame, settings[i] )
            elif i in fades:
                movs.animate_transition( start_view, end_view, frames[i], first_frame, fades[i] )
            else:
                movs.animate_transition( start_view, end_view, frames[i], first_frame )
            #add an action
            if i in actions:
                mdo_cmd = actions[i]#+";set_view ("+str( views[i] )+")"
                print mdo_cmd
                cmd.mdo(first_frame, mdo_cmd)
            first_frame += frames[i]
            start_view = end_view[:]
        cmd.frame(1)
    
    ## views = readViews( "viewfile.txt" )
    ## frames = readFrames( "viewfile.frm" )
    ## actions = readActions( "viewfile.act" )
    ##print "Length ",len(views)
    #for v in views:
     #   print v
    #cur_view = views[0]
    views = []
    frames = []
    models = {}
    actions = {}
    scenes = {}
    settings = {}
    fades = {}
    scene_counter = 0
    cur_index = -1
    cmd.set( "scene_animation_duration","0" )
    #cmd.set_view( cur_view )
    
    cmd.extend("sn", show_next )
    cmd.extend("sp", show_prev )
    cmd.extend("sc", show_cur )
    cmd.extend("sinsert", insert_current )
    cmd.extend("sdelete", delete_current )
    cmd.extend("sreplace", replace_current )
    cmd.extend("sappend", append_view )
    cmd.extend("sscene", insert_scene )
    cmd.extend("sgo", go_to_view )
    cmd.extend("sreadold", read_all )
    cmd.extend("swriteold", write_views )
    cmd.extend("slist", list_frames )
    cmd.extend("ssetf", set_frames_current )
    cmd.extend("sshow", show_transition )
    cmd.extend("srecord", make_all )
    cmd.extend("saction", add_action_current )
    cmd.extend("sdelaction", clear_action_current )
    cmd.extend("sdumpactions", list_actions )
    cmd.extend("sappendaction", append_action_current )
    cmd.extend("smorph", morph_models )
    cmd.extend("sinterpsetting", interpolate_settings )
    cmd.extend("sdeletesetting", delete_settings )
    cmd.extend("scrossfade", crossfade )
    cmd.extend("swrite", writeKeyViewFile )
    cmd.extend("sread", readKeyViewFile )
    cmd.extend("ssetupview", setup_view )
    

### movs.py

    Math and animation routines for slerpy
    
    
    ##################################################################################
    #movs.py - Math and animation routines for slerpy
    #Copyright (C) 2006 Joel Bard
    #
    #This program is free software; you can redistribute it and/or
    #modify it under the terms of the GNU General Public License
    #as published by the Free Software Foundation; either version 2
    #of the License, or (at your option) any later version.
    #
    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.
    #
    #You should have received a copy of the GNU General Public License
    #along with this program; if not, write to the Free Software
    #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
    #################################################################################
    
    from pymol import cmd,stored
    from math import *
    
    def rmat2quat( M ):
        #M is a list of 9 values being the elements of the rotation matrix in row order
        T = M[0] + M[4] + M[8] + 1
        print "Trace ",T
        if T>0:
            S = 0.5 / sqrt(T)
            W = 0.25/S
            X = (M[7] - M[5])*S
            Y = (M[2] - M[6])*S
            Z = (M[3] - M[1])*S
        elif M[0] > M[4] and M[0] > M[8]:
            S = sqrt( 1.0 + M[0] - M[4] - M[8]) * 2
            X = 0.25 * S
            Y = (M[1] + M[3])/S
            Z = (M[2] + M[6])/S
            W = (M[5] - M[7])/S
        elif M[4] > M[8]:
            S = sqrt( 1.0 + M[4] - M[0] - M[8] ) * 2
            X = (M[1] + M[3])/S
            Y = 0.25 * S
            Z = (M[5] + M[7])/S
            W = (M[2] - M[6])/S
        else:
            S = sqrt( 1.0 + M[8] - M[0] - M[4]) * 2
            X = (M[2] + M[6])/S
            Y = (M[5] + M[7])/S
            Z = 0.25 * S
            W = (M[1] - M[3])/S
        return [X,Y,Z,W]
    
    def quat2rmat( Q ):
        #Q is a list of 4 values being the quaternion X Y Z W
        X=Q[0]
        Y=Q[1]
        Z=Q[2]
        W=Q[3]
        xx = X*X
        xy = X*Y
        xz = X*Z
        xw = X*W
        yy = Y*Y
        yz = Y*Z
        yw = Y*W
        zz = Z*Z
        zw = Z*W
    
        M= [1.0]*9
        M[0] = 1 - 2 * ( yy + zz )
        M[1] = 2 * ( xy - zw )
        M[2] = 2 * ( xz + yw )
        M[3] = 2 * ( xy + zw )
        M[4] = 1 - 2 * ( xx + zz )
        M[5] = 2 * ( yz - xw )
        M[6] = 2 * ( xz - yw )
        M[7] = 2 * ( yz + xw )
        M[8] = 1 - 2 * ( xx + yy )
        return M
    
    def quatconj( Q ):
        return [-Q[0],-Q[1],-Q[2],Q[3]]
    
    def quatmag( Q ):
        s = 0.0
        QC = quatconj(Q)
        for x in range(4):
            s += Q[x]*Q[x]
        print s
        return sqrt(s)
    
    def quatnorm( Q ):
        m = quatmag( Q )
        return [q/m for q in Q]
    
    def quatdotprod( q1, q2 ):
        dp = 0
        for i in range(4):
            dp += q1[i]*q2[i]
        return dp
    
    def vectnorm( V ):
        mag = 0.0
        for x in V:
            mag += x*x
        mag = sqrt(mag)
        return [x/mag for x in V]
    
    def quat2axisangle( Q ):
        #returns list where 0..2 are rot axis and 3 is angle
        qn = quatnorm( Q )
        cos_a = Q[3]
        angle = acos( cos_a ) * 2
        sin_a = sqrt( 1.0 - cos_a * cos_a )
        if fabs( sin_a ) < 0.000005:
            sin_a = 1
        ax_an = [ q/sin_a for q in Q[0:3] ]
        ax_an.append( angle )
        return ax_an
    
    def axisangle2quat( ax_an ):
        #ax_an is a list with axis coordinates followed by rotation angle
        axn = vectnorm( ax_an[:3] )
        angle = ax_an[3]
        sin_a = sin( angle / 2 )
        cos_a = cos( angle / 2 )
        Q = [ x * sin_a for x in axn ]
        Q.append( cos_a )
        return Q
        
    def rmat2axisangle( M ):
        q = rmat2quat( M )
        return quat2axisangle( q )
    
    def axisangle2rmat( a ):
        q = axisangle2quat( a )
        return quat2rmat( q )
    
    def animate_transition( start_view, end_view, nframes, first_frame, settings = [] ):
        #print "Views"
        #print start_view,'\n',end_view
    
        cview = start_view[:]
        cmd.set_view( start_view )
    
        #get initial and final quaternions for interpolation
        #print start_view[0:9]
        #get quaternions
        qstart = rmat2quat( start_view[0:9] )
        qend = rmat2quat( end_view[0:9] )
        
        #test for long way vs. short way
        if quatdotprod( qstart,qend ) < 0:
            qstart = [-q for q in qstart]
    
        axan_start = quat2axisangle( qstart )
        axan_end = quat2axisangle( qend )
            
        axan_cur = axan_start[:]
        frame_start = first_frame
        frame_end = frame_start + nframes
        doFade = 0
        doSetting = 0
        if len( settings ) == 4:
            settingName, selection, startVal, endVal = settings
            settingStep = (endVal-startVal)/float(nframes)
            print "Setting step ", settingStep
            doSetting = 1
        elif len( settings ) == 3:
            startVisSelection, endVisSelection, sticksOnly = settings
            settingStep = 1.0/float(nframes)
            doFade = 1
        for f in range( frame_start , frame_end):
            #get rotmat
            #using angle axis
    
            for i in range(4):
                axan_cur[i] = axan_cur[i] + (axan_end[i]-axan_start[i])/nframes
            newmat = axisangle2rmat( axan_cur )
            #print cview
            for i in range(9):
                cview[i] = newmat[i]
    
            mdo_cmd = "set_view (["
            for i in range(18):
                if i>8:
                    cview[i] = cview[i]+(end_view[i]-start_view[i])/nframes
                mdo_cmd += "%12.7f,"% cview[i]
            mdo_cmd = mdo_cmd[:-1]+"])"
            if doSetting:       
                val = float(f-frame_start)*settingStep + startVal 
                print val;
                mdo_cmd += "; set %s, %f, %s" % (settingName, val, selection)
                print mdo_cmd;
            #print "mdo ", mdo_cmd
            if doFade:
                val = float(f-frame_start)*settingStep
                otherVal = 1.0 - val
                mdo_cmd += "; set stick_transparency, %f, %s; set stick_transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
                if not sticksOnly:
                    #comment out surface transparency interpolation due to problem with transparent sticks in front of 
                    #transparent surfaces (get black holes)
                   # mdo_cmd += "; set transparency, %f, %s; set transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
                    mdo_cmd += "; set cartoon_transparency, %f, %s; set cartoon_transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
            cmd.mdo(f,mdo_cmd)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Slerpy&oldid=8080](https://pymolwiki.org/index.php?title=Slerpy&oldid=8080)"


---

## Supercell

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**[Included in psico](/index.php/Psico "Psico")**  
This command or function is available from [psico](/index.php/Psico "Psico"), which is a PyMOL extension.   
---  
Module  | [psico.xtal](http://pymol.org/psicoredirect.php?psico.xtal)  
  
[![](/images/4/48/SupercellExample.png)](/index.php/File:SupercellExample.png)

[](/index.php/File:SupercellExample.png "Enlarge")

Example with 2 unit cells in c-direction, created with: supercell 1,1,2,2x19

**supercell** can display multiple copies of the unit cell. Can also fill the unit cell (and its copies) with symmetry mates. 

See [thread on pymol-users mailing list](http://sourceforge.net/mailarchive/forum.php?thread_name=l2vdcf611bd1004140816zeca28714mf76b9f72008099ab%40mail.gmail.com&forum_name=pymol-users). 

Requires [numpy](http://numpy.scipy.org). 

## Example
    
    
    # import psico.xtal, or run the standalone script:
    run https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/xtal.py
    
    fetch 2x19, async=0
    supercell 2,1,1, 2x19, green
    supercell 1,1,2, 2x19, orange, name=super2
    

## See Also

  * [SuperSym](/index.php/SuperSym "SuperSym")
  * [symexp](/index.php/Symexp "Symexp")



Retrieved from "[https://pymolwiki.org/index.php?title=Supercell&oldid=13313](https://pymolwiki.org/index.php?title=Supercell&oldid=13313)"


---

## SuperSym

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/SuperSymPlugin.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/SuperSymPlugin.py)  
Author(s)  | [Stuart Ballard](/index.php?title=User:Srballard&action=edit&redlink=1 "User:Srballard \(page does not exist\)")  
License  |   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
[![](/images/3/34/SuperSymExample.png)](/index.php/File:SuperSymExample.png)

[](/index.php/File:SuperSymExample.png "Enlarge")

Symmetry partners for 1hpv showing 6-1 screw axis

[![](/images/d/da/SuperSymExample2.png)](/index.php/File:SuperSymExample2.png)

[](/index.php/File:SuperSymExample2.png "Enlarge")

Full cell of symmetry partners with symmetry axes displayed

SuperSym is a PyMOL plugin providing a large number of tools for visualization of space groups; unit cells; and symmetry axes, operators, and partners. 

The original source code is available from <https://sourceforge.net/projects/supersym/>

## Contents

  * 1 Dependencies
  * 2 Bugs
  * 3 Acknowledgments
  * 4 Installing SuperSym
  * 5 Feedback
  * 6 The Menu



## Dependencies

  * [cctbx](/index.php/CCTBX "CCTBX")
  * numpy



**PyMOL, cctbx and numpy must all be compiled with the same Python distribution!** See [CCTBX](/index.php/CCTBX "CCTBX"). 

This plugin was developed in 2010 for PyMOL version 1.2r1 and has not been updated since. 

## Bugs

Symmetry axes are not defined for all space groups, and do not display properly for some. 

## Acknowledgments

Primary coding and development was done by Stuart Ballard. All comments, questions, and issues should be directed to him at srballard@wisc.edu. 

Code for unit cell and symmetry axis building is borrowed from scripts created by Robert Campbell and Ralf W. Grosse-Kunstleve, available at <http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/>. Some of this code has been modified for use in SuperSym. 

[FindSurfaceResidues](/index.php/FindSurfaceResidues "FindSurfaceResidues") is utilized for some of SuperSym's graphics generation, with some modifications. 

## Installing SuperSym

To install SuperSym v1.2, download SuperSymPlugin12.py from <https://sourceforge.net/projects/supersym/>. In PyMOL, go to: 

  * Plugin > Manage Plugins > Install...



A file selector dialog will appear. Select SuperSymPlugin12.py. PyMOL will direct you to restart, and upon doing so SuperSym will be accessible through the Plugin menu. 

To use functions of SuperSym directly, without creating a drop-down menu, use the run command in PyMOL on SuperSymPlugin12.py. 

Note: previous errors resulting from incorrect naming of the plugin file have been resolved in v1.2. 

## Feedback

Please post any comments, complaints, bug fix requests, useful tricks, or cool adaptations of SuperSym here. 

## The Menu

  * **Default Symmetry Partner Set**
    * See **Build Symmetry Partners > Cell [0,0,0] (default)**
  * **Draw Unit Cell**
    * Creates a cgo object with unit cell axes as cylinders. This functions similarly to _show cell_ , but the cell axes are cylinders instead of lines, allowing for printing.
  * **Build Symmetry Partners >**
    * All options in this submenu generate sets of symmetry partners
    * **Cell [0,0,0] (default)**
      * Generates a suite of symmetry partners for a given object for the default unit cell, which is lattice position [0,0,0]
    * **Cell [x,y,z] (custom)**
      * Generates a suite of symmetry partners for a given object for a lattice position which you specify
    * **2x2x2 Block**
      * Generates 8 sets of symmetry partners for a given object, filling lattice positions [0,0,0] through [1,1,1]
    * **3x3x3 Block**
      * Generates 27 sets of symmetry partners for a given object, filling lattice positions [-1,-1,-1] through [1,1,1]. This option may take a long time to execute
    * **By Partner**
      * Generates only those symmetry partners which the user specifies by their defining symmetry operators
  * **Coloring >**
    * **Default Rainbow**
      * Colors all symmetry objects with a specified by their symmetry operations automatically
    * **Select color for each operation**
      * Select symmetry partners to color by their defining symmetry operation and select the color for each
    * **Select one color for custom set of operations**
      * Select a set of symmetry partners defined by symmetry operations and select one color for all of them
  * **Graphics >**
    * **Lines**
      * Convenience function to display symmetry partners as lines
    * **Ribbon**
    * Convenience function to display symmetry partners as ribbons
    * **Cartoon**
      * Convenience function to display symmetry partners as cartoons
    * **Sphere Surface (best for printing)**
      * Uses the findSurfaceResidues function and shows surface residues as spheres. If printing, this option saves at least 60% of materials relative to regular surfaces, with minimal loss in resolution
    * **Surface (high load render)**
      * Displays symmetry partners as surfaces. This option may take a very long time to execute
  * **Symmetry Axes >**
    * **Build Axes**
      * Builds all symmetry axes for the given object. This functionality will be customizable and extended in future versions
  * **Move symmetry partners**
    * Merely displays instructions for using built in hotkeys to move symmetry partners
  * **About**
    * Developer info
  * **Help**
    * Reference to this page



Retrieved from "[https://pymolwiki.org/index.php?title=SuperSym&oldid=13315](https://pymolwiki.org/index.php?title=SuperSym&oldid=13315)"


---

## Symmetry Axis

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.
    
    
    ###########################################################
    #
    #  Pymol script copyright Matthew O'Meara and Xavier Ambroggio 2007
    #
    #  Last updated Nov 29, 2007
    #
    #  Draw an axis given a point and a direction.  Optionally give color,
    #  length and width.
    #
    #  Usage: from the pymol menu click file->run...  then find this file.
    #  Then at the prompt, type
    #
    #           draw_axis x,y,z, i,k,j
    #
    #  where (x,y,z) is the point and (i,k,j) is the direction
    #
    #
    #  Also one can run the script as follows
    #
    #
    #           draw_axis x,y,z, i,k,j, length r,g,b, width,
    #
    #  where (r,g,b) is the color for the line in red, green, colors,
    #  width is the thickness of the line and length is the length of the
    #  line.
    #
    #
    #  For a fun example of the script run from the command line after the
    #  script is loaded
    #
    #           draw_axis_example
    #
    #
    
    from pymol.cgo import *    # get constants
    from pymol import cmd
    
    import math
    
    class Counter:
       def __init__(self):
           self.state = 1
    counter = Counter()
    
    def draw_axis(x=None, y=None, z=None, i=None, j=None, k=None, length=20.0, r=1.0, g=1.0, b=1.0, width=1.0 ):
       if x == None or y == None or z == None or i == None or j == None or k== None :
           print 'Usage: draw_axis x,y,z, i,k,j, length, r,g,b, width'
           print 'draw a line centered at (x,y,z) with the direction vector (i,j,k)'
           print 'length, color (r,g,b), and width arguments are optional'
    #        print 'For a fun example of the command, run draw_axis_example'
       else :
           x,y,z = float(x), float(y), float(z)
           i,j,k = float(i), float(j), float(k)
           r,g,b = float(r), float(g), float(b)
           width = float(width)
           length = float(length) / 2.0
    
           x1,y1,z1 = (x+i*length,y+j*length,z+k*length)
           x2,y2,z2 = (x-i*length,y-j*length,z-k*length)
    
           obj = [
               LINEWIDTH, width,
               BEGIN, LINES,
    
               COLOR,   r,  g,  b,
               VERTEX, x1, y1, z1,
               VERTEX, x2, y2, z2,
    
               END
               ]
    
           cmd.load_cgo(obj,'axis'+str(counter.state))
           counter.state += 1
    
    cmd.extend("draw_axis", draw_axis)
    
    # a simple example
    #draw_line(x=18.232,  y=17.150,  z=9.488,
    #          i=-.226639,j=0.708772,k=-.668039,
    #          r=1,       b=1,       g=1,
    #          width=1,   length=1)
    
    
    
    
    # a more complex example
    
    #import random
    #def example1(n, f):
    #    """draw a gradient field with n segments with the function f(x,y,z)=(i,j,k)"""
    #    for i in range(n):
    #        scale = 4
    #        x,y,z = [random.random()*scale for i in range(3)]
    #        i,j,k = f(x,y,z)
    
    #        draw_axis(x,y,z,i,j,k,abs(i),abs(j),abs(k))
    
    
    #def f(x,y,z):
    #    return (2*x,pow(z,2)+x,y-z)
    
    #cmd.extend("draw_axis_example", lambda :example1(1000,f))
    

Retrieved from "[https://pymolwiki.org/index.php?title=Symmetry_Axis&oldid=6392](https://pymolwiki.org/index.php?title=Symmetry_Axis&oldid=6392)"


---

## TransformSelectionByCameraView

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This little script was posted to the PyMol list. It will orient the molucule like cmd.orient() does, but does so by the camera view. 

<https://sourceforge.net/p/pymol/mailman/message/10097639/>
    
    
    # transform selection coordinates by the camera view
    #
    # The script answers this:
    #   Thanks!
    #   But translate[x,y,z] only translate the molecule.
    #   What I want  is to put longest length of molecule in the X axes, the 
    #   second Y axes, the third z axes.
    #   Just like what orient command does which change the view of camera but 
    #   not the coordinates.
    #   Now I want the coordinates also change after orient it.
    #
    cv=list(cmd.get_view())
    
    cmd.transform_selection("all", \
      cv[0:3]+[0.0]+ \
      cv[3:6]+[0.0]+ \
      cv[6:9]+[0.0]+ \
      cv[12:15]+[1.0], transpose=1)
    
    cmd.reset()
    

Retrieved from "[https://pymolwiki.org/index.php?title=TransformSelectionByCameraView&oldid=12248](https://pymolwiki.org/index.php?title=TransformSelectionByCameraView&oldid=12248)"


---

