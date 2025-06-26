# Category: CGO

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

## CGO Shapes

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A Compiled Graphics Object (CGO) is a Python list of floating point numbers, which is compiled by PyMOL into an geometric representation and associated with a given state. CGOs can be created in the PyMOL scene with the [Load_CGO](/index.php/Load_CGO "Load CGO") command. 

### CGO Commands

Based on an old version of the User Manual: 
    
    
    Lowercase names below should be replaced with floating-point numbers. Generally, the TRIANGLE primitive should only be used only as a last restore since it is much less effective to render than using a series of vertices with a BEGIN/END group.
    
    BEGIN, { POINTS | LINES | LINE_LOOP | LINE_STRIP | TRIANGLES | TRIANGLE_STRIP | TRIANGLE_FAN },
    
    VERTEX, x,  y,  z,
    
    COLOR,  red, green, blue, 
    
    NORMAL, normal-x, normal-y,  normal-z, 
    
    END,
    
    LINEWIDTH, line-width, 
    
    WIDTHSCALE, width-scale,   # for ray-tracing
    
    SPHERE, x, y, z,  radius    # uses the current color
    
    
    [[Category:CGO]]
    
    CYLINDER, x1, y1, z1, x2, y2, z2, radius,
              red1, green1, blue1, red2, green2, blue2,
    
    TRIANGLE,  x1, y1, z1, 
               x2, y2, z2,
               x3, y3, z3,
               normal-x1, normal-y1, normal-z1,
               normal-x2, normal-y2, normal-z2,
               normal-x3, normal-y3, normal-z3,
               red1, green1, blue1,          
               red2, green2, blue2,          
               red3, green3, blue3,          
    
    CONE,      x1, y1, z1,
               x2, y2, z2,
               r1, r2,
               red1, green1, blue1,          
               red2, green2, blue2,          
               cap1, cap2   # should the ends be solid (1) or open (0)?           
    

  


### EXAMPLES
    
    
    from pymol import cgo
    
    # Unit sphere at the origin
    plain_sphere = [cgo.SPHERE, 0, 0, 0, 1]
    cmd.load_cgo(plain_sphere, "plain_sphere")
    
    # Red unit sphere at (3, 4, 5)
    red_sphere = [cgo.COLOR, 1, 0, 0, cgo.SPHERE, 3, 4, 5, 1]
    cmd.load_cgo(red_sphere, "red_sphere")
    
    # Narrow yellow cylinder
    pos0 = [-1, -1, -1]
    pos1 = [-1, -1, 1]
    r = 0.1
    yellow = [1, 1, 0]
    yellow_cyl = [cgo.CYLINDER, *pos0, *pos1, r, *yellow, *yellow]
    cmd.load_cgo(yellow_cyl, "yellow_cylinder")
    
    # Magenta and green cone, open on one endray
    pos0 = [2, 2, 2]
    pos1 = [3, 3, 3]
    radii = [1, 0.5]
    color0 = [1, 0, 1]
    color1 = [0, 1, 0]
    caps = [0, 1]
    cone = [cgo.CONE] + pos0 + pos1 + radii + color0 + color1 + caps
    cmd.load_cgo(cone, "cone")
    
    # adjust the camera position
    set_view (\
         0.022470249,    0.398223877,   -0.917012811,\
         0.188309744,    0.899140060,    0.395077318,\
         0.981853366,   -0.181560174,   -0.054785311,\
         0.000002960,   -0.000002829,  -28.313688278,\
         1.387665749,    1.797374249,    2.392468214,\
        21.653684616,   34.973686218,  -20.000000000 )
    
    # force a black background and save a PNG image
    set ray_opaque_background, 1
    png cgo_examples.png, ray=1
    

  * [![Figure produced by the above example code.](/images/e/e6/Cgo_examples.jpg)](/index.php/File:Cgo_examples.jpg "Figure produced by the above example code.")

Figure produced by the above example code. 




  


### SEE ALSO

  * [Load_CGO](/index.php/Load_CGO "Load CGO")
  * [cgo_transparency](/index.php/Cgo_transparency "Cgo transparency")



Retrieved from "[https://pymolwiki.org/index.php?title=CGO_Shapes&oldid=13484](https://pymolwiki.org/index.php?title=CGO_Shapes&oldid=13484)"


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

## Cgo transparency

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This is an object-level setting to control transparency of Compiled Graphics Objects. 

Usage: 
    
    
    # PyMOL command
    set cgo_transparency, value [, object]
    
    # Python API
    cmd.set("cgo_transparency", value [, object])
    

Arguments 

  * value (float) number between 0 and 1 (with 0 meaning the object is opaque, and 1 == invisible).
  * object (str) name of the object to which to apply the setting (default: all CGO objects)



Retrieved from "[https://pymolwiki.org/index.php?title=Cgo_transparency&oldid=13483](https://pymolwiki.org/index.php?title=Cgo_transparency&oldid=13483)"


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

## CGOCylinder

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
  * 2 Examples
  * 3 Usage
  * 4 References



# Overview

[CGO](/index.php/Category:CGO "Category:CGO") cylinders are compiled graphics cylinder objects. If you set the cylinder height to something very small, you can mimic a circle (see examples). 

  * [![Example cylinder and label. The label was added as a label on a Pseudoatom.](/images/1/16/CylCGO.png)](/index.php/File:CylCGO.png "Example cylinder and label. The label was added as a label on a Pseudoatom.")

Example cylinder and label. The label was added as a label on a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). 

  * [![Cylinder with a very thin cyilnder, here 0.2 Angstroms tall to mimic a circle. The label was added as a label on a Pseudoatom.](/images/3/35/Cylsphere.png)](/index.php/File:Cylsphere.png "Cylinder with a very thin cyilnder, here 0.2 Angstroms tall to mimic a circle. The label was added as a label on a Pseudoatom.")

Cylinder with a very thin cyilnder, here 0.2 Angstroms tall to mimic a circle. The label was added as a label on a [Pseudoatom](/index.php/Pseudoatom "Pseudoatom"). 




# Examples
    
    
    # Cylinder #1 on left.
    x1,y1,z1 = 10, 0, 0 # start point
    r1,g1,b1 = 1,0,0 # color (red)
    x2,y2,z2 = 0.1, 0, 0 # end point
    r2,g2,b2 = 1,1,0 # color (yellow)
    radius = 10
    cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder1" )
    
    # Cylinder #2 (the circle)
    x1,y1,z1 = -0.1, 0, 0 # start point
    r1,g1,b1 = 1,0,0 # color (red)
    x2,y2,z2 = 0.1, 0, 0 # end point
    r2,g2,b2 = 1,1,0 # color (yellow)
    radius = 10
    cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder2" )
    
    # to make the cylinder transparent, use
    cmd.load_cgo( [ 25.0, 0.25, 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder1" )
    

# Usage
    
    
    cmd.load_cgo( [ a1, a2, cylSpec, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], objName )
    

The parameters are: 

**a1, a2**

    

    a1=25 means 'ALPHA' (if you do **from pymol.cgo import *** you can then use CGO keywords, like ALPHA, CYLINDER, etc. So, if a1="25.0" then set a2 to some value in the range 0-1.

**cylSpec**

    

    Cylinder specification: set this to 9.0 or if you did **from pymol.cgo import *** you can then use "CYLINDER", here.

**x1, y1, z1**

    

    coordinates for one end of the cylinder

**x2, y2, z2**

    

    coordinates for one end of the cylinder

**radius**

    

    radius of the cylinder

**r1, b1, g1**

    

    color one, (a triplet of values in the range 0-1

**r2, b2, g2**

    

    color two, (a triplet of values in the range 0-1; this allows a gradient

**objName**

    

    Name of cylinder object in PyMOL

# References

PyMOL list. Thanks to Tsjerk. 

Retrieved from "[https://pymolwiki.org/index.php?title=CGOCylinder&oldid=7043](https://pymolwiki.org/index.php?title=CGOCylinder&oldid=7043)"


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

## DrawGridBox

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [scripts/drawgridbox.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/scripts/drawgridbox.py)  
Author(s)  | [Cunliang Geng](/index.php?title=User:Clgeng&action=edit&redlink=1 "User:Clgeng \(page does not exist\)")  
License  | [BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)  
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
This script adds the [drawgridbox](/index.php?title=Drawgridbox&action=edit&redlink=1 "Drawgridbox \(page does not exist\)") command, which draw grid boxes around a selection. 

## Usage
    
    
     drawgridbox [ selection,  [nx, [ny, [nz, [padding, [lw, [r, [g, [b]]]]]]]]]
    

  


## Prameters

  * **selection** = str: atom selection {default: all}
  * **nx** = int: number of grids on axis X {default: 10}
  * **ny** = int: number of grids on axis Y {default: 10}
  * **nz** = int: number of grids on axis Z {default: 10}
  * **padding** = float: default 0.0
  * **lw** = float: line width {default: 2.0}
  * **r** = float: red color component, valid range is [0.0, 1.0] {default 1.0}
  * **g** = float: green color component, valid range is [0.0, 1.0] {default 1.0}
  * **b** = float: blue color component, valid range is [0.0, 1.0] {default 1.0}



  


## Examples
    
    
    load drawgridbox.py
    fetch 1cbh
    show surface
    drawgridbox 1cbh, nx=5, ny=5, nz=5,  lw=1, g=0, b=0
    

[![Drawgridbox.png](/images/5/5a/Drawgridbox.png)](/index.php/File:Drawgridbox.png)

  

    
    
    drawgridbox 1cbh, nx=5, ny=5, nz=5, padding=5, lw=1, r=0, g=0
    

[![Drawgridbox padding.png](/images/3/32/Drawgridbox_padding.png)](/index.php/File:Drawgridbox_padding.png)

Retrieved from "[https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904](https://pymolwiki.org/index.php?title=DrawGridBox&oldid=13904)"


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

## Load CGO

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**cmd.load_cgo()** loads Compiled Graphics Objects into the PyMOL session. CGOs are defined by lists of floating point numbers that define either a series of triangles/vertices and colors, or one of several predefined sequences that describe mathematical shapes. 

### PYMOL API
    
    
    cmd.load_cgo(object, name)
    

### ARGUMENTS

  * **object : list of floats** Coordinate data for the CGO.
  * **name : string** Name of the object (for the internal GUI)



### SEE ALSO

[CGO Shapes](/index.php/CGO_Shapes "CGO Shapes")

Retrieved from "[https://pymolwiki.org/index.php?title=Load_CGO&oldid=13474](https://pymolwiki.org/index.php?title=Load_CGO&oldid=13474)"


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

