# Category: Third Party Software

## Making Movies

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Movie Creation in PyMOL
  * 2 Encoding video files
    * 2.1 Animated GIF
    * 2.2 Avidemux
    * 2.3 Mencoder/Mplayer
      * 2.3.1 Troubleshooting
    * 2.4 Virtual Dub
  * 3 Keywords
    * 3.1 Other Options
  * 4 Movies on a Mac
    * 4.1 QT Pro
    * 4.2 Alternative
    * 4.3 In a Presentation
  * 5 Troubleshooting
    * 5.1 Quicktime



## Movie Creation in PyMOL

This page talks about how to take your export movie frames and turn them into an actual movie. For learning how to make a movie in the first place, see [MovieSchool](/index.php/MovieSchool "MovieSchool"). 

## Encoding video files

### Animated GIF

Quality note: Converting to video generally leads to a loss of quality and inflated file sizes. Gif animation gives a good alternative, retaining the original quality and maintaining small file sizes, however it lacks play controls. With skillful use of frame delays during the generation of the Gif it is possible to create pauses to emphasize certain frames. Ulead Gif Animator works well although there are many programs to choose from. To create files with many frames in Ulead add around 150 files as frames at a time. 

### Avidemux

[Avidemux](http://avidemux.sourceforge.net/) is a great tool to stich your image files together as a movie. It can read the PNG output stack from Pymol and encode a movie using almost all codecs currently available. The program has an easy to use graphical interface, making the conversion simple compared to some of the other options available to you. 

### Mencoder/Mplayer

[Mplayer](http://www.mplayerhq.hu) is an award-winning open source movie player. Mencoder (which comes packaged with Mplayer) is its movie encoder. Mencoder can take in various file formats (png,gif,jpg) and convert them to movies. 

Assuming you have created a lot of .png files and would like to encode a .mpeg, .avi or other video format, a number of solutions are known: 

  * The DiVX encoder using mplayer and mencoder? There's binaries for Unix and Windows. It makes rather nice compression on a 800x600 (probably higher). It doesn't take too long to produce the nicer quailty movies, but much longer than simply


    
    
    mencoder "mf://*.png" -mf type=png:fps=18 -ovc lavc -o output.avi
    

namely something like (this command worked BEST for my case; it's all one line connect it where the backslashes are). 
    
    
    mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts \
    vcodec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \
    "mf://*.png" -mf type=png:fps=18 -o output.avi
    

The mpeg4 codec requires a DivX plugin which is not a part of the default installation on some operation systems. The codec msmpeg4v2 makes movies which are more likely to be playable on standard Windows players and can be used with mencoder e.g. 
    
    
    mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq "mf://*.png" -mf type=png:fps=18 -o output.avi
    

This encoder line is however not optimised (yet), and the codec also produces a bit larger files than the mpeg4 at the same visual quality according to the Mplayer homepage. 

#### Troubleshooting

  1. If your movie shows up in PowerPoint with only the first frame showing, or is just a black square, try adding **-of asf** or **-of avi**.
  2. If the above line complains about the **codec** part, then replace the word **codec** in the command line with **vcodec**.



### Virtual Dub

[Virtual Dub](http://virtualdub.sourceforge.net/) is an open-source, robust piece of software that allows easy creation of .AVI files from image sequences (such as *.PNG created by PyMol). After opening Virtual Dub, select **File** Menu --> **Open Video File** \--> under "files of type" select **image sequence (*.png, etc)**. 

After viewing (and editing) the movie, choose **save as AVI (F7)** from **File** menu. This generates high quality, yet very large, AVI files; in order to down-size the files into formats which fit PowerPoint, free software such as [FormatFactory](http://www.formatoz.com/download.html) (**FF**) is available for download and easy use. 

## Keywords

mencoder, mplayer, movie 

### Other Options

  * For those with Photoshop CS3 installed: Load files using File - Scripts - Load Files into Stack. This loads the files as layers. Turn on the animation window at Window - Animation. You may need to reverse the frames for the correct direction - this is done at the animation window using the drop menu above the frames. You will also want to set a delay on the frames, select all frames using the dropdown menu and and then adjust to 0.2s on one of the frames. Then export the film using File - Export - Render Video. There are plenty of options for format and compression types. I found a quicktime movie with no compression at 10f/s to be right for my needs. The video could get quite big, but it will look as good as pymol. Alternatively you can also create animated GIFs in Photoshop. Once you have the frames of your clip as above, choose File - Save for Web and Devices and then save the GIF. The limit of 200 pictures is a problem in Photoshop.
  * Another good program for converting images into movies of different formats is VideoMach : <http://gromada.com/VideoMach.html>
  * TMPGEnc from <http://www.tmpgenc.net> is very fast, easy to use, and produces very nice ouput (MPEG-2). Unfortunately, it does not handle images larger than 720 x 576 pixels.
  * Adobe Premiere recipe, using Microsoft's MPEG4 V2, 960x720 @ 30 fps, which PowerPoint automatically treats as full-screen (due it's wacky metrics). Using this codec, a recent 24-second movie consumed only 4.5 MB of space, but looks much better than a 640x480 Cinepak-based movie with a file size of around ~40 MB. It definitely pays to use the latest technology.
  * A freeware jiffy to convert png files to an animation is imgcon, which proved to be very useful:<http://www.fmrib.ox.ac.uk/~yongyue/imgcondl.html>
  * The program convert, part of the ImageMagick suite of programs can be used provided the program mpeg2encode is in your path. mpeg2encode can be found here <http://www.mpeg.org/MSSG/>. To make an mpeg file, go


    
    
    convert *.png movie.mpeg
    

If you want to make an animated gif, do 
    
    
    convert *.png movie.gif
    

## Movies on a Mac

### QT Pro

On a Mac, it's easy enough to just fire up QuickTime Player and select 
    
    
    File->Open Image Sequence
    

This will prompt you to select the first png file in a folder and load all other pngs with the same base name in that folder. Then simply export these frames as a movie with any of the available quicktime codecs at a frame rate you like. 

### Alternative

The above requires QT Pro ($25). For free you can use this: [[1]](http://www.likelysoft.com/framed/). 

### In a Presentation

You don't need to convert the PyMol scene to a movie, just save session as a .psw (pymol show file) with the scenes embedded in there. Then from within Power|Point just set up a hyperlink (Insert hyperlink) and point it to the .psw file. In presentation mode when you click on the hyperlink pymol automatically boots and you can scroll thru your scens like a full screen powerpoint presentation. You have to make sure that you have pymol installed on your presentation computer (and you have your psw file) and just click OK when powerpoint warns of the perils of non-microsoft products. the final scene will return to your powerpoint presentation without any effort. 

# Troubleshooting

## Quicktime

Recent attention to transparencies and QuickTime show that one typically wants to set 
    
    
    set ray_opaque_background, off
    

to tell PyMol to stop using Alpha-channel for transparencies and switch to blending. 

If you're having transparency problems this could be a fix. 

Retrieved from "[https://pymolwiki.org/index.php?title=Making_Movies&oldid=8676](https://pymolwiki.org/index.php?title=Making_Movies&oldid=8676)"


---

## O

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

Python programmers, who wish to extract information from O files, may find the OdbParser module useful. The module reads both binary and formatted O files, and is very simple to use. Try it! Simply import the module, and use the get() method: 
    
    
    $ python
    Python 2.3.4 (#1, Mar 10 2006, 06:12:09)
    [GCC 3.4.5 20051201 (Red Hat 3.4.5-2)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
     >>>
     >>> import odbparser
     >>> db = odbparser.get("binary.o")
    

Now you have the O database as a Python dictionary 'db', with the datablock names as keys. For example, print out some information on the di-peptide from Baton: 
    
    
     >>> print db["di_atom_name"]
    ('N', 'CA', 'C', 'O', 'CB', 'N', 'CA', 'C', 'O', 'CB')
     >>> print x["di_residue_type"]
    ('ALA', 'ALA')
    

Go to [OdbParser](http://www.bioxray.dk/~mok/odbparser) for information about download, installation, etc. 

Retrieved from "[https://pymolwiki.org/index.php?title=O&oldid=4408](https://pymolwiki.org/index.php?title=O&oldid=4408)"


---

## PovRay

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL can export input files for [POV-Ray](http://www.povray.org/) with the ".pov" file extension: 
    
    
    set stick_ball
    save input.pov
    

It can also use POV-Ray directly for rendering with the [ray](/index.php/Ray#Renderer "Ray") command: 
    
    
    ray renderer=1
    

Since PyMOL 1.7.4, round stick caps are only exported correctly with [stick_ball](/index.php/Stick_ball "Stick ball")=on. 

## Nice PovRay settings

I typically use the make_pov.py script and "run" it from pymol once to load the function, and then I do _make_pov('povray.inp')_ to create the povray.inp file. Then I edit that file to insert some lines like: 
    
    
    fog {
      distance 10
      fog_type 2
      fog_alt 10.
      fog_offset -160.
      up <0.,1.,.4>
    colour rgbt<1.0, 1.0, 1.0, 0.1>
    turbulence 0.8
    }
    

In this case I'm not really doing depth-cueing but adding fog at the lower background edge (there were two planes defining the background and a surface below the molecule) rising up towards the front upper edge of the scene. 

"fog_type 2" means a "rising fog" along the "up" vector. fog_type 1 is a constant fog. To get pure depth cueing, you would want "up" to be along the <0., 0., 1.> vector (I think!). You'll need to play around with the distance and fog_offset parameters. You wouldn't necessarily want the "turbulence" parameter in there either. Check out "Atmospheric Effects" in the povray documentation for many more details: <http://www.povray.org/documentation/view/201/>

## make_pov.py v1
    
    
    # make_pov.py
    # Do "run make_pov.py" from within pymol and then execute the script
    # with "make_pov('povray.inp')" to create the povray.inp file.
    #
    # written by Robert Campbell 2003
    #
    from pymol import cmd
    
    def make_pov(file):
    	(header,data) = cmd.get_povray()
    	povfile=open(file,'w')
    	povfile.write(header)
    	povfile.write(data)
    	povfile.close()
    

## make_pov.py v2

This is a more extended version of an earlier extension of the version by Robert Campbell. The scene is written in two parts, a .pov file containing all meta data, such as the lights, camera and #defaults, and an include file (.inc) which contains the structure. In this way you have maximum control over your scene without having to edit a huge povray file. You may even want to consider splitting your scene up in separate parts (taken from the same perspective), which you combine in a global .pov file using #include statements. This will give even more control with regards to modifications to the scene. If 'clip' is set to near|far|both, then the corresponding clipping plane(s) is/are included in a CSG difference object. Note that the result may increase the render time significantly unless the scene is simple. 

Once you run **run make_pov.py** , run **make_pov** to execute the script. 

NB. the .pov file contains a commented statement with regards to a povray macro file, which allows transforming scenes and objects from model space to camera space and vice versa. The macro file is given below. 
    
    
    # make_pov.py
    # Do "run make_pov.py" from within pymol and then execute the script
    # with "make_pov('povray.inp')" to create the povray.inp file.
    #                                                                                                   
    # Original script written by Robert Campbell
    # Modified by Tsjerk A. Wassenaar
    #
    
    from pymol import cmd
    
    def make_pov(file, name="PymolObject", meta=True, clip=False ):
            f1, f2 = file, file[:-4] + '.inc'
    
            (header,data) = cmd.get_povray()
            povfile = open(f1,'w')
            if meta: povfile.write(header)
            povview = cmd.get_view()
    
            if clip:
                    objtype = "difference"
                    objclip = ""
                    if clip in ["near","both"]:
                            objclip = objclip + "plane { z, -%f }" % povview[15] 
                    if clip in ["far","both"]:
                            objclip = objclip + "plane { z, -%f }" % povview[16] 
            else:
                    objtype = "object"
                    objclip = ""
    
            povfile.write("""\n
    // Uncomment the following lines if you have the pymolmacro.inc include file and want to use it.
    /*
    #include \"pymolmacro.inc\"
    PYMOL_VIEW( %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f,
                %10.5f, %10.5f, %10.5f )
    */
    
    """ % povview)
    
            povfile.write("""
    #declare %s = union { #include "%s" }
    %s { %s %s }
    """ % (name, f2, objtype, name, objclip ) )
    
            povfile.close()
            povfile = open(f2,'w')
            povfile.write(data)
            povfile.close()
    
    cmd.extend('make_pov',make_pov)
    
    
    
    //
    //  PYMOLMACRO.INC v0.2 
    //
    //  (c)2005 Tsjerk Wassenaar, University of Groningen
    //
    //  This include file for Povray contains
    //  just a few macros which together allow
    //  the conversion between the model space
    //  (cartesian coordinates) and the Pymol
    //  camera space.
    //
    //  With these macros one can easily combine
    //  a Pymol scene with objects defined in the
    //  coordinate space of the original
    //  structure file.
    //
    //  The input consists of the output of the
    //  get_view() command in Pymol. This output
    //  consists of 18 floating point numbers
    //  defining a rotation matrix and shift
    //  vectors for the origin of rotation and
    //  for the camera position.
    //
    //  The macro PYMOL_VIEW loads a
    //  view obtained from Pymol.
    //
    //  It #declares two transformation statements:
    //
    //  FROM_PYMOL_VIEW
    //  TO_PYMOL_VIEW
    //
    //  The first can be used to transform the Pymol
    //  scene back to model (normal) space, the latter
    //  is used to transform other objects to appear in
    //  the scene on the correct position.
    //
    //  Additionally four macros are defined to transform
    //  vectors (points) from one space to another:
    //
    //  VEC2PYMOLSPACE( <x, y, z> )
    //  VEC2CARTSPACE( <x, y, z> )
    //  VEC2PYMOLVEC( <x, y, z> )
    //  VEC2CARTVEC( <x, y, z> ) 
    //
    //  *NEW*
    //
    //  If the view from pymol is stored as an array:
    //
    //  #declare M = array[18] {...}
    //
    //  then the macros
    //
    //  SET_PYMOL_VIEW 
    //    and 
    //  UNSET_PYMOL_VIEW
    //
    //  can be used directly to transform objects to and from that view:
    //  object { ... SET_PYMOL_VIEW( M ) }
    //
    //  This is especially useful if multiple views are defined 
    //  and the scene was set in one:
    //
    //  #declare VIEW1 = M;
    //  #declare VIEW2 = N;
    //  union { #include "file.inc" UNSET_PYMOL_VIEW( M ) SET_PYMOL_VIEW( N ) }
    //
    //  NOTE: transform statements are combined by POV-Ray prior to 
    //  transformations of objects, so there's little need to implement a macro
    //  SWITCH_PYMOL_VIEW( M, N )
    //  although that would appear simpler in the scenes  
    
    //  Tsjerk A. Wassenaar
    //  February 16, 2005
    //  April 5, 2005
    //  September 2, 2009
    //
    
    // Determinant of a matrix
    //------------------------
    #macro DET( M )
    
      #local a = M[0] * ( M[4]*M[8] - M[5]*M[7] ); 
      #local b = M[1] * ( M[3]*M[8] - M[5]*M[6] ); 
      #local c = M[2] * ( M[3]*M[7] - M[4]*M[6] );
    
      (a - b + c)
    
    #end // of DET()
    
    
    // The inverse of a matrix
    //------------------------
    #macro INV( m11, m12, m13, m21, m22, m23, m31, m32, m33 )
    
      #local M = array[9] { m11, m12, m13, m21, m22, m23, m31, m32, m33 };
      #local invdet = 1/DET( M );
    	
      #local t11 = invdet * ( m22*m33 - m23*m32 ); 
      #local t12 = invdet * ( m13*m32 - m12*m33 ); 
      #local t13 = invdet * ( m12*m23 - m13*m22 ); 
      #local t21 = invdet * ( m23*m31 - m21*m33 ); 
      #local t22 = invdet * ( m11*m33 - m13*m31 ); 
      #local t23 = invdet * ( m13*m21 - m11*m23 );
      #local t31 = invdet * ( m21*m32 - m22*m31 );
      #local t32 = invdet * ( m12*m31 - m11*m32 );
      #local t33 = invdet * ( m11*m22 - m12*m21 );
    
      t11, t12, t13, t21, t22, t23, t31, t32, t33, 0, 0, 0
    
    #end // of INV()
    
    #macro M_INV( M )
      #local invdet = 1/DET( M );
    	
      #local t11 = invdet * ( M[4]*M[8] - M[5]*M[7] ); 
      #local t21 = invdet * ( M[2]*M[7] - M[1]*M[8] ); 
      #local t31 = invdet * ( M[1]*M[5] - M[2]*M[4] ); 
    
      #local t12 = invdet * ( M[5]*M[6] - M[3]*M[8] ); 
      #local t22 = invdet * ( M[0]*M[8] - M[2]*M[6] ); 
      #local t32 = invdet * ( M[2]*M[3] - M[0]*M[5] );
    
      #local t13 = invdet * ( M[3]*M[7] - M[4]*M[6] );
      #local t23 = invdet * ( M[1]*M[6] - M[0]*M[7] );
      #local t33 = invdet * ( M[0]*M[4] - M[1]*M[3] );
    
      array[9] {t11, t12, t13, t21, t22, t23, t31, t32, t33}
    #end
    
    #macro MV_MUL( M, V )
        < M[0]*V.x + M[1]*V.y + M[2]*V.z,
          M[3]*V.x + M[4]*V.y + M[5]*V.z,
          M[6]*V.x + M[7]*V.y + M[8]*V.z >
    #end
    
    #macro SET_PYMOL_VIEW( M )
      transform {
        translate -< M[12], M[13], M[14] >
        matrix < M[0], M[1],  M[2],
    	     M[3], M[4],  M[5], 
    	     M[6], M[7],  M[8], 
    	     M[9], M[10], M[11] >
      } 
    #end // of SET_PYMOL_VIEW
    
    #macro UNSET_PYMOL_VIEW( M )
      transform {
        translate -< M[9], M[10], M[11] >
        matrix < INV( M[0], M[1], M[2], M[3], M[4], M[5], M[6], M[7], M[8] ) > 
        translate < M[12], M[13], M[14] >
      } 
    #end // of UNSET_PYMOL_VIEW
    
    #macro C2P_VEC( M, vec)
      #local nvec = vec - <M[12],M[13],M[14]>;
      #local nvec =
        < M[0]*nvec.x + M[1]*nvec.y + M[2]*nvec.z,
          M[3]*nvec.x + M[4]*nvec.y + M[5]*nvec.z,
          M[6]*nvec.x + M[7]*nvec.y + M[8]*nvec.z >; 
      nvec + <M[9],M[10],M[11]>
    #end
    
    #macro P2C_VEC( M, vec)
      MV_MUL( M_INV(M), vec - <M[9],M[10],M[11]> ) + <M[12],M[13],M[14]>
      //#local nvec = vec - <M[9],M[10],M[11]>;
      //#local N = M_INV( M ) ;
      //#local nvec =
      //  < N[0]*nvec.x + N[1]*nvec.y + N[2]*nvec.z,
      //    N[3]*nvec.x + N[4]*nvec.y + N[5]*nvec.z,
      //    N[6]*nvec.x + N[7]*nvec.y + N[8]*nvec.z >; 
      //nvec
    #end
    
    
    #macro PYMOL_VIEW( r11, r12, r13,     // 3x3 Rotation matrix ( Model space to Camera space )
    		   r21, r22, r23, 
    		   r31, r32, r33,
    		    c1,  c2,  c3,     // Camera position ( Model space )
    		    o1,  o2,  o3,     // Origin of rotation ( Model space )
    		    s1,  s2,  or)     // Slab near and far, orthoscopic flag ( discarded )
    
      #declare PYMOLVIEW_RMATRIX = array[9] { r11, r12, r13, 
    					  r21, r22, r23, 
    					  r31, r32, r33 }
      #declare PYMOLVIEW_CAMPOS  = < c1, c2, c3 >;
      #declare PYMOLVIEW_ORGPOS  = < o1, o2, o3 >;
    
      #declare TO_PYMOL_VIEW = transform {
        translate -< o1, o2, o3 >
        matrix < r11, r12, r13,
    	     r21, r22, r23, 
    	     r31, r32, r33, 
    	      c1,  c2,  c3 >
      }
    
      #declare FROM_PYMOL_VIEW = transform {
        translate -< c1, c2, c3>
        matrix < INV( r11, r12, r13, r21, r22, r23, r31, r32, r33 ) >
        translate  < o1, o2, o3>
      }
    
      #macro VEC2PYMOLSPACE(vec)
        #local nvec = vec - PYMOLVIEW_ORGPOS;
        #local nvec =
          < PYMOLVIEW_RMATRIX[0]*nvec.x + PYMOLVIEW_RMATRIX[3]*nvec.y + PYMOLVIEW_RMATRIX[6]*nvec.z,
            PYMOLVIEW_RMATRIX[1]*nvec.x + PYMOLVIEW_RMATRIX[4]*nvec.y + PYMOLVIEW_RMATRIX[7]*nvec.z,
            PYMOLVIEW_RMATRIX[2]*nvec.x + PYMOLVIEW_RMATRIX[5]*nvec.y + PYMOLVIEW_RMATRIX[8]*nvec.z >; 
        nvec + PYMOLVIEW_CAMPOS
      #end
    
      #macro VEC2CARTSPACE(vec)
    
        #local nvec = vec - PYMOLVIEW_CAMPOS;
    
        #local R = PYMOLVIEW_RMATRIX;
        #local invdet = 1/DET( R );
    
        #local T = array[9];
    
        #local T[0] = invdet * ( R[4]*R[8] - R[5]*R[7] ); 
        #local T[1] = invdet * ( R[2]*R[7] - R[1]*R[8] ); 
        #local T[2] = invdet * ( R[1]*R[5] - R[2]*R[4] ); 
        #local T[3] = invdet * ( R[5]*R[6] - R[3]*R[8] ); 
        #local T[4] = invdet * ( R[0]*R[8] - R[2]*R[6] ); 
        #local T[5] = invdet * ( R[2]*R[3] - R[0]*R[5] );
        #local T[6] = invdet * ( R[3]*R[7] - R[4]*R[6] );
        #local T[7] = invdet * ( R[1]*R[6] - R[0]*R[7] );
        #local T[8] = invdet * ( R[0]*R[4] - R[1]*R[3] );
    
        < T[0]*nvec.x + T[3]*nvec.y + T[6]*nvec.z + PYMOLVIEW_ORGPOS.x,
          T[1]*nvec.x + T[4]*nvec.y + T[7]*nvec.z + PYMOLVIEW_ORGPOS.y,
          T[2]*nvec.x + T[5]*nvec.y + T[8]*nvec.z + PYMOLVIEW_ORGPOS.z >
      #end
    
      #macro VEC2PYMOLVEC(vec)
        < PYMOLVIEW_RMATRIX[0]*vec.x + PYMOLVIEW_RMATRIX[3]*vec.y + PYMOLVIEW_RMATRIX[6]*vec.z,
          PYMOLVIEW_RMATRIX[1]*vec.x + PYMOLVIEW_RMATRIX[4]*vec.y + PYMOLVIEW_RMATRIX[7]*vec.z,
          PYMOLVIEW_RMATRIX[2]*vec.x + PYMOLVIEW_RMATRIX[5]*vec.y + PYMOLVIEW_RMATRIX[8]*vec.z >
      #end
    
      #macro VEC2CARTVEC(vec)
        #local nvec = vec - PYMOLVIEW_CAMPOS;
      #end
    
      #macro CAM2PYMOLCAM()
    
      #end
    
      #macro CAM2CARTCAM()
    
      #end
    #end
    

Retrieved from "[https://pymolwiki.org/index.php?title=PovRay&oldid=12676](https://pymolwiki.org/index.php?title=PovRay&oldid=12676)"


---

## Software Codecs

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Codecs

I've made movies with the DivX codec and those, in my experience, make the best movies. They take a while to encode, and are large, but have the best output quality. 

### Mplayer/Mencoder

I used **mencoder** to make a movie. The best encoding command line was, 
    
    
    mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts \
    codec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \
    mf://*.png -mf type=png:fps=18 -o output.avi
    

This gave me an 18 frame-per-second divx/mpeg4 movie with superb quality. The movie is large and is at [Arrestin Vs Rhodopsin (DivX)](http://vertrees.org/~tree/movie). The input are all the files ending in ".png" in the current directory; the output is output.avi, the movie. The three lines are actually on long one; concatenate them at the black slashes **\**. 

### Quicktime

Recent attention to transparencies and QuickTime show that one typically wants to set 
    
    
    set ray_opaque_background, off
    

to tell PyMol to stop using Alpha-channel for transparencies and switch to blending. 

If you're having transparency problems this could be a fix. 

Retrieved from "[https://pymolwiki.org/index.php?title=Software_Codecs&oldid=6524](https://pymolwiki.org/index.php?title=Software_Codecs&oldid=6524)"


---

## Surfaces and Voids

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Help with integrating output from programs that generate/analyse surfaces and voids in molecules. 

See also: 

  1. [http://sourceforge.net/search/?group_id=4546&words=cavity&type_of_search=mlists&limit=30](http://sourceforge.net/search/?group_id=4546&words=cavity&type_of_search=mlists&limit=30)
  2. <http://loschmidt.chemi.muni.cz/caver/>
  3. <http://alpha2.bmc.uu.se/usf/voidoo.html>
  4. <http://sts.bioengr.uic.edu/castp/>
  5. <http://www.ccl.net/cca/software/UNIX/pass/overview.shtml>



Retrieved from "[https://pymolwiki.org/index.php?title=Surfaces_and_Voids&oldid=5396](https://pymolwiki.org/index.php?title=Surfaces_and_Voids&oldid=5396)"


---

## SURFNET

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

A recipe for reading surfaces from [Roman Laskowski](http://www.ebi.ac.uk/~roman/)'s [SURFNET](https://www.ebi.ac.uk/thornton-srv/software/SURFNET/) program (for finding cavities in macromolecules) into PyMol for visualisation. 

1\. Create your surfaces in "CCP4" format in SURFNET. 
    
    
    Asides: 
    A. The "endianness" of SURFNET is set to big endian by default - see the 
    remarks about the SGI flag.  Change this if you're on a little endian machine, 
    e.g. LINUX/i386.
    B.SURFNET can be compiled against ccp4 version 5 and 6 libraries 
    by following the instructions in the SURFNET distribution and modifiying 
    the link lines at the end of ccp4link.scr to replace 
    
    $CLIB/libccp4.a 
    
    with 
    
    $CLIB/libccp4f.a $CLIB/libccp4c.a
    

2\. Use [Gerard Kleywegt](http://xray.bmc.uu.se/gerard)'s mapman from the [USF](http://structure.usc.edu/usf/) [RAVE](http://structure.usc.edu/usf/rave.html) package to convert the CCP4 density map to XPLOR format 

e.g. in a shell on LINUX: 
    
    
    $ lx_mapman
    
    MAPMAN > READ map1 gaps.den
    
    MAPMAN > WRITE map1 gaps.xplor XPLOR 
    

3\. Open the XPLOR map in PyMol 

4\. Generate a mesh or surface object from the map using isomesh or isosurface. 

e.g. on the PyMol command line: 
    
    
    isomesh gaps_mesh, gaps, 100.0
    

Retrieved from "[https://pymolwiki.org/index.php?title=SURFNET&oldid=12721](https://pymolwiki.org/index.php?title=SURFNET&oldid=12721)"


---

