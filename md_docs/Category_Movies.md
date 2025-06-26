# Category: Movies

## Backward

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**backward** moves the movie back one frame. 

### USAGE
    
    
    backward
    

### PYMOL API
    
    
    cmd.backward()
    

### SEE ALSO
    
    
    [Mset](/index.php/Mset "Mset"), [Forward](/index.php/Forward "Forward"), [Rewind](/index.php/Rewind "Rewind"), [Category:Movies](/index.php/Category:Movies "Category:Movies")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Backward&oldid=7436](https://pymolwiki.org/index.php?title=Backward&oldid=7436)"


---

## EMovie

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [PyMOL Plugin](/index.php/Plugins "Plugins")  
---|---  
Download  | [plugins/emovie.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/emovie.py)  
Author(s)  | [Israel Structural Proteomics Center](http://www.weizmann.ac.il/ISPC/home.html)  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 eMovie
    * 1.1 Program Features
    * 1.2 Screenshots
      * 1.2.1 eMovie in use
      * 1.2.2 eMovie menu bar
    * 1.3 Availability
    * 1.4 Authors
    * 1.5 Reference



# eMovie

eMovie is a plug-in for PyMOL that makes the creation of molecular movies both easy and intuitive via a breakthrough storyboard interface, similar in nature to what is used in the creation of traditional movies. The eMovie homepage is accessible at [www.weizmann.ac.il/ISPC/eMovie.html](http://www.weizmann.ac.il/ISPC/eMovie.html). 

## Program Features

eMovie is currently the most user-friendly way for users to create movies in PyMOL (even inexperienced users). Users interact with a user-friendly eMovie GUI that does not require typing commands into PyMOL. 

Modular actions such as zooms, rotations, fadings, and morphs (morphs require incentive PyMOL) can be inserted to any frame in the movie and the actions comprising the movie can be reviewed in list-format by viewing the eMovie storyboard. The storyboard also allows for deletion and reinsertion of actions. 

Movies can be saved, loaded, and exported as a series of image files to be later merged into a traditional movie format such as .mov using an external program like GraphicConverter. 

## Screenshots

### eMovie in use

[![EMovie in use.jpg](/images/3/38/EMovie_in_use.jpg)](/index.php/File:EMovie_in_use.jpg)

### eMovie menu bar

[![EMovie menubar.jpg](/images/3/3d/EMovie_menubar.jpg)](/index.php/File:EMovie_menubar.jpg)

## Availability

eMovie is freely available for download at the [eMovie homepage](http://www.weizmann.ac.il/ISPC/eMovie.html). 

## Authors

eMovie was created at the [Israel Structural Proteomics Center](http://www.weizmann.ac.il/ISPC/home.html) (ISPC) at the Weizmann Institute of Science. 

## Reference

Hodis, E., Schreiber, G., Rother, K., Sussman, J.L., eMovie: a storyboard-based tool for making molecular movies, Trends in Biochemical Sciences 32, 199-204 (2007). 

Retrieved from "[https://pymolwiki.org/index.php?title=EMovie&oldid=10436](https://pymolwiki.org/index.php?title=EMovie&oldid=10436)"


---

## Forward

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**forward** moves the movie one frame forward. 

### USAGE
    
    
    forward
    

### PYMOL API
    
    
    cmd.forward()
    

### SEE ALSO

[Mset](/index.php/Mset "Mset"), [Backward](/index.php/Backward "Backward"), [Rewind](/index.php/Rewind "Rewind"), [Category:Movies](/index.php/Category:Movies "Category:Movies")

Retrieved from "[https://pymolwiki.org/index.php?title=Forward&oldid=7475](https://pymolwiki.org/index.php?title=Forward&oldid=7475)"


---

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

## Movie fps

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The movie_fps setting only affects the playback speed within PyMOL and not the final exported movie. 

## See Also

  * [mplay](/index.php/Mplay "Mplay")
  * [movie.produce](/index.php/Movie.produce "Movie.produce")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_fps&oldid=10639](https://pymolwiki.org/index.php?title=Movie_fps&oldid=10639)"


---

## Movie loop

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Overview

This setting controls whether or not PyMOL loops your movies. 

## Syntax
    
    
    # turn off looping
    set movie_loop, 0
    
    # turn on looping
    set movie_loop, 1
    

## Examples

[movie_delay](/index.php?title=Movie_delay&action=edit&redlink=1 "Movie delay \(page does not exist\)")

Retrieved from "[https://pymolwiki.org/index.php?title=Movie_loop&oldid=6525](https://pymolwiki.org/index.php?title=Movie_loop&oldid=6525)"


---

## Movie panel

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

The Movie_panel is a time-line view (placed at the bottom of the screen when activated) that the user can interact with to more easily produce movies. A frame counter is shown as a small cursor. This can be dragged along the timeline to quickly access other frames in the movie. Whenever a camera setting is stored on a frame a blue colored cursor appears. 

[![](/images/a/a0/Showing_movie_panel.jpeg)](/index.php/File:Showing_movie_panel.jpeg)

[](/index.php/File:Showing_movie_panel.jpeg "Enlarge")

Notice the strip at the bottom of the screen showing the current frame. Notice also "Frame 50" in the right menu.

# Syntax
    
    
    # turn on the movie panel
    set movie_panel, 1
    

# See Also

[Category:Movies](/index.php/Category:Movies "Category:Movies")

Retrieved from "[https://pymolwiki.org/index.php?title=Movie_panel&oldid=6579](https://pymolwiki.org/index.php?title=Movie_panel&oldid=6579)"


---

## Movie.produce

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**movie.produce** can be used to generate a movie file (mpg) after using the movie commands, see [Movie_school](/index.php/Movie_school "Movie school")

  


## Contents

  * 1 USAGE
  * 2 ARGUMENTS
  * 3 EXAMPLES
  * 4 SEE ALSO



### USAGE
    
    
    movie.produce filename [, mode [, first [, last [, preserve [, encoder [, quality [, quiet ]]]]]]]
    

### ARGUMENTS

**filename**

    

    output file

**mode**

    

    normal {default}
    draw
    ray

**first**

    

    start frame {default: first frame}

**last**

    

    end frame {default: last frame}

**preserve**

    

    = 0: cleanup frames from temp directory {default}
    = 1: leave temp directory & frames

**encoder**

    

    which encode to use {default freemol's mpeg_encode}

**quality**

    

    integer 0..100: for visual quality of final movie; corresponds to movie_quality setting

**quiet**

    

    verbosity

### EXAMPLES
    
    
    # generate myfile.mpg
    movie.produce myfile.mpg
    
    # generate movie myfile.mpg ray'ing frames 20-100 while preserving the ppm images of each frame in the myfile.tmp directory
    movie.produce myfile.mpg, ray, 20, 100, 1
    

### SEE ALSO

  * [Movie_school](/index.php/Movie_school "Movie school")
  * [Movie_from_scenes](/index.php/Movie_from_scenes "Movie from scenes")
  * <http://freemol.org>



Retrieved from "[https://pymolwiki.org/index.php?title=Movie.produce&oldid=9160](https://pymolwiki.org/index.php?title=Movie.produce&oldid=9160)"


---

## Movie.roll

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This function will perform a movie roll--a 360 rotation of the camera about the scene--from the start frame to the end frame specified. So, if you specify fewer frames it'll roll faster than in you specified more frames. 

# Syntax
    
    
    # roll the scene from state 'first' to 'last' rolling over axis 'axis' and looping (1) or not (0).
    movie.roll [first=1 [, last=-1 [, loop=1 [, axis='y' ]]]]
    
    # roll the scene over states 1 through 50.
    movie.roll 1, 50
    

Retrieved from "[https://pymolwiki.org/index.php?title=Movie.roll&oldid=6534](https://pymolwiki.org/index.php?title=Movie.roll&oldid=6534)"


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

## Mview

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The mview command can store and delete movie keyframes. 

Keyframes store a view (camera or object position) and optionally the object [state](/index.php/State "State") and/or a [scene](/index.php/Scene "Scene"). Between keyframes, PyMOL will interpolate views and states, allowing for smooth animations. 

Before using mview, the movie timeline has to be set up with [mset](/index.php/Mset "Mset"). 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
    * 3.1 Ligand zoom
    * 3.2 360° rotation
    * 3.3 360° rotation of a single object
    * 3.4 Object-level state-sweep
    * 3.5 Ligand binding
    * 3.6 Scene based movie
  * 4 See Also



## Usage
    
    
    mview [ action [, first [, last [, power [, bias
        [, simple [, linear [, object [, wrap [, hand
        [, window [, cycles [, scene [, cut [, quiet
        [, auto [, state [, freeze ]]]]]]]]]]]]]]]]]]
    

## Arguments

  * **action** = str: one of store, clear, reset, purge, interpolate, uninterpolate, reinterpolate, toggle, toggle_interp, smooth {default: store}
  * **first** = int: frame number or 0 for current frame {default: 0}
  * **power** = float: slow down animation at keyframe (0.0) or not (1.0) {default: 0.0}
  * **object** = str: name of object for object keyframes, or empty for global (camera) keyframes {default: }
  * **scene** = str: name of scene to store scene with key frame {default: }
  * **cut** = float 0.0-1.0: scene switch moment (0.0: beginning of transition, 1.0: end of transition) {default: 0.5}
  * **auto** = -1/0/1: if freeze=0, then auto reinterpolate after store, clear, or toggle {default: -1 = use [movie_auto_interpolate](/index.php?title=Movie_auto_interpolate&action=edit&redlink=1 "Movie auto interpolate \(page does not exist\)")}
  * **state** = int: if > 0, then store object state {default: 0}
  * **freeze** = 0/1: never auto reinterpolate {default: 0}



## Examples

### Ligand zoom
    
    
    fetch 1rx1, async=0
    as cartoon
    as sticks, organic
    mset 1x70
    orient
    mview store, 1
    mview store, 70
    orient organic
    mview store, 30
    mview store, 40
    mplay
    

### 360° rotation
    
    
    fragment ala
    mset 1x90
    mview store, 1
    mview store, 90
    turn y, 120
    mview store, 30, power=1.0
    turn y, 120
    mview store, 60, power=1.0
    mplay
    

### 360° rotation of a single object
    
    
    set movie_auto_store, 0
    fragment ala
    fragment his
    translate [10, 0, 0], his
    zoom
    mset 1x90
    mview store, 1, object=his
    mview store, 90, object=his
    rotate y, 120, object=his  # keyword >>object<< is absolutely necessary! Otherwise movie will not work.
    mview store, 30, power=1.0, object=his
    rotate y, 120, object=his
    mview store, 60, power=1.0, object=his
    mplay
    

### Object-level state-sweep
    
    
    load <http://pymol.org/tmp/morph.pdb.gz>
    dss
    as cartoon
    mset 1x80
    mview store, 1, object=m, state=1
    mview store, 30, object=m, state=30
    mview store, 50, object=m, state=30
    mview store, 80, object=m, state=1
    mplay
    

### Ligand binding
    
    
    set movie_auto_store, 0
    fetch 1rx1, async=0
    extract ligand, organic
    as cartoon, 1rx1
    as sticks, ligand
    set_view (\
      0.527486444, -0.761115909, -0.377440333,\
      0.736519873, 0.631122172, -0.243357301,\
      0.423434794, -0.149625391, 0.893482506,\
      0.000059791, -0.000049331, -140.287048340,\
      34.670463562, 51.407436371, 17.568315506,\
      111.284034729, 169.290832520, -19.999998093 )
    mset 1x60
    mview store, 60, object=ligand
    translate [10, 0, 0], object=ligand
    mview store, 1, object=ligand
    mplay
    

### Scene based movie
    
    
    fragment ala
    as sticks
    color blue
    scene bluesticks, store
    as spheres
    color red
    turn y, 180
    scene redspheres, store
    mset 1x60
    mview store, 1, scene=bluesticks
    mview store, 30, scene=redspheres
    mplay
    

## See Also

  * [mset](/index.php/Mset "Mset")
  * [frame](/index.php/Frame "Frame")
  * [mdo](/index.php/Mdo "Mdo")
  * [mpng](/index.php/Mpng "Mpng")
  * [movie_loop](/index.php/Movie_loop "Movie loop") setting
  * [movie_auto_interpolate](/index.php?title=Movie_auto_interpolate&action=edit&redlink=1 "Movie auto interpolate \(page does not exist\)") setting
  * [movie_auto_store](/index.php?title=Movie_auto_store&action=edit&redlink=1 "Movie auto store \(page does not exist\)") setting



Retrieved from "[https://pymolwiki.org/index.php?title=Mview&oldid=12641](https://pymolwiki.org/index.php?title=Mview&oldid=12641)"


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

## Ray

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**ray** creates a ray-traced image of the current frame. 

[![](/images/b/bc/Gslike.png)](/index.php/File:Gslike.png)

[](/index.php/File:Gslike.png "Enlarge")

Varying settings to play with rendering options

## Contents

  * 1 Details
    * 1.1 Usage
    * 1.2 PyMol API
    * 1.3 Settings
      * 1.3.1 Modes
      * 1.3.2 Perspective
        * 1.3.2.1 Perspective Example Images
          * 1.3.2.1.1 Notes
      * 1.3.3 Renderer
    * 1.4 Performance
      * 1.4.1 Memory
    * 1.5 Examples
      * 1.5.1 Simple
      * 1.5.2 Specify Image Size
      * 1.5.3 Specify Renderer
      * 1.5.4 High Quality B&W Rendering
      * 1.5.5 High Quality Color
      * 1.5.6 Ray Tracing Stereo Images
    * 1.6 See also
    * 1.7 User comments



# Details

This command is used to make high-resolution photos fit for publication and formal movies. Please note, the **ray** command can take some time (up to several minutes, depending on image complexity and size). 

For those who are making movies with PyMOL, **Ray** is one of the most commonly used last steps before stitching the frames together to compile the movie. **Ray** has many powerful features such as setting the size of the image -- and it still works even if the [Viewport](/index.php/Viewport "Viewport") or screen is smaller than the requested output file dimensions. 

[![](/images/b/b9/No_ray_trace.png)](/index.php/File:No_ray_trace.png)

[](/index.php/File:No_ray_trace.png "Enlarge")

Image, not ray traced.

[![](/images/3/30/Ray_traced.png)](/index.php/File:Ray_traced.png)

[](/index.php/File:Ray_traced.png "Enlarge")

Image, ray traced.

## Usage
    
    
    ray [width,height [,renderer [,angle [,shift ]]]
    

**angle** and **shift** can be used to generate matched stereo pairs 

**width** and **height** can be set to any non-negative integer. If both are set to zero than the current window size is used and is equivalent to just using **ray** with no arguments. If one is set to zero (or missing) while the other is a positive integer, then the argument set to zero (or missing) will be scaled to preserve the current aspect ratio. 

## PyMol API
    
    
    cmd.ray(int width,int height,int renderer=-1,float shift=0)
    

## Settings

### Modes

Setting the **[Ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")** variable in PyMOL changes the way PyMOL's internal renderer represents proteins in the final output. New modes were recently added to give the user more options of molecular representation. New modes are: normal rendering, but with a black outline (nice for presentations); black and white only; quantized color with black outline (also, very nice for presentations; more of a _cartoony_ appearance). 

**Note:** Mode 3, the quantized color one, sort of **burns** the background if you're using this setting. This will make a pure white background somewhat "offwhite"; thus, a poster would look poor because you could see the border for the image. If you'll be using this mode, try the [ray_opaque_background](/index.php/Ray_opaque_background "Ray opaque background") setting. 
    
    
    # normal color
    set ray_trace_mode, 0
    
    # normal color + black outline
    set ray_trace_mode,  1
    
    # black outline only
    set ray_trace_mode,  2
    
    # quantized color + black outline
    set ray_trace_mode,  3
    
    set ray_trace_mode, 1 # (or 2 or 3; best with "bg_color white;set antialias,2")
    # These two new modes -- 2 and 3 -- are cool cartoon looking modes.
    
    # change the color of the outline to a named color, or a hex-code
    set ray_trace_color, magenta
    set ray_trace_color, 0x0033ff
    

Here are the example images for the new modes 

  * [![set ray_trace_mode,1](/images/a/a8/Ray_mode_1_ex.png)](/index.php/File:Ray_mode_1_ex.png "set ray_trace_mode,1")

set ray_trace_mode,1 

  * [![set ray_trace_mode,2](/images/9/95/Ray_mode_2_ex.png)](/index.php/File:Ray_mode_2_ex.png "set ray_trace_mode,2")

set ray_trace_mode,2 

  * [![set ray_trace_mode,3](/images/5/51/Ray_mode_3_ex.png)](/index.php/File:Ray_mode_3_ex.png "set ray_trace_mode,3")

set ray_trace_mode,3 




### Perspective

#### Perspective Example Images

  * [![Example with Perspective Turned Off](/images/3/31/No_persp.png)](/index.php/File:No_persp.png "Example with Perspective Turned Off")

Example with Perspective Turned Off 

  * [![Example with Perspective Turned On](/images/c/cc/Persp1.png)](/index.php/File:Persp1.png "Example with Perspective Turned On")

Example with Perspective Turned On 

  * [![Example with Perspective Turned On and Field of View Set High \(70\).](/images/8/88/Persp2.png)](/index.php/File:Persp2.png "Example with Perspective Turned On and Field of View Set High \(70\).")

Example with Perspective Turned On and Field of View Set High (70). 




##### Notes

PyMol 0.97 and prior used **orthoscopic** rendering -- that is, no perspective. Upon the arrival of 0.98 and later, we get perspective based rendering at a cost of a 4x decrease in render speed. If you want perspective 
    
    
    set orthoscopic, off
    

Otherwise 
    
    
    set orthoscopic, on
    

To magnify the effect of perspective on the scene, 
    
    
    set field_of_view, X
    

where 50<X<70\. Default is 20. 50-70 gives a very strong perspective effect. Nb. the field of view is in Y, not X as one would expect. 

  


### Renderer

**renderer = -1** is default (use value in ray_default_renderer) 

**renderer = 0** uses PyMOL's internal renderer 

**renderer = 1** uses PovRay's renderer. This is Unix-only and you must have "povray" in your path. It utilizes two temporary files: "tmp_pymol.pov" and "tmp_pymol.png". Also works on Mac via Povray37UnofficialMacCmd but it needs to be in your path as "povray". 

## Performance

  * The ray performance depends on distance between camera and molecule.



If the distance is big rendering takes much time. If the distance is too small distant parts of molecule dissolve. 

  * [![Too close to molecule](/images/7/70/Close_ray.png)](/index.php/File:Close_ray.png "Too close to molecule")

Too close to molecule 

  * [![Normal distance](/images/1/1c/Middle_ray.png)](/index.php/File:Middle_ray.png "Normal distance")

Normal distance 



  * Tip: If you have a rather complicated scene that is zoomed into only a part of the molecule, you can speed up the ray tracing by hiding everything else outside of a certain range of the zoomed-on point. For example, if I have a large molecule and I'm looking only at the 30-atom ligand bound to it, then I can do something like the following:


    
    
    # setup your complex scene
    ...
    
    # zoom on the hetero atom (ligand and not water) within 5 Angstroms
    select hh, het and not resn HOH
    zoom hh, 5
    
    # turn on depth cueing
    set depth_cue, 1
    
    # now, select stuff to hide; we select everything that is 
    # farther than 8 Ang from our main selection
    select th, (all) and not ( (all) within 8 of hh) )
    
    hide everything, th
    
    # any additional commands you want
    ...
    
    ray
    

As an example of the efficacy of this method, I ray traced a rather complex scene with all the atoms visible here's the output of ray: 
    
    
    PyMOL>ray
     Ray: render time: 24.50 sec. = 146.9 frames/hour (941.88 sec. accum.).
    

and here is the result when I soft-clipped everything else using the above method: 
    
    
    PyMOL>ray
     Ray: render time: 47.93 sec. = 75.1 frames/hour (989.80 sec. accum.).
    

The two images in the following gallery show the results of the rendering. 

  * [![normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI.](/images/1/1a/Ray_method_off.png)](/index.php/File:Ray_method_off.png "normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI.")

normal ray tracing. This took twice as long to make as the image to the right. Same size, and DPI. 

  * [![manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left.](/images/7/7d/Ray_method_on.png)](/index.php/File:Ray_method_on.png "manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left.")

manually hiding things you won't see anyway. This took 1/2 the time to render as compared to the same sized & DPId image at left. 




### Memory

If memory is an issue for you in PyMOL, try executing your rendering from a script rather than a PyMOL session file. An unfortunate unavoidable consequence of the fact that we use Python's portable, platform-independent "pickle" machinery for PyMOL session files. Packing or unpacking a Session or Scene file thus requires that there be two simultanous copies of the information to reside in RAM simultaneously: one native and a second in Python itself. 

So when memory is a limiting factor, scripts are recommended over sessions. 

## Examples

### Simple
    
    
    # ray trace the current scene using the default size of the viewport
    ray
    

### Specify Image Size
    
    
    # ray trace the current scene, but scaled to 1024x768 pixels
    ray 1024,768
    

### Specify Renderer
    
    
    # ray trace with an external renderer.
    ray renderer=0
    

### High Quality B&W Rendering

[![](/images/6/6c/1l9l.png)](/index.php/File:1l9l.png)

[](/index.php/File:1l9l.png "Enlarge")

Black and White (ray_trace_mode,2); click to see full image
    
    
    # Black and White Script
    load /tmp/3fib.pdb;
    show cartoon;
    set ray_trace_mode, 2;  # black and white cartoon
    bg_color white;
    set antialias, 2;
    ray 600,600
    png /tmp/1l9l.png
    

### High Quality Color

[![](/images/7/7e/1l9l_2.png)](/index.php/File:1l9l_2.png)

[](/index.php/File:1l9l_2.png "Enlarge")

Color mode (ray_trace_mode,3); click to see full image
    
    
    # Color Script
    load /tmp/thy_model/1l9l.pdb;
    hide lines;
    show cartoon;
    set ray_trace_mode, 3; # color
    bg_color white;
    set antialias, 2;
    remove resn HOH
    remove resn HET
    ray 600,600
    png /tmp/1l9l.png
    

### Ray Tracing Stereo Images

    _See[Stereo_Ray](/index.php/Stereo_Ray "Stereo Ray")_

## See also

  1. "help faster" for optimization tips with the builtin renderer. "help povray" for how to use PovRay instead of PyMOL's built-in ray-tracing engine. For high-quality photos, please also see the [Antialias](/index.php/Antialias "Antialias") command. [Ray shadows](/index.php/Ray_shadows "Ray shadows") for controlling shadows.
  2. See also [Ray Tracing](/index.php/Ray_Tracing "Ray Tracing").
  3. [Desaturation Tutorial](http://www.gimp.org/tutorials/Color2BW) \-- A good resource for making nice B&W images from color images (desaturation).
  4. [Ray Trace Gain](/index.php/Ray_Trace_Gain "Ray Trace Gain")



## User comments

How do I ray trace a publication-ready (~300dpi) image using PyMol?
    This answer is in the [Advanced Issues](/index.php/Category:Advanced_Issues "Category:Advanced Issues") (Image Manipulation Section).

Retrieved from "[https://pymolwiki.org/index.php?title=Ray&oldid=11725](https://pymolwiki.org/index.php?title=Ray&oldid=11725)"


---

## Rewind

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**rewind** goes to the beginning of the movie. 

### USAGE
    
    
    rewind
    

### PYMOL API
    
    
    cmd.rewind()
    

# See Also
    
    
    [Forward](/index.php/Forward "Forward"), [Backward](/index.php/Backward "Backward"), [Category:Movies](/index.php/Category:Movies "Category:Movies")
    

Retrieved from "[https://pymolwiki.org/index.php?title=Rewind&oldid=7597](https://pymolwiki.org/index.php?title=Rewind&oldid=7597)"


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

## States

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Introduction to States

Working with States in PyMOL is a very common task. A State is one particular conformation of an object. For example, one could load an NMR ensemble and set [all_states](/index.php/All_states "All states") or use the command [Split_States](/index.php/Split_States "Split States") to see all entries in the NMR file. Another example could be the set of states from a molecular dynamic (MD) simulation. Each conformer in the MD ensemble could be viewed as a state in PyMOL. 

  * [![All states in an ensemble are shown.](/images/7/76/All_states_on.png)](/index.php/File:All_states_on.png "All states in an ensemble are shown.")

All states in an ensemble are shown. 

  * [![The states are hidden.](/images/6/6f/All_states_off.png)](/index.php/File:All_states_off.png "The states are hidden.")

The states are hidden. 




If you are making a movie of a static coordinate set (such as a single crystal structure) then you have only one state. All objects in PyMOL can potentially consist of multiple states. For movies, see [Frames](/index.php/Frame "Frame"). 

# Using States in PyMOL

  * Upon loading an object you can separate each state into its own object with [Split_States](/index.php/Split_States "Split States").
  * States can be [colored](/index.php/Color "Color") individually.
  * One can show all states with the [all_states](/index.php/All_states "All states") setting.
  * One can iterate over states using [Iterate_State](/index.php/Iterate_State "Iterate State") or change properties with states using the [Alter_State](/index.php/Alter_State "Alter State") function.



# See Also

[PyMOL States Related Pages](/index.php/Category:States "Category:States")

Retrieved from "[https://pymolwiki.org/index.php?title=States&oldid=6528](https://pymolwiki.org/index.php?title=States&oldid=6528)"


---

