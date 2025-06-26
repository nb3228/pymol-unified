# Category: Scenes

## Scene

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**scene** makes it possible to save and restore multiple scenes scene within a single session. A scene consists of the view, all object activity information, all atom-wise visibility, color, representations, and the global frame index. 

## Contents

  * 1 Usage
    * 1.1 Arguments
  * 2 Using Scene
    * 2.1 Storing scenes
    * 2.2 Scenes as Movies
    * 2.3 Auto-play through Scenes
  * 3 Examples
  * 4 PyMOL API
  * 5 Notes
  * 6 See Also
  * 7 DEVELOPMENT TO DO



## Usage
    
    
    scene [ key [, action [, message [, view [, color [, active [, rep
        [, frame [, animate [, new_key ]]]]]]]]]]
    

### Arguments

  * **key** = string, new, auto, or *: use new for an automatically numbered new scene, use auto for the current scene (if one exists), and use * for all scenes (clear and recall actions only).
  * **action** = store, recall, insert_after, insert_before, next, previous, update, rename, clear or append: (default = recall). If rename, then a new_key argument must be explicitly defined.
  * **message** = string: a text message to display with the scene.
  * **view** = 1 or 0: controls whether the view is stored {default: 1}
  * **color** = 1 or 0: controls whether colors are stored {default: 1}
  * **active** = 1 or 0: controls whether activity (objects enabled/disabled) is stored {default: 1}
  * **rep** = 1 or 0: controls whether the representations are stored {default: 1}
  * **frame** = 1 or 0: controls whether the frame is stored {default: 1}
  * **animate** = float: animation duration in seconds {default: _scene_animation_duration_}
  * **new_key** = string: the new name for the scene



## Using Scene

The Scene command has quite a few actions/options that can be enabled by using the mouse and the keyboard through the usual Scene command or hot-keys. Also, you can shift the scenes around using the new [Scene_buttons](/index.php/Scene_buttons "Scene buttons") and just dragging the scene names. 

### Storing scenes
    
    
    # store this scene in the next spot, giving it the default name.
    scene auto, store
    

has the hot-key equivalent of **CTRL-PageDown** (FN+CTRL+DownArrow on the Mac). Try turning on [Scene_Buttons](/index.php?title=Scene_Buttons&action=edit&redlink=1 "Scene Buttons \(page does not exist\)") and then doing CTRL-PageDown; see the scene buttons popping up? 

### Scenes as Movies

If you desire to make a movie that only has camera changes or representation changes, then scenes are your friend. Simply setup each view and then when ready you can do Scene->Store from the PyMOL menus (or _scene auto, store_ on the command line or the third method Ctrl+PgDn (Fn+Ctrl+DownArrow on the Mac)). Do this for each view you setup. Once done, you can scroll through your scenes by pushing PgUp/PgDn. PyMOL automatically interpolates when you use the PgUp/PgDn buttons, so you get the desired smooth transitions. Mix this with [AxPyMOL](http://www.pymol.org/ax/) and you have movies in PowerPoint with very little work. 

### Auto-play through Scenes

With this simple trick you can auto-play through scenes. This is similar to "Movie > Program > Scene Loop" but uses only a single frame. 
    
    
    cmd.mset('1x1')
    cmd.set('scene_loop')
    cmd.set('movie_fps', 1.0 / 5.0)
    cmd.mdo(1, 'scene auto, next')
    cmd.mplay()
    

## Examples

Simple Examples. 
    
    
    scene F1, store
    scene F2, store, This view shows you the critical hydrogen bond.
     
    scene F1
    scene F2
    
    scene *
    

This example shows how to use scenes in a movie! 
    
    
    # SUMMARY
    #
    
    # This script demonstrates one way of creating a movie from scenes.
    # It assumes that we have three scenes, each running for 10 seconds
    # (300 frames apiece) including 2-second transitions.
    
    # 1) Load or create content for three scenes (this could just as easily
    #    come from a session file).
    
    load $TUT/1hpv.pdb
    util.cbc
    turn x,180
    orient
    as cartoon
    scene 001, store
    
    show sticks, organic
    orient organic
    scene 002, store
    
    hide cartoon
    show lines, byres organic expand 5
    turn x,45
    turn y,45
    scene 003, store
    
    # 2) Specify a 30-second movie -- state 1, 900 frames at 30 frames per second.
    
    mset 1 x900
    
    # 3) Program scene matrices as movie views at appopriate frames
    #    and also add y-axis rocking between scenes.
    
    scene 001, animate=0
    mview store, 1
    mview store, 240
    
    turn y,-30
    mview store, 70
    turn y,60
    mview store, 170
    
    scene 002, animate=0
    mview store, 300
    mview store, 540
    
    turn y,-30
    mview store, 370
    turn y,60
    mview store, 470
    
    scene 003, animate=0
    mview store, 600
    mview store, 840
    
    turn y,-30
    mview store, 670
    turn y,60
    mview store, 770
    
    # 4) Now interpolate the movie camera.
    
    mview interpolate
    mview smooth
    mview smooth
    
    # 5) Activate scene content at the appropriate movie frames.
     
    mdo 1: scene 001, view=0, quiet=1
    mdo 240: scene 002, view=0, quiet=1
    mdo 540: scene 003, view=0, quiet=1
    mdo 840: scene 001, view=0, quiet=1
    
    # 6) Force frame 1 content to load.
    
    rewind
    
    # 6) And play the movie.
    
    mplay
    

## PyMOL API
    
    
    cmd.scene(str key='auto', str action='recall', str-or-list message=None, bool view=1, bool color=1,
        bool active=1, bool rep=1, bool frame=1, float animate=-1, str new_key=None)
    

## Notes

  * To scroll through your frames, as in a presentation, just use the PG-UP and PG-DN keys. Very handy.
  * Scenes F1 through F12 are automatically bound to function keys provided that "set_key" hasn't been used to redefine the behaviour of the respective key.
  * If you have a script that modifies the representation of the molecules and stores them, quickly, then the stored frames may not be up to date. I suggest calling "refresh" between the commands.



## See Also

[View](/index.php/View "View"), [Set_View](/index.php/Set_View "Set View"), [Get_View](/index.php/Get_View "Get View"), [Movie_from_scenes](/index.php/Movie_from_scenes "Movie from scenes")

## DEVELOPMENT TO DO

Add support for save/restore of a certain global and object-and-state specific settings, such as: state, surface_color, ribbon_color, stick_color, transparency, sphere_transparency, etc. This would probably best be done by defining a class of "scene" settings which are treated in this manner. The current workaround is to create separate objects which are enabled/disabled differentially. 

Retrieved from "[https://pymolwiki.org/index.php?title=Scene&oldid=12266](https://pymolwiki.org/index.php?title=Scene&oldid=12266)"


---

## Scene buttons

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

# Overview

This setting controls whether or not PyMOL displays Scene_buttons or openGL buttons in the lower left hand corner that allows interaction with scenes. 

[![](/images/d/db/Sb1.png)](/index.php/File:Sb1.png)

[](/index.php/File:Sb1.png "Enlarge")

frame

# Syntax
    
    
    # turn on scene buttons
    set scene_buttons, 1
    
    # turn off scene buttons
    set scene_buttons, 0
    

# Usage

Once these are enabled you can interact with them by using the mouse: 

  * left-button will select a scene
  * middle-button will rove over the scenes (to see this, load a bunch of scenes and then middle-click on the first and then drag down)
  * right-button will move a scene in the stack (drag-and-drop).



Retrieved from "[https://pymolwiki.org/index.php?title=Scene_buttons&oldid=6623](https://pymolwiki.org/index.php?title=Scene_buttons&oldid=6623)"


---

