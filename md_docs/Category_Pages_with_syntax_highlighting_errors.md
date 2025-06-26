# Category: Pages with syntax highlighting errors

## ColorByDistance

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [Zhenting Gao](/index.php/User:Zhentg "User:Zhentg")  
License  |   
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Examples
  * 6 The Code



## Introduction

Color the surface of a binding site by the distance to a specific atom. 

## Usage

  * Open PyMOL
  * Load PDB files
  * run this Python script inside PyMOL
  * call the function 
    * colorByDistance anchor, bindingSite



## Required Arguments

Text 

## Optional Arguments

Text 

## Examples

Text 

## The Code
    
    
    def colorByDistance(anchor='anchor', site='site'):
    	#Based on: https://pymolwiki.org/index.php/Spectrum#Intermediate
    	#Author: Zhenting Gao (zhentgpicasa@gmail.com)
    	#Update: 11/12/2018
    	#Aim: Color atoms of a binding site based on their distance from a point
    	#May also refer to https://pymolwiki.org/index.php/Ramp_New#Ramp_.2B_Distance_Measure
    	# returns the length of the distance between atom A and atom B
    
    	diff_len = lambda x,y : math.sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]))
    
    	# fetch site from the PDB
    
    	#fetch site, async=0
    
    	# show it as surface
    
    	#as surface
    
    	# create the pseudoatom at the origin
    
    	cmd.pseudoatom("pOrig", anchor, label="origin")
    
    	# these are special PyMOL variables that will hold # the coordinates of 
    	# the atoms and the  pseudoatom
    
    	stored.origCoord = []
    	stored.distCoord = []
    
    	# copy the coordinates into those special variables 
    
    	cmd.iterate_state(1, "pOrig", 'stored.origCoord.append((x,y,z))')
    	cmd.iterate_state(1, site, 'stored.distCoord.append((x,y,z))')
    
    	# extend origCoord to be the same length as the other
    
    	stored.origCoord *= len(stored.distCoord)
    
    	# calculate the distances
    
    	stored.newB = map(lambda x,y: diff_len(x,y), stored.distCoord, stored.origCoord)
    	#print(stored.newB)
    	# put them into the b-factor of the protein
    
    	cmd.alter( site, "b=stored.newB.pop(0)")
    
    	# color by rainbow_rev or any other
    	# palette listed in "help spectrum"
    
    	cmd.spectrum(expression="b", palette="rainbow", selection=site)
    	cmd.set("surface_color","-1",site) #color the surface of the binding site by corresponding atom colors
    cmd.extend('colorByDistance', colorByDistance)

Retrieved from "[https://pymolwiki.org/index.php?title=ColorByDistance&oldid=12817](https://pymolwiki.org/index.php?title=ColorByDistance&oldid=12817)"


---

## ColorByGroup

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  |   
Author(s)  | [Zhenting Gao](/index.php/User:Zhentg "User:Zhentg")  
License  |   
  
## Contents

  * 1 Introduction
  * 2 Usage
  * 3 Required Arguments
  * 4 Optional Arguments
  * 5 Examples
  * 6 The Code



## Introduction

Color the objects by their groups. 

## Usage

  * Open PyMOL
  * Load PDB files, create groups
  * run this Python script inside PyMOL
  * call the function 
    * colorByGroup



## Required Arguments

Text 

## Optional Arguments

Text 

## Examples

Text 

## The Code
    
    
    def colorByGroup(palette="rainbow"):
    	#Based on: https://pymolwiki.org/index.php/Get_Names#ARGUMENTS
    	#Author: Zhenting Gao (zhentgpicasa@gmail.com)
    	#Update: 11/12/2018
    	#Aim: Color the objects by their groups.
    	# returns the length of the distance between atom A and atom B
    	groups=cmd.get_names("group_objects")
    	numClusters=len(groups)
    	import math
    	for x in range(numClusters):
    		cmd.alter(groups[x],"b="+str(x))
    	cmd.spectrum( "b", palette, " ".join(groups))
    	#https://pymolwiki.org/index.php/Spectrum
    cmd.extend('colorByGroup', colorByGroup)

Retrieved from "[https://pymolwiki.org/index.php?title=ColorByGroup&oldid=12816](https://pymolwiki.org/index.php?title=ColorByGroup&oldid=12816)"


---

## Get color index

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**get_color_index** in combination with **get_color_tuple** will retrieve the RGB values for a single color. 
    
    
    PyMOL>print cmd.get_color_index('black')
    1
    PyMOL>print cmd.get_color_tuple(1)
    (0.0, 0.0, 0.0)

  


## See Also

  * [Get_Color_Indices](/index.php/Get_Color_Indices "Get Color Indices")
  * [Get_Color_Tuples](/index.php/Get_Color_Tuples "Get Color Tuples")



Retrieved from "[https://pymolwiki.org/index.php?title=Get_color_index&oldid=12665](https://pymolwiki.org/index.php?title=Get_color_index&oldid=12665)"


---

## Timeline Python API

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Concepts & Classes
    * 1.1 Keyframe
    * 1.2 Clip
    * 1.3 Program
    * 1.4 Track
    * 1.5 Composition
    * 1.6 Sequence Composition
  * 2 Timeline Python API
    * 2.1 General comments of the API
    * 2.2 Explanation of terms
    * 2.3 API



## Concepts & Classes

**Note:** This is a incentive-only 3.0+ feature 

### Keyframe

Represents a state of an entity at a given time. Here, state does not necessarily pertain to the PyMOL concept of state but any property (transform, PyMOL state, setting, etc…) for a given entity (PyMOL Object or Camera). 

### Clip

Represents a collection of keyframes that describes an animation over time. 

### Program

Represents an animation of an object as a function of t. Where t is the interpolation factor from [0, 1] (0 represents the beginning of the animation and 1 represents the end). Visually, they are represented as clips with a keyframe count of two. PyMOL supplies several builtin programs including: nutate, rock, roll, state loop, etc… See Creating Custom Programs below to add complex, custom animations in the timeline. 

### Track

Represents a set of related keyframes and clips. These keyframes and clips are interpolated to and through each other within the same track and not outside. For each animated property, there shall be one Timeline track. 

### Composition

Represents a set of tracks. Allows for multiple objects and/or properties to be animated simultaneously. 

### Sequence Composition

A special composition that only stores clips of “normal” compositions. This is useful to combine multiple compositions into a single animation and blend between them seamlessly. Currently, the blending effect is not yet implemented. 

## Timeline Python API

### General comments of the API

All timeline functions force named-arguments for all parameters of each function. Thus, you must specify the name of the parameter followed by the argument value you want to pass. Examples are shown below. Also, uncommon from the rest of most of the other PyMOL API, builtin types are not commonly passed nor returned. Parameter and return types are specified for each function and a full description of their meaning is explained above. Examples of movies using the Python API can be found at: www.pymol.org. 

### Explanation of terms

  * _**Arguments**_ : a sequence of named parameters. For the timeline API, all arguments are named; positional arguments are disallowed. No arguments needed are shown by an empty sequence []. Defaulted arguments are shown by a proceeding assignment with value (ex: myParam = “myValue”). All non-defaulted parameters in this sequence must be provided. Sequence representation is for the purpose of documentation only. A sequence is not expected as input for the API.
  * _**Return**_ : a sequence of returned value(s). No returned values are depicted by empty sequence []. Sequence representation is for the purpose of documentation only. A sequence is not expected as input for the API.
  * _**Description**_ : explanation of the use of the API function call.
  * _**Example**_ : brief example of how to use the API.
  * _**Query**_ : As these functions mutate the timeline state, developers may often need to retrieve the current state of a timeline composition or element. These properties can be fetched from a variable that is synchronized with the Timeline backend and doesn’t require a call via `cmd`.
  * _**See Also**_ : Related timeline API calls. Many of these examples could be paired with the presented function.



### API

timeline_add_composition

_Signature_
    
    
    def timeline_add_composition() -> Composition
    

_Description_

    Appends a new composition to the timeline.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    

_See Also_

    

timeline_delete_composition

* * *

timeline_delete_composition

_Signature_
    
    
    def timeline_delete_composition(comp: Composition) -> None
    

_Description_

    Deletes the composition from the timeline. The composition and all track, clip, and keyframe instances managed under the deleted composition are thus invalidated and should not be used in the API again.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    cmd.timeline_delete_composition(comp=myComp)
    # myComp now invalidated
    

_See Also_

    

_timeline_add_composition_

* * *

timeline_add_sequence_composition

_Signature_
    
    
    def timeline_add_sequence_composition() -> Composition
    

_Description_

    Appends a new sequence composition to the Timeline.

_Example_
    
    
    mySeqComp = cmd.timeline_add_sequence_composition()
    

_See Also_

    

_timeline_add_composition_

* * *

timeline_set_composition_duration

_Signature_
    
    
    def timeline_set_composition_duration(comp: Composition, duration: int) -> None:
    

_Description_

    Sets the duration of comp to duration number of frames.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    cmd.timeline_set_composition_duration(comp=myComp, duration=120)
    # myComp is now 2 seconds (assuming myComp is set to at 60 fps (default))
    

_Query_

    Duration can be queried by duration.

_Example_
    
    
    cmd.timeline_set_composition_duration(comp=myComp, duration=240)
    print(myComp.duration) # prints “240”
    

_See Also_

    

timeline_add_composition

* * *

timeline_scrub_to

_Signature_
    
    
    def timeline_scrub_to(comp: Composition, frame_num: int)
    

_Description_

    Sets the current frame number of the composition to the given frame number. The argument frame_num may be outside the bounds of the composition’s duration but will not be reached during the animation.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    cmd.timeline_scrub_to(comp=myComp, frame_num=50)
    

_See Also_

    

timeline_set_composition_duration

* * *

timeline_add_object_track

_Signature_
    
    
    def timeline_add_object_track(comp: Composition, obj_name: str, track_type: TrackType=TrackType.TRANSFORM, track_type_detail: str="") -> Track
    

_Description_

    Appends a new object track onto comp. obj_name must be a valid name of an object currently managed in PyMOL.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    foo_track = cmd.timeline_add_object_track(comp=myComp, obj_name="foo")
    

_See Also_

    

timeline_delete_track

* * *

timeline_delete_track

_Signature_
    
    
    def timeline_delete_track(track: Track) -> None:
    

  
_Description_

    Deletes track from its composition. The track and all clip and keyframe instances managed under the deleted track are thus invalidated and should not be used in the API again.

_Example_
    
    
    myComp = cmd.timeline_add_composition()
    foo_track = cmd.timeline_add_object_track(comp=myComp, obj_name="foo")
    

_See Also_

    

timeline_add_object_track

* * *

timeline_store_scene_at

_Signature_
    
    
    def timeline_store_scene_at(comp: Composition, frame_num: int, scene_name: str) -> Keyframe
    

_Description_

    Appends a new scene keyframe onto comp at frame frame_num with the scene of name scene_name.

_Example_
    
    
    scene_kf = cmd.timeline_store_scene_at(comp=myComp, frame_num=42, scene_name="foo")
    cmd.timeline_delete_keyframe(keyframe=scene_kf)

_See Also_

    

timeline_add_keyframe

* * *

timeline_add_keyframe

_Signature_
    
    
    def timeline_add_keyframe(track: Track, frame_num: int, obj_state: int = -1) -> Keyframe
    

_Description_

    Appends a new object keyframe onto track at frame frame_num.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    

_See Also_

    

timeline_add_keyframe

* * *

timeline_delete_keyframe

_Signature_
    
    
    def timeline_delete_keyframe(keyframe: Keyframe) -> None
    

_Description_

    Deletes a keyframe

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_delete_keyframe(keyframe=bar_kf)
    

* * *

timeline_set_keyframe_interpolation

_Signature_
    
    
    def timeline_set_keyframe_interpolation(keyframe: Keyframe, interp_mode, KeyframeInterpolationMode) -> None
    

_Description_

    Describes how the properties between keyframe and the following keyframe are interpolated. By default, virtually all keyframes interpolate linearly (each time step between two keyframe changes in property by the same amount). Values of _KeyframeInterpolationMode_ include: 

    _CONSTANT_ ,
    _LINEAR_ ,
    _EASEIN_ ,
    _EASEOUT_ ,
    _EASEINOUT_

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_set_keyframe_interpolation(keyframe=bar_kf, interp_mode=KeyframeInterpolationMode.LINEAR)
    

* * *

timeline_move_keyframe_to

_Signature_
    
    
    def timeline_move_keyframe_to(keyframe: Track, frame_num: int) -> None
    

_Description_

    Moves keyframe to frame frame_num.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_move_keyframe_to(keyframe=bar_kf, frame_num=24)
    

_See Also_

    

timeline_move_clip_to

* * *

timeline_move_clip_to

_Signature_
    
    
    def timeline_move_clip_to(clip: Clip, frame_num: int)
    

_Description_

    Moves clip to frame frame_num.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_delete_keyframe(keyframe=bar_kf)
    

_See Also_

    

timeline_move_keyframe_to

* * *

timeline_set_clip_duration

_Signature_
    
    
    def timeline_set_clip_duration(clip: Clip, duration: int) -> None
    

_Description_

    Sets the duration of clip to duration frames.

_Example_
    
    
    bar_kf = cmd.timeline_add_keyframe(track=foo_track, frame_num=42)
    cmd.timeline_delete_keyframe(keyframe=bar_kf)
    

_See Also_

    

timeline_set_clip_speed, timeline_move_clip_to

* * *

timeline_set_clip_speed

_Signature_
    
    
    def timeline_set_clip_speed(clip: Clip, duration: int) -> None
    

_Description_

    Sets the duration of clip to duration frames.

_Example_
    
    
    bar_kf = cmd.timeline_set_clip_speed(track=foo_track, frame_num=42)
    

_See Also_

    

timeline_set_clip_duration, timeline_move_clip_to

* * *

timeline_register_program

_Signature_
    
    
    def timeline_register_program(program_name: str, program_type: object) -> None
    

_Description_

    Registers the program with program_type type and program_name name to an internal program database. The types are registered so that instances of it can be generated and represented as animation clips on the timeline. After a program is registered, you can store a program using timeline_store_program_at using program_name (case-sensitive) to identify the program.

_Example_
    
    
    class Drift:
        def animate(t: float) -> None:
            # ...
    
    cmd.timeline_register_program(program_name="drift", program_type=Drift)
    

_See Also_

    

timeline_store_program_at

* * *

timeline_store_program_at

_Signature_
    
    
    def timeline_store_program_at(track: Track, frame_num: int, duration: int, program: Union[object, str]) -> Clip
    

_Description_

    Adds an animation clip described by program onto track at frame frame_num for duration frames.

_Example_
    
    
    cmd.timeline_register_program(program_name="drift", program_type=Drift)
    cam_track = myComp.get_main_camera_track()
    cmd.timeline_store_program_at(track=cam_track, frame_num=60, duration=360, program="drift")
    

_See Also_

    

timeline_register_program

_Notes_

    To customize the properties of a program, instantiate a program, edit its properties, and provide it to the function. Example:
    
    
    roll_prg= pymol.timeline_programs.Roll(cam_track)
    roll_prg.num_loops = 2
    clip = cmd.timeline_store_program_at(track=cam_track, frame_num=60, duration=360, program=roll_prg)
    

    Since the program itself is agnostic to the duration of the clip that it’s stored in, in order to change the duration of the program, change the clip’s duration instead via timeline_set_clip_duration.

* * *

timeline_produce (NOT YET SUPPORTED)

_Signature_
    
    
    def timeline_produce(comp: Composition, render_mode: str = "draw", first_frame: int = None, last_frame: int = None, encoder: str = None, width: int = None, height: int = None, quality: int=100, quiet: bool=False)
    

[comp: Composition, render_mode: str = "draw", first_frame: int = None, last_frame: int = None, encoder: str = None, width: int = None, height: int = None, quality: int=100, quiet=False] 

_Description_

    Exports a comp to a playable movie (file format determined by encoder).

_Example_
    
    
    cmd.movie.timeline_produce(comp=comp)
    

Retrieved from "[https://pymolwiki.org/index.php?title=Timeline_Python_API&oldid=13874](https://pymolwiki.org/index.php?title=Timeline_Python_API&oldid=13874)"


---

