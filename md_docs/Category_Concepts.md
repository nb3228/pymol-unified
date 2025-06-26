# Category: Concepts

## Model Space and Camera Space

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Objects are defined in a Cartesian coordinate system, i.e. as xyz coordinates for each atom (or surface point or corner or...). This coordinate system is model space. To make objects visible, PyMOL places a "camera" such that it sees most of the molecule. The picture that is shown on the screen is the image that is "taken" by the camera. Initially, after loading an object, the camera is placed such that it looks on the object parallel to the model's z axis. The model's x and y axes are horizontal and vertical, respectively. So, initially, the model axes correspond to the physical directions of the screen: x and y are horizontal and vertical, z is perpendicular to the screen. 

In normal operation, to change the view, the camera is moved and the view of all objects is changed simultaneously; that is, the objects are not moved relative to one other. After changing the view, the x, y and z axes of model space no longer correspond to the directions of the screen. The commands [turn](/index.php/Turn "Turn"), [move](/index.php/Move "Move"), [center](/index.php/Center "Center"), [zoom](/index.php/Zoom "Zoom"), and [orient](/index.php/Orient "Orient"), as well as mouse action in viewing mode, can be used to move the camera in object space. 

For easy description of object movement relative to the camera, an additional coordinate system is defined: Camera space. The camera is situated in the origin of camera space, x and y are horizontal respectively vertical and z is the viewing direction. Camera space always corresponds to the physical directions of the screen--x and y are horizontal and vertical, and z is perpendicular to the screen. The relation between camera space and model space is described in the view matrix (see [get_view](/index.php/Get_view "Get view")). 

Changing the coordinates in model space is only necessary when objects must be moved relative to each other or when the transformed coordinates are to be written into a file. Model space coordinates can be changed with the commands [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), and [transform_selection](/index.php/Transform_selection "Transform selection"), or with the mouse in editing mode. 

  


## Space and File saving

When model coordinates are [saved](/index.php/Save "Save") (e.g. as a pdb file), PyMOL uses model space coordinates. In contrast, when a specific view is saved as [PovRay](/index.php/PovRay "PovRay") or [vrml](/index.php/Vrml "Vrml"), camera space coordinates are used (unless [geometry_export_mode](/index.php/Geometry_export_mode "Geometry export mode") is turned on, in which case no camera or lighting information is included in the output file). 

When image files in model space coordinates are needed, first choose colors and representations as needed, then align camera space to model space: 
    
    
    set_view (\
         1,    0,    0,\
         0,    1,    0,\
         0,    0,    1,\
         0,    0,    0,\
         0,    0,    0,\
         0,    300,  1 )
    

After setting the view, the molecule will probably be out of sight but will be saved to the file nonetheless. This can be used e.g. to receive a list of surface points in model space. 

## See Also

[turn](/index.php/Turn "Turn"), [move](/index.php/Move "Move"), [center](/index.php/Center "Center"), [zoom](/index.php/Zoom "Zoom"), [orient](/index.php/Orient "Orient"), [rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [Transform_selection](/index.php/Transform_selection "Transform selection"), [Get_View](/index.php/Get_View "Get View")

Retrieved from "[https://pymolwiki.org/index.php?title=Model_Space_and_Camera_Space&oldid=11529](https://pymolwiki.org/index.php?title=Model_Space_and_Camera_Space&oldid=11529)"


---

## Object Matrix

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Every object has an object matrix associated with it. The object matrix tells PyMOL whether the object coordinates are to be transformed before the object is displayed. The matrix is not a standard matrix in mathematical sense, it is something PyMOL-specific, also called TTT matrix: It is 4X4, with the upper left 3x3 forming a rotation matrix, the fourth column and row representing pre-rotation and post-rotation translation vectors respectively, and the 16th element always being 1.0. TTT stands for translation - transformation - translation. 

Initially, the object matrix is the identity matrix, so no transformation takes place. It is changed during structure alignment (e.g. align and super). By hand, it can be changed with the commands rotate, translate, and origin. 

It can be copied to other objects ([matrix_copy](/index.php/Matrix_copy "Matrix copy")) and reset to identity ([matrix_reset](/index.php/Matrix_reset "Matrix reset")). With [get_object_matrix](/index.php/Get_object_matrix "Get object matrix"), the current matrix can be read out. 

## See also

[rotate](/index.php/Rotate "Rotate"), [translate](/index.php/Translate "Translate"), [origin](/index.php/Origin "Origin"), [get_object_matrix](/index.php/Get_object_matrix "Get object matrix"), [matrix_copy](/index.php/Matrix_copy "Matrix copy"), [matrix_reset](/index.php/Matrix_reset "Matrix reset"), [align](/index.php/Align "Align"), [super](/index.php/Super "Super")

Retrieved from "[https://pymolwiki.org/index.php?title=Object_Matrix&oldid=8447](https://pymolwiki.org/index.php?title=Object_Matrix&oldid=8447)"


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

