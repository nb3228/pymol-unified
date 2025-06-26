# Category: OutputFormats

## COLLADA

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

**COLLADA** is an XML-based open 3D graphics format, designed to be used as an interchange format among various graphics programs. Beginning with version 1.7.3.2 (SVN rev4097) in the Open Source codebase, PyMOL is able to export COLLADA files, which have a ".dae" extension. 

## Contents

  * 1 Details
    * 1.1 USAGE
    * 1.2 Features
    * 1.3 Related Settings
    * 1.4 Limitations
    * 1.5 Increasing the Quality



# Details

COLLADA format defines a scene, with a camera and various geometric elements, which are colored and lit according to the properties of the material the object is made from and any effects (e.g. shaders) applied to it. PyMOL exports the scene exactly as it looks in the viewport when the save command is executed, so the orientations and colors of the objects should be the same. It currently saves every enabled object in the scene, regardless of any selection passed to the `save` command. 

### USAGE
    
    
    save file.dae
    

  


## Features

  * Full support of all standard representations (cartoon, stick, line, sphere, surface, etc.)
  * Customizable output quality and file size, based on quality settings for active representations (e.g. [surface_quality](/index.php/Surface_quality "Surface quality"))
  * Transparency support
  * As of PyMOL 2.4: support for exporting glTF files using the COLLADA2GLTF tool



## Related Settings

  * [collada_export_lighting](/index.php?title=Collada_export_lighting&action=edit&redlink=1 "Collada export lighting \(page does not exist\)") \- Optionally include PyMOL's lighting information.
  * [collada_geometry_mode](/index.php/Collada_geometry_mode "Collada geometry mode") \- Specify how geometry elements are presented.
  * [geometry_export_mode](/index.php/Geometry_export_mode "Geometry export mode") \- Exclude camera as well as lighting information.



  


## Limitations

  * **Shaders.** In version 1.4 of the COLLADA specification, only vertex-based meshes are defined, so all regular geometric shapes (e.g. spheres, cylinders, cones) have to be converted into a series of triangle meshes. This produces output similar to what you would see in PyMOL with `use_shaders` turned off. (Version 1.5 of the COLLADA spec includes "boundary representation" (<brep>) elements that define such geometric shapes as you might expect: a sphere is defined by a point and a radius, a cylinder by two points and a radius, a cone by two points and two radii. However, apparently due to the complexity of the physics portion of the 1.5 specification, few programs have even attempted to adopt it yet, so we're stuck with old-school triangle meshes.)


  * **Lighting.** Whereas in PyMOL, the lighting of the object appears similar no matter how the object is rotated--that is, the lights move with the camera or viewport--in COLLADA, the lights are fixed in the scene with the objects, and therefore, with PyMOL `[light_count](/index.php/Light_count "Light count")` of 2 or higher, only one side of the objects is lit with the directional lights; the other side appears dark. For this reason, PyMOL's lights are excluded from COLLADA output files by default. To include lighting in the output file, turn on the `[collada_export_lighting](/index.php?title=Collada_export_lighting&action=edit&redlink=1 "Collada export lighting \(page does not exist\)")` setting.


  * **Labels.** Like lights, labels would be fixed with the objects, so from all but a narrow range of orientations, the text would be unreadable.


  * **Volumes.** PyMOL's volume representation consists of a series of evenly spaced planes perpendicular to the viewport, the orientation of which is changed whenever the view is changed. Unfortunately, this makes them essentially useless in a COLLADA scene as well, as, once the original orientation is changed, the planes would look less like a 3D volume projection, and more like a stack of flat images slicing through the rest of the scene.



## Increasing the Quality

You can increase the number of triangles exported, and therefore the smoothness of the resulting objects, by increasing the following settings. 

  * surface_quality
  * sphere_quality
  * stick_quality
  * cartoon_sampling



Conversely, to decrease the file size, decrease the settings. 

Retrieved from "[https://pymolwiki.org/index.php?title=COLLADA&oldid=13178](https://pymolwiki.org/index.php?title=COLLADA&oldid=13178)"


---

## Vrml

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

PyMOL can export wireframe/VRML models to file. You can most easily do this through [exporting the mesh/surface](/index.php/Surface#Exporting_Surface.2FMesh_Coordinates_to_File "Surface"). 

  


## Increasing the Quality

You can increase the number of triangles exported by: 
    
    
    set surface_quality, 1
    set cartoon_sampling, 20
    

Retrieved from "[https://pymolwiki.org/index.php?title=Vrml&oldid=5094](https://pymolwiki.org/index.php?title=Vrml&oldid=5094)"


---

