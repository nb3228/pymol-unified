# Category: 3D

## Stereo mode

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/9/93/Anaglyph1.png)](/index.php/File:Anaglyph1.png)

[](/index.php/File:Anaglyph1.png "Enlarge")

Example of Anaglyph 3D in PyMOL

[![](/images/6/63/Ana4_focal_ray.png)](/index.php/File:Ana4_focal_ray.png)

[](/index.php/File:Ana4_focal_ray.png "Enlarge")

Example of a ray traced, focal blurred, anaglyph image.

[![](/images/8/82/3D-B-DNA.png)](/index.php/File:3D-B-DNA.png)

[](/index.php/File:3D-B-DNA.png "Enlarge")

3D B-DNA Example

[![](/images/2/21/Amyloid.png)](/index.php/File:Amyloid.png)

[](/index.php/File:Amyloid.png "Enlarge")

Amyloid beta example

The [stereo_mode](/index.php/Stereo_mode "Stereo mode") setting sets the type of stereo mode, if the [stereo](/index.php/Stereo "Stereo") setting is enabled. 

_You can also control both settings with the[stereo](/index.php/Stereo "Stereo") command, which is more convenient!_

## Contents

  * 1 Syntax
  * 2 Supported Stereo Modes
  * 3 Notes
    * 3.1 Anaglyph Color Quality and Ghosting
  * 4 See Also



## Syntax
    
    
    set stereo_mode, integer
    

Valid values for the **integer** argument are listed in the following table. 

## Supported Stereo Modes

Corresponding keyword arguments (instead of numeric values) can be passed to the [stereo](/index.php/Stereo "Stereo") command. 

value | description   
---|---  
1 | **quad-buffered**  
2 | **cross-eyed**  
3 | **walleye**  
4 | **geowall**  
5 | **sidebyside**  
6 | **stencil by row** , Zalman   
7 | **stencil by col**  
8 | **stencil checkerboard**  
9 | **stencil custom** for developers   
10 | **anaglyph** (requires green/magenta glasses)   
11 | **dynamic polarization**  
12 | **clone dynamic**  
  
## Notes

### Anaglyph Color Quality and Ghosting

To test the quality of your glasses and coloring, you can test for "ghosting". Choose a part of the structure where the channels for the left and right eye are nicely separated in space. Hold one lens over the part of the structure and see you see one image or if there is a faint second image for the other channel. Then hold the other other lens over the image and see if there is a faint second image. 

## See Also

  * [Stereo](/index.php/Stereo "Stereo")
  * [Stereo_angle](/index.php/Stereo_angle "Stereo angle")



Retrieved from "[https://pymolwiki.org/index.php?title=Stereo_Mode&oldid=11520](https://pymolwiki.org/index.php?title=Stereo_Mode&oldid=11520)"


---

