# Category: Volume

## Volume

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

[![](/images/b/b1/1okyVol.png)](/index.php/File:1okyVol.png)

[](/index.php/File:1okyVol.png "Enlarge")

Volume visualization of electron density for PDB 1oky

[![](/images/5/5d/1okyVolPanel.png)](/index.php/File:1okyVolPanel.png)

[](/index.php/File:1okyVolPanel.png "Enlarge")

Volume panel for the 1oky volume example. It has the iso-levels on the x-axis and the opacity (1.0 - transparency) on the y-axis.

  
Volume creates a new volume object from a map object. The data (3D scalar fields) are shown as a true 3D object using coloring and transparencies defined by the user to illustrate the data values. This technique supports single and multiple isosurfaces. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Example
  * 4 Screencasts
  * 5 Changes with PyMOL Version
  * 6 Ray Tracing
  * 7 Known Limitations
  * 8 See Also



## Usage
    
    
    volume name, map [, ramp [, selection [, buffer [, state [, carve ]]]]]
    

## Arguments

  * name = the name for the new volume object.
  * map = the name of the map object to use for computing the volume.
  * ramp = str: named color ramp {default: }
  * selection = an atom selection about which to display the mesh with an additional "buffer" (if provided).
  * carve = a radius about each atom in the selection for which to include density. If "carve" is not provided, then the whole brick is displayed.



## Example
    
    
    fetch 1oky, type=2fofc, async=0
    volume 1okyVol, 1oky_2fofc
    

## Screencasts

  * [Silent demo movie](http://www.youtube.com/watch?v=tuAo_8-_HIc) showing the basics of loading and using a volume in PyMOL. There are more capabilities, but this is the basic functionality.



## Changes with PyMOL Version

  * 1.4.0: first version with volume support
  * 1.7.2: 
    * pre-integrated volume rendering (volume_mode=1) as Incentive-PyMOL-only feature.
    * scripting support with custom color ramp ([volume_color](/index.php/Volume_color "Volume color"), [volume_ramp_new](/index.php?title=Volume_ramp_new&action=edit&redlink=1 "Volume ramp new \(page does not exist\)"))
    * improved volume panel, panel can be opened from the object menu ("C > panel")
    * lots of bugs fixed



## Ray Tracing

**There is no actual ray tracing support**. The volume rendering is implemented exclusively with OpenGL shaders. The recommended way to render a high resolution image is to use the [draw](/index.php/Draw "Draw") command. Example: 
    
    
    # render high resolution image on screen
    draw 4000, 3000, antialias=2
    png highres.png
    

Ray trace specific features like shadows or outlines ([ray_trace_mode](/index.php/Ray_trace_mode "Ray trace mode")) are not available with this approach. If such features are needed, the [ray_volume](/index.php?title=Ray_volume&action=edit&redlink=1 "Ray volume \(page does not exist\)") setting activates a hybrid solution, which blends the OpenGL rendered volume with the ray traced non-volume objects. Image sizes other than the current window size are not possible. 
    
    
    # compose on-screen volume with ray traced image
    set ray_volume
    ray
    png composed.png
    

Neither of these two solutions work with [headless (batch) mode](/index.php/Launching_PyMOL#Running_PyMOL_in_batch_mode "Launching PyMOL"). 

## Known Limitations

  * No real ray-tracing support yet
  * Multiple volume objects don't blend properly



## See Also

  * <http://pymol.org/volume>
  * <http://pymol.org/d/media:volumevisualization>
  * [volume_color](/index.php/Volume_color "Volume color")
  * [volume_ramp_new](/index.php?title=Volume_ramp_new&action=edit&redlink=1 "Volume ramp new \(page does not exist\)")
  * [map_new](/index.php/Map_new "Map new")
  * [isomesh](/index.php/Isomesh "Isomesh")
  * [isosurface](/index.php/Isosurface "Isosurface")



Retrieved from "[https://pymolwiki.org/index.php?title=Volume&oldid=12536](https://pymolwiki.org/index.php?title=Volume&oldid=12536)"


---

## Volume color

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

volume_color set or get the [volume](/index.php/Volume "Volume") colors. 

## Contents

  * 1 Usage
  * 2 Arguments
  * 3 Examples
  * 4 See Also



## Usage
    
    
    volume_color name [, ramp ]
    

## Arguments

  * name = str: volume object name
  * ramp = str, list or empty: named ramp, space delimited string or list with (x, color, alpha, ...) or (x, r, g, b, alpha, ...) values. If empty, get the current volume colors.



## Examples

Setting volume colors with one line: 
    
    
    fetch 1a00, map, type=2fofc, async=0
    volume vol, map
    volume_color vol, .8 cyan 0. 1. blue .3 2. yellow .3 
    

Using a named color ramp: 
    
    
    volume_ramp_new cyanblueyellow, \
        .8 cyan 0. \
        1. blue .3 \
        2. yellow .3
    volume_color vol, cyanblueyellow
    

Getting the current volume ramp: 
    
    
    PyMOL>volume_color vol
    ### cut below here and paste into script ###
    cmd.volume_ramp_new('ramp399', [\
         0.80, 0.00, 1.00, 1.00, 0.00, \
         1.00, 0.00, 0.00, 1.00, 0.30, \
         2.00, 1.00, 1.00, 0.00, 0.30, \
       ])
    ### cut above here and paste into script ###
    

## See Also

  * [volume](/index.php/Volume "Volume")



Retrieved from "[https://pymolwiki.org/index.php?title=Volume_color&oldid=12729](https://pymolwiki.org/index.php?title=Volume_color&oldid=12729)"


---

