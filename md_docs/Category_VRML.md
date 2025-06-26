# Category: VRML

## Surface

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

The surface representation of a protein, in PyMol, shows the ["Connolly" surface](http://en.wikipedia.org/wiki/Connolly_surface) or the surface that would be traced out by the **surfaces** of waters in contact with the protein at all possible positions. 

[![](/images/7/79/Surface_ex.png)](/index.php/File:Surface_ex.png)

[](/index.php/File:Surface_ex.png "Enlarge")

Surface Representation Example

  


## Contents

  * 1 Enabling
  * 2 Settings
    * 2.1 Examples
      * 2.1.1 Transparency
      * 2.1.2 Quality
      * 2.1.3 Probe Radius
  * 3 Tips
    * 3.1 Exporting Surface/Mesh Coordinates to File
      * 3.1.1 Older PyMOL Versions
      * 3.1.2 Newer PyMOL Versions
    * 3.2 Representation-independent Color Control
    * 3.3 Displaying a protein as surface with a ligand as sticks
    * 3.4 Calculating a partial surface
    * 3.5 Displaying surface inside a molecule
    * 3.6 Creating a Closed Surface
    * 3.7 Smooth surface quick (blob)
    * 3.8 Smooth surface accurate (blob)
    * 3.9 Huge Surfaces
  * 4 Performance



## Enabling

To enable the surface representation do 
    
    
    show surface, SEL
    

for any proper selection SEL. 

## Settings

  * [cavity_cull](/index.php/Cavity_cull "Cavity cull")
  * [surface_best](/index.php?title=Surface_best&action=edit&redlink=1 "Surface best \(page does not exist\)")
  * [surface_negative_color](/index.php/Surface_negative_color "Surface negative color")
  * [surface_carve_cutoff](/index.php/Surface_carve_cutoff "Surface carve cutoff")
  * [surface_negative_visible](/index.php/Surface_negative_visible "Surface negative visible")
  * [surface_carve_normal_cutoff](/index.php/Surface_carve_normal_cutoff "Surface carve normal cutoff")
  * [surface_normal](/index.php?title=Surface_normal&action=edit&redlink=1 "Surface normal \(page does not exist\)")
  * [surface_carve_selection](/index.php/Surface_carve_selection "Surface carve selection")
  * [surface_optimize_subsets](/index.php?title=Surface_optimize_subsets&action=edit&redlink=1 "Surface optimize subsets \(page does not exist\)")
  * [surface_carve_state](/index.php?title=Surface_carve_state&action=edit&redlink=1 "Surface carve state \(page does not exist\)")
  * [surface_poor](/index.php?title=Surface_poor&action=edit&redlink=1 "Surface poor \(page does not exist\)")
  * [surface_circumscribe](/index.php?title=Surface_circumscribe&action=edit&redlink=1 "Surface circumscribe \(page does not exist\)")
  * [surface_proximity](/index.php/Surface_proximity "Surface proximity")
  * [surface_clear_cutoff](/index.php?title=Surface_clear_cutoff&action=edit&redlink=1 "Surface clear cutoff \(page does not exist\)")
  * [surface_quality](/index.php/Surface_quality "Surface quality")
  * [surface_clear_selection](/index.php?title=Surface_clear_selection&action=edit&redlink=1 "Surface clear selection \(page does not exist\)")
  * [surface_ramp_above_mode](/index.php/Surface_ramp_above_mode "Surface ramp above mode")
  * [surface_clear_state](/index.php?title=Surface_clear_state&action=edit&redlink=1 "Surface clear state \(page does not exist\)")
  * [surface_solvent](/index.php/Surface_solvent "Surface solvent")
  * [surface_color](/index.php/Surface_color "Surface color")
  * [surface_trim_cutoff](/index.php?title=Surface_trim_cutoff&action=edit&redlink=1 "Surface trim cutoff \(page does not exist\)")
  * [surface_debug](/index.php?title=Surface_debug&action=edit&redlink=1 "Surface debug \(page does not exist\)")
  * [surface_trim_factor](/index.php?title=Surface_trim_factor&action=edit&redlink=1 "Surface trim factor \(page does not exist\)")
  * [surface_miserable](/index.php?title=Surface_miserable&action=edit&redlink=1 "Surface miserable \(page does not exist\)")
  * [surface_type](/index.php/Surface_type "Surface type")
  * [surface_mode](/index.php/Surface_mode "Surface mode")



### Examples

#### Transparency

To adjust the transparency of surfaces try: 
    
    
    set transparency, 0.5
    

Where 1.0 will be an invisible and 0.0 a completely solid surface. 

#### Quality

To smooth your surface representation try: 
    
    
    set surface_quality, 1
    

or higher if you wish, though it will take longer and might look odd. 

#### Probe Radius

To change the probe radius other than default 1.4 Å, you need to change the solvent radius, say, 1.6 Å: 
    
    
    set solvent_radius, 1.6
    

If the surface does not change correspondingly, use: 
    
    
    rebuild
    

## Tips

### Exporting Surface/Mesh Coordinates to File

PyMOL can export its coordinates as WRL wireframe model files for VRML input. 

#### Older PyMOL Versions
    
    
    # export the coordinates to povray
    open("surface.inc","w").write(cmd.get_povray()[1])
    

#### Newer PyMOL Versions
    
    
    # export the coordinates to .wrl file
    save myscene.wrl
    

or 
    
    
    # export the coordinates to .obj file. Only surface representation can be saved as .obj.
    # NOTE: the coordinates are saved in camera coordinate system.
    save myscene.obj
    

### Representation-independent Color Control

To color the surface representation a different color than the underlying cartoon or ligand representations, simply duplicate the object, show only the surface in the duplicate, and show only the cartoon and/or ligands in the original object. 

Or use the [surface_color](/index.php/Surface_color "Surface color") setting that is available. 

  


[![](/images/5/51/Representation_independent_color_control.jpg)](/index.php/File:Representation_independent_color_control.jpg)

[](/index.php/File:Representation_independent_color_control.jpg "Enlarge")

Representation-independent Color Control Example

### Displaying a protein as surface with a ligand as sticks

An easy way to do this is to create separate objects for each type of display. 

1 Load your protein 

2 Select the ligand 

3 Create a separate object for the ligand 

4 Remove ligand atoms from the protein 

5 Display both objects separately 

Example: 
    
    
    load prot.ent,protein
    select ligand,resn FAD
    create lig_sticks,ligand
    remove ligand
    show sticks,lig_sticks
    show surface,protein
    

  
Even easier is to: 

1 Load the protein 

2 S (Show) > organic > stick 

3 S (Show) > surface 

### Calculating a partial surface

There is, until now, an undocumented way to calculate a surface for only a part of an object without creating a new one: 
    
    
    flag ignore, not A/49-63/, set
    delete indicate
    show surface
    

If the surface was already computed, then you'll also need to issue the command: 
    
    
    rebuild
    

See [Get_Area](/index.php/Get_Area "Get Area") for more information on surface area calculations. 

### Displaying surface inside a molecule

As far as I can tell, setting ambient to zero alone doesn't quite do the job, since some triangles still get lit by the light source. The best combination I can find is: 
    
    
    set ambient=0
    set direct=0.7
    set reflect=0.0
    set backface_cull=0
    

Which gives no shadows and only a few artifacts. 

As an alternative, you might just consider showing the inside of the surface directly...that will create less visual artifacts, and so long as ambient and direct are sufficiently low, it will look reasonable in "ray". 
    
    
    util.ray_shadows("heavy")
    set two_sided_lighting=1
    set backface_cull=0
    

### Creating a Closed Surface

  * [![Example OPEN Surface](/images/6/6a/Surface_open.png)](/index.php/File:Surface_open.png "Example OPEN Surface")

Example OPEN Surface 

  * [![Example CLOSED Surface](/images/f/f5/Surface_closed.png)](/index.php/File:Surface_closed.png "Example CLOSED Surface")

Example CLOSED Surface 




To create what I'll call a **closed surface** (see images), you need to first make your atom selections, then create a new object for that selection then show the surface for that object. Here's an example. 
    
    
     sel A, id 1-100
     create B, A
     show surface, B
    

### Smooth surface quick (blob)

To get a quick blob type surface (not as accurate): 
    
    
    set solvent_radius, 4   
    alter all, vdw=4 
    sort
    set surface_quality, 1
    

### Smooth surface accurate (blob)

To get an accurate blob type surface: 
    
    
    set surface_quality, 1
    alter all, b=50
    alter all, q=1
    set gaussian_resolution,5
    map_new mapA, gaussian, 1, sele or pdb, 6
    isosurface surfA, mapA
    

**Notes:** Set gaussian resolution is variable with a larger number causing a more smooth surface (4 is medium and 8 is very smooth). The temperature factor field (b) has at least as much impact as the resolution setting, so increasing b factors is the more computationally efficient way of increasing the blur effect. If you are displaying more then one surface in a .pse file you must create a new map for each one (if you have three you will create mapA for the first, mapB for the second and mapC for the third), then you apply an isosurface to each map (isosurface surfA, mapA - isosurface surfB, mapB - isosurface surfC, mapC). 

### Huge Surfaces

If your protein or complex is too large to render ([ray](/index.php/Ray "Ray") runs out of RAM, for example) then check out the [tip for huge surfaces](/index.php/Huge_surfaces "Huge surfaces"). 

## Performance

To optimize performance and responsiveness, PyMOL tends to defer compute-intensive tasks until their results are actually needed. Thus, 
    
    
    cmd.show("surface")
    

doesn't actually show a surface, it only sets the surface visibility flag on the atoms present (for future reference). An actual surface won't be computed until PyMOL is asked to refresh or render the display. When running a script, you can force an update by calling: 
    
    
    cmd.refresh()
    

after cmd.show. 

Retrieved from "[https://pymolwiki.org/index.php?title=Surface&oldid=12463](https://pymolwiki.org/index.php?title=Surface&oldid=12463)"


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

