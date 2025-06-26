# Category: Using Pymol

## Mouse Controls

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## MacPyMol And PyMolX11 Hybrid

We strongly recommend (require?) use of a three-button mouse with MacPyMOL and PyMOLX11Hybrid. However, it is possible to use MacPyMOL in a limited way on Macs that have a single button mouse thanks to some built-in mouse remapping in the Mac OS X GLUT implementation. Her is how that works... 

If the Mac is hooked up to a 1-button mouse (only) when MacPyMOL or HybridX11PyMOL are launched, then Mac OS X itself will furnish translations as follows: 
    
    
    Rotate:       Click & Drag
    XY-Translate: Option-Click & Dreg
    Zoom:         Control-Click & Drag
    Select:       Click & Release
    Box-Select:   Shift-Click & Drag
    Box-Deselect: Shift-Option-Click & Drag
    Clipping:    Control-Shift-Click & Drag with...
     * near plane controlled by vertical motion
     * far plane controlled by horizontal motion
    

Note that not all mouse actions are possible with a one-button mouse. For example, I don't think molecular editing is possible with a one-button mouse. 

# See Also

[Configuring the Mouse](/index.php/Config_Mouse "Config Mouse")

Retrieved from "[https://pymolwiki.org/index.php?title=Mouse_Controls&oldid=8663](https://pymolwiki.org/index.php?title=Mouse_Controls&oldid=8663)"


---

## Mouse Settings

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Retrieved from "[https://pymolwiki.org/index.php?title=Mouse_Settings&oldid=6520](https://pymolwiki.org/index.php?title=Mouse_Settings&oldid=6520)"


---

## Ray Tracing

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Adjusting ray trace-image size

Much, much more at [Ray](/index.php/Ray "Ray"). 

The pymol ray tracer can generate an image of any size. 
    
    
    ray height,width
    

Example: 
    
    
    ray 3000,2400
    png filename.png
    

## Ray tracing maps

For better quality maps with a white background. 
    
    
    set ray_trace_fog,0
    set ray_shadows,0
    set antialias,1
    ray 1600,1200
    png img.png
    

(This will take quite a while...) 

## CGO label orientation

You could use the cmd.rotate and cmd.translate to position the labels, but it is likely to be somewhat painful. If I'm not mistaken, the rotation will always be about and axis through the origin and so you may need to translate the label into the correct position. 

Thus if you have your label in an object called 'text', you could do, 
    
    
      cmd.rotate(axis='x',angle=20.,object='text')
    

and repeat this with different angles, until you get the orientation correct. Then use: 
    
    
      cmd.translate(vector='[1.,2.,3.]',object='text')
    

(using the appropriate vector, of course!) to position the label. Not ideal, but if it is sufficiently important, it can be done! 

Retrieved from "[https://pymolwiki.org/index.php?title=Ray_Tracing&oldid=5549](https://pymolwiki.org/index.php?title=Ray_Tracing&oldid=5549)"


---

## Category:Objects and Selections

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Information About Objects and Selections in PyMol. Easily one of the most important sections of information. 

## Hiding Selection Names

To hide a selection name in the GUI on the right (the one that responds to "set internal_gui") just name your selections with leading underscores. Example: ` name ` versus ` _name `

## Subcategories

This category has the following 4 subcategories, out of 4 total. 

### R

  * [Representations](/index.php/Category:Representations "Category:Representations")



### S

  * [Selector Quick Reference](/index.php/Category:Selector_Quick_Reference "Category:Selector Quick Reference")



### W

  * [Working with Objects](/index.php/Category:Working_with_Objects "Category:Working with Objects")
  * [Working with Selections](/index.php/Category:Working_with_Selections "Category:Working with Selections")



## Pages in category "Objects and Selections"

The following 4 pages are in this category, out of 4 total. 

### C

  * [Color](/index.php/Color "Color")



### D

  * [Displaying Biochemical Properties](/index.php/Displaying_Biochemical_Properties "Displaying Biochemical Properties")
  * [Symmetry](/index.php/Symmetry "Symmetry")



### N

  * [Nonstandard Amino Acids](/index.php/Nonstandard_Amino_Acids "Nonstandard Amino Acids")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Objects_and_Selections&oldid=3156](https://pymolwiki.org/index.php?title=Category:Objects_and_Selections&oldid=3156)"


---

