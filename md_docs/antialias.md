# antialias

antialias

### Table of Contents

  * Antialiasing

    *       * Antialiasing in Ray-tracing

      * Real-time Antialiasing




# Antialiasing

Anti-aliasing is a graphical technique for smoothing pixel-based images to reduce discrete edges. 

[![No Anti-aliasing](/dokuwiki/lib/exe/fetch.php?w=200&h=200&tok=a3a51a&media=no_aa.png)](/dokuwiki/lib/exe/detail.php?id=antialias&media=no_aa.png "no_aa.png") |  [![Anti-aliased](/dokuwiki/lib/exe/fetch.php?w=200&h=200&tok=9860d9&media=aa.png)](/dokuwiki/lib/exe/detail.php?id=antialias&media=aa.png "aa.png")  
---|---  
No Antialiasing (rough edges) |  Anti-aliased Image (smooth edges)   
  
These settings (both **antialias** for ray-tracing and **antialias_shader** for real-time antialiasing) can be changed using the menu : 

Setting → Rendering → Antialias 

[![AA Menu](/dokuwiki/lib/exe/fetch.php?w=250&h=250&tok=94a91e&media=aa_menu.png)](/dokuwiki/lib/exe/detail.php?id=antialias&media=aa_menu.png "aa_menu.png")

or setting these settings: 
    
    
    # set real-time antialiasing on, using SMAA  (off by default)
    set antialias_shader, 2
    
    # set ray-tracing antialias level to 2 (1 by default)
    set antialias, 2

### Antialiasing in Ray-tracing

Antialiasing is enabled by default for ray-tracing. The **antialias** setting controls the level of antialiasing. Higher numbers take longer to ray trace but provide much smoother and better looking images (range: 0-4, 0 is off, 4 is a higher setting, better image quality. default: 1). 

### Real-time Antialiasing

As of PyMOL 1.7.2, we now provide real-time rendering when shaders are available and your graphics card handle these shaders, using the **antialias_shader** setting: 

  * 0 = no anti-aliasing (default)

  * 1 = fast approximate anti-aliasing (FXAA 1-pass)

  * 2 = subpixel morphological anti-aliasing (SMAA 3-passes)

  * 3-8 = debugging/different parameters for SMAA 

    * 3 = edges/1st-pass

    * 4 = blending values/2nd-pass

    * 5,6,7,2/8 = SMAA presets/parameters (low,medium,high,ultra)




There are differences between FXAA and SMAA, most notably that SMAA has 3-passes and FXAA has only 1-pass, depending on your graphics card, SMAA might perform slower. Also, currently we do these Anti-alias passes after everything in the scene is rendered (after the labels but before the indicators). FXAA degrades labels slightly worse than SMAA. In the next version of PyMOL, we will make sure these anti-alias passes happen before labels are rendered as well, so this will not be an issue. 

[![NoAntialiasing](/dokuwiki/lib/exe/fetch.php?w=200&h=200&tok=254bfd&media=no_aa_labels.png)](/dokuwiki/lib/exe/detail.php?id=antialias&media=no_aa_labels.png "no_aa_labels.png")| [![FXAA](/dokuwiki/lib/exe/fetch.php?w=200&h=200&tok=eb8d1e&media=fxaa_labels.png)](/dokuwiki/lib/exe/detail.php?id=antialias&media=fxaa_labels.png "fxaa_labels.png")| [![SMAA](/dokuwiki/lib/exe/fetch.php?w=200&h=200&tok=265971&media=smaa_labels.png)](/dokuwiki/lib/exe/detail.php?id=antialias&media=smaa_labels.png "smaa_labels.png")  
---|---|---  
No Antialiasing  
(Good labels, bad edges) |  FXAA   
(Better edges, worse labels) |  SMAA  
(Good edges, good labels)   
  
There are subtle differences between these pictures, however, when you look closely, they are clearly different. Look at the pixels in “LEU”: the FXAA image has the label “bleeding” into other pixels, while the SMAA labels look much closer to the original without AA. In this case, the edges in SMAA seem notably better as well. 

For most cases, you should use SMAA (i.e., **antialias_shader** =2), unless you are not using labels and are seeing a notable performance difference. 

antialias.txt · Last modified: 2014/06/27 20:39 by holder
