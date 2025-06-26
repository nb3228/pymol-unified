# rendering

rendering

### Table of Contents

  * Rendering

    * Draw vs. Ray

    * General Topics

    * Related Settings

    * Special Cases




# Rendering

This section covers techniques for generating high-quality image output. 

## Draw vs. Ray

The [draw](/dokuwiki/doku.php?id=command:draw "command:draw") command enables you to quickly create an antialiased image of the scene with current [representations](/dokuwiki/doku.php?id=representation "representation") using OpenGL. However, when using [draw](/dokuwiki/doku.php?id=command:draw "command:draw"), lighting and colors are interpolated, analytical geometries are tessellated as triangles, there is no way of displaying shadows, and you cannot render cross-sections. 

The [ray](/dokuwiki/doku.php?id=command:ray "command:ray") command use a raytracer to perform per-pixel lighting and coloring, to render analytical geometries with per-pixel accuracy, to display shadows, and to render cross sections. However, such calculations can take a considerable amount of time. 

Although the image quality obtained from [draw](/dokuwiki/doku.php?id=command:draw "command:draw") can be [improved](/dokuwiki/doku.php?id=representation:improve "representation:improve") by adjusting various settings, it is usually easier to just use [ray](/dokuwiki/doku.php?id=command:ray "command:ray"). 

## General Topics

[Image Size](/dokuwiki/doku.php?id=rendering:size "rendering:size")

## Related Settings

[Camera](/dokuwiki/doku.php?id=setting:animation "setting:animation") | [Lighting](/dokuwiki/doku.php?id=setting:light "setting:light") | [Background](/dokuwiki/doku.php?id=setting:bg "setting:bg") | [Ray](/dokuwiki/doku.php?id=setting:ray "setting:ray") | [Transparency](/dokuwiki/doku.php?id=setting:transparency "setting:transparency")

## Special Cases

[Large Structures](/dokuwiki/doku.php?id=rendering:large "rendering:large")

rendering.txt · Last modified: 2013/08/19 21:00 (external edit)
