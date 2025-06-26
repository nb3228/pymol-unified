# labels

labels

### Table of Contents

  * Labels

    * See Also

    * Fonts

    * Special Characters




# Labels

Currently, per-atom labels are the easiest the create: 
    
    
    label elem n, "Nitrogen"
    
    
    label name ca, resn+resi+chain

## See Also

[label](/dokuwiki/doku.php?id=command:label "command:label")

## Fonts

The fonts included with PyMOL 1.x are listed here with their identifiers. 

  * 5 - DejaVuSans (regular)

  * 6 - DejaVuSans (italic)

  * 7 - DejaVuSans (bold)

  * 8 - DejaVuSans (bold italic)

  * 9 - DejaVuSerif (regular)

  * 10 - DejaVuSerif (bold)

  * 11 - DejaVuMono (regular)

  * 12 - DejaVuMono (italic)

  * 13 - DejaVuMono (bold)

  * 14 - DejaVuMono (bold italic)

  * 15 - Gentium (regular)

  * 16 - Gentium (italic) 

  * 17 - DejaVuSerif (italic)

  * 18 - DejaVuSerif (bold italic)




## Special Characters

PyMOL uses the UTF8 unicode encodings for special characters. 

As of PyMOL 1.x, all scalable fonts (5-18) can display special characters such as Greek “alpha” and “beta” via unicode. However, the Angstrom symbol is missing from fonts 15 & 16\. 
    
    
    set label_font_id, 15
    
    # Alpha-Helix using greek Alpha character
    
    label 12/ca, "\316\261-Helix"
    
    # Beta-Sheet using greek Beta character
    
    label 20/ca, "\316\362-Sheet"
    

[![:greek.png](/dokuwiki/lib/exe/fetch.php?media=greek.png)](/dokuwiki/lib/exe/detail.php?id=labels&media=greek.png "greek.png")
    
    
    set label_font_id, 17
    
    label 6/ca, "5.0 \342\204\253"

[![:angstrom.png](/dokuwiki/lib/exe/fetch.php?media=angstrom.png)](/dokuwiki/lib/exe/detail.php?id=labels&media=angstrom.png "angstrom.png")

labels.txt · Last modified: 2013/08/19 21:01 (external edit)
