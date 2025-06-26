# selection

selection

### Table of Contents

  * Atom Selections

    * Selection by Matched Identifiers

    * Selection by Properties

    * Logical and Similarity Operations

    * Completing Groups

    * Geometric and Chemical Proximity Operations

    * Chemical Category

    * Sequence and Position




# Atom Selections

Conceptually, PyMOL atom selections are simply lists of atoms. These lists are used as input to many commands, as [selection expression](/dokuwiki/doku.php?id=concept:selection_expression "concept:selection_expression") arguments 

If you wish to use a selection expression more than once, then you can use the [select](/dokuwiki/doku.php?id=command:select "command:select") command to associate it with a [unique name](/dokuwiki/doku.php?id=concept:unique_name "concept:unique_name") that can be used repeatedly in other selection expressions. 

An [alphabetical list](/dokuwiki/doku.php?id=selection:alpha "selection:alpha") of selection operators is available along with a list of those that are [deprecated](/dokuwiki/doku.php?id=selection:deprecated "selection:deprecated"). 

Selection [macros](/dokuwiki/doku.php?id=selection:macros "selection:macros") provide a quickhand way of specifying atoms without having to provide explicit operators. 

## Selection by Matched Identifiers

[name](/dokuwiki/doku.php?id=selection:name "selection:name") | [resname](/dokuwiki/doku.php?id=selection:resname "selection:resname") | [resident](/dokuwiki/doku.php?id=selection:resident "selection:resident") | [chain](/dokuwiki/doku.php?id=selection:chain "selection:chain") | [segment](/dokuwiki/doku.php?id=selection:segment "selection:segment") | [object](/dokuwiki/doku.php?id=selection:object "selection:object") | [altloc](/dokuwiki/doku.php?id=selection:altloc "selection:altloc") | [id](/dokuwiki/doku.php?id=selection:id "selection:id") | [index](/dokuwiki/doku.php?id=selection:index "selection:index") | [rank](/dokuwiki/doku.php?id=selection:rank "selection:rank") | [element](/dokuwiki/doku.php?id=selection:element "selection:element")

## Selection by Properties

[b](/dokuwiki/doku.php?id=selection:b "selection:b") | [q](/dokuwiki/doku.php?id=selection:q "selection:q") | [partial_charge](/dokuwiki/doku.php?id=selection:partial_charge "selection:partial_charge") | [formal_charge](/dokuwiki/doku.php?id=selection:formal_charge "selection:formal_charge") | [ss](/dokuwiki/doku.php?id=selection:ss "selection:ss") | [hetatm](/dokuwiki/doku.php?id=selection:hetatm "selection:hetatm") | [bonded](/dokuwiki/doku.php?id=selection:bonded "selection:bonded") | [enabled](/dokuwiki/doku.php?id=selection:enabled "selection:enabled") | [present](/dokuwiki/doku.php?id=selection:present "selection:present") | [state](/dokuwiki/doku.php?id=selection:state "selection:state") | [visible](/dokuwiki/doku.php?id=selection:visible "selection:visible") | [color](/dokuwiki/doku.php?id=selection:color "selection:color") | [cartoon_color](/dokuwiki/doku.php?id=selection:cartoon_color "selection:cartoon_color") | [ribbon_color](/dokuwiki/doku.php?id=selection:ribbon_color "selection:ribbon_color") | [flag](/dokuwiki/doku.php?id=selection:flag "selection:flag") | [numeric_type](/dokuwiki/doku.php?id=selection:numeric_type "selection:numeric_type") | [text_type](/dokuwiki/doku.php?id=selection:text_type "selection:text_type") | [User-Defined Properties (New in 1.6.1)](/dokuwiki/doku.php?id=selection:user-defined_properties "selection:user-defined_properties")

## Logical and Similarity Operations

[not](/dokuwiki/doku.php?id=selection:not "selection:not") | [and](/dokuwiki/doku.php?id=selection:and "selection:and") | [or](/dokuwiki/doku.php?id=selection:or "selection:or") | [in](/dokuwiki/doku.php?id=selection:in "selection:in") | [like](/dokuwiki/doku.php?id=selection:like "selection:like")

## Completing Groups

[byresidue](/dokuwiki/doku.php?id=selection:byresidue "selection:byresidue") | [bychain](/dokuwiki/doku.php?id=selection:bychain "selection:bychain") | [bysegment](/dokuwiki/doku.php?id=selection:bysegment "selection:bysegment") | [byobject](/dokuwiki/doku.php?id=selection:byobject "selection:byobject") | [bycalpha](/dokuwiki/doku.php?id=selection:bycalpha "selection:bycalpha") | [byfragment](/dokuwiki/doku.php?id=selection:byfragment "selection:byfragment") | [bymolecule](/dokuwiki/doku.php?id=selection:bymolecule "selection:bymolecule")

## Geometric and Chemical Proximity Operations

[within](/dokuwiki/doku.php?id=selection:within "selection:within") | [near_to](/dokuwiki/doku.php?id=selection:near_to "selection:near_to") | [beyond](/dokuwiki/doku.php?id=selection:beyond "selection:beyond") | [around](/dokuwiki/doku.php?id=selection:around "selection:around") | [expand](/dokuwiki/doku.php?id=selection:expand "selection:expand") | [neighbor](/dokuwiki/doku.php?id=selection:neighbor "selection:neighbor") | [bound_to](/dokuwiki/doku.php?id=selection:bound_to "selection:bound_to") | [extend](/dokuwiki/doku.php?id=selection:extend "selection:extend") | [gap](/dokuwiki/doku.php?id=selection:gap "selection:gap")

## Chemical Category

[polymer](/dokuwiki/doku.php?id=selection:polymer "selection:polymer") | [solvent](/dokuwiki/doku.php?id=selection:solvent "selection:solvent") | [organic](/dokuwiki/doku.php?id=selection:organic "selection:organic") | [inorganic](/dokuwiki/doku.php?id=selection:inorganic "selection:inorganic") | [hydrogens](/dokuwiki/doku.php?id=selection:hydrogens "selection:hydrogens") | [donors](/dokuwiki/doku.php?id=selection:donors "selection:donors") | [acceptors](/dokuwiki/doku.php?id=selection:acceptors "selection:acceptors")

## Sequence and Position

[pepseq](/dokuwiki/doku.php?id=selection "selection") | [first](/dokuwiki/doku.php?id=selection:first "selection:first") | [last](/dokuwiki/doku.php?id=selection:last "selection:last")

selection.txt · Last modified: 2013/09/17 18:54 by bell
