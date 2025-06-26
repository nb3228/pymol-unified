# properties

properties

### Table of Contents

  * Example

  * Object Level Properties

  * Atom Level Properties

  * Selection Language




Beginning with incentive PyMOL 1.6, arbitrary object level and atom level properties are supported. Now, for each object and each atom you can possibly have an arbitrary dictionary associated with each. 

## Example
    
    
    cmd.set("load_object_props_default", "*")
    cmd.load("http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=2519&disopt=3DSaveSDF", "caffeine", format="sdf")
    cmd.wizard("message", "Vol: %.1f" % cmd.get_property("PUBCHEM_SHAPE_VOLUME", "caffeine"))

## Object Level Properties

Loading from files: 
    
    
    load abc.sdf, object_props=*
    
    
    set load_object_props_default, *
    load abc.sdf

Querying from the command line: 
    
    
    get_property_list abc
    get_property some_name, abc

## Atom Level Properties

Loading from files: 
    
    
    load abc.mae, atom_props=*
    
    
    set load_atom_props_default, *
    load abc.mae

Querying from the command line: 
    
    
    iterate abc, print(properties["some_name"])  # dictionary notation
    iterate abc, print(properties.some_name)     # dot notation
    iterate abc, print(p["some_name"])           # shortcut
    iterate abc, print(p.some_name)
    iterate abc, print(properties.all)           # dictionary of all properties

Coloring by numeric property value: 
    
    
    spectrum properties.some_name, blue yellow red

Scaling a setting by an atom property: 
    
    
    stored.values = []
    iterate all, stored.values.append(float(p["some_value"]))
    stored.v_min = min(stored.values)
    stored.v_range = max(stored.values) - stored.v_min
    alter all, s["transparency"] = (float(p["some_value"]) - stored.v_min) / stored.v_range

## Selection Language

Numeric: 
    
    
    select mysele, p.some_value = 3
    select mysele, p.some_value > 2.4

String, supports matching against sets with “+”: 
    
    
    select myseleA,    p.some_value in A
    select myseleAorB, p.some_value in A+B

properties.txt · Last modified: 2020/05/19 10:57 by holder
