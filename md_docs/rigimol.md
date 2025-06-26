# rigimol

rigimol

### Table of Contents

  * Molecular Morphing and RigiMOL

  * RigiMOL

    * The Basics

    * The Morph Command

    * Hints

  * Morphing Large Structures

  * Creating a Linear Morph

    *       * More on Smooth

  * Refining a Morph

    *       * Troubleshooting




# Molecular Morphing and RigiMOL

RigiMOL is a molecular morphing program to create aesthetically pleasing trajectories from a starting conformation to and ending conformation. 

Since PyMOL 1.6, standard morphing operations are available directly from the PyMOL user interface. 

  * Internal GUI: A > generate > morph

  * command: [morph](/dokuwiki/doku.php?id=command:morph "command:morph")




In previous version (<1.6), RigiMOL was loosely tied into PyMOL, so the typical scenario required a little programming, but we have prepackaged scripts (see the Morphing section of the Moviemaking Tutorial) to help you achieve your goals. 

# RigiMOL

## The Basics

The following code creates a morph similar to [this animated gif.](http://pymol.org/dsc/dokuwiki/lib/exe/fetch.php?media=morph.gif "http://pymol.org/dsc/dokuwiki/lib/exe/fetch.php?media=morph.gif")

Because RigiMOL creates a trajectory from a start state and an end state, you need just that – a two state object for RigiMOL to morph. Also, RigiMOL works with absolute coordinates, so it's advisable to align your structures before morphing them. We also suggest removing non-1-to-1 mapped atoms: if a solvent atom appears in the start state, but not the end state, that will confuse RigiMOL (how does it morph to nothing?). Knowing this, let's look at typical RigiMOL usage: 
    
    
    # fetch two conformers
    fetch 2iww 2iwv, async=0
    #  
    # remove non aligned structure
    remove 2iww and not c. A
    remove 2iwv and not c. A
    remove solvent or het
    #
    # align the two structures
    align 2iww, 2iwv, object=aln, cycles=0
    #   
    # create the input structure from the aligned structures
    # we begin with the start state
    create rin, 2iww in aln, 1, 1
    #   
    # now add the end-state as state 2
    create rin, 2iwv in aln, 1, 2
    #   
    # now we import the epymol module, this is where
    # the rigimol code can be found
    from epymol import rigimol
    #   
    # call RigiMOL.  Here we ask it to take our 2-state
    # object "rin", create the morph object "rout",
    # use low refinement (2), and update the PyMOL GUI
    # while it thinks (async_=1)
    rigimol.morph("rin", "rout", refinement=2)

This step will take a few minutes. PyMOL will update your screen as it calculates the morph. The morph will be a multi-state object. To save the morph type: 
    
    
    save rigimol_morph.pdb, rout, state=0

Now, if you want to see the opening-closing of this protein, you can type something like: 
    
    
    # position the protein
    orient
    turn y, -90
    #
    # zoom on the opening/closing
    zoom rout and i. 222, 10
    #
    # show as nicely rendered spheres
    set sphere_mode, 5
    as spheres
    mplay

## The Morph Command

To access the actual morph command you can type: 
    
    
    from epymol import rigimol

and then make calls to “rigimol.morph”. 

The usage for the morph command is: 
    
    
    morph(source, target, first=0, last=-1, refinement=10, async_=-1, steps=30, quiet=1)

  * source = name of input object

  * target = name of rigimol output

  * first = which is first state?

  * last = which is last state?

  * refinement = number of steps of refinement PyMOL does; high = slower but more appealing

  * async_ = if async_=0 then RigiMOL will do the work in the background and update PyMOL when it's done; if async_=1 then RigiMOL will update PyMOL on the fly as it calculate the morph.

  * steps = number of states generated in the transition from conformation 1 to conformation 2. (Coming soon in version 1.4.2.)




## Hints

  * PyMOL has a bug where it can sometimes destroy memory after the multi-state morph is returned. Save your data immediately after RigiMOL returns.




# Morphing Large Structures

To morph large structures, first make a linear morph and then use RigiMOL to refine the morph. Both steps are discussed below. 

# Creating a Linear Morph

PyMOL can quickly create an unrefined linear morph. This just creates the linear combination of coordinates of the N-states from the start state to the end states. That is, for the current state S and total number of states N: (x,y,z) = (S/N) * (x,y,z) + ((N-S) / N) * (x,y,z) 

For example, let's create a linear morph of PDB 2kyv from it's first state to it's fourth state: 
    
    
    # fetch the 20-state PDB, we will just use states 1 and 4.
    
    fetch 2kyx, async=0
    
    # create the input object with 100 states.  The first 50
    # will be state 1 from 2kyx, and 51-100 will be 
    # state 4 from 2kyx.  The new object's name will be "m_out"
    # for "morph_out"
    
    for x in range(51): cmd.create("m_out", "2kyx", 1, x)
    for x in range(51,101): cmd.create("m_out", "2kyx", 4, x)
    
    # perform the linear smoothing across all states
    
    smooth m_out, 1, 100, 1, 100

### More on Smooth

The 'smooth' command takes the following parameters: 

  * object/selection name

  * passes

  * window size

  * first state

  * last state




PyMOL will, for 'passes' number of times, smooth all atoms in the selection from state 'first' to state 'last' using a window size of 'window size'. PyMOL centers the smoothing window on the current state and smooths +/- 'window size' / 2 in each direction. Thus, if we have 100 states, and specify a window size of 100, we get the following: 

  * state1 = state1 + state2 + … + state50 / 50; 

  * state2 = (state1) + state2 + … + state50 / 51;

  * state3 = (state1 + state2) + state3 + … + state50 / 52;

  * …

  * state50 = (state1 + state2 + … state49) + state50 + state51 + … + state100 / 100




where all “backward” states are in parenthesis. 

# Refining a Morph

PyMOL can now refine the linear morph using RigiMOL. Continuing from the above example using object “m_in” we now call rigimol.refine: 
    
    
    # import rigimol so can we can use it:
      
    from epymol import rigimol
    
    # use rigimol to refine the linear morph
    
    rigimol.refine(2, "m_out")

This “refine” command takes the linear morph and uses PyMOL's internal sculpting and bump checking to improve the morph. Here, 25 is the number of iterations to apply when adjusting the parameters. For a massive structure, like the ribosome, try smaller numbers like those < 10\. 

### Troubleshooting

For best results: 

  * remove any segment identifiers: alter *, segi=''

  * ensure chain identifiers match between structures

  * ensure residue numbers are the same between structures

  * remove solvent

  * remove heteroatoms




rigimol.txt · Last modified: 2020/04/21 06:10 by holder
  *[GUI]: Graphical User Interface
