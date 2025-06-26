# Category: UI Script

## Jump

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Type  | [Python Script](/index.php/Running_Scripts "Running Scripts")  
---|---  
Download  | [jump.py](https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/jump.py)  
Author(s)  | [Sean M. Law](/index.php/User:Slaw "User:Slaw")  
License  | \-   
This code has been put under version control in the project [Pymol-script-repo](/index.php/Git_intro "Git intro")  
  
## Contents

  * 1 Description
  * 2 Usage
  * 3 Example
  * 4 PyMOL API
  * 5 See Also



### Description

In the case of movies or simulation trajectories, it is often useful to quickly jump between different frames but not necessarily in any order. This script allows the user to specify a list of frames and flipping through the set of frames is controlled by a user-defined keystroke (limited to those available through the [Set Key](/index.php/Set_Key "Set Key") command). This script is a nice tool for identifying different protein states from a simulation and can be used in conjunction with the [Modevectors](/index.php/Modevectors "Modevectors") script in a very useful way by first characterizing the directionality of the movements (using [Modevectors](/index.php/Modevectors "Modevectors")) and then visualizing the change directly using this script. 

Note: This script also requires the [Check Key](/index.php/Check_Key "Check Key") script! 

  


### Usage

load the script using the [run](/index.php/Run "Run") command 
    
    
    jump [frame1 [,...frameN [,sleep=0.1 [,x=1 [,key='pgdn']]]]]
    

### Example
    
    
    jump 1, 10, 100
    
    #This causes the movie to jump from frame 1 to frame 10 and finally to frame 100 each time the 'pgdn' key is pressed
    #Note that there is no limit to the number of frames that can be specified
    
    jump 1, 10, 100, 1000, key='F1'
    #This sets the hotkey to F1.  When pressed, it will jump through the frames 1, 10, 100, 1000 respectively
    
    jump 1, 10, 100, sleep=0.01
    #This does the same thing as the first example except that the pause between frames is 10x faster
    #One could consider the "sleep" option to be similar to the frame rate
    
    jump 1, 10, 100, x=10
    #This does the same thing as the first example except that each press of the hotkey will 
    #cause the script to run through the specified frames 10 times before stopping.
    #Note that you will not be able to control or exit the script during the 10 iterations!
    

### PyMOL API
    
    
    from pymol import cmd
    import time
    
    def jump (*args, **kwargs):
        
        """
        Author Sean M. Law
        University of Michigan
        seanlaw_(at)_umich_dot_edu
        """
    
        keystroke='pgdn'
        
        for key in kwargs:
            if (key == "key"):
                keystroke=kwargs["key"]
                keystroke=check_key(keystroke)
                if (keystroke == None):
                    return
            else:
                continue
    
        cmd.set_key(keystroke, jumpit, args, kwargs)
    
        return
    
    cmd.extend("jump", jump)
    
    def jumpit (*args, **kwargs):
    
        x=1
        sleep=0.1
    
        for key in kwargs:
            if (key == "sleep"):
                sleep=kwargs["sleep"]
                sleep=float(sleep)
            elif (key == "x"):
                x=kwargs["x"]
                x=int(x)
            else:
                continue
    
        args=list(args)
        
        for i in range (x):
            for j in range (len(args)):
                cmd.frame(int(args[j]))
                cmd.refresh()
                time.sleep(sleep)
    
        return
    cmd.extend("jumpit", jumpit)
    

  


### See Also

[Check Key](/index.php/Check_Key "Check Key")

Retrieved from "[https://pymolwiki.org/index.php?title=Jump&oldid=11148](https://pymolwiki.org/index.php?title=Jump&oldid=11148)"


---

