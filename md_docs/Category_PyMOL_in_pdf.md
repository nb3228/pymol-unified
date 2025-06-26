# Category: PyMOL in pdf

## 3d pdf

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
    * 1.1 Requirements
    * 1.2 Get IDFT file from PyMOL
  * 2 Windows - Convert the IDTF to U3D
    * 2.1 Quick
    * 2.2 General
  * 3 Linux install - Convert the IDTF to U3D
  * 4 LaTeX code
    * 4.1 LaTeX code pdf
    * 4.2 LaTeX beamer presentation
  * 5 More on 3D PDFs
    * 5.1 See Also



# Overview

PyMOL can convert to formats (vrml2 and idtf) that can be converted to a [3D PDF](/index.php/File:Pymol.pdf "File:Pymol.pdf") (will not work with most PDF browser plugins; must be downloaded and viewed with certain viewers liked Adobe Acrobat 9.2+). 

**Note**. The figure is only represented in cartoon format. No sticks yet... 

## Requirements

  * PyMOL
  * [Universal 3D Sample Software - u3d converter - IDTF to U3D](http://sourceforge.net/projects/u3d/files/Universal%203D%20Sample%20Software/)
  * LaTeX (pdflatex)



## Get IDFT file from PyMOL

  * Save your molecule to an IDTF file in PyMOL:


    
    
    save pymol.idtf, *
    

PyMOL will print a line that looks like: 
    
    
     3Daac=20.0, 3Droll=0, 3Dc2c=0 0 1, 3Droo=62.45, 3Dcoo=0 0 -62.45
    

copy this line into the pymol.tex file overwriting the same line in the file. 

# Windows - Convert the IDTF to U3D

## Quick

  1. Download **Universal 3D Sample Software - u3d converter - IDTF to U3D** , and extract to Desktop  

  2. Navigate to the **bin** folder: \Desktop\U3D_A_061228_5\Bin\Win32\Release   

  3. Copy **pymol.idtf** in here.
  4. Hold **shift** key, right click in folder, click **Open command window here**.



Then copy this inot command window: 
    
    
    IDTFConverter -input pymol.idtf -output pymol.u3d
    

Copy the pymol.u3d into your LaTeX folder 

## General

  1. Download **Universal 3D Sample Software - u3d converter - IDTF to U3D** , and extract to C:\Program Files (x86)\U3D_A_061228_5  

  2. Make a **IDTF2U3D.bat** and afterwards copy to **C:\Windows**.


    
    
    @echo off
    rem -- Run IDTFConverter --
    set IDTF_EXE_DIR=C:\Program Files (x86)\U3D_A_061228_5\Bin\Win32\Release
    
    if exist "%IDTF_EXE_DIR%\IDTFConverter.exe" goto haveProg
    echo "%IDTF_EXE_DIR%\IDTFConverter.exe" not found
    goto eof
    
    :haveProg
    if .%1==. goto noarg
    
    "%IDTF_EXE_DIR%\IDTFConverter.exe" -input %1 -output %~n1.u3d 
    goto eof
    
    :noarg
    echo no arguments given.
    echo IDTF2U3D pymol.idtf
    
    :eof
    echo Done
    

  1. Navigate to **pymol.idtf**. Right click, Open with..., Browse to: **C:\Windows** and select **IDTF2U3D.bat**



Copy the pymol.u3d into your LaTeX folder 

# Linux install - Convert the IDTF to U3D

  * Currently you have to compile the u3d converter on Linux. I did that with:
  * Some versions of Acrobat on Linux incorrectly parse the 3D data. Adobe knows about this and plans to fix it. Ironically, I created a 3D PDF on Linux but could only view it on Mac OS X.



See current versions here: <http://www2.iaas.msu.ru/tmp/u3d>
    
    
    cd ~/Downloads/ ;
    wget http://www2.iaas.msu.ru/tmp/u3d/u3d-1.4.3.tar.gz ;
    tar -xzf u3d-1.4.3.tar.gz ;
    cd u3d-1.4.3 ;
    cmake . ;
    make ;
    chmod g+wx . ;
    set IDDIR=$PWD ;
    
    cat > IDTF2U3D << EOF
    #!/bin/tcsh
    set FOLD=\$PWD
    set IN=\$PWD/\$1
    cp -f \$IN $IDDIR
    cd $IDDIR
    ./IDTFConverter -input \$1 -output \$1:r:t.u3d
    cp -f \$1:r:t.u3d \$FOLD/
    EOF
    
    chmod +x IDTF2U3D
    
    ln -s $PWD/IDTF2U3D ~/bin/IDTF2U3D
    

Then write in terminal 
    
    
    IDTF2U3D pymol.idtf
    

# LaTeX code

## LaTeX code pdf

Example [3D PDF](/index.php/File:Pymol.pdf "File:Pymol.pdf")  
Download Example, and open in acrobat reader: [1zqa PDF](http://www.fys.ku.dk/~tlinnet/1zqa.pdf)

  * The following LaTeX code saved as "pymol.tex":


    
    
    \documentclass[a4paper]{article}
    \usepackage[3D]{movie15}
    \usepackage[UKenglish]{babel}
    \usepackage[colorlinks=true]{hyperref} 
    \begin{document}
    \title{PyMOL 3D Objects in PDF}
    \author{Jason Vertrees}
    \maketitle
    \begin{figure}[!htb]
    \centering
    \includemovie[
    poster,
    toolbar, %same as `controls'
    label=pymol.ud3
    text=(pymol.u3d),
    3Dlights=CAD,
    % replace the next line with what PyMOL output
    3Daac=20.0, 3Droll=0, 3Dc2c=0 0 1, 3Droo=243.39, 3Dcoo=0 0 -243.39
    ]{\linewidth}{\linewidth}{pymol.u3d}
    \caption{A PyMOL object embedded in PDF, using U3D data format.}
    \label{pym:ex3d}
    \end{figure}
    \end{document}
    

## LaTeX beamer presentation

Download Example, and open in acrobat reader: [beamer PDF](http://www.fys.ku.dk/~tlinnet/beamer.pdf)   
Download Example, and open in acrobat reader: [beamer 1zqa PDF](http://www.fys.ku.dk/~tlinnet/beamer_1zqa.pdf)
    
    
    \documentclass{beamer}
    \usepackage[3D]{movie15}
    \usepackage[UKenglish]{babel}
    
    %% See: http://en.wikibooks.org/wiki/LaTeX/Presentations
    %% See: http://www.hartwork.org/beamer-theme-matrix/
    \usetheme{Copenhagen} 
    \usecolortheme{beaver}
    
    \begin{document}
    \title{Simple Beamer Class}   
    \author{Troels Linnet} 
    \date{\today} 
    \frame{\titlepage} 
    \frame{\frametitle{Table of contents}\tableofcontents} 
    
    \section{3D object} 
    \frame{\frametitle{Pymol 3D object}
    \begin{figure}[!htb]
    \centering
    \includemovie[
    poster,
    toolbar, %same as `controls'
    label=pymol.ud3
    text=(pymol.u3d),
    3Dlights=CAD,
    % replace the next line with what PyMOL output
    3Daac=20.0, 3Droll=0, 3Dc2c=0 0 1, 3Droo=243.39, 3Dcoo=0 0 -243.39
    ]{10cm}{6cm}{pymol.u3d}
    \caption{A PyMOL object embedded in PDF, using U3D data format.}
    \label{pym:ex3d}
    \end{figure}
    }
    
    \section{Lists}
    \subsection{Lists I}
    \frame{\frametitle{Unnumbered lists}
    \begin{block}{My bloc}
    test of bloc
    \end{block}
    
    \begin{itemize}
    \item Introduction to  \LaTeX  
    \item Beamer class
    \end{itemize} 
    }
    
    \end{document}
    

  * Create the PDF using LaTeX:


    
    
     pdflatex pymol.tex
    

# More on 3D PDFs

  * [3D PDFs at Adobe](http://www.adobe.com/manufacturing/3dpdfsamples/3dsolutions/)
  * [Future of scientific publishing](http://molecularmodelingbasics.blogspot.dk/2009/12/future-of-scientific-publishing-is-here.html)
  * [Instructions for Jmol](http://jpswalsh.com/how-to-embed-a-3d-molecule-in-a-pdf/)



## See Also

  * [ Movie in pdf](/index.php/Movie_pdf "Movie pdf")



Retrieved from "[https://pymolwiki.org/index.php?title=3d_pdf&oldid=11982](https://pymolwiki.org/index.php?title=3d_pdf&oldid=11982)"


---

## Movie pdf

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Contents

  * 1 Overview
    * 1.1 Requirements
    * 1.2 Get .mpg file from PyMOL
    * 1.3 LaTeX code
    * 1.4 LaTeX beamer presentation
    * 1.5 More on movies in PDFs
    * 1.6 See Also



# Overview

PyMOL can save videos to formats (.mpg) that can be converted to a [movie PDF](/index.php?title=File:Pymol_movie&action=edit&redlink=1 "File:Pymol movie \(page does not exist\)") (will not work with most PDF browser plugins; must be downloaded and viewed with certain viewers liked Adobe Acrobat 9.2+). 

## Requirements

  * PyMOL
  * LaTeX (pdflatex)



## Get .mpg file from PyMOL

  * Save a **.mpg** file from your video.   




You can try this tutorial:   
[ Example movie](/index.php/Biochemistry_student_intro#GUI_-_Scene_loop "Biochemistry student intro")

## LaTeX code

Download Example, and open in acrobat reader: [movie PDF](http://www.fys.ku.dk/~tlinnet/movie.pdf)   
Note, the resolution of the movie, has been cut intentionally to save filesize.  
Download Example, and open in acrobat reader: [movie PDF High Resolution](http://www.fys.ku.dk/~tlinnet/movie_highres.pdf)   


  * The following LaTeX code saved as "pymol.tex":


    
    
    \documentclass[a4paper]{article}
    \usepackage[3D]{movie15}
    \usepackage[UKenglish]{babel}
    \usepackage[colorlinks=true]{hyperref} 
    \begin{document}
    \title{PyMOL 3D Objects in PDF}
    \author{Troels Linnet}
    \maketitle
    
    \begin{figure}[!htb]
    \centering
    \includemovie[
    poster,
    toolbar, %same as `controls'
    text={\small(Click to play movie)}
    ]{320px}{216px}{movie.mpg}
    \caption{A PyMOL movie object embedded in PDF, using mpg format.}
    \label{mov:ex3d}
    \end{figure}
    
    See movie \ref{mov:ex3d}
    
    \end{document}
    

## LaTeX beamer presentation

Download Example, and open in acrobat reader: [beamer movie PDF](http://www.fys.ku.dk/~tlinnet/beamer_movie.pdf)   
Note, the resolution of the movie, has been cut intentionally to save filesize.   
Download Example, and open in acrobat reader: [beamer movie PDF High Resolution](http://www.fys.ku.dk/~tlinnet/beamer_movie_highres.pdf)   

    
    
    \documentclass{beamer}
    \usepackage[3D]{movie15}
    \usepackage[UKenglish]{babel}
    
    %% See: http://en.wikibooks.org/wiki/LaTeX/Presentations
    %% See: http://www.hartwork.org/beamer-theme-matrix/
    \usetheme{Copenhagen} 
    \usecolortheme{beaver}
    
    \begin{document}
    \title{Simple Beamer Class}   
    \author{Troels Linnet} 
    \date{\today} 
    \frame{\titlepage} 
    \frame{\frametitle{Table of contents}\tableofcontents} 
    
    \section{Movie} 
    \frame{\frametitle{Pymol Movie object}
    
    \begin{figure}[!htb]
    \centering
    \includemovie[
    poster,
    toolbar, %same as `controls'
    text={\small(Click to play movie)}
    ]{10cm}{6cm}{movie.mpg}
    \caption{A PyMOL movie object embedded in PDF, using mpg format.}
    \label{mov:ex3d}
    \end{figure}
    
    }
    
    \section{Lists}
    \subsection{Lists I}
    \frame{\frametitle{Unnumbered lists}
    \begin{block}{My bloc}
    test of bloc
    \end{block}
    
    \begin{itemize}
    \item Introduction to  \LaTeX  
    \item Beamer class
    \end{itemize} 
    }
    
    \end{document}
    

  * Create the PDF using LaTeX:


    
    
     pdflatex pymol.tex
    

## More on movies in PDFs

  * [PDFmovie](http://pages.uoregon.edu/noeckel/PDFmovie.html)



## See Also

  * [ 3d protein image in pdf](/index.php/3d_pdf "3d pdf")
  * [ Example Scene_loop](/index.php/Biochemistry_student_intro#GUI_-_Scene_loop "Biochemistry student intro")



Retrieved from "[https://pymolwiki.org/index.php?title=Movie_pdf&oldid=11101](https://pymolwiki.org/index.php?title=Movie_pdf&oldid=11101)"


---

