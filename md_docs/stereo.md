# stereo

stereo

### Table of Contents

  * NVidia Stereo 3D NVision

  * iZ3D for Open-Source PyMOL

  * iZ3D for Incentive PyMOL

  * Zalman




### NVidia Stereo 3D NVision

PyMOL does work with the NVidia 3D NVision setup, but needs the right hardware. For example, PyMOL does work with a FX 580. Please see the following links for more help on setting up the 3D NVision System on your machine with a Quadro card. 

Please note that the Quadro FX1400 does not support stereographic viewing in Win7. The response from NVidia is: most stereo features on Windows 7 are only support[ed] on Quadro *600 series(G80) and beyond. I believe clone mode stereo is support on Quadro *500 series(G70) and beyond. 

#### NVidia 3D Support Drivers

  * [List of approved 3D Display Devices](http://www.nvidia.com/object/3d-vision-pro-requirements.html "http://www.nvidia.com/object/3d-vision-pro-requirements.html")

  * [Latest Graphics Drivers Vista/Win7](http://www.nvidia.com/object/quadro-win7-winvista-64bit-259.12-whql-driver.html "http://www.nvidia.com/object/quadro-win7-winvista-64bit-259.12-whql-driver.html")

  * [Latest Graphics Drivers WinXP x64](http://www.nvidia.com/object/quadro-winxp-x64-259.12-whql-driver.html "http://www.nvidia.com/object/quadro-winxp-x64-259.12-whql-driver.html")

  * [Latest Quadro Drivers](http://www.nvidia.com/object/quadro-3d-vision-usbdriver-258.49-driver.html "http://www.nvidia.com/object/quadro-3d-vision-usbdriver-258.49-driver.html")




#### NVidia Setup Documentation

  * [NVIDIA 3D Vision on Quadro Professional Graphics Boards](http://pymol.org/sites/default/files/configuring nvidia 3d vision on quadro_0.pdf "http://pymol.org/sites/default/files/configuring nvidia 3d vision on quadro_0.pdf")




### iZ3D for Open-Source PyMOL

If you build PyMOL from the open-source code and want to use it with an iZ3D, use these instructions: 

  * For clone mode (Win: nView Clone, or Linux: TwinView Option Stereo 4), launch as: ./pymol -t 12 

  * For dual-screen mode, launch as: ./pymol -t 11 and then spread the window over both screens, edge to edge. 

  * Then, issue 'stereo on' to activate stereo.




### iZ3D for Incentive PyMOL

If you are using incentive PyMOL, use these instructions; 

  * Launch PyMOL in stereo mode for iz3D using: [pathname]\PymolWin.exe“ -J -x -S -t 11 -X 0 -Y 0 -W 3141 -H 990 

  * First, set the convergence based on the central object on the screen. Then, set the separation to optimize your 3D image with minimal ghosting 




### Zalman

Zalman monitors use interlaced technology to create the 3D effect. All odd numbered lines go to the right eye (for example) and all even numbered lines can go to the left eye. This creates two separate images for the viewer, thus creating the 3D effect. 

This solution is cheaper than the NVidia NVision systems as all you need is a Zalman monitor (no Quadro card or 3D NVision kit necessary). However, the Zalman 3D does have limitations: 

  * The resolution of your image is halved;

  * Non-OpenGL Menu text is nearly unreadable

  * The viewing angle is very small (about 10 degrees)




#### Supported Zalman Monitors

  * [AOC E2352PHZ](http://us.aoc.com/monitor_displays/e2352phz "http://us.aoc.com/monitor_displays/e2352phz")

  * [Asus VG23AH](https://www.asus.com/Monitors/VG23AH/ "https://www.asus.com/Monitors/VG23AH/")

  * [Asus VG27AH](https://www.asus.com/us/Monitors/VG27AH/ "https://www.asus.com/us/Monitors/VG27AH/")

  * [Asus VG278HR](https://www.asus.com/Monitors/VG278HR/ "https://www.asus.com/Monitors/VG278HR/")

  * [LG D2342P](http://www.lg.com/hk_en/monitors/lg-D2343P-home "http://www.lg.com/hk_en/monitors/lg-D2343P-home")

  * [LG D2343PB-BN](http://www.lg.com/us/commercial/lcd-computer-monitors/lg-D2343PB-BN "http://www.lg.com/us/commercial/lcd-computer-monitors/lg-D2343PB-BN")

  * [LG D2342P-PN](http://www.lg.com/us/computer-products/monitors/LG-led-monitor-D2342P-PN.jsp "http://www.lg.com/us/computer-products/monitors/LG-led-monitor-D2342P-PN.jsp")




stereo.txt · Last modified: 2016/01/22 21:54 by holder
