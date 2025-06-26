# Category: Stereo 3D Display

## Stereo 3D Display Options

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

This page is for aggregating the latest know-how and links to current Stereo 3D display options best suited for molecular graphics applications like PyMOL. Please strive to provide objective factual information based on first-hand experiences while using the displays for real work and teaching. 

Let's figure out together what stereo solutions work well in this brave new post-CRT world! 

## Contents

  * 1 Categories
  * 2 Active Stereo 3D (High-Refresh) Displays
    * 2.1 LCD Displays (120 Hz)
      * 2.1.1 NVidia NVision 3D Setup
        * 2.1.1.1 Necessary Hardware
        * 2.1.1.2 Necessary Software
        * 2.1.1.3 Installation Instructions
          * 2.1.1.3.1 System Setup
          * 2.1.1.3.2 Running PyMOL
    * 2.2 DLP Projection Televisions
  * 3 Passive Stereo 3D Displays
    * 3.1 One Piece Multi-layer LCD Displays
      * 3.1.1 Zalman
      * 3.1.2 LG
        * 3.1.2.1 On Linux
        * 3.1.2.2 On Windows
        * 3.1.2.3 On Mac
      * 3.1.3 iZ3D
    * 3.2 Mirror-based Multi-LCD Solutions
  * 4 Autostereoscopic LCD Displays
  * 5 Stereo 3D Projectors
    * 5.1 Active Stereo 3D DLP Projectors
    * 5.2 Passive Stereo 3D Adaptor Hardware for Active Stereo 3D Projectors
  * 6 Who Says What?



## Categories

  * **Active Stereo 3D** \-- requires expensive and/or bulky shutter glasses. For decades, this has been the standard solution for stereo 3D molecular visualization on the desktop.
  * **Passive Stereo 3D** \-- requires inexpensive lightweight polarized glasses. This is the standard solution for delivering stereo 3D to audiences of more than a small group of people.
  * **Autostereoscopic 3D** \-- means that no glasses are required. However, current autostereoscopic displays tend do not to work well for molecular graphics due to their inability to faithfully represent fine detail such as thin bonds and wire meshes.



## Active Stereo 3D (High-Refresh) Displays

This solution would be analogous to using desktop CRT monitors with shutter glasses. 

### LCD Displays (120 Hz)

  * [ACER GD235HZ](http://us.acer.com/ac/en/US/content/model/ET.UG5HP.001). 1920x1080 120Hz 2ms .


  * [Planar SA2311W](http://www.planar3d.com/3d-products/sa2311w/). This is a high-end 23" 3D-ready monitor. The resolution is 1900x1280 and has a 2ms refresh time. This worked well for us.


  * [ViewSonic VX2265wm](http://www.engadget.com/2008/08/26/viewsonic-shows-off-a-120hz-lcd-display-for-computers) (VX2268wm in Europe) - [On the market](http://www.google.com/products?q=ViewSonic%20VX2265wm&btnG=Search+Froogle&lmode=unknown) and [verified working under FC 12](http://sourceforge.net/mailarchive/forum.php?thread_name=DED5C399-7534-4D3D-8B19-E3676C4F1867%40weizmann.ac.il&forum_name=pymol-users)


  * [Samsung 2233RZ](http://www.nvidia.com/object/product_GeForce_3D_VisionBundle_us.html) \- On the market. Quad buffered stereo in Linux works with a [G8x based graphics core](http://en.wikipedia.org/wiki/Nvidia_Quadro) or better Quadro FX card with the 3 pin mini din stereo connector (currently, the cheapest card that works in Linux is the Quadro FX 3700), 195.22 (or newer) nvidia linux binary driver, and the Nvidia 3d vision kit. Even though the Quadro FX 1400/3450/4000 cards have a 3 pin stereo connector, these will not work with Nvidia 3D vision since these have core versions less than G8x. For more information see this forum post [at the Nvidia Forums](http://forums.nvidia.com/index.php?showtopic=91072&view=findpost&p=968627). - SP


  * USB only based stereo with the 3D vision kit works only in MS Windows (e.g. with a low end Quadro FX 370 that has no 3 pin mini din stereo connector). For more information see this forum post [at the Nvidia Forums](http://forums.nvidia.com/index.php?showtopic=91072&view=findpost&p=968627). - SP


  * The 195.22 Nvidia linux drivers do not support the Samsung 2233RZ in Stereo mode 3 or 10 for quad buffered stereo with other stereo kits, emitters, or goggles such as those purchased from NuVision, Stereographics, or Edimensional. You cannot use NuVision, Stereographics, or Edimensional goggles with the Nvidia 3D Vision emitters. - SP


  * NVidia 3D NVision kit only supports DirectX software for GeForce (gaming cards) on Windows; users are reporting that they are not able to run PyMOL with NVision with these cards. Get a newer model low end quadro (> G8x graphics core) without the 3 pin mini din (e.g. Quadro 370) or with the 3 pin mini din (e.g. Quadro 3700) for Windows.



#### NVidia NVision 3D Setup

The NVidia 3D NVision setup provides a very nice 3D experience. You need the following to enable PyMOL to show NVision 3D on Windows. Please review the hardware and software requirements before moving on to the installation and setup. 

##### Necessary Hardware

  * Monitor: 120 Hz LCD: a [2233RZ](http://www.samsung.com/us/consumer/office/monitors/specialty/LS22CMFKFV/ZA/index.idx?pagetype=prd_detail&returnurl=%7CSamsung) or a [Fuhzion vx2265wm](http://www.viewsonic.com/products/desktop-monitors/lcd/x-series/vx2265wm-fuhzion-lcd.htm%7CViewSonic)
  * Cable: [Dual Link DVI cable](http://images.google.com/imgres?imgurl=http://www.logicsupply.com/images/dvi_connector_types.gif&imgrefurl=http://www.logicsupply.com/faq&usg=__G2BLaVTqBN4ie8fz_LJR1zc3zBc=&h=261&w=440&sz=15&hl=en&start=0&sig2=_hFM6ICIsxPq5WIAv8BCqg&zoom=1&tbnid=NIcKIs_BW_2rmM:&tbnh=135&tbnw=228&ei=KHN2TL-UC8P_lgfr44nsCw&prev=/images%3Fq%3Ddual%2Blink%2Bdvi%26hl%3Den%26biw%3D1475%26bih%3D1042%26gbv%3D2%26tbs%3Disch:1&itbs=1&iact=hc&vpx=136&vpy=323&dur=3153&hovh=173&hovw=292&tx=227&ty=74&oei=KHN2TL-UC8P_lgfr44nsCw&esq=1&page=1&ndsp=30&ved=1t:429,r:6,s:0); most 120Hz monitors will come with this cable--regardless, the cable is necessary
  * Quadro Card: recent [Quadro](http://www.nvidia.com/page/quadrofx_family.html) series graphics card (not a GeForce card) such as an FX 380 or 570 or later. The GeForce cards do not support windowed openGL stereo, so we do not support these series of cards for the NVision 3D solution. For linux, you must have a quadro card that has a 3 pin mini din connector. The cheapest/oldest card that will work with linux is the Quadro 3700. 
    * **WARNING** : The Quadro FX1400 does not support 3d vision stereo on Windows7 or Linux.
  * Emitter: [3D Vision](http://www.nvidia.com/object/3d-vision-main.html%7CGeForce) hardware kit (an emitter with 3D shutter glasses). For Linux, make sure your kit comes with the 3 pin mini din "VESA" to 2.5mm stereo cable to connect from the stereo output on the video card into the emitter. See the [3 pin Mini-DIN connector](/index.php/3_pin_Mini-DIN_connector "3 pin Mini-DIN connector") article for tips on how to make one of these cables if yours is missing.
  * GeForce Cards from series 400 onward have gained OpenGl support in recent Nvidia driver iterations (314+). This allows Pymol to be viewed in 3D using the quad buffered stereo setting with a GeForce card, 120Hz screen and 3D Vision kit.



##### Necessary Software

  * Windows XP 32 bit (testing other OSs soon!), Windows Vista
  * Latest Quadro [Drivers from NVidia](http://www.nvidia.com/Download/index.aspx?lang=en-us%7CGraphics).
  * Latest [Graphics drivers for the NVision system](http://www.nvidia.com/Download/index.aspx?lang=en-us%7C3D)\--under **Product Type** choose **3D Vision**.



##### Installation Instructions

###### System Setup

  1. Install the Quadro **Graphics Drivers** and reboot your machines
  2. Install the NVision Installation, hooking up the 3D emitter and glasses as directed in the instructions 
     1. Make sure the 3D demos work
     2. Complete the **3D Vision Drivers** install (I had errors/warnings about old drivers but this did't matter)
  3. Specify how to drive the 3D by, click on



    

    

    **Windows Start Button** > **Control Panel** > **NVidia Control Panel** > **Manage 3D Settings** (tab) > **Global Settings** (tab on the right) > **Base Profile** (tab). Then, under **Settings** choose **Stereo - Display Mode**. Next, select **Generic Active Stereo (with NVidia IR Emitter)**. If you have a DLP monitor/TV choose the corresponding DLP option. You **must** also set **Stereo - Enable** to **on**.

###### Running PyMOL

That's it! PyMOL should now work in Quad Buffered 3D Stereo using the NVidia 3D NVision system. To run PyMOL in 3D mode on: 

  * Windows



    

    **Start > PyMOL > PyMOL > PyMOL 3D Launch (last menu option) > PyMOL Stereo (Quad Buffered 3D)**

  * Linux



    

    pymol -S -t 1
    _Note that hardware stereo may not work in Xorg unless window compositing is turned off; Gnome3, Unity, etc all use window compositing as part of their eye candy. A window manager that should work by default is the MATE desktop. Install this if you have trouble_

  * Mac



    

    Sorry, at this time the NVision system is not known to work on Macs.

### DLP Projection Televisions

Projection televisions tend to be too large and fuzzy for desktop use. Also, a band of about 20 pixels around on the edge of the display are invisible, and this limitation cannot be eliminated through overscan since the image must be scanned at native resolution in order to support stereo 3D. The workaround is to shrink the PyMOL window to cover the visible portion of the screen. It is worth noting that true 3D-capable LCDs (as distinct from 3D-capable HDTVs) do not suffer from this problem. 

Aside from the above concerns, the quality of the DLP stereo 3D effect is exceptional: there is absolutely no ghosting or cross-talk between the two images. 

  * [Samsung 3D-Ready DLP HDTVs](http://pages.samsung.com/us/dlp3d) \- work with PyMOL 1.2b3 & later without any special drivers. Quadro driver support is still lacking as of Feb. 1st, 2009 - WLD


  * [Mitsubishi 3D-Ready DLP HDTVs](http://www.mitsubishi-tv.com/) \- not yet tested, but are expected to work with PyMOL 1.2b3 & later without any special drivers. - WLD



See [The 3D HDTV List](http://www.3dmovielist.com/3dhdtvs.html) for more 3D-capable HDTV options. 

## Passive Stereo 3D Displays

### One Piece Multi-layer LCD Displays

Affordable! 

#### Zalman

**iZ3D, the original supplier of Zalman display drivers has ceased operation and support as of 31 July 2012. DO NOT PURCHASE THESE MONITORS WITHOUT FURTHER CONFIRMATION of display support, the iZ3D support (required drivers, etc) is not activatable. If you do have further information, please post it here.** [Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)") 12:21, 12 September 2012 (CDT) 

  * [Zalman 22-inch 3D LCD monitor](http://www.zalman.co.kr/eng/product/Product_read.asp?Idx=219) \- works with PyMOL 1.2b3 & later without any special drivers. Great stereo quality provided that all drawn lines are at least 2 pixels thick. Menus are a bit awkward to use while in stereo mode, but even so, this 650 USD display provides excellent 3D molecular visualization in both full-screen in windowed modes. - WLD (**The Zalman ZM-M220W is DeLano Scientific's RECOMMENDED SOLUTION as of Feb 11, 2009!**).
  * [Zalman 24-inch 3D LCD monitor](http://www.zalman.co.kr/Eng/product/Product_Read.asp?idx=391) \- also works with PyMOL 1.2b3 & later under LINUX (Centos 5 x86_64 plain kernel + NVidia driver from ELRepo). I'm using an NVidia Quadro FX 580 (G96GL) graphics card (£125). Monitor cost around £350. PyMOL automagically detects that quad buffered stereo is available on startup.--[Bosmith](/index.php?title=User:Bosmith&action=edit&redlink=1 "User:Bosmith \(page does not exist\)") 16:32, 2 December 2010 (UTC)



#### LG

  * [LG D2342P-PN](http://www.lg.com/us/computer-products/monitors/LG-led-monitor-D2342P-PN.jsp)
  * [LG DM2752D](http://www.lg.com/uk/support-product/lg-DM2752D-PZ) This (and other LG passive 3D monitors/TVs) work with PyMOL. I'm driving them using NVidia Quadro 600 graphics cards under LINUX (CentOSes 5, 6 & 7) with the NVidia driver from ELRepo.



##### On Linux

  * setup by editing the xorg.conf file:



in the Device section of xorg.conf add: 
    
    
       Option "Stereo" "7"
    

in the Screen section of xorg.conf and an additional: 
    
    
       Section "Extensions"
           Option         "Composite" "Disable"
       EndSection
    

N.B. the current Gnome 3 (gnome-shell) in RHEL 7 derivatives (Scientific, CentOS, etc.) is a compositing window manager and is not properly stereo aware, so you will need to use a different window manager. See e.g. [[1]](https://sbgrid.org/wiki/stereo) for an alternative. 

  * launch using:


    
    
       pymol -S -t 6
    

##### On Windows

  1. From the Start Menu: In the "PyMOL" folder, go into the "Stereo 3D Launch" subfolder, and select "PyMOL Zalman 3D (By Row)". You might want to control-drag a copy of that shortcut on to your desktop in order to drag & drop content files onto it for stereo 3D visualization
  2. From the Command Line:


    
    
       "C:\Program Files\DeLano Scientific\PyMOL\PyMOLWin.exe" -S -t 6
    

##### On Mac

  1. MacPyMOL: Copy and rename the "MacPyMOL" application bundle to "MacPyMOLZalman". You can then double-click on the MacPyMOLZalman icon or drop data files directly onto it to visualize content in the Zalman stereo 3D mode.
  2. PyMOL X11 Hybrid Mode: Copy and rename the "MacPyMOL" appplication bundle to "PyMOLX11Zalman". After launching X11, you can then double-click on the PyMOLX11Zalman icon or drop data files directly onto it to visualize that content in the Zalman stereo 3D mode.



_For all platforms, remember to toggle stereo on and off using the "set stereo" command:_
    
    
       set stereo, on
    

#### iZ3D

  * [IZ3D](http://www.iz3d.com) \- works with PyMOL 1.2b3 & later without any special drivers. However, this display exhibits far too much cross-talk and interference between the two stereo images. Not suitable for professional use. - WLD


  * IZ3D is closed as of 31 July 2012, and will not offer support to their products. [Jedgold](/index.php?title=User:Jedgold&action=edit&redlink=1 "User:Jedgold \(page does not exist\)") 12:26, 12 September 2012 (CDT)



### Mirror-based Multi-LCD Solutions

Expensive! 

  * [Planar3D](http://www.planar3d.com) "I have used these displays with nVidia Quadro graphics cards under both Windows and Linux running both PyMOL and Maestro. They work well, and the stereo quality is excellent!" - WLD.
  * [Omnia MIMO](http://www.inition.co.uk/inition/product.php?URL_=product_stereovis_omnia_mimo&SubCatID_=3)



## Autostereoscopic LCD Displays

Some autostereoscopic displays have the ability to switch between 2D and 3D display modes. Others are built for 3D only. 

  * [Dimension Technologies Inc.](http://www.dti3d.com)
  * [SeeReal Technologies](http://www.seereal.com)
  * [NewSight Corp.](http://www.newsight.com/3d-products/displays.html)



## Stereo 3D Projectors

Although these displays require shutter glasses out of the box, when combined with the adapters below and a special "silvered" screen, they can be used to project Passive Stereo 3D to a large audience. 

### Active Stereo 3D DLP Projectors

  * [DepthQ Stereoscopic](http://www.depthq.com) "The original DepthQ gave a very good stereo 3D effort with PyMOL, but I haven't seen their latest products." - WLD.
  * [Christie MIRAGE S+4K SXGA+ 6500 LUMEN DLP™ STEREOSCOPIC PROJECTOR](http://www.christiedigital.com/AMEN/Products/christieMirageS4K.htm) "I have been very impressed with the stereo 3D effect produced by MIRAGE projectors equipped with StereoGraphic ZScreens running PyMOL under Windows with a high-end nVidia Quadro card." - WLD.



### Passive Stereo 3D Adaptor Hardware for Active Stereo 3D Projectors

These devices make it possible for a large audience to see projected stereo 3D using inexpensive polarized glasses. 

  * [RealD StereoGraphics Projection ZScreen](http://reald-corporate.com/scientific/projectorzscreen.asp)



## Who Says What?

If you provide a specific quote or endorsement above, please say who you are so that everyone can know the source of the information. 

  * WLD = Warren L. DeLano of DeLano Scientific LLC
  * SP = Sabuj Pattanayek of the Center For Structural Biology, Vanderbilt University



Retrieved from "[https://pymolwiki.org/index.php?title=Stereo_3D_Display_Options&oldid=12808](https://pymolwiki.org/index.php?title=Stereo_3D_Display_Options&oldid=12808)"


---

## XFree86 Configuration

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

## Iiyama HM204DT + Nvidia Quadro FX 1100 on Fedora Core 2 (updated as of 2005-02-17)

` `

`
    
    
    Section "Module"
           Load  "dbe"
           Load  "extmod"
           Load  "fbdevhw"
           Load  "glx"
           Load  "record"
           Load  "freetype"
           Load  "type1"
           #Load  "dri"    # THIS MUST BE COMMENTED OUT
    EndSection
    

``
    
    
    Section "Monitor"
      Identifier   "hm204dt"
      VendorName   "Iiyama"
      ModelName    "HM204DT"
      HorizSync    30 - 142
      VertRefresh  50 - 200
      Option      "dpms"
    
      # 1400x1100 @ 120.00 Hz (GTF) hsync: 141.48 kHz; pclk: 275.04 MHz
      Modeline "1400x1100_120"  275.04  1400 1520 1672 1944  1100 1101 1104 1179  -HSync +Vsync
    EndSection
    

``
    
    
    Section "Device"
           Identifier  "quadrofx"
           Driver      "nvidia"
           VendorName  "PNY"
           BoardName   "NVIDIA Quadro FX 1100"
           Option "AllowDFPStereo" "true"
           Option "Stereo" "3"
    EndSection
    

``
    
    
    Section "Screen"
           Identifier "Screen0"
           Device     "quadrofx"
           Monitor    "hm204dt"
           DefaultDepth     24
           SubSection "Display"
                   Viewport   0 0
                   Depth     24
                   Modes   "1400x1100_120" "1280x1024"
           EndSubSection
    EndSection
    

```

``

Retrieved from "[https://pymolwiki.org/index.php?title=XFree86_Configuration&oldid=4398](https://pymolwiki.org/index.php?title=XFree86_Configuration&oldid=4398)"


---

## Category:Hardware Options

The printable version is no longer supported and may have rendering errors. Please update your browser bookmarks and please use the default browser print function instead.

Hardware Options Here 

## Pages in category "Hardware Options"

The following 2 pages are in this category, out of 2 total. 

### M

  * [Monitors Hardware Options](/index.php/Monitors_Hardware_Options "Monitors Hardware Options")



### S

  * [Stereo 3D Display Options](/index.php/Stereo_3D_Display_Options "Stereo 3D Display Options")



Retrieved from "[https://pymolwiki.org/index.php?title=Category:Hardware_Options&oldid=3063](https://pymolwiki.org/index.php?title=Category:Hardware_Options&oldid=3063)"


---

