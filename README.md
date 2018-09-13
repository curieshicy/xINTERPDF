# xINTERPDF 

(A video demo about installation and usage is available at https://www.youtube.com/watch?v=lAFZ5VYEH1g).

Python GUI program for analyzing organic pair distribution function (PDF) data collected at synchrotron and/or laboratory X-ray sources. It uses <b>DiffPy-CMI</b> (http://www.diffpy.org/products/diffpycmi/index.html) as a backend for simulation of PDFs. Currently it supports (1) The study of intermolecular interaction (e.g. hydrogen bonds) by subtracting out the scattering signal of single molecule(s) in real space. (2) The PDF model fit of the crystalline organic compound using the method proposed by Prill et al. (<i>J. Appl. Cryst.</i>, 48, 171-178, 2015.) (3) The phase quantification of physical mixtures of organics. (4) Generate Score/Scree plots based on Principle Component Analysis. 

## Citation

If you like the program and use it in your own work, it would be appreciated if you cite the following papers.

(1) Chenyang Shi, “xINTERPDF: a GUI program for analyzing intermolecular pair distribution functions of organic compounds from X-ray total scattering data”, <i>J. Appl. Cryst.</i>, 51, (2018), https://doi.org/10.1107/S1600576718012359.
<br>
(2) Pavol Juhás, Christopher L. Farrow, Xiaohao Yang, Kevin R. Knox, and Simon J. L. Billinge, “Complex modeling: a strategy and software program for combining multiple information sources to solve ill posed structure and nanostructure inverse problems”, <i>Acta Crystallogr. A</i>, 71, 562-568, 2015. 

## Overview
A quick view of the program is shown below.
### The Main Window

![overview](https://user-images.githubusercontent.com/8492535/41437338-bbb61878-6fe9-11e8-9dca-3556858c97ed.png)

### The interface for study of intermolecular PDF

![interpdf](https://user-images.githubusercontent.com/8492535/35756647-b06a457c-0831-11e8-82b4-7d6ef6c39178.png)

### The window for visualization of intermolecular PDF 

![screen shot 2018-06-16 at 3 14 36 pm](https://user-images.githubusercontent.com/8492535/41502096-15a816f8-7178-11e8-9ea3-3e4843fbed82.png)

### The interface for PDF model fit of measured organic crystalline PDF

![screen shot 2018-02-13 at 11 07 39 am](https://user-images.githubusercontent.com/8492535/36172498-afd203dc-10cb-11e8-8514-1654050e5d32.png)

### A breakdown of the total fit to the organic crystalline PDF

![e3a](https://user-images.githubusercontent.com/8492535/41437344-be14db18-6fe9-11e8-9df0-4793eda34031.png)

### The interface for phase quantification

![e4b](https://user-images.githubusercontent.com/8492535/41437457-2132cd0e-6fea-11e8-9e44-0400aa0f72bd.png)

### Generation of Score/Scree plots from Principle Component Analysis

![e5a](https://user-images.githubusercontent.com/8492535/41437348-c10cd4d8-6fe9-11e8-8261-2f6fbd832a6c.png)

## Installation

(See also Slides 11-14 in User Guide available at https://github.com/curieshicy/xINTERPDF/blob/master/xINTERPDF_User_Guide_20180615.pdf)

xINTERPDF can be installed on Linux and macOS computers. The easiest way to install it is through <b>conda</b>. Here is an example of installing it on macOS 10.10.3. 

(1) Download Anaconda Distribution for macOS at https://www.anaconda.com/download/?lang=en-us#macos. Select Python 2.7 version to install.

(2) Invoke a terminal, type <b>conda config --get channels</b> to check any channels that have been added. <b>diffpy</b> is required. If you don’t see it, type <b>conda config --add channels diffpy</b> to add it.  

(3) Type <b>conda create –c curieshicy –n xinterpdf xinterpdf</b> to install the xINTERPDF program.

(4) Once the installation is complete. Type <b>source activate xinterpdf</b> to start the virtual environment and <b>xinterpdf</b> to invoke the main window of xINTERPDF

(Alternatively) If conda install fails, one may download the raw files (Logo.gif and cli.py) at https://github.com/curieshicy/xINTERPDF/tree/master/Conda_Recipe_macOS_Linux/Conda_Build_Recipe_macOS/xinterpdf (macOS) or https://github.com/curieshicy/xINTERPDF/tree/master/Conda_Recipe_macOS_Linux/Conda_Build_Recipe_Linux/xinterpdf (Linux). To start the program, in a terminal, navigate to the folder where you put both files, and type <b>python cli.py</b> to invoke the main window. Make sure you have installed Diffpy-CMI, matplotlib (2.0.2) and Scikit-Learn (0.19.1). Follow http://www.diffpy.org/products/diffpycmi/index.html to install DiffPy-CMI. If you have conda, matplotlib can be installed by <b>conda install matplotlib=2.0.2</b>; Scikit-Learn can be installed by <b>conda install scikit-learn=0.19.1</b>  

## Manual

A user guide for xINTERPDF can be found at https://github.com/curieshicy/xINTERPDF/blob/master/xINTERPDF_User_Guide_20180615.pdf. It contains an overview of the program, a step-by-step installation guide, an explaination of techincal terms and a demonstration on three examples. The example files are available at https://github.com/curieshicy/xINTERPDF/tree/master/Examples. 

## Questions and suggestions

Questions and suggestions are welcome. You may contact me at cs3000@columbia.edu. Thanks for considering using xINTERPDF.









