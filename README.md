# xINTERPDF
Python GUI program for analyzing organic pair distribution function (PDF) data collected at synchrotron and laboratory X-ray sources. It uses <b>DiffPy-CMI</b> (http://www.diffpy.org/products/diffpycmi/index.html) as a backend for simulation of PDFs. 

## Citation

If you like the program and use it in your own work, please cite the following papers.

(1) Chenyang Shi, “xINTERPDF: a GUI program for analyzing intermolecular pair distribution functions in organic compounds from X-ray total scattering data”, under review, <i>J. Appl. Cryst.</i>, 2018.

(2) Pavol Juhás, Christopher L. Farrow, Xiaohao Yang, Kevin R. Knox, and Simon J. L. Billinge, <i>Acta Crystallogr., A</i> 71, 562-568, 2015. 


A quick view of the program is shown below.
### The Main Window

<img width="977" alt="main" src="https://user-images.githubusercontent.com/8492535/35762628-89f81c6a-085f-11e8-9841-d87f83060de4.png">

### The interface for study of intermolecular PDF

![interpdf](https://user-images.githubusercontent.com/8492535/35756647-b06a457c-0831-11e8-82b4-7d6ef6c39178.png)

### The window for visualization of intermolecular PDF 

![visual](https://user-images.githubusercontent.com/8492535/35756649-b0909b3c-0831-11e8-8d73-3369eef25364.png)

### The interface for PDF model fit of measured organic crystalline PDF

![screen shot 2018-02-13 at 11 07 39 am](https://user-images.githubusercontent.com/8492535/36172498-afd203dc-10cb-11e8-8514-1654050e5d32.png)

### The interface for phase quantification

![screen shot 2018-04-04 at 4 01 20 pm](https://user-images.githubusercontent.com/8492535/38377317-432135ba-38c0-11e8-823f-12a84ffde7c5.png)

## Installation

(See also Slides 9-13 in User Guide available at https://github.com/curieshicy/xINTERPDF/blob/master/xINTERPDF_User_Guide_0.1.0.pdf)

xINTERPDF can be installed on Linux and macOS computers. The easiest way to install it is through <b>conda</b>. Here is an example of installing it on macOS 10.10.3. 

(1) Download Anaconda Distribution for macOS at https://www.anaconda.com/download/?lang=en-us#macos. Select Python 2.7 version to install.

(2) Invoke a terminal, type <b>conda config --get channels</b> to check any channels that have been added. diffpy and xinterpdf are required. If you don’t see both, type <b>conda config --add channels diffpy</b> and <b>conda config --add channels xinterpdf</b> to add them.  

(3) Type <b>conda create –c curieshicy –n xinterpdf xinterpdf</b> to install it.

(4) Once the installation is complete. Type <b>source activate xinterpdf</b> to start the virtual environment and <b>xinterpdf</b> to invoke the main window of xINTERPDF

(Alternatively) If conda install fails, one may download the raw files (Logo.gif and cli.py) at https://github.com/curieshicy/xINTERPDF/tree/master/Conda_Build_Recipe/xinterpdf. To start the program, in a terminal, navigate the folder where you put both files, and type <b>python cli.py</b> to invoke the main window. Make sure you have installed Diffpy-CMI and matplotlib (2.0.2). Follow http://www.diffpy.org/products/diffpycmi/index.html to install DiffPy-CMI. If you have conda, matplotlib can be installed by <b>conda install matplotlib=2.0.2</b>.  

## Manual

A user guide for xINTERPDF can be found at https://github.com/curieshicy/xINTERPDF/blob/master/xINTERPDF_User_Guide_0.1.0.pdf. It contains an overview of the program, a step-by-step installation guide, an explaination of techincal terms and a demonstration on three examples. The example files are available at https://github.com/curieshicy/xINTERPDF/tree/master/Examples. 

## Questions and suggestions

Questions and suggestions are welcome. You may contact me at cs3000@columbia.edu. Thanks for considering using xINTERPDF.









