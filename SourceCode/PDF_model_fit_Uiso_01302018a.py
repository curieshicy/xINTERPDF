""""
This is the code for OrganicPDF GUI program. It is used to extract structure information from a measured Organic PDF.
The code is written by Chenyang Shi at AbbVie.
Send him an email at chenyang.shi@abbvie.com for questions and suggestions. 

"""

from Tkinter import *
import Tkinter as tk
import ttk
from math import exp
import os  # for loading files or exporting files
import tkFileDialog
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

##DiffPy-CMI loading
from pyobjcryst.crystal import CreateCrystalFromCIF
from diffpy.srfit.pdf import DebyePDFGenerator, PDFParser, PDFGenerator, PDFContribution
from diffpy.srfit.fitbase import Profile, Calculator
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.Structure import Structure, loadStructure
from pyobjcryst.molecule import Molecule
from scipy.optimize import leastsq


import numpy

class quatzero(Calculator):
    """Return the zero-th component of normalized quaternion.
        """
    
    def __call__(self, q1, q2, q3):
        ssq = q1**2 + q2**2 + q3**2
        q0 = numpy.sqrt(1.0 - ssq) if ssq < 1 else 0.0
        return q0


class PDF_model_fit_Uiso():
    def __init__(self, master):
        self.master = master
        self.master.title("Model Fit to Organic PDF")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(900, 700) # width + height
        self.master.resizable(False, False)
        
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Arial", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')

        ttk.Label(self.master, text = "To perform a fit to crystalline organic PDF, one needs to prepare a structure file for the single molecule.\n"
            "In addition, a crystal structure for compound is needed as input.", font = ("Arial", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()
        
        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Molecule Structure", command = self.load_molecule,
                             style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Crystal Structure",command = self.load_crystal,
                             style = "TButton").grid(row = 1, column = 3, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Load PDF data",command = self.load_pdf_data,
                                        style = "TButton").grid(row = 1, column = 6, columnspan = 3, padx = 5)

    
        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ## these are the entries
        self.Qdamp = StringVar()
        self.Qbroad = StringVar()
        self.Uiso_intra = StringVar()
        self.Uiso_inter = StringVar()
        self.Qmin = StringVar()
        self.Qmax = StringVar()
        self.lat_a = StringVar()
        self.lat_b = StringVar()
        self.lat_c = StringVar()
        self.alpha = StringVar()
        self.beta  = StringVar()
        self.gamma  = StringVar()
        self.scale = StringVar()
        self.zoom_mol = StringVar()
        self.zoom_intra = StringVar()
        self.delta2_cryst = StringVar()
        self.delta2_mol  = StringVar()
        self.delta2_intra  = StringVar()
        self.rmin  = StringVar()
        self.rmax  = StringVar()
        self.optimizers = StringVar()
        
        ##set default values
        self.Qmin.set(0.0)
        self.zoom_mol.set(1.0)
        self.zoom_intra.set(1.0)
        self.delta2_cryst.set(0.0)
        self.delta2_mol.set(0.0)
        self.delta2_intra.set(0.0)
        self.rmin.set(1.0)
        self.rmax.set(20.0)
        
        ##set other default values just for test, remove later
        self.Qdamp.set(0.04)
        self.Qbroad.set(0.02)
        self.Uiso_intra.set(0.005)
        self.Uiso_inter.set(0.05)
        self.Qmax.set(24.0)  ##this Qmax is set to be same for both molecules

        self.scale.set(1)
        
        #these are the checkbuttons
        self.var_0  = IntVar()
        self.var_1  = IntVar()
        self.var_2  = IntVar()
        self.var_3  = IntVar()
        self.var_4  = IntVar()
        self.var_5  = IntVar()
        self.var_6  = IntVar()
        self.var_7  = IntVar()
        self.var_8  = IntVar()
        self.var_9  = IntVar()
        self.var_10 = IntVar()
        self.var_11 = IntVar()
        self.var_12 = IntVar()
        self.var_13 = IntVar()
        self.var_14 = IntVar()
        self.var_15 = IntVar()
        self.var_16 = IntVar()
        self.var_17 = IntVar()
        
        ##if values are fixed, then they are not refineable in the fit. Qdamp, Qbroad, Qmin, Qmax, three delta values should be fixed.
        self.var_0.set(1)
        self.var_1.set(1)
        self.var_4.set(1)
        self.var_5.set(1)
        self.var_13.set(1)
        self.var_14.set(1)
        self.var_15.set(1)
        self.var_16.set(1)
        self.var_17.set(1)

        #################first layer Qdamp, Qbroad, Uiso_inter, Uiso_intra, Qmin, Qmax##########################
        #Qdamp
        ttk.Label(self.top_frame, text = "Qdamp", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Qdamp).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_0,
                        style = "TCheckbutton").grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #Qbroad
        ttk.Label(self.top_frame, text = "Qbroad", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Arial", 12),
                     textvariable = self.Qbroad).grid(row = 3, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_1,
                        style = "TCheckbutton").grid(row = 3, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso intra
        ttk.Label(self.top_frame, text = "Uiso_intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Uiso_intra).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_2,
                        style = "TCheckbutton").grid(row = 3, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso inter
        ttk.Label(self.top_frame, text = "Uiso_inter", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Uiso_inter).grid(row = 3, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_3,
                        style = "TCheckbutton").grid(row = 3, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #Qmin
        ttk.Label(self.top_frame, text = "Qmin", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Qmin).grid(row = 3, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_4,
                        style = "TCheckbutton").grid(row = 3, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #Qmax
        ttk.Label(self.top_frame, text = "Qmax", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Qmax).grid(row = 3, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_5,
                        style = "TCheckbutton").grid(row = 3, column = 17, columnspan = 1, padx = 5, sticky = "sw")

        #################second layer lat_a, lat_b. lat_c, alpha, beta, gamma##########################        
        #lat_a
        ttk.Label(self.top_frame, text = "Lat_a", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.lat_a).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_6,
                        style = "TCheckbutton").grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #lat_b
        ttk.Label(self.top_frame, text = "Lat_b", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Arial", 12),
                     textvariable = self.lat_b).grid(row = 4, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_7,
                        style = "TCheckbutton").grid(row = 4, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #lat_c
        ttk.Label(self.top_frame, text = "Lat_c", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.lat_c).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_8,
                        style = "TCheckbutton").grid(row = 4, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #alpha
        ttk.Label(self.top_frame, text = "Alpha", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.alpha).grid(row = 4, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_9,
                        style = "TCheckbutton").grid(row = 4, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #beta
        ttk.Label(self.top_frame, text = "Beta", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.beta).grid(row = 4, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_10,
                        style = "TCheckbutton").grid(row = 4, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #gamma
        ttk.Label(self.top_frame, text = "Gamma", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.gamma).grid(row = 4, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_11,
                        style = "TCheckbutton").grid(row = 4, column = 17, columnspan = 1, padx = 5, sticky = "sw")


        #################third layer scale, zoom_mol, zoom_intra, delta2_cryst, delta2_mole, delta2_intra##########################        
        #scale
        ttk.Label(self.top_frame, text = "Scale", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.scale).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_12,
                        style = "TCheckbutton").grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #zoom_mol
        ttk.Label(self.top_frame, text = "Zoom_Mol", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Arial", 12),
                     textvariable = self.zoom_mol).grid(row = 5, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_13,
                        style = "TCheckbutton").grid(row = 5, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #zoom_intra
        ttk.Label(self.top_frame, text = "Zoom_Intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.zoom_intra).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_14,
                        style = "TCheckbutton").grid(row = 5, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_cryst
        ttk.Label(self.top_frame, text = "Delta2_Cryst", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.delta2_cryst).grid(row = 5, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_15,
                        style = "TCheckbutton").grid(row = 5, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_mol
        ttk.Label(self.top_frame, text = "Delta2_Mol", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.delta2_mol).grid(row = 5, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_16,
                        style = "TCheckbutton").grid(row = 5, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_intra
        ttk.Label(self.top_frame, text = "Delta2_Intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.delta2_intra).grid(row = 5, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_17,
                        style = "TCheckbutton").grid(row = 5, column = 17, columnspan = 1, padx = 5, sticky = "sw")

        ## layout for step 3, Select Optimizers and run

        ttk.Label(self.top_frame, text = "Step 3: Select Optimzers", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        #rmin
        ttk.Label(self.top_frame, text = "Rmin", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                      textvariable = self.rmin).grid(row = 7, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                  
        #rmax
        ttk.Label(self.top_frame, text = "Rmax", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                      textvariable = self.rmax).grid(row = 7, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        self.combobox = ttk.Combobox(self.top_frame, font = ("Arial", 12), textvariable = self.optimizers, values = ["Least Squares", "Fmin",
                              "Basin Hopping"])
        self.combobox.current(0)
        self.combobox.grid(row = 7, column = 4, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.top_frame, text = "Run The FIT", command = self.run_the_fit,
                             style = "TButton").grid(row = 7, column = 6, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Reset To Default", command = self.reset_to_default,
                             style = "TButton").grid(row = 7, column = 9, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "***Fit Use Rigid Body***", command = self.rigid_body_fit,
                       style = "TButton").grid(row = 7, column = 11, columnspan = 3, padx = 5, sticky = "sw")

    ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.bottom_frame = ttk.Frame(self.master, padding = (10, 10))
        self.bottom_frame.pack()

        self.fig = plt.figure(figsize=(12, 5), dpi=100) ##create a figure; modify the size here
        self.fig.add_subplot()
        
        plt.title("PDF Model Fit To Crystalline Organic PDF")
        plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        
        #        plt.xticks([])
        #        plt.yticks([])
        
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.bottom_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.bottom_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    
    def load_molecule(self):
        self.input_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                      filetypes = [("Text Documents", "*.xyz")])

    def load_crystal(self):
        self.input_cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".cif",
                                                        filetypes = [("Text Documents", "*.cif")])
                                                        
        self.lat_a.set(loadStructure(self.input_cryst_struct).lattice.a)
        self.lat_b.set(loadStructure(self.input_cryst_struct).lattice.b)
        self.lat_c.set(loadStructure(self.input_cryst_struct).lattice.c)
        self.alpha.set(loadStructure(self.input_cryst_struct).lattice.alpha)
        self.beta.set(loadStructure(self.input_cryst_struct).lattice.beta)
        self.gamma.set(loadStructure(self.input_cryst_struct).lattice.gamma)

    
    def load_pdf_data(self):
        self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                         filetypes = [("Text Documents", "*.gr")])
    
    def reset_to_default(self):
        ##checkbox
        self.var_0.set(1)
        self.var_1.set(1)
        self.var_4.set(1)
        self.var_5.set(1)
        self.var_13.set(1)
        self.var_14.set(1)
        self.var_15.set(1)
        self.var_16.set(1)
        self.var_17.set(1)
        ##entry
        self.Qmin.set(0.0)
        self.zoom_mol.set(1.0)
        self.zoom_intra.set(1.0)
        self.delta2_cryst.set(0.0)
        self.delta2_mol.set(0.0)
        self.delta2_intra.set(0.0)
        self.rmin.set(1.0)
        self.rmax.set(20.0)
        self.combobox.current(0)

        self.lat_a.set(loadStructure(self.input_cryst_struct).lattice.a)
        self.lat_b.set(loadStructure(self.input_cryst_struct).lattice.b)
        self.lat_c.set(loadStructure(self.input_cryst_struct).lattice.c)
        self.alpha.set(loadStructure(self.input_cryst_struct).lattice.alpha)
        self.beta.set(loadStructure(self.input_cryst_struct).lattice.beta)
        self.gamma.set(loadStructure(self.input_cryst_struct).lattice.gamma)
    
    ###########This is an attempt to perform a fit using Acta A paper approach##########
    def rigid_body_fit(self):

        def getCrystalWithMolecule(ciffilename):
            """
            Load the CIF file and group all scatterers as one Molecule.

            This assumes the asymmetric unit is a single whole Molecule,
            i.e., there are no other molecules or independent atoms in
            the crystal and the molecule does not need to be
            symmetry-expanded.
            """
            # center of the molecule in fractional coordinates
            xyzmol = numpy.array([0.0, 0.0, 0.0])
            with open(ciffilename) as fp:
                crst = CreateCrystalFromCIF(fp, oneScatteringPowerPerElement=True)
            mol = Molecule(crst, "mol")
            for i in reversed(range(crst.GetNbScatterer())):
                a = crst.GetScatterer(i)
                xyzm = numpy.array([a.X, a.Y, a.Z]) - xyzmol
                xc, yc, zc = crst.FractionalToOrthonormalCoords(*xyzm)
                mol.AddAtom(xc, yc, zc, a.GetScatteringPower(), a.GetName())
                crst.RemoveScatterer(a)
            crst.AddScatterer(mol)
            mol.X, mol.Y, mol.Z = xyzmol
            return crst

#        def plotrecipe(f, what='obs, calc, diff'):
#            from matplotlib.pyplot import plot
#            choice = set(what.replace(',', ' ').split())
#            
#            r = f.pcnt.r.value
#            gobs = f.pcnt.y.value
#            gcalc = f.pcnt.evaluate()
#            gdiff = gobs - gcalc
#            yall = [('obs', gobs),
#            ('calc', gcalc),
#            ('diff', gdiff + g_diff_baseline)]
#            pargs = sum([(r, y) for n, y in yall if n in choice], ())
#            rv = plot(*pargs)
#            if 'calc' in choice and 'diff' in choice:
#                rv[-1].set_color(rv[-2].get_color())
#            return rv

        def plotrecipe(f, c1 = "r-", c2 = "g-", c3 = "k-"):
    
            self.fig.clf()
            r = f.pcnt.r.value
            g = f.pcnt.y.value
            gcalc = f.pcnt.evaluate()
            diffzero = -0.8 * max(g) * numpy.ones_like(g)
            diff = g - gcalc + diffzero
            
            self.fig.add_subplot(111)
            plt.plot(r,g,'bo',markersize = 5, label="G(r) Data")
            plt.plot(r, gcalc,c1,label="G(r) Fit", lw= 2)
            plt.plot(r,diff,c2,label="G(r) diff", lw = 2)
            plt.plot(r,diffzero,c3, lw = 2)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.legend(loc=1)
            
            plt.title("PDF Model Fit Using a Rigid Body Approach")
            plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.canvas.draw()
            
            return



        # define PDFContribution compososed of inter and intra-molecular components.
        pcnt = PDFContribution('pcnt')
        pcnt.loadData(self.input_pdf_gr)
        pcnt.qdamp = float(self.Qdamp.get())
        pcnt.qbroad = float(self.Qbroad.get())
        pcnt.setCalculationRange(xmin= float(self.rmin.get()), xmax = float(self.rmax.get()), dx = 0.01)
        # intra-molecular contribution from one molecule:
        crst = getCrystalWithMolecule(self.input_cryst_struct)
        mol = crst.GetScatterer("mol")
        pcnt.addStructure('mol', mol, periodic=False)
        # inter-molecular contributions (wide)
        wcrst = getCrystalWithMolecule(self.input_cryst_struct)
        wmol = wcrst.GetScatterer('mol')
        pcnt.addStructure('wcrst', wcrst, periodic=True)
        pcnt.addStructure('wmol', wmol, periodic=False)
        pcnt.setEquation('scale * (mol + wcrst - wmol)')

        # speed up simulation by using parallel jobs.
        pcnt.wcrst.parallel(4)

        # define FitRecipe for this organic crystal
        ocfit = FitRecipe()
        ocfit.clearFitHooks()
        ocfit.addContribution(pcnt)

        # expose the overall data scale
        ocfit.addVar(ocfit.pcnt.scale, name='scale', value=0.1)
        ocfit.addVar(ocfit.pcnt.qdamp, name='qdamp')

        # expose unit cell parameters
        pwcrst = ocfit.pcnt.wcrst.phase
        ocfit.addVar(pwcrst.a)
        ocfit.addVar(pwcrst.b)
        ocfit.addVar(pwcrst.c)
        # objcryst uses radians for cell angle beta, but let us refine
        # in degrees.
        ocfit.newVar('beta', value=numpy.degrees(pwcrst.beta.value))
        ocfit.constrain(pwcrst.beta, 'radians(beta)')

        # objcryst links all atoms of the same type to one scattering power
        # with common Biso.  It is thus sufficient to constrain just one of
        # those.

        # first let's ensure there are only 2 independent scattering powers
        # for "C" and "H" in the CIF file.
        assert 2 == crst.GetScatteringPowerRegistry().GetNb(), \
        "unexpected number of atom species"

        # expose isotropic displacement parameters for carbon and hydrogen
        # in the molecule.  CIF file has sites labeled "C1", "H1" so we can
        # refer to them by those names.
        pmol = ocfit.pcnt.mol.phase
        ocfit.addVar(pmol.C1.Biso, name='bisoC', value=1.0)
        ocfit.restrain('bisoC', lb=0.02, sig=1e-4)
        # hydrogen gives negligible contribution, so we keep them constant
        ocfit.addVar(pmol.H1.Biso, name='bisoH', value=1, fixed=True)

        # for intra-molecular component we use the same binter value for all atoms:
        ocfit.addVar(pwcrst.mol.C1.Biso, name='binter', value=10)
        ocfit.constrain(pwcrst.mol.H1.Biso, 'binter')

        # expose rotation quaternions for the molecule
        pwmol = ocfit.pcnt.wcrst.phase.mol
        ocfit.addVar(pwmol.q1, name='wq1')
        ocfit.addVar(pwmol.q2, name='wq2')
        ocfit.addVar(pwmol.q3, name='wq3')

        # constrain q0 so we have normalized quaternion.
        # first define a function "fixq0" that calculates q0 from [q1, q2, q3]
        ocfit.registerCalculator(quatzero('fixq0'), argnames=[])
        ocfit.constrain(pwmol.q0, 'fixq0(wq1, wq2, wq3)')

        ocfit.fix('all')
        ocfit.free('scale')
        ocfit.scalarResidual()
        print("== INITIAL ==\n")
        print(FitResults(ocfit))

        ocfit.free('a', 'b', 'c')
        ocfit.free('bisoC', 'binter')
        leastsq(ocfit.residual, ocfit.values)
        print("== INTERMEDIATE ==\n")
        print(FitResults(ocfit))

        ocfit.free('wq1', 'wq2', 'wq3')
        leastsq(ocfit.residual, ocfit.values)
        print("== REFINED ==\n")
        print(FitResults(ocfit))

        # rotate molecule away from the minimum
        ocfit.wq1 << 0.1
        ocfit.fix('a', 'b', 'c', 'binter', 'wq1', 'wq2', 'wq3')
        leastsq(ocfit.residual, ocfit.values)
        print("== MISORIENTED ==\n")
        print(FitResults(ocfit))
        plotrecipe(ocfit)

        # and finally refine its orientation
        ocfit.free('wq1', 'wq2', 'wq3')
        leastsq(ocfit.residual, ocfit.values)
        print("== ORIENTED BACK ==\n")
        print(FitResults(ocfit))
        plotrecipe(ocfit, c1 = "m-", c2 = "y-", c3 = "c-")
        plt.show()

    def run_the_fit(self):
        
        ####codes for PDF fits
        def makeRecipe(stru1, stru2, stru3, datname):
            ## The Profile
            profile = Profile()
            
            # Load data and add it to the profile
            parser = PDFParser()
            parser.parseFile(datname)
            profile.loadParsedData(parser)
            #profile.setScatteringType("N")
            profile.setCalculationRange(xmin= float(self.rmin.get()), xmax = float(self.rmax.get()), dx = 0.01)
            
            ## The ProfileGenerator
            generator_MEF_Cryst_B = PDFGenerator("G_MEF_Cryst_B")
            generator_MEF_Cryst_B.setStructure(stru1, periodic = True)
            generator_MEF_Mole_B = DebyePDFGenerator("G_MEF_Mole_B")
            generator_MEF_Mole_B.setStructure(stru2, periodic = False)
            generator_MEF_Intra = DebyePDFGenerator("G_MEF_Intra")
            generator_MEF_Intra.setStructure(stru3, periodic = False)
            
            ## The FitContribution
            # Add both generators to the FitContribution. Add the Profile. This will
            # send the metadata to the generators.
            contribution = FitContribution("MEF")
            contribution.addProfileGenerator(generator_MEF_Cryst_B)
            contribution.addProfileGenerator(generator_MEF_Mole_B)
            contribution.addProfileGenerator(generator_MEF_Intra)
            contribution.setProfile(profile, xname = "r")
            #write down the fit equation:
            #(G_MEF_Cryst_B - G_MEF_Mole_B) gives the intermolecular PDF, using a larger atomic displacement parameter
            #G_MEF_Intra gives intramolecular PDF, using a smaller atomic displacement parameter.
            #The sum of both parts gives the total PDF.
            contribution.setEquation("scale * (G_MEF_Cryst_B - G_MEF_Mole_B + G_MEF_Intra)")
            
            # Make the FitRecipe and add the FitContribution.
            recipe = FitRecipe()
            recipe.addContribution(contribution)

            ##############self.var_0 is Qdamp#####
            recipe.newVar("Qdamp", float(self.Qdamp.get()))
            recipe.newVar("Qbroad", float(self.Qbroad.get()))

            if int(self.var_0.get()) == 1: ##1 is fixed; 0 to free parameters
                generator_MEF_Cryst_B.qdamp.value = float(self.Qdamp.get())
                generator_MEF_Mole_B.qdamp.value = float(self.Qdamp.get())
                generator_MEF_Intra.qdamp.value = float(self.Qdamp.get())

            elif int(self.var_0.get()) == 0:
                recipe.constrain(generator_MEF_Cryst_B.qdamp, "Qdamp")
                recipe.constrain(generator_MEF_Mole_B.qdamp, "Qdamp")
                recipe.constrain(generator_MEF_Intra.qdamp, "Qdamp")
           
           ##############self.var_1 is Qbroad#####

            if int(self.var_1.get()) == 1: ##1 is fixed; 0 to free parameters
                generator_MEF_Cryst_B.qbroad.value = float(self.Qbroad.get())
                generator_MEF_Mole_B.qbroad.value = float(self.Qbroad.get())
                generator_MEF_Intra.qbroad.value = float(self.Qbroad.get())
            
            elif int(self.var_1.get()) == 0:
                recipe.constrain(generator_MEF_Cryst_B.qbroad, "Qbroad")
                recipe.constrain(generator_MEF_Mole_B.qbroad, "Qbroad")
                recipe.constrain(generator_MEF_Intra.qbroad, "Qbroad")
            
            # Vary the gloabal scale as well.
            
            recipe.newVar("Global_Scale", float(self.scale.get()))
            
            if int(self.var_12.get()) == 1: ##1 is fixed; 0 to free parameters
                contribution.scale = float(self.scale.get())
            
            elif int(self.var_12.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(contribution.scale, "Global_Scale")
            

            #############################################################################################
            ############### First the MEF_Cryst_B parameters ############################################
            #############################################################################################
            phase_MEF_Cryst_B = generator_MEF_Cryst_B.phase
            
            lat = phase_MEF_Cryst_B.getLattice()
            atoms = phase_MEF_Cryst_B.getScatterers()
            
            recipe.newVar("Uiso_Inter", float(self.Uiso_inter.get()))
            recipe.newVar("lat_a", float(self.lat_a.get()))
            recipe.newVar("lat_b", float(self.lat_b.get()))
            recipe.newVar("lat_c", float(self.lat_c.get()))
            recipe.newVar("alpha", float(self.alpha.get()))
            recipe.newVar("beta", float(self.beta.get()))
            recipe.newVar("gamma", float(self.gamma.get()))

            if int(self.var_6.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.a = float(self.lat_a.get())
            elif int(self.var_6.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.a, "lat_a")

            if int(self.var_7.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.b = float(self.lat_b.get())
            elif int(self.var_7.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.b, "lat_b")

            if int(self.var_8.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.c = float(self.lat_c.get())
            elif int(self.var_8.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.c, "lat_c")
            
            if int(self.var_9.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.alpha = float(self.alpha.get())
            elif int(self.var_9.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.alpha, "alpha")

            if int(self.var_10.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.beta = float(self.beta.get())
            elif int(self.var_10.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.beta, "beta")
            
            if int(self.var_11.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.gamma = float(self.gamma.get())
            elif int(self.var_11.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.gamma, "gamma")

            
            if int(self.var_3.get()) == 1: ##1 is fixed; 0 to free parameters
                for atom in atoms:
                    atom.Uiso = float(self.Uiso_inter.get())
            #atom.Uiso = float(self.Uiso_inter.get())
            
            elif int(self.var_3.get()) == 0: ##1 is fixed; 0 to free parameters
                for atom in atoms:
                    recipe.constrain(atom.Uiso, "Uiso_Inter")

            recipe.newVar("delta2_cryst", float(self.delta2_cryst.get()))

            if int(self.var_15.get()) == 1: ##1 is fixed; 0 to free parameters
                generator_MEF_Cryst_B.delta2 = float(self.delta2_cryst.get())
                #recipe.newVar(generator_MEF_Cryst_B.delta2, float(self.delta2_cryst.get()), fixed = True)

            elif int(self.var_15.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(generator_MEF_Cryst_B.delta2, "delta2_cryst")
            
            
            #############################################################################################
            ############### Second the MEF_Mole_B parameters ############################################
            #############################################################################################
            phase_MEF_Mole_B = generator_MEF_Mole_B.phase
            generator_MEF_Mole_B.setQmin(float(self.Qmin.get()))
            generator_MEF_Mole_B.setQmax(float(self.Qmax.get()))

            lat = phase_MEF_Mole_B.getLattice()
            recipe.newVar("zoom_Mole_B", float(self.zoom_mol.get()))

            if int(self.var_13.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.a = float(self.zoom_mol.get())
                lat.b = float(self.zoom_mol.get())
                lat.c = float(self.zoom_mol.get())

            elif int(self.var_13.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.a, "zoom_Mole_B")
                recipe.constrain(lat.b, "zoom_Mole_B")
                recipe.constrain(lat.c, "zoom_Mole_B")

            # Constrain fractional xyz parameters
            atoms = phase_MEF_Mole_B.getScatterers()
            # Constrain ADPs
            
            if int(self.var_3.get()) == 1: ##1 is fixed; 0 to free parameters
                for atom in atoms:
                    atom.Uiso = float(self.Uiso_inter.get())

            elif int(self.var_3.get()) == 0: ##1 is fixed; 0 to free parameters
                for atom in atoms:
                    recipe.constrain(atom.Uiso, "Uiso_Inter")

            recipe.newVar("delta2_mol", float(self.delta2_mol.get()))

            if int(self.var_16.get()) == 1: ##1 is fixed; 0 to free parameters
                generator_MEF_Mole_B.delta2 = float(self.delta2_mol.get())

            elif int(self.var_16.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(generator_MEF_Mole_B.delta2, "delta2_mol")

            #############################################################################################
            ############### Third the intra molecule parameters##########################################
            #############################################################################################
            phase_MEF_Intra = generator_MEF_Intra.phase
            generator_MEF_Intra.setQmin(float(self.Qmin.get()))
            generator_MEF_Intra.setQmax(float(self.Qmax.get()))

            lat = phase_MEF_Intra.getLattice()
            recipe.newVar("zoom_Intra", float(self.zoom_intra.get()))


            if int(self.var_14.get()) == 1: ##1 is fixed; 0 to free parameters
                lat.a = float(self.zoom_intra.get())
                lat.b = float(self.zoom_intra.get())
                lat.c = float(self.zoom_intra.get())

            elif int(self.var_14.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(lat.a, "zoom_Intra")
                recipe.constrain(lat.b, "zoom_Intra")
                recipe.constrain(lat.c, "zoom_Intra")

            # Constrain fractional xyz parameters
            atoms = phase_MEF_Intra.getScatterers()
            # Constrain ADPs
            recipe.newVar("Uiso_Intra", float(self.Uiso_intra.get()))

            
            if int(self.var_2.get()) == 1: ##1 is fixed; 0 to free parameters
                for atom in atoms:
                    #recipe.addVar(atom.Uiso, float(self.Uiso_intra.get()), fixed = True)
                    atom.Uiso = float(self.Uiso_intra.get())
                        
            elif int(self.var_2.get()) == 0:
                for atom in atoms:
                    recipe.constrain(atom.Uiso, "Uiso_Intra")

            recipe.newVar("delta2_intra", float(self.delta2_intra.get()))

            if int(self.var_17.get()) == 1: ##1 is fixed; 0 to free parameters
                generator_MEF_Intra.delta2 = float(self.delta2_intra.get())
            elif int(self.var_17.get()) == 0: ##1 is fixed; 0 to free parameters
                recipe.constrain(generator_MEF_Intra.delta2, "delta2_intra")

            # Give the recipe away so it can be used!
            return recipe
                        
                        
        def plotResults(recipe):
 
            self.fig.clf()
            r = recipe.MEF.profile.x
            g = recipe.MEF.profile.y
            gcalc = recipe.MEF.profile.ycalc
            diffzero = -0.8 * max(g) * numpy.ones_like(g)
            diff = g - gcalc + diffzero

            self.fig.add_subplot(111)
            plt.plot(r,g,'bo',markersize = 5, label="G(r) Data")
            plt.plot(r, gcalc,'r-',label="G(r) Fit", lw= 2)
            plt.plot(r,diff,'g-',label="G(r) diff", lw = 2)
            plt.plot(r,diffzero,'k-', lw = 2)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.legend(loc=1)

            plt.title("PDF Model Fit To Crystalline Organic PDF")
            plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.canvas.draw()

            return

        
        data  = self.input_pdf_gr
        stru1 = Structure(filename = self.input_cryst_struct)
        stru2 = Structure(filename = self.input_mole_struct)
        stru3 = Structure(filename = self.input_mole_struct)
        recipe = makeRecipe(stru1, stru2, stru3, data)
        
#        _plot = True
#        if _plot:
#            from diffpy.srfit.fitbase.fithook import PlotFitHook
#            recipe.pushFitHook(PlotFitHook())
#        recipe.fithooks[0].verbose = 3

        from diffpy.srfit.fitbase.fithook import PlotFitHook
        recipe.pushFitHook(PlotFitHook())
        recipe.fithooks[0].verbose = 3

        ###different optimizers###
        from scipy.optimize import leastsq
        from scipy.optimize import fmin
        from scipy.optimize import basinhopping
        from scipy.optimize import differential_evolution
        
        if self.optimizers.get() == "Least Squares":
            leastsq(recipe.residual, recipe.values)
        
        elif self.optimizers.get() == "Fmin":
            fmin(recipe.scalarResidual, recipe.getValues())

        elif self.optimizers.get() == "Basin Hopping":
            basinhopping(recipe.scalarResidual, recipe.getValues())

#        elif self.optimizers.get() == "Differential Evolution":
#            differential_evolution(recipe.scalarResidual, recipe.getValues())

        plotResults(recipe)

        basename = str(self.input_pdf_gr)[:-3]
        stru1.write(basename + "_cryst.cif", "cif")
        stru2.write(basename +  "_mole.xyz", "xyz")
        stru3.write(basename +  "_intra.xyz", "xyz")

        profile = recipe.MEF.profile
        
        profile.savetxt(basename + ".fit")

        res = FitResults(recipe)
        res.printResults()

        res.saveResults(basename + ".res")

#        if _plot:
#            plotResults(recipe)

        plotResults(recipe)

#        p2, C, info, msg, success = leastsq(recipe.residual, recipe.values, full_output = 1)
        
#        while info["nfev"]< 2000:
#            res = FitResults(recipe)
#            res.printResults()
#
#            plotResults(recipe)


def main():
    root = Tk()
    GUI = PDF_model_fit_Uiso(root)
    root.mainloop()

if __name__ == "__main__": main()
