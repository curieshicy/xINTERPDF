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
from diffpy.srfit.pdf import DebyePDFGenerator, PDFParser, PDFGenerator
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.Structure import Structure

import numpy


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
        
        ##set default values
        self.Qmin.set(0.0)
        self.zoom_mol.set(1.0)
        self.zoom_intra.set(1.0)
        self.delta2_cryst.set(0.0)
        self.delta2_mol.set(0.0)
        self.delta2_intra.set(0.0)
        self.rmin.set(1.0)
        self.rmax.set(15.0)
        
        ##set other default values just for test, remove later
        self.Qdamp.set(0.02902)
        self.Uiso_intra.set(0.005)
        self.Uiso_inter.set(0.05)
        self.Qmax.set(24.0)
        self.lat_a.set(14.556)
        self.lat_b.set(6.811)
        self.lat_c.set(7.657)
        self.alpha.set(119.57)
        self.beta.set(103.93)
        self.gamma.set(91.30)
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
        #Uiso inter
        ttk.Label(self.top_frame, text = "Uiso_intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Uiso_inter).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_2,
                        style = "TCheckbutton").grid(row = 3, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso intra
        ttk.Label(self.top_frame, text = "Uiso_inter", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 12),
                     textvariable = self.Uiso_intra).grid(row = 3, column = 10, columnspan = 1, padx = 5, sticky = "sw")
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

        self.combobox = ttk.Combobox(self.top_frame, font = ("Arial", 12), values = ["Least Squares", "fmin",
                              "Basin Hopping"])
        self.combobox.current(0)
        self.combobox.grid(row = 7, column = 4, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.top_frame, text = "Run The FIT", command = self.run_the_fit,
                             style = "TButton").grid(row = 7, column = 6, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Reset To Default", command = self.reset_to_default,
                             style = "TButton").grid(row = 7, column = 9, columnspan = 3, padx = 5, sticky = "sw")

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
    
    def load_pdf_data(self):
        self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                         filetypes = [("Text Documents", "*.gr")])
    
    def reset_to_default(self):
        ##checkbox
        self.var_0.set(1)
        self.var_1.set(1)
        self.var_4.set(1)
        self.var_5.set(1)
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


    def run_the_fit(self):
        
        ####codes for PDF fits
        def makeRecipe(stru1, stru2, stru3, datname):
            ## The Profile
            profile = Profile()
            
            # Load data and add it to the profile
            parser = PDFParser()
            parser.parseFile(datname)
            profile.loadParsedData(parser)
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
            
            generator_MEF_Cryst_B.qdamp.value = float(self.Qdamp.get())
            generator_MEF_Mole_B.qdamp.value = float(self.Qdamp.get())
            generator_MEF_Intra.qdamp.value = float(self.Qdamp.get())
            
            # Vary the gloabal scale as well.
            recipe.addVar(contribution.scale, float(self.scale.get()))
            
            #############################################################################################
            ############### First the MEF_Cryst_B parameters ############################################
            #############################################################################################
            phase_MEF_Cryst_B = generator_MEF_Cryst_B.phase
            
            lat = phase_MEF_Cryst_B.getLattice()
            atoms = phase_MEF_Cryst_B.getScatterers()
            
            recipe.newVar("Uiso_Inter", float(self.Uiso_inter.get()), tag = "T1")
            recipe.newVar("lat_a", float(self.lat_a.get()), tag = "lat")
            recipe.newVar("lat_b", float(self.lat_b.get()), tag = "lat")
            recipe.newVar("lat_c", float(self.lat_c.get()), tag = "lat")
            recipe.newVar("alpha", float(self.alpha.get()), tag = "lat")
            recipe.newVar("beta",  float(self.beta.get()), tag = "lat")
            recipe.newVar("gamma", float(self.gamma.get()), tag = "lat")
            
            recipe.constrain(lat.a, "lat_a")
            recipe.constrain(lat.b, "lat_b")
            recipe.constrain(lat.c, "lat_c")
            recipe.constrain(lat.alpha, "alpha")
            recipe.constrain(lat.beta, "beta")
            recipe.constrain(lat.gamma, "gamma")
            
            for atom in atoms:
                if atom.element.title() == "N":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
                
                elif atom.element.title() == "O":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
                
                elif atom.element.title() == "C":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
                
                elif atom.element.title() == "H":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
            
            generator_MEF_Cryst_B.delta2.value = float(self.delta2_cryst.get())
            
            
            #############################################################################################
            ############### Second the MEF_Mole_B parameters ############################################
            #############################################################################################
            phase_MEF_Mole_B = generator_MEF_Mole_B.phase
            generator_MEF_Mole_B.setQmin(float(self.Qmin.get()))
            generator_MEF_Mole_B.setQmax(float(self.Qmax.get()))
            recipe.newVar("zoom_Mole_B", float(self.zoom_mol.get()), tag = "lat2")
            
            lat = phase_MEF_Mole_B.getLattice()
            recipe.constrain(lat.a, "zoom_Mole_B")
            recipe.constrain(lat.b, "zoom_Mole_B")
            recipe.constrain(lat.c, "zoom_Mole_B")
            # Constrain fractional xyz parameters
            atoms = phase_MEF_Mole_B.getScatterers()
            # Constrain ADPs
            
            for atom in atoms:
                if atom.element.title() == "C":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
                
                elif atom.element.title() == "O":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
                
                elif atom.element.title() == "N":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
                
                elif atom.element.title() == "H":
                    recipe.constrain(atom.Uiso, "Uiso_Inter")
            
            generator_MEF_Mole_B.delta2.value = float(self.delta2_mol.get())
            
            #############################################################################################
            ############### Third the intra molecule parameters##########################################
            #############################################################################################
            phase_MEF_Intra = generator_MEF_Intra.phase
            generator_MEF_Intra.setQmin(float(self.Qmin.get()))
            generator_MEF_Intra.setQmax(float(self.Qmax.get()))
            recipe.newVar("zoom_Intra", float(self.zoom_intra.get()), tag = "lat3")
            
            lat = phase_MEF_Intra.getLattice()
            recipe.constrain(lat.a, "zoom_Intra")
            recipe.constrain(lat.b, "zoom_Intra")
            recipe.constrain(lat.c, "zoom_Intra")
            # Constrain fractional xyz parameters
            atoms = phase_MEF_Intra.getScatterers()
            # Constrain ADPs
            recipe.newVar("Uiso_Intra", float(self.Uiso_intra.get()), tag = "T2")
            
            for atom in atoms:
                if atom.element.title() == "C":
                    recipe.constrain(atom.Uiso, "Uiso_Intra")
                
                elif atom.element.title() == "O":
                    recipe.constrain(atom.Uiso, "Uiso_Intra")
                
                elif atom.element.title() == "N":
                    recipe.constrain(atom.Uiso, "Uiso_Intra")
                
                elif atom.element.title() == "H":
                    recipe.constrain(atom.Uiso, "Uiso_Intra")
            
            generator_MEF_Intra.delta2.value = float(self.delta2_intra.get())
            
            # Give the recipe away so it can be used!
            return recipe
                        
                        
        def plotResults(recipe):
     
            
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


        stru1 = Structure(filename = self.input_cryst_struct)
        stru2 = Structure(filename = self.input_mole_struct)
        stru3 = Structure(filename = self.input_mole_struct)
        
        data  = self.input_pdf_gr
        
        print data

        recipe = makeRecipe(stru1, stru2, stru3, data)
        
        from scipy.optimize import leastsq
        leastsq(recipe.residual, recipe.values)
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
