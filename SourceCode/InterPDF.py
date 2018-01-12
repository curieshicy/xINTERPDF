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
import matplotlib.gridspec as gridspec
##DiffPy-CMI loading
from pyobjcryst.crystal import CreateCrystalFromCIF
from diffpy.srfit.pdf import DebyePDFGenerator, PDFParser, PDFGenerator
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.Structure import Structure
from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator

import numpy

class InterPDF():
    def __init__(self, master):
        self.master = master
        self.master.title("Extract Intermolecular PDF")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 300) # width + height
        self.master.resizable(False, False)
        
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Arial", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        
        ttk.Label(self.master, text = "The study of intermolecular contribution via subtracting PDF of a single molecule from the total PDF.", font = ("Arial", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()

        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures", justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 0, column = 3, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Molecule Structure", command = self.load_molecule,
                 style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Crystal Structure",command = self.load_crystal,
                 style = "TButton").grid(row = 1, column = 3, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Load PDF data",command = self.load_pdf_data,
                 style = "TButton").grid(row = 1, column = 6, columnspan = 3, padx = 5)

        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 2, column = 3, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        self.middle_frame = ttk.Frame(self.master)
        self.middle_frame.pack()
        
        ######adding tabs##############
        self.nb = ttk.Notebook(self.middle_frame)
        self.nb.pack(fill = BOTH, expand = "yes")
        
        self.tab_1 = ttk.Frame(self.middle_frame)
        self.tab_2 = ttk.Frame(self.middle_frame)
        
        self.nb.add(self.tab_1, text = "Molecule Parameters")
        self.nb.add(self.tab_2, text = "Crystal Parameters")
        

        ## these are the entries for tab 1
        self.qdamp_1 = StringVar()
        self.qbroad_1 = StringVar()
        self.qmin_1 = StringVar()
        self.qmax_1 = StringVar()
        self.rmin_1  = StringVar()
        self.rmax_1  = StringVar()
        self.rstep_1 = StringVar()
        self.delta2_mol = StringVar()
        self.mol_pc_dpc = StringVar()

        ## these are the entries for tab 2
        self.qdamp_2 = StringVar()
        self.qbroad_2 = StringVar()
        self.qmin_2 = StringVar()
        self.qmax_2 = StringVar()
        self.rmin_2  = StringVar()
        self.rmax_2  = StringVar()
        self.rstep_2 = StringVar()
        self.delta2_cryst = StringVar()
        self.cryst_pc_dpc = StringVar()

        ##set default values
        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin_1.set(0.0)
        self.rmax_1.set(20.0)
        self.rstep_1.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
        
        self.qdamp_2.set(0.04)
        self.qbroad_2.set(0.02)
        self.qmin_2.set(0.0)
        self.qmax_2.set(20.0)
        self.rmin_2.set(0.0)
        self.rmax_2.set(20.0)
        self.rstep_2.set(0.01)
        self.delta2_cryst.set(0.0)
        self.cryst_pc_dpc.set("PDFCalculator")

        #######################################Tab 1#########################################################
        #Qdamp
        ttk.Label(self.tab_1, text = "Qdamp", justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = self.qdamp_1).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_1, text = "Qbroad", justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8,font = ("Arial", 12),
                textvariable = self.qbroad_1).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_1, text = "Qmin", justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = self.qmin_1).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_1, text = "Qmax", justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = self.qmax_1).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                
        #rmin
        ttk.Label(self.tab_1, text = "Rmin", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                  textvariable = self.rmin_1).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #rmax
        ttk.Label(self.tab_1, text = "Rmax", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                  textvariable = self.rmax_1).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_1, text = "Rstep", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                  textvariable = self.rstep_1).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #delta2
        ttk.Label(self.tab_1, text = "Delta2", justify = CENTER,
              font = ("Arial", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = self.delta2_mol).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                            
        ttk.Label(self.tab_1, text = "Calculator", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")

        self.combobox_1 = ttk.Combobox(self.tab_1, font = ("Arial", 12),
                            textvariable = self.mol_pc_dpc, values = ["DebyePDFCalculator","PDFCalculator"])
        self.combobox_1.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.tab_1, text = "Reset To Default", command = self.reset_to_default_mol,
                     style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_1, text = "Expand Thermals", command = self.load_mol_uiso_elements,
                    style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                 
    #######################################Tab 2#########################################################
        #Qdamp
        ttk.Label(self.tab_2, text = "Qdamp", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.qdamp_2).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_2, text = "Qbroad", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8,font = ("Arial", 12),
        textvariable = self.qbroad_2).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_2, text = "Qmin", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.qmin_2).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_2, text = "Qmax", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.qmax_2).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rmin
        ttk.Label(self.tab_2, text = "Rmin", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.rmin_2).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #rmax
        ttk.Label(self.tab_2, text = "Rmax", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.rmax_2).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_2, text = "Rstep", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.rstep_2).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #delta2
        ttk.Label(self.tab_2, text = "Delta2", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
        textvariable = self.delta2_cryst).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Label(self.tab_2, text = "Calculator", justify = CENTER,
            font = ("Arial", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")
    
        self.combobox_2 = ttk.Combobox(self.tab_2, font = ("Arial", 12), textvariable = self.cryst_pc_dpc, values = ["PDFCalculator", "DebyePDFCalculator"])
        self.combobox_2.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Button(self.tab_2, text = "Reset To Default", command = self.reset_to_default_cryst,
           style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_2, text = "Expand Thermals", command = self.load_cryst_uiso_elements,
                      style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        self.bottom_frame = ttk.Frame(self.master)
        self.bottom_frame.pack()

        ttk.Label(self.bottom_frame, text = "Step 3: Run The Simulation",justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 0, column = 1, columnspan = 3, padx = 5, pady = 5)

        ttk.Button(self.bottom_frame, text = "Visualize", command = self.visualize,
                 style = "TButton").grid(row = 0, column = 4, columnspan = 3, padx = 5)


    def load_molecule(self):
        self.input_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                          filetypes = [("Text Documents", "*.xyz")])
                                                          
    def load_crystal(self):
        self.input_cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".cif",
                                                        filetypes = [("Text Documents", "*.cif")])
                                                          
    def load_pdf_data(self):
        self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                     filetypes = [("Text Documents", "*.gr")])

    def load_mol_uiso_elements(self):
    
        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
            
        ##the labels for Us, the thermals
        ttk.Label(self.tab_1, text = "U11", justify = CENTER,
                  font = ("Arial", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "U22", justify = CENTER,
                font = ("Arial", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "U33", justify = CENTER,
                font = ("Arial", 14, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "U12", justify = CENTER,
                font = ("Arial", 14, "bold")).grid(row = 2, column = 8, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "U13", justify = CENTER,
                font = ("Arial", 14, "bold")).grid(row = 2, column = 9, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "U23", justify = CENTER,
                font = ("Arial", 14, "bold")).grid(row = 2, column = 10, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "Occ", justify = CENTER,
                          font = ("Arial", 14, "bold")).grid(row = 2, column = 11, columnspan = 1, padx = 5)
                

        for i, j in enumerate(self.mole_elements_list):
            
            setattr(self, "mol_{}_U11".format(j), StringVar())
            setattr(self, "mol_{}_U22".format(j), StringVar())
            setattr(self, "mol_{}_U33".format(j), StringVar())
            setattr(self, "mol_{}_U12".format(j), StringVar())
            setattr(self, "mol_{}_U13".format(j), StringVar())
            setattr(self, "mol_{}_U23".format(j), StringVar())
            setattr(self, "mol_{}_Occ".format(j), StringVar())
            getattr(self, "mol_{}_U11".format(j)).set(0.005)
            getattr(self, "mol_{}_U22".format(j)).set(0.005)
            getattr(self, "mol_{}_U33".format(j)).set(0.005)
            getattr(self, "mol_{}_U12".format(j)).set(0.000)
            getattr(self, "mol_{}_U13".format(j)).set(0.000)
            getattr(self, "mol_{}_U23".format(j)).set(0.000)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)

            ttk.Label(self.tab_1, text = self.mole_elements_list[i], justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")

            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_U11".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_U22".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_U33".format(j))).grid(row = 3+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_U12".format(j))).grid(row = 3+i, column = 8, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_U13".format(j))).grid(row = 3+i, column = 9, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_U23".format(j))).grid(row = 3+i, column = 10, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "mol_{}_Occ".format(j))).grid(row = 3+i, column = 11, columnspan = 1, padx = 5, sticky = "sw")

    def load_cryst_uiso_elements(self):
        
        self.cryst_elements = set(Structure(filename = self.input_cryst_struct).element) #return a dictionary
            
            #put elements into a list
        self.cryst_elements_list = []
        for i in self.cryst_elements:
            self.cryst_elements_list.append(i)
            ##the labels for Us
        ttk.Label(self.tab_2, text = "U11", justify = CENTER,
                      font = ("Arial", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "U22", justify = CENTER,
                      font = ("Arial", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "U33", justify = CENTER,
                      font = ("Arial", 14, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "U12", justify = CENTER,
                      font = ("Arial", 14, "bold")).grid(row = 2, column = 8, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "U13", justify = CENTER,
                      font = ("Arial", 14, "bold")).grid(row = 2, column = 9, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "U23", justify = CENTER,
                      font = ("Arial", 14, "bold")).grid(row = 2, column = 10, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "Occ", justify = CENTER,
                    font = ("Arial", 14, "bold")).grid(row = 2, column = 11, columnspan = 1, padx = 5)

        for i, j in enumerate(self.cryst_elements):

            setattr(self, "cryst_{}_U11".format(j), StringVar())
            setattr(self, "cryst_{}_U22".format(j), StringVar())
            setattr(self, "cryst_{}_U33".format(j), StringVar())
            setattr(self, "cryst_{}_U12".format(j), StringVar())
            setattr(self, "cryst_{}_U13".format(j), StringVar())
            setattr(self, "cryst_{}_U23".format(j), StringVar())
            setattr(self, "cryst_{}_Occ".format(j), StringVar())
            getattr(self, "cryst_{}_U11".format(j)).set(0.005)
            getattr(self, "cryst_{}_U22".format(j)).set(0.005)
            getattr(self, "cryst_{}_U33".format(j)).set(0.005)
            getattr(self, "cryst_{}_U12".format(j)).set(0.000)
            getattr(self, "cryst_{}_U13".format(j)).set(0.000)
            getattr(self, "cryst_{}_U23".format(j)).set(0.000)
            getattr(self, "cryst_{}_Occ".format(j)).set(1.000)
    
            ttk.Label(self.tab_2, text = self.cryst_elements_list[i], justify = CENTER,
                font = ("Arial", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_U11".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_U22".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_U33".format(j))).grid(row = 3+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_U12".format(j))).grid(row = 3+i, column = 8, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_U13".format(j))).grid(row = 3+i, column = 9, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_U23".format(j))).grid(row = 3+i, column = 10, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Arial", 12),
                textvariable = getattr(self, "cryst_{}_Occ".format(j))).grid(row = 3+i, column = 11, columnspan = 1, padx = 5, sticky = "sw")

    def reset_to_default_mol(self):

        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin_1.set(0.0)
        self.rmax_1.set(20.0)
        self.rstep_1.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
    
    def reset_to_default_cryst(self):
        
        self.qdamp_2.set(0.04)
        self.qbroad_2.set(0.02)
        self.qmin_2.set(0.0)
        self.qmax_2.set(20.0)
        self.rmin_2.set(0.0)
        self.rmax_2.set(20.0)
        self.rstep_2.set(0.01)
        self.delta2_cryst.set(0.0)
        self.cryst_pc_dpc.set("PDFCalculator")
    
    def visualize(self):  ###create a two-panel figure, no curves ###

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.top = tk.Toplevel()
        self.top.title("Visualize the plots")
      
        self.top_frame = ttk.Frame(self.top, padding = (10, 10))
        self.top_frame.pack()
        self.fig = plt.figure(figsize=(10, 6), dpi=100) ##create a figure; modify the size here
        
        self.fig.add_subplot(211)
        
        plt.title("Intermolecular PDFs")
        plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        
        self.fig.add_subplot(212)


        plt.title("Difference PDFs")
        plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        
        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.top_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.top_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.bottom_frame = ttk.Frame(self.top, padding = (10, 10))
        self.bottom_frame.pack()

        ttk.Button(self.bottom_frame, text = "Make The Plot", command = self.make_the_plot,
           style = "TButton").grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.bottom_frame, text = "Fine Tune Parameters", command = self.fine_tune_parameters,
                      style = "TButton").grid(row = 0, column = 2, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.bottom_frame, text = "Clear The Plot", command = self.clear_the_plot,
               style = "TButton").grid(row = 0, column = 4, columnspan = 2, padx = 5, sticky = "sw")
    
    def plot_with_cfgs(self, mol_cfg_given, cryst_cfg_given):
    
    ###this is to generate the plots when the configuration parameters are given for molecule and crystals###
    ###There are four combinations of PC, DPC for molecules and crystals###

        if self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            r1, g1 = mol_pc(self.mol_struc)
            cryst_pc = PDFCalculator(**cryst_cfg_given)
            r2, g2 = cryst_pc(self.cryst_struc)
            
            self.fig.add_subplot(211)
            plt.title("Intermolecular PDFs")
            #plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r1, g1, "r-", lw=2, label = "mol_pc")
            plt.plot(r2, g2, "b-", lw=2, label = "cryst_pc")
            plt.legend(loc=0)
            
            self.fig.add_subplot(212)
            plt.title("Difference PDFs")
            plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r1, g2 - g1, "g-", lw=2, label = "cryst_pc - mol_pc")
            plt.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()

        elif self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            r3, g3 = mol_pc(self.mol_struc)
            cryst_pc = DebyePDFCalculator(**cryst_cfg_given)
            r4, g4 = cryst_pc(self.cryst_struc)
            
            self.fig.add_subplot(211)
            plt.title("Intermolecular PDFs")
            #plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r3, g3, "r-", lw=2, label = "mol_pc")
            plt.plot(r4, g4, "b-", lw=2, label = "cryst_dpc")
            plt.legend(loc=0)
            
            self.fig.add_subplot(212)
            plt.title("Difference PDFs")
            plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r3, g4 - g3, "g-", lw=2, label = "cryst_dpc - mol_pc")
            plt.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()

        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            r5, g5 = mol_dpc(self.mol_struc)
            cryst_pc = PDFCalculator(**cryst_cfg_given)
            r6, g6 = cryst_pc(self.cryst_struc)
            
            self.fig.add_subplot(211)
            plt.title("Intermolecular PDFs")
            #plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r5, g5, "r-", lw=2, label = "mol_dpc")
            plt.plot(r6, g6, "b-", lw=2, label = "cryst_pc")
            plt.legend(loc=0)
            
            self.fig.add_subplot(212)
            plt.title("Difference PDFs")
            plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r5, g6 - g5, "g-", lw=2, label = "cryst_pc - mol_dpc")
            plt.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()

        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "DeybePDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            r7, g7 = mol_dpc(self.mol_struc)
            cryst_dpc = DebyePDFCalculator(**cryst_cfg_given)
            r8, g8 = cryst_dpc(self.cryst_struc)
            
            self.fig.add_subplot(211)
            plt.title("Intermolecular PDFs")
            #plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r7, g7, "r-", lw=2, label = "mol_dpc")
            plt.plot(r8, g8, "b-", lw=2, label = "cryst_dpc")
            plt.legend(loc=0)
            
            self.fig.add_subplot(212)
            plt.title("Difference PDFs")
            plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            plt.xticks(fontsize = 11)
            plt.yticks(fontsize = 11)
            plt.plot(r7, g8 - g7, "g-", lw=2, label = "cryst_dpc - mol_dpc")
            plt.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()
    
    def make_the_plot(self):
    ###first the molecular PDF ###
        self.mol_cfg ={
        "qmax": float(self.qmax_1.get()),
        "qmin": float(self.qmin_1.get()),
        "rmin": float(self.rmin_1.get()),
        "rmax": float(self.rmax_1.get()),
        "qdamp": float(self.qdamp_1.get()),
        "qbroad": float(self.qbroad_1.get()),
        "delta2": float(self.delta2_mol.get()),
        "rstep": float(self.rstep_1.get())
                       }

        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        self.mol_struc = Structure(filename =  self.input_mole_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
        
        #print self.mole_elements_list

        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].U11 = float(getattr(self, "mol_{}_U11".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].U22 = float(getattr(self, "mol_{}_U22".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].U33 = float(getattr(self, "mol_{}_U33".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].U12 = float(getattr(self, "mol_{}_U12".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].U13 = float(getattr(self, "mol_{}_U13".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].U23 = float(getattr(self, "mol_{}_U23".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].occupancy = float(getattr(self, "mol_{}_Occ".format(j)).get())


        ###second the crystalline PDF ###
        self.cryst_cfg ={
        "qmax": float(self.qmax_2.get()),
        "qmin": float(self.qmin_2.get()),
        "rmin": float(self.rmin_2.get()),
        "rmax": float(self.rmax_2.get()),
        "qdamp": float(self.qdamp_2.get()),
        "qbroad": float(self.qbroad_2.get()),
        "delta2": float(self.delta2_cryst.get()),
        "rstep": float(self.rstep_2.get())
        }
        
        self.cryst_elements = set(Structure(filename = self.input_cryst_struct).element) #return a dictionary
        self.cryst_struc = Structure(filename =  self.input_cryst_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.cryst_elements_list = []
        for i in self.cryst_elements:
            self.cryst_elements_list.append(i)
        
        #print self.mole_elements_list
        
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].U11 = float(getattr(self, "cryst_{}_U11".format(j)).get())
            self.cryst_struc[self.cryst_struc.element == j].U22 = float(getattr(self, "cryst_{}_U22".format(j)).get())
            self.cryst_struc[self.cryst_struc.element == j].U33 = float(getattr(self, "cryst_{}_U33".format(j)).get())
            self.cryst_struc[self.cryst_struc.element == j].U12 = float(getattr(self, "cryst_{}_U12".format(j)).get())
            self.cryst_struc[self.cryst_struc.element == j].U13 = float(getattr(self, "cryst_{}_U13".format(j)).get())
            self.cryst_struc[self.cryst_struc.element == j].U23 = float(getattr(self, "cryst_{}_U23".format(j)).get())
            self.cryst_struc[self.cryst_struc.element == j].occupancy = float(getattr(self, "cryst_{}_Occ".format(j)).get())

        self.plot_with_cfgs(self.mol_cfg, self.cryst_cfg)



    def fine_tune_parameters(self):
        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.panel = tk.Toplevel()
        self.panel.title("Fine Tune Parameters")
        
        self.panel_frame = ttk.Frame(self.panel, padding = (10, 10))
        self.panel_frame.pack()
        ######adding tabs##############
        self.panel_tab = ttk.Notebook(self.panel_frame)
        self.panel_tab.pack(fill = BOTH, expand = "yes")
        
        self.panel_tab_1 = ttk.Frame(self.panel_frame)
        self.panel_tab_2 = ttk.Frame(self.panel_frame)
        
        self.panel_tab.add(self.panel_tab_1, text = "Molecule Parameters")
        self.panel_tab.add(self.panel_tab_2, text = "Crystal Parameters")
    
    #########These are parameters for molecules###############
        ttk.Label(self.panel_tab_1, text = "Qdamp", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qbroad", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmin", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmax", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmin", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmax", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rstep", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Delta2", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")


        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        
        self.qdamp_1_entry = StringVar()
        self.qbroad_1_entry = StringVar()
        self.qmin_1_entry = StringVar()
        self.qmax_1_entry = StringVar()
        self.rmin_entry  = StringVar()
        self.rmax_entry  = StringVar()
        self.rstep_entry = StringVar()
        self.delta2_mol_entry = StringVar()

        self.m_qdamp_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 0.5, orient = HORIZONTAL, variable = self.qdamp_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qdamp_scale.grid(row = 0, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_qbroad_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 0.5, orient = HORIZONTAL, variable = self.qbroad_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qbroad_scale.grid(row = 1, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_qmin_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 2.0, orient = HORIZONTAL, variable = self.qmin_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qmin_scale.grid(row = 2, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_qmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 30.0, orient = HORIZONTAL, variable = self.qmax_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmin_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.m_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.m_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rstep_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.m_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_delta2_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_mol_entry, command = self.update_plot_fine_tune_para)
        self.m_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")
    
        self.qdamp_1_entry.set(0.04)
        self.qbroad_1_entry.set(0.02)
        self.qmin_1_entry.set(0.0)
        self.qmax_1_entry.set(20.0)
        self.rmin_entry.set(0.0)
        self.rmax_entry.set(20.0)
        self.rstep_entry.set(0.01)
        self.delta2_mol_entry.set(0.0)

        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "2.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "30.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "100.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "200.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "10.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        
        em1 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.qdamp_1_entry)
        em1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em1.bind("<Return>", self.update_plot_fine_tune_para)
        em2 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.qbroad_1_entry)
        em2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em2.bind("<Return>", self.update_plot_fine_tune_para)
        em3 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.qmin_1_entry)
        em3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em3.bind("<Return>", self.update_plot_fine_tune_para)
        em4 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.qmax_1_entry)
        em4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em4.bind("<Return>", self.update_plot_fine_tune_para)
        em5 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.rmin_entry)
        em5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em5.bind("<Return>", self.update_plot_fine_tune_para)
        em6 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.rmax_entry)
        em6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em6.bind("<Return>", self.update_plot_fine_tune_para)
        em7 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.rstep_entry)
        em7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em7.bind("<Return>", self.update_plot_fine_tune_para)
        em8 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Arial", 12),textvariable = self.delta2_mol_entry)
        em8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em8.bind("<Return>", self.update_plot_fine_tune_para)
        
        ttk.Button(self.panel_tab_1, text = "Reset All", command = self.Reset_m_c_Entry,
                   style = "TButton").grid(row = 8, column = 0, columnspan = 2, padx = 5)

        #########These are parameters for crystals###############
        ttk.Label(self.panel_tab_2, text = "Qdamp", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qbroad", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qmin", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qmax", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rmin", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rmax", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rstep", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Delta2", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")

        self.qdamp_2_entry = StringVar()
        self.qbroad_2_entry = StringVar()
        self.qmin_2_entry = StringVar()
        self.qmax_2_entry = StringVar()
#        self.rmin_2_entry  = StringVar()
#        self.rmax_2_entry  = StringVar()
#        self.rstep_2_entry = StringVar()
        self.delta2_cryst_entry = StringVar()

        self.c_qdamp_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 0.5, orient = HORIZONTAL, variable = self.qdamp_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qdamp_scale.grid(row = 0, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_qbroad_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 0.5, orient = HORIZONTAL, variable = self.qbroad_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qbroad_scale.grid(row = 1, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_qmin_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 2.0, orient = HORIZONTAL, variable = self.qmin_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qmin_scale.grid(row = 2, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_qmax_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 30.0, orient = HORIZONTAL, variable = self.qmax_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_rmin_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.c_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        self.c_rmax_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.c_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_rstep_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.c_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_delta2_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_cryst_entry, command = self.update_plot_fine_tune_para)
        self.c_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.qdamp_2_entry.set(0.04)
        self.qbroad_2_entry.set(0.02)
        self.qmin_2_entry.set(0.0)
        self.qmax_2_entry.set(20.0)
#        self.rmin_2_entry.set(0.0)
#        self.rmax_2_entry.set(20.0)
#        self.rstep_2_entry.set(0.01)
        self.delta2_cryst_entry.set(0.0)

        ttk.Label(self.panel_tab_2, text = "0.5", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.5", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "2.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "30.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "100.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "200.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "1.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "10.0", justify = CENTER,
        font = ("Arial", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        

        ec1 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.qdamp_2_entry)
        ec1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec1.bind("<Return>", self.update_plot_fine_tune_para)
        ec2 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.qbroad_2_entry)
        ec2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec2.bind("<Return>", self.update_plot_fine_tune_para)
        ec3 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.qmin_2_entry)
        ec3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec3.bind("<Return>", self.update_plot_fine_tune_para)
        ec4 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.qmax_2_entry)
        ec4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec4.bind("<Return>", self.update_plot_fine_tune_para)
        ec5 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.rmin_entry)
        ec5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec5.bind("<Return>", self.update_plot_fine_tune_para)
        ec6 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.rmax_entry)
        ec6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec6.bind("<Return>", self.update_plot_fine_tune_para)
        ec7 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.rstep_entry)
        ec7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec7.bind("<Return>", self.update_plot_fine_tune_para)
        ec8 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Arial", 12),textvariable = self.delta2_cryst_entry)
        ec8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec8.bind("<Return>", self.update_plot_fine_tune_para)

        ttk.Button(self.panel_tab_2, text = "Reset All", command = self.Reset_m_c_Entry,
           style = "TButton").grid(row = 8, column = 0, columnspan = 2, padx = 5)
    
    
    def update_plot_fine_tune_para(self, event):
    
        self.fig.clf()
        ###first the molecular PDF ###
        self.mol_fine_tune_cfg ={
        "qmax": float(self.qmax_1_entry.get()),
        "qmin": float(self.qmin_1_entry.get()),
        "rmin": float(self.rmin_entry.get()),
        "rmax": float(self.rmax_entry.get()),
        "qdamp": float(self.qdamp_1_entry.get()),
        "qbroad": float(self.qbroad_1_entry.get()),
        "delta2": float(self.delta2_mol_entry.get()),
        "rstep": float(self.rstep_entry.get())
        }

        self.cryst_fine_tune_cfg ={
        "qmax": float(self.qmax_2_entry.get()),
        "qmin": float(self.qmin_2_entry.get()),
        "rmin": float(self.rmin_entry.get()),
        "rmax": float(self.rmax_entry.get()),
        "qdamp": float(self.qdamp_2_entry.get()),
        "qbroad": float(self.qbroad_2_entry.get()),
        "delta2": float(self.delta2_cryst_entry.get()),
        "rstep": float(self.rstep_entry.get())
        }
        
        self.plot_with_cfgs(self.mol_fine_tune_cfg, self.cryst_fine_tune_cfg)

    def Reset_m_c_Entry(self):

        self.fig.clf()

        self.qdamp_1_entry.set(0.04)
        self.qbroad_1_entry.set(0.02)
        self.qmin_1_entry.set(0.0)
        self.qmax_1_entry.set(20.0)
        self.rmin_entry.set(0.0)
        self.rmax_entry.set(20.0)
        self.rstep_entry.set(0.01)
        self.delta2_mol_entry.set(0.0)
        
        self.qdamp_2_entry.set(0.04)
        self.qbroad_2_entry.set(0.02)
        self.qmin_2_entry.set(0.0)
        self.qmax_2_entry.set(20.0)
        self.delta2_cryst_entry.set(0.0)
        
        
        self.mol_reset_cfg ={
        "qmax": 20.0,
        "qmin": 0.0,
        "rmin": 0.0,
        "rmax": 20.0,
        "qdamp": 0.04,
        "qbroad": 0.02,
        "delta2": 0.0,
        "rstep": 0.01
        }
        
        self.cryst_reset_cfg ={
        "qmax": 20.0,
        "qmin": 0.0,
        "rmin": 0.0,
        "rmax": 20.0,
        "qdamp": 0.04,
        "qbroad": 0.02,
        "delta2": 0.0,
        "rstep": 0.01
        }
        
        
        self.plot_with_cfgs(self.mol_reset_cfg, self.cryst_reset_cfg)


    def clear_the_plot(self):
        self.fig.clf()

        self.fig.add_subplot(211)
        
        plt.title("Intermolecular PDFs")
        #plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        
        self.fig.add_subplot(212)
        plt.title("Difference PDFs")
        plt.xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        plt.ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        
        self.fig.tight_layout()
        self.canvas.draw()

def main():
    root = Tk()
    GUI = InterPDF(root)
    root.mainloop()

if __name__ == "__main__": main()