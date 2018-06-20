""""
This is the code for xINTERPDF program. It is used to extract structural information from a measured X-ray pair distribution function for organic compounds.
The code is written by Chenyang Shi at AbbVie.
Send him an email at cs3000@columbia.edu for questions and suggestions.

"""
#loading Tkinter modules
from Tkinter import *
import Tkinter as tk
import ttk
import tkFileDialog
import tkMessageBox

import numpy
from math import exp
import webbrowser
import os  # for loading files or exporting files

##loading matplotlib modules
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
##loading DiffPy-CMI modules
from pyobjcryst.crystal import CreateCrystalFromCIF
from diffpy.srfit.pdf import DebyePDFGenerator, PDFParser, PDFGenerator, PDFContribution
from diffpy.srfit.fitbase import Profile, Calculator
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.Structure import Structure, loadStructure
from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator
from pyobjcryst.molecule import Molecule
from scipy.optimize import minimize
from scipy.stats import pearsonr
from sklearn.decomposition import PCA as SKPCA
import itertools


class Overall_Look:
    def __init__(self, master):
        self.master = master
        self.master.title("xINTERPDF")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(600, 250) # width + height
        self.master.resizable(False, False)
        
        ## add menu bar
        self.menubar = Menu(self.master)
        
        ##Dropdown menu 2
        self.interpdf_menu = Menu(self.menubar, tearoff = 0)
        self.interpdf_menu.add_command(label = "Crystalline", font = ("Times", 16), command = self.c_inter_pdf)
        self.interpdf_menu.add_command(label = "Amorphous", font = ("Times", 16), command = self.a_inter_pdf)
        self.interpdf_menu.add_command(label = "ASD", font = ("Times", 16), command = self.asd_inter_pdf)
        self.menubar.add_cascade(label = "Inter PDF", font = ("Times", 16), menu = self.interpdf_menu)

        ##Dropdown menu 3
        self.pdffit_menu = Menu(self.menubar, tearoff = 0)
        self.pdffit_menu.add_command(label = "PDF Fit", font = ("Times", 16), command = self.pdf_model_uiso)
        self.pdffit_menu.add_command(label = "Breakdown of Fit", font = ("Times", 16), command = self.breakdown_fit)
        self.menubar.add_cascade(label = "PDF Fit", font = ("Times", 16), menu = self.pdffit_menu)
        
        ##Dropdown menu 4
        self.utility_menu = Menu(self.menubar, tearoff = 0)
        self.utility_menu.add_command(label = "Linear Combination", font = ("Times", 16), command = self.linear_comb)
        self.utility_menu.add_command(label = "PCA", font = ("Times", 16), command = self.pca)
        self.menubar.add_cascade(label = "Utilities", font = ("Times", 16), menu = self.utility_menu)

        ##Dropdown menu 5
        self.help_menu = Menu(self.menubar, tearoff = 0)
        self.homepage = self.help_menu.add_command(label = "Homepage", font = ("Times", 16), command = self.callback)
        self.menubar.add_cascade(label = "Help", font = ("Times", 16), menu = self.help_menu)
        self.master.config(menu = self.menubar)

        ##define style in ttk##
        self.style = ttk.Style()
        self.style.theme_use('alt')
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91', font = ('Times', 16))
        self.style.configure('Header.TLabel', font = ('Times', 23, 'bold'))
        
        ##load image##
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        self.logo = PhotoImage(file = os.path.join(__location__, "Logo.gif")).subsample(3,3)

        ##this is the top frame, we will write a few sentences about this program##
        self.frame_header = ttk.Frame(master, relief = RIDGE, padding = (30, 15))
        self.frame_header.pack()

        ttk.Label(self.frame_header, image = self.logo).grid(row = 0, column = 0, rowspan = 2, columnspan = 2, sticky = "w")
        ttk.Label(self.frame_header, text = 'Welcome to use xINTERPDF program!',
                  style = 'Header.TLabel').grid(row = 0, column = 2, columnspan = 4)
        ttk.Label(self.frame_header, wraplength = 600,
                  text = ("This program uses DiffPy-CMI as a backend engine to simulate pair distribution functions (PDFs). It is designed to extract structural information from the measured X-ray PDF data for organic materials. Currently it supports (1) The study of intermolecular interaction (e.g. hydrogen bonds) by subtracting out the scattering signal of single molecule(s) in real space. (2) The PDF model fit of the crystalline organic compound using the method proposed by Prill et al. (J. Appl. Cryst., 2015, 48, 171-178.). (3) The phase quantification of physical mixtures of organics. (4) Generate Score/Scree plots based on Principle Component Analysis. Homepage: http://www.diffpy.org/products/xinterpdf.html."
                          )).grid(row = 1, column = 2, columnspan = 4)

        self.bottom_frame = ttk.Frame(master, relief = FLAT, padding = (30, 15)) 
        self.bottom_frame.pack()

    def pdf_model_uiso(self):
        self.pdf_model_uiso = tk.Toplevel(self.master)
        self.GUI = PDF_model_fit_Uiso(self.pdf_model_uiso)
    
    def breakdown_fit(self):
        self.breakdown_fit = tk.Toplevel(self.master)
        self.GUI = Break_down_fit(self.breakdown_fit)

    def c_inter_pdf(self):
        self.c_inter_pdf = tk.Toplevel(self.master)
        self.GUI = C_InterPDF(self.c_inter_pdf)
    
    def a_inter_pdf(self):
        self.a_inter_pdf = tk.Toplevel(self.master)
        self.GUI = A_InterPDF(self.a_inter_pdf)
    
    def asd_inter_pdf(self):
        self.asd_inter_pdf = tk.Toplevel(self.master)
        self.GUI = ASD_InterPDF(self.asd_inter_pdf)
    
    def linear_comb(self):
        self.linear_comb = tk.Toplevel(self.master)
        self.GUI = Linear_Comb(self.linear_comb)
    
    def pca(self):
        self.pca = tk.Toplevel(self.master)
        self.GUI = PCA(self.pca)

    def callback(self):
        webbrowser.open_new(r"http://www.diffpy.org/products/xinterpdf.html")

class Linear_Comb():
    def __init__(self, master):
        self.master = master
        self.master.title("Linear Combination of two (three) PDFs to match the third (fourth) one.")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 600) # width + height
        self.master.resizable(False, False)
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        
        ttk.Label(self.master, text = "Find scale factors that minimize difference/maximize similarity between s0*[s1*PDF-A + s2*PDF-B + (1-s1-s2)*PDF-C] and PDF-D", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()
                  
        ##here are the layout for step 1, load structure files
                  
        ttk.Label(self.top_frame, text = "Step 1: Load PDF Data", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
                  
        ttk.Button(self.top_frame, text = "Load PDF Data 1", command = self.load_pdf_1,
                             style = "TButton").grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load PDF Data 2",command = self.load_pdf_2,
                             style = "TButton").grid(row = 1, column = 2, columnspan = 2, padx = 5)
        ttk.Button(self.top_frame, text = "Load PDF Data 3 (Optional)",command = self.load_pdf_3,
                             style = "TButton").grid(row = 1, column = 4, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Load target PDF",command = self.load_pdf_target,
                             style = "TButton").grid(row = 1, column = 7, columnspan = 2, padx = 5)
                             
        self.optimizers = StringVar()
        self.rmin = StringVar()
        self.rmax = StringVar()
        self.rstep = StringVar()
        self.rmin.set(0)
        self.rmax.set(15)
        self.rstep.set(0.01)

        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Label(self.top_frame, text = "Rmin", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                  textvariable = self.rmin).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                  
        ttk.Label(self.top_frame, text = "Rmax", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                textvariable = self.rmax).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                
        ttk.Label(self.top_frame, text = "Rstep", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 3, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                  textvariable = self.rstep).grid(row = 3, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        
        ttk.Label(self.top_frame, text = "Select Optimizers", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 3, column = 6, columnspan = 1, padx = 5, sticky = "sw")

        self.combobox = ttk.Combobox(self.top_frame, font = ("Times", 12), textvariable = self.optimizers, values = ["Least Squares", "Pearson",
                                                                                                               "Least Absolute Deviation"])
        self.combobox.current(0)
        self.combobox.grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        
        
        ttk.Label(self.top_frame, text = "Step 3: Run The Fit", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Run The FIT", command = self.run_the_fit, style = "TButton").grid(row = 5, column = 0, columnspan = 2, padx = 5)
        
        ttk.Button(self.top_frame, text = "Save The FIT", command = self.save_the_fit, style = "TButton").grid(row = 5, column = 2, columnspan = 2, padx = 5)

        ttk.Button(self.top_frame, text = "Start Another Fit", command = self.start_new_fit, style = "TButton").grid(row = 5, column = 4, columnspan = 2, padx = 5)
        

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.bottom_frame = ttk.Frame(self.master, padding = (10, 10))
        self.bottom_frame.pack()

        self.fig = plt.figure(figsize=(9, 6), dpi=100) ##create a figure; modify the size here
        
        self.ax1 = self.fig.add_subplot(211)
        self.ax1.set_title("Plots of individual PDFs")
        self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
#        self.ax1.set_xticks(fontsize = 11)
#        self.ax1.set_yticks(fontsize = 11)

        self.ax2 = self.fig.add_subplot(212)
        self.ax2.set_title("Optimization of PDFs to match the target one")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.bottom_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.bottom_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)


    def load_pdf_1(self):
        try:
            self.input_pdf_gr_1 = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                             filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr_1)
                    
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def load_pdf_2(self):
        try:
            self.input_pdf_gr_2 = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                             filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr_2)
                    
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def load_pdf_3(self):
        try:
            self.input_pdf_gr_3 = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                           filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr_3)
                
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def load_pdf_target(self):
        try:
            self.input_pdf_gr_t = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                     filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr_t)
            
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def start_new_fit(self):
        self.rmin.set(0)
        self.rmax.set(15)
        self.rstep.set(0.01)
        self.combobox.current(0)
        
        self.fig.clf()
        self.ax1 = self.fig.add_subplot(211)
        self.ax1.set_title("Plots of individual PDFs")
        self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.set_title("Optimization of PDFs to match the target one")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()
        self.canvas.draw()

#        self.input_pdf_gr_1 = None
#        self.input_pdf_gr_2 = None
#        self.input_pdf_gr_3 = None
#        self.input_pdf_gr_t = None

    def run_the_fit(self):
    
        #make sure the PDF files are loaded
        self.pdf_1_ = 1 # 1 is success; 0 is failure
        self.pdf_2_ = 1
        self.pdf_3_ = 1
        self.pdf_t_ = 1
        self.global_scale_ = 1.0
        self.scale_1_ = 0.5
        
        ##these scale factors for matching three PDFs against one
        self._scale_1 = 0.33
        self._scale_2 = 0.33
        
        try:
            self.input_pdf_gr_1
        except:
            self.pdf_1_ = 0
        
        try:
            self.input_pdf_gr_2
        except:
            self.pdf_2_ = 0

        try:
            self.input_pdf_gr_3
        except:
            self.pdf_3_ = 0
        
        try:
            self.input_pdf_gr_t
        except:
            self.pdf_t_ = 0
    
        def return_x_y_array(input_gr, rmin, rmax, rstep):
            profile = Profile()
            parser = PDFParser()
            parser.parseFile(input_gr)
            profile.loadParsedData(parser)
            profile.setCalculationRange(xmin = rmin, xmax = rmax, dx = rstep)

            r = profile.x
            g = profile.y
            
            return (r, g)
        
        #Scenario 1: Two PDFs against a target PDF
        if self.pdf_1_ == 1 and self.pdf_2_ == 1 and self.pdf_3_ == 0 and self.pdf_t_ == 1:

            self.r1, self.g1 = return_x_y_array(self.input_pdf_gr_1, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            self.r2, self.g2 = return_x_y_array(self.input_pdf_gr_2, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            self.rt, self.gt = return_x_y_array(self.input_pdf_gr_t, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            
            def object_func_LS(para_list):
                scale = para_list[0]
                GSF = para_list[1]
                g_1_2 = (self.g1*scale + (1 - scale)*self.g2)*GSF
                target = sum((self.gt - g_1_2)**2)
                return target

            def object_func_LAD(para_list):
                scale = para_list[0]
                GSF = para_list[1]
                g_1_2 = (self.g1*scale + (1 - scale)*self.g2)*GSF
                target = sum(numpy.abs(g_1_2 - self.gt))
                return target

            def object_Pearson(para_list):
                scale = para_list[0]
                GSF = para_list[1]
                g_1_2 = (self.g1*scale + (1 - scale)*self.g2)*GSF
                target = pearsonr(g_1_2, self.gt)[0]
                if target > 0:
                    return (-1)*target
                else:
                    return target
        
            self.x0_ = [self.scale_1_, self.global_scale_]
        
            if self.optimizers.get() == "Least Squares":
                self.x0_ = minimize(object_func_LS, self.x0_, method = "nelder-mead").x
            
            elif self.optimizers.get() == "Pearson":
                self.x0_ = minimize(object_Pearson, self.x0_, method = "nelder-mead").x

            elif self.optimizers.get() == "Least Absolute Deviation":
                self.x0_ = minimize(object_func_LAD, self.x0_, method = "nelder-mead").x


            self.sum_fit = (self.g1*self.x0_[0] + self.g2*(1 - self.x0_[0]))*self.x0_[1]
            diffzero = -0.8*max(self.gt)*numpy.ones_like(self.gt)
            diff = self.gt - self.sum_fit + diffzero


            self.fig.clf()
            self.ax1 = self.fig.add_subplot(211)
        
            self.ax1.set_title("Plots of individual PDFs")
            self.ax1.plot(self.r1, self.g1, "k-", lw=2, label = "PDF 1")
            self.ax1.plot(self.r2, self.g2, "m-", lw=2, label = "PDF 2")
            self.ax1.plot(self.rt, self.gt, "b-", lw=2, label = "PDF target")

            self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax1.xaxis.set_tick_params(labelsize=11)
            self.ax1.yaxis.set_tick_params(labelsize=11)
            
            self.ax1.legend(loc = 0)

            self.ax2 = self.fig.add_subplot(212)
            self.ax2.set_title("Optimization of PDFs to match the target one")
            self.ax2.plot(self.r2, self.sum_fit, "r-", lw=2, label = "PDF fit, scale = {0:.4f}, global scale = {1:.4f}".format(self.x0_[0], self.x0_[1]))
            self.ax2.plot(self.rt, self.gt, "b-", lw=2, label = "PDF target")
            self.ax2.plot(self.r2, diff, "g-", lw=2, label = "PDF difference")
            self.ax2.plot(self.r2, diffzero, "k--", lw=1)

            self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax2.xaxis.set_tick_params(labelsize=11)
            self.ax2.yaxis.set_tick_params(labelsize=11)
            
            self.ax2.legend(loc = 0)
    
            self.fig.tight_layout()
            self.canvas.draw()

        #Scenario 2: Three PDFs against a target PDF
        elif self.pdf_1_ == 1 and self.pdf_2_ == 1 and self.pdf_3_ == 1 and self.pdf_t_ == 1:
            
            self.r1, self.g1 = return_x_y_array(self.input_pdf_gr_1, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            self.r2, self.g2 = return_x_y_array(self.input_pdf_gr_2, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            self.r3, self.g3 = return_x_y_array(self.input_pdf_gr_3, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            self.rt, self.gt = return_x_y_array(self.input_pdf_gr_t, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            
            def object_func_LS(para_list):
                scale1 = para_list[0]
                scale2 = para_list[1]
                scale3 = 1 - scale1 - scale2
                GSF = para_list[2]
                g_1_2_3 = (self.g1*scale1 + scale2*self.g2 + self.g3*scale3)*GSF
                target = sum((self.gt - g_1_2_3)**2)
                return target

            def object_func_LAD(para_list):
                scale1 = para_list[0]
                scale2 = para_list[1]
                scale3 = 1 - scale1 - scale2
                GSF = para_list[2]
                g_1_2_3 = (self.g1*scale1 + scale2*self.g2 + self.g3*scale3)*GSF
                target = sum(numpy.abs(g_1_2_3 - self.gt))
                return target

            def object_Pearson(para_list):
                scale1 = para_list[0]
                scale2 = para_list[1]
                scale3 = 1 - scale1 - scale2
                GSF = para_list[2]
                g_1_2_3 = (self.g1*scale1 + scale2*self.g2 + self.g3*scale3)*GSF
                target = pearsonr(g_1_2_3, self.gt)[0]
                if target > 0:
                    return (-1)*target
                else:
                    return target
        
            self.x0_ = [self._scale_1, self._scale_2, self.global_scale_]
            
            if self.optimizers.get() == "Least Squares":
                self.x0_ = minimize(object_func_LS, self.x0_, method = "nelder-mead").x
        
            elif self.optimizers.get() == "Pearson":
                self.x0_ = minimize(object_Pearson, self.x0_, method = "nelder-mead").x

            elif self.optimizers.get() == "Least Absolute Deviation":
                self.x0_ = minimize(object_func_LAD, self.x0_, method = "nelder-mead").x
        
        
            self.sum_fit = (self.g1*self.x0_[0] + self.g2*self.x0_[1] + self.g3*(1 - self.x0_[0] - self.x0_[1]))*self.x0_[2]
            diffzero = -0.8*max(self.gt)*numpy.ones_like(self.gt)
            diff = self.gt - self.sum_fit + diffzero
            
            
            self.fig.clf()
            self.ax1 = self.fig.add_subplot(211)
            
            self.ax1.set_title("Plots of individual PDFs")
            self.ax1.plot(self.r1, self.g1, "k-", lw=2, label = "PDF 1")
            self.ax1.plot(self.r2, self.g2, "m-", lw=2, label = "PDF 2")
            self.ax1.plot(self.r3, self.g3, "y-", lw=2, label = "PDF 3")
            self.ax1.plot(self.rt, self.gt, "b-", lw=2, label = "PDF target")
            
            self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax1.xaxis.set_tick_params(labelsize=11)
            self.ax1.yaxis.set_tick_params(labelsize=11)
            
            self.ax1.legend(loc = 0)
            
            self.ax2 = self.fig.add_subplot(212)
            self.ax2.set_title("Optimization of PDFs to match the target one")
            self.ax2.plot(self.r2, self.sum_fit, "r-", lw=2, label = "PDF fit, scale1 = {0:.4f}, scale2 = {1:.4f}, global scale = {2:.4f}".format(self.x0_[0], self.x0_[1], self.x0_[2]))
            self.ax2.plot(self.rt, self.gt, "b-", lw=2, label = "PDF target")
            self.ax2.plot(self.r2, diff, "g-", lw=2, label = "PDF difference")
            self.ax2.plot(self.r2, diffzero, "k--", lw=1)
            
            self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax2.xaxis.set_tick_params(labelsize=11)
            self.ax2.yaxis.set_tick_params(labelsize=11)
            
            self.ax2.legend(loc = 0)
            
            self.fig.tight_layout()
            self.canvas.draw()

    def save_the_fit(self):
        self.out_file_name = tkFileDialog.asksaveasfilename(defaultextension = ".txt", filetypes = [("Text Documents", "*.txt")])
        if self.pdf_1_ == 1 and self.pdf_2_ == 1 and self.pdf_3_ == 0 and self.pdf_t_ == 1:
            numpy.savetxt(self.out_file_name,
                          zip(self.r1,
                              self.g1,
                              self.g2,
                              self.gt,
                              self.sum_fit,
                              self.gt - self.sum_fit),
                          
                          header =
                          "****************************************************************************************\n"
                          "****************************************************************************************\n"
                          "******Here are the raw data for the plotting the curves.********************************\n"
                          "******Best fit scales PDF 1 by a factor of {0:.4f} and PDF 2 by a factor of {1:.4f}*****\n"
                          "*************************The global scale factor is {2:.4f} ****************************\n"
                          "******From left to right, the data correspond to ***************************************\n"
                          "******(1) r (2)PDF 1 (3) PDF 2 (4) PDF target (5)Fit (6) Diff****** ********************\n"
                          "****************************************************************************************\n".format(self.x0_[0], 1-self.x0_[0], self.x0_[1]))

        elif self.pdf_1_ == 1 and self.pdf_2_ == 1 and self.pdf_3_ == 1 and self.pdf_t_ == 1:
            numpy.savetxt(self.out_file_name,
                  zip(self.r1,
                      self.g1,
                      self.g2,
                      self.g3,
                      self.gt,
                      self.sum_fit,
                      self.gt - self.sum_fit),
                  
                  header =
                  "****************************************************************************************\n"
                  "****************************************************************************************\n"
                  "******Here are the raw data for the plotting the curves.********************************\n"
                  "******Best fit scales PDF 1 by a factor of {0:.4f} and PDF 2 by a factor of {1:.4f}*****\n"
                  "*************************The global scale factor is {2:.4f} ****************************\n"
                  "******From left to right, the data correspond to ***************************************\n"
                  "******(1) r (2)PDF 1 (3) PDF 2 (4) PDF 3 (5) PDF target (6)Fit (7) Diff****** ********************\n"
                  "****************************************************************************************\n".format(self.x0_[0], self.x0_[1], self.x0_[2]))


class C_InterPDF():
    def __init__(self, master):
        self.master = master
        self.master.title("Extract Intermolecular PDF of Crystalline Compounds")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 300) # width + height
        self.master.resizable(False, False)
        
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        self.style.configure('TNotebook.Tab', background = 'grey91', font = ("Times", 16, "bold"))
        
        ttk.Label(self.master, text = "The study of intermolecular contribution via subtracting PDF of a single molecule from the total PDF.", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()

        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures and Data", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 0, column = 3, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Molecule Structure", command = self.load_molecule,
                 style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Crystal Structure",command = self.load_crystal,
                 style = "TButton").grid(row = 1, column = 3, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Load PDF data",command = self.load_pdf_data,
                 style = "TButton").grid(row = 1, column = 6, columnspan = 3, padx = 5)

        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 2, column = 3, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

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
        self.rmin  = StringVar()
        self.rmax  = StringVar()
        self.rstep = StringVar()
        self.delta2_mol = StringVar()
        self.mol_pc_dpc = StringVar()

        ## these are the entries for tab 2
        self.qdamp_2 = StringVar()
        self.qbroad_2 = StringVar()
        self.qmin_2 = StringVar()
        self.qmax_2 = StringVar()
        self.delta2_cryst = StringVar()
        self.cryst_pc_dpc = StringVar()

        ##set default values
        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin.set(0.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
        
        self.qdamp_2.set(0.04)
        self.qbroad_2.set(0.02)
        self.qmin_2.set(0.0)
        self.qmax_2.set(20.0)
        self.delta2_cryst.set(0.0)
        self.cryst_pc_dpc.set("PDFCalculator")
        
        self.atom_den = StringVar()
        
        self.Y_scale_entry = StringVar()
        self.Y_scale_entry.set(1.0)


        #######################################Tab 1#########################################################
        #Qdamp
        ttk.Label(self.tab_1, text = "Qdamp", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qdamp_1).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_1, text = "Qbroad", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8,font = ("Times", 12),
                textvariable = self.qbroad_1).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_1, text = "Qmin", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qmin_1).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_1, text = "Qmax", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qmax_1).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                
        #rmin
        ttk.Label(self.tab_1, text = "Rmin", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rmin).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #rmax
        ttk.Label(self.tab_1, text = "Rmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rmax).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_1, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rstep).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #delta2
        ttk.Label(self.tab_1, text = "Delta2", justify = CENTER,
              font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.delta2_mol).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                            
        ttk.Label(self.tab_1, text = "Calculator", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")

        self.combobox_1 = ttk.Combobox(self.tab_1, font = ("Times", 12),
                            textvariable = self.mol_pc_dpc, values = ["DebyePDFCalculator","PDFCalculator"])
        self.combobox_1.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.tab_1, text = "Reset To Default", command = self.reset_to_default_mol,
                     style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_1, text = "Expand Thermals", command = self.load_mol_uiso_elements,
                    style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                 
    #######################################Tab 2#########################################################
        #Qdamp
        ttk.Label(self.tab_2, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.qdamp_2).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_2, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8,font = ("Times", 12),
        textvariable = self.qbroad_2).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_2, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.qmin_2).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_2, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.qmax_2).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rmin
        ttk.Label(self.tab_2, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.rmin).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #rmax
        ttk.Label(self.tab_2, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.rmax).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_2, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.rstep).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #delta2
        ttk.Label(self.tab_2, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.delta2_cryst).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Label(self.tab_2, text = "Calculator", justify = CENTER,
            font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")
    
        self.combobox_2 = ttk.Combobox(self.tab_2, font = ("Times", 12), textvariable = self.cryst_pc_dpc, values = ["PDFCalculator", "DebyePDFCalculator"])
        self.combobox_2.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Label(self.tab_2, text = "Atom Density", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 7, column = 2,  columnspan = 1, padx = 5, sticky = "sw")
        
        
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
                  textvariable = self.atom_den).grid(row = 7, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Button(self.tab_2, text = "Reset To Default", command = self.reset_to_default_cryst,
           style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_2, text = "Expand Thermals", command = self.load_cryst_uiso_elements,
                      style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        
        self.bottom_frame = ttk.Frame(self.master)
        self.bottom_frame.pack()

        ttk.Label(self.bottom_frame, text = "Step 3: Run The Simulation",justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 1, columnspan = 3, padx = 5, pady = 5)

        ttk.Button(self.bottom_frame, text = "Visualize", command = self.visualize,
                 style = "TButton").grid(row = 0, column = 4, columnspan = 3, padx = 5)


    def load_molecule(self):
        try:
            self.input_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                          filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.input_mole_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")


    def load_crystal(self):
        try:
            self.input_cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".cif",
                                                        filetypes = [("Text Documents", "*.cif")])
            Structure(filename = self.input_cryst_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The crystal structure file cannot be parsed properly by the program! Please check...")
    
    def load_pdf_data(self):
        try:
            self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                     filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
        
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")


    def load_mol_uiso_elements(self):
        
    
        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
            
        ##the labels for Us, the thermals
        ttk.Label(self.tab_1, text = "Uiso", justify = CENTER,
                          font = ("Times", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "Occ", justify = CENTER,
                          font = ("Times", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)
                

        for i, j in enumerate(self.mole_elements_list):
            
            setattr(self, "mol_{}_Uiso".format(j), StringVar())
            setattr(self, "mol_{}_Occ".format(j), StringVar())
            getattr(self, "mol_{}_Uiso".format(j)).set(0.005)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)

            ttk.Label(self.tab_1, text = self.mole_elements_list[i], justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "mol_{}_Uiso".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "mol_{}_Occ".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")

    def load_cryst_uiso_elements(self):
        
        self.cryst_elements = set(Structure(filename = self.input_cryst_struct).element) #return a dictionary
            
            #put elements into a list
        self.cryst_elements_list = []
        for i in self.cryst_elements:
            self.cryst_elements_list.append(i)
            ##the labels for Us
        ttk.Label(self.tab_2, text = "Uiso", justify = CENTER,
                      font = ("Times", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "Occ", justify = CENTER,
                    font = ("Times", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)

        for i, j in enumerate(self.cryst_elements):

            setattr(self, "cryst_{}_Uiso".format(j), StringVar())
            setattr(self, "cryst_{}_Occ".format(j), StringVar())
            getattr(self, "cryst_{}_Uiso".format(j)).set(0.005)
            getattr(self, "cryst_{}_Occ".format(j)).set(1.000)
    
            ttk.Label(self.tab_2, text = self.cryst_elements_list[i], justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "cryst_{}_Uiso".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "cryst_{}_Occ".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")

    def reset_to_default_mol(self):

        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin.set(0.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
    
        for i, j in enumerate(self.mole_elements_list):
        
            getattr(self, "mol_{}_Uiso".format(j)).set(0.005)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)
    
    def reset_to_default_cryst(self):
        
        self.qdamp_2.set(0.04)
        self.qbroad_2.set(0.02)
        self.qmin_2.set(0.0)
        self.qmax_2.set(20.0)
        self.delta2_cryst.set(0.0)
        self.atom_den.set(loadStructure(self.input_cryst_struct).occupancy.sum()/loadStructure(self.input_cryst_struct).lattice.volume)
        self.cryst_pc_dpc.set("PDFCalculator")

        for i, j in enumerate(self.cryst_elements):
    
            getattr(self, "cryst_{}_Uiso".format(j)).set(0.005)
            getattr(self, "cryst_{}_Occ".format(j)).set(1.000)
    
    def visualize(self):  ###create a two-panel figure, no curves ###

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.top = tk.Toplevel()
        self.top.title("Intermolecular PDFs")
      
        self.vis_top_frame = ttk.Frame(self.top, padding = (10, 10))
        self.vis_top_frame.pack()
        self.fig = plt.figure(figsize=(10, 6), dpi=100) ##create a figure; modify the size here
        
        self.ax1 = self.fig.add_subplot(211)
        
        self.ax1.set_title("Individual PDFs")
        self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)


        self.ax2.set_title("Difference PDFs")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.vis_top_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        #self.canvas.draw()

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.vis_top_frame)
        #self.toolbar.pack()
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.vis_bottom_frame = ttk.Frame(self.top, padding = (10, 10))
        self.vis_bottom_frame.pack()

        ttk.Button(self.vis_bottom_frame, text = "Make The Plots", command = self.make_the_plot,
           style = "TButton").grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Fine Tune Parameters", command = self.fine_tune_parameters,
                      style = "TButton").grid(row = 0, column = 2, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Clear The Plots", command = self.clear_the_plot,
               style = "TButton").grid(row = 0, column = 4, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Save The Plots", command = self.save_the_plot,
               style = "TButton").grid(row = 0, column = 6, columnspan = 2, padx = 5, sticky = "sw")
    
    
    
    def plot_with_cfgs(self, mol_cfg_given, cryst_cfg_given):
        
        self.atom_den.set(loadStructure(self.input_cryst_struct).occupancy.sum()/loadStructure(self.input_cryst_struct).lattice.volume)


        if self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_pc(self.mol_struc)
            cryst_pc = PDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_pc(self.cryst_struc)

            try:

                self.fig.clf()
                profile = Profile()
                parser = PDFParser()
                parser.parseFile(self.input_pdf_gr)
                profile.loadParsedData(parser)
                
                self.r_exp = profile.x[:len(self.r1)]
                self.g_exp = profile.y[:len(self.r1)]
                
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_pc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_pc")
                self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_pc - mol_pc")
                self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1, "m-", lw=2, label = "exp - mol_pc")
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()
            
            except:
                self.fig.clf()

                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_pc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_pc")

                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_pc - mol_pc")
                self.ax2.legend(loc=0)
                self.fig.tight_layout()
                self.canvas.show()
    


        elif self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_pc(self.mol_struc)
            cryst_pc = DebyePDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_pc(self.cryst_struc)
            self.g2 = self.g2 - 4.0*numpy.pi*self.r2*float(self.atom_den.get())
            
            try:
                self.fig.clf()

                profile = Profile()
                parser = PDFParser()
                parser.parseFile(self.input_pdf_gr)
                profile.loadParsedData(parser)
                
                self.r_exp = profile.x[:len(self.r1)]
                self.g_exp = profile.y[:len(self.r1)]
                
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_pc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_dpc")
                self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")

                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_dpc - mol_pc")
                self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1, "k-", lw=2, label = "exp")
                
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()

            except:
    
                self.fig.clf()

                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_pc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_dpc")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_dpc - mol_pc")
                self.ax2.legend(loc=0)
                

                self.fig.tight_layout()
                self.canvas.show()


        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_dpc(self.mol_struc)
            cryst_pc = PDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_pc(self.cryst_struc)

            try:
                self.fig.clf()

                profile = Profile()
                parser = PDFParser()
                parser.parseFile(self.input_pdf_gr)
                profile.loadParsedData(parser)
                
                self.r_exp = profile.x[:len(self.r1)]
                self.g_exp = profile.y[:len(self.r1)]


                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_dpc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_pc")
                self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
                
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_pc - mol_dpc")
                self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1, "m-", lw=2, label = "exp - mol_dpc")
                
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()
            
            except:
                self.fig.clf()
        
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_dpc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_pc")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_pc - mol_dpc")
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()


        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_dpc(self.mol_struc)
            cryst_dpc = DebyePDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_dpc(self.cryst_struc)
            self.g2 = self.g2 - 4.0*numpy.pi*self.r2*float(self.atom_den.get())
            
            try:
                self.fig.clf()
     
                profile = Profile()
                parser = PDFParser()
                parser.parseFile(self.input_pdf_gr)
                profile.loadParsedData(parser)
                
                self.r_exp = profile.x[:len(self.r1)]
                self.g_exp = profile.y[:len(self.r1)]
                
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_dpc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_dpc")
                self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_dpc - mol_dpc")
                self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1, "k-", lw=2, label = "exp - mol_dpc")
                
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()
        
            except:
                self.fig.clf()
    
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_dpc")
                self.ax1.plot(self.r2, self.g2, "b-", lw=2, label = "cryst_dpc")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r1, self.g2 - self.g1, "g-", lw=2, label = "cryst_dpc - mol_dpc")
                self.ax2.legend(loc=0)
                self.fig.tight_layout()
                self.canvas.show()

    
    def make_the_plot(self):
    ###first the molecular PDF ###
        self.mol_cfg ={
        "qmax": float(self.qmax_1.get()),
        "qmin": float(self.qmin_1.get()),
        "rmin": float(self.rmin.get()),
        "rmax": float(self.rmax.get()),
        "qdamp": float(self.qdamp_1.get()),
        "qbroad": float(self.qbroad_1.get()),
        "delta2": float(self.delta2_mol.get()),
        "rstep": float(self.rstep.get())
                       }

        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        self.mol_struc = Structure(filename =  self.input_mole_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
        
        #print self.mole_elements_list

        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "mol_{}_Uiso".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].occupancy = float(getattr(self, "mol_{}_Occ".format(j)).get())


        ###second the crystalline PDF ###
        self.cryst_cfg ={
        "qmax": float(self.qmax_2.get()),
        "qmin": float(self.qmin_2.get()),
        "rmin": float(self.rmin.get()),
        "rmax": float(self.rmax.get()),
        "qdamp": float(self.qdamp_2.get()),
        "qbroad": float(self.qbroad_2.get()),
        "delta2": float(self.delta2_cryst.get()),
        "rstep": float(self.rstep.get())
        }
        
        self.cryst_elements = set(Structure(filename = self.input_cryst_struct).element) #return a dictionary
        self.cryst_struc = Structure(filename =  self.input_cryst_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.cryst_elements_list = []
        for i in self.cryst_elements:
            self.cryst_elements_list.append(i)
        
        #print self.mole_elements_list
        
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].Uisoequiv = float(getattr(self, "cryst_{}_Uiso".format(j)).get())
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
        self.panel_tab_3 = ttk.Frame(self.panel_frame)
        
        self.panel_tab.add(self.panel_tab_1, text = "Molecule Parameters")
        self.panel_tab.add(self.panel_tab_2, text = "Crystal Parameters")
        self.panel_tab.add(self.panel_tab_3, text = "Experimental PDF")
    
    #########These are parameters for molecules###############
        ttk.Label(self.panel_tab_1, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")


        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        
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
        
        self.m_qmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 40.0, orient = HORIZONTAL, variable = self.qmax_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmin_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.m_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.m_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rstep_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.m_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_delta2_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_mol_entry, command = self.update_plot_fine_tune_para)
        self.m_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")
    
        self.qdamp_1_entry.set(float(self.qdamp_1.get()))
        self.qbroad_1_entry.set(float(self.qbroad_1.get()))
        self.qmin_1_entry.set(float(self.qmin_1.get()))
        self.qmax_1_entry.set(float(self.qmax_1.get()))
        self.rmin_entry.set(float(self.rmin.get()))
        self.rmax_entry.set(float(self.rmax.get()))
        self.rstep_entry.set(float(self.rstep.get()))
        self.delta2_mol_entry.set(float(self.delta2_mol.get()))

        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "2.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "40.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "100.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "200.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "10.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        
        em1 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qdamp_1_entry)
        em1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em1.bind("<Return>", self.update_plot_fine_tune_para)
        em2 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qbroad_1_entry)
        em2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em2.bind("<Return>", self.update_plot_fine_tune_para)
        em3 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qmin_1_entry)
        em3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em3.bind("<Return>", self.update_plot_fine_tune_para)
        em4 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qmax_1_entry)
        em4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em4.bind("<Return>", self.update_plot_fine_tune_para)
        em5 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rmin_entry)
        em5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em5.bind("<Return>", self.update_plot_fine_tune_para)
        em6 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rmax_entry)
        em6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em6.bind("<Return>", self.update_plot_fine_tune_para)
        em7 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rstep_entry)
        em7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em7.bind("<Return>", self.update_plot_fine_tune_para)
        em8 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.delta2_mol_entry)
        em8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em8.bind("<Return>", self.update_plot_fine_tune_para)
        
        ##add fine tune options for Uisos for molecules
        for i, j in enumerate(self.mole_elements_list):

            setattr(self, "ft_mol_{}_Uiso".format(j), StringVar())
            getattr(self, "ft_mol_{}_Uiso".format(j)).set(float(getattr(self, "mol_{}_Uiso".format(j)).get()))

            ttk.Label(self.panel_tab_1, text = self.mole_elements_list[i] + "_Uiso", justify = CENTER,
                      font = ("Times", 16, "bold")).grid(row = 8+i, column = 0, columnspan = 2, padx = 5, sticky = "sw")
            ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
                      font = ("Times", 16, "bold")).grid(row = 8+i, column = 2, columnspan = 1, padx = 5, sticky = "sw")
            self.Uiso_m_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = getattr(self, "ft_mol_{}_Uiso".format(j)), command = self.update_plot_fine_tune_para)
            self.Uiso_m_scale.grid(row = 8+i, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
            ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            uiso_m_entry = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),
                                     textvariable = getattr(self, "ft_mol_{}_Uiso".format(j)))
            uiso_m_entry.grid(row = 8+i, column = 8, columnspan = 2, padx = 5, sticky = "sw")
            uiso_m_entry.bind("<Return>", self.update_plot_fine_tune_para)
    
        ttk.Button(self.panel_tab_1, text = "Reset All", command = self.Reset_m_c_e_Entry,
                   style = "TButton").grid(row = 8 + len(self.mole_elements_list), column = 0, columnspan = 2, padx = 5)

        #########These are parameters for crystals###############
        ttk.Label(self.panel_tab_2, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")

        self.qdamp_2_entry = StringVar()
        self.qbroad_2_entry = StringVar()
        self.qmin_2_entry = StringVar()
        self.qmax_2_entry = StringVar()
        self.delta2_cryst_entry = StringVar()

        self.c_qdamp_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 0.5, orient = HORIZONTAL, variable = self.qdamp_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qdamp_scale.grid(row = 0, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_qbroad_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 0.5, orient = HORIZONTAL, variable = self.qbroad_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qbroad_scale.grid(row = 1, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_qmin_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 2.0, orient = HORIZONTAL, variable = self.qmin_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qmin_scale.grid(row = 2, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_qmax_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 40.0, orient = HORIZONTAL, variable = self.qmax_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_rmin_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.c_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        self.c_rmax_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.c_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_rstep_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.c_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_delta2_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_cryst_entry, command = self.update_plot_fine_tune_para)
        self.c_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.qdamp_2_entry.set(float(self.qdamp_2.get()))
        self.qbroad_2_entry.set(float(self.qbroad_2.get()))
        self.qmin_2_entry.set(float(self.qmin_2.get()))
        self.qmax_2_entry.set(float(self.qmax_2.get()))
        self.delta2_cryst_entry.set(float(self.delta2_cryst.get()))

        ttk.Label(self.panel_tab_2, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "2.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "40.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "100.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "200.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "1.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "10.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        

        ec1 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qdamp_2_entry)
        ec1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec1.bind("<Return>", self.update_plot_fine_tune_para)
        ec2 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qbroad_2_entry)
        ec2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec2.bind("<Return>", self.update_plot_fine_tune_para)
        ec3 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qmin_2_entry)
        ec3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec3.bind("<Return>", self.update_plot_fine_tune_para)
        ec4 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qmax_2_entry)
        ec4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec4.bind("<Return>", self.update_plot_fine_tune_para)
        ec5 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.rmin_entry)
        ec5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec5.bind("<Return>", self.update_plot_fine_tune_para)
        ec6 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.rmax_entry)
        ec6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec6.bind("<Return>", self.update_plot_fine_tune_para)
        ec7 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.rstep_entry)
        ec7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec7.bind("<Return>", self.update_plot_fine_tune_para)
        ec8 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.delta2_cryst_entry)
        ec8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec8.bind("<Return>", self.update_plot_fine_tune_para)

        ##add fine tune options for Uisos for crystals
        for i, j in enumerate(self.cryst_elements_list):

            setattr(self, "ft_cryst_{}_Uiso".format(j), StringVar())
            getattr(self, "ft_cryst_{}_Uiso".format(j)).set(float(getattr(self, "cryst_{}_Uiso".format(j)).get()))
                    
            ttk.Label(self.panel_tab_2, text = self.cryst_elements_list[i] + "_Uiso", justify = CENTER,
                              font = ("Times", 16, "bold")).grid(row = 8+i, column = 0, columnspan = 2, padx = 5, sticky = "sw")
            ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 2, columnspan = 1, padx = 5, sticky = "sw")
            self.Uiso_c_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = getattr(self, "ft_cryst_{}_Uiso".format(j)), command = self.update_plot_fine_tune_para)
            self.Uiso_c_scale.grid(row = 8+i, column = 3, columnspan = 4, padx = 5, sticky = "sw")
                              
            ttk.Label(self.panel_tab_2, text = "1.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            uiso_c_entry = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12), textvariable = getattr(self, "ft_cryst_{}_Uiso".format(j)))
            uiso_c_entry.grid(row = 8+i, column = 8, columnspan = 2, padx = 5, sticky = "sw")
            uiso_c_entry.bind("<Return>", self.update_plot_fine_tune_para)

        ttk.Button(self.panel_tab_2, text = "Reset All", command = self.Reset_m_c_e_Entry,
           style = "TButton").grid(row = 8 + len(self.cryst_elements_list), column = 0, columnspan = 2, padx = 5)
        
        ##########there are the paramters for tuning experimental PDFs#############
        
        ttk.Label(self.panel_tab_3, text = "Y_scale", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "0.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        self.y_scale = ttk.Scale(self.panel_tab_3, from_ = 0, to = 10, orient = HORIZONTAL, variable = self.Y_scale_entry,
                                       command = self.update_plot_fine_tune_para)
        self.y_scale.grid(row = 0, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "10.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
                  
        ee1 = ttk.Entry(self.panel_tab_3, width = 4, font = ("Times", 12),textvariable = self.Y_scale_entry)
        ee1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ee1.bind("<Return>", self.update_plot_fine_tune_para)

        ttk.Button(self.panel_tab_3, text = "Reset All", command = self.Reset_m_c_e_Entry, style = "TButton").grid(row = 1, column = 0,
                   columnspan = 2, padx = 5)
    
    
    def update_plot_fine_tune_para(self, event):
    
        self.fig.clf()
        ###first the molecular PDF ###
        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "ft_mol_{}_Uiso".format(j)).get())
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].Uisoequiv = float(getattr(self, "ft_cryst_{}_Uiso".format(j)).get())

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

    def Reset_m_c_e_Entry(self):

        self.fig.clf()

        for j in self.mole_elements_list:
            getattr(self, "ft_mol_{}_Uiso".format(j)).set(0.005)
        
        for j in self.cryst_elements_list:
            getattr(self, "ft_cryst_{}_Uiso".format(j)).set(0.005)
        
        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "ft_mol_{}_Uiso".format(j)).get())
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].Uisoequiv = float(getattr(self, "ft_cryst_{}_Uiso".format(j)).get())

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
        
        self.Y_scale_entry.set(1.0)
        
        
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

        self.ax1 = self.fig.add_subplot(211)
        
        self.ax1.set_title("Individual PDFs")
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.set_title("Difference PDFs")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()
        self.canvas.draw()

    def save_the_plot(self):
        self.out_file_name = tkFileDialog.asksaveasfilename(defaultextension = ".txt", filetypes = [("Text Documents", "*.txt")])
        ##Scenario 1
        if self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator":
            
            try:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g_exp*float(self.Y_scale_entry.get()),
                    self.g2 - self.g1,
                    self.g_exp*float(self.Y_scale_entry.get()) - self.g1),
                header =
                      "***********************************************************************************\n"
                      "***********************************************************************************\n"
                      "******Here are the raw data for the plotting the curves.***************************\n"
                      "******You have chosen to use PDFCalculator for both molecule and crystal.*********\n"
                      "******You also have scaled the measured PDF by a factor of {0:.3f} ****************\n"
                      "******From left to right, the data correspond to (1)r (2)m (3)c (4)e (5)c-m (6)e-m \n"
                      "***********************************************************************************\n"
                      "***********************************************************************************\n".format(float(self.Y_scale_entry.get())))
        
            except:
                
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g2 - self.g1
                    ),
                header =
                      "*****************************************************************************\n"
                      "*****************************************************************************\n"
                      "*********Here are the raw data for the plotting the curves.******************\n"
                      "*********You have chosen to use PDFCalculator for both molecule and crystal.\n"
                      "*********From left to right, the data correspond to (1)r (2)m (3)c (4)c-m ***\n"
                      "*****************************************************************************\n"
                      "*****************************************************************************\n")



        elif self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC

            try:
                
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                  self.g1,
                  self.g2,
                  self.g_exp*float(self.Y_scale_entry.get()),
                  self.g2 - self.g1,
                  self.g_exp*float(self.Y_scale_entry.get()) - self.g1),
                
                header =
                  "****************************************************************************************\n"
                  "****************************************************************************************\n"
                  "******Here are the raw data for the plotting the curves.********************************\n"
                  "******You have chosen to use PDFCalculator for molecule; DebyePDFCalculator for crystal.\n"
                  "******You also have scaled the measured PDF by a factor of {0:.3f} *********************\n"
                  "******From left to right, the data correspond to (1)r (2)m (3)c (4)e (5)c-m (6)e-m *****\n"
                  "****************************************************************************************\n"
                  "****************************************************************************************\n".format(float(self.Y_scale_entry.get())))
        
            except:

                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g2 - self.g1
                    ),
                header =
                              "*********************************************************************************************\n"
                              "*********************************************************************************************\n"
                              "*********Here are the raw data for the plotting the curves.***********************************\n"
                              "*********You have chosen to use PDFCalculator for molecule and DebyePDFCalculator for crystal.\n"
                              "*********From left to right, the data correspond to (1)r (2)m (3)c (4)c-m ********************\n"
                              "**********************************************************************************************\n"
                              "**********************************************************************************************\n")


        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            try:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g_exp*float(self.Y_scale_entry.get()),
                    self.g2 - self.g1,
                    self.g_exp*float(self.Y_scale_entry.get()) - self.g1),
                header =
                  "************************************************************************************************\n"
                  "************************************************************************************************\n"
                  "******Here are the raw data for the plotting the curves.****************************************\n"
                  "******You have chosen to use DeybePDFCalculator for molecule; PDFCalculator for crystal.*********\n"
                  "******You also have scaled the measured PDF by a factor of {0:.3f} ******************************\n"
                  "******From left to right, the data correspond to (1)r (2)m (3)c (4)e (5)c-m (6)e-m **************\n"
                  "*************************************************************************************************\n"
                  "*************************************************************************************************\n".format(float(self.Y_scale_entry.get())))
        
            except:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g2 - self.g1
                    ),
                header =
                              "******************************************************************************************\n"
                              "******************************************************************************************\n"
                              "*********Here are the raw data for the plotting the curves.*******************************\n"
                              "*********You have chosen to use DebyePDFCalculator for molecule; PDFCalculator for crystal.\n"
                              "*********From left to right, the data correspond to (1)r (2)m (3)c (4)c-m *****************\n"
                              "*******************************************************************************************\n"
                              "*******************************************************************************************\n")



        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC
            try:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g_exp*float(self.Y_scale_entry.get()),
                    self.g2 - self.g1,
                    self.g_exp*float(self.Y_scale_entry.get()) - self.g1),
                header =
                      "****************************************************************************************\n"
                      "****************************************************************************************\n"
                      "******Here are the raw data for the plotting the curves.********************************\n"
                      "******You have chosen to use DebyePDFCalculator for both molecule and crystal.*********\n"
                      "******You also have scaled the measured PDF by a factor of {0:.3f} *********************\n"
                      "******From left to right, the data correspond to (1)r (2)m (3)c (4)e (5)c-m (6)e-m *****\n"
                      "****************************************************************************************\n"
                      "****************************************************************************************\n".format(float(self.Y_scale_entry.get())))
            
            except:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g2,
                    self.g2 - self.g1
                    ),
                header =
                              "**********************************************************************************\n"
                              "**********************************************************************************\n"
                              "*********Here are the raw data for the plotting the curves.***********************\n"
                              "*********You have chosen to use DebyePDFCalculator for both molecule and crystal.\n"
                              "*********From left to right, the data correspond to (1)r (2)m (3)c (4)c-m ********\n"
                              "**********************************************************************************\n"
                              "**********************************************************************************\n")

class A_InterPDF():
    def __init__(self, master):
        self.master = master
        self.master.title("Extract Intermolecular PDF of Amorphous Compounds")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 300) # width + height
        self.master.resizable(False, False)
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        
        ttk.Label(self.master, text = "The study of intermolecular contribution via subtracting PDF of a single molecule from the total PDF.", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()

        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 0, column = 1, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Molecule Structure", command = self.load_molecule,
                 style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load PDF data",command = self.load_pdf_data,
                 style = "TButton").grid(row = 1, column = 3, columnspan = 3, padx = 5)

        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 2, column = 1, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        self.middle_frame = ttk.Frame(self.master)
        self.middle_frame.pack()
        
        ######adding tabs##############
        self.nb = ttk.Notebook(self.middle_frame)
        self.nb.pack(fill = BOTH, expand = "yes")
        
        self.tab_1 = ttk.Frame(self.middle_frame)
        self.nb.add(self.tab_1, text = "Molecule Parameters")

        ## these are the entries for tab 1
        self.qdamp_1 = StringVar()
        self.qbroad_1 = StringVar()
        self.qmin_1 = StringVar()
        self.qmax_1 = StringVar()
        self.rmin  = StringVar()
        self.rmax  = StringVar()
        self.rstep = StringVar()
        self.delta2_mol = StringVar()
        self.mol_pc_dpc = StringVar()

        ##set default values
        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin.set(0.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
        
        self.Y_scale_entry = StringVar()
        self.Y_scale_entry.set(1.0)

        #######################################Tab 1#########################################################
        #Qdamp
        ttk.Label(self.tab_1, text = "Qdamp", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qdamp_1).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_1, text = "Qbroad", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8,font = ("Times", 12),
                textvariable = self.qbroad_1).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_1, text = "Qmin", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qmin_1).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_1, text = "Qmax", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qmax_1).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                
        #rmin
        ttk.Label(self.tab_1, text = "Rmin", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rmin).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #rmax
        ttk.Label(self.tab_1, text = "Rmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rmax).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_1, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rstep).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #delta2
        ttk.Label(self.tab_1, text = "Delta2", justify = CENTER,
              font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.delta2_mol).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                            
        ttk.Label(self.tab_1, text = "Calculator", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")

        self.combobox_1 = ttk.Combobox(self.tab_1, font = ("Times", 12),
                            textvariable = self.mol_pc_dpc, values = ["DebyePDFCalculator","PDFCalculator"])
        self.combobox_1.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.tab_1, text = "Reset To Default", command = self.reset_to_default_mol,
                     style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_1, text = "Expand Thermals", command = self.load_mol_uiso_elements,
                    style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                
        self.bottom_frame = ttk.Frame(self.master)
        self.bottom_frame.pack()

        ttk.Label(self.bottom_frame, text = "Step 3: Run The Simulation",justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 1, columnspan = 3, padx = 5, pady = 5)

        ttk.Button(self.bottom_frame, text = "Visualize", command = self.visualize,
                 style = "TButton").grid(row = 0, column = 4, columnspan = 3, padx = 5)


    def load_molecule(self):
        try:
            self.input_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                          filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.input_mole_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")


    def load_crystal(self):
        try:
            self.input_cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".cif",
                                                        filetypes = [("Text Documents", "*.cif")])
            Structure(filename = self.input_cryst_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The crystal structure file cannot be parsed properly by the program! Please check...")
    
    def load_pdf_data(self):
        try:
            self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                     filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
        
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def load_mol_uiso_elements(self):
        
    
        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
            
        ##the labels for Us, the thermals
        ttk.Label(self.tab_1, text = "Uiso", justify = CENTER,
                          font = ("Times", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "Occ", justify = CENTER,
                          font = ("Times", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)
                

        for i, j in enumerate(self.mole_elements_list):
            
            setattr(self, "mol_{}_Uiso".format(j), StringVar())
            setattr(self, "mol_{}_Occ".format(j), StringVar())
            getattr(self, "mol_{}_Uiso".format(j)).set(0.005)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)

            ttk.Label(self.tab_1, text = self.mole_elements_list[i], justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "mol_{}_Uiso".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "mol_{}_Occ".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")

    def reset_to_default_mol(self):

        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin.set(0.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
    
        for i, j in enumerate(self.mole_elements_list):
        
            getattr(self, "mol_{}_Uiso".format(j)).set(0.005)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)
    
    def visualize(self):  ###create a two-panel figure, no curves ###

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.top = tk.Toplevel()
        self.top.title("Intermolecular PDFs")
      
        self.vis_top_frame = ttk.Frame(self.top, padding = (10, 10))
        self.vis_top_frame.pack()
        self.fig = plt.figure(figsize=(10, 6), dpi=100) ##create a figure; modify the size here
        
        self.ax1 = self.fig.add_subplot(211)
        
        self.ax1.set_title("Individual PDFs")
        self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)

        self.ax2.set_title("Difference PDFs")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.vis_top_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        #self.canvas.draw()

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.vis_top_frame)
        #self.toolbar.pack()
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.vis_bottom_frame = ttk.Frame(self.top, padding = (10, 10))
        self.vis_bottom_frame.pack()

        ttk.Button(self.vis_bottom_frame, text = "Make The Plots", command = self.make_the_plot,
           style = "TButton").grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Fine Tune Parameters", command = self.fine_tune_parameters,
                      style = "TButton").grid(row = 0, column = 2, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Clear The Plots", command = self.clear_the_plot,
               style = "TButton").grid(row = 0, column = 4, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Save The Plots", command = self.save_the_plot,
               style = "TButton").grid(row = 0, column = 6, columnspan = 2, padx = 5, sticky = "sw")
    
    
    
    def plot_with_cfgs(self, mol_cfg_given):
        
        #self.atom_den.set(loadStructure(self.input_cryst_struct).occupancy.sum()/loadStructure(self.input_cryst_struct).lattice.volume)


        if self.mol_pc_dpc.get() == "PDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_pc(self.mol_struc)

            try:

                self.fig.clf()
                profile = Profile()
                parser = PDFParser()
                parser.parseFile(self.input_pdf_gr)
                profile.loadParsedData(parser)
                
                self.r_exp = profile.x[:len(self.r1)]
                self.g_exp = profile.y[:len(self.r1)]
                
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_pc")
                self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1, "m-", lw=2, label = "exp - mol_pc")
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()
            
            except:
                self.fig.clf()

                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_pc")

                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)

                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()


        elif self.mol_pc_dpc.get() == "DebyePDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_dpc(self.mol_struc)

            try:
                self.fig.clf()

                profile = Profile()
                parser = PDFParser()
                parser.parseFile(self.input_pdf_gr)
                profile.loadParsedData(parser)
                
                self.r_exp = profile.x[:len(self.r1)]
                self.g_exp = profile.y[:len(self.r1)]
    
    
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_dpc")
                self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
                
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1, "m-", lw=2, label = "exp - mol_dpc")
                
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()
            
            except:
                self.fig.clf()
        
                self.ax1 = self.fig.add_subplot(211)
                self.ax1.set_title("Individual PDFs")
                self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax1.xaxis.set_tick_params(labelsize=11)
                self.ax1.yaxis.set_tick_params(labelsize=11)
                self.ax1.plot(self.r1, self.g1, "r-", lw=2, label = "mol_dpc")
                
                self.ax1.legend(loc=0)
                
                self.ax2 = self.fig.add_subplot(212)
                self.ax2.set_title("Difference PDFs")
                self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
                self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
                self.ax2.xaxis.set_tick_params(labelsize=11)
                self.ax2.yaxis.set_tick_params(labelsize=11)
                self.ax2.legend(loc=0)
                
                self.fig.tight_layout()
                self.canvas.show()

    
    def make_the_plot(self):
    ###first the molecular PDF ###
        self.mol_cfg ={
        "qmax": float(self.qmax_1.get()),
        "qmin": float(self.qmin_1.get()),
        "rmin": float(self.rmin.get()),
        "rmax": float(self.rmax.get()),
        "qdamp": float(self.qdamp_1.get()),
        "qbroad": float(self.qbroad_1.get()),
        "delta2": float(self.delta2_mol.get()),
        "rstep": float(self.rstep.get())
                       }

        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        self.mol_struc = Structure(filename =  self.input_mole_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
        
        #print self.mole_elements_list

        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "mol_{}_Uiso".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].occupancy = float(getattr(self, "mol_{}_Occ".format(j)).get())

        self.plot_with_cfgs(self.mol_cfg)

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
        self.panel_tab.add(self.panel_tab_2, text = "Experimental PDF")
    
    #########These are parameters for molecules###############
        ttk.Label(self.panel_tab_1, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")


        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        
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
        
        self.m_qmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 40.0, orient = HORIZONTAL, variable = self.qmax_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmin_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.m_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.m_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rstep_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.m_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_delta2_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_mol_entry, command = self.update_plot_fine_tune_para)
        self.m_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")
    
        self.qdamp_1_entry.set(float(self.qdamp_1.get()))
        self.qbroad_1_entry.set(float(self.qbroad_1.get()))
        self.qmin_1_entry.set(float(self.qmin_1.get()))
        self.qmax_1_entry.set(float(self.qmax_1.get()))
        self.rmin_entry.set(float(self.rmin.get()))
        self.rmax_entry.set(float(self.rmax.get()))
        self.rstep_entry.set(float(self.rstep.get()))
        self.delta2_mol_entry.set(float(self.delta2_mol.get()))

        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "2.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "40.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "100.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "200.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "10.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        
        em1 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qdamp_1_entry)
        em1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em1.bind("<Return>", self.update_plot_fine_tune_para)
        em2 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qbroad_1_entry)
        em2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em2.bind("<Return>", self.update_plot_fine_tune_para)
        em3 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qmin_1_entry)
        em3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em3.bind("<Return>", self.update_plot_fine_tune_para)
        em4 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qmax_1_entry)
        em4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em4.bind("<Return>", self.update_plot_fine_tune_para)
        em5 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rmin_entry)
        em5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em5.bind("<Return>", self.update_plot_fine_tune_para)
        em6 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rmax_entry)
        em6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em6.bind("<Return>", self.update_plot_fine_tune_para)
        em7 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rstep_entry)
        em7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em7.bind("<Return>", self.update_plot_fine_tune_para)
        em8 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.delta2_mol_entry)
        em8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em8.bind("<Return>", self.update_plot_fine_tune_para)
        
        ##add fine tune options for Uisos for molecules
        for i, j in enumerate(self.mole_elements_list):

            setattr(self, "ft_mol_{}_Uiso".format(j), StringVar())
            getattr(self, "ft_mol_{}_Uiso".format(j)).set(float(getattr(self, "mol_{}_Uiso".format(j)).get()))

            ttk.Label(self.panel_tab_1, text = self.mole_elements_list[i] + "_Uiso", justify = CENTER,
                      font = ("Times", 16, "bold")).grid(row = 8+i, column = 0, columnspan = 2, padx = 5, sticky = "sw")
            ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
                      font = ("Times", 16, "bold")).grid(row = 8+i, column = 2, columnspan = 1, padx = 5, sticky = "sw")
            self.Uiso_m_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = getattr(self, "ft_mol_{}_Uiso".format(j)), command = self.update_plot_fine_tune_para)
            self.Uiso_m_scale.grid(row = 8+i, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
            ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            uiso_m_entry = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),
                                     textvariable = getattr(self, "ft_mol_{}_Uiso".format(j)))
            uiso_m_entry.grid(row = 8+i, column = 8, columnspan = 2, padx = 5, sticky = "sw")
            uiso_m_entry.bind("<Return>", self.update_plot_fine_tune_para)
    
        ttk.Button(self.panel_tab_1, text = "Reset All", command = self.Reset_m_e_Entry,
                   style = "TButton").grid(row = 8 + len(self.mole_elements_list), column = 0, columnspan = 2, padx = 5)
        
        ##########there are the paramters for tuning experimental PDFs#############
        
        ttk.Label(self.panel_tab_2, text = "Y_scale", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        self.y_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 10, orient = HORIZONTAL, variable = self.Y_scale_entry,
                                       command = self.update_plot_fine_tune_para)
        self.y_scale.grid(row = 0, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "10.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
                  
        ee1 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.Y_scale_entry)
        ee1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ee1.bind("<Return>", self.update_plot_fine_tune_para)

        ttk.Button(self.panel_tab_2, text = "Reset All", command = self.Reset_m_e_Entry, style = "TButton").grid(row = 1, column = 0,
                   columnspan = 2, padx = 5)
    
    
    def update_plot_fine_tune_para(self, event):
    
        self.fig.clf()
        ###first the molecular PDF ###
        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "ft_mol_{}_Uiso".format(j)).get())

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
        
        self.plot_with_cfgs(self.mol_fine_tune_cfg)

    def Reset_m_e_Entry(self):

        self.fig.clf()

        for j in self.mole_elements_list:
            getattr(self, "ft_mol_{}_Uiso".format(j)).set(0.005)
        
        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "ft_mol_{}_Uiso".format(j)).get())

        self.qdamp_1_entry.set(0.04)
        self.qbroad_1_entry.set(0.02)
        self.qmin_1_entry.set(0.0)
        self.qmax_1_entry.set(20.0)
        self.rmin_entry.set(0.0)
        self.rmax_entry.set(20.0)
        self.rstep_entry.set(0.01)
        self.delta2_mol_entry.set(0.0)
        
        self.Y_scale_entry.set(1.0)
        
        
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
        
        
        self.plot_with_cfgs(self.mol_reset_cfg)


    def clear_the_plot(self):
        self.fig.clf()

        self.ax1 = self.fig.add_subplot(211)
        
        self.ax1.set_title("Individual PDFs")
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.set_title("Difference PDFs")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()
        self.canvas.draw()

    def save_the_plot(self):
        self.out_file_name = tkFileDialog.asksaveasfilename(defaultextension = ".txt", filetypes = [("Text Documents", "*.txt")])
        ##Scenario 1
        if self.mol_pc_dpc.get() == "PDFCalculator":
            
            try:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g_exp*float(self.Y_scale_entry.get()),
                    self.g_exp*float(self.Y_scale_entry.get()) - self.g1),
                header =
                      "***********************************************************************************\n"
                      "***********************************************************************************\n"
                      "******Here are the raw data for the plotting the curves.***************************\n"
                      "******You have chosen to use PDFCalculator for molecule.***************************\n"
                      "******You also have scaled the measured PDF by a factor of {0:.3f} ****************\n"
                      "******From left to right, the data correspond to (1)r (2)m (3)e (4)e-m ************\n"
                      "***********************************************************************************\n"
                      "***********************************************************************************\n".format(float(self.Y_scale_entry.get())))
        
            except:
                
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1),
                              
                header =
                      "*****************************************************************************\n"
                      "*****************************************************************************\n"
                      "*********Here are the raw data for the plotting the curves.******************\n"
                      "*********You have chosen to use PDFCalculator for molecule.******************\n"
                      "*********From left to right, the data correspond to (1)r (2)m ***************\n"
                      "*****************************************************************************\n"
                      "*****************************************************************************\n")


        elif self.mol_pc_dpc.get() == "DebyePDFCalculator": # DPC
            try:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1,
                    self.g_exp*float(self.Y_scale_entry.get()),
                    self.g_exp*float(self.Y_scale_entry.get()) - self.g1),
                header =
                      "****************************************************************************************\n"
                      "****************************************************************************************\n"
                      "******Here are the raw data for the plotting the curves.********************************\n"
                      "******You have chosen to use DebyePDFCalculator for molecule.***************************\n"
                      "******You also have scaled the measured PDF by a factor of {0:.3f} *********************\n"
                      "******From left to right, the data correspond to (1)r (2)m (3)e (4)e-m *****************\n"
                      "****************************************************************************************\n"
                      "****************************************************************************************\n".format(float(self.Y_scale_entry.get())))
            
            except:
                numpy.savetxt(self.out_file_name,
                zip(self.r1,
                    self.g1
                    ),
                header =
                              "**********************************************************************************\n"
                              "**********************************************************************************\n"
                              "*********Here are the raw data for the plotting the curves.***********************\n"
                              "*********You have chosen to use DebyePDFCalculator for molecule.******************\n"
                              "*********From left to right, the data correspond to (1)r (2)m ********************\n"
                              "**********************************************************************************\n"
                              "**********************************************************************************\n")


class ASD_InterPDF():
    def __init__(self, master):
        self.master = master
        self.master.title("Extract Intermolecular PDF of ASDs")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 300) # width + height
        self.master.resizable(False, False)
        
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        
        ttk.Label(self.master, text = "The study of intermolecular contribution via subtracting PDF of a single molecule from the total PDF.", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()

        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 0, column = 3, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Molecule Structure 1", command = self.load_molecule,
                 style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Molecule Structure 2",command = self.load_crystal,
                 style = "TButton").grid(row = 1, column = 3, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Load PDF data",command = self.load_pdf_data,
                 style = "TButton").grid(row = 1, column = 6, columnspan = 3, padx = 5)

        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 2, column = 3, columnspan = 3, padx = 5, pady = 5, sticky = "sw")

        self.middle_frame = ttk.Frame(self.master)
        self.middle_frame.pack()
        
        ######adding tabs##############
        self.nb = ttk.Notebook(self.middle_frame)
        self.nb.pack(fill = BOTH, expand = "yes")
        
        self.tab_1 = ttk.Frame(self.middle_frame)
        self.tab_2 = ttk.Frame(self.middle_frame)
        
        self.nb.add(self.tab_1, text = "Molecule 1 Parameters")
        self.nb.add(self.tab_2, text = "Molecule 2 Parameters")

        ## these are the entries for tab 1
        self.qdamp_1 = StringVar()
        self.qbroad_1 = StringVar()
        self.qmin_1 = StringVar()
        self.qmax_1 = StringVar()
        self.rmin  = StringVar()
        self.rmax  = StringVar()
        self.rstep = StringVar()
        self.delta2_mol = StringVar()
        self.mol_pc_dpc = StringVar()

        ## these are the entries for tab 2
        self.qdamp_2 = StringVar()
        self.qbroad_2 = StringVar()
        self.qmin_2 = StringVar()
        self.qmax_2 = StringVar()
        self.delta2_cryst = StringVar()
        self.cryst_pc_dpc = StringVar()

        ##set default values
        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin.set(0.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
        
        self.qdamp_2.set(0.04)
        self.qbroad_2.set(0.02)
        self.qmin_2.set(0.0)
        self.qmax_2.set(20.0)
        self.delta2_cryst.set(0.0)
        self.cryst_pc_dpc.set("DebyePDFCalculator")
        
        self.Y_scale_entry = StringVar()
        self.Y_scale_entry.set(1.0)
        
        self.mol_1_Y_scale_entry = StringVar()
        self.mol_1_Y_scale_entry.set(1.0)
        
        self.mol_2_Y_scale_entry = StringVar()
        self.mol_2_Y_scale_entry.set(1.0)


        #######################################Tab 1#########################################################
        #Qdamp
        ttk.Label(self.tab_1, text = "Qdamp", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qdamp_1).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_1, text = "Qbroad", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8,font = ("Times", 12),
                textvariable = self.qbroad_1).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_1, text = "Qmin", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qmin_1).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_1, text = "Qmax", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.qmax_1).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                
        #rmin
        ttk.Label(self.tab_1, text = "Rmin", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rmin).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #rmax
        ttk.Label(self.tab_1, text = "Rmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rmax).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_1, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.rstep).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        #delta2
        ttk.Label(self.tab_1, text = "Delta2", justify = CENTER,
              font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = self.delta2_mol).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")
                            
        ttk.Label(self.tab_1, text = "Calculator", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")

        self.combobox_1 = ttk.Combobox(self.tab_1, font = ("Times", 12),
                            textvariable = self.mol_pc_dpc, values = ["DebyePDFCalculator","PDFCalculator"])
        self.combobox_1.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Label(self.tab_1, text = "Scale Factor", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 7, column = 2,  columnspan = 1, padx = 5, sticky = "sw")

        ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                  textvariable = self.mol_1_Y_scale_entry).grid(row = 7, column = 3, columnspan = 1, padx = 5, sticky = "sw")
    
        ttk.Button(self.tab_1, text = "Reset To Default", command = self.reset_to_default_mol,
                     style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_1, text = "Expand Thermals", command = self.load_mol_uiso_elements,
                    style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                 
    #######################################Tab 2#########################################################
        #Qdamp
        ttk.Label(self.tab_2, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.qdamp_2).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qbroad
        ttk.Label(self.tab_2, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8,font = ("Times", 12),
        textvariable = self.qbroad_2).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #Qmin
        ttk.Label(self.tab_2, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.qmin_2).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #Qmax
        ttk.Label(self.tab_2, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.qmax_2).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rmin
        ttk.Label(self.tab_2, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.rmin).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #rmax
        ttk.Label(self.tab_2, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.rmax).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.tab_2, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.rstep).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")

        #delta2
        ttk.Label(self.tab_2, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
        textvariable = self.delta2_cryst).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Label(self.tab_2, text = "Calculator", justify = CENTER,
            font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 1, padx = 5, sticky = "sw")
    
        self.combobox_2 = ttk.Combobox(self.tab_2, font = ("Times", 12), textvariable = self.cryst_pc_dpc, values = ["PDFCalculator", "DebyePDFCalculator"])
        self.combobox_2.grid(row = 7, column = 1, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Label(self.tab_2, text = "Scale Factor", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 7, column = 2,  columnspan = 1, padx = 5, sticky = "sw")
        
        ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
                  textvariable = self.mol_2_Y_scale_entry).grid(row = 7, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Button(self.tab_2, text = "Reset To Default", command = self.reset_to_default_cryst,
           style = "TButton").grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Button(self.tab_2, text = "Expand Thermals", command = self.load_cryst_uiso_elements,
                      style = "TButton").grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        
        
        self.bottom_frame = ttk.Frame(self.master)
        self.bottom_frame.pack()

        ttk.Label(self.bottom_frame, text = "Step 3: Run The Simulation",justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 1, columnspan = 3, padx = 5, pady = 5)

        ttk.Button(self.bottom_frame, text = "Visualize", command = self.visualize,
                 style = "TButton").grid(row = 0, column = 4, columnspan = 3, padx = 5)


    def load_molecule(self):
        try:
            self.input_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                          filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.input_mole_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")


    def load_crystal(self):
        try:
            self.input_cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                        filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.input_cryst_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The crystal structure file cannot be parsed properly by the program! Please check...")
    
    def load_pdf_data(self):
        try:
            self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr",
                                                     filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
        
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def load_mol_uiso_elements(self):
        
    
        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
            
        ##the labels for Us, the thermals
        ttk.Label(self.tab_1, text = "Uiso", justify = CENTER,
                          font = ("Times", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_1, text = "Occ", justify = CENTER,
                          font = ("Times", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)
                

        for i, j in enumerate(self.mole_elements_list):
            
            setattr(self, "mol_{}_Uiso".format(j), StringVar())
            setattr(self, "mol_{}_Occ".format(j), StringVar())
            getattr(self, "mol_{}_Uiso".format(j)).set(0.005)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)

            ttk.Label(self.tab_1, text = self.mole_elements_list[i], justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "mol_{}_Uiso".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_1, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "mol_{}_Occ".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")

    def load_cryst_uiso_elements(self):
        
        self.cryst_elements = set(Structure(filename = self.input_cryst_struct).element) #return a dictionary
            
            #put elements into a list
        self.cryst_elements_list = []
        for i in self.cryst_elements:
            self.cryst_elements_list.append(i)
            ##the labels for Us
        ttk.Label(self.tab_2, text = "Uiso", justify = CENTER,
                      font = ("Times", 14, "bold")).grid(row = 2, column = 5, columnspan = 1, padx = 5)
        ttk.Label(self.tab_2, text = "Occ", justify = CENTER,
                    font = ("Times", 14, "bold")).grid(row = 2, column = 6, columnspan = 1, padx = 5)

        for i, j in enumerate(self.cryst_elements):

            setattr(self, "cryst_{}_Uiso".format(j), StringVar())
            setattr(self, "cryst_{}_Occ".format(j), StringVar())
            getattr(self, "cryst_{}_Uiso".format(j)).set(0.005)
            getattr(self, "cryst_{}_Occ".format(j)).set(1.000)
    
            ttk.Label(self.tab_2, text = self.cryst_elements_list[i], justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 3+i, column = 4, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "cryst_{}_Uiso".format(j))).grid(row = 3+i, column = 5, columnspan = 1, padx = 5, sticky = "sw")
            ttk.Entry(self.tab_2, width = 8, font = ("Times", 12),
                textvariable = getattr(self, "cryst_{}_Occ".format(j))).grid(row = 3+i, column = 6, columnspan = 1, padx = 5, sticky = "sw")

    def reset_to_default_mol(self):

        self.qdamp_1.set(0.04)
        self.qbroad_1.set(0.02)
        self.qmin_1.set(0.0)
        self.qmax_1.set(20.0)
        self.rmin.set(0.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        self.delta2_mol.set(0.0)
        self.mol_pc_dpc.set("DebyePDFCalculator")
        self.mol_1_Y_scale_entry.set(1.0)
    
        for i, j in enumerate(self.mole_elements_list):
        
            getattr(self, "mol_{}_Uiso".format(j)).set(0.005)
            getattr(self, "mol_{}_Occ".format(j)).set(1.000)
    
    def reset_to_default_cryst(self):
        
        self.qdamp_2.set(0.04)
        self.qbroad_2.set(0.02)
        self.qmin_2.set(0.0)
        self.qmax_2.set(20.0)
        self.delta2_cryst.set(0.0)
        self.mol_2_Y_scale_entry.set(1.0)
        #self.atom_den.set(loadStructure(self.input_cryst_struct).occupancy.sum()/loadStructure(self.input_cryst_struct).lattice.volume)
        self.cryst_pc_dpc.set("DebyePDFCalculator")

        for i, j in enumerate(self.cryst_elements):
    
            getattr(self, "cryst_{}_Uiso".format(j)).set(0.005)
            getattr(self, "cryst_{}_Occ".format(j)).set(1.000)
    
    def visualize(self):  ###create a two-panel figure, no curves ###

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.top = tk.Toplevel()
        self.top.title("Intermolecular PDFs")
      
        self.vis_top_frame = ttk.Frame(self.top, padding = (10, 10))
        self.vis_top_frame.pack()
        self.fig = plt.figure(figsize=(10, 6), dpi=100) ##create a figure; modify the size here
        
        self.ax1 = self.fig.add_subplot(211)
        
        self.ax1.set_title("Individual PDFs")
        self.ax1.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)

        self.ax2.set_title("Difference PDFs")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.vis_top_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        #self.canvas.draw()

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.vis_top_frame)
        #self.toolbar.pack()
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.vis_bottom_frame = ttk.Frame(self.top, padding = (10, 10))
        self.vis_bottom_frame.pack()

        ttk.Button(self.vis_bottom_frame, text = "Make The Plots", command = self.make_the_plot,
           style = "TButton").grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Fine Tune Parameters", command = self.fine_tune_parameters,
                      style = "TButton").grid(row = 0, column = 2, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Clear The Plots", command = self.clear_the_plot,
               style = "TButton").grid(row = 0, column = 4, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.vis_bottom_frame, text = "Save The Plots", command = self.save_the_plot,
               style = "TButton").grid(row = 0, column = 6, columnspan = 2, padx = 5, sticky = "sw")
    
    
    
    def plot_with_cfgs(self, mol_cfg_given, cryst_cfg_given):

        if self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_pc(self.mol_struc)
            cryst_pc = PDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_pc(self.cryst_struc)

            self.fig.clf()
            profile = Profile()
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
            profile.loadParsedData(parser)
            
            self.r_exp = profile.x[:len(self.r1)]
            self.g_exp = profile.y[:len(self.r1)]
            
            self.ax1 = self.fig.add_subplot(211)
            self.ax1.set_title("Individual PDFs")
            self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax1.xaxis.set_tick_params(labelsize=11)
            self.ax1.yaxis.set_tick_params(labelsize=11)
            self.ax1.plot(self.r1, self.g1*float(self.mol_1_Y_scale_entry.get()), "r-", lw=2, label = "mol_1_pc")
            self.ax1.plot(self.r2, self.g2*float(self.mol_2_Y_scale_entry.get()), "b-", lw=2, label = "mol_2_pc")
            self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
            
            self.ax1.legend(loc=0)
            
            self.ax2 = self.fig.add_subplot(212)
            self.ax2.set_title("Difference PDFs")
            self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax2.xaxis.set_tick_params(labelsize=11)
            self.ax2.yaxis.set_tick_params(labelsize=11)
            self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1*float(self.mol_1_Y_scale_entry.get()) -self.g2*float(self.mol_2_Y_scale_entry.get()), "m-", lw=2, label = "exp - mol_1_pc - mol_2_pc")
            self.ax2.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()

        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_dpc(self.mol_struc)
            cryst_pc = PDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_pc(self.cryst_struc)
            
            self.fig.clf()
            profile = Profile()
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
            profile.loadParsedData(parser)
            
            self.r_exp = profile.x[:len(self.r1)]
            self.g_exp = profile.y[:len(self.r1)]
            
            self.ax1 = self.fig.add_subplot(211)
            self.ax1.set_title("Individual PDFs")
            self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax1.xaxis.set_tick_params(labelsize=11)
            self.ax1.yaxis.set_tick_params(labelsize=11)
            self.ax1.plot(self.r1, self.g1*float(self.mol_1_Y_scale_entry.get()), "r-", lw=2, label = "mol_1_dpc")
            self.ax1.plot(self.r2, self.g2*float(self.mol_2_Y_scale_entry.get()), "b-", lw=2, label = "mol_2_pc")
            self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
            
            self.ax1.legend(loc=0)
            
            self.ax2 = self.fig.add_subplot(212)
            self.ax2.set_title("Difference PDFs")
            self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax2.xaxis.set_tick_params(labelsize=11)
            self.ax2.yaxis.set_tick_params(labelsize=11)
            self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1*float(self.mol_1_Y_scale_entry.get()) -self.g2*float(self.mol_2_Y_scale_entry.get()), "m-", lw=2, label = "exp - mol_1_dpc - mol_2_pc")
            self.ax2.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()

        elif self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC
            mol_pc = PDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_pc(self.mol_struc)
            cryst_dpc = DebyePDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_dpc(self.cryst_struc)
            
            self.fig.clf()
            profile = Profile()
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
            profile.loadParsedData(parser)
            
            self.r_exp = profile.x[:len(self.r1)]
            self.g_exp = profile.y[:len(self.r1)]
            
            self.ax1 = self.fig.add_subplot(211)
            self.ax1.set_title("Individual PDFs")
            self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax1.xaxis.set_tick_params(labelsize=11)
            self.ax1.yaxis.set_tick_params(labelsize=11)
            self.ax1.plot(self.r1, self.g1*float(self.mol_1_Y_scale_entry.get()), "r-", lw=2, label = "mol_1_pc")
            self.ax1.plot(self.r2, self.g2*float(self.mol_2_Y_scale_entry.get()), "b-", lw=2, label = "mol_2_dpc")
            self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
            
            self.ax1.legend(loc=0)
            
            self.ax2 = self.fig.add_subplot(212)
            self.ax2.set_title("Difference PDFs")
            self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax2.xaxis.set_tick_params(labelsize=11)
            self.ax2.yaxis.set_tick_params(labelsize=11)
            self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1*float(self.mol_1_Y_scale_entry.get()) -self.g2*float(self.mol_2_Y_scale_entry.get()), "m-", lw=2, label = "exp - mol_1_pc - mol_2_dpc")
            self.ax2.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()

        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator": # PC
            mol_dpc = DebyePDFCalculator(**mol_cfg_given)
            self.r1, self.g1 = mol_dpc(self.mol_struc)
            cryst_dpc = DebyePDFCalculator(**cryst_cfg_given)
            self.r2, self.g2 = cryst_dpc(self.cryst_struc)
            
            self.fig.clf()
            profile = Profile()
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
            profile.loadParsedData(parser)
            
            self.r_exp = profile.x[:len(self.r1)]
            self.g_exp = profile.y[:len(self.r1)]
            
            self.ax1 = self.fig.add_subplot(211)
            self.ax1.set_title("Individual PDFs")
            self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax1.xaxis.set_tick_params(labelsize=11)
            self.ax1.yaxis.set_tick_params(labelsize=11)
            self.ax1.plot(self.r1, self.g1*float(self.mol_1_Y_scale_entry.get()), "r-", lw=2, label = "mol_1_dpc")
            self.ax1.plot(self.r2, self.g2*float(self.mol_2_Y_scale_entry.get()), "b-", lw=2, label = "mol_2_dpc")
            self.ax1.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()), "k-", lw=2, label = "exp")
            
            self.ax1.legend(loc=0)
            
            self.ax2 = self.fig.add_subplot(212)
            self.ax2.set_title("Difference PDFs")
            self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.ax2.xaxis.set_tick_params(labelsize=11)
            self.ax2.yaxis.set_tick_params(labelsize=11)
            self.ax2.plot(self.r_exp, self.g_exp*float(self.Y_scale_entry.get()) - self.g1*float(self.mol_1_Y_scale_entry.get()) -self.g2*float(self.mol_2_Y_scale_entry.get()), "m-", lw=2, label = "exp - mol_1_dpc - mol_2_dpc")
            self.ax2.legend(loc=0)
            
            self.fig.tight_layout()
            self.canvas.show()


    def make_the_plot(self):
    ###first the molecular PDF ###
        self.mol_cfg ={
        "qmax": float(self.qmax_1.get()),
        "qmin": float(self.qmin_1.get()),
        "rmin": float(self.rmin.get()),
        "rmax": float(self.rmax.get()),
        "qdamp": float(self.qdamp_1.get()),
        "qbroad": float(self.qbroad_1.get()),
        "delta2": float(self.delta2_mol.get()),
        "rstep": float(self.rstep.get())
                       }

        self.mole_elements = set(Structure(filename = self.input_mole_struct).element) #return a dictionary
        self.mol_struc = Structure(filename =  self.input_mole_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.mole_elements_list = []
        for i in self.mole_elements:
            self.mole_elements_list.append(i)
        
        #print self.mole_elements_list

        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "mol_{}_Uiso".format(j)).get())
            self.mol_struc[self.mol_struc.element == j].occupancy = float(getattr(self, "mol_{}_Occ".format(j)).get())


        ###second the crystalline PDF ###
        self.cryst_cfg ={
        "qmax": float(self.qmax_2.get()),
        "qmin": float(self.qmin_2.get()),
        "rmin": float(self.rmin.get()),
        "rmax": float(self.rmax.get()),
        "qdamp": float(self.qdamp_2.get()),
        "qbroad": float(self.qbroad_2.get()),
        "delta2": float(self.delta2_cryst.get()),
        "rstep": float(self.rstep.get())
        }
        
        self.cryst_elements = set(Structure(filename = self.input_cryst_struct).element) #return a dictionary
        self.cryst_struc = Structure(filename =  self.input_cryst_struct)
        #put elements into a list ["C", "H", "O", "N" ...]
        self.cryst_elements_list = []
        for i in self.cryst_elements:
            self.cryst_elements_list.append(i)
        
        #print self.mole_elements_list
        
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].Uisoequiv = float(getattr(self, "cryst_{}_Uiso".format(j)).get())
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
        self.panel_tab_3 = ttk.Frame(self.panel_frame)
        
        self.panel_tab.add(self.panel_tab_1, text = "Molecule 1 Parameters")
        self.panel_tab.add(self.panel_tab_2, text = "Molecule 2 Parameters")
        self.panel_tab.add(self.panel_tab_3, text = "Experimental PDF")
    
    #########These are parameters for molecules###############
        ttk.Label(self.panel_tab_1, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")


        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        
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
        
        self.m_qmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 40.0, orient = HORIZONTAL, variable = self.qmax_1_entry, command = self.update_plot_fine_tune_para)
        self.m_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmin_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.m_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rmax_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.m_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_rstep_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.m_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
        self.m_delta2_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_mol_entry, command = self.update_plot_fine_tune_para)
        self.m_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")
    
        self.qdamp_1_entry.set(float(self.qdamp_1.get()))
        self.qbroad_1_entry.set(float(self.qbroad_1.get()))
        self.qmin_1_entry.set(float(self.qmin_1.get()))
        self.qmax_1_entry.set(float(self.qmax_1.get()))
        self.rmin_entry.set(float(self.rmin.get()))
        self.rmax_entry.set(float(self.rmax.get()))
        self.rstep_entry.set(float(self.rstep.get()))
        self.delta2_mol_entry.set(float(self.delta2_mol.get()))

        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "2.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "40.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "100.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "200.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_1, text = "10.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        
        em1 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qdamp_1_entry)
        em1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em1.bind("<Return>", self.update_plot_fine_tune_para)
        em2 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qbroad_1_entry)
        em2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em2.bind("<Return>", self.update_plot_fine_tune_para)
        em3 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qmin_1_entry)
        em3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em3.bind("<Return>", self.update_plot_fine_tune_para)
        em4 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.qmax_1_entry)
        em4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em4.bind("<Return>", self.update_plot_fine_tune_para)
        em5 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rmin_entry)
        em5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em5.bind("<Return>", self.update_plot_fine_tune_para)
        em6 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rmax_entry)
        em6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em6.bind("<Return>", self.update_plot_fine_tune_para)
        em7 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.rstep_entry)
        em7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em7.bind("<Return>", self.update_plot_fine_tune_para)
        em8 = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),textvariable = self.delta2_mol_entry)
        em8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        em8.bind("<Return>", self.update_plot_fine_tune_para)
        
        ##add fine tune options for Uisos for molecules
        for i, j in enumerate(self.mole_elements_list):

            setattr(self, "ft_mol_{}_Uiso".format(j), StringVar())
            getattr(self, "ft_mol_{}_Uiso".format(j)).set(float(getattr(self, "mol_{}_Uiso".format(j)).get()))

            ttk.Label(self.panel_tab_1, text = self.mole_elements_list[i] + "_Uiso", justify = CENTER,
                      font = ("Times", 16, "bold")).grid(row = 8+i, column = 0, columnspan = 2, padx = 5, sticky = "sw")
            ttk.Label(self.panel_tab_1, text = "0.0", justify = CENTER,
                      font = ("Times", 16, "bold")).grid(row = 8+i, column = 2, columnspan = 1, padx = 5, sticky = "sw")
            self.Uiso_m_scale = ttk.Scale(self.panel_tab_1, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = getattr(self, "ft_mol_{}_Uiso".format(j)), command = self.update_plot_fine_tune_para)
            self.Uiso_m_scale.grid(row = 8+i, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        
            ttk.Label(self.panel_tab_1, text = "1.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            uiso_m_entry = ttk.Entry(self.panel_tab_1, width = 4, font = ("Times", 12),
                                     textvariable = getattr(self, "ft_mol_{}_Uiso".format(j)))
            uiso_m_entry.grid(row = 8+i, column = 8, columnspan = 2, padx = 5, sticky = "sw")
            uiso_m_entry.bind("<Return>", self.update_plot_fine_tune_para)
    
        ttk.Button(self.panel_tab_1, text = "Reset All", command = self.Reset_m_c_e_Entry,
                   style = "TButton").grid(row = 8 + len(self.mole_elements_list), column = 0, columnspan = 2, padx = 5)

        #########These are parameters for crystals###############
        ttk.Label(self.panel_tab_2, text = "Qdamp", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qbroad", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Qmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rmin", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rmax", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Rstep", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "Delta2", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 2, columnspan = 1, padx = 5, sticky = "sw")

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

        self.c_qmax_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 40.0, orient = HORIZONTAL, variable = self.qmax_2_entry, command = self.update_plot_fine_tune_para)
        self.c_qmax_scale.grid(row = 3, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_rmin_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 100.0, orient = HORIZONTAL, variable = self.rmin_entry, command = self.update_plot_fine_tune_para)
        self.c_rmin_scale.grid(row = 4, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        self.c_rmax_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 200.0, orient = HORIZONTAL, variable = self.rmax_entry, command = self.update_plot_fine_tune_para)
        self.c_rmax_scale.grid(row = 5, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_rstep_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = self.rstep_entry, command = self.update_plot_fine_tune_para)
        self.c_rstep_scale.grid(row = 6, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.c_delta2_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 10.0, orient = HORIZONTAL, variable = self.delta2_cryst_entry, command = self.update_plot_fine_tune_para)
        self.c_delta2_scale.grid(row = 7, column = 3, columnspan = 4, padx = 5, sticky = "sw")

        self.qdamp_2_entry.set(float(self.qdamp_2.get()))
        self.qbroad_2_entry.set(float(self.qbroad_2.get()))
        self.qmin_2_entry.set(float(self.qmin_2.get()))
        self.qmax_2_entry.set(float(self.qmax_2.get()))
#        self.rmin_2_entry.set(0.0)
#        self.rmax_2_entry.set(20.0)
#        self.rstep_2_entry.set(0.01)
        self.delta2_cryst_entry.set(float(self.delta2_cryst.get()))

        ttk.Label(self.panel_tab_2, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "0.5", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "2.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "40.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "100.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "200.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "1.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_2, text = "10.0", justify = CENTER,
        font = ("Times", 16, "bold")).grid(row = 7, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        

        ec1 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qdamp_2_entry)
        ec1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec1.bind("<Return>", self.update_plot_fine_tune_para)
        ec2 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qbroad_2_entry)
        ec2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec2.bind("<Return>", self.update_plot_fine_tune_para)
        ec3 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qmin_2_entry)
        ec3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec3.bind("<Return>", self.update_plot_fine_tune_para)
        ec4 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.qmax_2_entry)
        ec4.grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec4.bind("<Return>", self.update_plot_fine_tune_para)
        ec5 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.rmin_entry)
        ec5.grid(row = 4, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec5.bind("<Return>", self.update_plot_fine_tune_para)
        ec6 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.rmax_entry)
        ec6.grid(row = 5, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec6.bind("<Return>", self.update_plot_fine_tune_para)
        ec7 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.rstep_entry)
        ec7.grid(row = 6, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec7.bind("<Return>", self.update_plot_fine_tune_para)
        ec8 = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12),textvariable = self.delta2_cryst_entry)
        ec8.grid(row = 7, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ec8.bind("<Return>", self.update_plot_fine_tune_para)

        ##add fine tune options for Uisos for crystals
        for i, j in enumerate(self.cryst_elements_list):

            setattr(self, "ft_cryst_{}_Uiso".format(j), StringVar())
            getattr(self, "ft_cryst_{}_Uiso".format(j)).set(float(getattr(self, "cryst_{}_Uiso".format(j)).get()))
                    
            ttk.Label(self.panel_tab_2, text = self.cryst_elements_list[i] + "_Uiso", justify = CENTER,
                              font = ("Times", 16, "bold")).grid(row = 8+i, column = 0, columnspan = 2, padx = 5, sticky = "sw")
            ttk.Label(self.panel_tab_2, text = "0.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 2, columnspan = 1, padx = 5, sticky = "sw")
            self.Uiso_c_scale = ttk.Scale(self.panel_tab_2, from_ = 0, to = 1.0, orient = HORIZONTAL, variable = getattr(self, "ft_cryst_{}_Uiso".format(j)), command = self.update_plot_fine_tune_para)
            self.Uiso_c_scale.grid(row = 8+i, column = 3, columnspan = 4, padx = 5, sticky = "sw")
                              
            ttk.Label(self.panel_tab_2, text = "1.0", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 8+i, column = 7, columnspan = 1, padx = 5, sticky = "sw")
            uiso_c_entry = ttk.Entry(self.panel_tab_2, width = 4, font = ("Times", 12), textvariable = getattr(self, "ft_cryst_{}_Uiso".format(j)))
            uiso_c_entry.grid(row = 8+i, column = 8, columnspan = 2, padx = 5, sticky = "sw")
            uiso_c_entry.bind("<Return>", self.update_plot_fine_tune_para)

        ttk.Button(self.panel_tab_2, text = "Reset All", command = self.Reset_m_c_e_Entry,
           style = "TButton").grid(row = 8 + len(self.cryst_elements_list), column = 0, columnspan = 2, padx = 5)
        
        ##########there are the paramters for tuning experimental PDFs#############
        
        ttk.Label(self.panel_tab_3, text = "Y_scale_exp", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "0.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        self.y_scale = ttk.Scale(self.panel_tab_3, from_ = 0, to = 10, orient = HORIZONTAL, variable = self.Y_scale_entry,
                                       command = self.update_plot_fine_tune_para)
        self.y_scale.grid(row = 0, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "10.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 7, columnspan = 1, padx = 5, sticky = "sw")
                  
        ee1 = ttk.Entry(self.panel_tab_3, width = 4, font = ("Times", 12),textvariable = self.Y_scale_entry)
        ee1.grid(row = 0, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ee1.bind("<Return>", self.update_plot_fine_tune_para)
        

        ttk.Label(self.panel_tab_3, text = "Y_scale_mol_1", justify = CENTER,
                    font = ("Times", 16, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "0.0", justify = CENTER,
                    font = ("Times", 16, "bold")).grid(row = 1, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        self.y_scale = ttk.Scale(self.panel_tab_3, from_ = 0, to = 10, orient = HORIZONTAL, variable = self.mol_1_Y_scale_entry,
                    command = self.update_plot_fine_tune_para)
        self.y_scale.grid(row = 1, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "10.0", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")

        ee2 = ttk.Entry(self.panel_tab_3, width = 4, font = ("Times", 12),textvariable = self.mol_1_Y_scale_entry)
        ee2.grid(row = 1, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ee2.bind("<Return>", self.update_plot_fine_tune_para)

        ttk.Label(self.panel_tab_3, text = "Y_scale_mol_2", justify = CENTER,
               font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "0.0", justify = CENTER,
               font = ("Times", 16, "bold")).grid(row = 2, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        self.y_scale = ttk.Scale(self.panel_tab_3, from_ = 0, to = 10, orient = HORIZONTAL, variable = self.mol_2_Y_scale_entry,
                              command = self.update_plot_fine_tune_para)
        self.y_scale.grid(row = 2, column = 3, columnspan = 4, padx = 5, sticky = "sw")
        ttk.Label(self.panel_tab_3, text = "10.0", justify = CENTER,
               font = ("Times", 16, "bold")).grid(row = 2, column = 7, columnspan = 1, padx = 5, sticky = "sw")

        ee3 = ttk.Entry(self.panel_tab_3, width = 4, font = ("Times", 12),textvariable = self.mol_2_Y_scale_entry)
        ee3.grid(row = 2, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ee3.bind("<Return>", self.update_plot_fine_tune_para)


        ttk.Button(self.panel_tab_3, text = "Reset All", command = self.Reset_m_c_e_Entry, style = "TButton").grid(row = 3, column = 0,
                   columnspan = 2, padx = 5)
    
    
    def update_plot_fine_tune_para(self, event):
    
        self.fig.clf()
        ###first the molecular PDF ###
        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "ft_mol_{}_Uiso".format(j)).get())
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].Uisoequiv = float(getattr(self, "ft_cryst_{}_Uiso".format(j)).get())

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

    def Reset_m_c_e_Entry(self):

        self.fig.clf()

        for j in self.mole_elements_list:
            getattr(self, "ft_mol_{}_Uiso".format(j)).set(0.005)
        
        for j in self.cryst_elements_list:
            getattr(self, "ft_cryst_{}_Uiso".format(j)).set(0.005)
        
        for j in self.mole_elements_list:
            self.mol_struc[self.mol_struc.element == j].Uisoequiv = float(getattr(self, "ft_mol_{}_Uiso".format(j)).get())
        for j in self.cryst_elements_list:
            self.cryst_struc[self.cryst_struc.element == j].Uisoequiv = float(getattr(self, "ft_cryst_{}_Uiso".format(j)).get())

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
        
        self.Y_scale_entry.set(1.0)
        self.mol_1_Y_scale_entry.set(1.0)
        self.mol_2_Y_scale_entry.set(1.0)
        
        
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

        self.ax1 = self.fig.add_subplot(211)
        
        self.ax1.set_title("Individual PDFs")
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.set_title("Difference PDFs")
        self.ax2.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.xaxis.set_tick_params(labelsize=11)
        self.ax2.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()
        self.canvas.draw()

    def save_the_plot(self):
        self.out_file_name = tkFileDialog.asksaveasfilename(defaultextension = ".txt", filetypes = [("Text Documents", "*.txt")])
        
        s_m1_ = float(self.mol_1_Y_scale_entry.get())
        s_m2_ = float(self.mol_2_Y_scale_entry.get())
        s_e_ = float(self.Y_scale_entry.get())
        
        
        ##Scenario 1
        if self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator":
            numpy.savetxt(self.out_file_name,
            zip(self.r1,
                self.g1*s_m1_,
                self.g2*s_m2_,
                self.g_exp*s_e_,
                self.g_exp*s_e_ - self.g1*s_m1_ - self.g2*s_m2_),
            header =
                  "***********************************************************************************\n"
                  "***********************************************************************************\n"
                  "******Here are the raw data for the plotting the curves.***************************\n"
                  "******You have chosen to use PDFCalculator for both molecules.*********************\n"
                  "******You also have scaled the molecule 1 PDF by a factor of {0:.3f} **************\n"
                  "******You also have scaled the measured 2 PDF by a factor of {1:.3f} **************\n"
                  "******You also have scaled the measured PDF by a factor of {2:.3f} ****************\n"
                  "******From left to right, the data correspond to (1)r (2)m1 (3)m2 (4)e (5)e-m1-m2* \n"
                  "***********************************************************************************\n"
                  "***********************************************************************************\n".format(s_m1_, s_m2_, s_e_))


        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "PDFCalculator":
            numpy.savetxt(self.out_file_name,
              zip(self.r1,
                  self.g1*s_m1_,
                  self.g2*s_m2_,
                  self.g_exp*s_e_,
                  self.g_exp*s_e_ - self.g1*s_m1_ - self.g2*s_m2_),
              header =
              "***********************************************************************************\n"
              "***********************************************************************************\n"
              "******Here are the raw data for the plotting the curves.***************************\n"
              "******You have chosen to use DPC for mole 1; PC for mole 2.************************\n"
              "******You also have scaled the molecule 1 PDF by a factor of {0:.3f} **************\n"
              "******You also have scaled the measured 2 PDF by a factor of {1:.3f} **************\n"
              "******You also have scaled the measured PDF by a factor of {2:.3f} ****************\n"
              "******From left to right, the data correspond to (1)r (2)m1 (3)m2 (4)e (5)e-m1-m2* \n"
              "***********************************************************************************\n"
              "***********************************************************************************\n".format(s_m1_, s_m2_, s_e_))
        
        elif self.mol_pc_dpc.get() == "PDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator":
            numpy.savetxt(self.out_file_name,
              zip(self.r1,
                  self.g1*s_m1_,
                  self.g2*s_m2_,
                  self.g_exp*s_e_,
                  self.g_exp*s_e_ - self.g1*s_m1_ - self.g2*s_m2_),
              header =
              "***********************************************************************************\n"
              "***********************************************************************************\n"
              "******Here are the raw data for the plotting the curves.***************************\n"
              "******You have chosen to use PC for mole 1; DPC for mole 2.************************\n"
              "******You also have scaled the molecule 1 PDF by a factor of {0:.3f} **************\n"
              "******You also have scaled the measured 2 PDF by a factor of {1:.3f} **************\n"
              "******You also have scaled the measured PDF by a factor of {2:.3f} ****************\n"
              "******From left to right, the data correspond to (1)r (2)m1 (3)m2 (4)e (5)e-m1-m2* \n"
              "***********************************************************************************\n"
              "***********************************************************************************\n".format(s_m1_, s_m2_, s_e_))
        
        elif self.mol_pc_dpc.get() == "DebyePDFCalculator" and self.cryst_pc_dpc.get() == "DebyePDFCalculator":
            numpy.savetxt(self.out_file_name,
            zip(self.r1,
                self.g1*s_m1_,
                self.g2*s_m2_,
                self.g_exp*s_e_,
                self.g_exp*s_e_ - self.g1*s_m1_ - self.g2*s_m2_),
            header =
            "***********************************************************************************\n"
            "***********************************************************************************\n"
            "******Here are the raw data for the plotting the curves.***************************\n"
            "******You have chosen to use DebyePDFCalculator for both molecules.****************\n"
            "******You also have scaled the molecule 1 PDF by a factor of {0:.3f} **************\n"
            "******You also have scaled the measured 2 PDF by a factor of {1:.3f} **************\n"
            "******You also have scaled the measured PDF by a factor of {2:.3f} ****************\n"
            "******From left to right, the data correspond to (1)r (2)m1 (3)m2 (4)e (5)e-m1-m2* \n"
            "***********************************************************************************\n"
            "***********************************************************************************\n".format(s_m1_, s_m2_, s_e_))


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
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')

        ttk.Label(self.master, text = "To perform a fit to crystalline organic PDF, one needs to prepare a structure file for the single molecule.\n"
            "In addition, a crystal structure for compound is needed as input.", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()
        
        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Molecule Structure", command = self.load_molecule,
                             style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Crystal Structure",command = self.load_crystal,
                             style = "TButton").grid(row = 1, column = 3, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Load PDF data",command = self.load_pdf_data,
                                        style = "TButton").grid(row = 1, column = 6, columnspan = 3, padx = 5)

    
        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")

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
        self.rstep  = StringVar()
        self.optimizers = StringVar()
        
        self.lat_a_cons = StringVar()
        self.lat_b_cons = StringVar()
        self.lat_c_cons = StringVar()
        self.alpha_cons = StringVar()
        self.beta_cons = StringVar()
        self.gamma_cons = StringVar()

        ##set default values
        self.Qmin.set(0.0)
        self.zoom_mol.set(1.0)
        self.zoom_intra.set(1.0)
        self.delta2_cryst.set(0.0)
        self.delta2_mol.set(0.0)
        self.delta2_intra.set(0.0)
        self.rmin.set(1.0)
        self.rmax.set(20.0)
        self.rstep.set(0.01)
        
        ##set other default values just for test, remove later
        self.Qdamp.set(0.03)
        self.Qbroad.set(0.00)
        self.Uiso_intra.set(0.003)
        self.Uiso_inter.set(0.02)
        self.Qmax.set(24.0)  ##this Qmax is set to be same for both molecules

        self.scale.set(1)
        
        #these are the checkbuttons
        self.var_0  = IntVar()
        self.var_1  = IntVar()
        self.var_2  = IntVar()
        self.var_3  = IntVar()
        #self.var_4  = IntVar()
        #self.var_5  = IntVar()
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
        #self.var_4.set(1)
        #self.var_5.set(1)
        self.var_13.set(1)
        self.var_14.set(1)
        self.var_15.set(1)
        self.var_16.set(1)
        self.var_17.set(1)

        #################first layer Qdamp, Qbroad, Uiso_inter, Uiso_intra, Qmin, Qmax##########################
        #Qdamp
        ttk.Label(self.top_frame, text = "Qdamp", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.Qdamp).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_0,
                        style = "TCheckbutton").grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #Qbroad
        ttk.Label(self.top_frame, text = "Qbroad", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                     textvariable = self.Qbroad).grid(row = 3, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_1,
                        style = "TCheckbutton").grid(row = 3, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso intra
        ttk.Label(self.top_frame, text = "Uiso_intra", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.Uiso_intra).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_2,
                        style = "TCheckbutton").grid(row = 3, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso inter
        ttk.Label(self.top_frame, text = "Uiso_inter", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.Uiso_inter).grid(row = 3, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_3,
                        style = "TCheckbutton").grid(row = 3, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #Qmin
        ttk.Label(self.top_frame, text = "Qmin", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.Qmin).grid(row = 3, column = 13, columnspan = 1, padx = 5, sticky = "sw")
            #ttk.Checkbutton(self.top_frame, variable = self.var_4,
            #           style = "TCheckbutton").grid(row = 3, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #Qmax
        ttk.Label(self.top_frame, text = "Qmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.Qmax).grid(row = 3, column = 16, columnspan = 1, padx = 5, sticky = "sw")
            # ttk.Checkbutton(self.top_frame, variable = self.var_5,
            #           style = "TCheckbutton").grid(row = 3, column = 17, columnspan = 1, padx = 5, sticky = "sw")

        #################second layer lat_a, lat_b. lat_c, alpha, beta, gamma##########################        
        #lat_a
        ttk.Label(self.top_frame, text = "Lat_a", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.lat_a).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_6,
                        style = "TCheckbutton").grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #lat_b
        ttk.Label(self.top_frame, text = "Lat_b", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                     textvariable = self.lat_b).grid(row = 4, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_7,
                        style = "TCheckbutton").grid(row = 4, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #lat_c
        ttk.Label(self.top_frame, text = "Lat_c", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 4, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.lat_c).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_8,
                        style = "TCheckbutton").grid(row = 4, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #alpha
        ttk.Label(self.top_frame, text = "Alpha", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 4, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.alpha).grid(row = 4, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_9,
                        style = "TCheckbutton").grid(row = 4, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #beta
        ttk.Label(self.top_frame, text = "Beta", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 4, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.beta).grid(row = 4, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_10,
                        style = "TCheckbutton").grid(row = 4, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #gamma
        ttk.Label(self.top_frame, text = "Gamma", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 4, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.gamma).grid(row = 4, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_11,
                        style = "TCheckbutton").grid(row = 4, column = 17, columnspan = 1, padx = 5, sticky = "sw")
                        
        ###################new layer of constraints for a, b, c, alpha, beta, gamma################################################
        ttk.Label(self.top_frame, text = "Constraints", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        #lat_a
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                textvariable = self.lat_a_cons).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        #lat_b
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                textvariable = self.lat_b_cons).grid(row = 5, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        #lat_c
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                textvariable = self.lat_c_cons).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        #alpha
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                textvariable = self.alpha_cons).grid(row = 5, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        #beta
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                textvariable = self.beta_cons).grid(row = 5, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        #gamma
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                textvariable = self.gamma_cons).grid(row = 5, column = 16, columnspan = 1, padx = 5, sticky = "sw")

        #################third layer scale, zoom_mol, zoom_intra, delta2_cryst, delta2_mole, delta2_intra##########################        
        #scale
        ttk.Label(self.top_frame, text = "Scale", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 6, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.scale).grid(row = 6, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_12,
                        style = "TCheckbutton").grid(row = 6, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #zoom_mol
        ttk.Label(self.top_frame, text = "Zoom_Mol", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 6, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Times", 12),
                     textvariable = self.zoom_mol).grid(row = 6, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_13,
                        style = "TCheckbutton").grid(row = 6, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #zoom_intra
        ttk.Label(self.top_frame, text = "Zoom_Intra", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 6, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.zoom_intra).grid(row = 6, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_14,
                        style = "TCheckbutton").grid(row = 6, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_cryst
        ttk.Label(self.top_frame, text = "Delta2_Cryst", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 6, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.delta2_cryst).grid(row = 6, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_15,
                        style = "TCheckbutton").grid(row = 6, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_mol
        ttk.Label(self.top_frame, text = "Delta2_Mol", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 6, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.delta2_mol).grid(row = 6, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_16,
                        style = "TCheckbutton").grid(row = 6, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_intra
        ttk.Label(self.top_frame, text = "Delta2_Intra", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 6, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                     textvariable = self.delta2_intra).grid(row = 6, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_17,
                        style = "TCheckbutton").grid(row = 6, column = 17, columnspan = 1, padx = 5, sticky = "sw")

        ## layout for step 3, Select Optimizers and run

        ttk.Label(self.top_frame, text = "Step 3: Select Optimzers", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        #rmin
        ttk.Label(self.top_frame, text = "Rmin", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 8, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                      textvariable = self.rmin).grid(row = 8, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                  
        #rmax
        ttk.Label(self.top_frame, text = "Rmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 8, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                      textvariable = self.rmax).grid(row = 8, column = 3, columnspan = 1, padx = 5, sticky = "sw")

        #rstep
        ttk.Label(self.top_frame, text = "Rstep", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 8, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                  textvariable = self.rstep).grid(row = 8, column = 5, columnspan = 1, padx = 5, sticky = "sw")


        self.combobox = ttk.Combobox(self.top_frame, font = ("Times", 12), textvariable = self.optimizers, values = ["Least Squares", "Fmin",
                              "Basin Hopping"])
        self.combobox.current(0)
        self.combobox.grid(row = 8, column = 6, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.top_frame, text = "Run The FIT", command = self.run_the_fit,
                             style = "TButton").grid(row = 8, column = 8, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Reset To Default", command = self.reset_to_default,
                             style = "TButton").grid(row = 8, column = 11, columnspan = 3, padx = 5, sticky = "sw")


    ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.bottom_frame = ttk.Frame(self.master, padding = (10, 10))
        self.bottom_frame.pack()

        self.fig = plt.figure(figsize=(12, 5), dpi=100) ##create a figure; modify the size here
        
        self.ax = self.fig.add_subplot(111)
        
        self.ax.set_title("PDF Model Fit To Crystalline Organic PDF")
        self.ax.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax.xaxis.set_tick_params(labelsize=11)
        self.ax.yaxis.set_tick_params(labelsize=11)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.bottom_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.bottom_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    
    def load_molecule(self):
        try:
            self.input_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz",
                                                      filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.input_mole_struct)

        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")


    def load_crystal(self):
        try:
            self.input_cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".cif",
                                                        filetypes = [("Text Documents", "*.cif")])
            Structure(filename = self.input_cryst_struct)

            self.lat_a.set(loadStructure(self.input_cryst_struct).lattice.a)
            self.lat_b.set(loadStructure(self.input_cryst_struct).lattice.b)
            self.lat_c.set(loadStructure(self.input_cryst_struct).lattice.c)
            self.alpha.set(loadStructure(self.input_cryst_struct).lattice.alpha)
            self.beta.set(loadStructure(self.input_cryst_struct).lattice.beta)
            self.gamma.set(loadStructure(self.input_cryst_struct).lattice.gamma)
            
            a = loadStructure(self.input_cryst_struct).lattice.a
            b = loadStructure(self.input_cryst_struct).lattice.b
            c = loadStructure(self.input_cryst_struct).lattice.c
            alpha = loadStructure(self.input_cryst_struct).lattice.alpha
            beta = loadStructure(self.input_cryst_struct).lattice.beta
            gamma = loadStructure(self.input_cryst_struct).lattice.gamma
        
        ##first lattice parameters
            if a==b==c:
                self.lat_a_cons.set("@1")
                self.lat_b_cons.set("@1")
                self.lat_c_cons.set("@1")
            elif a == b != c:
                self.lat_a_cons.set("@1")
                self.lat_b_cons.set("@1")
                self.lat_c_cons.set("@3")
            elif a == c != b:
                self.lat_a_cons.set("@1")
                self.lat_b_cons.set("@2")
                self.lat_c_cons.set("@1")
            elif b == c != a:
                self.lat_a_cons.set("@1")
                self.lat_b_cons.set("@2")
                self.lat_c_cons.set("@2")
            elif a != b != c:
                self.lat_a_cons.set("@1")
                self.lat_b_cons.set("@2")
                self.lat_c_cons.set("@3")
        ##then angles
            if alpha == beta == gamma:
                self.alpha_cons.set("@4")
                self.beta_cons.set("@4")
                self.gamma_cons.set("@4")
            elif alpha == beta != gamma:
                self.alpha_cons.set("@4")
                self.beta_cons.set("@4")
                self.gamma_cons.set("@6")
            elif alpha == gamma != beta:
                self.alpha_cons.set("@4")
                self.beta_cons.set("@5")
                self.gamma_cons.set("@4")
            elif beta == gamma!=alpha:
                self.alpha_cons.set("@4")
                self.beta_cons.set("@5")
                self.gamma_cons.set("@5")
            elif alpha != beta != gamma:
                self.alpha_cons.set("@4")
                self.beta_cons.set("@5")
                self.gamma_cons.set("@6")
        except:

            tkMessageBox.showwarning("Warning!", "The crystal structure file cannot be parsed properly by the program! Please check...")


    def load_pdf_data(self):
        try:
            self.input_pdf_gr = tkFileDialog.askopenfilename(defaultextension = ".gr", filetypes = [("Text Documents", "*.gr")])
            parser = PDFParser()
            parser.parseFile(self.input_pdf_gr)
                
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")
    
    def reset_to_default(self):
        ##checkbox
        self.var_0.set(1)
        self.var_1.set(1)
        #self.var_4.set(1)
        #self.var_5.set(1)
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
        self.rstep.set(0.01)
        self.combobox.current(0)
        
        ##set other default values just for test, remove later
        self.Qdamp.set(0.03)
        self.Qbroad.set(0.00)
        self.Uiso_intra.set(0.003)
        self.Uiso_inter.set(0.02)
        self.Qmax.set(24.0)  ##this Qmax is set to be same for both molecules
        
        self.scale.set(1)

        self.lat_a.set(loadStructure(self.input_cryst_struct).lattice.a)
        self.lat_b.set(loadStructure(self.input_cryst_struct).lattice.b)
        self.lat_c.set(loadStructure(self.input_cryst_struct).lattice.c)
        self.alpha.set(loadStructure(self.input_cryst_struct).lattice.alpha)
        self.beta.set(loadStructure(self.input_cryst_struct).lattice.beta)
        self.gamma.set(loadStructure(self.input_cryst_struct).lattice.gamma)
    

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
            profile.setCalculationRange(xmin= float(self.rmin.get()), xmax = float(self.rmax.get()), dx = float(self.rstep.get()))
            
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


            a = self.lat_a_cons.get()
            b = self.lat_b_cons.get()
            c = self.lat_c_cons.get()
            ##first lattice parameters

            if int(self.var_6.get()) == 0 and int(self.var_7.get()) == 0 and int(self.var_8.get()) == 0:
                if a==b==c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_a")
                    recipe.constrain(lat.c, "lat_a")
                elif a == b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_c", float(self.lat_c.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_a")
                    recipe.constrain(lat.c, "lat_c")
                elif a == c != b:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_b")
                    recipe.constrain(lat.c, "lat_a")
                elif b == c != a:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_b")
                    recipe.constrain(lat.c, "lat_b")
                elif a != b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.newVar("lat_c", float(self.lat_c.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_b")
                    recipe.constrain(lat.c, "lat_c")

            elif int(self.var_6.get()) == 0 and int(self.var_7.get()) == 1 and int(self.var_8.get()) == 0:
                lat.b = float(self.lat_b.get())
                if a==b==c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.c, "lat_a")
                elif a == b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_c", float(self.lat_c.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.c, "lat_c")
                elif a == c != b:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.c, "lat_a")
                elif b == c != a:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.c, "lat_b")
                elif a != b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_c", float(self.lat_c.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.c, "lat_c")

            elif int(self.var_6.get()) == 0 and int(self.var_7.get()) == 0 and int(self.var_8.get()) == 1:
                lat.c = float(self.lat_c.get())
            
                if a==b==c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_a")
                elif a == b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_a")
                elif a == c != b:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_b")
                elif b == c != a:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_b")
                elif a != b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.a, "lat_a")
                    recipe.constrain(lat.b, "lat_b")

            elif int(self.var_6.get()) == 1 and int(self.var_7.get()) == 0 and int(self.var_8.get()) == 0:
                lat.a = float(self.lat_a.get())
                
                if a==b==c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.constrain(lat.b, "lat_a")
                    recipe.constrain(lat.c, "lat_a")
                elif a == b != c:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_c", float(self.lat_c.get()))
                    recipe.constrain(lat.b, "lat_a")
                    recipe.constrain(lat.c, "lat_c")
                elif a == c != b:
                    recipe.newVar("lat_a", float(self.lat_a.get()))
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.b, "lat_b")
                    recipe.constrain(lat.c, "lat_a")
                elif b == c != a:
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.constrain(lat.b, "lat_b")
                    recipe.constrain(lat.c, "lat_b")
                elif a != b != c:
                    recipe.newVar("lat_b", float(self.lat_b.get()))
                    recipe.newVar("lat_c", float(self.lat_c.get()))
                    recipe.constrain(lat.b, "lat_b")
                    recipe.constrain(lat.c, "lat_c")

            elif int(self.var_6.get()) == 0 and int(self.var_7.get()) == 1 and int(self.var_8.get()) == 1:
                lat.b = float(self.lat_b.get())
                lat.c = float(self.lat_c.get())
                recipe.newVar("lat_a", float(self.lat_a.get()))
                recipe.constrain(lat.a, "lat_a")

            elif int(self.var_6.get()) == 1 and int(self.var_7.get()) == 0 and int(self.var_8.get()) == 1:
                lat.a = float(self.lat_a.get())
                lat.c = float(self.lat_c.get())
                recipe.newVar("lat_b", float(self.lat_b.get()))
                recipe.constrain(lat.b, "lat_b")

            elif int(self.var_6.get()) == 1 and int(self.var_7.get()) == 1 and int(self.var_8.get()) == 0:
                lat.a = float(self.lat_a.get())
                lat.b = float(self.lat_b.get())
                recipe.newVar("lat_c", float(self.lat_c.get()))
                recipe.constrain(lat.c, "lat_c")

            elif int(self.var_6.get()) == 1 and int(self.var_7.get()) == 1 and int(self.var_8.get()) == 1:
                lat.a = float(self.lat_a.get())
                lat.b = float(self.lat_b.get())
                lat.c = float(self.lat_c.get())

            ###second, alpha beta gamma,

            alpha = self.alpha_cons.get()
            beta  = self.beta_cons.get()
            gamma = self.gamma_cons.get()

            if int(self.var_9.get()) == 0 and int(self.var_10.get()) == 0 and int(self.var_11.get()) == 0:
                if alpha==beta==gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "alpha")
                    recipe.constrain(lat.gamma, "alpha")
                elif alpha == beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("gamma", float(self.gamma.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "alpha")
                    recipe.constrain(lat.gamma, "gamma")
                elif alpha == gamma != beta:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "beta")
                    recipe.constrain(lat.gamma, "alpha")
                elif beta == gamma != alpha:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "beta")
                    recipe.constrain(lat.gamma, "beta")
                elif alpha != beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.newVar("gamma", float(self.gamma.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "beta")
                    recipe.constrain(lat.gamma, "gamma")
            
            elif int(self.var_9.get()) == 0 and int(self.var_10.get()) == 1 and int(self.var_11.get()) == 0:
                lat.b = float(self.lat_b.get())
                if alpha==beta==gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.gamma, "alpha")
                elif alpha == beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("gamma", float(self.gamma.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.gamma, "gamma")
                elif alpha == gamma != beta:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.gamma, "alpha")
                elif beta == gamma != alpha:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.gamma, "beta")
                elif alpha != beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("gamma", float(self.gamma.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.gamma, "gamma")

            elif int(self.var_9.get()) == 0 and int(self.var_10.get()) == 0 and int(self.var_11.get()) == 1:
                lat.c = float(self.lat_c.get())
        
                if alpha==beta==gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "alpha")
                elif alpha == beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "alpha")
                elif alpha == gamma != beta:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "beta")
                elif beta == gamma != alpha:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "beta")
                elif alpha != beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.alpha, "alpha")
                    recipe.constrain(lat.beta, "beta")
            
            elif int(self.var_9.get()) == 1 and int(self.var_10.get()) == 0 and int(self.var_11.get()) == 0:
                lat.a = float(self.lat_a.get())
                
                if alpha==beta==gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.constrain(lat.beta, "alpha")
                    recipe.constrain(lat.gamma, "alpha")
                elif alpha == beta != gamma:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("gamma", float(self.gamma.get()))
                    recipe.constrain(lat.beta, "alpha")
                    recipe.constrain(lat.gamma, "gamma")
                elif alpha == gamma != beta:
                    recipe.newVar("alpha", float(self.alpha.get()))
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.beta, "beta")
                    recipe.constrain(lat.gamma, "alpha")
                elif beta == gamma != alpha:
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.constrain(lat.beta, "beta")
                    recipe.constrain(lat.gamma, "beta")
                elif alpha != beta != gamma:
                    recipe.newVar("beta", float(self.beta.get()))
                    recipe.newVar("gamma", float(self.gamma.get()))
                    recipe.constrain(lat.beta, "beta")
                    recipe.constrain(lat.gamma, "gamma")

            elif int(self.var_9.get()) == 0 and int(self.var_10.get()) == 1 and int(self.var_11.get()) == 1:
                lat.beta = float(self.beta.get())
                lat.gamma = float(self.gamma.get())
                recipe.newVar("alpha", float(self.alpha.get()))
                recipe.constrain(lat.a, "alpha")
            
            elif int(self.var_9.get()) == 1 and int(self.var_10.get()) == 0 and int(self.var_11.get()) == 1:
                lat.alpha = float(self.alpha.get())
                lat.gamma = float(self.gamma.get())
                recipe.newVar("beta", float(self.beta.get()))
                recipe.constrain(lat.b, "beta")
            
            elif int(self.var_9.get()) == 1 and int(self.var_10.get()) == 1 and int(self.var_11.get()) == 0:
                lat.alpha = float(self.alpha.get())
                lat.beta = float(self.beta.get())
                recipe.newVar("gamma", float(self.gamma.get()))
                recipe.constrain(lat.gamma, "gamma")
            
            elif int(self.var_9.get()) == 1 and int(self.var_10.get()) == 1 and int(self.var_11.get()) == 1:
                lat.alpha = float(self.alpha.get())
                lat.beta = float(self.beta.get())
                lat.gamma = float(self.gamma.get())

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
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(r,g,'bo',markersize = 5, label="G(r) Data")
            self.ax.plot(r, gcalc,'r-',label="G(r) Fit", lw= 2)
            self.ax.plot(r,diff,'g-',label="G(r) diff", lw = 2)
            self.ax.plot(r,diffzero,'k-', lw = 2)
            self.ax.xaxis.set_tick_params(labelsize=11)
            self.ax.yaxis.set_tick_params(labelsize=11)
            self.ax.legend(loc=1)
            self.ax.set_title("PDF Model Fit To Crystalline Organic PDF")
            self.ax.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
            self.ax.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
            self.canvas.draw()

            return

        data  = self.input_pdf_gr
        stru1 = Structure(filename = self.input_cryst_struct)
        stru2 = Structure(filename = self.input_mole_struct)
        stru3 = Structure(filename = self.input_mole_struct)
        recipe = makeRecipe(stru1, stru2, stru3, data)
        
        recipe.fithooks[0].verbose = 3

        ###different optimizers###
        from scipy.optimize import leastsq
        from scipy.optimize import fmin
        from scipy.optimize import basinhopping
        
        if self.optimizers.get() == "Least Squares":
            leastsq(recipe.residual, recipe.values)
        
        elif self.optimizers.get() == "Fmin":
            fmin(recipe.scalarResidual, recipe.getValues())

        elif self.optimizers.get() == "Basin Hopping":
            basinhopping(recipe.scalarResidual, recipe.getValues())

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


class Break_down_fit():
    def __init__(self, master):
        self.master = master
        self.master.title("Breakdown the intra and intermolecular parts of the total fit")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 600) # width + height
        self.master.resizable(False, False)
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        
        ttk.Label(self.master, text = "Based on the fit result, calculate intra- and inter-part of the fit.", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()

        self.qmax = StringVar()
        self.qmin = StringVar()
        
        self.qmax.set(24.0)
        self.qmin.set(0.0)

        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1: Load Structures and Fits", justify = CENTER,
                font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 3, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load Intra Structure", command = self.load_intra_molecule,
                 style = "TButton").grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Mole Structure", command = self.load_inter_molecule,
                            style = "TButton").grid(row = 1, column = 2, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Load Crystal Structure",command = self.load_crystal,
                 style = "TButton").grid(row = 1, column = 4, columnspan = 2, padx = 5)
        ttk.Button(self.top_frame, text = "Load Fit Data",command = self.load_pdf_fit,
                 style = "TButton").grid(row = 1, column = 6 , columnspan = 2, padx = 5)
        ttk.Button(self.top_frame, text = "Load Fit Result",command = self.load_fit_result,
                            style = "TButton").grid(row = 1, column = 8, columnspan = 3, padx = 5)

        ttk.Label(self.top_frame, text = "Step 2: Make and Save the Plot", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        
        ttk.Label(self.top_frame, text = "Qmin", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                            textvariable = self.qmin).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")
                            
        ttk.Label(self.top_frame, text = "Qmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                  textvariable = self.qmax).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        
        ttk.Button(self.top_frame, text = "Make The Plot", command = self.make_the_plot,
                   style = "TButton").grid(row = 3, column = 4, columnspan = 2, padx = 5)
        ttk.Button(self.top_frame, text = "Clear The Plot", command = self.clear_the_plot,
                       style = "TButton").grid(row = 3, column = 6, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Save The Plot", command = self.save_the_plot,
                                  style = "TButton").grid(row = 3, column = 8, columnspan = 3, padx = 5, sticky = "sw")

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.bottom_frame = ttk.Frame(self.master, padding = (10, 10))
        self.bottom_frame.pack()
        
        self.fig = plt.figure(figsize=(9, 6), dpi=100) ##create a figure; modify the size here
        
        self.ax = self.fig.add_subplot(111)
        
        self.ax.set_title("Breakdown of PDF Fit")
        self.ax.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax.xaxis.set_tick_params(labelsize=11)
        self.ax.yaxis.set_tick_params(labelsize=11)
        
        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.bottom_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.bottom_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def load_intra_molecule(self):
        try:
            self.intra_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz", filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.intra_mole_struct)
        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")

    def load_inter_molecule(self):
        try:
            self.inter_mole_struct = tkFileDialog.askopenfilename(defaultextension = ".xyz", filetypes = [("Text Documents", "*.xyz")])
            Structure(filename = self.inter_mole_struct)
        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")

    def load_crystal(self):
        try:
            self.cryst_struct = tkFileDialog.askopenfilename(defaultextension = ".cif", filetypes = [("Text Documents", "*.cif"), ("Text Documents", "*.stru")])
            Structure(filename = self.cryst_struct)
        
        except:
            tkMessageBox.showwarning("Warning!", "The molecule structure file cannot be parsed properly by the program! Please check...")

    def load_pdf_fit(self):
        try:
            self.pdf_fit = tkFileDialog.askopenfilename(defaultextension = ".fit", filetypes = [("Text Documents", "*.fit")])
            numpy.loadtxt(self.pdf_fit)
            
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")

    def load_fit_result(self):
        try:
            self.fit_res = tkFileDialog.askopenfilename(defaultextension = ".res", filetypes = [("Text Documents", "*.res")])
        
        except:
            tkMessageBox.showwarning("Warning!", "The PDF data file cannot be parsed properly by the program! Please check...")


    def make_the_plot(self):
        
        # load fit data
        fit_data = numpy.loadtxt(self.pdf_fit).transpose()
        fit_x, fit_y, gr = fit_data[0], fit_data[1], fit_data[2]
        
        # retrieve data
        
        res_file = open(self.fit_res)
        split_res_file = []
        for line in res_file:
            split_res_file.append(line.split())
        
        retrieve_params = [] # in the order of ['Global Scale','Qdamp', 'Qbroad', 'Uiso_Inter', 'Uiso_Intra', 'delta2_intra', 'delta2_mol', 'delta2_cryst']
        
        for line in split_res_file:
            if len(line) != 0:
        
                if line[0] == 'Global_Scale':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'Qdamp':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'Qbroad':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'Uiso_Inter':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'Uiso_Intra':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'delta2_intra':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'delta2_mol':
                    retrieve_params.append(line[1])
                
                elif line[0] == 'delta2_cryst':
                    retrieve_params.append(line[1])
            
        retrieve_params = [float(i) for i in retrieve_params]
        
        scale = retrieve_params[0]
        qdamp = retrieve_params[1]
        qbroad = retrieve_params[2]
        Uisointer = retrieve_params[3]
        Uisointra = retrieve_params[4]
        delta2_intra = retrieve_params[5]
        delta2_mole = retrieve_params[6]
        delta2_cryst = retrieve_params[7]
        
        ##FIXME
        
        x_array = numpy.loadtxt(self.pdf_fit).transpose()[0]
        
        rmin = min(x_array)
        rmax = max(x_array)
        rstep = x_array[1] - x_array[0]

        qmin = float(self.qmin.get())
        qmax = float(self.qmax.get())

        S = Structure(filename = self.intra_mole_struct)
        duplicate_elements = [i for i in S.element]
        element_list = list(set(duplicate_elements))
        ##
        
        cfg_intra = {   'qmax' : qmax,
                        'qmin' : qmin,
                        'rmin' : rmin,
                        'rmax' : rmax,
                        'rstep': rstep,
                        'qdamp': qdamp,
                        'qbroad': qbroad,
                        'delta2':delta2_intra
        
        }

        cfg_mole = {        'qmax' : qmax,
                            'qmin' :qmin,
                            'rmin' : rmin,
                            'rmax' : rmax,
                            'rstep': rstep,
                            'qdamp': qdamp,
                            'qbroad': qbroad,
                            'delta2': delta2_mole
    
        }
        
        cfg_cryst = {       'qmax' : qmax,
                            'qmin' : qmin,
                            'rmin' : rmin,
                            'rmax' : rmax,
                            'rstep': rstep,
                            'qdamp': qdamp,
                            'qbroad': qbroad,
                            'delta2': delta2_cryst

                    }
    
    
        # Phase 1: xyz_intra
        intra = Structure(filename = self.intra_mole_struct)
        # set ADPs for all atoms
        intra_thermal = Uisointra
        for i in element_list:
            intra[intra.element== i].Uisoequiv = intra_thermal

        dpc_intra = DebyePDFCalculator(**cfg_intra)
        self.r1, self.g1 = dpc_intra(intra)

        # Phase 2: xyz_mole
        mole = Structure(filename = self.inter_mole_struct)
            # set ADPs for all atoms
        mole_thermal = Uisointer
        for i in element_list:
            mole[mole.element== i].Uisoequiv = mole_thermal

        dpc_mole = DebyePDFCalculator(**cfg_mole)
        self.r2, self.g2 = dpc_mole(mole)

        # Phase 3: cryst
        cryst = Structure(filename = self.cryst_struct)
        # set ADPs for all atoms
        cryst_thermal = Uisointer
        for i in element_list:
            cryst[cryst.element== i].Uisoequiv = cryst_thermal

        pc_cryst = PDFCalculator(**cfg_cryst)
        self.r3, self.g3 = pc_cryst(cryst)

        self.fit_sum = scale * (self.g1 + self.g3 - self.g2)
        
        self.fig.clf()
        
        self.ax1 = self.fig.add_subplot(321)
        self.ax1.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax1.xaxis.set_tick_params(labelsize=11)
        self.ax1.yaxis.set_tick_params(labelsize=11)
        self.ax1.plot(self.r1, scale*self.g1, 'r-', lw=2, label = 'caculated IntraPDF')
        self.ax1.legend(loc=0)

        self.ax2 = self.fig.add_subplot(323)
        self.ax2.plot(self.r1, scale*self.g3, 'b-', lw=2, label = 'calculated CrystalPDF')
        self.ax2.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax2.legend(loc=0)

        self.ax3 = self.fig.add_subplot(325)
        self.ax3.plot(self.r1, scale*self.g2, 'k-', lw=2, label = 'calculated MolePDF')
        self.ax3.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax3.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)

        self.ax3.legend(loc=0)
        
        self.ax4 = self.fig.add_subplot(322)
        self.ax4.plot(self.r1, scale*(self.g3 - self.g2), 'm-', lw=2, label = 'calculated InterPDF')
        self.ax4.legend(loc=0)

        self.ax5 = self.fig.add_subplot(324)
        self.ax5.plot(fit_x, fit_y, 'c-', lw=2, label = 'total fit')
        self.ax5.plot(self.r1, self.fit_sum, 'y-', lw=2, label = 'Intra + InterPDF')
        self.ax5.legend(loc=0)

        self.ax6 = self.fig.add_subplot(326)
        self.ax6.plot(fit_x, gr, 'g-', lw=2, label = 'experimental PDF')
        self.ax6.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)

        self.ax6.legend(loc=0)
        self.fig.tight_layout()
        self.canvas.draw()

    def clear_the_plot(self):
        
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Breakdown of PDF Fit")
        self.ax.set_xlabel(ur"r (\u00c5)", labelpad = 3, fontsize = 15)
        self.ax.set_ylabel(ur"PDF, G (\u00c5$^{-2})$", labelpad = 10, fontsize = 15)
        self.ax.xaxis.set_tick_params(labelsize=11)
        self.ax.yaxis.set_tick_params(labelsize=11)
        self.fig.tight_layout()
        self.canvas.draw()


    def save_the_plot(self):
        self.out_file_name = tkFileDialog.asksaveasfilename(defaultextension = ".txt", filetypes = [("Text Documents", "*.txt")])
        numpy.savetxt(self.out_file_name,
                          zip(self.r1,
                              self.g1,
                              self.g2,
                              self.g3,
                              self.g3 - self.g2,
                              self.fit_sum),
                          header =
                          "***********************************************************************************\n"
                          "***********************************************************************************\n"
                          "******Here are the raw data for the plotting the curves.***************************\n"
                          "****Each column from left to right, corresponds to, (1)r (2)simulated G_Intra******\n"
                          "****(3)simualted G_crystal (4) Simulated G_mole (5) Simulated G_Inter**************\n"
                          "****(6)simulated G_total***********************************************************\n"
                          "***********************************************************************************\n"
                            )
class PCA():
    def __init__(self, master):
        self.master = master
        self.master.title("Load PDFs and generate score/scree plot")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(800, 600) # width + height
        self.master.resizable(False, False)
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91', font = ("Times", 16, "bold"))
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91')
        
        ttk.Label(self.master, text = "Generate Score Plot and Scree Plot based on Principle Component Analysis", font = ("Times", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (10, 10))
        self.top_frame.pack()
        
        self.legend_marker_size = StringVar()
        self.rmin = StringVar()
        self.rmax = StringVar()
        self.rstep = StringVar()
        
        self.legend_marker_size.set(8)
        self.rmin.set(1.0)
        self.rmax.set(40.0)
        self.rstep.set(0.01)
        
        ##here are the layout for step 1, load structure files
        
        ttk.Label(self.top_frame, text = "Step 1: Load a series of PDFs", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 0, column = 0, columnspan = 3, padx = 5, sticky = "sw")
        
        ttk.Button(self.top_frame, text = "Load a series of PDFs", command = self.load_pdfs,
                 style = "TButton").grid(row = 1, column = 0, columnspan = 3, padx = 5, sticky = "sw")
                 
        ttk.Label(self.top_frame, text = "Rmin", justify = CENTER,
                   font = ("Times", 16, "bold")).grid(row = 1, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                   textvariable = self.rmin).grid(row = 1, column = 4, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Label(self.top_frame, text = "Rmax", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 1, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                  textvariable = self.rmax).grid(row = 1, column = 6, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Label(self.top_frame, text = "Rstep", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 1, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                  textvariable = self.rstep).grid(row = 1, column = 8, columnspan = 1, padx = 5, sticky = "sw")

        ttk.Label(self.top_frame, text = "Legend marker size", justify = CENTER,
                  font = ("Times", 16, "bold")).grid(row = 1, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Times", 12),
                  textvariable = self.legend_marker_size).grid(row = 1, column = 10, columnspan = 1, padx = 5, sticky = "sw")


        ttk.Label(self.top_frame, text = "Step 2: Make and Save the Plot", justify = CENTER, font = ("Times", 16, "bold")).grid(row = 2, column = 0, columnspan = 3, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Make The Plot", command = self.generate_PCA_plot,
                 style = "TButton").grid(row = 3, column = 0, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Clear The Plot", command = self.clear_the_pca_plot,
                 style = "TButton").grid(row = 3, column = 3, columnspan = 3, padx = 5, sticky = "sw")
        ttk.Button(self.top_frame, text = "Save The Plot", command = self.save_pca_results,
                 style = "TButton").grid(row = 3, column = 6, columnspan = 3, padx = 5, sticky = "sw")

        ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.bottom_frame = ttk.Frame(self.master, padding = (10, 10))
        self.bottom_frame.pack()

        self.fig = plt.figure(figsize=(9, 5), dpi=100) ##create a figure; modify the size here

        self.ax1 = self.fig.add_subplot(121)
        
        self.ax1.set_xlabel('Principal Component 1',fontsize = 12)
        self.ax1.set_ylabel('Principal Component 2', fontsize = 12)
        self.ax1.set_title('Plot of first two principle components')
        
        self.ax2 = self.fig.add_subplot(122)
        self.ax2.set_ylabel('Percentage of Explained Variance',fontsize = 12)
        self.ax2.set_xlabel('Principal Component',fontsize = 12)
        self.ax2.set_title('Scree Plot')

        self.fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.bottom_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.bottom_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    def load_pdfs(self):
        try:
            self.input_pdfs = list(tkFileDialog.askopenfilenames(defaultextension = ".gr", filetypes = [("Text Documents", "*.gr")]))
            for pdf in self.input_pdfs:
                parser = PDFParser()
                parser.parseFile(pdf)
        
        except:
            tkMessageBox.showwarning("Warning!", "One of more PDF data file cannot be parsed properly by the program! Please check...")
            

        if len(self.input_pdfs) == 1:
            tkMessageBox.showwarning("Warning!", "You need more than one PDF file for generating Score/Scree plots...")


    def return_x_y_array(self, input_gr, rmin, rmax, rstep):
        profile = Profile()
        parser = PDFParser()
        parser.parseFile(input_gr)
        profile.loadParsedData(parser)
        profile.setCalculationRange(xmin = rmin, xmax = rmax, dx = rstep)
        r = profile.x
        g = profile.y
        return (r, g)

    def generate_PCA_plot(self):
        
        ##store data
        data = []
        name_list = self.input_pdfs
        
        for name in name_list:
            r, g = self.return_x_y_array(name, float(self.rmin.get()), float(self.rmax.get()), float(self.rstep.get()))
            data.append(g)
        
        data = numpy.array(data)
        
        pca = SKPCA()
        pca.fit(data)
        self.data_transform = pca.transform(data)
        
        self.variance_ratio = pca.explained_variance_ratio_

        markers_database = itertools.cycle(['bo', 'bs', 'bh', 'b^', 'bv',
                                        'ko', 'ks', 'kh', 'k^', 'kv',
                                        'ro', 'rs', 'rh', 'r^', 'rv',
                                        'mo', 'ms', 'mh', 'm^', 'mv',
                                        'go', 'gs', 'gh', 'g^', 'gv',
                                        'b<', 'b>', 'bp', 'bD', 'b*',
                                        'k<', 'k>', 'kp', 'kD', 'k*',
                                        'r<', 'r>', 'rp', 'rD', 'r*',
                                        'm<', 'm>', 'mp', 'mD', 'm*',
                                        'g<', 'g>', 'gp', 'gD', 'g*'])
        
        self.ax1 = self.fig.add_subplot(121)
        
        self.pretty_name_list = []
        for name in name_list:
            self.pretty_name_list.append(name.split('/')[-1])

        for x_point, y_point, label in zip(self.data_transform[:,0], self.data_transform[:,1], self.pretty_name_list):
            self.ax1.plot(x_point, y_point, next(markers_database), markersize = int(self.legend_marker_size.get()), label = label)

        self.ax1.set_xlabel('Principal Component 1 ({:.3f})'.format(self.variance_ratio[0]), fontsize = 12)
        self.ax1.set_ylabel('Principal Component 2 ({:.3f})'.format(self.variance_ratio[1]), fontsize = 12)
        self.ax1.axvline(x=0, linestyle = 'dotted', color = 'r', lw=2)
        self.ax1.axhline(y=0, linestyle = 'dotted', color = 'r', lw=2)
        self.ax1.legend(loc=0, prop = {'size': int(self.legend_marker_size.get())})
        self.ax1.set_title('Plot of first two principle components')

        self.ax2 = self.fig.add_subplot(122)
        ybar = self.variance_ratio*100
        xbar = ['PC' + str(x) for x in range(1, len(ybar) + 1)]
        self.ax2.bar(range(1, len(ybar)+1), ybar, tick_label = xbar)
        self.ax2.set_ylabel('Percentage of Explained Variance',fontsize = 12)
        self.ax2.set_xlabel('Principal Component',fontsize = 12)
        self.ax2.set_title('Scree Plot')
        self.fig.tight_layout()
        self.canvas.draw()

    def save_pca_results(self):
        self.out_file_name = tkFileDialog.asksaveasfilename(defaultextension = ".txt", filetypes = [("Text Documents", "*.txt")])
        filenames = numpy.array(self.pretty_name_list)
        DAT = [filenames, numpy.around(self.variance_ratio, decimals = 4)]
        num_of_pcs = len(self.variance_ratio)
        for pc_value in self.data_transform.transpose():
            DAT.append(numpy.around(pc_value, decimals = 4))
        
        data = numpy.array(DAT).transpose()
        
        total_cols = num_of_pcs + 2

        numpy.savetxt(self.out_file_name, data, delimiter=" ", fmt="%-8s",
                      header =
                      "***********************************************************************************\n"
                      "***********************************************************************************\n"
                      "**************Here are the raw data for the scores/scree plots*********************\n"
                      "************************There are in total {} columns of data**********************\n"
                      "**************The first column corresponds to each file name***********************\n"
                      "********The second column corresponds to explained variance for scree plot*********\n"
                      "********From third column to the end are {} principle component values*************\n".format(total_cols, num_of_pcs)
                      )

    def clear_the_pca_plot(self):
        
        self.fig.clf()
        
        self.ax1 = self.fig.add_subplot(121)
        self.ax1.set_xlabel('Principal Component 1',fontsize = 12)
        self.ax1.set_ylabel('Principal Component 2', fontsize = 12)
        self.ax1.set_title('Plot of first two principle components')
        self.ax1.set_title('Plot of first two principle components')
        
        self.ax2 = self.fig.add_subplot(122)
        self.ax2.set_ylabel('Percentage of Explained Variance',fontsize = 12)
        self.ax2.set_xlabel('Principal Component',fontsize = 12)
        self.ax2.set_title('Scree Plot')
        
        self.fig.tight_layout()
        self.canvas.draw()


def main():
    root = Tk()
    GUI = Overall_Look(root)
    root.mainloop()

if __name__ == "__main__": main()
