""""
This is the code for OrganicPDF GUI program. It is used to extract information from a measured Organic PDF.
The code is written by Chenyang Shi at AbbVie.
Send him an email at chenyang.shi@abbvie.com for questions and suggestions. 

"""

from Tkinter import *
import Tkinter as tk
import ttk
from math import exp
import os  # for loading files or exporting files
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure



class PDF_model_fit():
    def __init__(self, master):
        self.master = master
        self.master.title("Model Fit to Organic PDF")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(900, 650) # width + height
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
    
        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2: Set Parameters", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")

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
        
        self.var_0  = StringVar()
        self.var_1 = StringVar()
        self.var_2 = StringVar()
        self.var_3 = StringVar()
        self.var_4 = StringVar()
        self.var_5 = StringVar()
        self.var_6 = StringVar()
        self.var_7 = StringVar()
        self.var_8 = StringVar()
        self.var_9 = StringVar()
        self.var_10 = StringVar()
        self.var_11 = StringVar()
        self.var_12 = StringVar()
        self.var_13 = StringVar()
        self.var_14 = StringVar()
        self.var_15 = StringVar()
        self.var_16 = StringVar()
        self.var_17 = StringVar()

        #################first layer Qdamp, Qbroad, Uiso_inter, Uiso_intra, Qmin, Qmax##########################
        #Qdamp
        ttk.Label(self.top_frame, text = "Qdamp", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.Qdamp).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_0,
                        style = "TCheckbutton").grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #Qbroad
        ttk.Label(self.top_frame, text = "Qbroad", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Arial", 10),
                     textvariable = self.Qbroad).grid(row = 3, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_1,
                        style = "TCheckbutton").grid(row = 3, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso inter
        ttk.Label(self.top_frame, text = "Uiso_intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.Uiso_inter).grid(row = 3, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_2,
                        style = "TCheckbutton").grid(row = 3, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso intra
        ttk.Label(self.top_frame, text = "Uiso_inter", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.Uiso_intra).grid(row = 3, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_3,
                        style = "TCheckbutton").grid(row = 3, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #Qmin
        ttk.Label(self.top_frame, text = "Qmin", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.Qmin).grid(row = 3, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_4,
                        style = "TCheckbutton").grid(row = 3, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #Qmax
        ttk.Label(self.top_frame, text = "Qmax", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 3, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.Qmax).grid(row = 3, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_5,
                        style = "TCheckbutton").grid(row = 3, column = 17, columnspan = 1, padx = 5, sticky = "sw")

        #################second layer lat_a, lat_b. lat_c, alpha, beta, gamma##########################        
        #lat_a
        ttk.Label(self.top_frame, text = "Lat_a", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.lat_a).grid(row = 4, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_6,
                        style = "TCheckbutton").grid(row = 4, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #lat_b
        ttk.Label(self.top_frame, text = "Lat_b", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Arial", 10),
                     textvariable = self.lat_b).grid(row = 4, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_7,
                        style = "TCheckbutton").grid(row = 4, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #lat_c
        ttk.Label(self.top_frame, text = "Lat_c", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.lat_c).grid(row = 4, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_8,
                        style = "TCheckbutton").grid(row = 4, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #alpha
        ttk.Label(self.top_frame, text = "Alpha", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.alpha).grid(row = 4, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_9,
                        style = "TCheckbutton").grid(row = 4, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #beta
        ttk.Label(self.top_frame, text = "Beta", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.beta).grid(row = 4, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_10,
                        style = "TCheckbutton").grid(row = 4, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #gamma
        ttk.Label(self.top_frame, text = "Gamma", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 4, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.gamma).grid(row = 4, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_11,
                        style = "TCheckbutton").grid(row = 4, column = 17, columnspan = 1, padx = 5, sticky = "sw")


        #################third layer scale, zoom_mol, zoom_intra, delta2_cryst, delta2_mole, delta2_intra##########################        
        #scale
        ttk.Label(self.top_frame, text = "Scale", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 0, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.scale).grid(row = 5, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_12,
                        style = "TCheckbutton").grid(row = 5, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        #zoom_mol
        ttk.Label(self.top_frame, text = "Zoom_Mol", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8,font = ("Arial", 10),
                     textvariable = self.zoom_mol).grid(row = 5, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_13,
                        style = "TCheckbutton").grid(row = 5, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        #zoom_intra
        ttk.Label(self.top_frame, text = "Zoom_Intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.zoom_intra).grid(row = 5, column = 7, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_14,
                        style = "TCheckbutton").grid(row = 5, column = 8, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_cryst
        ttk.Label(self.top_frame, text = "Delta2_Cryst", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 9, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.delta2_cryst).grid(row = 5, column = 10, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_15,
                        style = "TCheckbutton").grid(row = 5, column = 11, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_mol
        ttk.Label(self.top_frame, text = "Delta2_Mol", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 12, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.delta2_mol).grid(row = 5, column = 13, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_16,
                        style = "TCheckbutton").grid(row = 5, column = 14, columnspan = 1, padx = 5, sticky = "sw")
        #delta2_intra
        ttk.Label(self.top_frame, text = "Delta2_Intra", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 5, column = 15, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10),
                     textvariable = self.delta2_intra).grid(row = 5, column = 16, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_17,
                        style = "TCheckbutton").grid(row = 5, column = 17, columnspan = 1, padx = 5, sticky = "sw")
        
        ## layout for step 3, Select Optimizers and run
        ttk.Label(self.top_frame, text = "Step 3: Select Optimzers", justify = CENTER,
                  font = ("Arial", 16, "bold")).grid(row = 6, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        self.combobox = ttk.Combobox(self.top_frame, font = ("Arial", 12), values = ["Least Squares", "fmin",
                              "Basin Hopping"])
        self.combobox.current(0)
        self.combobox.grid(row = 7, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        
        ttk.Button(self.top_frame, text = "Run The FIT!", command = self.run_the_fit,
                             style = "TButton").grid(row = 7, column = 3, columnspan = 3, padx = 5)
        ttk.Button(self.top_frame, text = "Reset To Default", command = self.reset_to_default,
                             style = "TButton").grid(row = 7, column = 6, columnspan = 3, padx = 5, sticky = "sw")

    ##the bottom frame, the interactive matplotlib canvas for interactive plot/visualization
        self.bottom_frame = ttk.Frame(self.master, padding = (10, 10))
        self.bottom_frame.pack()
        self.fig = plt.figure(figsize=(13, 6), dpi=100) ##create a figure; modify the size here
        self.fig.add_subplot()
    
        plt.title("PDF Model Fit To Crystalline Organic PDF")
        plt.xlabel("r $(\AA)$", labelpad = 10, fontsize = 15)
        plt.ylabel("PDF, G $(\AA^{-2})$", labelpad = 10, fontsize = 15)

        plt.xticks([])
        plt.yticks([])

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.bottom_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.bottom_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    
    def load_molecule():
        pass

    def load_crystal():
        pass
    
    def run_the_fit():
        pass
    
    def reset_to_default():
        pass

def main():
    root = Tk()
    GUI = PDF_model_fit(root)
    root.mainloop()

if __name__ == "__main__": main()
