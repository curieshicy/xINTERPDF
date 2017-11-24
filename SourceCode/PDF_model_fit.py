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


class PDF_model_fit():
    def __init__(self, master):
        self.master = master
        self.master.title("Model fit to organic PDF")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(900, 650) # width + height
        self.master.resizable(False, False)
        
        
        ##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91')
        self.style.configure('TCheckbutton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91', font = ('Arial', 18))
        self.style.configure('Header.TLabel', font = ('Arial', 24, 'bold'))

        ttk.Label(self.master, text = "To perform a fit to crystalline organic PDF, one needs to prepare a strcuture file for single molecule.\n"
            "In addition, a crystal structure for compound is needed as input.", font = ("Arial", 16, "bold"), justify = CENTER).pack(side = TOP)
        self.top_frame = ttk.Frame(self.master, padding = (30, 15))
        self.top_frame.pack()
        
        ##here are the layout for step 1, load structure files

        ttk.Label(self.top_frame, text = "Step 1, Load Structure Files", justify = CENTER, font = ("Arial", 15, "bold")).grid(row = 0, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Label(self.top_frame, text = "Load molecule structure", justify = LEFT, font = ("Arial", 15, "bold")).grid(row = 1, column = 0, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load molecule structure", command = self.load_molecule, style = "TButton").grid(row = 1, column = 3, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Label(self.top_frame, text = "Load crystal structure", justify = RIGHT, font = ("Arial", 15, "bold")).grid(row = 1, column = 6, columnspan = 2, padx = 5, sticky = "sw")

        ttk.Button(self.top_frame, text = "Load crystal structure",command = self.load_crystal, style = "TButton").grid(row = 1, column = 12, columnspan = 2, padx = 5, sticky = "sw")
    
        ## layout for step 2, set parameters
        ttk.Label(self.top_frame, text = "Step 2, Set Parameters", justify = CENTER, font = ("Arial", 15, "bold")).grid(row = 2, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Label(self.top_frame, text = "Instrument/Thermals", justify = LEFT, font = ("Arial", 15, "bold")).grid(row = 3, column = 0, columnspan = 2, padx = 5, sticky = "sw")
        

        self.Qdamp = StringVar()
        self.Qbroad = StringVar()
        self.var_0  = StringVar()
        self.var_1 = StringVar()
        self.var_2 = StringVar()
        
        #Qdamp
        ttk.Label(self.top_frame, text = "Qdamp", justify = LEFT, font = ("Arial", 15, "bold")).grid(row = 3, column = 1, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 10, font = ("Arial", 10), textvariable = self.Qdamp).grid(row = 3, column = 2, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_0, style = "TCheckbutton").grid(row = 3, column = 3, columnspan = 1, padx = 5, sticky = "sw")
        #Qbroad
        ttk.Label(self.top_frame, text = "Qbroad", justify = LEFT, font = ("Arial", 15, "bold")).grid(row = 3, column = 4, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10), textvariable = self.Qdamp).grid(row = 3, column = 5, columnspan = 1, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_0, style = "TCheckbutton").grid(row = 3, column = 6, columnspan = 1, padx = 5, sticky = "sw")
        #Uiso inter
        ttk.Label(self.top_frame, text = "Uiso_intra", justify = LEFT, font = ("Arial", 15, "bold")).grid(row = 3, column = 7, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10), textvariable = self.Qdamp).grid(row = 3, column = 8, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_1, style = "TCheckbutton").grid(row = 3, column = 9, columnspan = 2, padx = 5, sticky = "sw")
        #Uiso intra
        ttk.Label(self.top_frame, text = "Uiso_inter", justify = LEFT, font = ("Arial", 15, "bold")).grid(row = 3, column = 10, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Entry(self.top_frame, width = 8, font = ("Arial", 10), textvariable = self.Qdamp).grid(row = 3, column = 11, columnspan = 2, padx = 5, sticky = "sw")
        ttk.Checkbutton(self.top_frame, variable = self.var_2, style = "TCheckbutton").grid(row = 3, column = 12, columnspan = 2, padx = 5, sticky = "sw")
    




    def load_molecule():
        pass

    def load_crystal():
        pass


def main():
    root = Tk()
    GUI = PDF_model_fit(root)
    root.mainloop()

if __name__ == "__main__": main()
