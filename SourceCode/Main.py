""""
This is the code for OrganicPDF GUI program. It is used to extract information from a measured Organic PDF.
The code is written by Chenyang Shi at AbbVie.
Send him an email at chenyang.shi@abbvie.com for questions and suggestions. 

"""

from Tkinter import *
import Tkinter as tk
import ttk
from math import exp
from PDF_model_fit import PDF_model_fit

import os  # for loading files or exporting files


class Overall_Look:
    def __init__(self, master):
        self.master = master
        self.master.title("Organic PDF GUI")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(600, 300) # width + height
        self.master.resizable(False, False)


##define style in ttk##
        self.style = ttk.Style()
        self.style.configure('TFrame', background = 'grey91')
        self.style.configure('TButton', background = 'grey91')
        self.style.configure('TLabel', background = 'grey91', font = ('Arial', 18))
        self.style.configure('Header.TLabel', font = ('Arial', 24, 'bold'))

        
##load image##
        self.logo = PhotoImage(file = "Logo.gif").subsample(4,4)

##this is the top frame, we will write a few sentences about this program##
        self.frame_header = ttk.Frame(master, relief = RIDGE, padding = (30, 15))
        self.frame_header.pack()

        ttk.Label(self.frame_header, image = self.logo).grid(row = 0, column = 0, rowspan = 2, columnspan = 2, sticky = "w")
        ttk.Label(self.frame_header, text = 'Welcome to use Organic PDF GUI program!',
                  style = 'Header.TLabel').grid(row = 0, column = 2, columnspan = 4)
        ttk.Label(self.frame_header, wraplength = 600,
                  text = ("This program is designed to extract structural information from a measured orgainc PDF from synchrotron X-ray or neutron sources or laboroatry X-ray facilities. Currently it supports (1) Study of intermolecular hydrogen bond by subtracting signal of a single molecule. (2) PDF model fit of a crystalline orgainc PDF using the method proposed by Prill et al. (J. Appl. Cryst., 2015, 48, 171-178. (3) Solve crystal structure from powder PDF data using the method proposed by Prill et al. (Acta Crystallogr. A.,2016, 72, 62-72))"
                          )).grid(row = 1, column = 2, columnspan = 4)
        

        self.bottom_frame = ttk.Frame(master, relief = FLAT, padding = (30, 15)) 
        self.bottom_frame.pack()
        
        Button(self.bottom_frame, text = 'Intermolecular PDF', relief = RAISED,
                    command = self.inter_pdf, font = ('Arial', 20, "bold")).grid(row = 0, column = 0, columnspan = 3, padx = 5, pady = 5)
        Button(self.bottom_frame, text = 'Crystalline PDF fit', relief = RAISED,
                    command = self.pdf_model, font = ('Arial', 20, "bold")).grid(row = 0, column = 5, columnspan = 3, padx = 5, pady = 5)
        Button(self.bottom_frame, text = 'Structure Solve', relief = RAISED,
                    command = self.stru_sol, font = ('Arial', 20, "bold")).grid(row = 0, column = 10, columnspan = 3, padx = 5, pady = 5)
    
    def inter_pdf(self):
        self.pdf_model = tk.Toplevel(self.master)
        self.GUI = Inter_PDF_study(self.inter_pdf)


    def pdf_model(self):
        self.pdf_model = tk.Toplevel(self.master)
        self.GUI = PDF_model_fit(self.pdf_model)
    

    def stru_sol(self):
        self.pdf_model = tk.Toplevel(self.master)
        self.GUI = Structure_Solve(self.stru_sol)

def main():
    root = Tk()
    GUI = Overall_Look(root)
    root.mainloop()

if __name__ == "__main__": main()
