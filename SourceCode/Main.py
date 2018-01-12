""""
This is the code for OrganicPDF GUI program. It is used to extract information from a measured Organic PDF.
The code is written by Chenyang Shi at AbbVie.
Send him an email at chenyang.shi@abbvie.com for questions and suggestions. 

"""

from Tkinter import *
import Tkinter as tk
import ttk
from math import exp
from PDF_model_fit_Uiso import PDF_model_fit_Uiso
from InterPDF import InterPDF

import os  # for loading files or exporting files


class Overall_Look:
    def __init__(self, master):
        self.master = master
        self.master.title("Organic PDFGUI")
        self.master.configure(background = "grey91") #the color will be changed later
        self.master.minsize(600, 250) # width + height
        self.master.resizable(False, False)
        
## add menu bar
        self.menubar = Menu(self.master)
        ## Dropdown menu 1
        self.file_menu = Menu(self.menubar, tearoff = 0)
        self.file_menu.add_command(label = "New")
        self.file_menu.add_command(label = "Open")
        self.file_menu.add_command(label = "Save")
        self.file_menu.add_command(label = "Save As")
        self.file_menu.add_command(label = "Exit")
        self.menubar.add_cascade(label = "File", menu = self.file_menu)
        
        ##Dropdown menu 2
        self.interpdf_menu = Menu(self.menubar, tearoff = 0)
        self.interpdf_menu.add_command(label = "Inter PDF", command = self.inter_pdf)
        self.menubar.add_cascade(label = "Inter PDF", menu = self.interpdf_menu)

        ##Dropdown menu 3
        self.pdffit_menu = Menu(self.menubar, tearoff = 0)
        self.pdffit_menu.add_command(label = "PDF Fit", command = self.pdf_model_uiso)
        self.menubar.add_cascade(label = "PDF Fit", menu = self.pdffit_menu)
        
#        ##Dropdown menu 4
#        self.pdfsolve_menu = Menu(self.menubar, tearoff = 0)
#        self.pdfsolve_menu.add_command(label = "PDF Solve")
#        self.menubar.add_cascade(label = "PDF_Solve", menu = self.pdfsolve_menu)

        ##Dropdown menu 5
        self.help_menu = Menu(self.menubar, tearoff = 0)
        self.help_menu.add_command(label = "Help")
        self.help_menu.add_command(label = "About")
        self.menubar.add_cascade(label = "Help", menu = self.help_menu)
        
        self.master.config(menu = self.menubar)


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
        ttk.Label(self.frame_header, text = 'Welcome to use Organic PDFGUI program!',
                  style = 'Header.TLabel').grid(row = 0, column = 2, columnspan = 4)
        ttk.Label(self.frame_header, wraplength = 600,
                  text = ("This program is designed to extract structural information from measured orgainc PDFs from synchrotron X-ray or neutron sources or laboratory X-ray facilities. Currently it supports (1) Study of intermolecular hydrogen bond by subtracting signal of a single molecule. (2) PDF model fit of a crystalline orgainc PDF using the method proposed by Prill et al. (J. Appl. Cryst., 2015, 48, 171-178. (3) Solve crystal structure from powder PDF data using the method proposed by Prill et al. (Acta Crystallogr. A., 2016, 72, 62-72))"
                          )).grid(row = 1, column = 2, columnspan = 4)
        

        self.bottom_frame = ttk.Frame(master, relief = FLAT, padding = (30, 15)) 
        self.bottom_frame.pack()

#    def inter_pdf(self):
#        self.inter_pdf = tk.Toplevel(self.master)
#        self.GUI = Inter_PDF_study(self.inter_pdf)

    def pdf_model_uiso(self):
        self.pdf_model_uiso = tk.Toplevel(self.master)
        self.GUI = PDF_model_fit_Uiso(self.pdf_model_uiso)
    
#    def stru_sol(self):
#        self.stru_sol = tk.Toplevel(self.master)
#        self.GUI = Structure_Solve(self.stru_sol)

    def inter_pdf(self):
        self.inter_pdf = tk.Toplevel(self.master)
        self.GUI = InterPDF(self.inter_pdf)


def main():
    root = Tk()
    GUI = Overall_Look(root)
    root.mainloop()

if __name__ == "__main__": main()
