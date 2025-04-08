#This script will print a pdf of the beam being investigated
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import pickle as pickle
import FreeSimpleGUI as sg
import csv
import steel_functions as st
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, MultiRow, MultiColumn
from pylatex.utils import italic, NoEscape, bold
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os
import platform
import subprocess

def Main(Beam):
    #Loop saves each figure which is stored within the object beam
    for i in Beam.Combined_Deflections:
        Beam.Combined_Deflections[i].savefig('Fig'+str(i)+'.png')

    #Beam.Combined_Deflections[0].savefig('Fig1.png')
    #Set Geometry options for the document
    geometry_options = {'tmargin': '1cm', 'lmargin': '1cm'}
    #Create the document
    doc = Document(Beam.values['Name'], geometry_options=geometry_options)

    #Begin the first section of the document
    with doc.create(Section('Design of '+Beam.values['Name'], numbering=False)):
        #Create the first sub heading
        doc.append(Subsection('Beam Properties',numbering=False))
        doc.append(NoEscape('Section Type:\t' + Beam.values['SectionType']))
        doc.append('\n\n')
        doc.append(NoEscape('Section Size:\t'+Beam.values['SectionSize']))
        doc.append('\n\n')
        doc.append(NoEscape('Length:\t' + str(Beam.values['Length'])))
        doc.append('\n\n')
        if Beam.values['Cantilever'] ==True:
            doc.append(NoEscape('Cantilever Length:\t' + str(Beam.values['CLength'])))
            doc.append('\n\n')
        doc.append(NoEscape('Modulus of Elasticity, E =\t' + str(Beam.E)))
        doc.append('\n\n')
        doc.append(NoEscape('Second Moment of Area, Iy =\t' + str(Beam.Iy)))
        doc.append('\n\n')
        doc.append(NoEscape('Second Moment of Area, Ix =\t' + str(Beam.Ix)))
        doc.append('\n\n')
        try:
            doc.append(NoEscape('Moment Capacity, $\phi M_{n}$ = \t' + str(round(Beam.section_properties[0][0]/1000,2)) + ' KNm'))
            doc.append('\n\n')
            #doc.append('\nMoment Capacity, \t' + str(Beam.section_properties[0][0]/1000))
        except:
            pass
        try:
            doc.append('Moment Capacity\t' + str(Beam.section_properties[0]['PhiMbx']))
            doc.append('\n\n')
        except:
            pass
        doc.append(Subsection('Deflection Checks', numbering=False))
        for i in Beam.Combined_Deflections:
            with doc.create(Figure(position='htbp')) as plot:
                plot.add_image('Fig'+str(i)+'.png', width=NoEscape(r'0.8\textwidth'))


    doc.generate_tex()
    subprocess.run(['pdflatex', '-interaction=nonstopmode', Beam.values['Name'] + '.tex'])
    os.remove("Fig1.png")
    os.remove(Beam.values['Name'] + '.tex')
    os.remove(Beam.values['Name'] + '.aux')
    os.remove(Beam.values['Name'] + '.log')