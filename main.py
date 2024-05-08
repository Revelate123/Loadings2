# This script will calculate design actions for use in simple beam problems

# It will determine the moment and shear of a simply supported beam
# It will determine the required Second moment of inertia to meet deflection criteria

import scipy.special as special

import os
import matplotlib.pyplot as plt
import matplotlib
import pickle as pickle
import PySimpleGUI as sg
import csv
import Output
import steel_functions as st
from flask import Flask, flash, redirect, render_template, request, session, jsonify
from flask_session import Session
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import layouts




def open_file(Name, Project):
    try:
        with open(Project + '\\' + Name + '.pkl', 'rb') as inp:
            foo = pickle.load(inp)
            return foo
    except (OSError, IOError) as e:
        print('No info')


def Method(Name, Length, sample_rate, E, I, Iy, Loadings, Cpe_V, Cpi_V, Cpe_H, Cpi_H, P, values, Cantilever):
    Name1 = Beam(Name, Length, sample_rate, E, I, Iy, Loadings, Cpe_V, Cpi_V, Cpe_H, Cpi_H, P, values, Cantilever)
    Name1.LoadCase(Loadings, Cantilever)

    with open(values['Project'] + '\\' + Name + '.pkl', 'wb') as outp:
        pickle.dump(Name1, outp, pickle.HIGHEST_PROTOCOL)


def variable_write(values, Name):
    # if values['Member Type'] == 'Beam':
    with open(Name + '.txt', 'w', newline='') as csv_file:
        data = [[str(i), values[i]] for i in values]
        print(data)
        writer = csv.writer(csv_file)
        writer.writerows(data)


def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
    if canvas.children:
        for child in canvas.winfo_children():
            child.destroy()
    if canvas_toolbar.children:
        for child in canvas_toolbar.winfo_children():
            child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
    figure_canvas_agg.draw()
    toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
    toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)


class Toolbar(NavigationToolbar2Tk):
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)


def variables(Name):
    Dictionary = {}
    with open(str(Name) + '.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            Dictionary[row[0]] = row[1]
    return Dictionary






def delete_figure_agg(figure_agg):
    figure_agg.get_tk_widget().forget()
    plt.close('all')


def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def draw_graph(Beam,i):
    # fig = plt.figure(figsize=(12,30),dpi=80)
    fig = matplotlib.figure.Figure(figsize=(12, 4),dpi = 80)
    try:
        # Length = np.arange(0,Beam.Length + Beam.values['CLength']+ 1/Beam.sample_rate,1/Beam.sample_rate)
        Length = range(round(Beam.Length * Beam.sample_rate + Beam.values['CLength'] * Beam.sample_rate + 1))
        Length = [i / Beam.sample_rate for i in Length]
        Zeros = [0] * (int(Beam.Length + Beam.values['CLength']) * Beam.sample_rate + 1)
    except:

        # Length =  np.arange(0,Beam.Length + 1/Beam.sample_rate,1/Beam.sample_rate)
        Length = range(round(Beam.Length * Beam.sample_rate + 1))
        Length = [i / Beam.sample_rate for i in Length]



    j = 0
    R1 = float(Beam.values['Load_CaseG1'+str(i)])* getattr(Beam,'G1').R1 + float(Beam.values['Load_CaseG2'+str(i)])* getattr(Beam,'G2').R1 + float(Beam.values['Load_CaseQ'+str(i)])* getattr(Beam,'Q').R1 - float(Beam.values['Load_CaseWup'+str(i)])* getattr(Beam,'Wup').R1 + float(Beam.values['Load_CaseWdown'+str(i)])* getattr(Beam,'Wdown').R1 + float(Beam.values['Load_CaseWOoP'+str(i)])* getattr(Beam,'WOoP').R1
    R2 = float(Beam.values['Load_CaseG1'+str(i)])* getattr(Beam,'G1').R2 + float(Beam.values['Load_CaseG2'+str(i)])* getattr(Beam,'G2').R2 + float(Beam.values['Load_CaseQ'+str(i)])* getattr(Beam,'Q').R2 - float(Beam.values['Load_CaseWup'+str(i)])* getattr(Beam,'Wup').R2 + float(Beam.values['Load_CaseWdown'+str(i)])* getattr(Beam,'Wdown').R2 + float(Beam.values['Load_CaseWOoP'+str(i)])* getattr(Beam,'WOoP').R2
    if Beam.values['Load_CaseM'+str(i)] ==True:
        j += 1
        Combined_Moment = [float(Beam.values['Load_CaseG1'+str(i)])* x1 + float(Beam.values['Load_CaseG2'+str(i)])* x2 + float(Beam.values['Load_CaseQ'+str(i)])* x3 - float(Beam.values['Load_CaseWup'+str(i)])* x4 + float(Beam.values['Load_CaseWdown'+str(i)])* x5 + float(Beam.values['Load_CaseWOoP'+str(i)])* x6 for x1, x2, x3, x4, x5, x6 in zip(getattr(Beam, 'G1').Moment,getattr(Beam, 'G2').Moment,getattr(Beam, 'Q').Moment,getattr(Beam, 'Wup').Moment,getattr(Beam, 'Wdown').Moment,getattr(Beam, 'WOoP').Moment)]


    if Beam.values['Load_CaseS' + str(i)] == True:
        j +=1
        Combined_Shear = [float(Beam.values['Load_CaseG1'+str(i)]) * x1 + float(Beam.values['Load_CaseG2'+str(i)]) * x2 + float(
            Beam.values['Load_CaseQ'+str(i)]) * x3 - float(Beam.values['Load_CaseWup'+str(i)]) * x4 + float(
            Beam.values['Load_CaseWdown'+str(i)]) * x5 + float(Beam.values['Load_CaseWOoP'+str(i)]) * x6 for x1, x2, x3, x4, x5, x6 in
                    zip(getattr(Beam, 'G1').Shear, getattr(Beam, 'G2').Shear,getattr(Beam, 'Q').Shear,
                        getattr(Beam, 'Wup').Shear, getattr(Beam, 'Wdown').Shear, getattr(Beam, 'WOoP').Shear)]
    if Beam.values['Load_CaseD' + str(i)] == True:
        j += 1
        Combined_Deflection = [float(Beam.values['Load_CaseG1'+str(i)]) * x1 + float(Beam.values['Load_CaseG2'+str(i)]) * x2 + float(
            Beam.values['Load_CaseQ'+str(i)]) * x3 - float(Beam.values['Load_CaseWup'+str(i)]) * x4 + float(
            Beam.values['Load_CaseWdown'+str(i)]) * x5 + float(Beam.values['Load_CaseWOoP'+str(i)]) * x6 for x1, x2, x3, x4, x5, x6 in
                    zip(getattr(Beam, 'G1').Deflection, getattr(Beam, 'G2').Deflection,getattr(Beam, 'Q').Deflection,
                        getattr(Beam, 'Wup').Deflection, getattr(Beam, 'Wdown').Deflection, getattr(Beam, 'WOoP').Deflection)]

        # Factor deflections for timber creep factors
        if Beam.values['SectionType'] == 'MGP 10' or Beam.values['SectionType'] == 'MGP 12':
            if Beam.values['j2'] == '<15%':
                Combined_Deflection = [2 * i for i in Combined_Deflection]

            else:
                Combined_Deflection = [3 * i for i in Combined_Deflection]
        # Find the Maximum deflection from the load case selected and add to the Beam object

        getattr(Beam, 'MaxDeflections')[i] = max(Combined_Deflection, key=abs)
    j1 = 1
    if Beam.values['Load_CaseM'+str(i)] ==True:
        ax1 = fig.add_subplot(1, j, j1)
        ax1.set_ylabel('Moment [KNm]')
        ax1.set_xlabel('Length [m]')
        ax1.set_title(Beam.values['Load_Case'+str(i)])
        ax1.plot(Length, Combined_Moment)
        ax1.text(0, -abs(1 / 10 * max(Combined_Moment, key=abs)), str(round(R1, 1)) + 'KN',
                 verticalalignment='top', horizontalalignment='center')
        ax1.text(Beam.Length, -abs(1 / 10 * max(Combined_Moment, key=abs)), str(round(R2, 1)) + 'KN',
                 verticalalignment='top', horizontalalignment='center')
        ax1.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
        try:
            ax1.plot([0, Beam.Length], [-Beam.section_properties[0][0]/1000, -Beam.section_properties[0][0]/1000], '--r')
        except:
            pass
        try:
            ax1.plot([0, Beam.Length], [-Beam.section_properties[0]['PhiMbx'], -Beam.section_properties[0]['PhiMbx']], '--r')
        except:
            pass
        try:
            ax1.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
        except:
            pass
        j1 += 1

    if Beam.values['Load_CaseS' + str(i)] == True:
        ax1 = fig.add_subplot(1, j, j1)
        ax1.set_ylabel('Shear [KN]')
        ax1.set_xlabel('Length [m]')
        ax1.set_title(Beam.values['Load_Case'+str(i)])
        ax1.plot(Length, Combined_Shear)
        ax1.text(0, -abs(1 / 10 * max(Combined_Shear, key=abs)), str(round(R1, 1)) + 'KN',
                 verticalalignment='top', horizontalalignment='center')
        ax1.text(Beam.Length, -abs(1 / 10 * max(Combined_Shear, key=abs)), str(round(R2, 1)) + 'KN',
                 verticalalignment='top', horizontalalignment='center')
        ax1.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
        try:
            ax1.plot([0, Beam.Length], [-Beam.section_properties[1][0]/1000, -Beam.section_properties[1][0]/1000], '--r')
        except:
            pass
        try:
            ax1.plot([0, Beam.Length], [-Beam.section_properties[0]['PhiVu'], -Beam.section_properties[0]['PhiVu']], '--r')
        except:
            pass
        try:
            ax1.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
        except:
            pass
        j1 +=1

    if Beam.values['Load_CaseD' + str(i)] == True:
        ax1 = fig.add_subplot(1, j, j1)
        ax1.set_ylabel('Deflection [mm]')
        ax1.set_xlabel('Length [m]')
        ax1.set_title(Beam.values['Load_Case'+str(i)])
        ax1.plot(Length, Combined_Deflection)
        ax1.text(0, -abs(1 / 10 * max(Combined_Deflection, key=abs)), str(round(R1, 1)) + 'KN',
                 verticalalignment='top', horizontalalignment='center')
        ax1.text(Beam.Length, -abs(1 / 10 * max(Combined_Deflection, key=abs)), str(round(R2, 1)) + 'KN',
                 verticalalignment='top', horizontalalignment='center')
        ax1.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
        ax1.plot([0, Beam.Length], [-max(Length)/Beam.values['Load_CaseD_L'+str(i)]*1000, -max(Length)/Beam.values['Load_CaseD_L'+str(i)]*1000], '--r')
        try:
            ax1.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
            ax1.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [-Beam.values['CLength'] / 0.180, -Beam.values['CLength'] / 0.180], '--r')
        except:
            pass
        j1 += 1

    fig.tight_layout()
    getattr(Beam,'Combined_Deflections')[i] = fig
    return fig

def PassFail(Beam):
    String = ''
    try:
        # Length = np.arange(0,Beam.Length + Beam.values['CLength']+ 1/Beam.sample_rate,1/Beam.sample_rate)
        Length = range(round(Beam.Length * Beam.sample_rate + Beam.values['CLength'] * Beam.sample_rate + 1))
        Length = [i / Beam.sample_rate for i in Length]
        Zeros = [0] * (int(Beam.Length + Beam.values['CLength']) * Beam.sample_rate + 1)
    except:

        # Length =  np.arange(0,Beam.Length + 1/Beam.sample_rate,1/Beam.sample_rate)
        Length = range(round(Beam.Length * Beam.sample_rate + 1))
        Length = [i / Beam.sample_rate for i in Length]
    for i in range(int(Beam.values['Load_Cases'])):
        try:
            if abs(max(Length) / Beam.values['Load_CaseD_L' + str(i)] * 1000) >= abs(Beam.MaxDeflections[i]):
                String += Beam.values['Load_Case'+str(i)] + ' = Pass | Limit = ' + str(round(abs(max(Length) / Beam.values['Load_CaseD_L' + str(i)] * 1000),1)) +'mm > '+ str(round(abs(Beam.MaxDeflections[i]),1)) +'mm\n'
            else:
                String += Beam.values['Load_Case'+str(i)] + ' = Fail | Limit = '+ str(round(abs(max(Length) / Beam.values['Load_CaseD_L' + str(i)] * 1000),1)) +'mm < '+ str(round(abs(Beam.MaxDeflections[i]),1)) + 'mm\n'
        except:
            continue
    return String
    # Create an event loop


variable = variables('Loading')
#Layouts(variable)




app = Flask(__name__)

# Configure session to use filesystem (instead of signed cookies)
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)



@app.route('/', methods=['GET', 'POST'])
def main():
    if request.method == "POST":
        pass
    else:
        return render_template("beam_calculator.html")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8080)

