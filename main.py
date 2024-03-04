# This script will calculate design actions for use in simple beam problems

# It will determine the moment and shear of a simply supported beam
# It will determine the required Second moment of inertia to meet deflection criteria
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import pickle as pickle
import PySimpleGUI as sg
import csv
import Output
import steel_functions as st
import Timber as Timber1
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, MultiRow, MultiColumn
from pylatex.utils import italic, NoEscape, bold
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


def Loads(G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit, Gpoint, Qpoint, Wpointsls, a):
    wdown12G15Q = 1.2 * G * Trib_width + 1.5 * Q * Trib_width
    wdown135G = 1.35 * G * Trib_width
    wdown12GWuls = 1.2 * G * Trib_width + Wsls * 1.48 * Trib_width
    wup = min(0.9 * G * Trib_width - Wsls * 1.48 * Trib_width, 0)
    wdown = max(wdown135G, wdown12G15Q, wdown12GWuls)
    Mulsdown = wdown * Length ** 2 / 8
    Vulsdown = wdown * Length / 2
    Mulsup = wup * Length ** 2 / 8
    Vulsup = wup * Length / 2

    w12G15Q = 1.2 * Gpoint + 1.5 * Qpoint
    w135G = 1.35 * Gpoint
    w12GWuls = 1.2 * Gpoint + 1.48 * Wpointsls
    wdownpoint = max(w12G15Q, w135G, w12GWuls)
    wuppoint = min(0.9 * Gpoint - 1.48 * Wpointsls, 0)

    Mulsdownpoint = wdownpoint * a * (Length - a) / Length
    Mulsuppoint = wuppoint * a * (Length - a) / Length
    x = min(a, (Length - a))
    y = max(wdownpoint, -1 * wuppoint)
    Vpoint = y * x / Length

    GI = 5 * G * Trib_width * 10 ** 3 * Length ** 4 / (
                384 * 200 * 10 ** 9 * Glimit * 10 ** -3) * 10 ** 12 + Gpoint * 10 ** 3 * Length * 10 ** 3 / (
                     48 * 200 ** 9 * Glimit * 10 ** -3)
    QI = 5 * 0.7 * Q * Trib_width * 10 ** 3 * Length ** 4 / (
                384 * 200 * 10 ** 9 * Qlimit * 10 ** -3) * 10 ** 12 + Qpoint * 10 ** 3 * Length * 10 ** 3 / (
                     48 * 200 ** 9 * Qlimit * 10 ** -3)
    WI = 5 * Wsls * Trib_width * 10 ** 3 * Length ** 4 / (
                384 * 200 * 10 ** 9 * Wlimit * 10 ** -3) * 10 ** 12 + Wpointsls * 10 ** 3 * Length * 10 ** 3 / (
                     48 * 200 ** 9 * Wlimit * 10 ** -3)

    return Mulsdown, max(Vulsdown,
                         Vulsup), Mulsup, GI, QI, WI, wdown, wup, wdown12G15Q, wdown135G, wdown12GWuls, Mulsdownpoint, Mulsuppoint, Vpoint, wdownpoint, wuppoint, Mulsdownpoint, Mulsuppoint, Vpoint


def Write_Loads(Name, G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit, Gpoint, Qpoint, Wpointsls, a):
    Load = Loads(G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit, Gpoint, Qpoint, Wpointsls, a)
    with open(Name + '.txt', 'w') as f:
        f.write(Name)
        f.write('\n')
        f.write('\n')
        f.write('Length:  ' + str(Length))
        f.write('\nG:  ' + str(G))
        f.write('\nQ:  ' + str(Q))
        f.write('\nWsls:  ' + str(Wsls))
        f.write('\nTrib Width:  ' + str(Trib_width))
        f.write('\nG Point:  ' + str(Gpoint))
        f.write('\nQ Point:  ' + str(Qpoint))
        f.write('\nWsls Point:  ' + str(Wpointsls))
        f.write('\n1.2G + 1.5Q:  ' + str(Load[8]) + ' kN/m')
        f.write('\n1.35G:  ' + str(Load[9]) + ' kN/m')
        f.write('\n1.2G + 1.48Wsls:  ' + str(Load[10]) + ' kN/m')
        f.write('\nMaximum UDL down:' + str(Load[6]) + ' kN/m')
        f.write('\nMaximum UDL up:' + str(Load[7]) + ' kN/m')
        f.write('\nMaximum point load down:  ' + str(Load[14]))
        f.write('\nMaximum point load up:  ' + str(Load[15]))
        f.write('\nMaximum positive Moment due to UDL:  ' + str(Load[0]) + ' kNm')
        f.write('\nMaximum negative Moment due to UDL:  ' + str(Load[2]) + ' kNm')
        f.write('\nMaximum shear due to UDL:  ' + str(Load[1]) + ' kN')
        f.write('\nMaximum positive moment due to point load:' + str(Load[16]))
        f.write('\nMaximum negative moment due to point load:' + str(Load[17]))
        f.write('\nMaximum Shear due to point load:' + str(Load[18]))
        f.write('\n\nRequired Ix for G:  ' + str(Load[3]) + ' mm4')
        f.write('\n\nRequired Ix for phiQ:  ' + str(Load[4]) + ' mm4')
        f.write('\n\nRequired Ix for Wsls:  ' + str(Load[5]) + ' mm4')


Name = 'xx'
Length = 7.4
Glimit = min(Length * 1000 / 360, 100)
Qlimit = min(Length * 1000 / 250, 100)
Wlimit = min(Length * 1000 / 150, 100)
G = 0.5
Q = 2
Wsls = 0
Trib_width = 2.5
Gpoint = 0
Qpoint = 0
Wpointsls = 0
a = Length / 2

Write_Loads(Name, G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit, Gpoint, Qpoint, Wpointsls, a)


class Load_Case:
    def __init__(self, Name, Length, sample_rate, E, I, values):
        self.Name = Name
        self.Length = Length  # In metres
        self.Moment = [0] * round(Length * sample_rate + 1)
        self.Shear = [0] * round(Length * sample_rate + 1)
        self.Deflection = [0] * round(Length * sample_rate + 1)
        try:
            if values['CLength'] > 0:
                self.Deflection = [0] * round(values['CLength'] * sample_rate + Length * sample_rate + 1)
                self.Moment = [0] * round(values['CLength'] * sample_rate + Length * sample_rate + 1)
                self.Shear = [0] * round(values['CLength'] * sample_rate + Length * sample_rate + 1)
        except:
            pass
        self.R1 = 0
        self.R2 = 0
        self.sample_rate = sample_rate
        self.Marea = [0] * round(Length * sample_rate + 1)
        self.E = E
        self.I = I
        self.values = values

    def __str__(self):
        return f"{self.R1} {self.R2} {self.sample_rate}"

    def UDL(self, a, b, w):
        chec = round(self.Length * self.sample_rate + 1)
        Moment = [0] * round(self.Length * self.sample_rate + 1)
        Shear = [0] * round(self.Length * self.sample_rate + 1)
        Deflection = [0] * round(self.Length * self.sample_rate + 1)
        Marea = [0] * round(self.Length * self.sample_rate + 1)
        c = self.Length - b
        b = b - a
        R1 = w * b / 2 / self.Length * (2 * c + b)
        R2 = w * b / 2 / self.Length * (2 * a + b)
        for i in range(round(self.Length * self.sample_rate + 1)):
            if i <= a * self.sample_rate:
                Shear[i] = R1
                Moment[i] = R1 * i / self.sample_rate
            elif i > a * self.sample_rate and i < (a + b) * self.sample_rate:
                Shear[i] = R1 - w * (i / self.sample_rate - a)
                Moment[i] = R1 * i / self.sample_rate - w / 2 * (i / self.sample_rate - a) ** 2
            elif i >= (a + b) * self.sample_rate:
                Shear[i] = -R2
                Moment[i] = R2 * (self.Length - i / self.sample_rate)
                # print(i)
        for i in range(round(self.Length * self.sample_rate)):
            Marea[i] = np.trapz(np.array(Moment[:i]), dx=1 / self.sample_rate)
        centroid = min(range(len(Marea)), key=lambda x: abs(Marea[x] - max(Marea) / 2))
        x1 = lambda x: R1 * x ** 2 / 2
        y1 = integrate.quad(x1, 0, a)
        x1a = lambda x: R1 * x
        y1a = integrate.quad(x1a, 0, a)
        x2 = lambda x: (R1 * x - w / 2 * (x - a) ** 2) * x
        x2a = lambda x: R1 * x - w / 2 * (x - a) ** 2
        y2 = integrate.quad(x2, a, a + b)
        y2a = integrate.quad(x2a, a, a + b)
        x3 = lambda x: R2 * (self.Length - x) * x
        x3a = lambda x: R2 * (self.Length - x)
        y3 = integrate.quad(x3, a + b, self.Length)
        y3a = integrate.quad(x3a, a + b, self.Length)
        try:
            End_deflection = (self.Length - (y1[0] + y2[0] + y3[0]) / (y1a[0] + y2a[0] + y3a[0])) * (
                        y1a[0] + y2a[0] + y3a[0])
            End_deflection1 = Marea[int(self.Length * self.sample_rate)] * (self.Length - centroid / self.sample_rate)
            print(End_deflection1, End_deflection)
            for i in range(round(self.Length * self.sample_rate)):
                if i == 0:
                    Deflection[i] = 0
                elif i < a * self.sample_rate:
                    y1 = integrate.quad(x1, 0, i / self.sample_rate)
                    y1a = integrate.quad(x1a, 0, i / self.sample_rate)

                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                                i / self.sample_rate - (y1[0]) / (y1a[0])) * (y1a[0])) / (self.E * self.I) * 10 ** 6

                elif i < a * self.sample_rate + b * self.sample_rate and i >= a * self.sample_rate:
                    y1 = integrate.quad(x1, 0, a)
                    y1a = integrate.quad(x1a, 0, a)
                    y2 = integrate.quad(x2, a, i / self.sample_rate)
                    y2a = integrate.quad(x2a, a, i / self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                                i / self.sample_rate - (y1[0] + y2[0]) / (y1a[0] + y2a[0])) * (y1a[0] + y2a[0])) / (
                                                self.E * self.I) * 10 ** 6
                elif i >= a * self.sample_rate + b * self.sample_rate:
                    y1 = integrate.quad(x1, 0, a)
                    y1a = integrate.quad(x1a, 0, a)
                    y2 = integrate.quad(x2, a, a + b)
                    y2a = integrate.quad(x2a, a, a + b)
                    y3 = integrate.quad(x3, a + b, i / self.sample_rate)
                    y3a = integrate.quad(x3a, a + b, i / self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                                i / self.sample_rate - (y1[0] + y2[0] + y3[0]) / (y1a[0] + y2a[0] + y3a[0])) * (
                                                 y1a[0] + y2a[0] + y3a[0])) / (self.E * self.I) * 10 ** 6
        except:
            pass
        for i in range(round(self.Length * self.sample_rate + 1)):
            self.Moment[i] += -Moment[i]
            self.Shear[i] += Shear[i]
            self.Deflection[i] += -Deflection[i]
            self.Marea[i] += Marea[i]
        try:

            for i in range(round(self.Length * self.sample_rate + 1),
                           round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):
                self.Deflection[i] += -(i / (self.Length * self.sample_rate) * End_deflection - (
                            i / self.sample_rate - (y1[0] + y2[0] + y3[0]) / (y1a[0] + y2a[0] + y3a[0])) * (
                                                    y1a[0] + y2a[0] + y3a[0])) / (self.E * self.I) * 10 ** 6
        except:
            pass
        self.R1 += R1
        self.R2 += R2

    def Cantilever_UDL(self, a, b, w):
        b = b - a
        Moment = [0] * round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        Shear = [0] * round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        Deflection = [0] * round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        R1 = w * b * (a + b / 2) / self.Length
        R2 = R1 + w * b
        x1a = lambda x: R1 * x
        x1 = lambda x: R1 * x ** 2
        y1 = integrate.quad(x1, 0, self.Length)
        y1a = integrate.quad(x1a, 0, self.Length)
        x2a = lambda x: R1 * x - R2 * (x - self.Length)
        x2 = lambda x: (R1 * x - R2 * (x - self.Length)) * x
        y2 = integrate.quad(x2, self.Length, self.Length + a)
        y2a = integrate.quad(x2a, self.Length, self.Length + a)
        x3a = lambda x: R1 * x - R2 * (x - self.Length) + w * (x - self.Length - a) ** 2 / 2
        x3 = lambda x: (R1 * x - R2 * (x - self.Length) + w * (x - self.Length - a) ** 2 / 2) * x
        y3 = integrate.quad(x3, self.Length + a, self.Length + a + b)
        y3a = integrate.quad(x3a, self.Length + a, self.Length + a + b)
        try:
            End_deflection = (self.Length - (y1[0]) / (y1a[0])) * (
                y1a[0])
            for i in range(round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):
                if i == 0:
                    Moment[i] = R1 * (i / self.sample_rate)
                    Shear[i] = R1
                elif i <= self.Length * self.sample_rate:
                    Moment[i] = R1 * (i / self.sample_rate)
                    Shear[i] = R1
                    y1 = integrate.quad(x1, 0, i / self.sample_rate)
                    y1a = integrate.quad(x1a, 0, i / self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                                i / self.sample_rate - (y1[0]) / (y1a[0])) * (y1a[0])) / (self.E * self.I) * 10 ** 6
                elif i > self.Length * self.sample_rate and i <= self.Length * self.sample_rate + a * self.sample_rate:
                    Moment[i] = R1 * (i / self.sample_rate) - R2 * (i / self.sample_rate - self.Length)
                    Shear[i] = R1 - R2
                    y1 = integrate.quad(x1, 0, self.Length)
                    y1a = integrate.quad(x1a, 0, self.Length)
                    y2 = integrate.quad(x2, self.Length, i / self.sample_rate)
                    y2a = integrate.quad(x2a, self.Length, i / self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                            i / self.sample_rate - (y1[0] + y2[0]) / (y1a[0] + y2a[0])) * (y1a[0] + y2a[0])) / (
                                            self.E * self.I) * 10 ** 6
                elif i > self.Length * self.sample_rate + a * self.sample_rate and i <= self.Length * self.sample_rate + a * self.sample_rate + b * self.sample_rate:
                    Moment[i] = R1 * (i / self.sample_rate) - R2 * (i / self.sample_rate - self.Length) + w * (
                                i / self.sample_rate - self.Length - a) ** 2 / 2
                    Shear[i] = R1 - R2 + w * (i / self.sample_rate - self.Length - a)
                    y1 = integrate.quad(x1, 0, self.Length)
                    y1a = integrate.quad(x1a, 0, self.Length)
                    y2 = integrate.quad(x2, self.Length, a + self.Length)
                    y2a = integrate.quad(x2a, self.Length, a + self.Length)
                    y3 = integrate.quad(x3, a + self.Length, i / self.sample_rate)
                    y3a = integrate.quad(x3a, a + self.Length, i / self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                            i / self.sample_rate - (y1[0] + y2[0] + y3[0]) / (y1a[0] + y2a[0] + y3a[0])) * (
                                             y1a[0] + y2a[0] + y3a[0])) / (self.E * self.I) * 10 ** 6
        except:
            pass
        for i in range(round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):
            self.Moment[i] += Moment[i]
            self.Shear[i] += -Shear[i]
            self.Deflection[i] += Deflection[i]
        self.R1 += -R1
        self.R2 += R2

    def Point_load(self, a, P):
        b = self.Length - a
        Moment = [0] * round(self.Length * self.sample_rate + 1)
        Shear = [0] * round(self.Length * self.sample_rate + 1)
        Deflection = [0] * round(self.Length * self.sample_rate + 1)
        Marea = [0] * round(self.Length * self.sample_rate + 1)
        R1 = P * b / self.Length
        R2 = P * a / self.Length
        for i in range(round(self.Length * self.sample_rate + 1)):
            if a == 0 or b == 0:
                Moment[i] = P * b * (i / self.sample_rate) / self.Length
                Shear[i] = R1
            elif b == 0:
                Moment[i] = P * a * (self.Length - i / self.sample_rate) / self.Length
                Shear[i] = -R2
            else:
                if i < a * self.sample_rate:
                    Moment[i] = P * b * (i / self.sample_rate) / self.Length
                    Shear[i] = R1
                    Deflection[i] = P * b * (i / self.sample_rate) / (6 * self.Length * self.E * self.I) * (
                                self.Length ** 2 - b ** 2 - (i / self.sample_rate) ** 2) * 10 ** 6
                else:
                    Moment[i] = P * a * (self.Length - i / self.sample_rate) / self.Length
                    Shear[i] = -R2
                    Deflection[i] = P * a * (self.Length - i / self.sample_rate) / (
                                6 * self.Length * self.E * self.I) * (
                                            self.Length * 2 * (i / self.sample_rate) - a ** 2 - (
                                                i / self.sample_rate) ** 2) * 10 ** 6
        for i in range(round(self.Length * self.sample_rate + 1)):
            self.Moment[i] += -Moment[i]
            self.Shear[i] += Shear[i]
            self.Deflection[i] += -Deflection[i]
            self.Marea[i] += Marea[i]
            try:

                for i in range(round(self.Length * self.sample_rate + 1),
                               round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):

                    self.Deflection[i] += (P*a*(self.Length - a)* (i/self.sample_rate - self.Length)) /(6*self.E * self.I * self.Length)*(self.Length + a)*10**3
            except:
                pass
        self.R1 += R1
        self.R2 += R2

    def Cantilever_Point_load(self, a, P):
        Moment = [0] * round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        Shear = [0] * round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        Deflection = [0] * round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        R1 = P * a / self.Length
        R2 = P / self.Length * (self.Length + a)
        try:
            for i in range(round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):
                if i <= self.Length * self.sample_rate:
                    Moment[i] = P * a * (i / self.sample_rate) / self.Length
                    Shear[i] = R1
                    Deflection[i] = (P * a * i / self.sample_rate) / (6 * self.Length * self.E * self.I) * (
                                self.Length ** 2 - (i / self.sample_rate) ** 2) * 10 ** 6
                elif i > self.Length * self.sample_rate and i <= self.Length * self.sample_rate + a * self.sample_rate:
                    Moment[i] = P * (a - (i / self.sample_rate - self.Length))
                    Shear[i] = -R2 + R1
                    Deflection[i] = -(P * (i / self.sample_rate - self.Length)) / (6 * self.E * self.I) * (
                                2 * a * self.Length + 3 * a * (i / self.sample_rate - self.Length) - (
                                    i / self.sample_rate - self.Length) ** 2) * 10 ** 6
                    Deflect = Deflection[i]
                elif i > self.Length * self.sample_rate + a * self.sample_rate:
                    Deflection[i] = (i / self.sample_rate - self.Length) / a * Deflect
        except:
            pass
        for i in range(round(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):
            self.Moment[i] += Moment[i]
            self.Shear[i] += -Shear[i]
            self.Deflection[i] += Deflection[i]
        self.R1 += -R1
        self.R2 += R2


class ORMType:
    def __init__(self, key, value):
        self.key = value

    def __repr__(self):
        return "ok"


class Beam:
    def __init__(self, Name, Length, sample_rate, E, Ix, Iy, Loadings, Cpe_V, Cpi_V, Cpe_H, Cpi_H, P, values,
                 Cantilever):
        self.Name = Name
        self.E = E
        self.Ix = Ix
        self.Iy = Iy
        self.Length = Length
        self.sample_rate = sample_rate
        self.Loadings = Loadings
        self.Cpe_V = Cpe_V
        self.Cpi_V = Cpi_V
        self.Cpe_H = Cpe_H
        self.Cpi_H = Cpi_H
        self.P = P
        self.values = values
        self.Cantilever = Cantilever
        self.MaxDeflections = {}
        self.Combined_Deflections ={}
        self.section_properties = {}

    def __str__(self):
        return self.Name

    def LoadCase(self, Loadings, Cantilever):

        for key in Loadings:
            if key == 'WOoP':
                setattr(self, key, Load_Case(self.Name, self.Length, self.sample_rate, self.E, self.Iy, self.values))
            else:
                setattr(self, key, Load_Case(self.Name, self.Length, self.sample_rate, self.E, self.Ix, self.values))
            # self.key = Load_Case(self.name, self.Length, self.sample_rate, self.E, self.I)
            # try:
            for j in range(len(Loadings[key][0])):
                getattr(self, key).UDL(Loadings[key][0][j][0], Loadings[key][0][j][1], Loadings[key][0][j][2])
                # self.key.UDL(Loadings[key][0][j][0],Loadings[key][0][j][1],Loadings[key][0][j][2])
            # except:
            # pass
            # try:
            for j in range(len(Loadings[key][1])):
                getattr(self, key).Point_load(Loadings[key][1][j][0], Loadings[key][1][j][1])
            # except:
            # pass
            if self.values['Cantilever'] == True:

                for j in range(len(Cantilever[key][0])):
                    getattr(self, key).Cantilever_UDL(Cantilever[key][0][j][0], Cantilever[key][0][j][1],
                                                      Cantilever[key][0][j][2])
                    # self.key.UDL(Loadings[key][0][j][0],Loadings[key][0][j][1],Loadings[key][0][j][2])
                # except:
                # pass
                # try:
                for j in range(len(Cantilever[key][1])):
                    getattr(self, key).Cantilever_Point_load(Cantilever[key][1][j][0], Cantilever[key][1][j][1])


def info(Name):
    try:
        with open(Name + '.pkl', 'rb') as inp:
            foo = pickle.load(inp)

    except (OSError, IOError) as e:
        print('No info')


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


def KeyCheck(Dictionary, string):
    try:
        output = Dictionary[string]
        if output == 'True':
            output = 1
        elif output == 'False':
            output = 0
        if output == '':
            output = 0
    except KeyError:
        output = 0
    return output
def KeyCheck1(Dictionary, string):
    try:
        output = Dictionary[string]
        if output == 'True':
            output = 1
        elif output == 'False':
            output = 0
        if output == '':
            output = 1
    except KeyError:
        output = 1
    return output

def Cantilever_Layout(variable):
    Layout = Layout = [[sg.Column([
        [sg.Text('Cantilever Length')],
    ]), sg.Column([
        [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'CLength'), enable_events=True, key='CLength'),sg.Button('Update', key = 'CLength_Update')],
    ])],
        [sg.Text('Number of UDL\'s'),
         sg.Input(default_text=int(float(KeyCheck(variable, 'CUDLS'))), key='CUDLS', enable_events=True, size=(5, 1))]
    ]
    try:
        CUDLS = int(float(KeyCheck(variable, 'CUDLS')))
    except:
        CUDLS = 1
    for i in range(CUDLS):
        i = str(i) + 'C'
        Layout += [[sg.Column([
            [sg.Text('Load Case')],
            [sg.Text('G1')],
            [sg.Text('G2')],
            [sg.Text('Q')],
            [sg.Text('Wup')],
            [sg.Text('Wdown')],
            [sg.Text('WOoP')],
        ]),
            sg.Column([
                [sg.Text('Magnitude KPa')],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1M' + str(i)), enable_events=True,
                          key='G1M' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2M' + str(i)), enable_events=True,
                          key='G2M' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QM' + str(i)), enable_events=True,
                          key='QM' + str(i))],
                [sg.Text()],
                [sg.Text()],
                [sg.Text()]
                # [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupM'+str(i)),enable_events=True,key='WupM'+str(i))],
                # [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownM'+str(i)),enable_events=True,key='WdownM'+str(i))],
                # [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPM'+str(i)),enable_events=True,key='WOoPM'+str(i))]
            ]),
            sg.Column([
                [sg.Text('Tributary Width m')],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1T' + str(i)), enable_events=True,
                          key='G1T' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2T' + str(i)), enable_events=True,
                          key='G2T' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QT' + str(i)), enable_events=True,
                          key='QT' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupT' + str(i)), enable_events=True,
                          key='WupT' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownT' + str(i)), enable_events=True,
                          key='WdownT' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPT' + str(i)), enable_events=True,
                          key='WOoPT' + str(i))]
            ]),
            sg.Column([
                [sg.Text('Start m')],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1a' + str(i)), enable_events=True,
                          key='G1a' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2a' + str(i)), enable_events=True,
                          key='G2a' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Qa' + str(i)), enable_events=True,
                          key='Qa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wupa' + str(i)), enable_events=True,
                          key='Wupa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wdowna' + str(i)), enable_events=True,
                          key='Wdowna' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPa' + str(i)), enable_events=True,
                          key='WOoPa' + str(i))],
            ]),
            sg.Column([
                [sg.Text('Stop m')],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1b' + str(i)), enable_events=True,
                          key='G1b' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2b' + str(i)), enable_events=True,
                          key='G2b' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Qb' + str(i)), enable_events=True,
                          key='Qb' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wupb' + str(i)), enable_events=True,
                          key='Wupb' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wdownb' + str(i)), enable_events=True,
                          key='Wdownb' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPb' + str(i)), enable_events=True,
                          key='WOoPb' + str(i))],
            ])]]
    Layout += [[sg.Text('New point Loads'),
                sg.Input(key='CPoint_loads', default_text=int(float(KeyCheck(variable, 'CPoint_loads'))),
                         enable_events=True, size=(5, 1))]]
    try:
        CPoint_loads = int(float(KeyCheck(variable, 'CPoint_loads')))
    except:
        CPoint_loads = 1
    for i in range(CPoint_loads):
        i = str(i) + 'C'
        Layout += [[sg.Column([
            [sg.Text('Load Case')],
            [sg.Text('G1')],
            [sg.Text('G2')],
            [sg.Text('Q')],
            [sg.Text('Wup')],
            [sg.Text('Wdown')],
            [sg.Text('WOoP')],
        ]),
            sg.Column([
                [sg.Text('Magnitude KN')],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1P' + str(i)), enable_events=True,
                          key='G1P' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2P' + str(i)), enable_events=True,
                          key='G2P' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QP' + str(i)), enable_events=True,
                          key='QP' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupP' + str(i)), enable_events=True,
                          key='WupP' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownP' + str(i)), enable_events=True,
                          key='WdownP' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPP' + str(i)), enable_events=True,
                          key='WOoPP' + str(i))]
            ]),
            sg.Column([
                [sg.Text('Location')],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1Pa' + str(i)), enable_events=True,
                          key='G1Pa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2Pa' + str(i)), enable_events=True,
                          key='G2Pa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QPa' + str(i)), enable_events=True,
                          key='QPa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupPa' + str(i)), enable_events=True,
                          key='WupPa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownPa' + str(i)), enable_events=True,
                          key='WdownPa' + str(i))],
                [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPPa' + str(i)), enable_events=True,
                          key='WOoPPa' + str(i))]
            ])]]

    Layout += [[sg.Text('Existing point loads'),
                sg.Input(key='CEx_Point_loads', default_text=int(float(KeyCheck(variable, 'CEx_Point_loads'))),
                         enable_events=True,
                         size=(5, 1))]]
    try:
        directory = KeyCheck(variable, 'Project')
    except:
        directory = os.fsencode(os.getcwd())
    existing = []
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".pkl"):
            filename = os.path.basename(filename)
            filename = os.path.splitext(filename)[0]
            if filename != KeyCheck(variable, 'Name'):
                existing += [filename]
            continue
        else:
            continue
    try:
        CEx_Point_loads = int(float(KeyCheck(variable, 'CEx_Point_loads')))
    except:
        CEx_Point_loads = 1
    for i in range(CEx_Point_loads):
        i = str(i) + 'C'
        try:
            Layout += [[sg.Combo(existing, size=(10, 1), key='Ex' + str(i),
                                 default_value=KeyCheck(variable, 'Ex' + str(i)), enable_events=True),
                        sg.Combo(['Left', 'Right'], size=(5, 1), key='L/R' + str(i),
                                 default_value=KeyCheck(variable, 'L/R' + str(i)), enable_events=True),
                        sg.Input(default_text=KeyCheck(variable, 'Exa' + str(i)), key='Exa' + str(i), size=(5, 1),
                                 enable_events=True)]]
        except:
            Ex_Point_loads = 0
    return Layout


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
def Layouts(variable):
    from openpyxl import load_workbook
    sg.theme('GrayGrayGray')
    wb = load_workbook(filename='Steel Design Calculator.xlsx')
    # SectionType = wb.sheetnames
    SectionSize = []
    # print(sheet_ranges['A6'].value)
    SectionType = ['Universal_Beam', 'Universal_Column', 'PFC', 'RHS', 'SHS', 'CHS', 'MGP 10', 'MGP 12']
    for j in wb['Universal_Beam']:
        SectionSize.append(j[0].value)
    SectionSize = list(filter(lambda item: item is not None, SectionSize))

    for j1 in range(4):
        SectionSize.pop(0)

    Layout = [[sg.Column([
        [sg.Text('Current Project'),
         sg.Input(key='Project', default_text=KeyCheck(variable, 'Project'), enable_events=True),
         sg.FolderBrowse('Change Project',enable_events=True, key='Project')],
        [sg.Column([
        [sg.Text('Name')],
        [sg.Text('Length')],
        [sg.Text('E')],
        [sg.Text('Choose Section Type:')],
        [sg.Text('Choose Section Size:')],
        [sg.Text('Seasoned')],
        [sg.Text('Duration of Load')],
        [sg.Text('Initial Moisture content:')],
        [sg.Text('Ix')],
        [sg.Text('Iy')]]),sg.Column([
        [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'Name'), enable_events=True, key='Name'),
         sg.Button('Reload Existing', key='Reload')],
        [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'Length'), enable_events=True, key='Length'),sg.Button('Update',key='Length_Update'),
         sg.Checkbox('Cantilever', default=KeyCheck(variable, 'Cantilever'), key='Cantilever', enable_events=True)],
        [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'E'), enable_events=True, key='E'),
         sg.Checkbox('Override for E', default=KeyCheck(variable, 'E_check'), key='E_check', enable_events=True)],
        [sg.Combo(SectionType, key='SectionType', enable_events=True, default_value=KeyCheck(variable, 'SectionType'),
                  size=(30, 1)),sg.Checkbox('Self Weight',default=KeyCheck(variable,'SelfWeight'),key='SelfWeight',enable_events=True)],
        # SectionType[0]
        [sg.Combo(SectionSize, key='SectionSize', enable_events=True, default_value=KeyCheck(variable, 'SectionSize'),
                  size=(30, 1))],  # SectionSize[0]
        [sg.Combo(['Seasoned', 'Unseasoned'], default_value='Seasoned', key='Seasoned')],
        [sg.Combo(['5 seconds', '5 minutes', '5 hours', '5 days', '5 months', '50+ years'],
                      default_value='5 months', key='load_duration')],
        [sg.Combo(['<15%', '>25%'], key='j2', default_value=KeyCheck(variable, 'j2'), size=(10, 1))],
        [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'Ix'), enable_events=True, key='Ix'),
         sg.Checkbox('Override for Ix', default=KeyCheck(variable, 'Ix_check'), key='Ix_check', enable_events=True)],
        [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'Iy'), enable_events=True, key='Iy'),
         sg.Checkbox('Override for Iy', default=KeyCheck(variable, 'Iy_check'), key='Iy_check', enable_events=True)]])]])]]
    Layout += steel_calculator()
    Layout += [[sg.Column([
        [sg.Text('Base Wind pressure (ULS):'),
         sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'P'), enable_events=True, key='P')],
        [sg.Text('Vertical Cp,e:'),
         sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpe_V'), enable_events=True, key='Cpe_V')],
        [sg.Text('Vertical Cp,i:'),
         sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpi_V'), enable_events=True, key='Cpi_V')],
        [sg.Text('Horizontal Cp,e:'),
         sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpe_H'), enable_events=True, key='Cpe_H')],
        [sg.Text('Horizontal Cp,i:'),
         sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpi_H'), enable_events=True, key='Cpi_H')],
        [sg.Text('Number of UDL\'s'),
         sg.Input(default_text=int(float(KeyCheck(variable, 'UDLS'))), key='UDLS', enable_events=True, size=(5, 1))]
    ])]]
    try:
        UDLS = int(float(KeyCheck(variable, 'UDLS')))
    except:
        UDLS = 1
    for i in range(UDLS):
        Layout += [[sg.Input(default_text=KeyCheck(variable, 'Load_case' + str(i)), enable_events=True,
                             key='Load_case' + str(i))], [sg.Column([
            [sg.Text('Load Case')],
            [sg.Text('G1')],
            [sg.Text('G2')],
            [sg.Text('Q')],
            [sg.Text('Wup')],
            [sg.Text('Wdown')],
            [sg.Text('WOoP')],
        ]),
                       sg.Column([
                           [sg.Text('Magnitude KPa')],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1M' + str(i)), enable_events=True,
                                     key='G1M' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2M' + str(i)), enable_events=True,
                                     key='G2M' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QM' + str(i)), enable_events=True,
                                     key='QM' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupM' + str(i)), enable_events=True,
                                     key='WupM' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownM' + str(i)),
                                     enable_events=True, key='WdownM' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPM' + str(i)), enable_events=True,
                                     key='WOoPM' + str(i))]
                           # [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupM'+str(i)),enable_events=True,key='WupM'+str(i))],
                           # [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownM'+str(i)),enable_events=True,key='WdownM'+str(i))],
                           # [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPM'+str(i)),enable_events=True,key='WOoPM'+str(i))]
                       ]),
                       sg.Column([
                           [sg.Text('Tributary Width m')],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1T' + str(i)), enable_events=True,
                                     key='G1T' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2T' + str(i)), enable_events=True,
                                     key='G2T' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QT' + str(i)), enable_events=True,
                                     key='QT' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupT' + str(i)), enable_events=True,
                                     key='WupT' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownT' + str(i)),
                                     enable_events=True, key='WdownT' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPT' + str(i)), enable_events=True,
                                     key='WOoPT' + str(i))]
                       ]),
                       sg.Column([
                           [sg.Text('Start m')],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1a' + str(i)), enable_events=True,
                                     key='G1a' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2a' + str(i)), enable_events=True,
                                     key='G2a' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Qa' + str(i)), enable_events=True,
                                     key='Qa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wupa' + str(i)), enable_events=True,
                                     key='Wupa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wdowna' + str(i)),
                                     enable_events=True, key='Wdowna' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPa' + str(i)), enable_events=True,
                                     key='WOoPa' + str(i))],
                       ]),
                       sg.Column([
                           [sg.Text('Stop m')],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1b' + str(i)), enable_events=True,
                                     key='G1b' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2b' + str(i)), enable_events=True,
                                     key='G2b' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Qb' + str(i)), enable_events=True,
                                     key='Qb' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wupb' + str(i)), enable_events=True,
                                     key='Wupb' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Wdownb' + str(i)),
                                     enable_events=True, key='Wdownb' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPb' + str(i)), enable_events=True,
                                     key='WOoPb' + str(i))],
                       ])]]
    Layout += [[sg.Text('New point Loads'),
                sg.Input(key='Point_loads', default_text=int(float(KeyCheck(variable, 'Point_loads'))),
                         enable_events=True, size=(5, 1))]]
    try:
        Point_loads = int(float(KeyCheck(variable, 'Point_loads')))
    except:
        Point_loads = 1
    for i in range(Point_loads):
        Layout += [[sg.Input(default_text=KeyCheck(variable, 'Load_caseP' + str(i)), enable_events=True,
                             key='Load_caseP' + str(i))], [sg.Column([
            [sg.Text('Load Case')],
            [sg.Text('G1')],
            [sg.Text('G2')],
            [sg.Text('Q')],
            [sg.Text('Wup')],
            [sg.Text('Wdown')],
            [sg.Text('WOoP')],
        ]),
                       sg.Column([
                           [sg.Text('Magnitude KN')],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1P' + str(i)), enable_events=True,
                                     key='G1P' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2P' + str(i)), enable_events=True,
                                     key='G2P' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QP' + str(i)), enable_events=True,
                                     key='QP' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupP' + str(i)), enable_events=True,
                                     key='WupP' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownP' + str(i)),
                                     enable_events=True, key='WdownP' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPP' + str(i)), enable_events=True,
                                     key='WOoPP' + str(i))]
                       ]),
                       sg.Column([
                           [sg.Text('Location')],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G1Pa' + str(i)), enable_events=True,
                                     key='G1Pa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'G2Pa' + str(i)), enable_events=True,
                                     key='G2Pa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'QPa' + str(i)), enable_events=True,
                                     key='QPa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WupPa' + str(i)), enable_events=True,
                                     key='WupPa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WdownPa' + str(i)),
                                     enable_events=True, key='WdownPa' + str(i))],
                           [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'WOoPPa' + str(i)),
                                     enable_events=True, key='WOoPPa' + str(i))]
                       ])]]

    Layout += [[sg.Text('Existing point loads'),
                sg.Input(key='Ex_Point_loads', default_text=int(float(KeyCheck(variable, 'Ex_Point_loads'))),
                         enable_events=True,
                         size=(5, 1))]]
    try:
        directory = KeyCheck(variable, 'Project')
        existing = []
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".pkl"):
                filename = os.path.basename(filename)
                filename = os.path.splitext(filename)[0]
                if filename != KeyCheck(variable, 'Name'):
                    existing += [filename]
                continue
            else:
                continue
    except:
        directory = os.fsencode(os.getcwd())
        existing = []
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".pkl"):
                filename = os.path.basename(filename)
                filename = os.path.splitext(filename)[0]
                if filename != KeyCheck(variable, 'Name'):
                    existing += [filename]
                continue
            else:
                continue
    try:
        Ex_Point_loads = int(float(KeyCheck(variable, 'Ex_Point_loads')))
    except:
        Ex_Point_loads = 1
    for i in range(Ex_Point_loads):
        try:
            Layout += [[sg.Combo(existing, size=(10, 1), key='Ex' + str(i),
                                 default_value=KeyCheck(variable, 'Ex' + str(i)), enable_events=True),
                        sg.Combo(['Left', 'Right'], size=(5, 1), key='L/R' + str(i),
                                 default_value=KeyCheck(variable, 'L/R' + str(i)), enable_events=True),
                        sg.Input(default_text=KeyCheck(variable, 'Exa' + str(i)), key='Exa' + str(i), size=(5, 1),
                                 enable_events=True)]]
        except:
            Ex_Point_loads = 0

    Layout += [[sg.Column([
        [sg.Button('Calculate', key='Calculate')],
        [sg.Text('Number of Load Cases:'),sg.Input(key='Load_Cases', default_text=int(float(KeyCheck(variable, 'Load_Cases'))),size=(5,1)),sg.Button('Default Load Cases',key='Default_Load_Cases')],
    ])]]
    Layout +=[[sg.Text('',key='P/F_Check')]]
    try:
        Load_Cases = int(float(KeyCheck(variable, 'Load_Cases')))
    except:
        Load_Cases = 1

    Canvas = []
    for i in range(Load_Cases):
        #try:
            COL1 = [[sg.Text('Load Case')],[sg.Input(key='Load_Case'+str(i), default_text=KeyCheck(variable, 'Load_Case'+str(i)),size=(20,1),enable_events=True)]]
            COL2 = [[sg.Text('G1')],
                [sg.Input(key='Load_CaseG1' + str(i), default_text=KeyCheck(variable, 'Load_CaseG1' + str(i)),size=(5,1),enable_events=True)]]
            COL3 = [[sg.Text('G2')],
                [sg.Input(key='Load_CaseG2' + str(i), default_text=KeyCheck(variable, 'Load_CaseG2' + str(i)),size=(5,1),enable_events=True)]]
            COL4 = [[sg.Text('Q')],
                [sg.Input(key='Load_CaseQ' + str(i), default_text=KeyCheck(variable, 'Load_CaseQ' + str(i)),size=(5,1),enable_events=True)]]
            COL5 = [[sg.Text('Wup')],
                [sg.Input(key='Load_CaseWup' + str(i), default_text=KeyCheck(variable, 'Load_CaseWup' + str(i)),size=(5,1),enable_events=True)]]
            COL6 = [[sg.Text('Wdown')],
                [sg.Input(key='Load_CaseWdown' + str(i), default_text=KeyCheck(variable, 'Load_CaseWdown' + str(i)),size=(5,1),enable_events=True)]]
            COL7 = [[sg.Text('WOoP')],
                [sg.Input(key='Load_CaseWOoP' + str(i), default_text=KeyCheck(variable, 'Load_CaseWOoP' + str(i)),size=(5,1),enable_events=True)]]
            COL8 = [[sg.Text('Moment')],
                [sg.Checkbox('Moment',key='Load_CaseM' + str(i), default=bool(float(KeyCheck(variable, 'Load_CaseM' + str(i)))),enable_events=True)]]
            COL9 = [[sg.Text('Shear')],
                [sg.Checkbox('Shear',key='Load_CaseS' + str(i), default=bool(float(KeyCheck(variable, 'Load_CaseS' + str(i)))),enable_events=True)]]
            COL10 = [[sg.Text('Deflection')],
                [sg.Checkbox('Deflection',key='Load_CaseD' + str(i), default=bool(float(KeyCheck(variable, 'Load_CaseD' + str(i)))),enable_events=True)]]
            COL11 = [[sg.Text('Deflection Limit L/')],[sg.Input(default_text=KeyCheck1(variable,'Load_CaseD_L'+str(i)),size=(5,1), key='Load_CaseD_L' + str(i), enable_events=True)]]

            Canvas += [
                [sg.Column(COL1), sg.Column(COL2), sg.Column(COL3), sg.Column(COL4), sg.Column(COL5), sg.Column(COL6),
                 sg.Column(COL7), sg.Column(COL8), sg.Column(COL9), sg.Column(COL10),sg.Column(COL11)]]
            Canvas += [[sg.Canvas(key='controls_cv'+str(i))],[sg.Canvas(key='Graph'+str(i))]]
        #except:
           # Load_Cases = 0;

    if KeyCheck(variable, 'Cantilever') == True:
        Layout = [[sg.Column([[sg.Column(Layout, vertical_alignment='t'), sg.VSeparator(),
                               sg.Column(Cantilever_Layout(variable), vertical_alignment='t'), sg.VSeparator(),
                               sg.Column(
                                   Canvas,
                                   vertical_alignment='t')]], scrollable=True, expand_x=True,
                             expand_y=True,key ='Column')]]
    else:
        Layout = [[sg.Column([[sg.Column(Layout, vertical_alignment='t'), sg.VSeparator(), sg.Column(
            Canvas, vertical_alignment='t')]], scrollable=True, expand_x=True,
                             expand_y=True,key ='Column')]]
    window = sg.Window('Window', Layout, resizable=True).finalize()
    window.Maximize()
    from openpyxl import load_workbook

    wb = load_workbook(filename='Steel Design Calculator.xlsx')
    # SectionType = wb.sheetnames
    SectionSize = []
    # print(sheet_ranges['A6'].value)
    SectionType = ['Universal_Beam', 'Universal_Column', 'PFC', 'RHS', 'SHS', 'CHS', 'T-section']
    for j in wb['Universal_Beam']:
        SectionSize.append(j[0].value)
    SectionSize = list(filter(lambda item: item is not None, SectionSize))
    print(SectionSize)
    for j1 in range(4):
        SectionSize.pop(0)
    fig_canvas_agg = None
    while True:
        event, values = window.read()
        if event == sg.WINDOW_CLOSED or event == 'Quit':
            break
        for i in values:
            if event == i:
                variable_write(values, 'Loading')
        if event == 'UDLS' or event == 'Point_loads' or event == 'Ex_Point_loads' or event == 'CUDLS' or event == 'CPoint_loads' or event == 'CEx_Point_loads' or event == 'Cantilever':
            window.close()
            variable = variables('Loading')
            Layouts(variable)
        if event == 'Default_Load_Cases':

            variable_write(values, 'Loading')
            window.close()

            variable = variables('Loading')
            for i in variable:
                if 'Load_Case' in i:
                    variable[i] = '0'
            variable['Load_Cases'] ='10'
            variable['Load_Case0'] = '1.2G2 + 1.5Q'
            variable['Load_CaseG20']='1.2'
            variable['Load_CaseQ0']='1.5'
            variable['Load_CaseM0'] = '1'
            variable['Load_CaseS0'] = '1'

            variable['Load_Case1'] = '1.35G2'
            variable['Load_CaseG21']='1.35'
            variable['Load_CaseM1'] = '1'
            variable['Load_CaseS1'] = '1'

            variable['Load_Case2'] = '0.9G1 + 1.5Wup'
            variable['Load_CaseG12']='0.9'
            variable['Load_CaseWup2']='1.5'
            variable['Load_CaseM2'] = '1'
            variable['Load_CaseS2'] = '1'

            variable['Load_Case3'] = '1.2G2 + 1.5Wdown'
            variable['Load_CaseG23'] = '1.2'
            variable['Load_CaseWdown3'] = '1.5'
            variable['Load_CaseM3'] = '1'
            variable['Load_CaseS3'] = '1'

            variable['Load_Case4'] = '1.5WOoP'
            variable['Load_CaseWOoP4']='1.5'
            variable['Load_CaseM4'] = '1'
            variable['Load_CaseS4'] = '1'

            variable['Load_Case5'] = 'G2'
            variable['Load_CaseG25']='1'
            variable['Load_CaseD5'] = '1'
            variable['Load_CaseD_L5'] = '360'

            variable['Load_Case6'] = 'Q'
            variable['Load_CaseQ6']='1'
            variable['Load_CaseD6'] = '1'
            variable['Load_CaseD_L6'] = '250'

            variable['Load_Case8'] = 'Wdown'
            variable['Load_CaseWdown8']='1'
            variable['Load_CaseD8'] = '1'
            variable['Load_CaseD_L8'] = '150'

            variable['Load_Case7'] = 'Wup'
            variable['Load_CaseWup7']='1'
            variable['Load_CaseD7'] = '1'
            variable['Load_CaseD_L7'] = '150'

            variable['Load_Case9'] = 'WOoP'
            variable['Load_CaseWOoP9']='1'
            variable['Load_CaseD9'] = '1'
            variable['Load_CaseD_L9'] = '150'
            Layouts(variable)
        if event == 'Length_Update':
            variable_write(values, 'Loading')
            for i in range(int(variable['UDLS'])):
                window['segment length'].update(values['Length'])
                window['G1b' + str(i)].update(values['Length'])
                window['G2b' + str(i)].update(values['Length'])
                window['Qb' + str(i)].update(values['Length'])
                window['Wupb' + str(i)].update(values['Length'])
                window['Wdownb' + str(i)].update(values['Length'])
                window['WOoPb' + str(i)].update(values['Length'])
        if event == 'CLength_Update':
            variable_write(values, 'Loading')
            for i in range(int(variable['CUDLS'])):
                i = str(i) + 'C'
                window['G1b' + str(i)].update(values['Length'])
                window['G2b' + str(i)].update(values['Length'])
                window['Qb' + str(i)].update(values['Length'])
                window['Wupb' + str(i)].update(values['Length'])
                window['Wdownb' + str(i)].update(values['Length'])
                window['WOoPb' + str(i)].update(values['Length'])
        if event == 'calculate' or event == 'Calculate':
           try:
                SectionType1 = values['SectionType']
                SectionSize1 = values['SectionSize']

                alpha_m = float(values['alpha_m'])
                restraint = values['restraint']
                load_height_position = values['load_height_position']
                longitudinal_position = values['longitudinal_position']
                ends_with_restraint = values['ends_with_restraint']
                section_properties = st.section_properties(SectionType1, SectionSize1, 0, 0, 0)
                st.compact(section_properties, SectionType1, float(values['segment length']), alpha_m, restraint,
                           load_height_position,
                           longitudinal_position, ends_with_restraint)
                st.shear(section_properties, SectionType1)
                st.shear_moment(section_properties, 0)
                st.axial_compression(section_properties, SectionType1, float(values['segment length']))

                print(section_properties)
                try:
                    del Timber12
                    del Timber13
                except:
                    pass
                try:
                    del section_properties1
                except:
                    pass
                section_properties1 = section_properties
                # End program if user closes window or
                # presses the OK button
                # if event == 'quit' or event == sg.WIN_CLOSED:

                window['PhiMsx'].update(str(round(section_properties['PhiMsx'], 1)) + 'kNm')

                window['PhiMbx'].update(str(round(section_properties['PhiMbx'], 1)) + 'kNm')

                window['PhiMsy'].update(str(round(section_properties['PhiMsy'], 1)) + 'kNm')
                window['PhiMby'].update(str(round(section_properties['PhiMby'], 1)) + 'kNm')
                window['PhiNsx'].update(str(round(section_properties['PhiNsx'], 1)) + 'kN')
                window['PhiNcx'].update(str(round(section_properties['PhiNcx'], 1)) + 'kN')
                window['PhiNsy'].update(str(round(section_properties['PhiNsy'], 1)) + 'kN')
                window['PhiNcy'].update(str(round(section_properties['PhiNcy'], 1)) + 'kN')
                window['PhiVu'].update(str(round(section_properties['PhiVu'], 2)) + ' KN')
           except:
            MGP10 = {'90x90': [90, 90], '140x90': [140, 90], '190x90': [190, 90], '240x90': [240, 90],
                        '300x90': [300, 90], '90x45': [90, 45], '140x45': [140, 45], '190x45': [190, 45],
                        '240x45': [240, 45]}
            Timber = MGP10[values['SectionSize']]
            if values['Seasoned'] == 'Seasoned':
                Seasoned = True
            else:
                Seasoned = False
            try:
                del section_properties1
            except:
                pass
            try:
                del Timber12
                del Timber13
            except:
                pass
            Timber12 = Timber1.beam_moment(Timber[1],Timber[0],values['load_duration'],Seasoned,1,1,1,float(values['segment length'])*1000,values['SectionType'],True)
            Timber13 = Timber1.beam_shear(Timber[1],Timber[0],values['load_duration'],Seasoned,1,1,1,float(values['segment length'])*1000,values['SectionType'])
        if event == 'Calculate':
            variable_write(values, 'Loading')
            if fig_canvas_agg:
                delete_figure_agg(fig_canvas_agg)
            plt.close('all')
            G1U = []
            G2U = []
            QU = []
            WupU = []
            WdownU = []
            WOoPU = []
            G1P = []
            G2P = []
            QP = []
            WupP = []
            WdownP = []
            WOoPP = []
            for i in values:
                try:
                    values[i] = float(values[i])
                except:
                    print(values[i])
            for i in range(UDLS):
                G1U += [
                    [values['G1a' + str(i)], values['G1b' + str(i)], values['G1M' + str(i)] * values['G1T' + str(i)]]]
                G2U += [
                    [values['G2a' + str(i)], values['G2b' + str(i)], values['G2M' + str(i)] * values['G2T' + str(i)]]]
                QU += [[values['Qa' + str(i)], values['Qb' + str(i)], values['QM' + str(i)] * values['QT' + str(i)]]]

                WupU += [[values['Wupa' + str(i)], values['Wupb' + str(i)], values['WupM' + str(i)] * values['WupT' + str(i)]]]

                WdownU += [[values['Wdowna' + str(i)], values['Wdownb' + str(i)], values['WdownM' + str(i)] * values['WdownT' + str(i)]]]

                WOoPU += [[values['WOoPa' + str(i)], values['WOoPb' + str(i)], values['WOoPM' + str(i)] * values['WOoPT' + str(i)]]]


            for i in range(Point_loads):
                G1P += [[values['G1Pa' + str(i)], values['G1P' + str(i)]]]
                G2P += [[values['G2Pa' + str(i)], values['G2P' + str(i)]]]
                QP += [[values['QPa' + str(i)], values['QP' + str(i)]]]
                WupP += [[values['WupPa' + str(i)], values['WupP' + str(i)]]]
                WdownP += [[values['WdownPa' + str(i)], values['WdownP' + str(i)]]]
                WOoPP += [[values['WOoPPa' + str(i)], values['WOoPP' + str(i)]]]
            for i in range(Ex_Point_loads):
                Beam = open_file(values['Ex' + str(i)], values['Project'])
                if values['L/R' + str(i)] == 'Left':
                    G1P += [[values['Exa' + str(i)], Beam.G1.R1]]
                    G2P += [[values['Exa' + str(i)], Beam.G2.R1]]
                    QP += [[values['Exa' + str(i)], Beam.Q.R1]]
                    WupP += [[values['Exa' + str(i)], Beam.Wup.R1]]
                    WdownP += [[values['Exa' + str(i)], Beam.Wdown.R1]]
                    WOoPP += [[values['Exa' + str(i)], Beam.WOoP.R1]]
                else:
                    G1P += [[values['Exa' + str(i)], Beam.G1.R2]]
                    G2P += [[values['Exa' + str(i)], Beam.G2.R2]]
                    QP += [[values['Exa' + str(i)], Beam.Q.R2]]
                    WupP += [[values['Exa' + str(i)], Beam.Wup.R2]]
                    WdownP += [[values['Exa' + str(i)], Beam.Wdown.R2]]
                    WOoPP += [[values['Exa' + str(i)], Beam.WOoP.R2]]
            Loadings = {'G1': [G1U, G1P], 'G2': [G2U, G2P], 'Q': [QU, QP], 'Wup': [WupU, WupP],
                        'Wdown': [WdownU, WdownP], 'WOoP': [WOoPU, WOoPP]}
            if values['Ix_check'] == False:
                try:
                    SectionType1 = values['SectionType']
                    SectionSize1 = values['SectionSize']
                    section_properties = st.section_properties(SectionType1, SectionSize1, 0, 0, 0)
                    Ix = section_properties['Ix'] * 10 ** -12
                except:
                    MGP10 = {'90x90': [90, 90], '140x90': [140, 90], '190x90': [190, 90], '240x90': [240, 90],
                             '300x90': [300, 90], '90x45': [90, 45], '140x45': [140, 45], '190x45': [190, 45],
                             '240x45': [240, 45]}
                    Timber = MGP10[values['SectionSize']]
                    Ix = Timber[1] * Timber[0] ** 3 / 12 * 10 ** -12
            else:
                Ix = values['Ix'] * 10 ** -12
            if values['Iy_check'] == False:
                try:
                    SectionType1 = values['SectionType']
                    SectionSize1 = values['SectionSize']
                    section_properties = st.section_properties(SectionType1, SectionSize1, 0, 0, 0)
                    Iy = section_properties['Iy'] * 10 ** -12
                except:
                    MGP10 = {'90x90': [90, 90], '140x90': [140, 90], '190x90': [190, 90], '240x90': [240, 90],
                             '300x90': [300, 90],'90x45':[90,45],'140x45':[140,45],'190x45':[190,45],'240x45':[240,45]}
                    Timber = MGP10[values['SectionSize']]
                    Iy = Timber[0] * Timber[1] ** 3 / 12 * 10 ** -12
            else:
                Iy = values['Iy'] * 10 ** -12

            CG1U = []
            CG2U = []
            CQU = []
            CWupU = []
            CWdownU = []
            CWOoPU = []
            CG1P = []
            CG2P = []
            CQP = []
            CWupP = []
            CWdownP = []
            CWOoPP = []
            try:
                CUDLS = int(float(KeyCheck(variable, 'CUDLS')))
            except:
                CUDLS = 1
            try:
                CPoint_loads = int(float(KeyCheck(variable, 'CPoint_loads')))
            except:
                CPoint_loads = 1
            try:
                CEx_Point_loads = int(float(KeyCheck(variable, 'CEx_Point_loads')))
            except:
                CEx_Point_loads = 1
            for i in values:
                try:
                    values[i] = float(values[i])
                except:
                    print(values[i])
            try:
                for i in range(CUDLS):
                    i = str(i) + 'C'
                    CG1U += [[values['G1a' + str(i)], values['G1b' + str(i)],
                              values['G1M' + str(i)] * values['G1T' + str(i)]]]
                    CG2U += [[values['G2a' + str(i)], values['G2b' + str(i)],
                              values['G2M' + str(i)] * values['G2T' + str(i)]]]
                    CQU += [
                        [values['Qa' + str(i)], values['Qb' + str(i)], values['QM' + str(i)] * values['QT' + str(i)]]]
                    CWupM = (values['Cpe_V'] + values['Cpi_V']) * values['P']
                    CWupU += [[values['Wupa' + str(i)], values['Wupb' + str(i)], CWupM * values['WupT' + str(i)]]]
                    CWdownM = (values['Cpe_V'] + values['Cpi_V']) * values['P']
                    CWdownU += [
                        [values['Wdowna' + str(i)], values['Wdownb' + str(i)], CWdownM * values['WdownT' + str(i)]]]
                    CWOoPM = (values['Cpe_H'] + values['Cpi_H']) * values['P']
                    CWOoPU += [[values['WOoPa' + str(i)], values['WOoPb' + str(i)], CWOoPM * values['WOoPT' + str(i)]]]
                for i in range(CPoint_loads):
                    i = str(i) + 'C'
                    CG1P += [[values['G1Pa' + str(i)], values['G1P' + str(i)]]]
                    CG2P += [[values['G2Pa' + str(i)], values['G2P' + str(i)]]]
                    CQP += [[values['QPa' + str(i)], values['QP' + str(i)]]]
                    CWupP += [[values['WupPa' + str(i)], values['WupP' + str(i)]]]
                    CWdownP += [[values['WdownPa' + str(i)], values['WdownP' + str(i)]]]
                    CWOoPP += [[values['WOoPPa' + str(i)], values['WOoPP' + str(i)]]]
                for i in range(CEx_Point_loads):
                    i = str(i) + 'C'
                    Beam = open_file(values['Ex' + str(i)], values['Project'])
                    if values['L/R' + str(i)] == 'Left':
                        CG1P += [[values['Exa' + str(i)], Beam.G1.R1]]
                        CG2P += [[values['Exa' + str(i)], Beam.G2.R1]]
                        CQP += [[values['Exa' + str(i)], Beam.Q.R1]]
                        CWupP += [[values['Exa' + str(i)], Beam.Wup.R1]]
                        CWdownP += [[values['Exa' + str(i)], Beam.Wdown.R1]]
                        CWOoPP += [[values['Exa' + str(i)], Beam.WOoP.R1]]
                    else:
                        CG1P += [[values['Exa' + str(i)], Beam.G1.R2]]
                        CG2P += [[values['Exa' + str(i)], Beam.G2.R2]]
                        CQP += [[values['Exa' + str(i)], Beam.Q.R2]]
                        CWupP += [[values['Exa' + str(i)], Beam.Wup.R2]]
                        CWdownP += [[values['Exa' + str(i)], Beam.Wdown.R2]]
                        CWOoPP += [[values['Exa' + str(i)], Beam.WOoP.R2]]
                Cantilever = {'G1': [CG1U, CG1P], 'G2': [CG2U, CG2P], 'Q': [CQU, CQP], 'Wup': [CWupU, CWupP],
                              'Wdown': [CWdownU, CWdownP], 'WOoP': [CWOoPU, CWOoPP]}
            except:
                Cantilever = 0
            if values['E_check'] == True:
                E = values['E'] * 10 ** 9
            else:
                if values['SectionType'] == 'MGP 10':
                    E = 10 * 10 ** 9
                elif values['SectionType'] == 'MGP 12':
                    E = 12.7 * 10 ** 9
                else:
                    E = 200 * 10 ** 9
            Method(values['Name'], values['Length'], 100, E, Ix, Iy, Loadings,
                   values['Cpe_V'], values['Cpi_V'], values['Cpe_H'], values['Cpi_H'], values['P'], values, Cantilever)

            Beam1 = open_file(values['Name'], values['Project'])
            try:
                getattr(Beam1,'section_properties')[0] = section_properties1
            except:
                pass
            try:
                getattr(Beam1,'section_properties')[0] = Timber12
                getattr(Beam1, 'section_properties')[1] = Timber13
            except:
                pass
            # fig_canvas_agg = draw_figure(window['ULSGraph'].TKCanvas,draw_graph(Beam1))
            for i in range(Load_Cases):

                draw_figure_w_toolbar(window['Graph'+str(i)].TKCanvas, draw_graph(Beam1,i), window['controls_cv'+str(i)].TKCanvas)

            window['P/F_Check'].update(PassFail(Beam1))
            window.refresh()
            window['Column'].contents_changed()
            Output.Main(Beam1)
        if event == 'Reload':
            try:
                directory = KeyCheck(variable, 'Project')
                existing = []
                for file in os.listdir(directory):
                    filename = os.fsdecode(file)
                    if filename.endswith(".pkl"):
                        filename = os.path.basename(filename)
                        filename = os.path.splitext(filename)[0]
                        existing += [filename]
                        continue
                    else:
                        continue
            except:
                directory = os.fsencode(os.getcwd())
            Layout1 = [
                [sg.Text('Select Existing Beam')],
                [sg.Combo(existing, key='Reload_Ex', enable_events=True)],
                [sg.Button('Continue', key='Reload_continue')]
            ]
            Project = values['Project']
            window1 = sg.Window('Select Existing Beam', Layout1, resizable=True).finalize()
            window.close()
            while True:
                event, values = window1.read()
                if event == sg.WINDOW_CLOSED or event == 'Quit':
                    break
                elif event == 'Reload_continue':
                    Beam = open_file(values['Reload_Ex'], KeyCheck(variable, 'Project'))
                    window1.close()
                    Beam.values['Project'] = Project
                    Layouts(Beam.values)
                    break

        if event == 'dLimit':
            for x in values['dLimit'].split('/'):
                if x.isdigit():
                    Limit = int(x)
            window['dLimit1'].update(str(round(float(values['Length']) * 1000 / Limit, 2)) + ' mm')
            if float(values['Length']) * 1000 / Limit > d[0] + d[1]:
                window['dCheck'].update('OK')
            else:
                window['dCheck'].update('NG')
        if event == 'SectionType':
            if values['SectionType'] == 'Universal_Beam' or values['SectionType'] == 'Universal_Column' or values[
                'SectionType'] == 'PFC' or values['SectionType'] == 'RHS' or values['SectionType'] == 'SHS' or values[
                'SectionType'] == 'CHS':
                item = values[event]
                matrix = []
                SectionNames = wb[item]
                for x in SectionNames:
                    matrix.append(x[0].value)
                for i in range(7):
                    matrix.pop(0)
                SectionSize = list(filter(lambda item: item is not None, SectionSize))
                window['SectionSize'].update(value=matrix[7], values=matrix)
            elif values['SectionType'] == 'MGP 10' or values['SectionType'] == 'MGP 12':
                SectionSize = ['90x45','90x90','140x45', '140x90', '190x45','190x90', '240x45','240x90']
                window['SectionSize'].update(value='90x90', values=SectionSize)

        # window['b1'].update('OK')
        elif event == 'b2' or event == sg.WIN_CLOSED:
            break
        elif event == 'restraint':
            if values['restraint'] == 'FU' or values['restraint'] == 'PU':
                window['ends_with_restraint'].update(values=['Any'], value='Any')
            elif values['restraint'] == 'FF' or values['restraint'] == 'FP' or values['restraint'] == 'PP':
                window['ends_with_restraint'].update(values=['None', 'One', 'Both'], value='None')
            elif values['restraint'] == 'FL' or values['restraint'] == 'PL' or values['restraint'] == 'LL':
                window['ends_with_restraint'].update(values=['None'], value='None')

        elif event == 'print_calcs':
            SectionType1 = values['SectionType']
            SectionSize1 = values['SectionSize']
            Length = float(values['Length'])
            alpha_m = float(values['alpha_m'])
            restraint = values['restraint']
            load_height_position = values['load_height_position']
            longitudinal_position = values['longitudinal_position']
            ends_with_restraint = values['ends_with_restraint']
            section_properties = st.section_properties(SectionType1, SectionSize1, 0, 0, 0)
            calculations = st.compact(section_properties, SectionType1, Length, alpha_m, restraint,
                                      load_height_position,
                                      longitudinal_position, ends_with_restraint)
            st.printcalcs(SectionSize1, section_properties, values['print_location'], values['job_name'], '', 0, 0, 0,
                          True)
    window.close()


def steel_calculator():
    # load in data to be read
    from openpyxl import load_workbook

    wb = load_workbook(filename='Steel Design Calculator.xlsx')
    # SectionType = wb.sheetnames
    SectionSize = []
    # print(sheet_ranges['A6'].value)
    SectionType = ['Universal_Beam', 'Universal_Column', 'PFC', 'RHS', 'SHS', 'CHS', 'T-section']
    for j in wb['Universal_Beam']:
        SectionSize.append(j[0].value)
    SectionSize = list(filter(lambda item: item is not None, SectionSize))
    print(SectionSize)
    for j1 in range(4):
        SectionSize.pop(0)
    print(SectionSize)
    # load in GUI library

    sg.theme('GrayGrayGray')
    # ttk_style = 'vista'
    layout = [
        [sg.Column([
            [sg.Text('Input segment Length in metres below:')],
            [sg.Input(default_text=KeyCheck(variable, 'segment length'), key='segment length', size=(5, 1))],
            [sg.Text('Input alpha m below:', key='t3')],
            [sg.Input(key='alpha_m', size=(5, 1), default_text=KeyCheck(variable, 'alpha_m'))],
            [sg.Combo(['FF', 'FP', 'FL', 'PP', 'PL', 'LL', 'FU', 'PU'], key='restraint', default_value='FF',
                      enable_events=True)],
            [sg.Combo(['Shear centre', 'Top flange'], key='load_height_position', default_value='Shear centre',
                      size=(15, 1))],
            [sg.Combo(['Within segment', 'At segment end'], key='longitudinal_position', default_value='Within segment',
                      size=(15, 1)), sg.Text('Load height')],
            [sg.Combo(['Any', 'None', 'One', 'Both'], key='ends_with_restraint', default_value='One'),
             sg.Text('Ends with Lateral restraint')],
            [sg.Button('Calculate Steel', key='calculate', use_ttk_buttons=True)],
            [sg.Column([
                [sg.Text('\u03A6Msx = '), sg.Text('', key='PhiMsx')],
                [sg.Text('\u03A6Mbx = '), sg.Text('', key='PhiMbx')],
                [sg.Text('\u03A6Msy = '), sg.Text('', key='PhiMsy')],
                [sg.Text('\u03A6Mby = '), sg.Text('', key='PhiMby')]
            ]), sg.Column([
                [sg.Text('\u03A6Nsx = '), sg.Text('', key='PhiNsx')],
                [sg.Text('\u03A6Ncx = '), sg.Text('', key='PhiNcx')],
                [sg.Text('\u03A6Nsy = '), sg.Text(key='PhiNsy')],
                [sg.Text('\u03A6Ncy = '), sg.Text(key='PhiNcy')]

            ]), sg.Column([
                [sg.Text('\u03A6Vu = '), sg.Text('', key='PhiVu')],
                [sg.Text('\u03A6Vvm = '), sg.Text('', key='PhiVvm')]
            ])]])]]

    # Create the window
    return layout

    # Create an event loop


variable = variables('Loading')
Layouts(variable)


x = np.linspace(0, float(Beam.Length + Beam.values['CLength']),
                int(Beam.Length * Beam.sample_rate + Beam.values['CLength'] * Beam.sample_rate + 1))
plt.plot(x, Beam.Q.Shear)
plt.show()
