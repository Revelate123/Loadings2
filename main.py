#This script will calculate design actions for use in simple beam problems

#It will determine the moment and shear of a simply supported beam
#It will determine the required Second moment of inertia to meet deflection criteria
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import pickle as pickle
import PySimpleGUI as sg
import csv
import steel_functions as st
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def Loads(G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit,Gpoint,Qpoint,Wpointsls, a):

    wdown12G15Q = 1.2 * G * Trib_width + 1.5 * Q * Trib_width
    wdown135G = 1.35*G*Trib_width
    wdown12GWuls = 1.2*G*Trib_width + Wsls*1.48*Trib_width
    wup = min(0.9* G * Trib_width - Wsls*1.48*Trib_width,0)
    wdown = max(wdown135G,wdown12G15Q,wdown12GWuls)
    Mulsdown = wdown*Length**2/8
    Vulsdown = wdown*Length/2
    Mulsup = wup*Length**2/8
    Vulsup = wup * Length/2

    w12G15Q = 1.2*Gpoint + 1.5*Qpoint
    w135G = 1.35*Gpoint
    w12GWuls = 1.2*Gpoint + 1.48*Wpointsls
    wdownpoint = max(w12G15Q,w135G,w12GWuls)
    wuppoint = min(0.9*Gpoint - 1.48*Wpointsls,0)

    Mulsdownpoint = wdownpoint*a*(Length-a)/Length
    Mulsuppoint = wuppoint*a*(Length-a)/Length
    x = min(a,(Length-a))
    y = max(wdownpoint,-1*wuppoint)
    Vpoint = y*x/Length

    GI = 5*G*Trib_width*10**3*Length**4 / (384*200*10**9*Glimit*10**-3)*10**12 + Gpoint*10**3*Length*10**3/(48*200**9*Glimit*10**-3)
    QI = 5 * 0.7 * Q * Trib_width * 10 ** 3 * Length ** 4 / (384 * 200 * 10 ** 9 * Qlimit * 10 ** -3) * 10 ** 12 + Qpoint*10**3*Length*10**3/(48*200**9*Qlimit*10**-3)
    WI = 5 * Wsls * Trib_width * 10 ** 3 * Length ** 4 / (384 * 200 * 10 ** 9 * Wlimit * 10 ** -3) * 10 ** 12 + Wpointsls*10**3*Length*10**3/(48*200**9*Wlimit*10**-3)

    return Mulsdown, max(Vulsdown,Vulsup), Mulsup, GI, QI, WI, wdown,wup,wdown12G15Q,wdown135G,wdown12GWuls,Mulsdownpoint,Mulsuppoint,Vpoint,wdownpoint,wuppoint,Mulsdownpoint,Mulsuppoint,Vpoint

def Write_Loads(Name,G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit,Gpoint,Qpoint,Wpointsls, a):
    Load = Loads(G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit,Gpoint,Qpoint,Wpointsls, a)
    with open(Name + '.txt', 'w') as f:
        f.write(Name)
        f.write('\n')
        f.write('\n')
        f.write('Length:  '+ str(Length))
        f.write('\nG:  ' + str(G))
        f.write('\nQ:  '+str(Q))
        f.write('\nWsls:  '+str(Wsls))
        f.write('\nTrib Width:  ' + str(Trib_width))
        f.write('\nG Point:  '+ str(Gpoint))
        f.write('\nQ Point:  ' + str(Qpoint))
        f.write('\nWsls Point:  ' + str(Wpointsls))
        f.write('\n1.2G + 1.5Q:  ' + str(Load[8]) +' kN/m')
        f.write('\n1.35G:  ' + str(Load[9])+' kN/m')
        f.write('\n1.2G + 1.48Wsls:  '+ str(Load[10])+' kN/m')
        f.write('\nMaximum UDL down:' + str(Load[6])+' kN/m')
        f.write('\nMaximum UDL up:' + str(Load[7]) + ' kN/m')
        f.write('\nMaximum point load down:  '+ str(Load[14]))
        f.write('\nMaximum point load up:  '+ str(Load[15]))
        f.write('\nMaximum positive Moment due to UDL:  '+ str(Load[0])+ ' kNm')
        f.write('\nMaximum negative Moment due to UDL:  ' + str(Load[2])+ ' kNm')
        f.write('\nMaximum shear due to UDL:  '+ str(Load[1])+ ' kN')
        f.write('\nMaximum positive moment due to point load:'  + str(Load[16]))
        f.write('\nMaximum negative moment due to point load:' + str(Load[17]))
        f.write('\nMaximum Shear due to point load:' + str(Load[18]))
        f.write('\n\nRequired Ix for G:  ' + str(Load[3])+ ' mm4')
        f.write('\n\nRequired Ix for phiQ:  ' + str(Load[4])+ ' mm4')
        f.write('\n\nRequired Ix for Wsls:  ' + str(Load[5])+ ' mm4')

Name = 'xx'
Length = 7.4
Glimit = min(Length*1000/360,100)
Qlimit = min(Length*1000/250,100)
Wlimit = min(Length*1000/150,100)
G = 0.5
Q = 2
Wsls = 0
Trib_width = 2.5
Gpoint = 0
Qpoint = 0
Wpointsls = 0
a = Length/2

Write_Loads(Name,G, Q, Wsls, Trib_width, Length, Glimit, Qlimit, Wlimit,Gpoint,Qpoint,Wpointsls, a)



class Load_Case:
    def __init__(self,Name,Length,sample_rate,E,I,values):
        self.Name = Name
        self.Length = Length # In metres
        self.Moment = [0] * int(Length * sample_rate + 1)
        self.Shear = [0] * int(Length * sample_rate + 1)
        self.Deflection = [0] * int(Length * sample_rate + 1)
        try:
            if values['CLength'] > 0:
                self.Deflection = [0]*int(values['CLength']*sample_rate+ Length * sample_rate + 1)
                self.Moment = [0] * int(values['CLength'] * sample_rate + Length * sample_rate + 1)
                self.Shear = [0] * int(values['CLength'] * sample_rate + Length * sample_rate + 1)
        except:
            pass
        self.R1 = 0
        self.R2 = 0
        self.sample_rate = sample_rate
        self.Marea = [0] * int(Length * sample_rate + 1)
        self.E = E
        self.I = I
        self.values =values
    def __str__(self):
        return f"{self.R1} {self.R2} {self.sample_rate}"
    def UDL(self,a,b,w):
        Moment = [0] * int(self.Length * self.sample_rate + 1)
        Shear = [0] * int(self.Length * self.sample_rate + 1)
        Deflection = [0] * int(self.Length * self.sample_rate + 1)
        Marea = [0] * int(self.Length * self.sample_rate + 1)
        c = self.Length - a - b
        R1 = w*b/2/self.Length*(2*c + b)
        R2 = w*b/2/self.Length*(2*a + b)
        for i in range(int(self.Length*self.sample_rate + 1)):
            if i <= a*self.sample_rate:
                Shear[i] = R1
                Moment[i] = R1*i/self.sample_rate
            elif i >a*self.sample_rate and i < (a + b)*self.sample_rate:
                Shear[i] = R1 - w*(i/self.sample_rate - a)
                Moment[i] = R1*i/self.sample_rate - w/2*(i/self.sample_rate - a)**2
            elif i >= (a + b)*self.sample_rate:
                Shear[i] = -R2
                Moment[i] = R2*(self.Length - i/self.sample_rate)
                #print(i)
        for i in range(int(self.Length*self.sample_rate)):
            Marea[i] = np.trapz(np.array(Moment[:i]),dx=1/self.sample_rate)
        centroid = min(range(len(Marea)),key=lambda x:abs(Marea[x]-max(Marea)/2))
        x1 = lambda x: R1*x**2/2
        y1 = integrate.quad(x1,0,a)
        x1a = lambda x: R1*x
        y1a = integrate.quad(x1a,0,a)
        x2 = lambda x: (R1*x - w/2*(x-a)**2)*x
        x2a = lambda x: R1*x - w/2*(x-a)**2
        y2 = integrate.quad(x2,a,a+b)
        y2a = integrate.quad(x2a, a, a + b)
        x3 = lambda x: R2*(self.Length - x)*x
        x3a = lambda x: R2*(self.Length - x)
        y3 = integrate.quad(x3,a+b,self.Length)
        y3a = integrate.quad(x3a, a + b, self.Length)
        try:
            End_deflection = (self.Length - (y1[0] + y2[0] + y3[0])/(y1a[0]+y2a[0] + y3a[0]))*(y1a[0]+y2a[0] + y3a[0])
            End_deflection1 = Marea[int(self.Length * self.sample_rate)]*(self.Length - centroid/self.sample_rate)
            print(End_deflection1,End_deflection)
            for i in range(int(self.Length*self.sample_rate)):
                if i == 0:
                    Deflection[i] = 0
                elif i < a*self.sample_rate:
                    y1 = integrate.quad(x1,0,i/self.sample_rate)
                    y1a = integrate.quad(x1a,0,i/self.sample_rate)

                    Deflection[i] = (i/(self.Length*self.sample_rate)*End_deflection - (i/self.sample_rate - (y1[0])/(y1a[0]))*(y1a[0]))/(self.E*self.I)*10**6

                elif i <a*self.sample_rate + b*self.sample_rate and i >= a*self.sample_rate:
                    y1 = integrate.quad(x1, 0, a)
                    y1a = integrate.quad(x1a, 0, a)
                    y2 = integrate.quad(x2, a,i/self.sample_rate)
                    y2a =  integrate.quad(x2a, a,i/self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (i/self.sample_rate  - (y1[0] + y2[0])/(y1a[0]+y2a[0]))*(y1a[0]+y2a[0]))/(self.E*self.I)*10**6
                elif i >= a*self.sample_rate + b*self.sample_rate:
                    y1 = integrate.quad(x1, 0, a)
                    y1a = integrate.quad(x1a, 0, a)
                    y2 = integrate.quad(x2, a, a + b)
                    y2a = integrate.quad(x2a, a, a+b)
                    y3 = integrate.quad(x3, a + b,i/self.sample_rate)
                    y3a = integrate.quad(x3a, a + b, i/self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (i/self.sample_rate  - (y1[0] + y2[0] + y3[0])/(y1a[0]+y2a[0] + y3a[0]))*(y1a[0]+y2a[0] + y3a[0]))/(self.E*self.I)*10**6
        except:
            pass
        for i in range(int(self.Length*self.sample_rate + 1)):
            self.Moment[i] += -Moment[i]
            self.Shear[i] += Shear[i]
            self.Deflection[i] += -Deflection[i]
            self.Marea[i] += Marea[i]
        try:

            for i in range(int(self.Length*self.sample_rate + 1),int(self.Length * self.sample_rate + self.values['CLength']*self.sample_rate +1)):
                self.Deflection[i] += -(i / (self.Length * self.sample_rate) * End_deflection - (i/self.sample_rate  - (y1[0] + y2[0] + y3[0])/(y1a[0]+y2a[0] + y3a[0]))*(y1a[0]+y2a[0] + y3a[0]))/(self.E*self.I)*10**6
        except:
            pass
        self.R1 += R1
        self.R2 += R2

    def Cantilever_UDL(self,a,b,w):
        Moment = [0] * int(self.Length * self.sample_rate + self.values['CLength']*self.sample_rate +1)
        Shear = [0] * int(self.Length * self.sample_rate + self.values['CLength']*self.sample_rate +1)
        Deflection = [0] * int(self.Length * self.sample_rate + self.values['CLength']*self.sample_rate +1)
        R1 = w*b*(a+b/2)/self.Length
        R2 = R1 + w*b
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
            for i in range(int(self.Length * self.sample_rate + self.values['CLength']*self.sample_rate +1)):
                if i == 0:
                    Moment[i] = R1 * (i / self.sample_rate)
                    Shear[i] = R1
                elif i <= self.Length*self.sample_rate:
                    Moment[i] = R1*(i/self.sample_rate)
                    Shear[i] = R1
                    y1 = integrate.quad(x1, 0, i/self.sample_rate)
                    y1a = integrate.quad(x1a, 0, i/self.sample_rate)
                    Deflection[i] = (i/(self.Length*self.sample_rate)*End_deflection - (i/self.sample_rate - (y1[0])/(y1a[0]))*(y1a[0]))/(self.E*self.I)*10**6
                elif i > self.Length*self.sample_rate and i <= self.Length*self.sample_rate + a*self.sample_rate:
                    Moment[i] = R1*(i/self.sample_rate) - R2*(i/self.sample_rate - self.Length)
                    Shear[i] = R1 -R2
                    y1 = integrate.quad(x1, 0, self.Length)
                    y1a = integrate.quad(x1a, 0, self.Length)
                    y2 = integrate.quad(x2, self.Length, i / self.sample_rate)
                    y2a = integrate.quad(x2a, self.Length, i / self.sample_rate)
                    Deflection[i] = (i / (self.Length * self.sample_rate) * End_deflection - (
                                i / self.sample_rate - (y1[0] + y2[0]) / (y1a[0] + y2a[0])) * (y1a[0] + y2a[0])) / (
                                                self.E * self.I) * 10 ** 6
                elif i > self.Length*self.sample_rate + a*self.sample_rate and i <= self.Length*self.sample_rate + a*self.sample_rate + b*self.sample_rate:
                    Moment[i] = R1*(i/self.sample_rate) - R2*(i/self.sample_rate - self.Length) + w*(i/self.sample_rate - self.Length - a)**2/2
                    Shear[i] = R1 - R2 + w*(i/self.sample_rate - self.Length - a)
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
        for i in range(int(self.Length*self.sample_rate + self.values['CLength']*self.sample_rate + 1)):
            self.Moment[i] += Moment[i]
            self.Shear[i] += -Shear[i]
            self.Deflection[i] += Deflection[i]
        self.R1 += -R1
        self.R2 += R2

    def Point_load(self,a,P):
        b = self.Length - a
        Moment = [0] * int(self.Length * self.sample_rate + 1)
        Shear = [0] * int(self.Length * self.sample_rate + 1)
        Deflection = [0] * int(self.Length * self.sample_rate + 1)
        Marea = [0] * int(self.Length * self.sample_rate + 1)
        R1 = P*b/self.Length
        R2 = P*a/self.Length
        for i in range(int(self.Length*self.sample_rate + 1)):
            if a == 0 or b == 0:
                Moment[i] = P * b * (i / self.sample_rate) / self.Length
                Shear[i] = R1
            elif b == 0:
                Moment[i] = P * a * (self.Length - i / self.sample_rate) / self.Length
                Shear[i] = -R2
            else:
                if i < a*self.sample_rate:
                    Moment[i] = P*b*(i/self.sample_rate)/self.Length
                    Shear[i] = R1
                    Deflection[i] = P*b*(i/self.sample_rate)/(6*self.Length*self.E*self.I)*(self.Length**2 - b**2 - (i/self.sample_rate)**2) *10**6
                else:
                    Moment[i] = P*a*(self.Length - i/self.sample_rate)/self.Length
                    Shear[i] = -R2
                    Deflection[i] = P * a * (self.Length - i / self.sample_rate) / (6 * self.Length * self.E * self.I) * (
                                self.Length*2*(i/self.sample_rate) - a ** 2 - (i / self.sample_rate) ** 2) *10**6
        for i in range(int(self.Length*self.sample_rate + 1)):
            self.Moment[i] += -Moment[i]
            self.Shear[i] += Shear[i]
            self.Deflection[i] += -Deflection[i]
            self.Marea[i] += Marea[i]
        self.R1 += R1
        self.R2 += R2
    def Cantilever_Point_load(self,a,P):
        Moment = [0] * int(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        Shear = [0] * int(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        Deflection = [0] * int(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)
        R1 = P*a/self.Length
        R2 = P/self.Length*(self.Length + a)
        try:
            for i in range(int(self.Length * self.sample_rate + self.values['CLength'] * self.sample_rate + 1)):
                if i <= self.Length*self.sample_rate:
                    Moment[i] = P*a*(i/self.sample_rate)/self.Length
                    Shear[i] = R1
                    Deflection[i] = (P*a*i/self.sample_rate)/(6*self.Length*self.E*self.I)*(self.Length**2 - (i/self.sample_rate)**2)*10**6
                elif i > self.Length*self.sample_rate and i <= self.Length*self.sample_rate + a *self.sample_rate:
                    Moment[i] = P*(a-(i/self.sample_rate - self.Length))
                    Shear[i] = -R2 + R1
                    Deflection[i] = -(P*(i/self.sample_rate - self.Length))/(6*self.E*self.I)*(2*a*self.Length + 3*a*(i/self.sample_rate - self.Length) - (i/self.sample_rate - self.Length)**2)*10**6
                    Deflect = Deflection[i]
                elif i > self.Length*self.sample_rate + a *self.sample_rate:
                    Deflection[i] = (i/self.sample_rate - self.Length)/a * Deflect
        except:
            pass
        for i in range(int(self.Length*self.sample_rate + self.values['CLength']*self.sample_rate + 1)):
            self.Moment[i] += Moment[i]
            self.Shear[i] += -Shear[i]
            self.Deflection[i] += Deflection[i]
        self.R1 += -R1
        self.R2 += R2

class ORMType:
    def __init__(self, key,value):
        self.key = value

    def __repr__(self):
        return "ok"

class Beam:
    def __init__(self,Name,Length,sample_rate,E,Ix,Iy,Loadings,Cpe_V,Cpi_V,Cpe_H,Cpi_H,P,values,Cantilever):
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
    def __str__(self):
        return self.Name

    def LoadCase(self,Loadings,Cantilever):

        for key in Loadings:
            if key == 'WOoP':
                setattr(self,key,Load_Case(self.Name, self.Length, self.sample_rate, self.E, self.Iy,self.values))
            else:
                setattr(self, key, Load_Case(self.Name, self.Length, self.sample_rate, self.E, self.Ix,self.values))
            #self.key = Load_Case(self.name, self.Length, self.sample_rate, self.E, self.I)
            #try:
            for j in range(len(Loadings[key][0])):
                    getattr(self,key).UDL(Loadings[key][0][j][0],Loadings[key][0][j][1],Loadings[key][0][j][2])
                        #self.key.UDL(Loadings[key][0][j][0],Loadings[key][0][j][1],Loadings[key][0][j][2])
            #except:
                #pass
            #try:
            for j in range(len(Loadings[key][1])):

                    getattr(self,key).Point_load(Loadings[key][1][j][0],Loadings[key][1][j][1])
            #except:
               # pass
            if self.values['Cantilever'] == True:

                for j in range(len(Cantilever[key][0])):
                        getattr(self,key).Cantilever_UDL(Cantilever[key][0][j][0],Cantilever[key][0][j][1],Cantilever[key][0][j][2])
                            #self.key.UDL(Loadings[key][0][j][0],Loadings[key][0][j][1],Loadings[key][0][j][2])
                #except:
                    #pass
                #try:
                for j in range(len(Cantilever[key][1])):
                        getattr(self,key).Cantilever_Point_load(Cantilever[key][1][j][0],Cantilever[key][1][j][1])
Loading_cases = ['G1','G2','Q','W']
def Checks(Beam):
    print('////\n', Beam.Name)
    #for i in range(len(Beam.sample_rate*Beam.Length + 1)):
    try:


        G_12_M = [1.2 * i for i in getattr(Beam, 'G2').Moment]
        Q_15_M = [1.5 * i for i in getattr(Beam, 'Q').Moment]
        ULS_M = [x + y for x, y in zip(G_12_M, Q_15_M)]  # 1.2G + 1.5Q

        G_12_S = [1.2 * i for i in getattr(Beam, 'G2').Shear]
        Q_15_S = [1.5 * i for i in getattr(Beam, 'Q').Shear]
        ULS_S = [x + y for x, y in zip(G_12_S, Q_15_S)]  # 1.2G + 1.5Q
        ULS_Wup_M = [0.9*x - y for x, y in zip(getattr(Beam, 'G1').Moment, getattr(Beam, 'Wup').Moment)]

        ULS_Wup_S = [0.9*x - y for x, y in zip(getattr(Beam, 'G1').Shear, getattr(Beam, 'Wup').Shear)]

        ULS_Wdown_M = [x + y for x, y in zip(G_12_M, getattr(Beam,'Wdown').Moment)]
        ULS_Wdown_S = [x + y for x, y in zip(G_12_S, getattr(Beam, 'Wdown').Shear)]
        G_135_M = [1.35 * i for i in getattr(Beam, 'G2').Moment]
        G_135_S = [1.35 * i for i in getattr(Beam, 'G2').Shear]
        print('In-Plane')
        print('The Maximum Moment is: ',
                  round(max(max(ULS_M, key=abs), max(G_135_M, key=abs), max(ULS_Wup_M, key=abs),max(ULS_Wdown_M,key=abs), key=abs), 2), ' KNm')
        print('The Maximum Shear is: ',
                  round(max(max(ULS_S, key=abs), max(G_135_S, key=abs), max(ULS_Wup_S, key=abs),max(ULS_Wdown_S,key=abs), key=abs), 2), ' KN')
        print('Maximum G1 R1 = ',round(getattr(Beam,'G1').R1,2))
        print('Maximum G1 R2 = ', round(getattr(Beam, 'G1').R2, 2))
        print('Maximum G2 R1 = ', round(getattr(Beam, 'G2').R1, 2))
        print('Maximum G2 R2 = ', round(getattr(Beam, 'G2').R2, 2))
        print('Maximum Q R1 = ', round(getattr(Beam, 'Q').R1, 2))
        print('Maximum Q R2 = ', round(getattr(Beam, 'Q').R2, 2))
        print('Maximum R1 = ', round(max(getattr(Beam,'G1').R1,getattr(Beam,'G2').R1,getattr(Beam,'Q').R1,getattr(Beam,'Wup').R1,getattr(Beam,'Wdown').R1,getattr(Beam,'WOoP').R1,key=abs),2))
        print('Maximum R2 = ', round(
            max(getattr(Beam, 'G1').R2, getattr(Beam, 'G2').R2, getattr(Beam, 'Q').R2, getattr(Beam, 'Wup').R2,
                getattr(Beam, 'Wdown').R2), 2))
        print('The maximum Deflection for G is: ', round(max(max(getattr(Beam, 'G2').Deflection, key=abs),max(getattr(Beam, 'G1').Deflection, key=abs),key=abs), 2),
                  'mm L/300 = ', Beam.Length / 0.3)
        print('The maximum Deflection for Q is: ', round(max(getattr(Beam, 'Q').Deflection, key=abs), 2),
                  'mm L/250 = ', Beam.Length / 0.250)
        print('The maximum Deflection for Wup is: ', round(0.7 * max(getattr(Beam, 'Wup').Deflection, key=abs), 2),
                  'mm L/150 = ', Beam.Length / 0.150)
        print('The maximum Deflection for Wdown is: ',
                  round(0.7 * max(getattr(Beam, 'Wdown').Deflection, key=abs), 2), 'mm L/150 = ', Beam.Length / 0.150)

        getattr(Beam, 'WOoP')
        print('\nOut-of-Plane')
        print('The maximum Moment Out of plane is: ', round(max(getattr(Beam, 'WOoP').Moment, key=abs), 2), 'KNm')
        print('The maximum Shear Out of plane is: ', round(max(getattr(Beam, 'WOoP').Shear, key=abs), 2), 'KN')
        print('The maximum Deflection Out of plane is: ',
              round(0.7 * max(getattr(Beam, 'WOoP').Deflection, key=abs), 2), 'mm L/150 =', Beam.Length / 0.150)

    except:
                print('ERROR')
Name = 'L3'






def info(Name):
    try:
        with open(Name+'.pkl', 'rb') as inp:
            foo = pickle.load(inp)
        Checks(foo)
    except (OSError, IOError) as e:
        print('No info')

def open_file(Name,Project):
    try:
        with open(Project + '\\'+ Name+'.pkl', 'rb') as inp:
            foo = pickle.load(inp)
            return foo
    except (OSError, IOError) as e:
        print('No info')


def Method(Name,Length,sample_rate,E,I,Iy,Loadings,Cpe_V,Cpi_V,Cpe_H,Cpi_H,P,values,Cantilever):
    Name1 =  Beam(Name,Length,sample_rate,E,I,Iy,Loadings,Cpe_V,Cpi_V,Cpe_H,Cpi_H,P,values,Cantilever)
    Name1.LoadCase(Loadings,Cantilever)
    Checks(Name1)
    with open(values['Project']+'\\'+Name+'.pkl','wb') as outp:
        pickle.dump(Name1,outp,pickle.HIGHEST_PROTOCOL)

Name = 'FB3'
Length = 2.6
Trib = 5/2
d = 240
b = 90
E = 12*10**9
Cpe = 0.9
Cpi = 0.3
p = 1.28
#L1 = open_file('L1',values['Project'])
#L2 = open_file('L2',values['Project'])
G1 = [[[0,Length,0.2*Trib]]]
G2 = [[[0,Length,0.5*Trib]]]
Q = [[[0,Length,1.5*Trib]]]
Wup = [[[0,Length,(Cpe+Cpi)*p*0]]]
Wdown = [[[0,Length,0*Trib]]]
WOoP = [[[0,Length,(Cpe+Cpi)*p*0]]]
Loadings = {'G1':G1,'G2':G2,'Q':Q,'Wup':Wup,'Wdown':Wdown,'WOoP':WOoP}
#Method(Name,Length,100,E,b*d**3/12*10**-12,d*b**3/12*10**-12,Loadings)

def variable_write(values,Name):
    #if values['Member Type'] == 'Beam':
        with open(Name +'.txt','w',newline='') as csv_file:
            data = [[str(i),values[i]] for i in values]
            print(data)
            writer = csv.writer(csv_file)
            writer.writerows(data)


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
    except KeyError:
        output = 0
    return output
def Cantilever_Layout(variable):
    Layout = Layout = [[sg.Column([
        [sg.Text('Cantilever Length')],
        ]),sg.Column([
            [sg.Input(size=(10,1),default_text=KeyCheck(variable,'CLength'),enable_events=True,key='CLength')],
        ])],
    [sg.Text('Number of UDL\'s'),sg.Input(default_text=int(float(KeyCheck(variable,'CUDLS'))),key='CUDLS',enable_events=True,size=(5,1))]
    ]
    try:
        CUDLS = int(float(KeyCheck(variable,'CUDLS')))
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
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1M'+str(i)),enable_events=True,key='G1M'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2M'+str(i)),enable_events=True,key='G2M'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QM'+str(i)),enable_events=True,key='QM'+str(i))],
                [sg.Text()],
                [sg.Text()],
                [sg.Text()]
                #[sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupM'+str(i)),enable_events=True,key='WupM'+str(i))],
                #[sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownM'+str(i)),enable_events=True,key='WdownM'+str(i))],
                #[sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPM'+str(i)),enable_events=True,key='WOoPM'+str(i))]
            ]),
        sg.Column([
            [sg.Text('Tributary Width m')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1T'+str(i)),enable_events=True,key='G1T'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2T'+str(i)),enable_events=True,key='G2T'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QT'+str(i)),enable_events=True,key='QT'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupT'+str(i)),enable_events=True,key='WupT'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownT'+str(i)),enable_events=True,key='WdownT'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPT'+str(i)),enable_events=True,key='WOoPT'+str(i))]
        ]),
        sg.Column([
            [sg.Text('Start m')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1a'+str(i)),enable_events=True,key='G1a'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2a'+str(i)),enable_events=True,key='G2a'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Qa'+str(i)),enable_events=True,key='Qa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wupa'+str(i)),enable_events=True,key='Wupa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wdowna'+str(i)),enable_events=True,key='Wdowna'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPa'+str(i)),enable_events=True,key='WOoPa'+str(i))],
        ]),
        sg.Column([
            [sg.Text('Extent m')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1b'+str(i)),enable_events=True,key='G1b'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2b'+str(i)),enable_events=True,key='G2b'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Qb'+str(i)),enable_events=True,key='Qb'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wupb'+str(i)),enable_events=True,key='Wupb'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wdownb'+str(i)),enable_events=True,key='Wdownb'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPb'+str(i)),enable_events=True,key='WOoPb'+str(i))],
        ])]]
    Layout += [[sg.Text('New point Loads'),sg.Input(key='CPoint_loads',default_text=int(float(KeyCheck(variable,'CPoint_loads'))),enable_events=True,size=(5,1))]]
    try:
        CPoint_loads = int(float(KeyCheck(variable,'CPoint_loads')))
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
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1P'+str(i)),enable_events=True,key='G1P'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2P'+str(i)),enable_events=True,key='G2P'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QP'+str(i)),enable_events=True,key='QP'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupP'+str(i)),enable_events=True,key='WupP'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownP'+str(i)),enable_events=True,key='WdownP'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPP'+str(i)),enable_events=True,key='WOoPP'+str(i))]
            ]),
        sg.Column([
            [sg.Text('Location')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1Pa'+str(i)),enable_events=True,key='G1Pa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2Pa'+str(i)),enable_events=True,key='G2Pa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QPa'+str(i)),enable_events=True,key='QPa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupPa'+str(i)),enable_events=True,key='WupPa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownPa'+str(i)),enable_events=True,key='WdownPa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPPa'+str(i)),enable_events=True,key='WOoPPa'+str(i))]
        ])]]

    Layout += [[sg.Text('Existing point loads'),
                sg.Input(key='CEx_Point_loads', default_text=int(float(KeyCheck(variable, 'CEx_Point_loads'))), enable_events=True,
                         size=(5, 1))]]
    try:
        directory = KeyCheck(variable,'Project')
    except:
        directory = os.fsencode(os.getcwd())
    existing = []
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".pkl"):
            filename = os.path.basename(filename)
            filename = os.path.splitext(filename)[0]
            if filename != KeyCheck(variable,'Name'):
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
            Layout += [[sg.Combo(existing,size=(10,1),key='Ex'+str(i),default_value= KeyCheck(variable,'Ex'+str(i)),enable_events=True),
                   sg.Combo(['Left','Right'],size=(5,1),key='L/R' + str(i),default_value=KeyCheck(variable,'L/R'+str(i)),enable_events=True),
                   sg.Input(default_text=KeyCheck(variable,'Exa'+str(i)),key = 'Exa'+str(i),size=(5,1),enable_events=True)]]
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

def draw_graph(Beam):

    #fig = plt.figure(figsize=(12,30),dpi=80)
    fig = matplotlib.figure.Figure(figsize=(12,30),dpi=80)
    try:
        Length = np.arange(0,Beam.Length + Beam.values['CLength']+ 1/Beam.sample_rate,1/Beam.sample_rate)
        Zeros = [0]*(int(Beam.Length + Beam.values['CLength'])*Beam.sample_rate+1)
    except:
        Length =  np.arange(0,Beam.Length + 1/Beam.sample_rate,1/Beam.sample_rate)
        Zeros = [0] * (int(Beam.Length) * Beam.sample_rate+1)

    G_12_M = [1.2 * i for i in getattr(Beam, 'G2').Moment]
    Q_15_M = [1.5 * i for i in getattr(Beam, 'Q').Moment]
    ULS_M = [x + y for x, y in zip(G_12_M, Q_15_M)]  # 1.2G + 1.5Q

    G_12_S = [1.2 * i for i in getattr(Beam, 'G2').Shear]
    Q_15_S = [1.5 * i for i in getattr(Beam, 'Q').Shear]
    ULS_S = [x + y for x, y in zip(G_12_S, Q_15_S)]  # 1.2G + 1.5Q
    ULS_Wup_M = [0.9 * x - y for x, y in zip(getattr(Beam, 'G1').Moment, getattr(Beam, 'Wup').Moment)]

    ULS_Wup_S = [0.9 * x - y for x, y in zip(getattr(Beam, 'G1').Shear, getattr(Beam, 'Wup').Shear)]

    ULS_Wdown_M = [x + y for x, y in zip(G_12_M, getattr(Beam, 'Wdown').Moment)]
    ULS_Wdown_S = [x + y for x, y in zip(G_12_S, getattr(Beam, 'Wdown').Shear)]
    G_135_M = [1.35 * i for i in getattr(Beam, 'G2').Moment]
    G_135_S = [1.35 * i for i in getattr(Beam, 'G2').Shear]

    ax1 = fig.add_subplot(10,2,1)
    ax1.set_ylabel('Moment [KNm]')
    ax1.set_xlabel('Length [m]')
    ax1.set_title('1.2G + 1.5Q MOMENT')
    ax1.plot(Length,ULS_M)
    ax1.text(0,-abs(1/10*max(ULS_M,key=abs)),str(round(1.2*Beam.G2.R1 + 1.5*Beam.Q.R1,1))+'KN',verticalalignment='top',horizontalalignment='center')
    ax1.text(Beam.Length, -abs(1 / 10 * max(ULS_M, key=abs)), str(round(1.2 * Beam.G2.R2 + 1.5 * Beam.Q.R2, 1)) + 'KN',
             verticalalignment='top', horizontalalignment='center')
    ax1.plot([0,Beam.Length],[0,0],'g',marker=6,markersize=15)
    try:
        ax1.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0],'g')
    except:
       pass

    ax2 = fig.add_subplot(10,2,2)
    ax2.set_ylabel('Shear [KNm]')
    ax2.set_xlabel('Length [m]')
    ax2.set_title('1.2G + 1.5Q Shear')
    ax2.plot(Length, ULS_S,'b')
    ax2.text(0, -abs(1 / 10 * max(ULS_S, key=abs)), str(round(1.2 * Beam.G2.R1 + 1.5 * Beam.Q.R1, 1)) + 'KN',
             verticalalignment='top', horizontalalignment='center')
    ax2.text(Beam.Length, -abs(1 / 10 * max(ULS_S, key=abs)), str(round(1.2 * Beam.G2.R2 + 1.5 * Beam.Q.R2, 1)) + 'KN',
             verticalalignment='top', horizontalalignment='center')
    ax2.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax2.plot([0,0],[0,ULS_S[0]],'b')
    ax2.plot([Length[-1],Length[-1]],[0,ULS_S[-1]],'b')
    try:
        ax2.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass

    ax1a = fig.add_subplot(10, 2, 3)
    ax1a.set_ylabel('Moment [KNm]')
    ax1a.set_xlabel('Length [m]')
    ax1a.set_title('1.2G + Wdown MOMENT')
    ax1a.plot(Length, ULS_Wdown_M)
    ax1a.text(0, -abs(1 / 10 * max(ULS_Wdown_M, key=abs)), str(round(1.2 * Beam.G2.R1 +  Beam.Wdown.R1, 1)) + 'KN',
             verticalalignment='top', horizontalalignment='center')
    ax1a.text(Beam.Length, -abs(1 / 10 * max(ULS_Wdown_M, key=abs)), str(round(1.2 * Beam.G2.R2 +  Beam.Wdown.R2, 1)) + 'KN',
             verticalalignment='top', horizontalalignment='center')
    ax1a.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    try:
        ax1a.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax2a = fig.add_subplot(10, 2, 4)
    ax2a.set_ylabel('Shear [KNm]')
    ax2a.set_xlabel('Length [m]')
    ax2a.set_title('1.2G + Wdown Shear')
    ax2a.plot(Length, ULS_Wdown_S,'b')
    ax2a.text(0, -abs(1 / 10 * max(ULS_Wdown_S, key=abs)), str(round(1.2 * Beam.G2.R1 + Beam.Wdown.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax2a.text(Beam.Length, -abs(1 / 10 * max(ULS_Wdown_S, key=abs)), str(round(1.2 * Beam.G2.R2 + Beam.Wdown.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax2a.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax2a.plot([0, 0], [0, ULS_Wdown_S[0]], 'b')
    ax2a.plot([Length[-1], Length[-1]], [0, ULS_Wdown_S[-1]], 'b')
    try:
        ax2a.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax1b = fig.add_subplot(10, 2, 5)
    ax1b.set_ylabel('Moment [KNm]')
    ax1b.set_xlabel('Length [m]')
    ax1b.set_title('1.35G MOMENT')
    ax1b.plot(Length, G_135_M)
    ax1b.text(0, -abs(1 / 10 * max(G_135_M, key=abs)), str(round(1.35 * Beam.G2.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax1b.text(Beam.Length, -abs(1 / 10 * max(G_135_M, key=abs)),
              str(round(1.35 * Beam.G2.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax1b.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    try:
        ax1b.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax2b = fig.add_subplot(10, 2, 6)
    ax2b.set_ylabel('Shear [KNm]')
    ax2b.set_xlabel('Length [m]')
    ax2b.set_title('1.35G Shear')
    ax2b.plot(Length,G_135_S,'b')
    ax2b.text(0, -abs(1 / 10 * max(G_135_S, key=abs)), str(round(1.35 * Beam.G2.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax2b.text(Beam.Length, -abs(1 / 10 * max(G_135_S, key=abs)),
              str(round(1.35 * Beam.G2.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax2b.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax2b.plot([0, 0], [0, G_135_S[0]], 'b')
    ax2b.plot([Length[-1], Length[-1]], [0, G_135_S[-1]], 'b')
    try:
        ax2b.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax3a = fig.add_subplot(10, 2, 7)
    ax3a.set_ylabel('Deflection [mm]')
    ax3a.set_xlabel('Length [m]')
    ax3a.set_title('G1 Deflection')
    ax3a.plot(Length, getattr(Beam.G1, 'Deflection'))
    ax3a.text(0, -abs(1 / 10 * max(getattr(Beam.G1, 'Deflection'), key=abs)), str(round(Beam.G1.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax3a.text(Beam.Length, -abs(1 / 10 * max(getattr(Beam.G1, 'Deflection'), key=abs)),
              str(round(Beam.G1.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax3a.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax3a.plot([0, Beam.Length], [Beam.Length/0.36, Beam.Length/0.36], linestyle=(0, (5, 10)),color = 'red')
    ax3a.plot([0, Beam.Length], [-Beam.Length / 0.36, -Beam.Length / 0.36], linestyle=(0, (5, 10)),color = 'red')
    ax3a.plot([Beam.Length, max(Length)], [(max(Length)-Beam.Length) / 0.18, (max(Length)-Beam.Length) / 0.18], linestyle=(0, (5, 10)),color = 'red')
    ax3a.plot([Beam.Length, max(Length)], [-(max(Length)-Beam.Length)/ 0.18, -(max(Length)-Beam.Length) / 0.18], linestyle=(0, (5, 10)),color = 'red')
    try:
        ax3a.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax3 = fig.add_subplot(10, 2, 8)
    ax3.set_ylabel('Deflection [mm]')
    ax3.set_xlabel('Length [m]')
    ax3.set_title('G2 Deflection')
    ax3.plot(Length, getattr(Beam.G2,'Deflection'))
    ax3.text(0, -abs(1 / 10 * max(getattr(Beam.G2, 'Deflection'), key=abs)), str(round(Beam.G2.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax3.text(Beam.Length, -abs(1 / 10 * max(getattr(Beam.G2, 'Deflection'), key=abs)),
              str(round(Beam.G2.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax3.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax3.plot([0, Beam.Length], [Beam.Length / 0.36, Beam.Length / 0.36],linestyle=(0, (5, 10)),color = 'red')
    ax3.plot([0, Beam.Length], [-Beam.Length / 0.36, -Beam.Length / 0.36], linestyle=(0, (5, 10)),color = 'red')
    ax3.plot([Beam.Length, max(Length)], [(max(Length)-Beam.Length)  / 0.18, (max(Length)-Beam.Length)  / 0.18], linestyle=(0, (5, 10)),color = 'red')
    ax3.plot([Beam.Length, max(Length)], [-(max(Length)-Beam.Length)  / 0.18, -(max(Length)-Beam.Length)  / 0.18], linestyle=(0, (5, 10)),color = 'red')
    try:
        ax3.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax4 = fig.add_subplot(10, 2, 9)
    ax4.set_ylabel('Deflection [mm]')
    ax4.set_xlabel('Length [m]')
    ax4.set_title('Q Deflection')
    ax4.plot(Length, getattr(Beam.Q, 'Deflection'))
    ax4.text(0, -abs(1 / 10 * max(getattr(Beam.Q, 'Deflection'), key=abs)), str(round(Beam.Q.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax4.text(Beam.Length, -abs(1 / 10 * max(getattr(Beam.Q, 'Deflection'), key=abs)),
              str(round(Beam.Q.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax4.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax4.plot([0, Beam.Length], [Beam.Length / 0.25, Beam.Length / 0.25], linestyle=(0, (5, 10)),color = 'red')
    ax4.plot([0, Beam.Length], [-Beam.Length / 0.25, -Beam.Length / 0.25], linestyle=(0, (5, 10)),color = 'red')
    ax4.plot([Beam.Length, max(Length)], [(max(Length)-Beam.Length)  / 0.125, (max(Length)-Beam.Length)  / 0.125], linestyle=(0, (5, 10)),color = 'red')
    ax4.plot([Beam.Length, max(Length)], [-(max(Length)-Beam.Length)  / 0.125, -(max(Length)-Beam.Length)  / 0.125], linestyle=(0, (5, 10)),color = 'red')
    try:
        ax4.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    ax5 = fig.add_subplot(10, 2, 10)
    ax5.set_ylabel('Deflection [mm]')
    ax5.set_xlabel('Length [m]')
    ax5.set_title('Wind Down Deflection')
    ax5.plot(Length, getattr(Beam.Wdown, 'Deflection'))
    ax5.text(0, -abs(1 / 10 * max(getattr(Beam.Wdown, 'Deflection'), key=abs)), str(round(Beam.Wdown.R1, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax5.text(Beam.Length, -abs(1 / 10 * max(getattr(Beam.Wdown, 'Deflection'), key=abs)),
              str(round(Beam.Wdown.R2, 1)) + 'KN',
              verticalalignment='top', horizontalalignment='center')
    ax5.plot([0, Beam.Length], [0, 0], 'g', marker=6, markersize=15)
    ax5.plot([0, Beam.Length], [Beam.Length / 0.15, Beam.Length / 0.15], linestyle=(0, (5, 10)),color = 'red')
    ax5.plot([0, Beam.Length], [-Beam.Length / 0.15, -Beam.Length / 0.15], linestyle=(0, (5, 10)),color = 'red')
    ax5.plot([Beam.Length, max(Length)], [(max(Length)-Beam.Length)  / 0.075, (max(Length)-Beam.Length)  / 0.075], linestyle=(0, (5, 10)),color = 'red')
    ax5.plot([Beam.Length, max(Length)], [-(max(Length)-Beam.Length)  / 0.075, -(max(Length)-Beam.Length)  / 0.075], linestyle=(0, (5, 10)),color = 'red')
    try:
        ax5.plot([Beam.Length, Beam.Length + Beam.values['CLength']], [0, 0], 'g')
    except:
        pass


    fig.tight_layout()
    return fig
def Layouts(variable):
    from openpyxl import load_workbook

    wb = load_workbook(filename='Steel Design Calculator.xlsx')
    # SectionType = wb.sheetnames
    SectionSize = []
    # print(sheet_ranges['A6'].value)
    SectionType = ['Universal_Beam', 'Universal_Column', 'PFC', 'RHS', 'SHS', 'CHS', 'MGP10','MGP12']
    for j in wb['Universal_Beam']:
        SectionSize.append(j[0].value)
    SectionSize = list(filter(lambda item: item is not None, SectionSize))

    for j1 in range(4):
        SectionSize.pop(0)

    Layout = [[sg.Column([
        [sg.Text('Name')],
        [sg.Text('Length')],
        [sg.Text('E')],
        [sg.Text('Choose Section Type:')],
        [sg.Text('Choose Section Size:')],
        #[sg.Text('b')],
        #[sg.Text('d')],
        [sg.Text('Ix')],
        [sg.Text('Iy')]
        ]),sg.Column([
            [sg.Input(size=(10,1),default_text=KeyCheck(variable,'Name'),enable_events=True,key='Name'),sg.Button('Reload Existing',key='Reload')],
            [sg.Input(size=(10,1),default_text=KeyCheck(variable,'Length'),enable_events=True,key='Length'),sg.Checkbox('Cantilever',default=KeyCheck(variable,'Cantilever'),key='Cantilever',enable_events=True)],
            [sg.Input(size=(10,1),default_text=KeyCheck(variable,'E'),enable_events=True,key='E'),sg.Checkbox('Override for E',default=KeyCheck(variable,'E_check'),key='E_check',enable_events=True)],
            [sg.Combo(SectionType, key='SectionType', enable_events=True, default_value=KeyCheck(variable,'SectionType'), size=(30, 1))],  #SectionType[0]
            [sg.Combo(SectionSize, key='SectionSize',enable_events=True, default_value=KeyCheck(variable, 'SectionSize'), size=(30, 1))],  # SectionSize[0]
            #[sg.Input(size=(10,1),default_text=KeyCheck(variable,'b'),enable_events=True,key='b')],
            #[sg.Input(size=(10,1),default_text=KeyCheck(variable,'d'),enable_events=True,key='d')],
            [sg.Input(size=(10,1),default_text=KeyCheck(variable,'Ix'),enable_events=True,key='Ix'),sg.Checkbox('Override for Ix',default=KeyCheck(variable,'Ix_check'),key = 'Ix_check',enable_events=True)],
            [sg.Input(size=(10, 1), default_text=KeyCheck(variable, 'Iy'), enable_events=True, key='Iy'),sg.Checkbox('Override for Iy',default=KeyCheck(variable,'Iy_check'),key = 'Iy_check',enable_events=True)]
        ]),sg.Column([[sg.Text('Current Project')],
        [sg.Text('Base Wind pressure:')],
        [sg.Text('Vertical Cp,e:')],
        [sg.Text('Vertical Cp,i:')],
        [sg.Text('Horizontal Cp,e:')],
        [sg.Text('Horizontal Cp,i:')]
    ],vertical_alignment='t'),
    sg.Column([[sg.Input(key='Project',default_text=KeyCheck(variable,'Project'),enable_events=True),sg.FolderBrowse('Change Project')],
        [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'P'), enable_events=True, key='P')],
        [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Cpe_V'),enable_events=True,key='Cpe_V')],
        [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpi_V'), enable_events=True, key='Cpi_V')],
        [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpe_H'), enable_events=True, key='Cpe_H')],
        [sg.Input(size=(5, 1), default_text=KeyCheck(variable, 'Cpi_H'), enable_events=True, key='Cpi_H')]

    ],vertical_alignment='t')],
    [sg.Text('Number of UDL\'s'),sg.Input(default_text=int(float(KeyCheck(variable,'UDLS'))),key='UDLS',enable_events=True,size=(5,1))]
    ]
    try:
        UDLS = int(float(KeyCheck(variable,'UDLS')))
    except:
        UDLS = 1
    for i in range(UDLS):
        Layout += [[sg.Input(default_text=KeyCheck(variable,'Load_case'+str(i)),enable_events=True,key='Load_case'+str(i))],[sg.Column([
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
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1M'+str(i)),enable_events=True,key='G1M'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2M'+str(i)),enable_events=True,key='G2M'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QM'+str(i)),enable_events=True,key='QM'+str(i))],
                [sg.Text()],
                [sg.Text()],
                [sg.Text()]
                #[sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupM'+str(i)),enable_events=True,key='WupM'+str(i))],
                #[sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownM'+str(i)),enable_events=True,key='WdownM'+str(i))],
                #[sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPM'+str(i)),enable_events=True,key='WOoPM'+str(i))]
            ]),
        sg.Column([
            [sg.Text('Tributary Width m')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1T'+str(i)),enable_events=True,key='G1T'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2T'+str(i)),enable_events=True,key='G2T'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QT'+str(i)),enable_events=True,key='QT'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupT'+str(i)),enable_events=True,key='WupT'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownT'+str(i)),enable_events=True,key='WdownT'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPT'+str(i)),enable_events=True,key='WOoPT'+str(i))]
        ]),
        sg.Column([
            [sg.Text('Start m')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1a'+str(i)),enable_events=True,key='G1a'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2a'+str(i)),enable_events=True,key='G2a'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Qa'+str(i)),enable_events=True,key='Qa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wupa'+str(i)),enable_events=True,key='Wupa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wdowna'+str(i)),enable_events=True,key='Wdowna'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPa'+str(i)),enable_events=True,key='WOoPa'+str(i))],
        ]),
        sg.Column([
            [sg.Text('Extent m')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1b'+str(i)),enable_events=True,key='G1b'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2b'+str(i)),enable_events=True,key='G2b'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Qb'+str(i)),enable_events=True,key='Qb'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wupb'+str(i)),enable_events=True,key='Wupb'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'Wdownb'+str(i)),enable_events=True,key='Wdownb'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPb'+str(i)),enable_events=True,key='WOoPb'+str(i))],
        ])]]
    Layout += [[sg.Text('New point Loads'),sg.Input(key='Point_loads',default_text=int(float(KeyCheck(variable,'Point_loads'))),enable_events=True,size=(5,1))]]
    try:
        Point_loads = int(float(KeyCheck(variable,'Point_loads')))
    except:
        Point_loads = 1
    for i in range(Point_loads):
        Layout += [[sg.Input(default_text=KeyCheck(variable,'Load_caseP'+str(i)),enable_events=True,key='Load_caseP'+str(i))],[sg.Column([
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
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1P'+str(i)),enable_events=True,key='G1P'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2P'+str(i)),enable_events=True,key='G2P'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QP'+str(i)),enable_events=True,key='QP'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupP'+str(i)),enable_events=True,key='WupP'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownP'+str(i)),enable_events=True,key='WdownP'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPP'+str(i)),enable_events=True,key='WOoPP'+str(i))]
            ]),
        sg.Column([
            [sg.Text('Location')],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G1Pa'+str(i)),enable_events=True,key='G1Pa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'G2Pa'+str(i)),enable_events=True,key='G2Pa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'QPa'+str(i)),enable_events=True,key='QPa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WupPa'+str(i)),enable_events=True,key='WupPa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WdownPa'+str(i)),enable_events=True,key='WdownPa'+str(i))],
                [sg.Input(size=(5,1),default_text=KeyCheck(variable,'WOoPPa'+str(i)),enable_events=True,key='WOoPPa'+str(i))]
        ])]]

    Layout += [[sg.Text('Existing point loads'),
                sg.Input(key='Ex_Point_loads', default_text=int(float(KeyCheck(variable, 'Ex_Point_loads'))), enable_events=True,
                         size=(5, 1))]]
    try:
        directory = KeyCheck(variable,'Project')
        existing = []
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".pkl"):
                filename = os.path.basename(filename)
                filename = os.path.splitext(filename)[0]
                if filename != KeyCheck(variable,'Name'):
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
            Layout += [[sg.Combo(existing,size=(10,1),key='Ex'+str(i),default_value= KeyCheck(variable,'Ex'+str(i)),enable_events=True),
                   sg.Combo(['Left','Right'],size=(5,1),key='L/R' + str(i),default_value=KeyCheck(variable,'L/R'+str(i)),enable_events=True),
                   sg.Input(default_text=KeyCheck(variable,'Exa'+str(i)),key = 'Exa'+str(i),size=(5,1),enable_events=True)]]
        except:
            Ex_Point_loads = 0

    Layout += [[sg.Column([
        [sg.Button('Calculate',key = 'Calculate')],
    [sg.Text(key='Output')]
    ])]]
    if KeyCheck(variable,'Cantilever') == True:
        Layout = [[sg.Column([[sg.Column(Layout,vertical_alignment='t'),sg.VSeparator(),sg.Column(Cantilever_Layout(variable),vertical_alignment='t'),sg.VSeparator(),sg.Column([[sg.Canvas(key='ULSGraph',size=(1200,2500))]],vertical_alignment='t'),sg.VSeparator(),sg.Column(steel_calculator(),vertical_alignment='t')]],scrollable=True,expand_x=True,expand_y=True)]]
    else:
        Layout = [[sg.Column([[sg.Column(Layout,vertical_alignment='t'), sg.VSeparator(),sg.Column([[sg.Canvas(key='ULSGraph',size=(1200,2500))]],vertical_alignment='t'),sg.VSeparator(),
                   sg.Column(steel_calculator(),vertical_alignment='t')]],scrollable=True,expand_x=True,expand_y=True)]]
    window = sg.Window('Window',Layout,resizable=True).finalize()
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
        event,values = window.read()
        if event == sg.WINDOW_CLOSED or event == 'Quit':
            break
        for i in values:
            if event == i:
                variable_write(values, 'Loading')
        if event == 'UDLS' or event == 'Point_loads' or event == 'Ex_Point_loads' or event == 'CUDLS' or event == 'CPoint_loads' or event == 'CEx_Point_loads' or event == 'Cantilever':
            window.close()
            variable = variables('Loading')
            Layouts(variable)
        if event == 'Length':
            window.close()
            variable = variables('Loading')
            for i in range(int(variable['UDLS'])):
                variable['G1b' + str(i)] = values['Length']
                variable['G2b' + str(i)] = values['Length']
                variable['Qb' + str(i)] = values['Length']
                variable['Wupb' + str(i)] = values['Length']
                variable['Wdownb' + str(i)] = values['Length']
                variable['WOoPb' + str(i)] = values['Length']
            Layouts(variable)
        if event == 'CLength':
            window.close()
            variable = variables('Loading')
            for i in range(int(variable['CUDLS'])):
                i = str(i)+'C'
                variable['G1b' + str(i)] = values['Length']
                variable['G2b' + str(i)] = values['Length']
                variable['Qb' + str(i)] = values['Length']
                variable['Wupb' + str(i)] = values['Length']
                variable['Wdownb' + str(i)] = values['Length']
                variable['WOoPb' + str(i)] = values['Length']
            Layouts(variable)
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
                G1U += [[values['G1a'+str(i)],values['G1b'+str(i)],values['G1M'+str(i)]*values['G1T'+str(i)]]]
                G2U += [[values['G2a'+str(i)],values['G2b'+str(i)],values['G2M'+str(i)]*values['G2T'+str(i)]]]
                QU += [[values['Qa' + str(i)], values['Qb' + str(i)], values['QM' + str(i)] * values['QT' + str(i)]]]
                WupM = (values['Cpe_V'] + values['Cpi_V'])*values['P']
                WupU += [[values['Wupa' + str(i)], values['Wupb' + str(i)], WupM * values['WupT' + str(i)]]]
                WdownM = (values['Cpe_V'] + values['Cpi_V']) * values['P']
                WdownU += [[values['Wdowna' + str(i)], values['Wdownb' + str(i)], WdownM * values['WdownT' + str(i)]]]
                WOoPM = (values['Cpe_H'] + values['Cpi_H'])*values['P']
                WOoPU += [[values['WOoPa' + str(i)], values['WOoPb' + str(i)], WOoPM * values['WOoPT' + str(i)]]]
            for i in range(Point_loads):
                G1P += [[values['G1Pa' + str(i)], values['G1P' + str(i)]]]
                G2P += [[values['G2Pa' + str(i)], values['G2P' + str(i)]]]
                QP += [[values['QPa' + str(i)], values['QP' + str(i)]]]
                WupP += [[values['WupPa' + str(i)], values['WupP' + str(i)]]]
                WdownP += [[values['WdownPa' + str(i)], values['WdownP' + str(i)]]]
                WOoPP += [[values['WOoPPa' + str(i)], values['WOoPP' + str(i)]]]
            for i in range(Ex_Point_loads):
                Beam = open_file(values['Ex'+str(i)],values['Project'])
                if values['L/R'+str(i)] == 'Left':
                    G1P += [[values['Exa' + str(i)],Beam.G1.R1]]
                    G2P += [[values['Exa' + str(i)],Beam.G2.R1]]
                    QP += [[values['Exa' + str(i)],Beam.Q.R1]]
                    WupP += [[values['Exa' + str(i)],Beam.Wup.R1]]
                    WdownP += [[values['Exa' + str(i)],Beam.Wdown.R1]]
                    WOoPP += [[ values['Exa' + str(i)],Beam.WOoP.R1]]
                else:
                    G1P += [[ values['Exa' + str(i)],Beam.G1.R2]]
                    G2P += [[values['Exa' + str(i)],Beam.G2.R2]]
                    QP += [[values['Exa' + str(i)],Beam.Q.R2]]
                    WupP += [[values['Exa' + str(i)],Beam.Wup.R2]]
                    WdownP += [[values['Exa' + str(i)],Beam.Wdown.R2]]
                    WOoPP += [[values['Exa' + str(i)],Beam.WOoP.R2]]
            Loadings = {'G1': [G1U,G1P], 'G2': [G2U,G2P], 'Q': [QU,QP], 'Wup': [WupU,WupP], 'Wdown': [WdownU,WdownP], 'WOoP': [WOoPU,WOoPP]}
            if values['Ix_check'] == False:
                try:
                    SectionType1 = values['SectionType']
                    SectionSize1 = values['SectionSize']
                    section_properties = st.section_properties(SectionType1, SectionSize1, 0, 0, 0)
                    Ix = section_properties['Ix'] * 10 **-12
                except:
                    MGP10 = {'90x90':[90,90],'140x90':[140,90],'190x90':[190,90],'240x90':[240,90],'300x90':[300,90]}
                    Timber = MGP10[values['SectionSize']]
                    Ix = Timber[0]*Timber[1]**3/12 * 10 **-12
            else:
                Ix = values['Ix']*10**-12
            if values['Iy_check'] == False:
                try:
                    SectionType1 = values['SectionType']
                    SectionSize1 = values['SectionSize']
                    section_properties = st.section_properties(SectionType1, SectionSize1, 0, 0, 0)
                    Iy = section_properties['Iy'] * 10 **-12
                except:
                    MGP10 = {'90x90': [90, 90], '140x90': [140, 90], '190x90': [190, 90], '240x90': [240, 90],
                             '300x90': [300, 90]}
                    Timber = MGP10[values['SectionSize']]
                    Iy = Timber[1] * Timber[0] ** 3 / 12 * 10 ** -12
            else:
                Iy = values['Iy']*10**-12

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
                    Beam = open_file(values['Ex' + str(i)],values['Project'])
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
                E = values['E']
            else:
                if values['SectionType'] == 'MGP10':
                    E = 10*10**9
                elif values['SectionType'] == 'MGP12':
                    E = 12.7 * 10 **9
                else:
                    E = 200*10**9
            Method(values['Name'], values['Length'], 100, E, Ix, Iy, Loadings,
                       values['Cpe_V'], values['Cpi_V'], values['Cpe_H'], values['Cpi_H'], values['P'], values,Cantilever)

            Beam1 = open_file(values['Name'],values['Project'])
            fig_canvas_agg = draw_figure(window['ULSGraph'].TKCanvas,draw_graph(Beam1))
            try:

                G_12_M = [1.2 * i for i in getattr(Beam1, 'G2').Moment]
                Q_15_M = [1.5 * i for i in getattr(Beam1, 'Q').Moment]
                ULS_M = [x + y for x, y in zip(G_12_M, Q_15_M)]  # 1.2G + 1.5Q

                G_12_S = [1.2 * i for i in getattr(Beam1, 'G2').Shear]
                Q_15_S = [1.5 * i for i in getattr(Beam1, 'Q').Shear]
                ULS_S = [x + y for x, y in zip(G_12_S, Q_15_S)]  # 1.2G + 1.5Q
                ULS_Wup_M = [0.9 * x - y for x, y in zip(getattr(Beam1, 'G1').Moment, getattr(Beam1, 'Wup').Moment)]

                ULS_Wup_S = [0.9 * x - y for x, y in zip(getattr(Beam1, 'G1').Shear, getattr(Beam1, 'Wup').Shear)]

                ULS_Wdown_M = [x + y for x, y in zip(G_12_M, getattr(Beam1, 'Wdown').Moment)]
                ULS_Wdown_S = [x + y for x, y in zip(G_12_S, getattr(Beam1, 'Wdown').Shear)]
                G_135_M = [1.35 * i for i in getattr(Beam1, 'G2').Moment]
                G_135_S = [1.35 * i for i in getattr(Beam1, 'G2').Shear]
                String = 'In-Plane\n' + 'The Maximum Moment is: ' + str(round(max(max(ULS_M, key=abs), max(G_135_M, key=abs), max(ULS_Wup_M, key=abs),
                                max(ULS_Wdown_M, key=abs), key=abs), 2)) + ' KNm\n'
                String += 'The Maximum Shear is: ' + str(round(max(max(ULS_S, key=abs), max(G_135_S, key=abs), max(ULS_Wup_S, key=abs),
                                max(ULS_Wdown_S, key=abs), key=abs), 2)) + ' KN\n'


                String += '\nMaximum R1 = ' + str(round(
                    max(getattr(Beam1, 'G1').R1, getattr(Beam1, 'G2').R1, getattr(Beam1, 'Q').R1, getattr(Beam1, 'Wup').R1,
                        getattr(Beam1, 'Wdown').R1, getattr(Beam1, 'WOoP').R1, key=abs), 2))
                String += '\nMaximum R2 = ' + str(round(
                    max(getattr(Beam1, 'G1').R2, getattr(Beam1, 'G2').R2, getattr(Beam1, 'Q').R2, getattr(Beam1, 'Wup').R2,
                        getattr(Beam1, 'Wdown').R2), 2))
                String += '\nThe maximum Deflection for G is: ' + str(round(
                    max(max(getattr(Beam1, 'G2').Deflection, key=abs), max(getattr(Beam1, 'G1').Deflection, key=abs),
                        key=abs), 2)) + 'mm L/300 = ' + str(Beam1.Length / 0.3)
                String += '\nThe maximum Deflection for Q is: ' + str(round(max(getattr(Beam1, 'Q').Deflection, key=abs), 2)) + 'mm L/250 = ' + str(Beam1.Length / 0.250)
                window['Output'].update(String)
                print('The maximum Deflection for Wup is: ',
                      round(0.7 * max(getattr(Beam, 'Wup').Deflection, key=abs), 2),
                      'mm L/150 = ', Beam.Length / 0.150)
                print('The maximum Deflection for Wdown is: ',
                      round(0.7 * max(getattr(Beam, 'Wdown').Deflection, key=abs), 2), 'mm L/150 = ',
                      Beam.Length / 0.150)

                getattr(Beam, 'WOoP')
                print('\nOut-of-Plane')
                print('The maximum Moment Out of plane is: ', round(max(getattr(Beam, 'WOoP').Moment, key=abs), 2),
                      'KNm')
                print('The maximum Shear Out of plane is: ', round(max(getattr(Beam, 'WOoP').Shear, key=abs), 2), 'KN')
                print('The maximum Deflection Out of plane is: ',
                      round(0.7 * max(getattr(Beam, 'WOoP').Deflection, key=abs), 2), 'mm L/150 =', Beam.Length / 0.150)

            except:
                print('ERROR')


            #print(Beam1.Loadings)
        if event == 'Reload':
            try:
                directory = KeyCheck(variable,'Project')
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
                [sg.Combo(existing,key='Reload_Ex',enable_events=True)],
                [sg.Button('Continue',key='Reload_continue')]
            ]
            Project = values['Project']
            window1 = sg.Window('Select Existing Beam',Layout1,resizable=True).finalize()
            window.close()
            while True:
                event, values = window1.read()
                if event == sg.WINDOW_CLOSED or event == 'Quit':
                    break
                elif event == 'Reload_continue':
                    Beam = open_file(values['Reload_Ex'],KeyCheck(variable,'Project'))
                    window1.close()
                    Beam.values['Project'] = Project
                    Layouts(Beam.values)
                    break
        if event == 'calculate':
            SectionType1 = values['SectionType']
            SectionSize1 = values['SectionSize']

            alpha_m = float(values['alpha_m'])
            restraint = values['restraint']
            load_height_position = values['load_height_position']
            longitudinal_position = values['longitudinal_position']
            ends_with_restraint = values['ends_with_restraint']
            section_properties = st.section_properties(SectionType1, SectionSize1,0,0,0)
            st.compact(section_properties, SectionType1, float(values['segment length']), alpha_m, restraint, load_height_position,
                                      longitudinal_position, ends_with_restraint)
            st.shear(section_properties,SectionType1)
            st.shear_moment(section_properties,0)
            st.axial_compression(section_properties,SectionType1,float(values['segment length']))

            print(section_properties)
            # End program if user closes window or
            # presses the OK button
            # if event == 'quit' or event == sg.WIN_CLOSED:
            window['PhiMsx'].update(str(round(section_properties['PhiMsx'], 1)) + 'kNm')

            window['PhiMbx'].update(str(round(section_properties['PhiMbx'], 1)) + 'kNm')
            window['PhiMsy'].update(str(round(section_properties['PhiMsy'],1)) + 'kNm')
            window['PhiMby'].update(str(round(section_properties['PhiMby'], 1)) + 'kNm')
            window['PhiNsx'].update(str(round(section_properties['PhiNsx'], 1)) + 'kN')
            window['PhiNcx'].update(str(round(section_properties['PhiNcx'], 1)) + 'kN')
            window['PhiNsy'].update(str(round(section_properties['PhiNsy'], 1)) + 'kN')
            window['PhiNcy'].update(str(round(section_properties['PhiNcy'], 1)) + 'kN')
            window['PhiVu'].update(str(round(section_properties['PhiVu'], 2)) + ' KN')
        elif event == 'calc_T_section':
            section_properties = st.section_properties(values['SectionType'],0,values['b'],values['d'],values['t'])
            calculations = st.compact(section_properties,values['SectionType'], float(values['Length']), float(values['alpha_m']),values['restraint'],
                                      values['load_height_position'], values['longitudinal_position'],values['ends_with_restraint'])
            window['result'].update('PhiMsx = ' + str(round(calculations['PhiMsx'], 1)) + 'kNm  Actual capacity ~90% of shown here')
            window['PhiMsx'].update('PhiMbx = ' + str(round(calculations['PhiMbx'], 1)) + 'kNm  Actual capacity ~90% of shown here')
        elif event == 'dLimit':
            for x in values['dLimit'].split('/'):
                if x.isdigit():
                    Limit = int(x)
            window['dLimit1'].update(str(round(float(values['Length'])*1000/Limit,2)) + ' mm')
            if float(values['Length'])*1000/Limit > d[0]+d[1]:
                window['dCheck'].update('OK')
            else:
                window['dCheck'].update('NG')
        if event == 'SectionType':
            if values['SectionType'] == 'Universal_Beam' or values['SectionType']=='Universal_Column' or values['SectionType'] == 'PFC' or values['SectionType'] == 'RHS' or values['SectionType']=='SHS' or values['SectionType']=='CHS':
                    item = values[event]
                    matrix = []
                    SectionNames = wb[item]
                    for x in SectionNames:
                        matrix.append(x[0].value)
                    for i in range(7):
                        matrix.pop(0)
                    SectionSize = list(filter(lambda item: item is not None, SectionSize))
                    window['SectionSize'].update(value=matrix[7], values=matrix)
            elif values['SectionType'] == 'MGP10':
                SectionSize = ['90x90','140x90','190x90','240x90']
                window['SectionSize'].update(value='90x90', values=SectionSize)

        # window['b1'].update('OK')
        elif event == 'b2' or event == sg.WIN_CLOSED:
            break
        elif event == 'restraint':
            if values['restraint'] == 'FU' or values['restraint'] == 'PU':
                window['ends_with_restraint'].update(values = ['Any'], value = 'Any')
            elif values['restraint'] == 'FF' or values['restraint'] == 'FP' or values['restraint'] == 'PP':
                window['ends_with_restraint'].update(values = ['None', 'One', 'Both'], value = 'None')
            elif values['restraint'] == 'FL' or values['restraint'] == 'PL' or values['restraint'] == 'LL':
                window['ends_with_restraint'].update(values = ['None'], value = 'None')

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
            st.printcalcs(SectionSize1,section_properties,values['print_location'],values['job_name'],'',0,0,0,True)
    window.close()


def steel_calculator():

    # load in data to be read
    from openpyxl import load_workbook

    wb = load_workbook(filename='Steel Design Calculator.xlsx')
    # SectionType = wb.sheetnames
    SectionSize = []
    # print(sheet_ranges['A6'].value)
    SectionType = ['Universal_Beam', 'Universal_Column', 'PFC', 'RHS','SHS','CHS', 'T-section']
    for j in wb['Universal_Beam']:
        SectionSize.append(j[0].value)
    SectionSize = list(filter(lambda item: item is not None, SectionSize))
    print(SectionSize)
    for j1 in range(4):
        SectionSize.pop(0)
    print(SectionSize)
    # load in GUI library

    # sg.theme('LightBrown11')
    # ttk_style = 'vista'
    layout = [
        [sg.Column([
            [sg.Text("Steel calculator", key='title')],
            [sg.Text('This calculator is not finished yet', key='t1')],
            [sg.Text('Available section types are Universal Columns, Universal Beams, and PFC\'s and RHS')],
            #[sg.Text('Input Total Length in metres below:', key='t2')],
            #[sg.Input(key='Length', size=(5, 1), default_text='1')],
            [sg.Text('Input segment Length in metres below:')],
            [sg.Input(default_text='1', key='segment length', size=(5, 1))],
            [sg.Text('Input alpha m below:', key='t3')],
            [sg.Input(key='alpha_m', size=(5, 1), default_text='1')],
            [sg.Combo(['FF', 'FP', 'FL', 'PP', 'PL', 'LL', 'FU', 'PU'], key='restraint', default_value='FP',
                      enable_events=True)],
            [sg.Combo(['Shear centre', 'Top flange'], key='load_height_position', default_value='Shear centre',
                      size=(15, 1))],
            [sg.Combo(['Within segment', 'At segment end'], key='longitudinal_position', default_value='Within segment',
                      size=(15, 1)),sg.Text('Load height')],
            [sg.Combo(['Any', 'None', 'One', 'Both'], key='ends_with_restraint', default_value='One'),sg.Text('Ends with Lateral restraint')],
            [sg.Button('Calculate', key='calculate', use_ttk_buttons=True)],
            [sg.Button('Back', key='back'), sg.Button('quit', key='b2', use_ttk_buttons=True)],
            [sg.Column([
                [sg.Text('\u03A6Msx = '),sg.Text('', key='PhiMsx')],
                [sg.Text('\u03A6Mbx = '),sg.Text('', key='PhiMbx')],
                [sg.Text('\u03A6Msy = '),sg.Text('', key='PhiMsy')],
                [sg.Text('\u03A6Mby = '),sg.Text('', key='PhiMby')]
            ]),sg.Column([
                [sg.Text('\u03A6Nsx = '),sg.Text('', key = 'PhiNsx')],
                [sg.Text('\u03A6Ncx = '),sg.Text('',key = 'PhiNcx')],
                [sg.Text('\u03A6Nsy = '),sg.Text(key = 'PhiNsy')],
                [sg.Text('\u03A6Ncy = '),sg.Text(key = 'PhiNcy')]

            ]),sg.Column([
                [sg.Text('\u03A6Vu = '), sg.Text('', key='PhiVu')],
                [sg.Text('\u03A6Vvm = '), sg.Text('', key='PhiVvm')]
            ])],
            [sg.Button('Print calculations', key='print_calcs')],
            [sg.Text('Type Job Name:'), sg.Input(default_text='610UB125', key='job_name')],
            [sg.Text('Choose destination:'),
             sg.Input(key='print_location', default_text=r'C:\Users\tduffett\PycharmProjects\pythonProject1'),
             sg.FolderBrowse()]]),
            sg.VSeparator()

           # ,sg.Column(
            #[
             #   [sg.Canvas(key='canvas')],
             #   [sg.Button('Calculate',key='deflection'),sg.Text('Point Loads:'),sg.Input(default_text='1 1',key = 'PL',size=(10,1))
             #       ,sg.Text('format:[Load,position,etc.]   Length'),sg.Input(default_text='2 2',key = 'L',size=(10,1)),
             #    sg.Text('format: [span1,span2,etc.]   UDL:'),sg.Input(default_text='1 2 3',key = 'UDL',size=(10,1)),sg.Text('format:[start,end,load]')]
            #]
            #)
        ]]


    # Create the window
    return layout

    # Create an event loop


variable = variables('Loading')
Layouts(variable)
Beam = open_file('FB6A')
print(getattr(Beam.G1,'R2'))
print(getattr(Beam.G2,'R2'))
print(getattr(Beam.Q,'R2'))


x = np.linspace(0, float(Beam.Length + Beam.values['CLength']), int(Beam.Length * Beam.sample_rate + Beam.values['CLength']*Beam.sample_rate + 1))
plt.plot(x, Beam.Q.Shear)
plt.show()