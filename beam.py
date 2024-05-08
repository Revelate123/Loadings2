import scipy.integrate as integrate
import numpy as np

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
            elif a * self.sample_rate < i < (a + b) * self.sample_rate:
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

