#!/usr/bin/env python
# encoding: utf-8


import numpy as np
import sympy as sp
import IPython.display as d


def F_eingesp(E, I, l, u, w, p0):
    """RSK eingespannt Linienlast

    :E: @todo
    :I: @todo
    :l: @todo
    :w: @todo
    :returns: tupel F1_ik, F2_ik, M3_ik

    """
    lam = l * (u * w**2 / E / I) ** (1/4)
    eps = l * np.sqrt(u * w**2 / E / A)

    o1 = (np.cosh(lam) + np.cos(lam)) / 2
    o2 = (np.sinh(lam) + np.sin(lam)) / 2
    o3 = (np.cosh(lam) - np.cos(lam)) / 2
    o4 = (np.sinh(lam) - np.sin(lam)) / 2

    F2_ik = p0*l / lam * (-o2 + o1*o2 - o3*o4) / (o3**2 - o2*o4)
    M3_ik = p0*l**2 / lam**2 * (o1*o3 - o3 - o4**2) / (o3**2 - o2*o4)

    return 0, F2_ik, M3_ik


def F_gelenk(E, I, l, u, w, p0):
    """RSK gelenkig Linienlast

    :E: @todo
    :I: @todo
    :l: @todo
    :w: @todo
    :returns: tupel F1_ki, F2_ki, M3_ki

    """
    lam = l * (u * w**2 / E / I) ** (1/4)
    eps = l * np.sqrt(u * w**2 / E / A)

    o1 = (np.cosh(lam) + np.cos(lam)) / 2
    o2 = (np.sinh(lam) + np.sin(lam)) / 2
    o3 = (np.cosh(lam) - np.cos(lam)) / 2
    o4 = (np.sinh(lam) - np.sin(lam)) / 2

    F2_ki = p0*l / lam * (o1**2 - o1 - o2*o4) / (o2*o3 - o1*o4)
    M3_ki = 0

    return 0, F2_ki, M3_ki


def Kw_eingesp_kik(E, A, I, l, u, w):
    """dyn. Steifigkeit

    :E: @todo
    :A: @todo
    :I: @todo
    :l: @todo
    :u: @todo
    :w: @todo
    :returns: @todo

    """
    lam = l * (u * w**2 / E / I) ** (1/4)
    eps = l * np.sqrt(u * w**2 / E / A)

    o1 = (np.cosh(lam) + np.cos(lam)) / 2
    o2 = (np.sinh(lam) + np.sin(lam)) / 2
    o3 = (np.cosh(lam) - np.cos(lam)) / 2
    o4 = (np.sinh(lam) - np.sin(lam)) / 2

    kik = np.matrix([[E*A/l * eps/np.tan(eps), 0, 0],
                     [0, E*I*lam**3/l**3 * (o1*o2 - o3*o4)/(o3**2 - o2*o4), E*I*lam**2/l**2 * (o4**2 - o1*o3)/(o3**2 - o2*o4)],
                     [0, E*I*lam**2/l**2 * (o4**2 - o1*o3)/(o3**2 - o2*o4), E*I*lam/l * (o2*o3 - o1*o4)/(o3**2 - o2*o4)]])

    return kik


def Kw_gelenk_kik(E, A, I, l, u, w):
    """dyn. Steifigkeit Gelenk an k

    :E: @todo
    :A: @todo
    :I: @todo
    :l: @todo
    :u: @todo
    :w: @todo
    :returns: @todo

    """
    lam = l * (u * w**2 / E / I) ** (1/4)
    eps = l * np.sqrt(u * w**2 / E / A)

    o1 = (np.cosh(lam) + np.cos(lam)) / 2
    o2 = (np.sinh(lam) + np.sin(lam)) / 2
    o3 = (np.cosh(lam) - np.cos(lam)) / 2
    o4 = (np.sinh(lam) - np.sin(lam)) / 2

    kik = np.matrix([[E*A/l * eps/np.tan(eps), 0, 0],
                     [0, E*I*lam**3/l**3 * (o1**2 - o2*o4)/(o2*o3 - o1*o4), 0],
                     [0, 0, 0]])

    return kik


def KT_eingesp(E, A, I, l, u, element='all'):
    """appr. Trennung, statischer Teil, eingesp Stab"""
    lam2 = l**2 / (I / A)
    K = np.matrix([[lam2, 0, 0, -lam2, 0, 0],
                  [0, 12, 6*l, 0, -12, 6*l],
                  [0, 6*l, 4*l**2, 0, -6*l, 2*l**2],
                  [-lam2, 0, 0, lam2, 0, 0],
                  [0, -12, -6*l, 0, 12, -6*l],
                  [0, 6*l, 2*l**2, 0, -6*l, 4*l**2]]) * E*I/l**3
    if element=='iki':
        return K[0:3, 0:3]
    elif element=='ikk':
        return K[0:3, 3:]
    elif element=='kii':
        return K[3:, 0:3]
    elif element=='kik':
        return K[3:, 3:]
    elif element=='all':
        return K


def MT_eingesp(E, A, I, l, u, element='all'):
    """appr. Trennung, dynamischer Teil, eingespannter Stab"""
    K = np.matrix([[140, 0, 0, 70, 0, 0],
                  [0, 156, 22*l, 0, 54, -13*l],
                  [0, 22*l, 4*l**2, 0, 13*l, -3*l**2],
                  [70, 0, 0, 140, 0, 0],
                  [0, 54, 13*l, 0, 156, -22*l],
                  [0, -13*l, -3*l**2, 0, -22*l, 4*l**2]]) * u*l/420
    if element=='iki':
        return K[0:3, 0:3]
    elif element=='ikk':
        return K[0:3, 3:]
    elif element=='kii':
        return K[3:, 0:3]
    elif element=='kik':
        return K[3:, 3:]
    elif element=='all':
        return K


def KT_gelenk(E, A, I, l, u, element='all'):
    """appr. Trennung, statischer Teil, gelenkiger Stab"""
    lam2 = l**2 / (I / A)
    K = np.matrix([[lam2, 0, 0, -lam2, 0, 0],
                  [0, 3, 3*l, 0, -3, 0],
                  [0, 3*l, 3*l**2, 0, -3*l, 0],
                  [-lam2, 0, 0, lam2, 0, 0],
                  [0, -3, -3*l, 0, 3, 0],
                  [0, 0, 0, 0, 0, 0]]) * E*I/l**3
    if element=='iki':
        return K[0:3, 0:3]
    elif element=='ikk':
        return K[0:3, 3:]
    elif element=='kii':
        return K[3:, 0:3]
    elif element=='kik':
        return K[3:, 3:]
    elif element=='all':
        return K


def MT_gelenk(E, A, I, l, u, element='all'):
    """appr. Trennung, dynamischer Teil, gelenkiger Stab"""
    K = np.matrix([[280, 0, 0, 140, 0, 0],
                  [0, 630, 72*l, 0, 117, 0],
                  [0, 72*l, 16*l**2, 0, 33*l, 0],
                  [140, 0, 0, 280, 0, 0],
                  [0, 117, 33*l, 0, 198, 0],
                  [0, 0, 0, 0, 0, 0]]) * u*l/840
    if element=='iki':
        return K[0:3, 0:3]
    elif element=='ikk':
        return K[0:3, 3:]
    elif element=='kii':
        return K[3:, 0:3]
    elif element=='kik':
        return K[3:, 3:]
    elif element=='all':
        return K


def K_eingesp(E, A, I, l, element='all'):
    """Stabsteifigkeitsmatrix für beidseitig eingespannten Stab"""
    K = np.matrix([[E*A/l, 0, 0, -E*A/l, 0, 0],
                  [0, 12*E*I/(l**3), 6*E*I/(l**2), 0, -12*E*I/(l**3), 6*E*I/(l**2)],
                  [0, 6*E*I/(l**2), 4*E*I/l, 0, -6*E*I/(l**2), 2*E*I/l],
                  [-E*A/l, 0, 0, E*A/l, 0, 0],
                  [0, -12*E*I/(l**3), -6*E*I/(l**2), 0, 12*E*I/(l**3), -6*E*I/(l**2)],
                  [0, 6*E*I/(l**2), 2*E*I/l, 0, -6*E*I/(l**2), 4*E*I/l]])
    if element=='iki':
        return K[0:3, 0:3]
    elif element=='ikk':
        return K[0:3, 3:]
    elif element=='kii':
        return K[3:, 0:3]
    elif element=='kik':
        return K[3:, 3:]
    elif element=='all':
        return K


def K_gelenk(E, A, I, l, element='all'):
    """Stabsteifigkeitsmatrix für rechtsseitigges Momentennullfeld"""
    K = np.matrix([[E*A/l, 0, 0, -E*A/l, 0, 0],
                  [0, 3*E*I/(l**3), 3*E*I/(l**2), 0, -3*E*I/(l**3), 0],
                  [0, 3*E*I/(l**2), 3*E*I/l, 0, -3*E*I/(l**2), 0],
                  [-E*A/l, 0, 0, E*A/l, 0, 0],
                  [0, -3*E*I/(l**3), -l*E*I/(l**2), 0, 3*E*I/(l**3), 0],
                  [0, 0, 0, 0, 0, 0]])
    if element=='iki':
        return K[0:3, 0:3]
    elif element=='ikk':
        return K[0:3, 3:]
    elif element=='kii':
        return K[3:, 0:3]
    elif element=='kik':
        return K[3:, 3:]
    elif element=='all':
        return K


def T(alpha):
    """Transformationsmatrix in Abhängigkeit von Drehwinkel Alpha (rad)"""
    c = np.round(np.cos(alpha), 8)
    s = np.round(np.sin(alpha), 8)
    T = np.matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    return T


def printmatrix(m, symbol=None):
    if len(m.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(m).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv += [r'\end{bmatrix}']
    r = '\n'.join(rv)
    if symbol:
        r = symbol + ' = ' + r
    return d.Latex('$$' + r + '$$')


def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)


class Output(object):

    """Docstring for Output. """

    def __init__(self, decimals=1):
        """@todo: to be defined1.

        :decimals: @todo

        """
        self.elements = list()
        self.decimals = decimals

    def addMatrix(self, value, name=None):
        """@todo: Docstring for addMatrix.

        :value: @todo
        :name: @todo
        :returns: @todo

        """
        if name:
            s = name + ' = '
        else:
            s = ''
        s += bmatrix(np.round(value, self.decimals))
        self.elements.append('$$' + s + '$$')

    def addSPMatrix(self, value, name=None):
        """@todo: Docstring for addMatrix.

        :value: @todo
        :name: @todo
        :returns: @todo

        """
        if name:
            s = name + ' = '
        else:
            s = ''
        s += bmatrix(value)
        self.elements.append('$$' + s + '$$')

    def addValue(self, value, name=None, unit=None):
        """@todo: Docstring for addValue.

        :value: @todo
        :name: @todo
        :unit: @todo
        :returns: @todo

        """
        if name:
            s = name + ' = '
        else:
            s = ''
        s += str(np.round(value, self.decimals))
        if unit:
            s += ' ' + unit
        self.elements.append('$$' + s + '$$')

    def addText(self, text):
        """@todo: Docstring for addText.

        :text: @todo
        :returns: @todo

        """
        self.elements.append('$$' + text + '$$')

    def latex(self):
        """@todo: Docstring for latex.
        :returns: @todo

        """
        return d.Latex('\n'.join(self.elements))
