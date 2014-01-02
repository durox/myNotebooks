#!/usr/bin/env python
# encoding: utf-8


import numpy as np
import sympy as sp
import IPython.display as d


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
