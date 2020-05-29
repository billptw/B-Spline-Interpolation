""""
CE7453 Assignment 1
Implementation of cubic B-spline interpolation algorithm
Input file contains a set of 2D points, outputs a cubic B-spline curve, and displays the B-spline curve

@author: Pung Tuck Weng, G1803000G
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import *


def read_input(filename):
    input_file = open(filename + '.txt', 'r').readlines()
    coord = list()
    for line in input_file:
        coord.append([float(i) for i in line.split()])
    return coord


def write_output(spline, filename):
    output_file = open(filename + '.txt', 'w')
    for line in spline:
        output_file.write(str(line))
    output_file.close()


def get_N_value(u, i, t, knots):
    def Nfirst(u, i, t, knots):
        num = (t[u] - knots[i]) ** 3
        denom = (knots[i + 1] - knots[i]) * (knots[i + 2] - knots[i]) * (knots[i + 3] - knots[i])
        if denom == 0:
            return Nsecond(u, i, t, knots)
        return num / denom

    def Nsecond(u, i, t, knots):
        num1 = ((t[u] - knots[i]) ** 2) * (knots[i + 2] - t[u])
        denom1 = (knots[i + 2] - knots[i + 1]) * (knots[i + 3] - knots[i]) * (knots[i + 2] - knots[i])
        num2 = (knots[i + 3] - t[u]) * (t[u] - knots[i]) * (t[u] - knots[i + 1])
        denom2 = (knots[i + 2] - knots[i + 1]) * (knots[i + 3] - knots[i + 1]) * (knots[i + 3] - knots[i])
        num3 = (knots[i + 4] - t[u]) * (t[u] - knots[i + 1]) ** 2
        denom3 = (knots[i + 2] - knots[i + 1]) * (knots[i + 4] - knots[i + 1]) * (knots[i + 3] - knots[i + 1])
        if denom1 == 0 or denom2 == 0 or denom3 == 0:
            return Nthird(u, i, t, knots)
        return num1 / denom1 + num2 / denom2 + num3 / denom3

    def Nthird(u, i, t, knots):
        num1 = (t[u] - knots[i]) * (knots[i + 3] - t[u]) ** 2
        denom1 = (knots[i + 3] - knots[i + 2]) * (knots[i + 3] - knots[i + 1]) * (knots[i + 3] - knots[i])
        num2 = (knots[i + 4] - t[u]) * (knots[i + 3] - t[u]) * (t[u] - knots[i + 1])
        denom2 = (knots[i + 3] - knots[i + 2]) * (knots[i + 4] - knots[i + 1]) * (knots[i + 3] - knots[i + 1])
        num3 = (knots[i + 4] - t[u]) ** 2 * (t[u] - knots[i + 2])
        denom3 = (knots[i + 3] - knots[i + 2]) * (knots[i + 4] - knots[i + 2]) * (knots[i + 4] - knots[i + 1])
        if denom1 == 0 or denom2 == 0 or denom3 == 0:
            return Nfourth(u, i, t, knots)
        return num1 / denom1 + num2 / denom2 + num3 / denom3

    def Nfourth(u, i, t, knots):
        num4 = (knots[i + 4] - t[u]) ** 3
        denom4 = (knots[i + 4] - knots[i + 3]) * (knots[i + 4] - knots[i + 2]) * (knots[i + 4] - knots[i + 1])
        if denom4 == 0.0:
            return 0.0
        return num4 / denom4

    if t[u] <= knots[i + 1] and t[u] >= knots[i]:
        return Nfirst(u, i, t, knots)
    elif t[u] <= knots[i + 2] and t[u] >= knots[i + 1]:
        return Nsecond(u, i, t, knots)
    elif t[u] <= knots[i + 3] and t[u] >= knots[i + 2]:
        return Nthird(u, i, t, knots)
    elif t[u] <= knots[i + 4] and t[u] >= knots[i + 3]:
        return Nfourth(u, i, t, knots)
    else:
        return 0.0


def get_spline_cp(coord, t, knots):
    n = len(coord) - 1

    n_matrix = np.zeros((n + 3, n + 3))
    i = 1
    j = 1
    while i < n:
        n_matrix[i + 1, j] = get_N_value(i, j, t, knots)
        n_matrix[i + 1, j + 1] = get_N_value(i, j + 1, t, knots)
        n_matrix[i + 1, j + 2] = get_N_value(i, j + 2, t, knots)
        # to check if input values are correct:
        # n_matrix[i+1, j] = i + j/10
        # n_matrix[i+1, j + 1] = i + (j+1)/10
        # n_matrix[i+1, j + 2] = i + (j+2)/10
        i += 1
        j += 1

    def diff_Nfirst(u, i, t, knots):
        x = Symbol('x')
        num = (x - knots[i]) ** 3
        denom = (knots[i + 1] - knots[i]) * (knots[i + 2] - knots[i]) * (knots[i + 3] - knots[i])
        if denom != 0 and num != 0:
            return diff(num, x, 2).subs(x, t[u]) / denom
        else:
            return 0.0

    def diff_Nsecond(u, i, t, knots):
        x = Symbol('x')
        num1 = ((x - knots[i]) ** 2) * (knots[i + 2] - x)
        denom1 = (knots[i + 2] - knots[i + 1]) * (knots[i + 3] - knots[i]) * (knots[i + 2] - knots[i])
        num2 = (knots[i + 3] - x) * (x - knots[i]) * (x - knots[i + 1])
        denom2 = (knots[i + 2] - knots[i + 1]) * (knots[i + 3] - knots[i + 1]) * (knots[i + 3] - knots[i])
        num3 = (knots[i + 4] - x) * (x - knots[i + 1]) ** 2
        denom3 = (knots[i + 2] - knots[i + 1]) * (knots[i + 4] - knots[i + 1]) * (knots[i + 3] - knots[i + 1])
        if denom1 != 0 and num1 != 0 or denom2 != 0 and num2 != 0 or denom3 != 0 and num3 != 0:
            return diff(num1, x, 2).subs(x, t[u]) / denom1 + diff(num2, x, 2).subs(x, t[u]) \
                   / denom2 + diff(num3, x, 2).subs(x, t[u]) / denom3
        else:
            return 0.0

    def diff_Nthird(u, i, t, knots):
        x = Symbol('x')
        num1 = (x - knots[i]) * (knots[i + 3] - x) ** 2
        denom1 = (knots[i + 3] - knots[i + 2]) * (knots[i + 3] - knots[i + 1]) * (knots[i + 3] - knots[i])
        num2 = (knots[i + 4] - x) * (knots[i + 3] - x) * (x - knots[i + 1])
        denom2 = (knots[i + 3] - knots[i + 2]) * (knots[i + 4] - knots[i + 1]) * (knots[i + 3] - knots[i + 1])
        num3 = (knots[i + 4] - x) ** 2 * (x - knots[i + 2])
        denom3 = (knots[i + 3] - knots[i + 2]) * (knots[i + 4] - knots[i + 2]) * (knots[i + 4] - knots[i + 1])
        if denom1 != 0 and num1 != 0 or denom2 != 0 and num2 != 0 or denom3 != 0 and num3 != 0:
            return diff(num1, x, 2).subs(x, t[u]) / denom1 + diff(num2, x, 2).subs(x, t[u]) \
                   / denom2 + diff(num3, x, 2).subs(x, t[u]) / denom3
        else:
            return 0.0

    def diff_Nfourth(u, i, t, knots):
        x = Symbol('x')
        num = (knots[i + 4] - x) ** 3
        denom = (knots[i + 4] - knots[i + 3]) * (knots[i + 4] - knots[i + 2]) * (knots[i + 4] - knots[i + 1])
        if denom != 0 and num != 0:
            return diff(num, x, 2).subs(x, t[u]) / denom
        else:
            return 0.0

    n_matrix[1, 0] = 1.0
    n_matrix[n + 1, n + 2] = 1.0
    n_matrix[0, 0] = diff_Nfourth(0, 0, t, knots)
    n_matrix[0, 1] = diff_Nthird(0, 1, t, knots)
    n_matrix[0, 2] = diff_Nsecond(0, 2, t, knots)
    n_matrix[n + 2, n] = diff_Nthird(n, n, t, knots)
    n_matrix[n + 2, n + 1] = diff_Nsecond(n, n + 1, t, knots)
    n_matrix[n + 2, n + 2] = diff_Nfirst(n, n + 2, t, knots)
    # to check if input values are correct:
    # n_matrix[0, 0] = 0.0
    # n_matrix[0, 1] = 0.1
    # n_matrix[0, 2] = 0.2
    # n_matrix[n + 2, n] = 4.4
    # n_matrix[n + 2, n + 1] = 4.5
    # n_matrix[n + 2, n + 2] = 4.6

    d_matrix = np.zeros((n + 3, 2))
    i = 1
    while i <= len(coord):
        d_matrix[i] = np.asarray(coord[i-1])
        i += 1

    p_matrix = np.linalg.solve(n_matrix, d_matrix)

    return p_matrix


def get_spline(coord):
    spline = list()
    k = 3   # degree k of cubic spline = 3
    spline.append(k)

    n = len(coord) - 1
    spline.append('\n')
    spline.append(n)

    def get_chord_parameterize(coord):
        t = np.zeros(len(coord))
        i = 1
        total = 0
        while i < len(coord):
            total += np.linalg.norm(np.asarray(coord[i]) - np.asarray(coord[i-1]))
            i += 1

        i = 1
        while i < len(coord):
            t[i] = np.linalg.norm(np.asarray(coord[i]) - np.asarray(coord[i-1])) / total + t[i-1]
            i += 1
        return t

    # optional: linear parameterization
    # t = np.linspace(0, 1, n + 1)

    t = get_chord_parameterize(coord)
    knots = np.append(np.zeros(k), t)
    knots = np.append(knots, np.ones(k))

    spline.append('\n\n')
    for val in knots:
        spline.append(val)
        spline.append('\t')

    spline_cp = get_spline_cp(coord, t, knots)

    spline.append('\n\n')
    i = 0
    while i < len(spline_cp):
        for val in spline_cp[i, :]:
            spline.append(val)
            spline.append('\t')
        spline.append('\n\n')
        i += 1

    # generating spline points for plotting
    spline_n = len(spline_cp) - 1
    points = np.linspace(knots[k], knots[spline_n + 1], n * 10)
    spline_points = []

    for point in points:
        value = 0.0
        i = 0
        while i <= spline_n:
            value += spline_cp[i] * get_N_value(0, i, [point], knots)
            i += 1
        spline_points.append(value.tolist())
    return [spline, spline_cp, spline_points]


def display_spline(input_coord, spline_cp, spline_points):
    input_coord_plot = np.array(input_coord).transpose()
    plt.plot(input_coord_plot[0], input_coord_plot[1], 'ro')
    spline_cp_plot = np.array(spline_points).transpose()
    plt.plot(spline_cp_plot[0], spline_cp_plot[1])
    spline_points_plot = np.array(spline_cp).transpose()
    plt.plot(spline_points_plot[0], spline_points_plot[1], 'gx')
    plt.legend(['Input Points', 'Interpolated Line', 'Control Points'])
    #plt.axis('equal')
    plt.show()

if __name__ == "__main__":
    input_coord = read_input('input')
    [spline, spline_cp, spline_points] = get_spline(input_coord)
    display_spline(input_coord, spline_cp, spline_points)
    write_output(spline, 'cubic')
