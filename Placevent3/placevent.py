#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Created by Daniel Sindhikara, sindhikara@gmail.com

Placevent

This program is designed to automatically place explicit solvent atoms/ions
based on 3D-RISM data. The 3D-RISM correlations should be in a DX file.

The details of the algorithm are described in: 

Placevent: An algorithm for prediction of explicit solvent atom distribution --
Application to HIV-1 protease and F-ATP synthase 

Daniel J. Sindhikara, Norio Yoshida, Fumio Hirata
http://dansindhikara.com/Software/Entries/2012/6/22_Placevent_New.html

This package requires Numpy.

g(r)_0 is printed in the occupation column
g(r)_i is printed in the beta column

    Copyright (C) 2012 Daniel J. Sindhikara

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

TODO: 
Add bulk water to avoid evacuation errors
Switch to Argparse-based method
'''

import sys
from math import *

import numpy as np

from Placevent3 import grid
from Placevent3 import pdb


def converttopop(distribution, delta, conc):
    '''Convert distribution fn to population fn using concentration.

    Returns: 
    1) Population function (numpy array)
    2) Expected total population in grid (float)
    3) Single voxel volume (float) 
    '''

    xlen = len(distribution)
    ylen = len(distribution[0])
    zlen = len(distribution[0][0])
    popzero = [[[0.0 for z in range(zlen)] for y in range(ylen)]
               for x in range(xlen)]
    gridvolume = delta[0] * delta[1] * delta[2]
    print('# Volume of each voxel is ', gridvolume, ' cubic angstroms')

    # z fast, y med, x slow

    totalpop = 0
    for i in range(xlen):
        for j in range(ylen):
            for k in range(zlen):
                popzero[i][j][k] = distribution[i][j][k] * conc \
                    * gridvolume
                totalpop += popzero[i][j][k]
    print('# The population is approximately ', totalpop, \
          '  within the given grid')
    return (popzero, totalpop, gridvolume)


#! modified by S.aono
#def doplacement(
#    popzero,
#    conc,
#    gridvolume,
#    origin,
#    delta,
#    shellindices,
#    grcutoff
#    ):
def doplacement(
    popzero,
    conc,
    gridvolume,
    origin,
    delta,
    shellindices,
    grcutoff,
    buffer_r,
    radius_lig,
    cen_lig
    ):
#! end modified by S.aono
    '''Does iterative portion of placement according to Placevent Algorithm

    Returns:
    NumPy array of "Atom" Class objects
    '''

    popi = np.array(popzero)  # popi (popzero) are instantaneous (initial) pops.
    bulkvoxelpop = gridvolume * conc
    max = 0
    topindices = [[[]]]
    print('# Initial 3D-RISM g(r) will be printed in occupancy column')
    placedcenters = []
    finished = 0
    print('# Doing placement...')
    index = 0
#! added by S.aono
    index_err = 0
    buffer_sq = buffer_r ** 2
#! end added by S.aono
    while finished == 0:
        index += 1
        mymax = popi.max()
        (maxi, maxj, maxk) = np.argwhere(popi == mymax)[0]
#! added by S.aono
        coord_current = np.array( [ maxi * delta[0] + origin[0],
                                    maxj * delta[1] + origin[1],
                                    maxk * delta[2] + origin[2] ] )
        ncen = len(placedcenters)
        coord_previous = []
        for icen in range(ncen):
            coord_previous.append( placedcenters[icen].coord )
        coord_previous = np.array( coord_previous )
        ierror = 0
        for icen in range(ncen):
            vec = coord_current - coord_previous[icen]
            dist = np.dot(vec, vec)
            if dist <= buffer_sq:
                ierror = -1
                break
        if not ierror == 0:
            popi[ maxi, maxj, maxk ] = popi[ maxi, maxj, maxk ] * 0.9
            index -= 1
            continue
#! end added by S.aono
        topindices.append([maxi, maxj, maxk])
        remainder = 1  # Remaining population to subtract
        indexradius = 0  # First shell radius
        while remainder > 0:
            availablepopulation = 0
            try:
                for indices in shellindices[indexradius]:  # Loop on shell.
                    try:
                        availablepopulation += popi[maxi + indices[0],
                                maxj + indices[1], maxk + indices[2]]
                    except IndexError:
                        availablepopulation += bulkvoxelpop # Steal from pseudobulk
                        #print('# Stopping! Population search went outside' \
                        #      ' grid.\n')
                        #print('# You may want to try again with a higher ' \
                        #      'concentration\n')
                        #remainder = -1  # Cue quitting search.
                        #finished = 1
                        #break
            except IndexError:
#! added by S.aono
                if not index == index_err:
                    index_err = index
                    index -= -1
                    popi[ maxi, maxj, maxk ] = 0.0
                    break
#! end added by S.aono
                print(' # Stopping!\n')
                print(' # Population search went outside maximum shell size.\n')
                print(' # This usually happens with a dilute concentration.\n')
                print(' # You may want to try again with a higher concentration.\n')
                print(' # Your results may still be reasonable.')
                remainder = -1  # Cue quitting search.
                finished = 1
                break
            if availablepopulation > remainder and remainder > -1:

                # There is enough in this shell to bring the remainder to zero.
                # Is this fraction right?
                fraction = (availablepopulation - remainder) / availablepopulation
                remainder = 0
                for indices in shellindices[indexradius]:
                    try:
                        popi[maxi + indices[0]][maxj + indices[1]][maxk
                            + indices[2]] = popi[maxi
                            + indices[0]][maxj + indices[1]][maxk
                            + indices[2]] * fraction
                    except IndexError:
                        pass # Don't need to adjust population outside grid
            elif remainder > -1:
                # Not enough remainder in this shell, set whole shell pop to zero.

                for indices in shellindices[indexradius]:
                    try:
                        popi[maxi + indices[0]][maxj + indices[1]][maxk
                            + indices[2]] = 0
                    except IndexError:
                        pass # Don't need to adjust population outside grid
                remainder -= availablepopulation
                indexradius += 1
                if indexradius > len(shellindices):  # Search is larger than maxshells
                    print(' # Stopping!\n')
                    print(' # Population search went outside grid.\n')
                    print(' # You may want to try again with a higher ' \
                          '  concentration.')
                    exit()

        myg0 = popzero[maxi][maxj][maxk] / (conc * gridvolume)  # Original g(r).
        mygi = mymax / (conc * gridvolume)  # g(r) when placed.
        placedcenters.append(pdb.Atom(serial=len(placedcenters),
                             resseq=index, coord=[maxi * delta[0]
                             + origin[0], maxj * delta[1] + origin[1],
                             maxk * delta[2] + origin[2]], occ=myg0,
                             tfac=mygi))
        if myg0 < grcutoff:
            finished = 1  # Already placed all atoms > grcutoff.
    return placedcenters


#! modified by S.aono
#def returncenters(guvfilename, molar, grcutoff):
def returncenters(guvfilename, molar, grcutoff, buffer_r, radius_lig, cen_lig):
#! end modified by S.aono

    '''Given distribution, returns placed atoms

    Input: guvfilename and molarity
    Returns: Placed centers as Atom objects
    '''

    conc = molar * 6.0221415E-4
    shellindices = grid.readshellindices()
    if guvfilename[-3:] == ".dx":
        (distribution, origin, delta, gridcount) = \
            grid.readdx(guvfilename)
    elif guvfilename[-3:] == ".h5":
        grids = grid.h5ToGrids(guvfilename)
        print("Warning, assuming target molecule is 'O', please contact the developer for more info.")
        if 'guv' not in grids['O'].keys():
            exit("Error, guv distribution not found in file")
        gGrid = grids['O']['guv']
        distribution = gGrid.distribution
        origin = gGrid.origin
        delta = gGrid.deltas
        gridcount = gGrid.gridcount
    else:
        exit("Error, incompatible file type! -> {0}".format(guvfilename))

#! added by S.aono
    if radius_lig > 0.0:
        nmin_i = int(np.floor((cen_lig[0] - radius_lig) / delta[0]))
        nmin_j = int(np.floor((cen_lig[1] - radius_lig) / delta[1]))
        nmin_k = int(np.floor((cen_lig[2] - radius_lig) / delta[2]))
        nmax_i = int(np.ceil((cen_lig[0] + radius_lig) / delta[0]))
        nmax_j = int(np.ceil((cen_lig[1] + radius_lig) / delta[1]))
        nmax_k = int(np.ceil((cen_lig[2] + radius_lig) / delta[2]))
        nmin_i = max( [nmin_i, 0] )
        nmin_j = max( [nmin_j, 0] )
        nmin_k = max( [nmin_k, 0] )
        nmax_i = min( [nmax_i, gridcount[0]-1] )
        nmax_j = min( [nmax_j, gridcount[1]-1] )
        nmax_k = min( [nmax_k, gridcount[2]-1] )
        gridcount[0] = nmax_i - nmin_i + 1
        gridcount[1] = nmax_j - nmin_j + 1
        gridcount[2] = nmax_k - nmin_k + 1
        origin[0] = delta[0] * nmin_i
        origin[1] = delta[1] * nmin_j
        origin[2] = delta[2] * nmin_k
        tmp = np.zeros( gridcount[0] *gridcount[1] *gridcount[2], \
                        dtype=float ).reshape(\
                       [gridcount[0], gridcount[1], gridcount[2]])
        for ix in range(gridcount[0]):
            for iy in range(gridcount[1]):
                tmp[ix][iy] = distribution[nmin_i+ix][nmin_j+iy][nmin_k:nmax_k+1]
        distribution = tmp
    print('# BOX_X = {:}, {:}'.format(origin[0], \
                                      (gridcount[0]-1) * delta[0] + origin[0]))
    print('# BOX_Y = {:}, {:}'.format(origin[1], \
                                      (gridcount[1]-1) * delta[1] + origin[1]))
    print('# BOX_Z = {:}, {:}'.format(origin[2], \
                                      (gridcount[2]-1) * delta[2] + origin[2]))
#! end added by S.aono
    (popzero, totalpop, gridvolume) = converttopop(distribution,
            delta, conc)
#! modified by S.aono
#    return doplacement(popzero, conc, gridvolume, origin, delta, shellindices,
#                       grcutoff)
    return doplacement(popzero, conc, gridvolume, origin, delta, shellindices,
                       grcutoff, buffer_r, radius_lig, cen_lig)
#! end modified by S.aono


def main():
    print('# placevent.py <.dx or MDF .h5 file> <concentration M (molar)> [cutoff g(r) ' \
          ' (default 1.5)]\n', flush=True)
    print('# Concentration (#/A^3) ~= molarity * 6.0221415E-4', flush=True)

    if len(sys.argv) < 3:
        print('Insufficient arguments need 2 : ', len(sys.argv) - 1, flush=True)
        exit()
#! modified by S.aono
#    if len(sys.argv) == 4:
    if len(sys.argv) >= 4:
        grcutoff = float(sys.argv[3])
    else:
        grcutoff = 1.5
    if len(sys.argv) >= 5:
        buffer_r = float(sys.argv[4])
    else:
        buffer_r = 0.0
    if len(sys.argv) == 9:
        radius_lig = float(sys.argv[5])
        cen_lig = np.array(
            [float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])]
        )
    else:
        radius_lig = -1.0
        cen_lig = np.zeros( 3, dtype=float )

    print('# grcutoff = {:}, buffer_r = {:}'.format(grcutoff, buffer_r), flush=True)
    print(
        '# radius_lig = {:}, cen_lig = {:}, {:}, {:}'\
        .format(radius_lig, cen_lig[0], cen_lig[1], cen_lig[2], flush=True)
    )
#! end modified by S.aono
    guvfilename = sys.argv[1]
    molar = float(sys.argv[2])
    print('# Your .dx file is', guvfilename, 'molarity is', molar, 'M', flush=True)
    placedcenters = returncenters(guvfilename, molar, grcutoff, buffer_r, radius_lig, cen_lig)
    for center in placedcenters:
        print(str(center)[:-2], flush=True)  # [:-2] is to get rid of the '\n'

if __name__ == '__main__':
    main()
