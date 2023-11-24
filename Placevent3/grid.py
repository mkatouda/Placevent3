# -*- coding: utf-8 -*-

#!/usr/bin/python

'''
Daniel J. Sindhikara
sindhikara@gmail.com
Handles 3D volumetric data
Copyright 2013 Daniel Jon Sindhikara
'''

from __future__ import division
import numpy as np
import os

class Grid:
    '''
Contains volumetric data
    '''

    def __init__(
        self,
        distribution,
        origin,
        gridcount,
        deltas,
        concentration=-1.0,
        ):
        if type(distribution) is list:
            self.distribution = np.array(distribution)
        elif type(distribution) is np.ndarray:
            self.distribution = distribution
        self.origin = origin
        self.gridcount = gridcount
        self.deltas = deltas
        self.concentration = concentration  # in molar

    def getvalue(self, coord): 
        return linearinterpolatevalue(self.distribution, self.origin,
                self.deltas, coord)

    def writedx(self, filename):
        printdxfrom3d(self.distribution, self.origin, self.deltas,
                      self.gridcount, filename)
    def writegrd(self, filename):
        printgrd(self.distribution, self.origin, self.deltas,
                      self.gridcount, filename)

    def nearestgridindices(self, coord):
        '''
        Given a 3D cartesian coordinate, return nearest grid indices
        '''

        gridindices = [int(round((coord[i] - self.origin[i])
                       / self.deltas[i])) for i in range(3)]
        return gridindices

    def coordfromindices(self, indices):
        return [indices[i] * self.deltas[i] + self.origin[i] for i in
                range(3)]

    def coarseRDF(self, coord):
        '''
        Given a 3D cartesian coordinate, return 1D distribution as list, spaced by gridspacing
        '''

        dist = []
        myindices = self.nearestgridindices(coord)
        for shellindex in shellindices:
            avg = 0.0
            try:
                for indexonshell in shellindex:
                    avg += self.distribution[myindices[0]
                            + indexonshell[0]][myindices[0]
                            + indexonshell[1]][myindices[2]
                            + indexonshell[2]]
                avg = avg / len(shellindex)
                dist.append(avg)
            except IndexError:
                break
        return dist

    def interpRDF(
        self,
        coord,
        delta,
        limit,
        concentration = False # If not False, print coordination number
        ):
        '''
        Given 3D cartesian coordinate (list-like object of floats),
        delta (float), and limit (float),
        
        return RDF

        by averaging linear interpolated g(r_vec) on points on spherical shell.
        '''

        spherefilename = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)),
                         'data', 'points', '200.pts')
        spherepoints = [[float(element) for element in line.split()]
                        for line in open(spherefilename).readlines()]
        rdf = []
        volints = [0.]
        for radius in np.arange(0, limit, delta):
            mysum = 0
            for point in spherepoints:
                mysum += self.getvalue([point[dim] * radius
                        + coord[dim] for dim in range(3)])
            rdf.append(float(mysum / len(spherepoints)))
            if radius > 0:
                volints.append(volints[-1] + 4./3. * 3.14159 * (radius**3 - (radius-delta)**3))
        if not concentration:
            return rdf
        else:
            coordnums = [volint * concentration for volint in volints]
            return rdf, coordnums


def rdf23dgrid(
    rdf,
    rdfdelta,
    griddelta,
    gridorigin,
    shellindices,
    ):
    '''
    Convert 1d RDF to 3D grid.
    '''

    gridcount = [len(shellindices) * 2 + 1] * 3  # + 1 is to make it odd, let
    mydist = np.array([[[1.0 for i in range(gridcount[0])] for j in
                      range(gridcount[1])] for k in
                      range(gridcount[2])])

    indexmultiplier = int(griddelta / rdfdelta)  # from 3D grid index to 1D rdf index
    for (i, shell) in enumerate(shellindices):  # loop through shells
        myrdf = rdf[i * indexmultiplier]
        for shellindex in shell:
            distindex = tuple([shellindex[dim] + int(gridcount[dim]
                              / 2) for dim in range(3)])

            # above will put shell 0 at the int(gridcount/2) aka numshells/2 + 1

            mydist[distindex] = myrdf
    return Grid(mydist, gridorigin, gridcount, [griddelta] * 3)


def dx2Grid(dxfilename):
    '''
    Reads a dx file into a Grid class object 
    '''

    # for now Ghetto rigged to use old readdx function

    (distribution, origin, deltas, gridcount) = readdx(dxfilename)
    return Grid(distribution, origin, gridcount, deltas)


def openreadfile(filename):
    '''
Open a file for reading.
Use gzip.open if it's a .gz file.
    '''
    from gzip import open as gopen
    if 'gz' in filename:
        f = gopen(filename,'rb')
    else:
        f = open(filename,'r')
    return(f)    


def readdx(filename):
    '''
    Reads dx files into memory.
    Returns grid data and single 4D array containing dx data [index][xindex][yindex][zindex]
    and the origin, deltas and gridcounts.
    '''
    dxfile = openreadfile(filename)
    dxlines = []
    for i in range(10):  # only need the first few lines to get grid data
        dxlines.append(dxfile.readline())
    dxfile.close()
    gridcounts = []
    origin = []
    deltas = [0, 0, 0]
    startline = 0
    for (i, line) in enumerate(dxlines):
        splitline = line.split()
        if len(splitline) > 2:
            if splitline[1] == '1':
                gridcounts.append(int(splitline[5]))
                gridcounts.append(int(splitline[6]))
                gridcounts.append(int(splitline[7]))

                # print("# gridcounts ",gridcount)

            if splitline[0] == 'origin':
                origin.append(float(splitline[1]))
                origin.append(float(splitline[2]))
                origin.append(float(splitline[3]))
            if splitline[0] == 'delta':
                if float(splitline[1]) > 0:
                    deltas[0] = float(splitline[1])
                if float(splitline[2]) > 0:
                    deltas[1] = float(splitline[2])
                if float(splitline[3]) > 0:
                    deltas[2] = float(splitline[3])
            if splitline[1] == '3':
                numpoints = int(splitline[9])

                # print("# Total number of gridpoints is ",numpoints)

                startline = i + 1
        if startline > 1:
            break
    
    # Read distribution Values
    dxfile = openreadfile(filename)
    dxtext = dxfile.read()
    dxfile.close()
    splittext = dxtext.split()
    del splittext[0:splittext.index('follows') + 1]  # get rid of header text, last word is "follows"
    floats = []
    for element in splittext:
        if len(element) > 0:
            try:
                floats.append(float(element))
            except ValueError:
                pass
    
    # Assign to 3D numpy array
    assert len(floats) == gridcounts[0]*gridcounts[1]*gridcounts[2]
    import numpy as np
    distribution = np.array(floats)
    distribution = np.reshape(distribution, gridcounts)

    return (distribution, origin, deltas, gridcounts)


def readUxDATA(filename, disttypes = ['g']):
    '''
    Reads 'UxDATA' (UVDATA or UUDATA) files into memory.
    Returns dictionary of 3D numpy arrays, and the origin, deltas and gridcounts.
    Dictionary keys refer to <speciesname>.<distributiontype>
    

    Since these files contain many distribution types, the type
    of distributions to keep should be specified. Default 'g' -> g(r)


    First line should look something like: 
    ## <version> <gridcounts> <gridspacing> <number of species>
    ## 20130619        64     0.5     2
   
    '''
    import numpy as np
    distfile = openreadfile(filename)
    filelines = distfile.readlines()
    linesplit = filelines[0].split()
    if len(linesplit) != 5:
        exit("Error! Cannot parse first line of UxDATA file")
    gridcounts = [int(linesplit[2]) for dim in range(3)]
    deltas = [float(linesplit[3]) for dim in range(3)]
    numspecies = int(linesplit[4])

    distcolumns = {}
    origin = False
    for line in filelines[:10]: # Should contain at least one label line
        if 'g(r)' in line:
            labels = line.split()
            for disttype in disttypes:
                for lindex, label in enumerate(labels):
                    if disttype in label[0] and 'k' not in label:
                        distcolumns[disttype] = lindex - 1
        if not origin and '#' not in line:
            origin = [float(element) for element in line.split()[:3]]
                        
    dists = {}
    for linenum, line in enumerate(filelines):
        splitline = line.split()
        if '##<' in line: #New species
            speciesname = splitline[1]
            print("Reading species:", speciesname)
            for key in distcolumns.keys():
                #print key
                distname = speciesname+'.'+key
                dists[distname] = []
                #print distname
        if "#" not in line and len(splitline) > 2:
            # Probably a data line
            for key, colnum in distcolumns.iteritems():
                distname = speciesname+'.'+key
                dists[distname].append(float(splitline[colnum]))
            
    for key, value in dists.iteritems():
        assert len(value) == gridcounts[0]*gridcounts[1]*gridcounts[2]
        myarray = np.array(value)
        dists[key] = np.reshape(myarray, gridcounts, order='F')
    return dists, origin, deltas, gridcounts

def data2Grids(uxdatafilename, disttypes=['g']):
    '''
    Reads a UxDATA into Grid class objects 
    '''
    (distributions, origin, deltas, gridcounts) = readUxDATA(uxdatafilename, disttypes=disttypes)
    grids = {}
    for name, dist in distributions.iteritems():
        grids[name] = Grid(dist, origin, gridcounts, deltas)
    return grids


def readTKRguv(filename):
    '''
    Reads TINKER 'guv' file format
    Returns lists of 3D numpy and the origin, deltas and gridcounts.



    First two lines should look something like:
    16.00000000000000         16.00000000000000         16.00000000000000     
           32           32           32
   
    '''
    print("Reading TINKER style guv file")
    print("Warning, origin set to default 0.0, 0.0, 0.0")
    print("Distributions are unlabeled.")
    print("Assuming distributions start on line 8")
    
    distfile = openreadfile(filename)
    filelines = distfile.readlines()
    
    sidelengths = [float(element) for element in filelines[0].split()]
    gridcounts = [int(element) for element in filelines[1].split()]
    deltas = [sidelengths[dim] / gridcounts[dim] for dim in range(3)]
    origin = [0.0]*3
    numspecies = int(filelines[6])
    print("Found %d species." % (numspecies))
    datalines = [line.split() for line in filelines[7:]]
    assert len(datalines) == gridcounts[0]*gridcounts[1]*gridcounts[2]
    dists = []

    for specnum in range(numspecies):
        dists.append(np.reshape(np.array([float(dataline[specnum]) for dataline in datalines]),
                                gridcounts))
    return dists, origin, deltas, gridcounts


def TKRguv2Grids(filename):
    '''
    Reads a TINKER guv into Grid class objects 
    '''

    (distributions, origin, deltas, gridcounts) = readTKRguv(filename)
    grids = []

    for dist in distributions:
        grids.append(Grid(dist, origin, gridcounts, deltas))
    return grids


def h5ToGrids(h5filename):
    '''
Reads distributions from an MDF hdf5 (h5) file into a dictionary
of Grid class objects.
Returns list of dictionaries:

gridlist[speciesindex]["distributiontype"] = GridObject

    '''
    import tables
    h5file = tables.openFile(h5filename)
    try: 
        offset = h5file.root.parameters_uv.coordinates_offset[:]
    except:
        offset = [0,0,0]
        print("No coordinates offset found. Box will be centered on [0.0, 0.0, 0.0]")
    

    #Get universal parameters
    gridcounts = h5file.root.parameters_uv.num_grid_points[:]
    deltas = [h5file.root.parameters_uv.grid_spacing.read()]*3
    origin = [-deltas[dim]*gridcounts[dim]*0.5+offset[dim] for dim in range(3)]

    speciesnames = [species for molecule in 
                    h5file.root.data_vv.parameters_vv.solvent 
                    for species in molecule.names]
    grids = dict((speciesname, {}) for speciesname in speciesnames)
    for key, value in h5file.root._v_children.iteritems():
        # 
        #print("key (in childen) is ",key)
        if "uv" in key[1:] and len(key) < 5:
            # This leave is an array of distributions           
            # key[1:] is to avoid uvv
            # len(key) < 5 is to avoid /parameters_uv/
            print("Reading:", key)
            value = value[:].T
            for specnum, dist in enumerate(value):
                assert len(dist) == gridcounts[0] * gridcounts[1] * gridcounts[2]
                threeddist = np.reshape(dist,[gridcounts[2], gridcounts[1], gridcounts[0]]).T
#                threeddist = np.reshape(dist,gridcounts)
                grids[speciesnames[specnum]][key] = Grid(threeddist, origin, gridcounts, deltas)
    return grids



def printdxfrom1dzfast(
    values,
    origin,
    delta,
    gridcounts,
    filename,
    ):
    ''' Print a dx file'''
    print("Deprecated Function")
    f = open(filename, 'w')
    f.write("#DX file from Dan's program\n")
    f.write('object 1 class gridpositions counts {0} {1} {2}\n'.format(gridcounts[0],
            gridcounts[1], gridcounts[2]))
    f.write('origin {0} {1} {2}\n'.format(origin[0], origin[1],
            origin[2]))
    f.write('delta {0} 0 0\n'.format(delta[0]))
    f.write('delta 0 {0} 0\n'.format(delta[1]))
    f.write('delta 0 0 {0}\n'.format(delta[2]))
    f.write('object 2 class gridconnections counts {0} {1} {2}\n'.format(gridcounts[0],
            gridcounts[1], gridcounts[2]))
    f.write('object 3 class array type double rank 0 items {0} data follows\n'.format(gridcounts[0]
            * gridcounts[1] * gridcounts[2]))
    for value in values:
        f.write('{0}\n'.format(value))
    f.write('object {0} class field\n'.format(filename))
    f.close()


def printdxfrom3d(
    distribution,
    origin,
    delta,
    gridcounts,
    filename,
    ):
    ''' Print a dx file given a 3d list
    This function is used by Grid objects.
    '''
    f = open(filename, 'w')
    f.write("#DX file from Dan Sindhikara's 'grid' program.\n")
    f.write('object 1 class gridpositions counts {0} {1} {2}\n'.format(gridcounts[0],
            gridcounts[1], gridcounts[2]))
    f.write('origin {0} {1} {2}\n'.format(origin[0], origin[1],
            origin[2]))
    f.write('delta {0} 0 0\n'.format(delta[0]))
    f.write('delta 0 {0} 0\n'.format(delta[1]))
    f.write('delta 0 0 {0}\n'.format(delta[2]))
    f.write('object 2 class gridconnections counts {0} {1} {2}\n'.format(gridcounts[0],
            gridcounts[1], gridcounts[2]))
    f.write('object 3 class array type double rank 0 items {0} data follows\n'.format(gridcounts[0]
            * gridcounts[1] * gridcounts[2]))
    for i in range(gridcounts[0]):
        for j in range(gridcounts[1]):
            for k in range(gridcounts[2]):
                f.write('{0}\n'.format(distribution[i][j][k]))
    f.write('object {0} class field\n'.format(filename))
    f.close()


def printgrd(
    distribution,
    origin,
    deltas,
    gridcounts,
    filename,
    ):
    ''' Print a .grd file readable by Accelrys DS given a 3d list.
    Note: This is not UHBD .grd format!!
    This function is used by Grid objects.
    '''
    boxdims = [(gridcounts[dim]-1)*deltas[dim] for dim in range(3)]
    f = open(filename, 'w')
    f.write("Volumetric data written by grid.py.\n")
    f.write('(1p,e12.5)\n')
    f.write('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n' % (boxdims[0],
                                    boxdims[1], boxdims[2],
                                    90.0, 90.0, 90.0))
    f.write('%5d %5d %5d\n' % (gridcounts[0]-1, gridcounts[1]-1, gridcounts[2]-1))
    
    #pleft and pright are number of grid points left or right of 0,0,0
    pleft = [ origin[dim] / deltas[dim] for dim in range(3)]
    pright = [gridcounts[dim] + pleft[dim] - 2 for dim in range(3)]

    # First number (1 or 3) represents fast variable (x or z)
    f.write('    1 %6d %6d %6d %6d %6d %6d\n' % (pleft[0], pright[0], pleft[1], pright[1], pleft[2], pright[2]))
    for z in range(gridcounts[2]):
        for y in range(gridcounts[1]):
            for x in range(gridcounts[0]):
                f.write('%12.5E\n' % (distribution[x][y][z]))
    f.close()

def getcoordfromindices(indices, origin, deltas):
    '''
Returns coordinates as a length 3 list of floats.
# Optimization on June 10, 2013
    '''
    return [float(indices[i]) * deltas[i] + origin[i] for i in range(3)]


def getindicesfromcoord(coord, origin, deltas):
    indices = []
    for i in range(3):
        indices.append(int((coord[i] - origin[i]) / deltas[i] + 0.5))
    return indices

def precomputeshellindices(maxindex):
    '''return a list of 3d lists containing the indices in successive search shells
i.e. 0 = [0,0,0]
     1 = [[-1,-1,-1],[-1,0,-1],..
    essentially how much to shift i,j,k indices from center to find grid point on shell
at index radius
    This will make evacuation phase faster
    '''

    from math import sqrt
    shellindices = [[[0, 0, 0]]]
    for index in range(1, maxindex):

        indicesinthisshell = []
        for i in range(-index, index + 1):
            for j in range(-index, index + 1):
                for k in range(-index, index + 1):
                    if int(sqrt(i * i + j * j + k * k)) == index:  # I think this will miss some
                        indicesinthisshell.append((i, j, k))
        shellindices.append(tuple(indicesinthisshell))
    return tuple(shellindices)


def createprecomputedindicesjson(numshells=40):
    '''
stores a local file called shells.json
    '''

    shellindices = precomputeshellindices(numshells)
    from json import dump
    import os
    outfile = os.path.join(os.path.dirname(__file__), "data_grid", 'shells.json')
    f = open(outfile, 'w')
    dump(shellindices, f)
    f.close()


def readshellindices():
    import os
    infile = os.path.join(os.path.dirname(__file__),"data_grid", 'shells.json')
    from json import load
    f = open(infile, 'rb')
    shellindices = load(f)
    return shellindices


def getlinearweightsandindices(origin, deltas, coord):
    '''
Given one 3D coordinate, return 8 corner indices and weights
that would allow direct linear interpolation
This subroutine is separated from linearinterpolatevalue to allow
precomputation of weights and indices.

This is using the formula for Trilinear interpolation
Significantly optimized (given typical Python constraints).
    '''

    
    x0i = int((coord[0] - origin[0]) / deltas[0])
    y0i = int((coord[1] - origin[1]) / deltas[1])
    z0i = int((coord[2] - origin[2]) / deltas[2])

    
    x0 = x0i * deltas[0] + origin[0]
    y0 = y0i * deltas[1] + origin[1]
    z0 = z0i * deltas[2] + origin[2]

    x1 = x0 + 1
    y1 = y0 + 1
    z1 = z0 + 1

    #Check if xyz1 is beyond array boundary:
    #x1[np.where(x1==input_array.shape[0])] = x0.max()
    #y1[np.where(y1==input_array.shape[1])] = y0.max()
    #z1[np.where(z1==input_array.shape[2])] = z0.max()

    dx = coord[0] - x0
    dy = coord[1] - y0
    dz = coord[2] - z0
    
    normweights = ((1-dx)*(1-dy)*(1-dz),\
                   dx*(1-dy)*(1-dz),\
                   (1-dx)*dy*(1-dz),\
                   (1-dx)*(1-dy)*dz,\
                   dx*(1-dy)*dz,\
                   (1-dx)*dy*dz,\
                   dx*dy*(1-dz),\
                   dx*dy*dz)

    cindices = ((x0i,y0i,z0i),
             (x0i+1,y0i,z0i),
             (x0i,y0i+1,z0i),
             (x0i,y0i,z0i+1),
             (x0i+1,y0i,z0i+1),
             (x0i,y0i+1,z0i+1),
             (x0i+1,y0i+1,z0i),
             (x0i+1,y0i+1,z0i+1))
    return (normweights, cindices)


def linearinterpolatevalue(
    distribution,
    origin,
    deltas,
    coord,
    ):
    '''given a 3d coordinate, using a linear interpolation from the 8 nearest gridpoints,
estimate the value at that coordinate
    '''

    (weights, cindices) = getlinearweightsandindices(origin, deltas,
            coord)
    value = 0
    for (i, mycindices) in enumerate(cindices):
        try:
            value += distribution[mycindices] * weights[i]
        except:
            print('Failed to find gridpoint at', mycindices)
            print('coordinate=', coord)
            return False
    return value


def calcrdf(
    distribution,
    origin,
    deltas,
    coord,
    deltar=0.1,
    maxradius=20.0,
    numsumgrids=20,
    ):
    '''
Calculates the radial distribution function about a point using the 3d distribution
    '''

    def shellintegral(radius, delta, center):
        sum = 0.0
        count = 0.0
        for i in range(int(2.0 * radius / delta)):
            for j in range(int(2.0 * radius / delta)):
                for k in range(int(2.0 * radius / delta)):
                    x = float(i) * delta - radius
                    y = float(j) * delta - radius
                    z = float(k) * delta - radius
                    thisrad = (x ** 2 + y ** 2 + z ** 2) ** 0.5
                    if thisrad > radius - delta / 2.0 and thisrad \
                        < radius + delta / 2.0:
                        mycoord = [center[0] + x, center[1] + y,
                                   center[2] + z]
                        sum += linearinterpolatevalue(distribution,
                                origin, deltas, mycoord)
                        count += 1.0
        sum = sum / count
        return sum

    radii = [(float(i) + 0.5) * deltar for i in range(int(maxradius
             / deltar))]
    gr = []
    for (i, rad) in enumerate(radii):
        subdelta = deltar * (float(i) + 0.5) / float(numsumgrids)
        gr.append(shellintegral(rad, subdelta, coord))

    return (radii, gr)



