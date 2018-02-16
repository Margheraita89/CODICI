#!/usr/bin/env python

import sys
import math
import vtk
import numpy as np
from vmtk import pypes
from vmtk import vmtkscripts
from vmtk import vtkvmtk
from vmtk import vmtkrenderer

computeareacells = 'ComputeAreaCells'

class ComputeAreaCells(pypes.pypeScript):
      
    def __init__(self):
        pypes.pypeScript.__init__(self)
        self.Zone = None

        self.SetScriptName('computeareacells')
        self.SetScriptDoc('')
        self.SetInputMembers([
            ['Zone0', 'zone0','vtkPolyData',1,'','the considered zone 0','vmtksurfacereader'],
            ['Zone1', 'zone1','vtkPolyData',1,'','the considered zone 1','vmtksurfacereader'],
            ['Zone2', 'zone2','vtkPolyData',1,'','the considered zone 2','vmtksurfacereader'],
            ['Zone3', 'zone3','vtkPolyData',1,'','the considered zone 3','vmtksurfacereader'],
            ])

        self.SetOutputMembers([
            ])

    def compute_area(self,zone):
        num_cells_zone = zone.GetNumberOfCells()
        print "Number of Cells = ", num_cells_zone
        
        num_tup = zone.GetPointData().GetArray("Principal Stress").GetNumberOfTuples()
        print "Number of Tuples = ", num_tup
        # for i in range(0,num_tup):
        #     # Ricavo il von Mises in generale
        #     # print i
        #     vM = zone.GetPointData().GetArray("Principal Stress").GetComponent(i,5)

        # Creo la matrice dove andro' ad immagazzinare tutti gli id delle celle nella zona da considerare
        mat_cells = []
        vM_mean = []
        somma = 0
        area_tot = 0
        vM_tot = 0

        ## Ricavo il von Mises per ogni cella
        for j in range(num_cells_zone):
            # Considero la j-esima cella
            cell = zone.GetCell(j)
            pointIds = cell.GetPointIds()
            area = cell.ComputeArea()
            row = []
            vM = []
            for i in range(pointIds.GetNumberOfIds()):
                # print pointIds.GetId(i)
                # print zone.GetPointData().GetArray("Principal Stress").GetComponent(pointIds.GetId(i),5)
                row.append(pointIds.GetId(i))
                vM.append(zone.GetPointData().GetArray("Principal Stress").GetComponent(pointIds.GetId(i),5))
            
            #print "Numero cella, Indici punti cella = ", j, row
            # print vM
            # print np.mean(vM)
            mat_cells.append(row)
            vM_mean.append(np.mean(vM))
            somma = somma + np.multiply(np.mean(vM),area)
            area_tot = area_tot + area

        print " "
        print "Valore della sommatoria = ", somma
        print "Valore dell'area totale = ", area_tot
        val_finale = np.divide(somma,area_tot)
        print "Valore finale dell'integrale per zona = ", val_finale

        return val_finale

    def Execute (self):

        # -----------------------------------------------------
        # load mesh and zone to be considered
        # -----------------------------------------------------
        zone0 = self.Zone0
        zone1 = self.Zone1
        zone2 = self.Zone2
        zone3 = self.Zone3

        vM_Z0 = compute_area(zone0)
        vM_Z1 = compute_area(zone1)
        vM_Z2 = compute_area(zone2)
        vM_Z3 = compute_area(zone3)

if __name__ == '__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()