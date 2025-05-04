import os
import vtk
import argparse
import numpy as np
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from MBFTools.PrePostComparison import PrePostMBFMap

class CompareMorphology(PrePostMBFMap):
    def __init__(self, args):
        super().__init__(args)
        InputFolderA = args.InputFolder
        InputFolderB = f"{InputFolderA[:-1]}B"
        
        self.CavityCapped_A = self.ReadVTPFile(f"{InputFolderA}/Morphology/{args.CavityCapped}")
        self.Endocardium_A = self.ReadVTPFile(f"{InputFolderA}/Morphology/{args.Endocardium}")
        self.Epicardium_A = self.ReadVTPFile(f"{InputFolderA}/Morphology/{args.Epicardium}")
        self.MBFTerritories_A = self.ReadVTUFile(f"{InputFolderA}/{args.MBFTerritories}")

        self.CavityCapped_B = self.ReadVTPFile(f"{InputFolderB}/Morphology/{args.CavityCapped}")
        self.Endocardium_B = self.ReadVTPFile(f"{InputFolderB}/Morphology/{args.Endocardium}")
        self.Epicardium_B = self.ReadVTPFile(f"{InputFolderB}/Morphology/{args.Epicardium}")
        self.MBFTerritories_B = self.ReadVTUFile(f"{InputFolderB}/{args.MBFTerritories}")

    def ComputeVolume(self, ClosedSurface):
        MassProp = vtk.vtkMassProperties()
        MassProp.SetInputData(ClosedSurface)
        MassProp.Update()
        
        return MassProp.GetVolume()
    
    def ComputeSurface(self, Surface):
        MassProp = vtk.vtkMassProperties()
        MassProp.SetInputData(Surface)
        MassProp.Update()

        return MassProp.GetSurfaceArea()

    def ComputeWallThickness(self, Endocardium, Epicardium):
        distancefilter = vtk.vtkDistancePolyDataFilter()
        distancefilter.SetInputData(1, Endocardium)
        distancefilter.SetInputData(0, Epicardium)
        distancefilter.SignedDistanceOff()
        distancefilter.Update()

        return distancefilter.GetOutput()
    
    def ProjectTerritories(self, Surface, Volume):
        TerritoyProfile_Array = vtk.vtkFloatArray()
        TerritoyProfile_Array.SetName("TerritoryMaps")
        TerritoyProfile_Array.SetNumberOfComponents(1)
        TerritoyProfile_Array.SetNumberOfTuples(Surface.GetNumberOfPoints())

        TerritoryMap = Volume.GetPointData().GetArray("TerritoryMaps")

        Locator = vtk.vtkPointLocator()
        Locator.SetDataSet(Volume)
        Locator.BuildLocator()

        for i in range(Surface.GetNumberOfPoints()):
            point = Surface.GetPoint(i)
            closest_point_id = Locator.FindClosestPoint(point)

            TerritoyProfile_Array.SetValue(i, TerritoryMap.GetValue(closest_point_id))

        Surface.GetPointData().AddArray(TerritoyProfile_Array)

        return Surface

    def ExtractWallThicknessInTerritory(self, Surface):
        MBF_Labels = super().ReadMBFLabels()
        WallThickness_data = dict()
        for key in MBF_Labels.keys():
            WallThickness_data[key] = np.array([])
            for i in MBF_Labels[key]:
                territory_ = self.ThresholdInBetweenPoints(Surface, "TerritoryMaps", i, i+1)
                wall_thickness = vtk_to_numpy(territory_.GetPointData().GetArray("TerritoryMaps"))
                WallThickness_data[key] = np.append(wall_thickness, WallThickness_data[key])

        return WallThickness_data

    def PlotResults(self):
        pass

    def main(self):
        pass

    def ThresholdInBetweenPoints(self, Surface, arrayname, value1, value2):
        Threshold=vtk.vtkThresholdPoints()
        Threshold.SetInputData(Surface)
        Threshold.ThresholdBetween(value1,value2)
        Threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,arrayname)
        Threshold.Update()
        return Threshold.GetOutput()

    def ReadVTPFile(self, FileName):
        reader=vtk.vtkXMLPolyDataReader()
        reader.SetFileName(FileName)
        reader.Update()
        
        return reader.GetOutput()
    
    def WriteVTPFile(self, FileName, Data):
        writer=vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(FileName)
        writer.SetInputData(Data)
        writer.Update()

    def ReadVTUFile(self, FileName):
        reader=vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(FileName)
        reader.Update()
        
        return reader.GetOutput()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type= str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-InputLabels", "--InputLabels", dest= "InputLabels", type= str, required= False, default= "MBF_Territories_Labels.dat")
    parser.add_argument("-CavityCapped", default= "CavityCapped.vtp", required= False, type= str, dest= "CavityCapped")
    parser.add_argument("-Endocardium", default= "Endocardium", required= False, type= str, dest= "Endocardium")
    parser.add_argument("-Epicardium", default= "Epicardium", required= False, type= str, dest= "Epicardium")
    parser.add_argument("-MBFTerritories", default= "MBF_Territories.vtu", required= False, type= str, dest= "MBFTerritories")

    args = parser.parse_args()
    CompareMorphology(args).main()
