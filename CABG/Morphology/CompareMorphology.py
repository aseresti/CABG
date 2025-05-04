import os
import vtk
import argparse
import matplotlib.pyplot as plt

class CompareMorphology():
    def __init__(self, args):
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
        

    def PlotResults(self):
        pass

    def main(self):
        pass

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
    parser.add_argument("-CavityCapped", default= "CavityCapped.vtp", required= False, type= str, dest= "CavityCapped")
    parser.add_argument("-Endocardium", default= "Endocardium", required= False, type= str, dest= "Endocardium")
    parser.add_argument("-Epicardium", default= "Epicardium", required= False, type= str, dest= "Epicardium")
    parser.add_argument("-MBFTerritories", default= "MBF_Territories.vtu", required= False, type= str, dest= "MBFTerritories")

    args = parser.parse_args()
    CompareMorphology(args).main()
