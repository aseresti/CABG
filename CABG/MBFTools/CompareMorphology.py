import os
import vtk
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from PrePostComparison import PrePostMBFMap
from utilities import ReadVTPFile, ReadVTUFile, WriteVTPFile

class CompareMorphology(PrePostMBFMap):
    def __init__(self, args):
        super().__init__(args)
        InputFolderA = args.InputFolder
        InputFolderB = f"{InputFolderA[:-1]}B"
        
        self.CavityCapped_A = ReadVTPFile(f"{InputFolderA}/Morphology/{args.CavityCapped}")
        self.Endocardium_A = ReadVTPFile(f"{InputFolderA}/Morphology/{args.Endocardium}")
        self.Epicardium_A = ReadVTPFile(f"{InputFolderA}/Morphology/{args.Epicardium}")
        self.MBFTerritories_A = ReadVTUFile(f"{InputFolderA}/{args.MBFTerritories}")
        self.output_surface_A = os.path.join(f"{InputFolderA}/Morphology", os.path.splitext(args.Epicardium)[0] + "_WallThickness.vtp")

        self.CavityCapped_B = ReadVTPFile(f"{InputFolderB}/Morphology/{args.CavityCapped}")
        self.Endocardium_B = ReadVTPFile(f"{InputFolderB}/Morphology/{args.Endocardium}")
        self.Epicardium_B = ReadVTPFile(f"{InputFolderB}/Morphology/{args.Epicardium}")
        self.MBFTerritories_B = ReadVTUFile(f"{InputFolderB}/{args.MBFTerritories}")
        self.output_surface_B = os.path.join(f"{InputFolderB}/Morphology", os.path.splitext(args.Epicardium)[0] + "_WallThickness.vtp")

    def ComputeVolume(self, ClosedSurface):
        MassProp = vtk.vtkMassProperties()
        MassProp.SetInputData(ClosedSurface)
        MassProp.Update()
        
        return MassProp.GetVolume()
    
    def ComputeSurfaceArea(self, Surface):
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
            if len(MBF_Labels[key]) > 0: 
                WallThickness_data[key] = np.array([])
                for i in MBF_Labels[key]:
                    territory_ = self.ThresholdInBetweenPoints(Surface, "TerritoryMaps", i, i+1)
                    wall_thickness = vtk_to_numpy(territory_.GetPointData().GetArray("Distance"))
                    WallThickness_data[key] = np.append(wall_thickness, WallThickness_data[key])

        return WallThickness_data

    def PlotResults(self, data):
        df = pd.DataFrame(data)
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="parameter", y="Value", hue="Time", palette="Set2")

        plt.ylabel("Value (units vary)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()

    def main(self):
        Volume_A = self.ComputeVolume(self.CavityCapped_A)
        Volume_B = self.ComputeVolume(self.CavityCapped_B)

        Surface_A_Endo = self.ComputeSurfaceArea(self.Endocardium_A)
        Surface_B_Endo = self.ComputeSurfaceArea(self.Endocardium_B)

        Surface_A_Epi = self.ComputeSurfaceArea(self.Epicardium_A)
        Surface_B_Epi = self.ComputeSurfaceArea(self.Epicardium_B)

        data = {"parameter":[], "Time":[], "Value":[]}
        print("Cavity of Volume- PreCABG: ", Volume_A)
        print("Cavity of Volume- PostCABG: ", Volume_B)
        data["parameter"].extend(["Cavity Volume", "Cavity Volume"])
        data["Time"].extend(["PreCABG (mL)", "PostCABG (mL)"])
        data["Value"].extend([Volume_A, Volume_B])

        print("EndoCardium Surface- PreCABG: ", Surface_A_Endo)
        print("EndoCardium Surface- PostCABG: ", Surface_B_Endo)
        data["parameter"].extend(["Endo Surface", "Endo Surface"])
        data["Time"].extend(["PreCABG (cm^2)", "PostCABG (cm^2)"])
        data["Value"].extend([Surface_A_Endo, Surface_B_Endo])

        print("EndoCardium Surface- PreCABG: ", Surface_A_Epi)
        print("EndoCardium Surface- PostCABG: ", Surface_B_Epi)
        data["parameter"].extend(["Epi Surface", "Epi Surface"])
        data["Time"].extend(["PreCABG (cm^2)", "PostCABG (cm^2)"])
        data["Value"].extend([Surface_A_Epi, Surface_B_Epi])

        Epicardium_WT_A = self.ComputeWallThickness(self.Endocardium_A, self.Epicardium_A)
        Epicardium_WT_B = self.ComputeWallThickness(self.Endocardium_B, self.Epicardium_B)
        Epicardium_WT_Territory_A = self.ProjectTerritories(Epicardium_WT_A, self.MBFTerritories_A)
        Epicardium_WT_Territory_B = self.ProjectTerritories(Epicardium_WT_B, self.MBFTerritories_B)

        WriteVTPFile(self.output_surface_A, Epicardium_WT_Territory_A)
        WriteVTPFile(self.output_surface_B, Epicardium_WT_Territory_B)

        WallThickness_Data_A = self.ExtractWallThicknessInTerritory(Epicardium_WT_Territory_A)
        WallThickness_Data_B = self.ExtractWallThicknessInTerritory(Epicardium_WT_Territory_B)

        for key in WallThickness_Data_A.keys():
            print(key, np.mean(WallThickness_Data_A[key]), np.mean(WallThickness_Data_B[key]))
            data["parameter"].extend([f"{key}_WallThickness", f"{key}_WallThickness"])
            data["Time"].extend(["PreCABG (0.1mm)", "PostCABG (0.1mm)"])
            data["Value"].extend([np.mean(WallThickness_Data_A[key])*100, np.mean(WallThickness_Data_B[key])*100])

        #Epicardium_WT_Territory_A.GetPointData().GetArray("Distance")
        self.PlotResults(data)


    def ThresholdInBetweenPoints(self, Surface, arrayname, value1, value2):
        Threshold=vtk.vtkThresholdPoints()
        Threshold.SetInputData(Surface)
        Threshold.ThresholdBetween(value1,value2)
        Threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,arrayname)
        Threshold.Update()
        return Threshold.GetOutput()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type= str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-InputLabels", "--InputLabels", dest= "InputLabels", type= str, required= False, default= "MBF_Territories_Labels.dat")
    parser.add_argument("-CavityCapped", default= "CavityCapped.vtp", required= False, type= str, dest= "CavityCapped")
    parser.add_argument("-Endocardium", default= "Endocardium.vtp", required= False, type= str, dest= "Endocardium")
    parser.add_argument("-Epicardium", default= "Epicardium.vtp", required= False, type= str, dest= "Epicardium")
    parser.add_argument("-MBFTerritories", default= "MBF_Territories.vtu", required= False, type= str, dest= "MBFTerritories")

    args = parser.parse_args()
    CompareMorphology(args).main()
