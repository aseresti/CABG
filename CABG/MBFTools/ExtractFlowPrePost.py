import os
import vtk
import argparse
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from utilities import ReadVTUFile, ThresholdInBetween
from ExtractFlowInTerritories import ExtractSubtendedFlow
from NormalizeMBFMap import MBFNormalization

class ExtractFlowPrePost(ExtractSubtendedFlow, MBFNormalization):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    def ReadPrePostFiles(self):
        self.MBF_A = ReadVTUFile(f"{self.args.InputFolder}/{self.args.InputMBF}")
        self.MBF_B = ReadVTUFile(f"{self.args.InputFolder[:-1]}B/{self.args.InputMBF}")
        self.InputLabels = f"{self.args.InputFolder[:-1]}B/{self.args.InputLabels}"

    def ReadMBFLabels(self):
        MBF_Labels = {"post_LAD": [], "post_LCx": [], "post_RCA": [], "NonIschemic": []}

        with open(self.InputLabels, "r") as ifile:
            for LINE in ifile:
                line = LINE.split()
                for key in MBF_Labels.keys()[:-1]:
                    if line[1].find(key)>=0: MBF_Labels[key].append(int(line[0]))
                    else: MBF_Labels["NonIschemic"].append(int(line[0]))

        return MBF_Labels

    def ReadTerritoryMBF(self, MBFMap, MBF_Labels, ArrayName = 0):
        MBF_data = {}
        Territories = {}
        for key in MBF_Labels.keys():
            MBF_data[key] = np.array([])
            AppendTerritory = vtk.vtkAppendFilter()
            for i in MBF_Labels[key]:
                territory_ = ThresholdInBetween(MBFMap, "TerritoryMaps", i, i)
                AppendTerritory.AddInputData(territory_)
                MBF_ = vtk_to_numpy(territory_.GetPointData().GetArray(ArrayName))
                MBF_data[key] = np.append(MBF_, MBF_data[key])
            AppendTerritory.Update()
            Territories[key] = AppendTerritory.GetOutput()
        
        return MBF_data, Territories

    def CollectFlowData(self, Territories):
        pass

    def CollectPrePostData(self):
        pass

    def BarPlot(self, BarData):
        pass

    def BoxPlot(self, BoxData):
        pass

    def Statistics(self, MBF_Data):
        pass

    def main(self):
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= True)
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = str, required = False, default = "ImageScalars")
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=True, dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="mm", required=False)
    args = parser.parse_args()