import os
import vtk
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from utilities import ReadVTUFile, ThresholdInBetween
from ExtractFlowInTerritories import ExtractSubtendedFlow

class ExtractFlowPrePost(ExtractSubtendedFlow):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    def ReadPrePostFiles(self):
        self.MBF_A = ReadVTUFile(os.path.join(self.args.InputFolder,self.args.InputMBF))
        self.MBF_B = ReadVTUFile(os.path.join(f"{self.args.InputFolder[:-1]}B",self.args.InputMBF))
        self.InputLabels = os.path.join(f"{self.args.InputFolder[:-1]}B",self.args.InputLabels)

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
        SubtendedFlow = {key: super().CalculateFlowInVoluem(item) for (key, item) in Territories.items()}            
        return SubtendedFlow

    def BarPlot(self, BarData):
        df = pd.DataFrame(BarData)
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Territory", y="Value", hue="Time", palette="pastel")

        plt.ylabel("Flow (mL/min)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()

    def BoxPlot(self, BoxData):
        pass

    def Statistics(self, MBF_Data):
        statistics_dict = {"Mean":0, "Stdev": 0, "Median":0, "IQR":0}
        MBFStatistics = {key: statistics_dict for key in MBF_Data.keys()}
        for (key, item) in MBF_Data.items():
            MBFStatistics[key]["Mean"] = np.mean(item)
            MBFStatistics[key]["Stdev"] = np.std(item)
            MBFStatistics[key]["Median"] = np.median(item)
            MBFStatistics[key]["IQR"] = np.percentile(item, 75) - np.percentile(item, 25)

        return MBFStatistics
    
    def Normalize(self, MBF, ArrayName):
        ScalarArray = MBF.GetPointData().GetArray(ArrayName)
        per_75th = np.percentile(vtk_to_numpy(ScalarArray), 75)
        IndexMBFArray = ScalarArray/per_75th
        IndexMBF = numpy_to_vtk(IndexMBFArray)
        IndexMBF.SetName("IndexMBF")
        MBF.GetPointData().AddArray(IndexMBF)

        return per_75th, MBF
    
    def main(self):
        self.ReadPrePostFiles()
        MBFLabels = self.ReadMBFLabels()
        MBFData_A, Territories_A = self.ReadTerritoryMBF(self.MBF_A, MBFLabels)
        MBFData_B, Territories_B = self.ReadTerritoryMBF(self.MBF_B, MBFLabels)
        Flow_A = self.CollectFlowData(Territories_A)
        Flow_B = self.CollectFlowData(Territories_B)

        Bardata = {"Territory": [], "Time": [], "Value": []}
        for key in Flow_A.keys():
            Bardata["Territory"].extend([key, key])
            Bardata["Time"].extend(["PreCABG", "PostCABG"])
            Bardata["Value"].extend([Flow_A[key], Flow_B[key]])
        
        MBFStat_A = self.Statistics(MBFData_A)
        MBFStat_B = self.Statistics(MBFData_B)
        perc75_A, IndexMBF_A = self.Normalize(MBFData_A, 0)
        perc75_B, IndexMBF_B = self.Normalize(MBFData_B, 0)
        IndexMBFData_A, _ = self.ReadTerritoryMBF(IndexMBF_A, MBFLabels)
        IndexMBFData_B, _ = self.ReadTerritoryMBF(IndexMBF_B, MBFLabels)
        IndexMBFStat_A = self.Statistics(IndexMBFData_A)
        IndexMBFStat_B = self.Statistics(IndexMBFData_B)
        #todo: Add Output .dat file

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type=str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = int, required = False, default = 0)
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=False, dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="cm", required=False)
    args = parser.parse_args()