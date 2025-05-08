import os
import vtk
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from utilities import ReadVTUFile, ThresholdInBetween, ExtractSurface
from NormalizeMBFMap import MBFNormalization

class PrePostMBFMap(MBFNormalization):
    def __init__(self, args):
        args.InputMBFMap = f"{args.InputFolder}/{args.InputMBF}"
        args.ArrayName = 0
        self.args = args
        super().__init__(args)
        self.MBF_A = ReadVTUFile(f"{args.InputFolder}/{args.InputMBF}")
        self.MBF_B = ReadVTUFile(f"{args.InputFolder[:-1]}B/{args.InputMBF}")
        self.InputLabels = f"{args.InputFolder[:-1]}B/{args.InputLabels}"

    def ReadTerritoryMBF(self, MBFMap, MBF_Labels, ArrayName = 0):
        MBF_data = {}
        Territories = {}
        for key in MBF_Labels.keys():
            MBF_data[key] = np.array([])
            AppendTerritory = vtk.vtkAppendFilter()
            for i in MBF_Labels[key]:
                territory_ = ThresholdInBetween(MBFMap, "TerritoryMaps", i, i)
                AppendTerritory.AddInputData(territory_)
                AppendTerritory.Update()
                MBF_ = vtk_to_numpy(territory_.GetPointData().GetArray(ArrayName))
                MBF_data[key] = np.append(MBF_, MBF_data[key])
            Territories[key] = AppendTerritory.GetOutput()
        
        return MBF_data, Territories


    def PlotBox(self, MBF_Labels, MBF_data_pre, MBF_data_post, ylabel = "MBF (ml/min/100g)"):
        color_list = ['aquamarine', 'sandybrown', 'palegreen', 'lightcyan', 'thistle', 'lavender', 'salmon', 'peachpuff']
        
        Labels = []
        IndexMBF = []
        colors = []
        for (i, key) in enumerate(MBF_Labels.keys()):
            if len(MBF_data_pre[key]) > 0:
                Labels.append(f'{key}_pre')
                Labels.append(f'{key}_post')
                IndexMBF.append(MBF_data_pre[key])
                IndexMBF.append(MBF_data_post[key])
                colors.append(color_list[i])
                colors.append(color_list[i])
            
        _, ax = plt.subplots()
        ax.set_ylabel(ylabel, fontdict={'size':20})
        bplot = ax.boxplot(IndexMBF, patch_artist=True, labels=Labels, showfliers= False)

        for patch, color in zip(bplot['boxes'], colors):
            patch.set_color(color)

        for median in bplot['medians']:
            median.set_color('darkblue')
            median.set_linewidth(2)

        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize = 13)
        plt.tight_layout()
        plt.show()

    def ReadMBFLabels(self):
        MBF_Labels = {"LAD": [], "LCx":[], "Intermedius":[], "Diag1":[], "Diag2":[], "PDA":[], "PL":[]}

        with open(self.InputLabels, "r") as ifile:
            for LINE in ifile:
                line = LINE.split()
                for key in MBF_Labels.keys():
                    if line[1].find(key)>=0: MBF_Labels[key].append(int(line[0]))

        return MBF_Labels

    def main(self):
        MBF_Labels = self.ReadMBFLabels()
        MBF_data_A, Territories_A = self.ReadTerritoryMBF(self.MBF_A, MBF_Labels)
        MBF_data_B, Territories_B = self.ReadTerritoryMBF(self.MBF_B, MBF_Labels)

        self.PlotBox(MBF_Labels, MBF_data_A, MBF_data_B)

        super().Normalize()
        self.args.InputMBFMap = f"{self.args.InputFolder[:-1]}B/{self.args.InputMBF}"
        super().__init__(self.args)
        super().Normalize()

        IndexMBF_A = ReadVTUFile(f"{args.InputFolder}/{os.path.splitext(args.InputMBF)[0]}_Normalized.vtu")
        IndexMBF_B = ReadVTUFile(f"{args.InputFolder[:-1]}B/{os.path.splitext(args.InputMBF)[0]}_Normalized.vtu")
        MBF_data_A, _ = self.ReadTerritoryMBF(IndexMBF_A, MBF_Labels, "IndexMBF")
        MBF_data_B, _ = self.ReadTerritoryMBF(IndexMBF_B, MBF_Labels, "IndexMBF")

        self.PlotBox(MBF_Labels, MBF_data_A, MBF_data_B, "IndexMBF")

        self.BarPlot(self.ProcessVolumeData(Territories_A, Territories_B))

    def ComputeVolume(self, ClosedSurface):
        tri_filter = vtk.vtkTriangleFilter()
        tri_filter.SetInputData(ClosedSurface)
        tri_filter.Update()

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(tri_filter.GetOutput())
        cleaner.Update()

        Mass = vtk.vtkMassProperties()
        Mass.SetInputData(cleaner.GetOutput())
        Mass.Update()

        return Mass.GetVolume()

    def ComputeTerritoryVolume(self, Territory):
        Volume_data = {}
        for (key, item) in Territory.items():
            if item.GetNumberOfPoints() > 0:
                Volume_data[key] = self.ComputeVolume(ExtractSurface(item))

        return Volume_data

    def ProcessVolumeData(self, Territories_A, Territories_B):
        Volume_MBF_A = self.ComputeVolume(ExtractSurface(self.MBF_A))
        Volume_MBF_B = self.ComputeVolume(ExtractSurface(self.MBF_B))
        Volume_data_A = self.ComputeTerritoryVolume(Territories_A)
        Volume_data_B = self.ComputeTerritoryVolume(Territories_B)

        data = {"Territory": [], "Time": [], "Value": []}
        data["Territory"].extend(["Myocardium", "Myocardium"])
        data["Time"].extend(["PreCABG", "PostCABG"])
        data["Value"].extend([Volume_MBF_A, Volume_MBF_B])

        for key in Volume_data_A.keys():
            data["Territory"].extend([key, key])
            data["Time"].extend(["PreCABG", "PostCABG"])
            data["Value"].extend([Volume_data_A[key], Volume_data_B[key]])

        return data

    def BarPlot(self, data):
        df = pd.DataFrame(data)
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Territory", y="Value", hue="Time", palette="pastel")

        plt.ylabel("Volume (mL)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type=str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-InputLabels", "--InputLabels", dest= "InputLabels", type= str, required= False, default= "MBF_Territories_Labels.dat")
    args = parser.parse_args()

    PrePostMBFMap(args).main()