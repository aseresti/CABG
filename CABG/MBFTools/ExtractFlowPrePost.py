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
        Ischemic_Labels = {"post_LAD": [], "post_LCx":[], "post_PL":[], "NonIschemic": []}
        keys = list(Ischemic_Labels.keys())[:-1]
        with open(self.InputLabels, "r") as ifile:
            for i, LINE in enumerate(ifile):
                if i == 0: continue
                line = LINE.split()
                if line[1].find(keys[0])>=0: Ischemic_Labels[keys[0]].append(int(line[0]))
                elif line[1].find(keys[1])>=0: Ischemic_Labels[keys[1]].append(int(line[0]))
                elif line[1].find(keys[2])>=0: Ischemic_Labels[keys[2]].append(int(line[0]))
                else: Ischemic_Labels["NonIschemic"].append(int(line[0]))


        Ischemic_Labels = {k:v for k, v in Ischemic_Labels.items() if len(v)>0}
        
        return Ischemic_Labels

    def ReadTerritoryMBF(self, MBFMap, MBF_Labels, ArrayName):
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

    def CollectFlowData(self, Territories, Unit, ArrayName):
        SubtendedFlow = {key: ExtractSubtendedFlow(self.args).CalculateFlowInVoluem(item, Unit, ArrayName) for (key, item) in Territories.items()}            
        return SubtendedFlow

    def BarPlot(self, BarData):
        df = pd.DataFrame(BarData)
        plt.figure(figsize=(8, 5))
        pastel_colors = sns.color_palette("pastel")
        selected_colors = [pastel_colors[3], pastel_colors[2]]
        ax = sns.barplot(data=df, x="Territory", y="Value", hue="Time", palette=selected_colors)
        for p in ax.patches:
            height = p.get_height()
            if not pd.isna(height):  # Check for missing values
                ax.text(
                    p.get_x() + p.get_width() / 2,
                    height,
                    f"{height:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=9
                )

        plt.ylabel("Flow (mL/min)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()

    def BoxPlot(self, BoxData):
        pass

    def TerritoryStatistics(self, MBF_Data):
        statistics_dict = {"Mean":0, "Stdev": 0, "Median":0, "IQR":0}
        MBFStatistics = {key: statistics_dict for key in MBF_Data.keys()}
        for (key, item) in MBF_Data.items():
            item = np.array(item)
            MBFStatistics[key]["Mean"] = np.mean(item)
            MBFStatistics[key]["std"] = np.std(item)
            MBFStatistics[key]["Median"] = np.median(item)
            MBFStatistics[key]["IQR"] = np.percentile(item, 75) - np.percentile(item, 25)

        return MBFStatistics
    
    def MyocardiumStatistics(self, MBF, ArrayName):
        Array = vtk_to_numpy(MBF.GetPointData().GetArray(ArrayName))
        Myocardium_Statistics  = {
            "Mean": np.mean(Array),
            "std": np.std(Array),
            "Median": np.median(Array),
            "IQR": np.percentile(Array, 75) - np.percentile(Array, 25)
        }
        return Myocardium_Statistics

    
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
        MBFData_A, Territories_A = self.ReadTerritoryMBF(self.MBF_A, MBFLabels, "ImageScalars")
        MBFData_B, Territories_B = self.ReadTerritoryMBF(self.MBF_B, MBFLabels, "scalars")
        Flow_A = self.CollectFlowData(Territories_A, self.args.Unit, "ImageScalars")
        Flow_B = self.CollectFlowData(Territories_B, self.args.Unit, "scalars")

        Bardata = {"Territory": [], "Time": [], "Value": []}
        for key in Flow_A.keys():
            Bardata["Territory"].extend([key, key])
            Bardata["Time"].extend(["PreCABG", "PostCABG"])
            Bardata["Value"].extend([Flow_A[key], Flow_B[key]])

        self.BarPlot(Bardata)
        
        MBFStat_A = self.TerritoryStatistics(MBFData_A)
        MBFStat_B = self.TerritoryStatistics(MBFData_B)
        perc75_A, IndexMBF_A = self.Normalize(self.MBF_A, "ImageScalars")
        perc75_B, IndexMBF_B = self.Normalize(self.MBF_B, "scalars")
        IndexMBFData_A, _ = self.ReadTerritoryMBF(IndexMBF_A, MBFLabels, "IndexMBF")
        IndexMBFData_B, _ = self.ReadTerritoryMBF(IndexMBF_B, MBFLabels, "IndexMBF")

        IndexMBFStat_A = self.TerritoryStatistics(IndexMBFData_A)
        IndexMBFStat_B = self.TerritoryStatistics(IndexMBFData_B)



        stats_A = self.MyocardiumStatistics(self.MBF_A, "ImageScalars")
        stats_B = self.MyocardiumStatistics(self.MBF_B, "scalars")
        index_stats_A = self.MyocardiumStatistics(IndexMBF_A, 'IndexMBF')
        index_stats_B = self.MyocardiumStatistics(IndexMBF_B, 'IndexMBF')
        opath = os.path.join(self.args.InputFolder, "PrePostStatistics.dat")
        with open(opath, 'w') as ofile:
            ofile.writelines("Case, Statistics, MBF_A, MBF_B, IndexMBF_A, IndexMBF_B\n")
            ofile.writelines(f"Myocardium, Mean, {stats_A['Mean']}, {stats_B['Mean']}, {index_stats_A['Mean']}, {index_stats_B['Mean']}\n")
            ofile.writelines(f"Myocardium, std, {stats_A['std']}, {stats_B['std']}, {index_stats_A['std']}, {index_stats_B['std']}\n")
            ofile.writelines(f"Myocardium, Median, {stats_A['Median']}, {stats_B['Median']}, {index_stats_A['Median']}, {index_stats_B['Median']}\n")
            ofile.writelines(f"Myocardium, IQR, {stats_A['IQR']}, {stats_B['IQR']}, {index_stats_A['IQR']}, {index_stats_B['IQR']}\n")
            ofile.writelines(f"Myocardium, 75th Percentile, 0, 0, {perc75_A}, {perc75_B}\n")
            for key in MBFStat_A.keys():
                ofile.writelines(f"{key}, Mean, {MBFStat_A[key]['Mean']}, {MBFStat_B[key]['Mean']}, {IndexMBFStat_A[key]['Mean']}, {IndexMBFStat_B[key]['Mean']}\n")
                ofile.writelines(f"{key}, std, {MBFStat_A[key]['std']}, {MBFStat_B[key]['std']}, {IndexMBFStat_A[key]['std']}, {IndexMBFStat_B[key]['std']}\n")
                ofile.writelines(f"{key}, Median, {MBFStat_A[key]['Median']}, {MBFStat_B[key]['Median']}, {IndexMBFStat_A[key]['Median']}, {IndexMBFStat_B[key]['Median']}\n")
                ofile.writelines(f"{key}, IQR, {MBFStat_A[key]['IQR']}, {MBFStat_B[key]['IQR']}, {IndexMBFStat_A[key]['IQR']}, {IndexMBFStat_B[key]['IQR']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type=str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-InputLabels", "--InputLabels", dest= "InputLabels", type= str, required= False, default="MBF_Territories_Labels.dat")
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = int, required = False, default = 0)
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=False, dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="cm", required=False)
    args = parser.parse_args()

    ExtractFlowPrePost(args).main()