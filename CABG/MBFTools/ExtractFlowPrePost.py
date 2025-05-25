import os
import vtk
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from utilities import ReadVTUFile
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
        Ischemic_Labels = {}
        for tag in self.args.TerritoryTag:
            Ischemic_Labels[f"{tag}"] =  []
        Ischemic_Labels["NonIschemic"] = []
        keys = list(Ischemic_Labels.keys())[:-1]
        with open(self.InputLabels, "r") as ifile:
            for i, LINE in enumerate(ifile):
                if i == 0: 
                    continue
                line = LINE.split()
                found = False
                for key in keys:
                    if line[1].find(key)>=0: 
                        Ischemic_Labels[key].append(int(line[0]))
                        found = True
                        break
                if not found: 
                    Ischemic_Labels["NonIschemic"].append(int(line[0]))


        #Ischemic_Labels = {k:v for k, v in Ischemic_Labels.items() if len(v)>0}
        
        return Ischemic_Labels

    def ReadTerritoryMBF(self, MBFMap, MBF_Labels, ArrayName):
        MBF_data = {}
        Territories = {}
        for key in MBF_Labels.keys():
            MBF_data[key] = np.array([])
            AppendTerritory = vtk.vtkAppendFilter()
            for i in MBF_Labels[key]:
                territory_ = self.ThresholdInBetween(MBFMap, "TerritoryMaps", i, i)
                AppendTerritory.AddInputData(territory_)
                MBF_ = vtk_to_numpy(territory_.GetCellData().GetArray(ArrayName))
                MBF_data[key] = np.append(MBF_, MBF_data[key])
            AppendTerritory.Update()
            Territories[key] = AppendTerritory.GetOutput()
        
        return MBF_data, Territories

    def CollectFlowData(self, Territories, ArrayName):
        Flow_Territoris = {k:v for k, v in Territories.items()}
        Volume_Territories = {k:v for k, v in Territories.items()}
        Average_Flow = {k:v for k, v in Territories.items()}

        for (key, value) in Territories.items():
            flow_, territory_volume_, NCell = self.CalculateCellDataFlow(value, ArrayName)
            Flow_Territoris[key] = flow_
            Volume_Territories[key] = territory_volume_
            Average_Flow[key] = flow_/NCell
            voxel_size = territory_volume_/NCell
            
        return Flow_Territoris, Volume_Territories, Average_Flow, voxel_size
    
    def ThresholdInBetween(self, Volume, arrayname, value1, value2):
        Threshold=vtk.vtkThreshold()
        Threshold.SetInputData(Volume)
        Threshold.SetLowerThreshold(value1)
        Threshold.SetUpperThreshold(value2)
        Threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,arrayname)
        Threshold.Update()
        return Threshold.GetOutput()
    
    def ConvertPointDataToCellData(self, pointdata):
        PointToCell = vtk.vtkPointDataToCellData()
        PointToCell.SetInputData(pointdata)
        PointToCell.Update()

        return PointToCell.GetOutput()

    def CalculateCellDataFlow(self, Territory, ArrayName):
        CellData = Territory.GetCellData()
        ImageScalars = vtk_to_numpy(CellData.GetArray(ArrayName))
        NCells = Territory.GetNumberOfCells()
        TerritoryVolume = 0
        TerritoryFlow = []
        for i in range(NCells):
            cell = Territory.GetCell(i)
            cell_bounds = cell.GetBounds()
            cell_volume = abs(cell_bounds[0] - cell_bounds[1]) * abs(cell_bounds[2] - cell_bounds[3]) * abs(cell_bounds[4] - cell_bounds[5])
            TerritoryFlow.append(ImageScalars[i]*cell_volume)
            TerritoryVolume += cell_volume

        if self.args.Unit == 'mm':
            return np.array(TerritoryFlow)/1000/100, TerritoryVolume/1000, NCells
        elif self.args.Unit == 'cm':
            return np.array(TerritoryFlow)/100, TerritoryVolume, NCells


    def BarPlot(self, BarData, ylabel):
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

        plt.ylabel(ylabel)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()

    def BoxPlot(self, BoxData):
        pass

    def TerritoryStatistics(self, Territories,ArrayName):
        MBFStatistics = {}
        for (key, item) in Territories.items():
            item_ = vtk_to_numpy(item.GetCellData().GetArray(ArrayName))
            MBFStatistics[key] = {
            "Mean": np.mean(item_),
            "std": np.std(item_),
            "Median": np.median(item_),
            "IQR": np.percentile(item_, 75) - np.percentile(item_, 25)
            }

        return MBFStatistics
    
    def MyocardiumStatistics(self, MBF, ArrayName):
        Array = vtk_to_numpy(MBF.GetCellData().GetArray(ArrayName))
        Myocardium_Statistics  = {
            "Mean": np.mean(Array),
            "std": np.std(Array),
            "Median": np.median(Array),
            "IQR": np.percentile(Array, 75) - np.percentile(Array, 25)
        }
        return Myocardium_Statistics

    
    def Normalize(self, MBF, ArrayName):
        ScalarArray = MBF.GetCellData().GetArray(ArrayName)
        per_75th = np.percentile(vtk_to_numpy(ScalarArray), 75)
        IndexMBFArray = ScalarArray/per_75th
        IndexMBF = numpy_to_vtk(IndexMBFArray)
        IndexMBF.SetName("IndexMBF")
        MBF.GetCellData().AddArray(IndexMBF)

        return per_75th, MBF
    
    def main(self):
        self.ReadPrePostFiles()
        MBFLabels = self.ReadMBFLabels()
        self.MBF_A = self.ConvertPointDataToCellData(self.MBF_A)
        self.MBF_B = self.ConvertPointDataToCellData(self.MBF_B)
        _, Territories_A = self.ReadTerritoryMBF(self.MBF_A, MBFLabels, "ImageScalars")
        _, Territories_B = self.ReadTerritoryMBF(self.MBF_B, MBFLabels, "scalars")
        Flow_A, Volume_A, AverageFlow_A, VoxelSize_A = self.CollectFlowData(Territories_A, "ImageScalars")
        Flow_B, Volume_B, AverageFlow_B, VoxelSize_B = self.CollectFlowData(Territories_B, "scalars")

        Bardata = {"Territory": [], "Time": [], "Value": []}
        for key in Flow_A.keys():
            Bardata["Territory"].extend([key, key])
            Bardata["Time"].extend(["PreCABG", "PostCABG"])
            Bardata["Value"].extend([np.sum(Flow_A[key]), np.sum(Flow_B[key])])

        #self.BarPlot(Bardata)

        Bardata = {"Territory": [], "Time": [], "Value": []}
        for key in Flow_A.keys():
            Bardata["Territory"].extend([key, key])
            Bardata["Time"].extend(["PreCABG", "PostCABG"])
            Bardata["Value"].extend([np.sum(AverageFlow_A[key])*1000, np.sum(AverageFlow_B[key])*1000])

        #self.BarPlot(Bardata)
        
        MBFStat_A = self.TerritoryStatistics(Territories_A, "ImageScalars")
        MBFStat_B = self.TerritoryStatistics(Territories_B, "scalars")

        perc75_A, IndexMBF_A = self.Normalize(self.MBF_A, "ImageScalars")
        perc75_B, IndexMBF_B = self.Normalize(self.MBF_B, "scalars")
        _, ITerritories_A = self.ReadTerritoryMBF(IndexMBF_A, MBFLabels, "IndexMBF")
        _, ITerritories_B = self.ReadTerritoryMBF(IndexMBF_B, MBFLabels, "IndexMBF")

        IndexMBFStat_A = self.TerritoryStatistics(ITerritories_A, "IndexMBF")
        IndexMBFStat_B = self.TerritoryStatistics(ITerritories_B, "IndexMBF")

        IndexFlow_A, _, AverageIndexFlow_A, _ = self.CollectFlowData(ITerritories_A, "IndexMBF")
        IndexFlow_B, _, AverageIndexFlow_B, _ = self.CollectFlowData(ITerritories_B, "IndexMBF")

        Bardata = {"Territory": [], "Time": [], "Value": []}
        for key in IndexFlow_A.keys():
            Bardata["Territory"].extend([key, key])
            Bardata["Time"].extend(["PreCABG", "PostCABG"])
            Bardata["Value"].extend([np.sum(IndexFlow_A[key]), np.sum(IndexFlow_B[key])])

        self.BarPlot(Bardata, "relative Flow (1/min)")

        Bardata = {"Territory": [], "Time": [], "Value": []}
        for key in IndexFlow_A.keys():
            Bardata["Territory"].extend([key, key])
            Bardata["Time"].extend(["PreCABG", "PostCABG"])
            Bardata["Value"].extend([np.sum(AverageIndexFlow_A[key]), np.sum(AverageIndexFlow_B[key])])

        #self.BarPlot(Bardata, "Average relative Flow (\u00b5/min/Voxel)")


        stats_A = self.MyocardiumStatistics(self.MBF_A, "ImageScalars")
        stats_B = self.MyocardiumStatistics(self.MBF_B, "scalars")
        index_stats_A = self.MyocardiumStatistics(IndexMBF_A, 'IndexMBF')
        index_stats_B = self.MyocardiumStatistics(IndexMBF_B, 'IndexMBF')

        Bardata2 = {"Territory": [], "Time": [], "Value": []}
        for key in IndexMBFStat_A.keys():
            Bardata2["Territory"].extend([key, key])
            Bardata2["Time"].extend(["PreCABG", "PostCABG"])
            Bardata2["Value"].extend([IndexMBFStat_A[key]['Mean'], IndexMBFStat_B[key]['Mean']])

        self.BarPlot(Bardata2, "Average Index MBF (1/min/100mL)")

        opath = os.path.join(self.args.InputFolder, "PrePostStatistics.dat")
        with open(opath, 'w') as ofile:
            ofile.writelines("--- MBF (mL/min/100mL) and Index MBF")
            ofile.writelines("Case, Statistics, MBF_A, MBF_B, IndexMBF_A, IndexMBF_B\n")
            ofile.writelines(f"Myocardium, Mean, {stats_A['Mean']}, {stats_B['Mean']}, {index_stats_A['Mean']}, {index_stats_B['Mean']}\n")
            ofile.writelines(f"Myocardium, std, {stats_A['std']}, {stats_B['std']}, {index_stats_A['std']}, {index_stats_B['std']}\n")
            ofile.writelines(f"Myocardium, Median, {stats_A['Median']}, {stats_B['Median']}, {index_stats_A['Median']}, {index_stats_B['Median']}\n")
            ofile.writelines(f"Myocardium, IQR, {stats_A['IQR']}, {stats_B['IQR']}, {index_stats_A['IQR']}, {index_stats_B['IQR']}\n")
            ofile.writelines(f"Myocardium, 75th Percentile, 0, 0, {perc75_A}, {perc75_B}\n")
            ofile.writelines(f"Myocardium, VoxelSize, {VoxelSize_A}, {VoxelSize_B}, _, _\n")
            for key in MBFStat_A.keys():
                ofile.writelines(f"{key}, Mean, {MBFStat_A[key]['Mean']}, {MBFStat_B[key]['Mean']}, {IndexMBFStat_A[key]['Mean']}, {IndexMBFStat_B[key]['Mean']}\n")
                ofile.writelines(f"{key}, std, {MBFStat_A[key]['std']}, {MBFStat_B[key]['std']}, {IndexMBFStat_A[key]['std']}, {IndexMBFStat_B[key]['std']}\n")
                ofile.writelines(f"{key}, Median, {MBFStat_A[key]['Median']}, {MBFStat_B[key]['Median']}, {IndexMBFStat_A[key]['Median']}, {IndexMBFStat_B[key]['Median']}\n")
                ofile.writelines(f"{key}, IQR, {MBFStat_A[key]['IQR']}, {MBFStat_B[key]['IQR']}, {IndexMBFStat_A[key]['IQR']}, {IndexMBFStat_B[key]['IQR']}\n")
                ofile.writelines(f"{key}, Territory Volume (mL), {Volume_A[key]}, {Volume_B[key]}, _, _ \n")
            ofile.writelines("--- Flow (mL/min) and ralative Flow (1/min)")
            for key in Flow_A.keys():
                ofile.writelines(f"{key}, Mean, {np.mean(Flow_A[key])}, {np.mean(Flow_B[key])}, {np.mean(IndexFlow_A[key])}, {np.mean(IndexFlow_B[key])}\n")
                ofile.writelines(f"{key}, std, {np.std(Flow_A[key])}, {np.std(Flow_B[key])}, {np.std(IndexFlow_A[key])}, {np.std(IndexFlow_B[key])}\n")
                ofile.writelines(f"{key}, Median, {np.median(Flow_A[key])}, {np.median(Flow_B[key])}, {np.median(IndexFlow_A[key])}, {np.median(IndexFlow_B[key])}\n")
                ofile.writelines(f"{key}, IQR, {np.percentile(Flow_A[key], 75) - np.percentile(Flow_A[key], 25)}, {np.percentile(Flow_B[key], 75) - np.percentile(Flow_B[key], 25)}, {np.percentile(IndexFlow_A[key], 75) - np.percentile(IndexFlow_A[key], 25)}, {np.percentile(IndexFlow_B[key], 75) - np.percentile(IndexFlow_B[key], 25)}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type=str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-InputLabels", "--InputLabels", dest= "InputLabels", type= str, required= False, default="MBF_Territories_Labels.dat")
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = int, required = False, default = 0)
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=True, nargs= "+", dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="cm", required=False)
    args = parser.parse_args()

    ExtractFlowPrePost(args).main()