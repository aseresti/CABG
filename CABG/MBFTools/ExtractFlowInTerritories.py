import os
import vtk
import argparse
import numpy as np
from utilities import ReadVTUFile, ThresholdInBetween
from vtk.util.numpy_support import vtk_to_numpy

class ExtractSubtendedFlow():
    def __init__(self, args):
        super().__init__()
        self.args = args

    def ReadMBFFiles(self):
        self.MBF = ReadVTUFile(args.InputMBF)
        labels = f"{os.path.splitext(self.args.InputMBF)[0]}_Labels.dat"
        self.Labels = {}
        with open(labels, "r") as ifile:
            for LINE in ifile:
                line = LINE.split()
                self.Labels[line[1]] = line[0]


    def CalculateVoxelFlow(self, Array, voxel, Unit):
        id_list = voxel.GetPointIds()
        cell_bounds = voxel.GetBounds()
        cell_volume = abs(cell_bounds[0] - cell_bounds[1]) * abs(cell_bounds[2] - cell_bounds[3]) * abs(cell_bounds[4] - cell_bounds[5])

        average_cell_mbf = 0
        for i in range(id_list.GetNumberOfIds()):
            average_cell_mbf += Array.GetValue(id_list.GetId(i))
        
        rho = 1.05
        if Unit == 'mm':
            return rho*average_cell_mbf/id_list.GetNumberOfIds()*cell_volume/1000/100
        elif Unit == 'cm':
            return rho*average_cell_mbf/id_list.GetNumberOfIds()*cell_volume/100

    def CalculateFlowInVoluem(self, Volume, Unit, ArrayName):
        MBFScalarArray = Volume.GetPointData().GetArray(ArrayName)
        NCells = Volume.GetNumberOfCells()
        Flow = 0
        for i in range(NCells):
            voxel = Volume.GetCell(i)
            Flow += self.CalculateVoxelFlow(MBFScalarArray, voxel, Unit)
        
        return Flow
    
    def ExtractSubtendedTerritory(self, TerritoryTag):
        self.ReadMBFFiles()
        SubtendedFlow = 0
        self.TerritoryTags = ""
        for (key, item) in self.Labels.items():
            if TerritoryTag in key:
                self.TerritoryTags += os.path.splitext(key)[0] + "+"
                territory_ = ThresholdInBetween(self.MBF, "TerritoryMaps", int(item), int(item))
                SubtendedFlow += self.CalculateFlowInVoluem(territory_, self.args.Unit)
        
        return SubtendedFlow
    
    def main(self):
        SubtendedFlow = self.ExtractSubtendedTerritory(self.args.TerritoryTag)
        print("Flow = ", int(SubtendedFlow*100)/100, "mL/min")
        ofile_path = f"./{os.path.splitext(os.path.basename(self.args.InputMBF))[0]}_MBFxVolume_{self.args.TerritoryTag}.dat"
        with open(ofile_path, "w") as ofile:
            ofile.writelines("Territory Tags:\n")
            ofile.writelines(f"{self.TerritoryTags}\n")
            ofile.writelines(f"Territory Flow: {int(SubtendedFlow*100)/100} mL/min")

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
            return np.array(TerritoryFlow)/1000/100, TerritoryVolume/1000, NCells, ImageScalars
        elif self.args.Unit == 'cm':
            return np.array(TerritoryFlow)/100, TerritoryVolume, NCells, ImageScalars

    def ExtractCellDataSubtendedTerritory(self, MBF, ArrayName):
        SubtendedFlow = 0
        self.TerritoryTags = ""
        NCells = 0
        TerritoryVolume = 0
        MBFScalarArray = np.array([])
        for (key, item) in self.Labels.items():
            if self.args.TerritoryTag in key:
                self.TerritoryTags += os.path.splitext(key)[0] + "+"
                territory_ = self.ThresholdInBetween(MBF, "TerritoryMaps", int(item), int(item))
                flow_, volume_, nCells, scalararray = self.CalculateCellDataFlow(territory_, ArrayName)
                SubtendedFlow += np.sum(flow_)
                NCells += nCells
                TerritoryVolume += volume_
                MBFScalarArray = np.append(scalararray, MBFScalarArray)
        
        return SubtendedFlow, TerritoryVolume/NCells*1000, TerritoryVolume, MBFScalarArray

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= True)
    #parser.add_argument("-InputLabel", "--InputLabel", dest= "InputLabel", type= str, required= True)
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = str, required = False, default = "ImageScalars")
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=True, dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="mm", required=False)
    args = parser.parse_args()
    
    ExtractSubtendedFlow(args).main()
