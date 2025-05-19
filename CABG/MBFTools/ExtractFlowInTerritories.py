import os
import vtk
import numpy as np
import argparse
from utilities import ReadVTUFile, ThresholdInBetween

class ExtractSubtendedFlow():
    def __init__(self, args):
        self.MBF = ReadVTUFile(args.InputMBF)
        self.args = args
        labels = f"{os.path.splitext(self.args.InputMBF)[0]}_Labels.dat"
        self.Labels = {}
        with open(labels, "r") as ifile:
            for LINE in ifile:
                line = LINE.split()
                self.Labels[line[1]] = line[0]


    def CalculateVoxelFlow(self, Array, voxel):
        id_list = voxel.GetPointIds()
        cell_bounds = voxel.GetBounds()
        cell_volume = abs(cell_bounds[0] - cell_bounds[1]) * abs(cell_bounds[2] - cell_bounds[3]) * abs(cell_bounds[4] - cell_bounds[5])

        average_cell_mbf = 0
        for i in range(id_list.GetNumberOfIds()):
            average_cell_mbf += Array.GetValue(id_list.GetId(i))
        
        if self.args.Unit == 'mm':
            return average_cell_mbf/id_list.GetNumberOfIds()*cell_volume/1000/100
        elif self.args.Unit == 'cm':
            return average_cell_mbf/id_list.GetNumberOfIds()*cell_volume/100

    def CalculateFlowInVoluem(self, Volume):
        MBFScalarArray = Volume.GetPointData().GetArray(self.args.ArrayName)
        NCells = Volume.GetNumberOfCells()
        Flow = 0
        for i in range(NCells):
            voxel = Volume.GetCell(i)
            Flow += self.CalculateVoxelFlow(MBFScalarArray, voxel)
        
        return Flow
    
    def ExtractSubtendedTerritory(self):
        SubtendedFlow = 0
        self.TerritoryTags = ""
        for (key, item) in self.Labels.items():
            if self.args.TerritoryTag in key:
                self.TerritoryTags += os.path.splitext(key)[0] + "+"
                territory_ = ThresholdInBetween(self.MBF, "TerritoryMaps", int(item), int(item))
                SubtendedFlow += self.CalculateFlowInVoluem(territory_)
        
        return SubtendedFlow
    
    def main(self):
        SubtendedFlow = self.ExtractSubtendedTerritory()
        print("Flow = ", int(SubtendedFlow*100)/100, "mL/min")
        ofile_path = f"./{os.path.splitext(os.path.basename(self.args.InputMBF))[0]}_MBFxVolume_{self.args.TerritoryTag}.dat"
        with open(ofile_path, "w") as ofile:
            ofile.writelines("Territory Tags:\n")
            ofile.writelines(f"{self.TerritoryTags}\n")
            ofile.writelines(f"Territory Flow: {int(SubtendedFlow*100)/100} mL/min")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= True)
    #parser.add_argument("-InputLabel", "--InputLabel", dest= "InputLabel", type= str, required= True)
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = str, required = False, default = "ImageScalars")
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=True, dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="mm", required=False)
    args = parser.parse_args()
    
    ExtractSubtendedFlow(args).main()
