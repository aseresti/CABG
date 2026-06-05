import os
import numpy as np
import argparse
from utilities import ReadVTUFile, WriteVTUFile
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class MBFNormalization():
    def __init__(self, args):
        super().__init__()
        self.args = args

    def Normalize(self, MBF):
        ImageScalars = MBF.GetPointData().GetArray(self.args.ArrayName)
        per_75th = np.percentile(vtk_to_numpy(ImageScalars), 75)
        IndexMBFArray = ImageScalars/per_75th
        IndexMBF = numpy_to_vtk(IndexMBFArray)
        IndexMBF.SetName("IndexMBF")
        MBF.GetPointData().AddArray(IndexMBF)

        return per_75th, MBF


    def main(self):
        MBF = ReadVTUFile(self.args.InputMBFMap)
        _, IndexMBF = self.Normalize(MBF)
        OPath = f"{os.path.splitext(self.args.InputMBFMap)[0]}_Normalized.vtu"
        WriteVTUFile(OPath, IndexMBF)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputMBFMap", "--InputMBFMap", dest = "InputMBFMap", type = str, required = True)
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = int, required = False, default= 0)
    args = parser.parse_args()
    
    MBFNormalization(args).main()