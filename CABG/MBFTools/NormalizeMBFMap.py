import os
import numpy as np
import argparse
from utilities import ReadVTUFile, WriteVTUFile
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class MBFNormalization():
    def __init__(self, args):
        self.args = args
        self.MBF = ReadVTUFile(self.args.InputMBFMap)

    def Normalize(self):
        ImageScalars = self.MBF.GetPointData().GetArray(self.args.ArrayName)
        per_75th = np.percentile(vtk_to_numpy(ImageScalars), 75)
        IndexMBF = numpy_to_vtk(ImageScalars/per_75th)
        IndexMBF.SetName("IndexMBF")
        self.MBF.GetPointData().AddArray(IndexMBF)
        OPath = f"{os.path.splitext(self.args.InputMBFMap)[0]}_Normalized.vtu"
        WriteVTUFile(OPath, self.MBF)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputMBFMap", "--InputMBFMap", dest = "InputMBFMap", type = str, required = True)
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = str, required = False, default = "ImageScalars")
    args = parser.parse_args()
    
    MBFNormalization(args).Normalize()