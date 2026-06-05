import vtk
import argparse
from utilities import ReadVTPFile

def ComputeVolume(Model):
    tri_filter = vtk.vtkTriangleFilter()
    tri_filter.SetInputData(Model)
    tri_filter.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(tri_filter.GetOutput())
    cleaner.Update()
    
    mass = vtk.vtkMassProperties()
    mass.SetInputData(cleaner.GetOutput())
    mass.Update()

    return mass.GetVolume()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputModel", "--InputModel", type=str, required=True, dest = "InputModel")
    args = parser.parse_args()

    Model = ReadVTPFile(args.InputModel)
    Volume = ComputeVolume(Model)
    print("The Volume is: ", Volume)
