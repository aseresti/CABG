import vtk
import os
import argparse
import glob
import xml.etree.ElementTree as ET
#from TubeFilter import *


class ConvertPath2VTP():
    def __init__(self,Args):
        self.Args = Args
        self.filenames = glob.glob(os.path.join(self.Args.InputFolder,"*.pth"))

    def main(self):
        for file in self.filenames:
            lumen = os.path.splitext(os.path.basename(file))[0]
            output_path = os.path.join(f"{self.Args.InputFolder}",f"{lumen}.vtp")
            points, _ = self.pth_to_points(file)
            polydata = self.points_to_vtp(points)
            self.WriteVTPFile(output_path, polydata)
        

    def pth_to_points(self,PathFile):
        with open(PathFile, "r") as path:
            tree = ET.parse(path)
        root = tree.getroot()

        # Extract points
        points = []
        for path_point in root.findall(".//path_point/pos"):
            x = float(path_point.attrib['x'])
            y = float(path_point.attrib['y'])
            z = float(path_point.attrib['z'])
            points.append([x, y, z])
        
        direction = []
        for tangent in root.findall(".//path_point/tangent"):
            x = float(tangent.attrib['x'])
            y = float(tangent.attrib['y'])
            z = float(tangent.attrib['z'])
            direction.append([x, y, z])

        return points, direction
        
    def points_to_vtp(self, points):
        # Create VTK points
        vtk_points = vtk.vtkPoints()
        for point in points:
            vtk_points.InsertNextPoint(point)
    
        # Create a polyline
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(points))
        for i in range(len(points)):
            polyline.GetPointIds().SetId(i, i)
    
        # Create a cell array to store the polyline
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyline)
    
        # Create a polydata object
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        polydata.SetLines(cells)

        return polydata
    
    def WriteVTPFile(self, output_vtp, polydata):
        # Write to a .vtp file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_vtp)
        writer.SetInputData(polydata)
        writer.Write()


if __name__=="__main__":
    Parser = argparse.ArgumentParser()
    Parser.add_argument("-InputFolder", "--InputFolder", dest="InputFolder", type=str, required=True)
    args = Parser.parse_args()

    ConvertPath2VTP(args).main()