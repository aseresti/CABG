import os
import vtk
import argparse

def convert_vtu_to_vtk(input_file, output_file):
    # Read the .vtu file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_file)
    reader.Update()

    # Get the unstructured grid data
    unstructured_grid = reader.GetOutput()

    # Write the data to a .vtk file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(unstructured_grid)
    writer.Write()

def main(args):
    """
    Takes the folder containing the .vtu files and converts them to .vtk files in the output folder.
    """

    # Iterate over all .vtu files in the input folder
    for filename in os.listdir(args.input_folder):
        if filename.endswith(".vtu"):
            input_file = os.path.join(args.input_folder, filename)
            output_file = os.path.join(args.output_folder, filename.replace(".vtu", ".vtk"))
            convert_vtu_to_vtk(input_file, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .vtu files to .vtk files.")
    parser.add_argument("-input_folder", type=str, required=True, help="Path to the folder containing .vtu files.")
    parser.add_argument("-output_folder", type=str, required=True, help="Path to the folder where .vtk files will be saved.")
    
    args = parser.parse_args()
    main(args)