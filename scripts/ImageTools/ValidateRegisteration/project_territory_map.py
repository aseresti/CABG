"""
This script takes the Elastix transform .txt file and inverts it.
Then, it applies the inverted transform to the regitered MBF map to get the original map in the original space.
Next, it uses a prob filter to project the territory labels in the re-registered map to the original MBF map.
Finally, the projected file is saved as a .vtu file.
"""

import os
import vtk
import argparse
import numpy as np
import SimpleITK as sitk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

def read_h5_transform(transform_file):
    # Read the H5 transform file and extract the parameters
    transform = sitk.ReadTransform(transform_file)
    inverted_transform = transform.GetInverse()
    return inverted_transform

def read_vtu_file(vtu_file):
    # Read the .vtu file using VTK
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()
    return reader.GetOutput()

def write_vtu_file(vtu_data, output_file):
    # Write the .vtu file using VTK
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(vtu_data)
    writer.Write()

def apply_transform_to_vtu(vtu_data, transform):
    # Convert sitk transform to VTK transform
    matrix = transform.GetParameters()[:9]
    t = transform.GetParameters()[9:12]
    c = transform.GetFixedParameters()

    m = [matrix[0:3], matrix[3:6], matrix[6:9]]
    
    # T_adjusted = T - c + np.dot(m, c)
    # T_adjusted = np.array(t) - np.array(c) + np.dot(np.array(m), np.array(c))
    # vt_t_adjusted = numpy_to_vtk(T_adjusted, deep=True)
    # vtk_transform = vtk.vtkTransform()
    # vtk_transform.SetMatrix([m[0][0], m[0][1], m[0][2], 0,
    #                          m[1][0], m[1][1], m[1][2], 0,
    #                          m[2][0], m[2][1], m[2][2], 0,
    #                          0, 0, 0, 1])
    # vtk_transform.Translate(T_adjusted)


    vtk_mat = vtk.vtkMatrix4x4()
    vtk_mat.Identity()

    for i in range(3):
        for j in range(3):
            vtk_mat.SetElement(i, j, m[i][j])

    
    t_adj = [0.0, 0.0, 0.0]
    for i in range(3):
        t_adj[i] = t[i] + c[i] - (m[i][0]*c[0] + m[i][1]*c[1] + m[i][2]*c[2])

    vtk_mat.SetElement(0, 3, t_adj[0])
    vtk_mat.SetElement(1, 3, t_adj[1])
    vtk_mat.SetElement(2, 3, t_adj[2])

    # sitk and vtk use different coordinate systems
    lps_to_ras = vtk.vtkMatrix4x4()
    lps_to_ras.Identity()
    lps_to_ras.SetElement(0, 0, 1.0)  
    lps_to_ras.SetElement(-1, 1, 1.0)

    final_matrix = vtk.vtkMatrix4x4()
    vtk.vtkMatrix4x4.Multiply4x4(lps_to_ras, vtk_mat, final_matrix)
    vtk.vtkMatrix4x4.Multiply4x4(final_matrix, lps_to_ras, final_matrix)

    vtk_transform = vtk.vtkTransform()
    vtk_transform.SetMatrix(final_matrix)

    # Apply the given transform using VTK's vtkTransform and vtkTransformFilter
    transform_filter = vtk.vtkTransformFilter()
    transform_filter.SetInputData(vtu_data)
    transform_filter.SetTransform(vtk_transform)
    transform_filter.Update()
    return transform_filter.GetOutput()

def project_labels_to_original_map(re_registered_mbf_map, original_mbf_map):
    # todo: use a different method
    source_points = re_registered_mbf_map.GetPoints()
    kd_tree = vtk.vtkKdTreePointLocator()
    kd_tree.SetDataSet(re_registered_mbf_map)
    kd_tree.BuildLocator()

    source_array = re_registered_mbf_map.GetPointData().GetArray("TerritoryMaps")
    target_array = vtk.vtkDoubleArray()
    target_array.SetName("ProjectedTerritoryMaps")
    target_array.SetNumberOfComponents(1)
    target_array.SetNumberOfTuples(original_mbf_map.GetNumberOfPoints())

    target_points = original_mbf_map.GetPoints()
    num_target_points = target_points.GetNumberOfPoints()
    for i in range(num_target_points):
        pt = target_points.GetPoint(i)
        # Find the closest point ID in the source mesh
        closest_id = kd_tree.FindClosestPoint(pt)
        
        # Extract data value and insert it into the target array
        val = source_array.GetTuple(closest_id)
        target_array.InsertNextTuple(val)

    original_mbf_map.GetPointData().AddArray(target_array)
    return original_mbf_map

def main(args):
    # Read the transform parameters
    inverted_transform = read_h5_transform(args.transform_file)
    registered_mbf_map = read_vtu_file(args.registered_mbf_map)
    original_mbf_map = read_vtu_file(args.original_mbf_map)

    
    # Apply the inverted transform to the registered MBF map (this is a placeholder, actual application logic will depend on the type of transform and data)
    # For example, you might use vtkTransform to apply the transformation to the data
    re_registered_mbf_map = apply_transform_to_vtu(registered_mbf_map, inverted_transform)
    output_folder = os.path.dirname(args.registered_mbf_map)
    write_vtu_file(re_registered_mbf_map, os.path.join(output_folder, "re_registered_mbf_map.vtu"))

    
    # Use a prob filter to project the territory labels in the re-registered map to the original MBF map (this is a placeholder, actual projection logic will depend on your specific requirements)
    projected_map = project_labels_to_original_map(re_registered_mbf_map, original_mbf_map)
    
    # Save the projected file as a .vtu file (this is a placeholder, actual saving logic will depend on your specific requirements)
    write_vtu_file(projected_map, os.path.join(output_folder, "projected_territory_map.vtu"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Invert Elastix transform and apply it to the registered MBF map.")
    parser.add_argument("-transform_file", type=str, required=True, help="Path to the Elastix transform .txt file.")
    parser.add_argument("-registered_mbf_map", type=str, required=True, help="Path to the registered MBF map file.")
    parser.add_argument("-original_mbf_map", type=str, required=True, help="Path to the original MBF map file.")

    
    args = parser.parse_args()
    main(args)