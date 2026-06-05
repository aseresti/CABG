#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:22 2023
Updated on Wed Jun 12 2024
Updated on Thu Oct 31 2024

The purpose of this script is to take the dicom file of the CT images and convert them into vti
images. User needs to specify the number of volumes in the stack of images refered to as NofCycle.
If the image stack is CTA the NofCycel is 0.

@author: aseresti@github.com
"""

import glob
import os
import argparse
import vtk
#import SimpleITK as sitk
#from vtk.util.numpy_support import numpy_to_vtk
#import numpy as np


class ConvertDicomtoVTI():
    def __init__(self,Args):
        self.Args = Args
    
    def convert(self,DCM1:str,OutputPath:str)->None:
        os.system(f'vmtkimagereader -ifile {DCM1} --pipe vmtkimagewriter -ofile {OutputPath}')
    
    """
    def ReadDicomSeries(self, data_directory:str):
    
        series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(data_directory)
        if series_IDs:
            print("---------Reading Dicom file")

        print("---------Executing SITK image")
        reader = sitk.ImageSeriesReader()
        names = reader.GetGDCMSeriesFileNames(data_directory, series_IDs[0])
        reader.SetFileNames(names)
        return reader.Execute()
    
    def sitk2vtk(self, sitk_image:sitk.Image):
        image_array = sitk.GetArrayFromImage(sitk_image)
        print("---------Converting sitk image into vtk")
        
        vtk_image = vtk.vtkImageData()
        vtk_image.SetDimensions(image_array.shape)
        vtk_image.SetSpacing(sitk_image.GetSpacing())
        vtk_image.SetOrigin(sitk_image.GetOrigin())

        vtk_array = numpy_to_vtk(image_array.ravel(), deep=True, array_type = vtk.VTK_FLOAT)
        vtk_image.GetPointData().SetScalars(vtk_array)

        return vtk_image
    
    def WriteVTK(self, Image):
        print("---------Writing VTK image")
        writer = vtk.vtkDataSetWriter()#vtkXMLImageDataWriter()
        writer.SetFileName(self.Args.OutputFileName)
        writer.SetInputData(Image)
        writer.Write()
    """

    def main(self):
        filenames = glob.glob(f'{self.Args.InputFolderName}/*.dcm')
        filenames = sorted(filenames)
        
        if self.Args.NumberOfCycles == 0:
            print("Converting CTA DCM data into vti")
            self.convert(filenames[0],f"./{self.Args.OutputFileName}")

            #sitk_image = self.ReadDicomSeries(self.Args.InputFolderName)
            #vtk_image = self.sitk2vtk(sitk_image)
            #self.WriteVTK(vtk_image)

            
        
        else:
            self.N_file_per_cycle = int(len(filenames)/self.Args.NumberOfCycles)
            for i in range(0,self.Args.NumberOfCycles):
                print("Converting CT-MPI DCM data into vtk")

                directoryname = f'perfusion_image_cycle_{i+1}'
                pathDicom = f'{self.Args.InputFolderName}/{directoryname}'
                os.system(f"mkdir {pathDicom}")
                for j in range((i)*self.N_file_per_cycle,(i+1)*self.N_file_per_cycle-1):
                    os.system(f'cp {filenames[j]} {pathDicom}')
            
                print(f'--- Looping over cycle: {i}')
                filenames_ = glob.glob(f'{pathDicom}/*.dcm')
                filenames_ = sorted(filenames_)
                self.convert(filenames_[0], f'./CTMPI_Image_{i+1}.vtk')
                
                #sitk_image = self.ReadDicomSeries(pathDicom)
                #vtk_image = self.sitk2vtk(sitk_image)

                self.output_file_path = f"./CTMPI_Image_{i}.vtk"
                
                #class_arguments = argparse.Namespace(InputSurface = None, InputImage = None)                
                #ConvertModel2LabelMap(class_arguments).WriteVTK(vtk_image)
                #self.WriteVTK(vtk_image)

                os.system(f"rm -rf {pathDicom}")
                
            

if __name__ == '__main__':
    #descreption
    parser = argparse.ArgumentParser(description='Thsi script takes a dicom folder with N cycles and outputs an averaged vti image')
    #Input
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = 'InputFolderName', help = 'The name of the folder with all of the dicom files')
    #NumberOfCycles
    parser.add_argument('-NofCycle', '--NofCycle', type = int, required = True, dest = 'NumberOfCycles', help = 'The number of perfusion images that are in the dicom folder. if CTA put 0')

    parser.add_argument('-OutputFileName', '--OutputFileName', type = str, required=True, dest = "OutputFileName")
    args = parser.parse_args()
    ConvertDicomtoVTI(args).main()
