"""
"""
import os
import sys
import vtk
vtk.vtkObject.GlobalWarningDisplayOff()
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from project_territory_map import read_vtu_file
sys.path.append(os.path.join(os.path.curdir, "scripts"))
from utilities import ThresholdInBetween
from vtk.util.numpy_support import vtk_to_numpy
sns.set_context("notebook")

def ReadLabels(InputLabels, TerritoryTag):
    MBF_Labels = {}
    for tag in TerritoryTag:
        MBF_Labels[tag] =  []
    MBF_Labels["Non-grafted\nstenosis>50%"] = []
    MBF_Labels["Non-grafted\nstenosis<50%"] = []
    with open(InputLabels, "r") as ifile:
        for i, LINE in enumerate(ifile):
            if i == 0: 
                continue
            line = LINE.strip().split()
            label = line[1]
            id_value = int(line[0])
            found = False

            for key in TerritoryTag:
                if key in label: 
                    MBF_Labels[key].append(id_value)
                    found = True
            
            if not found:
                if "NG" in label:
                    MBF_Labels["Non-grafted\nstenosis>50%"].append(id_value)
                    found = True
                else:
                    MBF_Labels["Non-grafted\nstenosis<50%"].append(id_value)

    MBF_Labels = {k.replace('post_', ''):v for k,v in MBF_Labels.items() if len(v)>0}
    
    return MBF_Labels

def CollectMBFData(MBF, Labels):
    for i in range(MBF.GetPointData().GetNumberOfArrays()):
            arrayname_ = MBF.GetPointData().GetArrayName(i)
            if 'scalars' in arrayname_.lower():
                ScalarArray = arrayname_

    AbsMBFData = {}
    
    for key in Labels.keys():
        AbsMBFData[key] = np.array([])
    
        for i in Labels[key]:
            territory_ = ThresholdInBetween(MBF, "TerritoryMaps", i, i)
            MBF_ = vtk_to_numpy(territory_.GetPointData().GetArray(ScalarArray))
            
            AbsMBFData[key] = np.append(AbsMBFData[key], MBF_)
            

    return AbsMBFData

def GetScalarArrayName (Data):
    for i in range(Data.GetPointData().GetNumberOfArrays()):
        array_name_ = Data.GetPointData().GetArrayName(i)
        if 'scalar' in array_name_.lower():
            return array_name_

def GetMBFArray(VTUVolume):
    ScalaraArrayName = GetScalarArrayName(VTUVolume)
    MBF_array = vtk_to_numpy(VTUVolume.GetPointData().GetArray(ScalaraArrayName))
    return MBF_array

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, dest="case", type=str)
    parser.add_argument("--tags", required=True, dest="tags", type= str, nargs="+")

    args = parser.parse_args()
    
    dir_ = "/Volumes/T7_Research/CABG/05_PrePostCABG/Validate_Registeration"
    pre_CABG_MBF = read_vtu_file(os.path.join(dir_, "pre-CABG", f"{args.case}A.vtu"))
    post_CABG_MBF = read_vtu_file(os.path.join(dir_, "re-registered", f"{args.case}B.vtu"))
    registered_mbf_map = read_vtu_file(os.path.join(dir_, "registered", f"{args.case}B.vtu"))
    label_file = os.path.join(dir_, "label-file", f"{args.case}.dat")

    tags = args.tags[0].split("+")
    tags = [f'post_{tag}' for tag in tags]
    MBF_Labels = ReadLabels(label_file, tags)
    print("MBF Labels:", MBF_Labels)
    pre_AbsMBFData = CollectMBFData(pre_CABG_MBF, MBF_Labels)
    print("Data pre-CABG:", pre_AbsMBFData)
    post_AbsMBFData = CollectMBFData(post_CABG_MBF, MBF_Labels)
    registered_AbsMBFData = CollectMBFData(registered_mbf_map, MBF_Labels)

    df_A = pd.DataFrame({
        'MBF': np.concatenate([pre_AbsMBFData[key] for key in pre_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] * len(pre_AbsMBFData[key]) for key in pre_AbsMBFData.keys()]),
        'Condition': 'Pre-CABG'
    })

    df_B = pd.DataFrame({
        'MBF': np.concatenate([post_AbsMBFData[key] for key in post_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] * len(post_AbsMBFData[key]) for key in post_AbsMBFData.keys()]),
        'Condition': 'Post-CABG'
    })

    df_C = pd.DataFrame({
        'MBF': np.concatenate([registered_AbsMBFData[key] for key in registered_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] * len(registered_AbsMBFData[key]) for key in registered_AbsMBFData.keys()]),
        'Condition': 'Registered'
    })

    df_long = pd.concat([df_A, df_B, df_C], ignore_index=True)

    
    g = sns.catplot(
    data=df_long, 
    x='Territory', 
    y='MBF', 
    hue='Condition', 
    # inner=None, 
    color=".8", 
    kind="violin",
    palette="Set2"
    )
    g.fig.set_size_inches(10, 6.5) 
    ax = g.axes[0, 0]
    ax.tick_params(axis='x', labelrotation=45)
    

    df_A = pd.DataFrame({
        'Mean': np.array([np.mean(pre_AbsMBFData[key]) for key in pre_AbsMBFData.keys()]),
        'Median MBF': np.array([np.median(pre_AbsMBFData[key]) for key in pre_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] for key in pre_AbsMBFData.keys()]),
        'Condition': 'Pre-CABG'
    })

    df_B = pd.DataFrame({
        'Mean': np.array([np.mean(post_AbsMBFData[key]) for key in post_AbsMBFData.keys()]),
        'Median MBF': np.array([np.median(post_AbsMBFData[key]) for key in post_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] for key in post_AbsMBFData.keys()]),
        'Condition': 'Post-CABG'
    })

    df_C = pd.DataFrame({
        'Mean': np.array([np.mean(registered_AbsMBFData[key]) for key in registered_AbsMBFData.keys()]),
        'Median MBF': np.array([np.median(registered_AbsMBFData[key]) for key in registered_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] for key in registered_AbsMBFData.keys()]),
        'Condition': 'Registered'
    })

    df_long_median = pd.concat([df_A, df_B, df_C], ignore_index=True)

    
    g = sns.relplot(
        data=df_long_median, 
        y='Territory', 
        x='Median MBF', 
        hue='Condition', 
        kind='scatter',
        # color=".8"
        palette="Set2" 
        )
    ax = g.axes[0, 0]
    ax.set_xlim(50, 250)
    g.fig.set_size_inches(7, 5) 
    plt.tight_layout()


    MBF_array_A = GetMBFArray(pre_CABG_MBF)
    ref_A = np.percentile(MBF_array_A, 75)
    MBF_array_B = GetMBFArray(post_CABG_MBF)
    ref_B = np.percentile(MBF_array_B, 75)
    MBF_array_C = GetMBFArray(registered_mbf_map)
    ref_C = np.percentile(MBF_array_C, 75)

    pre_AbsMBFData_index = {key:mbf/ref_A for key, mbf in pre_AbsMBFData.items()}
    post_AbsMBFData_index = {key:mbf/ref_B for key, mbf in post_AbsMBFData.items()}
    registered_AbsMBFData_index = {key:mbf/ref_C for key, mbf in registered_AbsMBFData.items()}
    

    df_A = pd.DataFrame({
        'Index MBF': np.concatenate([pre_AbsMBFData_index[key] for key in pre_AbsMBFData_index.keys()]),
        'Territory': np.concatenate([[key] * len(pre_AbsMBFData_index[key]) for key in pre_AbsMBFData_index.keys()]),
        'Condition': 'Pre-CABG'
    })

    df_B = pd.DataFrame({
        'Index MBF': np.concatenate([post_AbsMBFData_index[key] for key in post_AbsMBFData_index.keys()]),
        'Territory': np.concatenate([[key] * len(post_AbsMBFData[key]) for key in post_AbsMBFData.keys()]),
        'Condition': 'Post-CABG'
    })

    df_C = pd.DataFrame({
        'Index MBF': np.concatenate([registered_AbsMBFData_index[key] for key in registered_AbsMBFData_index.keys()]),
        'Territory': np.concatenate([[key] * len(registered_AbsMBFData[key]) for key in registered_AbsMBFData.keys()]),
        'Condition': 'Registered'
    })

    df_long = pd.concat([df_A, df_B, df_C], ignore_index=True)

    
    g = sns.catplot(
    data=df_long, 
    x='Territory', 
    y='Index MBF', 
    hue='Condition', 
    # inner=None, 
    color=".8", 
    kind="violin",
    palette="Set2"
    )
    g.fig.set_size_inches(10, 6.5) 
    ax = g.axes[0, 0]
    ax.tick_params(axis='x', labelrotation=45)


    df_A = pd.DataFrame({
        'Mean': np.array([np.mean(pre_AbsMBFData_index[key]) for key in pre_AbsMBFData.keys()]),
        'Median Index MBF': np.array([np.median(pre_AbsMBFData_index[key]) for key in pre_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] for key in pre_AbsMBFData.keys()]),
        'Condition': 'Pre-CABG'
    })

    df_B = pd.DataFrame({
        'Mean': np.array([np.mean(post_AbsMBFData_index[key]) for key in post_AbsMBFData.keys()]),
        'Median Index MBF': np.array([np.median(post_AbsMBFData_index[key]) for key in post_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] for key in post_AbsMBFData.keys()]),
        'Condition': 'Post-CABG'
    })

    df_C = pd.DataFrame({
        'Mean': np.array([np.mean(registered_AbsMBFData_index[key]) for key in registered_AbsMBFData.keys()]),
        'Median Index MBF': np.array([np.median(registered_AbsMBFData_index[key]) for key in registered_AbsMBFData.keys()]),
        'Territory': np.concatenate([[key] for key in registered_AbsMBFData.keys()]),
        'Condition': 'Registered'
    })

    df_long_median = pd.concat([df_A, df_B, df_C], ignore_index=True)

    g = sns.relplot(
        data=df_long_median, 
        y='Territory', 
        x='Median Index MBF', 
        hue='Condition', 
        kind='scatter',
        # color=".8"
        palette="Set2" 
        )
    ax = g.axes[0, 0]
    ax.set_xlim(0, 2)
    g.fig.set_size_inches(5, 5) 
    plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    main()