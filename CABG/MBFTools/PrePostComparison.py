import argparse
import numpy as np
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from utilities import ReadVTUFile, ThresholdInBetween

class PrePostMBFMap():
    def __init__(self, args):
        self.MBF_A = ReadVTUFile(f"{args.InputFolder}/{args.InputMBF}")
        self.MBF_B = ReadVTUFile(f"{args.InputFolder[:-1]}B/{args.InputMBF}")
        self.InputLabels = f"{args.InputFolder[:-1]}B/{args.InputLabels}"

    def ReadTerritoryMBF(self, MBFMap, MBF_Labels):
        MBF_data = {"LAD": np.array([]), "LCx": np.array([]), "Intermedius": np.array([]), "Diag1": np.array([]), "Diag2": np.array([]), "PDA": np.array([]), "PL": np.array([])}
        for key in MBF_Labels.keys():
            for i in MBF_Labels[key]:
                territory_b = ThresholdInBetween(MBFMap, "TerritoryMaps", i, i+1)
                indexMBF_ = vtk_to_numpy(territory_b.GetPointData().GetArray(0))
                MBF_data[key] = np.append(indexMBF_, MBF_data[key])
        
        return MBF_data


    def PlotBox(self, MBF_Labels, MBF_data_pre, MBF_data_post):
        color_list = ['aquamarine', 'sandybrown', 'palegreen', 'lightcyan', 'thistle', 'lavender', 'salmon', 'peachpuff']
        
        Labels = []
        IndexMBF = []
        colors = []
        for (i, key) in enumerate(MBF_Labels.keys()):
            if len(MBF_data_pre[key]) > 0:
                Labels.append(f'{key}_pre')
                Labels.append(f'{key}_post')
                IndexMBF.append(MBF_data_pre[key])
                IndexMBF.append(MBF_data_post[key])
                colors.append(color_list[i])
                colors.append(color_list[i])
            
        _, ax = plt.subplots()
        ax.set_ylabel("MBF (ml/min/100g)", fontdict={'color':'blue','size':20})
        bplot = ax.boxplot(IndexMBF, patch_artist=True, labels=Labels, showfliers= False)

        for patch, color in zip(bplot['boxes'], colors):
            patch.set_color(color)

        for median in bplot['medians']:
            median.set_color('darkblue')
            median.set_linewidth(2)

        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize = 13)
        plt.tight_layout()
        plt.show()



    def main(self):
        MBF_Labels = {"LAD": [], "LCx":[], "Intermedius":[], "Diag1":[], "Diag2":[], "PDA":[], "PL":[]}

        with open(self.InputLabels, "r") as ifile:
            for LINE in ifile:
                line = LINE.split()
                for key in MBF_Labels.keys():
                    if line[1].find(key)>=0: MBF_Labels[key].append(int(line[0]))

        MBF_data_pre = self.ReadTerritoryMBF(self.MBF_A, MBF_Labels)
        MBF_data_post = self.ReadTerritoryMBF(self.MBF_B, MBF_Labels)

        self.PlotBox(MBF_Labels, MBF_data_pre, MBF_data_post)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputFolderPre", "--InputFolderPre", type=str, required= True, dest= "InputFolder")
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= False, default= "MBF_Territories.vtu")
    parser.add_argument("-InputLabels", "--InputLabels", dest= "InputLabels", type= str, required= False, default= "MBF_Territories_Labels.dat")
    args = parser.parse_args()

    PrePostMBFMap(args).main()