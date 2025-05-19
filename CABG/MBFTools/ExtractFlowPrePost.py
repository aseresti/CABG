import os
import argparse
from utilities import ReadVTUFile, ThresholdInBetween
from ExtractFlowInTerritories import ExtractSubtendedFlow

class ExtractFlowPrePost(ExtractSubtendedFlow):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    def ReadPrePostFiles(self):
        pass

    def MergeSubtendedTerritories(self):
        pass

    def CollectFlowData(self, MBFMap):
        pass

    def CollectPrePostData(self):
        pass

    def BarPlot(self, BarData):
        pass

    def BoxPlot(self, BoxData):
        pass

    def Statistics(self, MBF_Data):
        pass

    def main(self):
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-InputMBF", "--InputMBF", dest= "InputMBF", type= str, required= True)
    parser.add_argument("-ArrayName", "--ArrayName", dest = "ArrayName", type = str, required = False, default = "ImageScalars")
    parser.add_argument("-TerritoryTag", "--TerritoryTag", type= str, required=True, dest = "TerritoryTag")
    parser.add_argument("-Unit", "--Unit", type= str, dest= "Unit", default="mm", required=False)
    args = parser.parse_args()