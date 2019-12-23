import platform
import os
from os import listdir
from os.path import isfile, join

scriptLocation = os.path.dirname(os.path.realpath(__file__))

microarrayFolder = os.path.join(scriptLocation, "data", "raw", "allen_HBA") 
DebNIFTIfile = os.path.join(scriptLocation, "data", "PET_meta_analyses_images", "DebMRIfiles", "FDG-PET Coordinates Final_ALE.nii")
