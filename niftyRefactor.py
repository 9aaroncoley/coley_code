import numpy as np
import nibabel as nib
import os

# class for parsing nifty files and normalizing them
class NiftyParser:
    threshold = 500
    mask_num = 1
    dataPathName = 'data/'
    resultsPathName = 'results/'
    normalized = "_norm"
    masked = "_mask"
    fileID = '.nii.gz'
    
    fileName = ''
    patientDir = ''

    def __init__(self, threshold):
        self.threshold = threshold

    def refactorFile(self, fileName, patientDir):
        self.fileName = fileName
        self.patientDir = patientDir

        #Load the nifty file data into usable data structures
        brain_vol = nib.load(self.dataPathName+self.patientDir+self.fileName)
        mask = nib.load(self.dataPathName+self.patientDir+self.fileName)
        self.brain_vol_data = brain_vol.get_fdata()
        self.mask_data = mask.get_fdata()
        
        numScans, rows, cols = self.brain_vol_data.shape

        #Loop over the 3D brain scan going one 2D layer at a time
        for i in range(numScans):
            for j in range(cols):
                for k in range(rows):
                    if self.brain_vol_data[i][k][j] > self.threshold:  #see three pixels in a row or something
                        break
                    self.brain_vol_data[i][k][j] = 0
                    self.mask_data[i][k][j] = 0

            for j in range(cols-1, -1, -1):
                for k in range(rows-1, -1, -1):
                    if self.brain_vol_data[i][k][j] > self.threshold:
                        break
                    self.brain_vol_data[i][k][j] = 0
                    self.mask_data[i][k][j] = 0

        self.createMask(numScans, rows, cols)
        self.saveFile()
    
    def createMask(self, numScans, rows, cols):
        for i in range(numScans):
            for j in range(cols):
                for k in range(rows):
                    if self.mask_data[i][k][j] != 0:
                        self.mask_data[i][k][j] = self.mask_num

    def saveFile(self):
        nii_img = nib.Nifti1Image(self.brain_vol_data, affine=np.eye(4))
        mask_img = nib.Nifti1Image(self.mask_data, affine=np.eye(4))
        caseName = self.fileName[:-len(self.fileID)]
        nii_location = self.resultsPathName+self.patientDir+caseName+self.normalized+self.fileID
        mask_location = self.resultsPathName+self.patientDir+caseName+self.masked+self.fileID
        nib.save(nii_img, nii_location)
        nib.save(mask_img, mask_location)

if __name__ == "__main__":
    parser = NiftyParser(threshold=500)

    patientFolders = [f+'/' for f in os.listdir("data/") if os.path.isdir("data/"+f)]
    for patientFolder in patientFolders:
        files = [f for f in os.listdir("data/"+patientFolder) if os.path.isfile("data/"+patientFolder+f)]
        for file in files:
            parser.refactorFile(fileName=file, patientDir=patientFolder)
            # print("Currently working on --> Patient: %s" % patientFolder, "File: %s" % file



# TO DO
# add superimposing on normalized files