import numpy as np
import nibabel as nib
import pandas
import math
import random  
import os
import csv
from glob import glob
from bs4 import BeautifulSoup

random.seed(1321)
np.random.seed(1321)

dir_lidc 	= './output_nii_mask/'
dir_info 	= './info/'
dir_label 	= './label/'
dir_mal     = './output_info_malignancy/'

with open('LIDC-IDRI_MetaData_ct.csv', 'rb') as csvfile:
    lidc_info = csv.DictReader(csvfile)
    PatientID = []
    SeriesUID = []
    for row in lidc_info:
        PatientID.append(row['PatientId'])
        SeriesUID.append(row['SeriesUID'])

def read_reading_session(reading_sessions, info):

    n_reading = 0
    for reading_session in reading_sessions:
        n_reading += 1

        Mal         = np.int16([])
        nodules 	= reading_session.find_all("unblindedReadNodule")
        label_no 	= 0
        for nodule in nodules:
            nodule_id = nodule.noduleID.text
            rois = nodule.find_all("roi")
            x_min = y_min = z_min = 999999
            x_max = y_max = z_max = -999999
            if len(rois) < 2:
                continue

            for roi in rois:
                z_pos = float(roi.imageZposition.text)
                z_min = min(z_min, z_pos)
                z_max = max(z_max, z_pos)
                edge_maps = roi.find_all("edgeMap")
                for edge_map in edge_maps:
                    x = int(edge_map.xCoord.text)
                    y = int(edge_map.yCoord.text)
                    x_min = min(x_min, x)
                    y_min = min(y_min, y)
                    x_max = max(x_max, x)
                    y_max = max(y_max, y)
                if x_max == x_min:
                    continue
                if y_max == y_min:
                    continue

            x_diameter = x_max - x_min
            x_center = x_min + x_diameter / 2
            y_diameter = y_max - y_min
            y_center = y_min + y_diameter / 2
            z_diameter = z_max - z_min
            z_center = z_min + z_diameter / 2
            z_center -= info[1]
            z_center /= info[0]
            z_center = int(z_center)

            if nodule.characteristics is None:
                print("!!!!Nodule:", nodule_id, " has no charecteristics")
                continue
            if nodule.characteristics.malignancy is None:
                print("!!!!Nodule:", nodule_id, " has no malignacy")
                continue

            malignacy          	= nodule.characteristics.malignancy.text
            sphericiy          	= nodule.characteristics.sphericity.text
            margin             	= nodule.characteristics.margin.text
            spiculation        	= nodule.characteristics.spiculation.text
            texture            	= nodule.characteristics.texture.text
            calcification      	= nodule.characteristics.calcification.text
            internal_structure 	= nodule.characteristics.internalStructure.text
            lobulation 			= nodule.characteristics.lobulation.text
            subtlety 			= nodule.characteristics.subtlety.text
            label_no += 1
            mal = [label_no, 1, x_center, y_center, z_center,
                                malignacy, sphericiy, margin, spiculation, texture, calcification, internal_structure, lobulation, subtlety]
            mal = np.int16(mal)
            Mal = np.append(Mal, mal, axis=0)
            print idno+' '+str(n_reading)+' nodule'+str(label_no)+' '+str(x_center)+' '+str(y_center)+' '+str(z_center)+' '+str(malignacy)

        nonNodules 	= reading_session.find_all("nonNodule")
        for nonNodule in nonNodules:
            z_center = float(nonNodule.imageZposition.text)
            z_center -= info[1]
            z_center /= info[0]
            z_center = int(z_center)
            x_center = int(nonNodule.locus.xCoord.text)
            y_center = int(nonNodule.locus.yCoord.text)
            label_no += 1
            mal = [label_no, 0, x_center, y_center, z_center,
                   0, 0, 0, 0, 0, 0, 0, 0, 0]
            mal = np.int16(mal)
            Mal = np.append(Mal, mal, axis=0)
            print idno+' '+str(n_reading)+' non-nodule '+str(x_center)+' '+str(y_center)+' '+str(z_center)

        Mal = np.reshape(Mal, [int(Mal.shape[0]/14),14])
        np.save(dir_mal + idno + '_' + str(n_reading) + '.npy', Mal)

xml_all = sorted(glob(dir_label+'*/*/*.xml'))
for xml_path in xml_all:
    with open(xml_path, 'r') as xml_file:
        markup = xml_file.read()
        xml = BeautifulSoup(markup, features="xml")
        if xml.LidcReadMessage is not None:
            patient_id = xml.LidcReadMessage.ResponseHeader.SeriesInstanceUid.text
            if patient_id in SeriesUID:
                idno    = str(PatientID[SeriesUID.index(patient_id)])
                if os.path.exists(dir_info+idno+'.npy'):
                    info    = np.load(dir_info+idno+'.npy')
                    reading_sessions = xml.LidcReadMessage.find_all("readingSession")
                    read_reading_session(reading_sessions, info)
