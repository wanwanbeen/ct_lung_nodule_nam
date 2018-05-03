import numpy as np
import nibabel as nib
import pandas
import math
import random  
import os
import csv
from glob import glob
from bs4 import BeautifulSoup
from scipy import ndimage

random.seed(1321)
np.random.seed(1321)

dir_lidc 	= './output_nii_mask/'
dir_info 	= './info/'
dir_label 	= './label/'
dir_mal     = './output_info_malignancy/'
dir_out     = './output_info_mask/'

with open('LIDC-IDRI_MetaData_ct.csv', 'rb') as csvfile:
    lidc_info = csv.DictReader(csvfile)
    PatientID = []
    SeriesUID = []
    for row in lidc_info:
        PatientID.append(row['PatientId'])
        SeriesUID.append(row['SeriesUID'])

def read_reading_session(reading_sessions, M, M_aff, info):

    n_reading = 0

    for reading_session in reading_sessions:
        n_reading += 1
        print idno + ' ' + str(n_reading)

        Mal         = np.int16([])
        Mtmp 		= np.zeros(M.shape)
        nodules 	= reading_session.find_all("unblindedReadNodule")
        label_no 	= 1
        for nodule in nodules:
            Mtmptmp   = np.zeros(M.shape)
            rois = nodule.find_all("roi")
            x_min = y_min = z_min = 999999
            x_max = y_max = z_max = -999999
            if len(rois) < 2:
                continue

            for roi in rois:
                z_pos = float(roi.imageZposition.text)
                z_min = min(z_min, z_pos)
                z_max = max(z_max, z_pos)
                z     = int((z_pos - info[1])/info[0])

                edge_maps = roi.find_all("edgeMap")
                for edge_map in edge_maps:
                    x = 512-int(edge_map.xCoord.text)
                    y = 512-int(edge_map.yCoord.text)

                    if x >= M.shape[0] or x <= 0 \
                            or y >= M.shape[1] or y <= 0 \
                            or z >= M.shape[2] or z <= 0:
                        print "wrong continue!!!"
                        return

                    Mtmptmp[x, y, z] = 1
                    x_min = min(x_min, x)
                    y_min = min(y_min, y)
                    x_max = max(x_max, x)
                    y_max = max(y_max, y)

                Mtmptmp[:,:,z] = ndimage.binary_fill_holes(Mtmptmp[:,:,z])

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
                malignacy = sphericiy = margin = spiculation = texture = calcification = internal_structure = lobulation = subtlety = 100
            elif nodule.characteristics.malignancy is None:
                malignacy = sphericiy = margin = spiculation = texture = calcification = internal_structure = lobulation = subtlety = 100
            else:
                malignacy          	= nodule.characteristics.malignancy.text
                sphericiy          	= nodule.characteristics.sphericity.text
                margin             	= nodule.characteristics.margin.text
                spiculation        	= nodule.characteristics.spiculation.text
                texture            	= nodule.characteristics.texture.text
                calcification      	= nodule.characteristics.calcification.text
                internal_structure 	= nodule.characteristics.internalStructure.text
                lobulation 			= nodule.characteristics.lobulation.text
                subtlety 			= nodule.characteristics.subtlety.text

            mal = [label_no, 1, x_center, y_center, z_center,
                                malignacy, sphericiy, margin, spiculation, texture, calcification, internal_structure, lobulation, subtlety]
            mal = np.int16(mal)
            Mal = np.append(Mal, mal, axis=0)
            print idno+' '+str(n_reading)+' nodule'+str(label_no)+' '+str(x_center)+' '+str(y_center)+' '+str(z_center)+' '+str(malignacy)
            Mtmp = Mtmp + Mtmptmp * label_no
            label_no += 1

        nonNodules 	= reading_session.find_all("nonNodule")

        for nonNodule in nonNodules:
            Mtmptmp = np.zeros(M.shape)
            z_center = float(nonNodule.imageZposition.text)
            z_center -= info[1]
            z_center /= info[0]
            z_center = int(z_center)
            x_center = 512-int(nonNodule.locus.xCoord.text)
            y_center = 512-int(nonNodule.locus.yCoord.text)

            if x_center >= M.shape[0] or x_center <= 0 \
                    or y_center >= M.shape[1] or y_center <= 0 \
                    or z_center >= M.shape[2] or z_center <= 0:
                print "wrong continue!!!"
                return

            Mtmptmp[x_center, y_center, z_center] = 1
            Mtmp = Mtmp + Mtmptmp * 100

            mal = [label_no, 0, x_center, y_center, z_center,
                   0, 0, 0, 0, 0, 0, 0, 0, 0]
            mal = np.int16(mal)
            Mal = np.append(Mal, mal, axis=0)
            print idno+' '+str(n_reading)+' non-nodule '+str(x_center)+' '+str(y_center)+' '+str(z_center)
            label_no += 1

        Mal = np.reshape(Mal, [int(Mal.shape[0]/14),14])
        np.save(dir_mal + idno + '_' + str(n_reading) + '.npy', np.int16(Mal))
        nib.Nifti1Image(np.int8(Mtmp), M_aff).to_filename(dir_out+idno + '_' + str(n_reading) + '.nii.gz')

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
                    M 		= nib.load(dir_lidc+idno+'_MASK.nii.gz')
                    M_aff	= M.affine
                    M 		= M.get_data()
                    info    = np.load(dir_info+idno+'.npy')
                    reading_sessions = xml.LidcReadMessage.find_all("readingSession")
                    read_reading_session(reading_sessions, M, M_aff, info)
