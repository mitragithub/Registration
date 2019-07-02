# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 17:27:48 2016

@author: balamurali
"""
import jsonontology_parse as jp
import SimpleITK as sitk
from skimage import measure
import numpy as np
import json
import os,sys
import h5py
import re
import glob
from skimage.util import view_as_blocks
from multiprocessing import Pool

def Down_Sample(image, block_size, func=np.sum, cval=0):

    if len(block_size) != image.ndim:
        raise ValueError("`block_size` must have the same length "
                         "as `image.shape`.")

    pad_width = []
    for i in range(len(block_size)):
        if block_size[i] < 1:
            raise ValueError("Down-sampling factors must be >= 1. Use "
                             "`skimage.transform.resize` to up-sample an "
                             "image.")
        if image.shape[i] % block_size[i] != 0:
            after_width = block_size[i] - (image.shape[i] % block_size[i])
        else:
            after_width = 0
        pad_width.append((0, after_width))

    image = np.pad(image, pad_width=pad_width, mode='constant',
                   constant_values=cval)

    out = view_as_blocks(image, block_size)

    for i in range(len(out.shape) // 2):
        out = func(out, axis=-1)

    return out

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def coordinate_change(pt,siz):
    height, width = siz
    #pt_x = pt[0]/width * 24000
    #pt_y = -pt[1]/height * 18000
    pt_x = pt[0] * 4
    pt_y = -pt[1] * 4
    return [pt_x, pt_y]

def getcomponents(imgslice):
    siz = imgslice.shape

    colors = np.union1d(np.array(imgslice).ravel(),np.array([]))
  
    json_data = """{"type":"FeatureCollection","features":["""
    unknown_colors = []
    
    for ci in colors[1:]: 
        count = True
        bwarr = imgslice==ci;
        con = measure.find_contours(bwarr,0.5)
        #con = measure.subdivide_polygon(con, degree=2)
        sum_1 = 0
        (br_reg_name,depth) = jp.ontology_find(ci)

	if br_reg_name is None:
	    #continue
	    br_reg_name = {'name': '##', 'acronym': '00'}
	    unknown_colors.append(ci)

        #print br_reg_name
        mul_pol = []
        for coni in con:
            #print coni
            conxy = np.array(coni)
            index_x = np.argmin((conxy[:,1]))
            slope = (conxy[:,0][index_x+1]-conxy[:,0][index_x])/(conxy[:,1][index_x+1]-conxy[:,1][index_x]+0.00001)
            if slope < 0:
                x_coord = list(conxy[:,1])
                y_coord = list(conxy[:,0])
                temp_arr = []
                for j in range(len(x_coord)):
                    temp_arr.append(coordinate_change([x_coord[j],y_coord[j]], siz))

                if len(con)==1:
                    coord = """\n{"type":"Feature","id":"%d","properties":{"name":"%s","acronym": "%s" },"geometry":{"type":"Polygon","coordinates":["""%(ci,br_reg_name["name"], br_reg_name["acronym"])+str(temp_arr)+"""]}},"""
                    count = False

                else:
                    mul_pol.append(temp_arr)

        if count is True:
            coord = """\n{"type":"Feature","id":"%d","properties":{"name":"%s","acronym":"%s"},"geometry":{"type":"MultiPolygon","coordinates":["""%(ci,br_reg_name["name"], br_reg_name["acronym"])+str(mul_pol)+"""]}},"""


        #file_json.write(coord)
	json_data = json_data + coord

    return json_data[:-1] + "]}", set(unknown_colors)


# outdir = 'jsondata'    

def single(element):
    unknown_cols = []

    print element
    count = element[-7:-4]
    data = h5py.File(element, 'r')
    data = np.asarray(data['seg'])

    data = Down_Sample(data, ((2,2)), func = np.max)
    data = Down_Sample(data, ((2,2)), func = np.max)
    #data = Down_Sample(data, ((2,2)), func = np.max)


    json_data, ukc = getcomponents(data.T)
    unknown_cols = np.union1d(unknown_cols,list(ukc))
    file_json = open('%s/atlas_%s_%s.json' % (outdir,pmdno, count),'w')
    file_json.write(json_data)
    file_json.close()
if __name__ =="__main__":
    pmdno = sys.argv[1]
    # outdir = '/sonas-hs/mitra/hpc/home/xli/snRNA_registration_pipeline/FinalTransformation/' + pmdno + '/JSON'
    outdir = '/sonas-hs/mitra/hpc/home/xli/RegistrationPipelineV3/Data/' + pmdno + '/Transformation_OUTPUT/JSON/'
    # os.system("mkdir " + outdir + '/' + pmdno)
    os.system('mkdir '+outdir )

    sliceno = None

    # imglist = natural_sort(glob.glob('/sonas-hs/mitra/hpc/home/xli/snRNA_registration_pipeline/FinalTransformation/' + pmdno + '/reg_high_seg_pad/*.mat'))
    imglist = natural_sort(glob.glob('/sonas-hs/mitra/hpc/home/xli/RegistrationPipelineV3/Data/' + pmdno + '/Transformation_OUTPUT/reg_high_seg_pad/*N*.mat'))

    p = Pool(10)
    p.map(single, imglist)


    # unknown_cols = []
    # # count = 1
    # for element in imglist:
    #     print element
    #     count = element[-8:-4]
    #     data = h5py.File(element, 'r')
    #     data = np.asarray(data['seg'])

    #     data = Down_Sample(data, ((2,2)), func = np.max)
    #     data = Down_Sample(data, ((2,2)), func = np.max)
    #     #data = Down_Sample(data, ((2,2)), func = np.max)


    #     json_data, ukc = getcomponents(data.T)
    #     unknown_cols = np.union1d(unknown_cols,list(ukc))
    #     file_json = open('%s/%s/atlas_%s_%s.json' % (outdir,pmdno,pmdno, count),'w')
    #     file_json.write(json_data)
    #     file_json.close()
    #     # count = count + 1


