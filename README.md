**A low resolution example has been prepared which can be run with registration_pipeline_example.m**


# Registration



Registration code involves modules:

1. Atlas free rigid slice alignment

2. 3D to 3D deformation and affine

3. 3D to 2D slices

4. 2D to 2D rigid slice alignment

These can be combined into various pipelines.  For example, 3D STPT brains wil use 2 only. RNAseq Nissl sections will use 1 for initialization, followed by 3.  Alternating nissl fluoro will use 1 for initialization, followed by 3 for nissl slices, and then 4 for other fluoro slices.



The code in this repository will be matlab based.

To run the registration:
1. make sure your down-sampled images for one brain are in the same folder.
2. run create_location_csv_MBA.m with correct meta information, this will create a geometry file saving location metadata for each sections. The geometry.csv will be save in the same folder.
3. run registration_pipeline_MBA.m with the correct atlas and segmentation file, with the input folder being the folder in step 1. The output will be saved in 'output_dir'
