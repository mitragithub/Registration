# Registration

Registration code will involve four modules:

1. Atlas free rigid slice alignment

2. 3D to 3D deformation and affine

3. 3D to 2D slices

4. 2D to 2D rigid slice alignment

These can be combined into various pipelines.  For example, 3D STPT brains wil use 2 only. RNAseq Nissl sections will use 1 for initialization, followed by 3.  Alternating nissl fluoro will use 1 for initialization, followed by 3 for nissl slices, and then 4 for other fluoro slices.

The code in this repository will be matlab based.