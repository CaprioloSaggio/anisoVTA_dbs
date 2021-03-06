#!/bin/bash

echo ""
echo "Splitting 4D image in 3D volumes"
fslsplit dti.nii

mv vol0000.nii.gz b0.nii.gz
mrconvert b0.nii.gz b0.nii
rm b0.nii.gz

echo ""
echo "Selecting a subsample of the 3D volumes, including a b0"
for i in {2..7}
do
VOLIN="vol000${i}.nii.gz";
VOLOUT="dti${i}.nii.gz";
NII="dti${i}.nii";
mv $VOLIN $VOLOUT;
mrconvert $VOLOUT $NII
done
rm vol*
rm *.gz

echo ""
echo "Coregistering b0 to anatomical image and saving the transformation"
flirt -in b0.nii -ref anat_t1.nii -out r_b0.nii -omat dti2anat.mat -dof 12
mrconvert r_b0.nii.gz r_b0.nii

echo ""
echo "Applying coregistration to all the considered diffusion-weighted volumes"
for i in {2..7}
do
NII="dti${i}.nii";
NIIOUT="r_${NII}";
flirt -in $NII -ref anat_t1.nii -out $NIIOUT -init dti2anat.mat -applyxfm
mrconvert "${NIIOUT}.gz" $NIIOUT 
done
rm *.gz

echo ""
echo "Merging the volumes"
fslmerge -t r_dti_biased r_b0.nii r_dti*
mrconvert r_dti_biased.nii.gz r_dti_biased.nii

echo ""
echo "Removing bias"
fslmaths r_dti_biased.nii -Tmean r_dti_m.nii
bet r_dti_m.nii r_dti_brain -m
mrconvert r_dti_brain_mask.nii.gz r_dti_brain_mask.mif
mrconvert r_dti_biased.nii r_dti_biased.mif
# dwibiascorrect ants r_dti_biased.mif r_dti.mif -mask r_dti_brain_mask.mif -fslgrad bvecs_short bvals_short -ants.b [100,3] -ants.c [1000,0.0] -ants.s 4
dwibiascorrect fsl r_dti_biased.mif r_dti.mif -fslgrad bvecs_short bvals_short -mask r_dti_brain_mask.mif
mrconvert r_dti.mif r_dti.nii

echo ""
echo "Computing diffusion tensor"
dtifit --data=r_dti.nii --mask=r_dti_brain_mask --bvecs=bvecs_short --bvals=bvals_short --save_tensor --out=dti
mrconvert dti_tensor.nii.gz dti_tensor.nii

rm *.gz

