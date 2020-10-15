mv *crd 4eyd_Hpp.crd
mv *top 4eyd_Hpp.top
python GAFF/acpype.py -x 4eyd_Hpp.crd -p 4eyd_Hpp.top
mv 4eyd_Hpp_GMX.top 4eyd.top
