#!/bin/bash
echo Glycine fragment_CO2
java -jar OMG.jar -mf C2H5NO2 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_CO2.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Glycine_fragment_CO2.sdf

echo dCysteine fragment_CO2
java -jar OMG.jar -mf C3H7NO2S -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_CO2.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_dCysteine_fragment_CO2.sdf

echo Phenylalanine fragment_benzene
java -jar OMG.jar -mf C9H11NO2 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Phenylalanine_fragment_benzene.sdf

echo Phenylalanine fragment_benzene_CO2_bis
java -jar OMG.jar -mf C9H11NO2 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene_CO2_bis.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Phenylalanine_fragment_benzene_CO2_bis.sdf

echo Phenylalanine fragment_benzene_CO2_CN
java -jar OMG.jar -mf C9H11NO2 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene_CO2_CN.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Phenylalanine_fragment_benzene_CO2_CN.sdf

echo Phenylalanine fragment_benzene_C2NO2
java -jar OMG.jar -mf C9H11NO2 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene_C2NO2.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Phenylalanine_fragment_benzene_C2NO2.sdf

echo Malic Acid fragment_CO2
java -jar OMG.jar -mf C4H6O5 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_CO2.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_malicacid_fragment_CO2.sdf

echo Uric Acid fragment_C3N2O
java -jar OMG.jar -mf C5H4N4O3 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_C3N2O.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_uricacid_fragment_C3N2O.sdf

echo Uric Acid fragment_CN2O
java -jar OMG.jar -mf C5H4N4O3 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_CN2O.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_uricacid_fragment_CN2O.sdf

echo Phenyllactic fragment_benzene
java -jar OMG.jar -mf C9H10O3 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Phenyllactic_fragment_benzene.sdf

echo Phenyllactic fragment_benzene_CO2_bis
java -jar OMG.jar -mf C9H10O3 -method ca -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene_CO2_bis.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_Phenyllactic_fragment_benzene_CO2_bis.sdf

echo p-Cresol Sulfate fragment_sulfate
java -jar OMG.jar -mf C7H8O4S -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_sulfate.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_p-Cresol_fragment_sulfate.sdf

echo p-Cresol Sulfate fragment_benzene
java -jar OMG.jar -mf C7H8O4S -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_benzene.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_p-Cresol_fragment_benzene.sdf

echo cholic fragment_cholic_core_C4O2_bis2
java -Xms256m -Xmx1024m -XX:MaxPermSize=1024m -jar OMG.jar -mf C24H40O5 -fr /media/sf_PhD/data/data_current/Structgen/fragments/fragment_cholic_core_C4O2_bis2.sdf -o /media/sf_PhD/data/data_current/Structgen/results/out_cholic_fragment_cholic_core_C4O2_bis2.sdf