#!/bin/bash
echo Glycine
java -jar OMG.jar -mf C2H5NO2 -o /media/sf_PhD/data/data_current/Structgen/results/out_Glycine.sdf

echo Acetyl-Glycine
java -jar OMG.jar -mf C4H7NO3 -o /media/sf_PhD/data/data_current/Structgen/results/out_Acetyl_Glycine.sdf

echo Glutamic Acid
java -jar OMG.jar -mf C5H9NO4 -o /media/sf_PhD/data/data_current/Structgen/results/out_Glutamic_Acid.sdf

echo Pyruvic Acid
java -jar OMG.jar -mf C3H4O3 -o /media/sf_PhD/data/data_current/Structgen/results/out_Pyruvic_Acid.sdf

echo Malic Acid
java -jar OMG.jar -mf C4H6O5 -o /media/sf_PhD/data/data_current/Structgen/results/out_Malic_Acid.sdf

echo Phosphoenolpyruvic Acid
java -Xms256m -Xmx1024m -XX:MaxPermSize=1024m -jar OMG.jar -mf C3H5O6P -o /media/sf_PhD/data/data_current/Structgen/results/out_Phosphoenolpyruvic_Acid.sdf

echo Creatinine
java -jar OMG.jar -mf C4H7N3O -o /media/sf_PhD/data/data_current/Structgen/results/out_Creatinine.sdf

echo Guanidinoacetic Acid
java -jar OMG.jar -mf C3H7N3O2 -o /media/sf_PhD/data/data_current/Structgen/results/out_Guanidinoacetic_Acid.sdf

echo Cytosine
java -jar OMG.jar -mf C4H5N3O -o /media/sf_PhD/data/data_current/Structgen/results/out_Cytosine.sdf

echo Histamine
java -jar OMG.jar -mf C5H9N3 -o /media/sf_PhD/data/data_current/Structgen/results/out_Histamine.sdf

echo D-Cysteine
java -jar OMG.jar -mf C3H7NO2S -o /media/sf_PhD/data/data_current/Structgen/results/out_DCysteine.sdf

echo P-cresol sulfate
java -Xms128m -Xmx1024m -XX:MaxPermSize=256m -jar OMG.jar -mf C7H8O4S -o /media/sf_PhD/data/data_current/Structgen/results/out_pcresol_sulfate.sdf





