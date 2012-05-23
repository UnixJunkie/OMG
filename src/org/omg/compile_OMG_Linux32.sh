#!/bin/bash
cp nautygetcan.c nautygetcan_backup.c
javac OMGJNI.java
cd ../..

javah -jni org.omg.OMGJNI 

mv org/omg/org_omg_OMGJNI.h org/omg/org_omg_OMGJNI_old.h
mv org_omg_OMGJNI.h org/omg/org_omg_OMGJNI.h

cd org
cd omg

gcc -o libnautygetcanLinux32.so -shared -I/usr/lib/jvm/java-6-openjdk/include -I/usr/lib/jvm/java-6-openjdk/include/linux nautygetcan.c nauty.c nautil.c naututil.c naugraph.c rng.c 

