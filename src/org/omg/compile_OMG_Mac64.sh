#!/bin/bash
cp nautygetcan.c nautygetcan_backup.c
javac OMGJNI.java
cd ../..

javah -jni org.omg.OMGJNI 

mv org/omg/org_omg_OMGJNI.h org/omg/org_omg_OMGJNI_old.h
mv org_omg_OMGJNI.h org/omg/org_omg_OMGJNI.h

cd org
cd omg

cc -c -I/System/Library/Frameworks/JavaVM.framework/Headers nautygetcan.c nauty.c nautil.c naututil.c naugraph.c rng.c 

cc -dynamiclib -o libnautygetcanMac64.jnilib nautygetcan.o nauty.o nautil.o naututil.o naugraph.c rng.o 

