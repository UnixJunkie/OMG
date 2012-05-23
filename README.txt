OMG - Open Molecule Generator
 
Copyright 2011 The Netherland Metabolomics Center Team License: GPL v3, see
doc/License-gpl-3.txt

1. Introduction
----------------------------------------------------------------------------------
You are currently reading the README file for the OMG Project. This project is
hosted under http://sourceforge.net/p/openmg

OMG is an open-source tool for the generation of chemical structures, implemented 
in the programming language Java(tm). The library is published under terms of the 
standard The GNU General Public License (GPL) v3. This has implications on what 
you can do with sources and binaries of the OMG library. 
For details, please refer to the file LICENSE, which is provided with this 
distribution.

PLEASE NOTE: OMG is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
A PARTICULAR PURPOSE.


2. System Requirements
----------------------------------------------------------------------------------
OMG.jar runs in the following OS;
- Ubuntu 32 bits
- Ubuntu 64 bits
- Mac OS X 64 bits


3. Using OMG tool
----------------------------------------------------------------------------------
In order to use the OMG tool in your program, you need to run the jar file via 
command line and provide some arguments. 

-ec:  elemental composition of the molecules to be generated.
-o:   SDF file where to store the molecules.  
-fr:  SDF file containing prescribed one or multiple substructures. In the case
	  of multiple substructures, they have to be non-overlapping. 

Here there some examples of how to run OMG using the command line: 

- Generating molecules
-- Generate molecules for the elemental composition C6H6
java -jar OMG.jar -ec C6H6

-- Generate molecules for the elemental composition C6H6 and store them in 
-- out_C6H6.sdf
java -jar OMG.jar -ec C6H6 -o out_C6H6.sdf

- Generating molecules with prescribed substructure(s)
-- Generate molecules for the elemental composition C2H5NO2 (glycine) using the 
-- prescribed substructure in fragment_CO2.sdf
java -jar OMG.jar -ec C2H5NO2 -fr fragment_CO2.sdf



4. Source Code
----------------------------------------------------------------------------------
You can download the source code at
http://sourceforge.net/p/openmg


5. Help
----------------------------------------------------------------------------------
If you need help don't hesitate to contact us (jpeironcely@gmail.com)

----------------------------------------------------------------------------------
Enjoy!  Comments and feedback are appreciated!

Julio E. Peironcely
jpeironcely@gmail.com
Leiden, Netherlands
http://juliopeironcely.com
