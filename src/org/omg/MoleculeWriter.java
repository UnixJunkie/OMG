package org.omg;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringWriter;

import org.omg.tools.Atom;


public class MoleculeWriter implements Runnable{
	final BufferedWriter outFile;
	final String sdfStr;
	
	MoleculeWriter(String str, BufferedWriter outFile) {
		this.outFile = outFile;
		sdfStr = str;
	}



	/**
	 * Writes a molecule to a file in SDF format. This is more or less equivalent to writeTo method but it does not depend on CDK. 
	 * The good is that it does not have the restrictions of having a "valid" atom container and is expected to work faster.
	 * The bad side is that it has less features and cannot handle all complicated chemical structures. 
	 * @param outFile
	 * @param mol_counter
	 * @param canString 
	 * @param adjacency 
	 */
	public void run() {
		try {
			outFile.write(sdfStr);
		} catch (IOException e) {
			System.err.println("Could not write molecule to output file.");
			e.printStackTrace();
		} 
	}
	
}
