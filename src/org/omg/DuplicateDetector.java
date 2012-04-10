package org.omg;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import fi.tkk.ics.jbliss.Graph;

public class DuplicateDetector {
	private HashMap<String, String> store = new HashMap<>();
	
	public static void main (String args[]){
		if (args.length < 1) {
			System.err.println("Input file name is missing. Quitting...");
			return;
		}
		try {
			DuplicateDetector dd = new DuplicateDetector();
			FileReader fileReader = new FileReader(args[0]);
			dd.detect(new BufferedReader(fileReader));
			fileReader.close();
		} catch (FileNotFoundException e) {
			System.err.println("File "+args[0]+" could not be opened.");
		} catch (IOException e) {
			System.err.println("Error reading from file.");
		}
	}

	/**
	 * @param args
	 * @param dd
	 * @throws IOException 
	 */
	private void detect(BufferedReader input) throws IOException {
			int dup = 0,total=0;
			String formula = input.readLine();
			System.out.println(formula);
			int atomCount = formula.length();
			char[] atoms = new char[atomCount];
			formula.getChars(0, atomCount, atoms, 0);
			String molStr;
			while ((molStr = input.readLine()) != null) {
				char[] bonds = new char[atomCount*atomCount];
				molStr.getChars(0, atomCount*atomCount, bonds, 0);
				String key = new Graph().canonize(atoms, bonds);
				total++;
				if (!store.containsKey(key)) {
					store.put(key, molStr);
				} else {
					System.out.println("Duplicate: "+molStr+" : "+store.get(key));
					System.out.println("Canonized: "+key);
					dup++;
				}
			}
			System.out.println("Duplicate count out of "+total+" for "+formula+" = "+dup);
	}
}
