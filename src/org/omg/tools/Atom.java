package org.omg.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class Atom {
	static HashMap<String,  List<Integer>> valenceTable;
	public final String symbol;
	public final int maxValence;
	public final List<Integer> valenceList;
	public boolean flag = false;	// can be used for different purposes, e.g., to mark the atom being in a fragment
	
	public Atom (String s, int id) {
		symbol = s;
		valenceList = valenceTable.get(symbol);
		if (valenceList == null) {
			System.err.println("The atom type "+symbol+" is not supported.");
			System.exit(200);
		}
		maxValence = valenceList.get(0);
	}

	static {
		// initialize the table
		valenceTable = new HashMap<String,  List<Integer>>();
		// TODO: read atom symbols from CDK?
		valenceTable.put("H", Arrays.asList(1));
		valenceTable.put("C", Arrays.asList(4));
		valenceTable.put("N", Arrays.asList(5,3));	// Remember to put the biggest number first!?
		valenceTable.put("O", Arrays.asList(2));
		valenceTable.put("S", Arrays.asList(6,4,2));
		valenceTable.put("P", Arrays.asList(5,3));
	}

}
