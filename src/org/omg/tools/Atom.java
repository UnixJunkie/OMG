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
	int ID;
	
	public Atom (String s, int id) {
		symbol = s;
		valenceList = valenceTable.get(symbol);
		maxValence = valenceList.get(0);
		ID = id;
	}

/*
	public static List<Atom> create (String s) {
		List<Integer> vlist = valenceTable.get(s);
		List<Atom> atomList = new ArrayList<Atom>(vlist.size());
		for (Integer v:vlist) {
			atomList.add(new Atom(s,v));
		}
		return atomList;
	}
*/
	static {
		// initialize the table
		valenceTable = new HashMap<>();
		// TODO: read atom symbols from CDK
		valenceTable.put("H", Arrays.asList(1));
		valenceTable.put("C", Arrays.asList(4));
		valenceTable.put("N", Arrays.asList(5,3));	// Remember to put the biggest number first!?
		valenceTable.put("O", Arrays.asList(2));
		valenceTable.put("S", Arrays.asList(6,4,2));
		valenceTable.put("P", Arrays.asList(5,3));
	}

}
