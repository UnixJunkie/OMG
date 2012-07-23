package org.omg.tools;

import java.util.ArrayList;

public class Util {
	public static ArrayList<String> parseFormula (String formula){
		boolean success = false;
		ArrayList<String> atoms = new ArrayList<String>();
		String sym = "";
		while (formula.length() > 0) {
			success = false;
			sym += formula.charAt(0);
			formula = formula.substring(1); 
			if (!Atom.valenceTable.containsKey(sym)) continue;
			success = true;
			int count=0;
			char ch;
			while (formula.length() > 0) {
				ch = formula.charAt(0);
				if ('0' <= ch && ch <= '9') {
					count = count*10 + (ch-'0');
					formula = formula.substring(1); 
				} else {
					break;
				}
			}
			if (count == 0) count = 1;
			for (int n=0; n<count; n++) {
				atoms.add(sym);
			}
			sym = "";
		}
		if (!success) {
			System.err.println("Could not parse the sub-expression: "+sym);
			return null;
		}
		return atoms;
	}
	
}
