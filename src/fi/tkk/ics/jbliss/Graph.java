/* 
 * @(#)Graph.java
 *
 * Copyright 2007-2010 by Tommi Junttila.
 * Released under the GNU General Public License version 3.
 */

package fi.tkk.ics.jbliss;

import java.io.PrintStream;
import java.util.*;

import org.omg.MolHelper;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;


public class Graph {
	/* The internal JNI interface to true bliss */
	private native long create();
	private native void destroy(long true_bliss);
	protected native int _add_vertex(long true_bliss, int color);
	protected native void _add_edge(long true_bliss, int v1, int v2);
	protected native int[] _canonical_labeling(long true_bliss, Object r);
	
	private static Map<String, Integer> colorTable; 
	private static final int bondColor = 1000;

	/**
	 * Find the canonical labeling and the automorphism group of the graph.
	 *
	 * @return           A canonical labeling permutation
	 */
	public int[] canonize(IAtomContainer atomContainer) {
		long bliss = create();
		assert bliss != 0;
		
		// Add each atom as a vertex with a color corresponding to the atom type
		for (IAtom atom : atomContainer.atoms()){
			_add_vertex(bliss, colorTable.get(atom.getSymbol()));
		}
		// Add each bond as a vertex connected to the adjacent atoms (turning a multi-graph to a simple one)
		for (IBond bond : atomContainer.bonds()) {
			int atom0 = atomContainer.getAtomNumber(bond.getAtom(0));
			int atom1 = atomContainer.getAtomNumber(bond.getAtom(1));
			int orderNumber = MolHelper.orderNumber(bond.getOrder());
			int vid=0;
			for (; 0<orderNumber; orderNumber--){
				vid = _add_vertex(bliss, bondColor);
				_add_edge(bliss, atom0, vid);
				_add_edge(bliss, atom1, vid);
			}
			bond.setID(""+vid);
		}
		int[] cf = _canonical_labeling(bliss, null); 
		destroy(bliss);

		return cf;
	}
	
	/**
	 * @param atomContainer
	 * @param cf
	 * @return
	 * @throws CloneNotSupportedException
	 */
	public static IAtomContainer relabel(IAtomContainer atomContainer, int[] cf)
			throws CloneNotSupportedException {
		IAtomContainer canonM_ext = (IAtomContainer) atomContainer.clone();
		int vid=0;
		for (IAtom a : canonM_ext.atoms()) {
			a.setID(""+(cf[vid++]));	
		}
//		int n = atomContainer.getAtomCount();
//		for (vid=0; vid<n; vid++)
//			atomContainer.getAtom(cf[vid]).setID(""+(vid));
		return canonM_ext;
	}
	

	static {
		/* Load the C++ library including the true bliss and
		 * the JNI interface code */
		System.loadLibrary("jbliss");
		
		// initialize the colorTable
		colorTable = new HashMap<>();
		// TODO: read atom symbols from CDK
		colorTable.put("C", 1);
		colorTable.put("N", 2);
		colorTable.put("O", 3);
		colorTable.put("Br", 4);
	}
}
