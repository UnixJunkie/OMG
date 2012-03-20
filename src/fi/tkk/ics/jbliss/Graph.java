/* 
 * @(#)Graph.java
 *
 * Copyright 2007-2010 by Tommi Junttila.
 * Released under the GNU General Public License version 3.
 */

package fi.tkk.ics.jbliss;

import java.io.PrintStream;
import java.util.*;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;


public class Graph {
	/* Intermediate translator stuff for mapping
       bliss automorphisms back to jbliss automorphisms */
	protected Map<Integer,Integer> _bliss_map;
	protected Map<Integer,Integer> _bliss_map_inv;
	protected Object         _reporter_param;

	/* The internal JNI interface to true bliss */
	private native long create();
	private native void destroy(long true_bliss);
	protected native int _add_vertex(long true_bliss, int color);
	protected native void _add_edge(long true_bliss, int v1, int v2);
	protected native int[] _canonical_labeling(long true_bliss, Object r);
	
	private static Map<String, Integer> colorTable; 


	/**
	 * Find the canonical labeling and the automorphism group of the graph.
	 * If the argument reporter is non-null,
	 * then a generating set of automorphisms is reported by calling its
	 * {@link Reporter#report} method for each generator.
	 *
	 * @return           A canonical labeling permutation
	 */
	public Map<Integer, Integer> canonical_labeling(IAtomContainer atomContainer) {
		long bliss = create();
		assert bliss != 0;
		_bliss_map     = new TreeMap<Integer,Integer>();
		_bliss_map_inv = new TreeMap<Integer,Integer>();
		
		int vertexID = 0;	// used to give a unique ID to the vertices corresponding to bonds
		// Add each atom as a vertex with a color corresponding to the atom type
		for(IAtom atom: atomContainer.atoms()){
			addVertex(bliss, Integer.parseInt(atom.getID()), colorTable.get(atom.getSymbol()));
			vertexID++;
		}
		// Add each bond as a vertex connected to the adjacent atoms (turning a multi-graph to a simple one)
		for (IBond bond : atomContainer.bonds()) {
			addVertex(bliss, vertexID, 0);
			_add_edge(bliss, _bliss_map.get(Integer.parseInt(bond.getAtom(0).getID())), _bliss_map.get(vertexID));
			_add_edge(bliss, _bliss_map.get(Integer.parseInt(bond.getAtom(1).getID())), _bliss_map.get(vertexID));
			vertexID++;
		}
		int[] cf = _canonical_labeling(bliss, null);
		destroy(bliss);
		TreeMap<Integer,Integer> labeling = new TreeMap<Integer,Integer>();
		for(Map.Entry<Integer,Integer> e : _bliss_map.entrySet())
			labeling.put(e.getKey(), cf[e.getValue()]);
		_bliss_map = null;
		_bliss_map_inv = null;
		return labeling;
	}
	
	/**
	 * @param bliss
	 * @param vid
	 * @param color
	 */
	private void addVertex(long bliss, int vid, int color) {
		int bliss_vertex = _add_vertex(bliss, color);
		_bliss_map.put(vid, bliss_vertex);
		_bliss_map_inv.put(bliss_vertex, vid);
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
