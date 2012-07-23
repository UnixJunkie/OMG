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
import org.omg.MolProcessor;
import org.omg.tools.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;


public class Graph {
	public int[] orbitRep;
	private int atomCount;

	public Graph() {
		throw new UnsupportedOperationException("Use the Graph constructor that takes atomCount as parameter.");
	}
	
	public Graph(int atomCount) {
		super();
		this.atomCount = atomCount;
		orbitRep = new int[atomCount];
	}
	/* The internal JNI interface to true bliss */
	private native long create();
	private native void destroy(long true_bliss);
	protected native int _add_vertex(long true_bliss, int color);
	protected native void _add_edge(long true_bliss, int v1, int v2);
	protected native int[] _canonical_labeling(long true_bliss, Object r);
	
	/**
	 * Upon receiving a new automorphism (permutation), we join the orbits accordingly.
	 * @param aut
	 */
	protected void _report(int[] aut)
	{
//		System.out.print("Aut: ");
		for (int i=0; i<atomCount; i++){
//			System.out.print(""+aut[i]+" ");
			if (orbitRep[orbitRep[i]] < orbitRep[orbitRep[aut[i]]])
				orbitRep[aut[i]] = orbitRep[orbitRep[i]];
			else
				orbitRep[i] = orbitRep[orbitRep[aut[i]]];
		}
//		System.out.print(";;\tOrbit:");
//		for (int i:orbitRep) System.out.print(" "+orbitRep[i]);
//		System.out.println();
	}

	private static Map<String, Integer> colorTable; 
	private static final int bondColor = 1000;

	
	public String canonize(char[] atoms, char[] bonds) {
		System.out.println("Whatever");
		// create a bliss instance for calculating the canonical labeling
		long bliss = create();
		assert bliss != 0;
		
		int atomCount=0;
		// Add each atom as a vertex with a color corresponding to the atom type
		for (char atom : atoms){
			_add_vertex(bliss, colorTable.get(""+atom));
			atomCount++;
		}
		int edges[][] = new int[atomCount][atomCount];
		int vid=0;
		// Add each bond as a vertex connected to the adjacent atoms (turning a multi-graph to a simple one)
		for (int atom0=0; atom0<atomCount; atom0++)
			for (int atom1=0; atom1<atomCount; atom1++) {
				char bond = bonds[atom0 * atomCount + atom1];
				int orderNumber = Integer.parseInt(""+bond);
				edges[atom0][atom1] = orderNumber;
				for (; 0<orderNumber; orderNumber--){
					vid = _add_vertex(bliss, bondColor);
					_add_edge(bliss, atom0, vid);
					_add_edge(bliss, atom1, vid);
				}
		}
		int[] cf = _canonical_labeling(bliss, null); 
		destroy(bliss);

		int canonicalBonds[] = new int[atomCount*atomCount];
		for (int left=0; left<atomCount; left++)
			for(int right=0; right<atomCount; right++){
				canonicalBonds[cf[left]*atomCount+cf[right]] = edges[left][right]+edges[right][left];
			}
		String str = "";
		for (int n:canonicalBonds) str+=n;
		return str;
	}

	
	
	/**
	 * Find the canonical labeling and the automorphism group of the graph.
	 *
	 * @return           A canonical labeling permutation
	 */
	public int[] canonize(IAtomContainer atomContainer, boolean report) {
		// initialize the orbit calculation
		atomCount = atomContainer.getAtomCount();
		orbitRep = new int[atomCount ];
		for (int i=0; i<atomCount; i++){
			orbitRep[i] = i;	// start with no symmetry
		}
		// create a bliss instance for calculating the canonical labeling
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
		int[] cf = _canonical_labeling(bliss, report?1:null); 
		destroy(bliss);

		return cf;
	}
	
	public int[] canonize(MolProcessor mol, boolean report) {
		for (int i=0; i<atomCount; i++){
			orbitRep[i] = i;	// start with no symmetry
		}
		// create a bliss instance for calculating the canonical labeling
		long bliss = create();
		assert bliss != 0;
		
		for (Atom a:mol.atoms){
			_add_vertex(bliss, colorTable.get(a.symbol));
		}
		for (int l=0; l<atomCount; l++)
			for (int r=l+1; r<atomCount; r++) 
				for (int o=0; o<mol.adjacency[l*atomCount+r]; o++) {
					int vid = _add_vertex(bliss, bondColor);
					_add_edge(bliss, l, vid);
					_add_edge(bliss, r, vid);				
				}
		int[] cf = _canonical_labeling(bliss, report?1:null); 
		destroy(bliss);

		return cf;
	}

		
	
	/**
	 * @param atomContainer
	 * @param perm
	 * @return
	 * @throws CloneNotSupportedException
	 */
	public static IAtomContainer relabel(IAtomContainer atomContainer, int[] perm)
			throws CloneNotSupportedException {
		IAtomContainer canonM_ext = (IAtomContainer) atomContainer.clone();
		int vid=0;
		for (IAtom a : canonM_ext.atoms()) {
			a.setID(""+(perm[vid++]));	
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
		colorTable.put("S", 4);
		colorTable.put("P", 5);
	}
}
