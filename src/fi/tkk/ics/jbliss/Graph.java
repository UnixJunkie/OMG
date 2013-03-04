
package fi.tkk.ics.jbliss;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.omg.MolProcessor;
import org.omg.tools.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * To convert a multi-graph to a simple graph, we create N vertexes for each edge of order N.
 * These vertexes are connected to the actual vertexes at the two ends of the edge.
 * These vertexes have a different color than the actual ones.
 *  
 * @author mmajid
 *
 */
public class Graph {
	public int[] orbitRep;
	private int atomCount;
	private static Map<String, Integer> colorTable;
	private static final int bondColor = 1000;

	public static boolean blissFound; 


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

	
	public int[] canonize(final MolProcessor mol) {
		// create a bliss instance for calculating the canonical labeling
		final long bliss = create();
		assert bliss != 0;
		
		for (Atom a:mol.atoms){
			_add_vertex(bliss, colorTable.get(a.symbol));
		}
		for (int l=0; l<atomCount; l++)
			for (int r=l+1; r<atomCount; r++) {
				int bondOrder = mol.getBondOrder(l, r);
				for (int o=0; o<bondOrder; o++) {
					final int vid = _add_vertex(bliss, bondColor);
					_add_edge(bliss, l, vid);
					_add_edge(bliss, r, vid);
				}
			}
		final int[] cf = _canonical_labeling(bliss, null); // no need to report the automorphisms back
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
		// initialize the colorTable
		colorTable = new HashMap<String, Integer>();
		// TODO: read atom symbols from CDK
		colorTable.put("C", 1);
		colorTable.put("N", 2);
		colorTable.put("O", 3);
		colorTable.put("S", 4);
		colorTable.put("P", 5);
		colorTable.put("H", 6);
		
		/* Load the C++ library including the true bliss and
		 * the JNI interface code */   

		String osName = System.getProperty("os.name");
		String osArch = System.getProperty("os.arch");
		String blissName="";
		
		String dirName = "", fullPath="";
		try {
			dirName = new File(".").getCanonicalPath();
			blissFound = true;
			if (osName.contains("Windows")){
				blissName = "bliss."+osArch+".dll";
				fullPath = dirName + "\\" + blissName;
				if (osArch.equals("amd64")) {
					blissFound = false;
				} else {
					System.load(fullPath);
				}
			} else {
				blissName = "bliss"+osName+osArch+".so";
				fullPath = dirName+"/"+blissName;
				System.load(fullPath);
			}
		} catch (IOException e) {
			System.err.println("Could not get the current directory.");
			e.printStackTrace();
			System.exit(11);
		} catch (UnsatisfiedLinkError ule2){
			//ule2.printStackTrace();
			//System.loadLibrary(blissName);
			System.err.println("Problems loading the jbliss library for '"+osName+":"+osArch+"' at "+fullPath);
			System.err.println("Make sure the file is correctly named and placed. If the file is not available, then unfortunately you cannot run the program on this platform.");
			blissFound = false;
		} finally {
			if (!blissFound) { 
				System.out.println("Prescribed fragments and canonical augmentation are disabled because bliss is not available on your platform.");
				System.out.println("For these features to be available, you need Java 1.6.0_41, Java 1.7 or later (under windows use 'Java x86 32bits').");
			}
		}
	}
}
