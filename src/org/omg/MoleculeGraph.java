package org.omg;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import fi.tkk.ics.jbliss.Graph;
import fi.tkk.ics.jbliss.Graph.Vertex;


/**
 * Changing a multi-graph to a single graph (that can be used by existing graph isomorphism checkers like: bliss or nauty
 * To do so, each "edge" in the multi-graph is converted to a vertex in the normal graph. All edges will have the same color
 * which will be different from the colors of other vertices.
 * 
 * @author mmajid
 *
 */
public class MoleculeGraph {
	IAtomContainer acontainer;
	Graph<Integer> atomGraph;
	
	public MoleculeGraph() {
	}


	public MoleculeGraph(IAtomContainer acontainer, Graph<Integer> aGraph) {
		super();
		this.acontainer = acontainer;
		this.atomGraph = aGraph;
	}


	public int initialize (String formula) throws CloneNotSupportedException{
		atomGraph = new Graph<>();
		acontainer = MolecularFormulaManipulator.getAtomContainer(
				MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));

		// Count and remove the Hydrogens
		int nH = 0;
		List<IAtom> listcont = new ArrayList<IAtom>();
		for(IAtom atom: acontainer.atoms()){
			String symbol = atom.getSymbol();
			if(symbol.equals("H")){
				nH++;
				listcont.add(atom);
			}
		}
		for(IAtom atom: listcont){
			acontainer.removeAtom(atom);
		}
		
		// set up other atoms and build the corresponding graph
		int colorCounter = 0;	// we will start from color 1 (color 0 is reserved for bonds)
		String prevSymbol = "";
		int na = acontainer.getAtomCount();
		for (int i=0; i<na; i++){
			IAtom atom = acontainer.getAtom(i);
			String symbol = atom.getSymbol();
			atom.setID(""+i);
			atom.setFlag(1, false);
			// Change the color for different atom types
			if (!prevSymbol.equals(symbol)) {	// TODO: Assuming that atoms are sorted by their symbol
				prevSymbol = symbol;
				colorCounter ++;
			}
			// Add a vertex for each atom (initially there are no bonds between the atoms)
			atomGraph.add_vertex(i, colorCounter);
		}
		return nH;
	}
	
	public int initialize (String formula, String fragments) throws CloneNotSupportedException, CDKException, FileNotFoundException{
		// generate the atom as in the above constructor
		int nH = this.initialize(formula);

		// Assuming that fragments is not null, add the fragments to the atom structure and canonicalize
		InputStream ins = new BufferedInputStream(new FileInputStream(fragments));
		MDLV2000Reader reader = new MDLV2000Reader(ins);
		ChemFile fileContents = (ChemFile)reader.read(new ChemFile());		        
		IChemSequence sequence = fileContents.getChemSequence(0);
		
		for (int i=0; i<sequence.getChemModelCount(); i++) {
			IAtomContainer frag = sequence.getChemModel(i).getMoleculeSet().getAtomContainer(0);
			
			try {
				acontainer = MolManipulator.buildFromFragment(acontainer,frag);
			} catch (CloneNotSupportedException e1) {
				e1.printStackTrace();
			}
		}
		// TODO: Add the edges corresponding to the fragment
		
		acontainer = MolManipulator.getcanonical(acontainer);
		return nH;
	}
		
	/**
	 * Checks if a molecule is complete, including checking the number of hydrogens against nH, and saturation.
	 * @param satCheck
	 * @param nH
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	public boolean isComplete(SaturationChecker satCheck, int nH) throws CloneNotSupportedException, CDKException{
		IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();

		for (IAtom atom : acprotonate.atoms()) {
			IAtomType type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acprotonate, atom);
			
			AtomTypeManipulator.configure(atom, type);
		}
		CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
		hAdder.addImplicitHydrogens(acprotonate);
		return (satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH));
	}
	
	
	/**
	 * Checks if all atoms in the molecule are connected.
	 * @return
	 */
	public boolean isConnected() {
		return (ConnectivityChecker.partitionIntoMolecules(acontainer).getAtomContainerCount() == 1);
	}
	
	
	/**
	 * Computes the canonical IAtomContainer based on the relabeling provided.
	 * @param labeling
	 * @param m_ext
	 * @return
	 * @throws CloneNotSupportedException 
	 */
	public static IAtomContainer relabel(Map<Integer, Integer> labeling, IAtomContainer m_ext) throws CloneNotSupportedException {
		assert labeling != null;

		IAtomContainer canonM_ext = (IAtomContainer) m_ext.clone();
		for (IAtom a : canonM_ext.atoms()) {
			a.setID(""+labeling.get(Integer.parseInt(a.getID())));
		}
		
		return canonM_ext;
	}

	
	/**
	 * It checks all possible ways to add one bond to a molecule, while considering the canonical augmentation of the corresponding graph.
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	ArrayList<MoleculeGraph> addOneBond() throws CloneNotSupportedException, CDKException{
    	Set<Graph<Integer>> visited = new TreeSet<Graph<Integer>>();
		ArrayList<MoleculeGraph> extMolList = new ArrayList<MoleculeGraph>();
		
		ArrayList<int[]> extBondlist = MolManipulator.extendMol(acontainer);
		
		Integer vertexID = atomGraph.nof_vertices();
		for(int[] bond : extBondlist){
			Graph<Integer> extGraph = atomGraph; //.copy();
			
			// add the bond as a new vertex and the corresponding edges to the graph representation
			extGraph.add_vertex(vertexID , 0);
			extGraph.add_edge(vertexID, Integer.parseInt(acontainer.getAtom(bond[0]).getID()));
			extGraph.add_edge(vertexID, Integer.parseInt(acontainer.getAtom(bond[1]).getID()));
			
			Map<Integer, Integer> canonical_labeling = extGraph.canonical_labeling();
			Graph<Integer> canonicalExtGraph = extGraph.relabel(canonical_labeling);
			extGraph.del_vertex(vertexID);

			if (visited.add(canonicalExtGraph) == false) continue;
			
			if (isCanonicallyAugmented(canonicalExtGraph)){
				IAtomContainer relabeledMolecule = relabel(canonical_labeling, acontainer);
				incBond(bond[0], bond[1], relabeledMolecule);
				extMolList.add(new MoleculeGraph (relabeledMolecule, canonicalExtGraph)); 
			}
		}
		return extMolList;
	}

	/**
	 * Checks using bliss 
	 * @param cGraph
	 * @return
	 */
	private boolean isCanonicallyAugmented(Graph<Integer> cGraph) {
		Graph<Integer> checkGraph = cGraph.copy();
		int bondCount = acontainer.getBondCount();
		if (bondCount==0) return true;
		int d=0;	// Since the bond vertices always have 0 as color, the first vertex in a canonical graph corresponds to a bond
		for (d=checkGraph.nof_vertices()-1; checkGraph.getVertexColor(d) != 0; d--);	// get the (vertex corresponding to the) last bond
		checkGraph.del_vertex(d);
		checkGraph = checkGraph.relabel(checkGraph.canonical_labeling());
		return checkGraph.compareTo(atomGraph) == 0;
	}
	
	private void incBond(int leftAtom, int rightAtom, IAtomContainer m_ext) {
		//Only one object bond is used, and recycled at every time we use it
		IBond bondAdd = m_ext.getBond(m_ext.getAtom(leftAtom), m_ext.getAtom(rightAtom));
		if(bondAdd == null){					
			m_ext.addBond(leftAtom, rightAtom, IBond.Order.SINGLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.SINGLE){
			m_ext.getBond(m_ext.getAtom(leftAtom), m_ext.getAtom(rightAtom)).setOrder(IBond.Order.DOUBLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.DOUBLE){
			m_ext.getBond(m_ext.getAtom(leftAtom), m_ext.getAtom(rightAtom)).setOrder(IBond.Order.TRIPLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.TRIPLE){
			m_ext.getBond(m_ext.getAtom(leftAtom), m_ext.getAtom(rightAtom)).setOrder(IBond.Order.QUADRUPLE);					
		}
	}

}
