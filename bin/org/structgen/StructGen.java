package org.structgen;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLWriter;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 * The main class ..
 * 
 * @author julio
 */
public class StructGen{
	/**Output File containing the list of graph. */
	BufferedWriter outFile;
	BufferedWriter outFile_intermediates;
	int graph_counter;
	static int id_counter;
	private int nH;
	private SaturationChecker satCheck;
	public static void main(String[] args) throws IOException{
		
		StructGen gen = new StructGen();
		String typerun = "mol";
		String formula = "C4H10";
		String fragments = null;
		int nvertex = 0;
		String out = "default_out.sdf";
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-t")){
				typerun = args[i+1];
			}
			else if(args[i].equals("-fr")){
				fragments = args[i+1];
			}
			else if(args[i].equals("-mf")){
				formula = args[i+1];
			}
			else if(args[i].equals("-o")){
				out = args[i+1];
			}
			else if(args[i].equals("-n")){
				nvertex = Integer.parseInt(args[i+1]);
			}
		}
		if(typerun.equals("graph")){
			
			gen.initializeGraph(nvertex,out);
		}
		else if(typerun.equals("multi")){
				
				gen.initializeMultiGraph(formula,out);
		}
		else if(typerun.equals("mol")){
				try {
					gen.initializeMolecule(formula,fragments, out);
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		
	}
	public StructGen(){
		
	}
	private void initializeGraph(int nvertex, String output) {
		long before = System.currentTimeMillis();
		int n_vertex = nvertex;
		graph_counter = 0;
		try {
			outFile = new BufferedWriter(new FileWriter(output));
			outFile.write(n_vertex + " vertices\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Graph g = new Graph(n_vertex);
		generateGraph(g);
		
		long after = System.currentTimeMillis();
		try {
			outFile.write("Duration: " + (after - before) + " miliseconds\n");
			outFile.write("Number of graphs: " + graph_counter + "\n");
			outFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}
	/**
	 * Start generating ..
	 * 
	 * @param graph The Graph
	 * @throws IOException
	 */
	public void generateGraph(Graph graph){
		if((GraphManipulator.isSaturated(graph))&&(graph.isComplete())){
			//We return the graph without changes, so we can finish the recursion

			graph.setComplete(true);//??? TODO change
			graph_counter++;
			try {
				outFile.write(graph_counter + "  " +graph.getCanonstr());
			} catch (IOException e) {
				e.printStackTrace();
			}

			return;
		}
		else{
			Graph canong = GraphManipulator.getcanonical(graph);
			/*generate all possible edges*/			
			ArrayList<Graph> extgraphlist = GraphManipulator.extendGraph(graph);
			
			for(Graph g_ext : extgraphlist){
				/*check constraints*/

				/*Graph canong_ext is canon(G')*/
				Graph canong_ext = GraphManipulator.getcanonical(g_ext);

				int lastV[] = canong_ext.lastE();
				
				/*Graph g_ext_e is G'-e'*/
				Graph g_ext_e = g_ext.clone();
				/*we have to map the edge in canon(G') to the labeling we have tin G'
				if we want to remove the right edge*/
				g_ext_e.removeEdge(canong_ext.getLab(lastV[0]), canong_ext.getLab(lastV[1]));

				/*Graph canong_ext_e is canon(G'-e')*/
				Graph canong_ext_e = GraphManipulator.getcanonical(g_ext_e);

				if(GraphManipulator.aresame(canong, canong_ext_e)||(graph.geteCount()==0)){

					g_ext.setComplete(true);

					graph_counter++;
					try {
						outFile.write(graph_counter + "  \n" + g_ext.toString());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					generateGraph(g_ext);
				}				
			}			
			return;
		}		
	}
	private void initializeMultiGraph(String formula, String output) {
		long before = System.currentTimeMillis();
		int n_vertex = 2;
		int[] lab = new int[n_vertex]; 
		int[] ptn = new int[n_vertex];
		
		lab[0]=0;lab[1]=1;lab[2]=2;
		ptn[0]=1;ptn[1]=0;ptn[2]=0;
		
		graph_counter = 0;
		try {
			outFile = new BufferedWriter(new FileWriter(output));
			outFile.write(n_vertex + " vertices\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		MultiGraph g = new MultiGraph(n_vertex, lab, ptn);
		generateMultiGraph(g);
		
		long after = System.currentTimeMillis();
		try {
			outFile.write("Duration: " + (after - before) + " miliseconds\n");
			outFile.write("Number of graphs: " + graph_counter + "\n");
			outFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public void generateMultiGraph(MultiGraph graph){
		if((MultiGraphManipulator.isSaturated(graph))&&(graph.isComplete())){
			//We return the graph without changes, so we can finish the recursion

			graph.setComplete(true);//??? TODO change
			graph_counter++;
			return;
		}
		else{
			MultiGraph canong = MultiGraphManipulator.getcanonical(graph);
			System.out.println("graph_matrix    " + graph.matrixToString());

			/*generate all possible edges*/			
			ArrayList<MultiGraph> extgraphlist = MultiGraphManipulator.extendMultiGraph(graph);
			for(MultiGraph mg : extgraphlist){
				System.out.println("mg " + mg.matrixToString());
			}
			for(MultiGraph g_ext : extgraphlist){

				/*Graph canong_ext is canon(G')*/				
				MultiGraph canong_ext = MultiGraphManipulator.getcanonical(g_ext);
				int lastV[] = new int[2];
				lastV[0] = canong_ext.getLastV1();
				lastV[1] = canong_ext.getLastV2();
				/*Graph g_ext_e is G'-e'*/
				MultiGraph g_ext_e = g_ext.clone();
				/*we have to map the edge in canon(G') to the labeling we have tin G'
				if we want to remove the right edge*/
				g_ext_e.removeEdge(canong_ext.getLab(lastV[0]), canong_ext.getLab(lastV[1]));

				/*Graph canong_ext_e is canon(G'-e')*/
				MultiGraph canong_ext_e = MultiGraphManipulator.getcanonical(g_ext_e);
				if(MultiGraphManipulator.aresame(canong, canong_ext_e)||(graph.geteCount()==0)){
					g_ext.setComplete(true);
					graph_counter++;
					try {
						outFile.write(g_ext.toString());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					generateMultiGraph(g_ext);
				}				
			}			
			return;
		}		
	}

	private void initializeMolecule(String formula, String fragments, String output) throws CDKException, FileNotFoundException {
		long before = System.currentTimeMillis();

		IAtomContainer acontainer = MolecularFormulaManipulator.getAtomContainer(MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));
		
		nH = 0;
		List<IAtom> listcont = new ArrayList<IAtom>();
		int atom_counter = 0;
		int id_counter = 0;
		for(IAtom atom: acontainer.atoms()){
			if(atom.getSymbol().equals("H")){
				nH++;
				listcont.add(atom);
			}
			else{
				atom.setID(""+atom_counter);
				atom_counter++;
			}
		}
		for(IAtom atom: listcont){
			acontainer.removeAtom(atom);
		}
		if(fragments != null){
			InputStream ins = new BufferedInputStream(new FileInputStream(fragments));
			MDLV2000Reader reader = new MDLV2000Reader(ins);
			ChemFile fileContents = (ChemFile)reader.read(new ChemFile());		        

			IMoleculeSet som = fileContents.getChemSequence(0).getChemModel(0).getMoleculeSet();
			IMolecule frag = som.getMolecule(0);

			try {
				acontainer = MolManipulator.builFromFragment(acontainer,frag);
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		       
//		       SmilesGenerator sg = new SmilesGenerator();
//				String smiles2 = sg.createSMILES(new Molecule(acontainer));
				acontainer = MolManipulator.getcanonical(acontainer);
		}		
		
		graph_counter = 0;
		try {
			outFile = new BufferedWriter(new FileWriter(output));
			outFile_intermediates = new BufferedWriter(new FileWriter("output_intermediates2.txt"));
			outFile_intermediates.write( "smiles"  + "\t" + "Id" + "\t" + "FatherId" + "\t" +
					"Father" + "\t" + "FatherCan" + "\t" + "FatherRed" + "FatherRedCan" + "\t" +
					"Child" + "ChildCan" + "\t" + "ChildRed" + "ChildRedCan" + "\t" +"\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			satCheck = new SaturationChecker();
			generateMol(acontainer);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		try {
			outFile.close();
			outFile_intermediates.close();
			System.out.println("molecules " + graph_counter);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long after = System.currentTimeMillis();
	
			System.out.println("Duration: " + (after - before) + " miliseconds\n");
	}

	
	
	private void generateMol(IAtomContainer acontainer) throws Exception {
		/*We check that the atoms are connected in the same molecules
		 * this is, they are not separated fragments*/
		IMoleculeSet fragments = ConnectivityChecker.partitionIntoMolecules(acontainer);
		/*We add hydrogens in order to check if the molecule is saturated.
		 * We will accept the molecule if the number of hydrogens necessary to saturate 
		 * is the same as the hydrogens in the original formula*/
		IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder());    	
    	for (IAtom atom : acprotonate.atoms()) {
    		IAtomType type = matcher.findMatchingAtomType(acprotonate, atom);
    		AtomTypeManipulator.configure(atom, type);
    	}
        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
    	hAdder.addImplicitHydrogens(acprotonate);
		if(satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH)){
		
				if(fragments.getMoleculeCount() == 1){
//				SmilesGenerator sg = new SmilesGenerator();
//				String smiles2 = sg.createSMILES(new Molecule(acprotonate));
				StringWriter writer = new StringWriter();
				MDLWriter mdlWriter = new MDLWriter(writer);
				mdlWriter.write(new Molecule(acprotonate));
				outFile.write(writer.toString()+ "> <Id> \n "+ acprotonate.getID() + "\n \n" +"$$$$\n");
				graph_counter++;

			}
			return;
		}
		else{
//			System.out.println("atoms acon1");
//			for(IAtom a : acontainer.atoms()){
//				System.out.println(a);
//			}
//			IAtomContainer canonM = MolManipulator.getcanonical(acontainer);
			IAtomContainer canonM = MolManipulator.getcanonical(MolManipulator.getreduction(acontainer));
//			System.out.println("atoms acon2");
//			for(IAtom a : canonM.atoms()){
//				System.out.println(a);
//			}
			
			ArrayList<IAtomContainer> extMollist = MolManipulator.extendMol(acontainer);
			ArrayList<IAtomContainer> validMollist = new ArrayList<IAtomContainer>();
			System.out.println("validMollist.size(): "+extMollist.size());
			for(IAtomContainer m_ext : extMollist){
				/*Graph canong_ext is canon(G')*/
//				IAtomContainer canonM_ext = MolManipulator.getcanonical(m_ext);

				for(IAtom a1 : m_ext.atoms()){
					System.out.println(" a0 "+ a1.getID()+";"+m_ext.getConnectedAtomsList(a1).size());
				}
				System.out.println("m_ext red  " +  MolManipulator.matrixToString(MolManipulator.getMatrix(MolManipulator.getreduction(m_ext))));

				IAtomContainer canonM_ext = MolManipulator.getcanonical(MolManipulator.getreduction(m_ext));
				
				if(m_ext.getBondCount()==9)
				{
					SmilesGenerator sg2 = new SmilesGenerator();
					String smiles3 = sg2.createSMILES(new Molecule((IAtomContainer) m_ext.clone()));
					System.out.println(" finished molecule: " + smiles3 + " "+ + m_ext.getAtomCount());
				}
//				
//				System.out.println("2lalalal " + MolManipulator.matrixToString(MolManipulator.getMatrix(m_ext)));
//				System.out.println("2lalalal " + MolManipulator.matrixToString(MolManipulator.getMatrix(canonM_ext)));
				
//				break;
				System.out.println("canonM_ext " + MolManipulator.matrixToString(MolManipulator.getMatrix(canonM_ext)));
				String lastBondID[] = new String[2];
				lastBondID[0] = canonM_ext.getBond(canonM_ext.getBondCount()-1).getAtom(0).getID();
				lastBondID[1] = canonM_ext.getBond(canonM_ext.getBondCount()-1).getAtom(1).getID();
				
				/*Graph g_ext_e is G'-e'*/

				System.out.println("abans red " + MolManipulator.matrixToString(MolManipulator.getMatrix(m_ext)));
				IAtomContainer m_ext_e = (IAtomContainer) MolManipulator.getreduction(m_ext).clone();
				System.out.println("despres red " + MolManipulator.matrixToString(MolManipulator.getMatrix(m_ext_e)));
				System.out.println("n_bondsA " + m_ext_e.getBondCount());
				int[] lastBond = new int[2];
				for(IAtom atom : m_ext_e.atoms()){
					if(atom.getID().equals(lastBondID[0]))
						lastBond[0]=m_ext_e.getAtomNumber(atom);
					if(atom.getID().equals(lastBondID[1]))
						lastBond[1] = m_ext_e.getAtomNumber(atom);
				}
				IBond bondremove = m_ext_e.getBond(m_ext_e.getAtom(lastBond[0]),m_ext_e.getAtom(lastBond[1]));
				
				
				if(bondremove.getOrder() == IBond.Order.SINGLE){
					m_ext_e.removeBond(bondremove);
					System.out.println("aquiii");
				}
				else if(bondremove.getOrder() == IBond.Order.DOUBLE){
					System.out.println("aquiii2");
					bondremove.setOrder(IBond.Order.SINGLE);
				}
				else if(bondremove.getOrder() == IBond.Order.TRIPLE){
					bondremove.setOrder(IBond.Order.DOUBLE);
				}
				else if(bondremove.getOrder() == IBond.Order.QUADRUPLE){
					bondremove.setOrder(IBond.Order.TRIPLE);
				}
				System.out.println("n_bondsB " + m_ext_e.getBondCount());
				/*Graph canong_ext_e is canon(G'-e')*/
				IAtomContainer canonM_ext_e = MolManipulator.getcanonical(m_ext_e);
				System.out.println("n_bondsC " + canonM_ext_e.getBondCount());
//				if(MolManipulator.aresame(canonM, canonM_ext_e)||(acontainer.getBondCount()==0)){
//					generateMol(m_ext);					
//				}	

				System.out.println("canonM        " +  MolManipulator.matrixToString(MolManipulator.getMatrix(canonM)));
				System.out.println("canonM_ext_e  " +  MolManipulator.matrixToString(MolManipulator.getMatrix(canonM_ext_e)));
//				break;
				SmilesGenerator sg = new SmilesGenerator();
				String smiles2 = sg.createSMILES(new Molecule(m_ext));
//				System.out.println(" hola " + smiles2 + " "+ m_ext.getAtomCount());
				id_counter++;
				m_ext.setID(Integer.toString(id_counter));
//				outFile_intermediates.write( smiles2 + "\n");
				
				outFile_intermediates.write(smiles2  + "\t" + m_ext.getID() + "\t" + acontainer.getID() + "\t" + 
						MolManipulator.matrixToString(MolManipulator.getMatrix(acontainer))+ "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(MolManipulator.getcanonical(acontainer)))+ "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(MolManipulator.getreduction(acontainer)))+ "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(canonM)) + "\t" +	
						
						MolManipulator.matrixToString(MolManipulator.getMatrix(m_ext))+  "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(MolManipulator.getcanonical(m_ext)))+ "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(MolManipulator.getreduction(m_ext)))+ "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(canonM_ext)) + "\t" +
						MolManipulator.matrixToString(MolManipulator.getMatrix(canonM_ext_e)) +
						"\n");
				
				if(MolManipulator.aresame(canonM, canonM_ext_e)||(acontainer.getBondCount()==0)){
					System.out.println("valid  ");
//					SmilesGenerator sg = new SmilesGenerator();
//					String smiles2 = sg.createSMILES(new Molecule(m_ext));
////					System.out.println(" hola " + smiles2 + " "+ m_ext.getAtomCount());
//					id_counter++;
//					m_ext.setID(Integer.toString(id_counter));
////					outFile_intermediates.write( smiles2 + "\n");
//					
//					outFile_intermediates.write(m_ext.getID() + "\t" + acontainer.getID() + "\t" + smiles2  + "\t" +
//							MolManipulator.matrixToString(MolManipulator.getMatrix(canonM)) + "\t" +	MolManipulator.matrixToString(MolManipulator.getMatrix(canonM_ext_e))+ "\n");
					m_ext.setProperty("Id", m_ext.getID());
					validMollist.add(m_ext);	
				}else
					System.out.println("Not valid  ");
				
			}	
			for(IAtomContainer to_extend : validMollist){
				generateMol(to_extend);
			}
			return;				
		}
		
		
		
	}

}
