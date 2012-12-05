package org.omg.tools;

import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.smsd.Substructure;

public class CDKUtil {
	final static SaturationChecker satCheck = new SaturationChecker();
	public static boolean acceptedByCDK(IAtomContainer acprotonate, int nH) throws CDKException {
		try{
			CDKAtomTypeMatcher typeMatcher = CDKAtomTypeMatcher.getInstance(acprotonate.getBuilder());
			for (IAtom atom : acprotonate.atoms()) {
				IAtomType type;
				type = typeMatcher.findMatchingAtomType(acprotonate, atom);
				if (type == null) return false;
				AtomTypeManipulator.configure(atom, type);
			}
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
			hAdder.addImplicitHydrogens(acprotonate);
			return (satCheck.isSaturated(acprotonate) && AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH);
		} catch (NullPointerException nle) {
			return false;
		}

	}

	public static boolean checkGood (IAtomContainer acprotonate, IAtomContainerSet goodlistquery) throws CDKException {
		if(goodlistquery != null){
			Substructure substructure = null;
			for(int i = 0; i < goodlistquery.getAtomContainerCount(); i++){
				substructure = new Substructure(goodlistquery.getAtomContainer(i),acprotonate,true, false, false);
				if (!substructure.isSubgraph()) {
					return false;
				}
			}
		}
		return true;
	}

	public static boolean checkBad (IAtomContainer acprotonate, IAtomContainerSet badlistquery) throws CDKException {
		// output only molecules that DON'T contain the substructures in badlistquery
		if(badlistquery != null){
			Substructure substructure = null;
			for(int i = 0; i < badlistquery.getAtomContainerCount(); i++){
				substructure = new Substructure(badlistquery.getAtomContainer(i),acprotonate,true, false, false);
				if (substructure.isSubgraph()) {
					return false;
				}
			}
		}
		return true;
	}
}
