package org.omg;

import java.util.ArrayList;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;

class Generator implements Runnable {
		/**
		 * 
		 */
		private final PMG pmg;
		MolHelper2 mol;
		
		public Generator(PMG pmgExecutor, MolHelper2 mol) {
			super();
			pmg = pmgExecutor;
			this.mol = mol;
		}

		@Override
		public void run() {
			try {
//				if (depth == atomCount) 
//					return;
				if (mol.isComplete(/*pmg.satCheck,*/ pmg.nH)){
					if (mol.isConnected()) {
//						if (!molSet.add(mol.canString)) System.err.println("Duplicate");
						long count = pmg.mol_counter.incrementAndGet();
						if(PMG.wFile){
							mol.writeMol(pmg.outFile, count);
						}
					}	
				}
//				else{
					// get all possible ways to add one bond to the molecule
					ArrayList<MolHelper2> extMolList = mol.addOneBond();
					
					for (MolHelper2  molecule : extMolList) {
						if (pmg.executor.getQueueSize() > pmg.executorCount*2) {
							mol = molecule;
							run();	// continue sequentially
						} else {
							pmg.generateParallelTask(molecule, true);
						}
					}

//					// recursively process all extended molecules
//					if (pmg.parallelExecution)
//					{	// make a parallel call
//						if (pmg.taskQueue.size() > PMG_FixedSizeExecutor.executorCount*100) {
//							pmg.parallelExecution = false;
//							System.out.println("Disabling parallelism at task queue of "+pmg.taskQueue.size()+" with active count: "+pmg.executor.getActiveCount());
//						}
//						for (MolHelper2  molecule : extMolList) {
//							pmg.generateTaskMol(molecule);
//						}
//					} else
//					{ // do a recursive call on the same thread 
//						if (pmg.taskQueue.size() < PMG_FixedSizeExecutor.executorCount*2){
//							pmg.parallelExecution = true;
//							System.out.println("Enabling parallelism at task queue of "+pmg.taskQueue.size()+" with active count: "+pmg.executor.getActiveCount());
//						}
//						for (MolHelper2  molecule : extMolList) {
//							mol = molecule;
//							pmg.startedTasks.incrementAndGet();
//							run();
//						}
//					}
//				}
			} catch (CloneNotSupportedException e){
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			}
			mol = null;
		}
		
	}