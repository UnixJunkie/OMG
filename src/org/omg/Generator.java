package org.omg;

import java.util.ArrayList;

import org.openscience.cdk.exception.CDKException;

class Generator implements Runnable {
		/**
		 * 
		 */
		private final PMG pmg;
		MolHelper2 mol;
		
		public Generator(PMG pmg_FixedSizeExecutor, MolHelper2 mol) {
			super();
			pmg = pmg_FixedSizeExecutor;
			this.mol = mol;
		}

		@Override
		public void run() {
			try {
//				if (depth == atomCount) 
//					return;
				if (mol.isComplete(pmg.satCheck, pmg.nH)){
					if (mol.isConnected()) {
//						if (!molSet.add(mol.canString)) System.err.println("Duplicate");
						pmg.mol_counter.incrementAndGet();
						if(PMG.wFile){
							mol.writeTo(pmg.outFile, pmg.mol_counter.get());
						}
					}	
				}
				else{
					// get all possible ways to add one bond to the molecule
					ArrayList<MolHelper2> extMolList = mol.addOneBond();
					
					for (MolHelper2  molecule : extMolList) {
						if (!pmg.generateParallelTask(molecule)) {
							mol = molecule;
							run();	// continue sequentially
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
				}
			} catch (CloneNotSupportedException e){
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			}
			pmg.startedTasks.decrementAndGet();
			mol = null;
		}
		
	}