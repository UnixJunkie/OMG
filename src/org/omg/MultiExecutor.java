package org.omg;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class MultiExecutor implements MultiCoreExecutor {

	// int TODO 
	ThreadPoolExecutor executor[];
	private boolean parallelExecution[]; 
	LinkedBlockingQueue<Runnable> taskQueue[];
	private int executorCount;
	private final static int MAX_PARALLEL = 100;	// stop parallel tasks
	private final static int MIN_PARALLEL = 10;		// restart parallel tasks

	@SuppressWarnings("unchecked")
	MultiExecutor (int pDegree){
		executor = new ThreadPoolExecutor[pDegree];
		parallelExecution = new boolean[pDegree];
		taskQueue = new LinkedBlockingQueue[pDegree];	// cannot specify the specific <Runnable> type
		for (int i=0; i<pDegree; i++){
			taskQueue[i] = new LinkedBlockingQueue<Runnable>();
			executor[i] = new ThreadPoolExecutor(1,1, 0L, TimeUnit.MILLISECONDS, taskQueue[i]);
			parallelExecution[i] = true;
		}
		executorCount = pDegree;
	}
	
	@Override
	public void execute(Runnable command) {
		for (int i=0; i<executorCount; i++) {
			if ((parallelExecution[i] && taskQueue[i].size() < MAX_PARALLEL) || taskQueue[i].size() < MIN_PARALLEL) {
				parallelExecution[i] = true;
				executor[i].execute(command);
				return;
			} else {
				parallelExecution[i] = false;
			}
		}
		throw new RejectedExecutionException("Disabling parallelism due to overload. Should continue sequentially...");
	}

	@Override
	public void shutdown() {
		for (ThreadPoolExecutor ex : executor){
			ex.shutdown();
		}
	}

	@Override
	public int getQueueSize() {
		int i = 0;
		for (BlockingQueue queue : taskQueue) i += queue.size(); 
		return i;
	}

	@Override
	public void execute(Runnable command, boolean force) {
		if (force) {
			int i = (int) (Math.round(Math.random()) % executorCount);
			executor[i].execute(command);
		} else {
			execute(command);
		}
	}

	@Override
	public boolean busy() {
		throw new RuntimeException("Not yet implemented."); // TODO
	}

}
