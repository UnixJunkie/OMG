package org.omg;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.management.RuntimeErrorException;


public class SingleExecutor implements MultiCoreExecutor {
	private static final int MAX_PARALLEL = 100;
	private static final int MIN_PARALLEL = 10;
	ThreadPoolExecutor executor;
	LinkedBlockingQueue<Runnable> taskQueue;
	private boolean parallelExecution = true;

	SingleExecutor (int pDegree) {
		taskQueue = new LinkedBlockingQueue<Runnable>();
		executor = new ThreadPoolExecutor(pDegree, pDegree, 0L, TimeUnit.MILLISECONDS, taskQueue);
	}
	
	@Override
	public void execute(Runnable task) {
		if ((parallelExecution && taskQueue.size() < MAX_PARALLEL) || taskQueue.size() < MIN_PARALLEL) {
			parallelExecution = true;
			executor.execute(task);
		} else {
			parallelExecution  = false;
			throw new RejectedExecutionException("Disabling parallelism due to overload. Should continue sequentially...");
		}
	}
	
	public void shutdown() {
		executor.shutdown();
	}

	@Override
	public int getQueueSize() {
		return taskQueue.size();
	}

	@Override
	public void execute(Runnable command, boolean force) {
		if (force) {
			executor.execute(command);
		} else {
			this.execute(command);
		}
	}

	@Override
	public boolean busy() {
		throw new RuntimeException("Not yet implemented."); // TODO
	}

}
