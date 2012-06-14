package org.omg;

import java.util.concurrent.Executor;

public interface MultiCoreExecutor extends Executor {
	void shutdown();

	int getQueueSize();

	void execute(Generator generator, boolean force);
}
