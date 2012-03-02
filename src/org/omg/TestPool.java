package org.omg;


import java.util.concurrent.*;
import java.util.*;
 
class TestPool
{
    int poolSize = 3;
 
    int maxPoolSize = 3;
 
    long keepAliveTime = 10;
 
    ThreadPoolExecutor threadPool = null;
 
    final ArrayBlockingQueue<Runnable> queue = new ArrayBlockingQueue<Runnable>(
            5);
 
    public TestPool()
    {
        threadPool = new ThreadPoolExecutor(poolSize, maxPoolSize,
                keepAliveTime, TimeUnit.SECONDS, queue);
 
    }
 
    class MyTask implements Runnable {

    	public MyTask(int count) {
			super();
			this.count = count;
		}
    	
		private int count;
		@Override
		public void run() {
            for (int i = 0; i < count; i++)
            {
                try
                {
                    System.out.println("Task "+count);
                    Thread.sleep(1000);
                } catch (InterruptedException ie)
                {
                }
            }
		}
    	
    }
    
    public void runTask(Runnable task)
    {
        // System.out.println("Task count.."+threadPool.getTaskCount() );
        // System.out.println("Queue Size before assigning the
        // task.."+queue.size() );
        threadPool.execute(task);
        // System.out.println("Queue Size after assigning the
        // task.."+queue.size() );
        // System.out.println("Pool Size after assigning the
        // task.."+threadPool.getActiveCount() );
        // System.out.println("Task count.."+threadPool.getTaskCount() );
        System.out.println("Task count.." + queue.size());
 
    }
 
    public void shutDown()
    {
        threadPool.shutdown();
    }
 
    public static void main(String args[])
    {
        TestPool mtpe = new TestPool();
        // start first one
        mtpe.runTask(new Runnable()
        {
            public void run()
            {
                for (int i = 0; i < 10; i++)
                {
                    try
                    {
                        System.out.println("First Task");
                        Thread.sleep(1000);
                    } catch (InterruptedException ie)
                    {
                    }
                }
            }
        });
        
        mtpe.runTask(mtpe.new MyTask(3));
        mtpe.runTask(mtpe.new MyTask(4));
        mtpe.runTask(mtpe.new MyTask(5));
        mtpe.runTask(mtpe.new MyTask(6));
        mtpe.runTask(mtpe.new MyTask(7));
        
        
//        // start second one
//        /*
//         * try{ Thread.sleep(500); }catch(InterruptedException
//         * ie){}
//         */
//        mtpe.runTask(new Runnable()
//        {
//            public void run()
//            {
//                for (int i = 0; i < 10; i++)
//                {
//                    try
//                    {
//                        System.out.println("Second Task");
//                        Thread.sleep(1000);
//                    } catch (InterruptedException ie)
//                    {
//                    }
//                }
//            }
//        });
//        // start third one
//        /*
//         * try{ Thread.sleep(500); }catch(InterruptedException
//         * ie){}
//         */
//        mtpe.runTask(new Runnable()
//        {
//            public void run()
//            {
//                for (int i = 0; i < 10; i++)
//                {
//                    try
//                    {
//                        System.out.println("Third Task");
//                        Thread.sleep(1000);
//                    } catch (InterruptedException ie)
//                    {
//                    }
//                }
//            }
//        });
//        // start fourth one
//        /*
//         * try{ Thread.sleep(500); }catch(InterruptedException
//         * ie){}
//         */
//        mtpe.runTask(new Runnable()
//        {
//            public void run()
//            {
//                for (int i = 0; i < 10; i++)
//                {
//                    try
//                    {
//                        System.out.println("Fourth Task");
//                        Thread.sleep(1000);
//                    } catch (InterruptedException ie)
//                    {
//                    }
//                }
//            }
//        });
//        // start fifth one
//        /*
//         * try{ Thread.sleep(500); }catch(InterruptedException
//         * ie){}
//         */
//        mtpe.runTask(new Runnable()
//        {
//            public void run()
//            {
//                for (int i = 0; i < 10; i++)
//                {
//                    try
//                    {
//                        System.out.println("Fifth Task");
//                        Thread.sleep(1000);
//                    } catch (InterruptedException ie)
//                    {
//                    }
//                }
//            }
//        });
//        // start Sixth one
//        /*
//         * try{ Thread.sleep(500); }catch(InterruptedException
//         * ie){}
//         */
//        mtpe.runTask(new Runnable()
//        {
//            public void run()
//            {
//                for (int i = 0; i < 10; i++)
//                {
//                    try
//                    {
//                        System.out.println("Sixth Task");
//                        Thread.sleep(1000);
//                    } catch (InterruptedException ie)
//                    {
//                    }
//                }
//            }
//        });
//        mtpe.shutDown();
//        System.out.println("Shutdown complete");
    }
 
}
