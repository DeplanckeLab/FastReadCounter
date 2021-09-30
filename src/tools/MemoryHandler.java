package tools;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryType;

public class MemoryHandler
{
	public static long maxEden;
	public static long maxSurvivor;
	public static long maxOld;
	public static long maxHeap; // Max memory of heap in bytes
	public static long maxHeap_legacy; // Max memory of heap in bytes (does not account for Eden/Survivor)
	
	static
	{
        maxHeap = 0l;
        maxHeap_legacy = Runtime.getRuntime().maxMemory();
        for (MemoryPoolMXBean mp : ManagementFactory.getMemoryPoolMXBeans()) 
        {
        	if(mp.getType() == MemoryType.HEAP)
        	{
        		if(mp.getName().contains("Survivor")) 
        		{
        			maxSurvivor = mp.getUsage().getMax();
        			maxHeap += mp.getUsage().getMax(); // 2* survivor
        		}
        		else if(mp.getName().contains("Eden")) maxEden = mp.getUsage().getMax();
        		maxHeap += mp.getUsage().getMax();
        	}
        }
	}
	
	// New size
	public static long getNewSize()
	{
        return maxEden + maxSurvivor;
	}
	
	// Current size of heap in bytes
	public static long getCurrentHeapSize()
	{
        return Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
	}
	
	// Current size of heap in percent
	public static float getCurrentHeapSizePercent()
	{
		float pc = ((float)getCurrentHeapSize() / maxHeap_legacy) * 100;
		if(pc > 100) pc = 100; // It can be slightly over?
        return pc;
	}
	
	public static boolean isHeapFull()
	{
		// A bit arbitrary?
		if(getCurrentHeapSizePercent() >= 90) return true;
		return false;
	}
	
	public static String toReadableSize(long bytes)
	{
        if (bytes < 1024) return bytes + " b";
        int z = (63 - Long.numberOfLeadingZeros(bytes)) / 10;
        return String.format("%.1f %sb", (double)bytes / (1L << (z*10)), " KMGTPE".charAt(z));
	}	
}
