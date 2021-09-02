package model;

import java.util.HashMap;

public class Global 
{
	public static HashMap<String, String> mappingGeneIdGeneName = null;
	public static HashMap<String, Integer> geneIndex = null;
	public static HashMap<String, String> mappingBarcodeName = null;

	public static long unmapped = 0;
	public static long notUnique = 0;
	public static long ambiguous = 0;
	public static long mapped = 0;
	public static long noFeature = 0;
	public static long toolowAqual = 0;
	public static long notDemultiplexed = 0;
	public static long nbReads = 0;
	public static long duplicates = 0;
	public static int foundGXTag = 0;
	
	//paired-end
	public static int mateUnmapped = 0;
	public static int properlyPaired = 0;
	public static int countedPair = 0;
	public static int countedUnique = 0;
}
