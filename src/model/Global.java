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
}
