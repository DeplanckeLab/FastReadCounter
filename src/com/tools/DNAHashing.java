package com.tools;

import java.util.HashMap;

/**
 * This class was implemented to replace the String.hash() function because the latter does not provide sufficient diversity for DNA sequences.
 * Testing the String.hash() function sometimes produces the same hash for two different DNA strings (~4 over 10M sequences from BAM file)
 * This bijective function should produce a unique value for {A,C,G,T} sequences ({-,N} are taken as Gs), and runs approximately as fast as the .hash() function 
 * 
 * @author Vincent Gardeux vincent.gardeux@epfl.ch
 *
 */
public class DNAHashing 
{
	private static HashMap<Character, Integer> map = new HashMap<Character, Integer>();
	private static char[] nuclChar= new char[] {'A', 'C', 'G', 'T'};
	
	static // A map is faster than a switch:case
	{
		map.put('A', 0);
		map.put('a', 0);
		map.put('C', 1);
		map.put('c', 1);
		map.put('G', 2);
		map.put('g', 2);
		map.put('T', 3);
		map.put('t', 3);
		map.put('N', 2); // return 'G' as default
		map.put('-', 2); // return 'G' as default
	}
	
	/*
	 * Adapted to Java from https://github.com/alexdobin/STAR/
	 */
	public static int convertNuclStrToInt32(String S) // We can use int for up to length 16nt
	{
		int intOut = 0;
	    for(int ii = 0; ii < S.length(); ii++)
	    {
	        int nt = map.get(S.charAt(ii)); // Ns and - are transformed in G for avoiding to reject the read
	        intOut = intOut << 2;
	        intOut += nt;
	    }
	    return intOut;
	}
	
	public static String convertNuclInt32toString(int nuclNum, int length)
	{
		char[] nuclOut = new char[length];

		for(int ii=1; ii<=length; ii++)
	    {
	        nuclOut[length-ii] = nuclChar[nuclNum & 3];
	        nuclNum = nuclNum >> 2;
	    }

	    return new String(nuclOut);
	}
	
	public static long convertNuclStrToInt64(String S) // We can use long for up to length 32nt
	{
	    long intOut = 0;
	    for (int ii = 0; ii < S.length(); ii++)
	    {
	        long nt = map.get(S.charAt(ii));
	        intOut = intOut << 2;
	        intOut += nt;
	    }
	    return intOut;
	}
	
	public static String convertNuclInt64toString(long nuclNum, int length)
	{
		char[] nuclOut = new char[length];

	    for(int ii=1; ii<=length; ii++)
	    {
	        nuclOut[length-ii] = nuclChar[(int)(nuclNum & 3)];
	        nuclNum = nuclNum >> 2;
	    }

	    return new String(nuclOut);
	}
}
