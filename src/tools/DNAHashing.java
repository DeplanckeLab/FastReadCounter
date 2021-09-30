package tools;

import java.util.HashMap;

import model.ErrorMessage;

/**
 * This class was implemented to replace the String.hash() function because the latter does not provide sufficient diversity for DNA sequences.
 * Testing the String.hash() function sometimes produces the same hash for two different DNA strings (~4 over 10M sequences from BAM file)
 * This bijective function should produce a unique value for {A,C,G,T,N,-} sequences, and runs approximately as fast as the .hash() function 
 * 
 * @author Vincent Gardeux vincent.gardeux@epfl.ch
 *
 */
public class DNAHashing 
{
	private static long[] powers = new long[80]; 
	private static HashMap<Character, Integer> map = new HashMap<Character, Integer>();
	
	static
	{
		for(int i = 0; i < powers.length; i++) powers[i] = (long)Math.pow(6, i);
		map.put('A', 0);
		map.put('C', 1);
		map.put('G', 2);
		map.put('T', 3);
		map.put('N', 4);
		map.put('-', 5);
	}
	
	public static long hashDNASequence(String seq)
	{
		int l = seq.length();
		if(l > 80) new ErrorMessage("This sequence is too big to be hashed (greater value than maximum long)");
		long hash = 0;
		for(int i = 0; i < l; i++)
		{
			Integer m = map.get(seq.charAt(i));
			hash += m*powers[i];
		}
		return hash;
	}
}
