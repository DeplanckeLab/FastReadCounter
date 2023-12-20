package com.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import com.errors.ErrorMessage;
import com.frc.parameters.Global;
import com.frc.parameters.Parameters;
import com.intervals.Forest;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;

public class Utils 
{
	public static DecimalFormat pcFormatter = new DecimalFormat("##.##");
	public static Long maxHeap = null;
	/**
	 * Read VCF different types of VCF file
	 * @param VCF input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static BufferedReader readVCF(File vcf) throws Exception
	{
		if(vcf.getAbsolutePath().endsWith(".vcf"))
		{
			return new BufferedReader(new FileReader(vcf));
		}
		else if(vcf.getAbsolutePath().endsWith(".vcf.gz"))
		{
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcf))));
		}
		else
		{
			System.err.println("The extension of the VCF file is not recognized : " + vcf.getAbsolutePath());
			System.err.println("It should be '.vcf', or '.vcf.gz'");
			System.exit(-1);
		}
		return null;
	}

	/**
	 * 
	 * @param barcodeFile Required to map the barcode with the sample name
	 * @return Unique set of all barcodes
	 */
	public static HashSet<String> readBarcodeFile(String barcodeFile)
	{
		HashSet<String> res = new HashSet<String>();
		Global.mappingBarcodeName = new HashMap<String, String>();
		Parameters.barcodeLength = -1;
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(barcodeFile));
			String line = br.readLine();
			int l = 0;
			while(line != null)
			{
				l++;
				String[] tokens = line.split("\t");
				if(tokens.length > 2) new ErrorMessage("Error on l." + l + ": More than two values per line are not authorized");
				if(tokens.length >= 1) 
				{
					String barcode = tokens[0].trim();
					if(Parameters.barcodeLength == -1) Parameters.barcodeLength = barcode.length();
					if(Parameters.barcodeLength != barcode.length()) new ErrorMessage("Error on l." + l + ": All barcodes should have the same length. Check your barcode file.");
					String mappedSample = barcode;
					if(tokens.length > 1) mappedSample = tokens[1].trim();
					Global.mappingBarcodeName.put(barcode, mappedSample);
					res.add(barcode);
				}
				line = br.readLine();
			}
			br.close();
		}
		catch(Exception e)
		{
			new ErrorMessage(e.getMessage());
		}
		return res;
	}
	
	public static String saveBarcode(String barcode, ArrayList<String> barcodes, int allowedDiff)
	{
		int[] diff = new int[barcodes.size()];
		for (int i = 0; i < barcode.length(); i++) 
		{
			for (int j = 0; j < barcodes.size(); j++) 
				if(barcode.charAt(i) != barcodes.get(j).charAt(i)) diff[j]++;
		}
		String sav = "";
		for (int i = 0; i < diff.length; i++) 
		{
			if(diff[i] == 0) 
			{
				sav = barcodes.get(i);
				break;
			}
			if(diff[i] <= allowedDiff)
			{
				if(!sav.equals(""))
				{
					//System.out.println(barcode + " cannot be saved because several barcodes can match it with " + allowedDiff + " errors.");
					sav = "";
					break;
				}
				sav = barcodes.get(i);
			}
		}
		if(sav.equals("")) return null;
		return sav;
	}
	
	public static HashSet<String> getOverlappingFeatures(String chr, int start, int end, Cigar c, boolean readNegativeStrandFlag)
	{
		if(Parameters.is_paired) new ErrorMessage("You should use this function with the 'firstOfPair' argument with paired-end data");
		return getOverlappingFeatures(chr, start, end, c, readNegativeStrandFlag, false); // We put false, it has no impact anyways
	}
	
	public static HashSet<String> getOverlappingFeatures(String chr, int start, int end, Cigar c, boolean readNegativeStrandFlag, boolean firstOfPair)
	{
		HashSet<String> res = new HashSet<>();
		List<CigarElement> l = c.getCigarElements();
		int s = start;
		for(CigarElement cigar:l)
		{
			switch(cigar.getOperator())
			{
				case M:
					res.addAll(Forest.findOverlappingFeatures(chr, s, s + cigar.getLength() - 1, readNegativeStrandFlag, firstOfPair)); // -1 Because the last letter is at the index before
					s += cigar.getLength();
					break;
				case N:
					s += cigar.getLength();
					break;
				case D:
					s += cigar.getLength();
					break;
				case EQ:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case H:
					// Hard clipping. Do nothing ? (alignment start is after any H & alignment end before any H)
					break;
				case I:
					// Do nothing
					break;
				case P:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case S:
					// Soft clipping. Do nothing (alignment start is after any S & alignment end before any S)
					break;
				case X:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
			}
		}
		
		if(s != end + 1) new ErrorMessage("Error while reading CIGAR");
		return res;
	}
	
	public static String[] sortKeys(Map<String, Integer> map)
	{
		String[] keys = map.keySet().toArray(new String[map.keySet().size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static String[] sort(Set<String> map)
	{
		String[] keys = map.toArray(new String[map.size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static String[] sortKeysByValues(Map<String, String> map)
	{
		List<String> values = new ArrayList<>(map.values());
		Collections.sort(values);
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	/* 
	 * Nb "cores"
	 */
	public static int getNbLogicalThreads()
	{
		return Runtime.getRuntime().availableProcessors();
	}
	
	/** Reads len bytes in a loop.
     * 
     * This function was taken from the great implementation of jvarkit from Pierre Lindenbaum: https://github.com/lindenb
     * 
     * @param in The InputStream to read from
     * @param buf The buffer to fill
     * @param off offset from the buffer
     * @param len the length of bytes to read
     * @throws IOException if it could not read requested number of bytes 
     * for any reason (including EOF)
     */
    public static void readFully(final InputStream in, final byte buf[], int off, int len ) throws IOException 
    {
    	int toRead = len;
    	while ( toRead > 0 )
    	{
    		int ret = in.read( buf, off, toRead );
    		if( ret < 0 ) throw new IOException("Premature EOF from inputStream");
    		toRead -= ret;
    		off += ret;
    	}
    }
     
    /** Reads len bytes in a loop.
     * 
     * This function was taken from the great implementation of jvarkit from Pierre Lindenbaum: https://github.com/lindenb
     * @param in The InputStream to read from
     * @param buf The buffer to fill
     * @throws IOException if it could not read requested number of bytes 
     * for any reason (including EOF)
     */
    public static void readFully(final InputStream in, final byte buf[]) throws IOException 
    {
    	readFully(in, buf, 0, buf.length);
    }
    
    /**
     * This function was taken from the great implementation of jvarkit from Pierre Lindenbaum: https://github.com/lindenb
     */
    public static SAMSequenceDictionary extractRequired(final SAMFileHeader h)
    {
    	if(h == null) new ErrorMessage("Cannot extract dictionary because SAM header was not provided.");
    	final SAMSequenceDictionary dict = h.getSequenceDictionary();
    	if(dict==null || dict.isEmpty()) new ErrorMessage("<bam>");
    	return dict;
    }
    
    public static long getFileSize(String path)
    {
		try
		{
			return Files.size(Paths.get(path));
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		return -1;
    }
	
	public static String toReadableTime(long ms)
	{
		if(ms < 1000) return ""+ms+" ms";
		long s = ms / 1000;
		ms = ms % 1000;
		if(s < 60) return s+" s "+ms+" ms";
		long mn = s / 60;
		s = s % 60;
		if(mn < 60) return mn +" mn "+s+" s "+ms+" ms";
		long h = mn / 60;
		mn = mn % 60;
		if(h < 24) return h +" h "+ mn +" mn "+s+" s "+ms+" ms";
		long d = h / 24;
		h = h % 24;
		return d+ " d " + h +" h "+ mn +" mn "+s+" s "+ms+" ms";
	}
}
