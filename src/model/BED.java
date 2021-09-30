package model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

public class BED 
{
	public static void readBED() throws Exception
	{
		System.out.println("\nReading BED file provided: " + Parameters.inputBEDFile.getAbsolutePath());
		Global.geneIndex = new HashMap<String, Integer>(); // Fill this as we go
		Global.mappingGeneIdGeneName = new HashMap<String, String>(); // For filling final matrices
		BufferedReader br = openBED(Parameters.inputBEDFile);
		String line = br.readLine();
		int l = 1;
		while(line.startsWith("browser") || line.startsWith("track")) { line = br.readLine(); l++;} // Skip these
		if(!Parameters.use_bam_tags) Forest.init();
		int nbFeatures = 0;
		while(line != null)
		{
			String[] tokens = line.split("\t");
			if(tokens.length < 3) new ErrorMessage("BED entry should contain at least 3 values at l." + l);

			// Parse Line
			String chr = tokens[0];
			long start = Long.parseLong(tokens[1]);
			long end = Long.parseLong(tokens[2]);
			String feature_id = chr + ":" + start + ":" + end;
			String feature_name = (tokens.length >= 4)?tokens[3]:"";
			boolean strand = (tokens.length >= 6)?tokens[5].equals("+"):true;
			
			Integer index = Global.geneIndex.get(feature_id);
			if(index == null)
			{
				// Here we use start+1 because BED files are 0-indexed
				if(!Parameters.use_bam_tags) Forest.addToTree(chr, new IntervalLabelled((int)(start + 1), (int)end, feature_id, strand));
				Global.geneIndex.put(feature_id, nbFeatures);
				Global.mappingGeneIdGeneName.put(feature_id, feature_name);
				nbFeatures++;
			}
			else System.err.println("l. " + l + ": " + feature_id + " already appeared in your BED file. All extra occurrences are ignored.");
			
			line = br.readLine(); l++;
		}
		br.close();
		
		System.out.println("In total " + nbFeatures + " annotations/features are found in the BED file.");
		
		if(nbFeatures == 0) new ErrorMessage("We couldn't parse the BED file. Please report this problem if the BED is in standard format.");
	}
	
	private static BufferedReader openBED(File bed) throws Exception
	{
		if(bed.getAbsolutePath().endsWith(".bed"))
		{
			return new BufferedReader(new FileReader(bed));
		}
		else if(bed.getAbsolutePath().endsWith(".bed.gz"))
		{
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(bed))));
		}
		else
		{
			System.err.println("The extension of the BED file is not recognized : " + bed.getAbsolutePath());
			System.err.println("It should be '.bed', or '.bed.gz'");
			System.exit(-1);
		}
		return null;
	}
}
