package model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

public class GTF 
{
	public static void readGTF() throws Exception
	{
		System.out.println("\nReading GTF file provided: " + Parameters.inputGTFFile.getAbsolutePath());
		Global.geneIndex = new HashMap<String, Integer>(); // Fill this as we go
		Global.mappingGeneIdGeneName = new HashMap<String, String>(); // For filling final matrices
		BufferedReader br = openGTF(Parameters.inputGTFFile);
		String line = br.readLine();
		if(!Parameters.use_bam_tags) Forest.init();
		int nbExons = 0;
		int nbGenes = 0;
		HashSet<String> uniqueGeneId = new HashSet<String>();
		ArrayList<String> uniqueGeneName = new ArrayList<String>(); // It's actually not unique and should not be since two geneID can have the same gene Name => ArrayList
		while(line != null)
		{
	    	if(!line.startsWith("#"))
    		{
				String[] tokens = line.split("\t");
				// Parse Line
				String[] params = tokens[8].split(";");
				long start = Long.parseLong(tokens[3]);
				long end = Long.parseLong(tokens[4]);
				String gene_name = null;
				String gene_id = null;
				String chr = tokens[0];
				String type = tokens[2];
				boolean strand = tokens[6].equals("+");
				for(String param:params) 
				{
					String[] values = param.trim().split("\\s+");
					values[1] = values[1].replaceAll("\"", "");
					if(values[0].equals("gene_name")) gene_name = values[1];
					else if(values[0].equals("gene_id")) gene_id = values[1];
				}
				if(gene_name == null) gene_name = gene_id;
				// Which type is it?
				if(type.equals("exon")) 
				{
					nbExons++;
					if(!Parameters.use_bam_tags) Forest.addToTree(chr, new IntervalLabelled((int)start, (int)end, gene_id, strand));
					if(uniqueGeneId.add(gene_id)) uniqueGeneName.add(gene_name);
				}
				else if(type.equals("gene"))
				{
					Global.geneIndex.put(gene_id, nbGenes);
					Global.mappingGeneIdGeneName.put(gene_id, gene_name);
					nbGenes++;
				}
			}
			line = br.readLine();
		}
		br.close();
		
		if(nbGenes == 0) {
			System.out.println("No Genes were detected in the GTF file. Probably the \"gene\" annotation is missing from the GTF file 3rd column?");
			System.out.println("Trying to \"save the day\" by collapsing exons to their annotated gene_id");
			for(String gene_id:uniqueGeneId) {
				Global.geneIndex.put(gene_id, nbGenes);
				Global.mappingGeneIdGeneName.put(gene_id, uniqueGeneName.get(nbGenes));
				nbGenes++;
			}
		}

		System.out.println(nbExons + " 'exons' are annotating " + uniqueGeneId.size() + " unique genes in the provided GTF file. In total " + nbGenes + " 'gene' annotations are found in the GTF file.");
		
		if(nbGenes == 0) {
			System.err.println("We couldn't parse the GTF file. Please report this problem if the GTF is in standard format. Or use another GTF from another source.");
			System.exit(-1);
		}
	}
	
	private static BufferedReader openGTF(File gtf) throws Exception
	{
		if(gtf.getAbsolutePath().endsWith(".gtf"))
		{
			return new BufferedReader(new FileReader(gtf));
		}
		else if(gtf.getAbsolutePath().endsWith(".gtf.gz"))
		{
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gtf))));
		}
		else
		{
			System.err.println("The extension of the GTF file is not recognized : " + gtf.getAbsolutePath());
			System.err.println("It should be '.gtf', or '.gtf.gz'");
			System.exit(-1);
		}
		return null;
	}
}
