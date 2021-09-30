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
				GTFLine infos = parseGTFLine(line);
				if(infos != null)
				{
		    		if(infos.gene_name == null) infos.gene_name = infos.gene_id;
					// Which type is it?
					if(infos.type.equals("exon")) 
					{
						nbExons++;
						if(!Parameters.use_bam_tags) Forest.addToTree(infos.chr, new IntervalLabelled((int)infos.start, (int)infos.end, infos.gene_id, infos.strand));
						if(uniqueGeneId.add(infos.gene_id)) uniqueGeneName.add(infos.gene_name);
					}
					else if(infos.type.equals("gene"))
					{
						Global.geneIndex.put(infos.gene_id, nbGenes);
						Global.mappingGeneIdGeneName.put(infos.gene_id, infos.gene_name);
						nbGenes++;
					}
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
	
	/**
	 * Fastest parsing of a GTF line
	 * @param line
	 * @return
	 */
	private static GTFLine parseGTFLine(String line)
	{
		GTFLine res = new GTFLine();
		int pos = 0, end = 0;
		char c;
		
		// Chromosome
		while((c = line.charAt(end)) != '\t') end++; // or line.indexOf() should be equivalent
		res.chr = line.substring(pos, end);
		end++; pos = end;
		
		// Skip source
		while((c = line.charAt(end)) != '\t') end++;
		end++; pos = end;
		
		// Type
		while((c = line.charAt(end)) != '\t') end++;
		res.type = line.substring(pos, end);
		if(!res.type.equals("exon") && !res.type.equals("gene")) return null;
		end++; pos = end;
		
		// Start
		while((c = line.charAt(end)) != '\t') end++; // or line.indexOf() should be equivalent
		res.start = Long.parseLong(line.substring(pos, end));
		end++; pos = end;

		// End
		while((c = line.charAt(end)) != '\t') end++; // or line.indexOf() should be equivalent
		res.end = Long.parseLong(line.substring(pos, end));
		end++; pos = end;
		
		// Ignore
		while((c = line.charAt(end)) != '\t') end++;
		end++; pos = end;
		
		// Strand
		while((c = line.charAt(end)) != '\t') end++; // or line.indexOf() should be equivalent
		res.strand = line.substring(pos, end).equals("+") ;
		end++; pos = end;
		
		// Ignore
		while((c = line.charAt(end)) != '\t') end++;
		end++; pos = end;
		
		// Info
		int found = 0;
		// I don't know which one gene_id / gene_name comes first, so using indexOf() will run twice the sequence from the start 
		for(int i = pos; i < line.length(); i++) // So ugly... But should be optimized.... but sorry for that...
		{
			if(found == 2) break;
			if(line.charAt(i) == 'g')
			{
				i++;
				if(line.charAt(i) == 'e')
				{
					i++;
					if(line.charAt(i) == 'n')
					{
						i++;
						if(line.charAt(i) == 'e')
						{
							i++;
							if(line.charAt(i) == '_')
							{
								i++;
								if(line.charAt(i) == 'i')
								{
									i++;
									if(line.charAt(i) == 'd')
									{
										i++;				
										c = line.charAt(i);
										while(c == ' ' || c == '\"') { i++; c = line.charAt(i);}
										pos = i;
										while(c != ' ' && c != '\"') { i++; c = line.charAt(i); }
										res.gene_id = line.substring(pos, i);
										found++;
									}
								}
								else if(line.charAt(i) == 'n')
								{
									i++;
									if(line.charAt(i) == 'a')
									{
										i++;
										if(line.charAt(i) == 'm')
										{
											i++;
											if(line.charAt(i) == 'e')
											{
												i++;				
												c = line.charAt(i);
												while(c == ' ' || c == '\"') { i++; c = line.charAt(i);}
												pos = i;
												while(c != ' ' && c != '\"') { i++; c = line.charAt(i); }
												res.gene_name = line.substring(pos, i);
												found++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		return res;
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

class GTFLine
{
	long start;
	long end;
	String chr;
	String type;
	boolean strand;
	String gene_id;
	String gene_name;
}
