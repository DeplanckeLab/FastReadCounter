package com.frc.vcf;

import java.io.BufferedReader;
import java.util.HashMap;
import java.util.HashSet;

import com.frc.parameters.Global;
import com.frc.parameters.Parameters;
import com.intervals.Forest;
import com.intervals.IntervalLabelled;
import com.tools.Utils;

public class VCF
{
	private static String[] knownCols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};
	private static HashSet<Character> authorized = new HashSet<Character>(); // Check if the variants are correctly formed
	
	static {
        authorized.add('a'); authorized.add('c'); authorized.add('g'); authorized.add('t');
        authorized.add('A'); authorized.add('C'); authorized.add('G'); authorized.add('T');
        authorized.add('N'); authorized.add('.');
	}
	
	public static void readVCF() throws Exception
	{
		System.out.println("\nReading VCF file provided: " + Parameters.inputVCFFile.getAbsolutePath());
		Global.geneIndex = new HashMap<String, Integer>(); // Fill this as we go
		Global.mappingGeneIdGeneName = new HashMap<String, String>(); // For filling final matrices
		BufferedReader br = Utils.readVCF(Parameters.inputVCFFile);
        
        // Reading Header
		int l = 0;
		Forest.init();
		HashMap<String, Integer> indexes = new HashMap<>();
        String line = br.readLine(); l++;
        while(!line.startsWith("#CHROM")) {line = br.readLine(); l++;}
    	String[] header = line.substring(1).split("\t"); // remove the # and split
    	for (int i = 0; i < header.length; i++) indexes.put(header[i], i);
    	System.out.println((indexes.size() - knownCols.length) + " samples found in VCF [Will not be used]");

    	long start = System.currentTimeMillis();
        line = br.readLine(); l++;
        int count = 0, nbFeatures = 0;
        while(line != null)
        {
        	count++;
            String[] tokens = line.split("\t");
        	SNP s = new SNP();
        	s.chr = tokens[indexes.get("CHROM")];
        	s.loc = Long.parseLong(tokens[indexes.get("POS")]);
            s.refAllele = tokens[indexes.get("REF")]; // May be different than in REF
            s.altAllele = tokens[indexes.get("ALT")]; // May be different than in REF
            s.rsId = tokens[indexes.get("ID")];

            if(checkIntegrity(s.refAllele, l) && checkIntegrity(s.altAllele, l))
            {
   				if(s.refAllele.length() == 1 && s.altAllele.length() == 1) //TODO consider indels
   				{
   					// Parse Line
   					String feature_id = s.chr + ":" + s.loc;
   					boolean strand = true;
   					s.isSNP = true;
   					
   					// TODO SPECIFIC CODE TO REMOVE
   					String geno1 = tokens[indexes.get("FVB_NJ")].split(":")[0];
   					int which1 = 0;
   					switch(geno1) {
   						case ".": which1 = 0; break;
   						case "1/1": which1 = 2; break;
   						case "1/0": which1 = 1; break;
   						case "0/1": which1 = 1; break;
   						case "0/0": which1 = 0; break;
   						default:
   							System.err.println("ERROR GENOTYPE: " + geno1);
   					}
   					
   					String geno2 = tokens[indexes.get("DBA_1J")].split(":")[0];
   					int which2 = 0;
   					switch(geno2) {
   						case ".": which2 = 0; break;
   						case "1/1": which2 = 2; break;
   						case "1/0": which2 = 1; break;
   						case "0/1": which2 = 1; break;
   						case "0/0": which2 = 0; break;
   						default:
   							System.err.println("ERROR GENOTYPE: " + geno2);
   					}
   					
   					if((which1 == 0 && which2 == 2) || (which1 == 2 && which2 == 0))
   					{
   						// TODO SPECIFIC CODE TO REMOVE
	   					Forest.addToTree(s.chr, new IntervalLabelled((int)s.loc - 1, (int)s.loc, feature_id, strand));
	   					Global.geneIndex.put(feature_id, nbFeatures);
	   					Global.mappingGeneIdGeneName.put(feature_id, s.rsId);
	   					nbFeatures++;
   					}
   				}
            	
   				if(nbFeatures % 1000000 == 0) System.out.println(nbFeatures + " SNPs read from VCF [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
	           	/*String geno = tokens[indexes.get(Parameters.sample1)];
	           	int pos = geno.indexOf(":");
	           	if(pos != - 1) geno = geno.substring(0, geno.indexOf(":"));
	           	int genotype1 = Utils.getGeno(geno);
	           	if(genotype1 == 0 || genotype1 == 2) // If Sample1 is homozygous
	           	{
	           		geno = tokens[indexes.get(Parameters.sample2)];
	           		pos = geno.indexOf(":");
	           		if(pos != - 1) geno = geno.substring(0, geno.indexOf(":"));
	           		int genotype2 = Utils.getGeno(geno);
	           		if((genotype1 == 0 && genotype2 == 2) || (genotype1 == 2 && genotype2 == 0)) 
	           		{
	           			het++;
	           			HashSet<IntervalLabelled> genes = GTF.findOverlappingGenes(s.chr, (int)s.loc, (int)(s.loc + 1));
	           			for(IntervalLabelled i:genes) 
	           			{
	           				if(s.refAllele.length() == 1 && s.altAllele.length() == 1) //TODO consider indels
	           				{
	           					Parameters.notIndels++;
	           					if(genotype1 == 0) s.isRef = 1;
	           					else s.isRef = 2;
	           					s.isSNP = true;
	           					i.variants.add(s);
	           				}
	           				inGenes++;
	           			}
	           		}
	           	}*/
            }
            line = br.readLine(); l++;
        }
        
        br.close();
        
        System.out.println(nbFeatures + " SNPs read from VCF [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
        System.out.println("The VCF File contained " + count + " variants");
        System.out.println(nbFeatures + " variants are SNPs (not INDELs) and thus will be used in the rest of the pipeline");
	}
		
	private static boolean checkIntegrity(String allele, int l)
	{
		String[] tokens = allele.split(",");
		HashSet<String> set = new HashSet<String>();
		for(String all:tokens)
		{
			if(set.contains(all))
			{
				System.out.println("[IGNORED] Malformed VCF at l. "+l+": Duplicate allele added to VariantContext: " + all); 
				return false; 
			}
			set.add(all);
			if(all.equals("")) 
			{ 
				System.out.println("[IGNORED] Malformed VCF at l. "+l+": empty alleles are not permitted in VCF records"); 
				return false; 
			}
	        for(int i = 0; i < all.length(); i++) 
	        {
	        	if(!authorized.contains(all.charAt(i))) 
	        	{ 
	        		System.out.println("[IGNORED] Malformed VCF at l. "+l+": unparsable vcf record with allele "+all.charAt(i)); 
	        		return false;
	        	}
	        }
		}
		return true;
	}
}
