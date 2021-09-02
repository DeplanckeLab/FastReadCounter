package model;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import tools.Utils;

public class ResultStruct 
{
	public String barcode;
	public long nbReads;
	public int unmapped;
	public int notUnique;
	public int ambiguous;
	public int mapped;
	public int noFeature;
	public int toolowAqual;
	public int[] counts;
	
	//paired-end
	public int mateUnmapped;
	public int properlyPaired;
	public int countedPair = 0;
	public int countedUnique = 0;
	
	public ResultStruct(int nbGenes)
	{
		this.counts = new int[nbGenes]; // genes 
	}
	
	public static void createOutputDGE(HashMap<String, ResultStruct> results)
	{
		try
		{
			// Create the read count and transcript count matrices
			BufferedWriter bw_reads = new BufferedWriter(new FileWriter(Parameters.outputFolder + "counts.txt"));
			BufferedWriter bw_reads_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "counts.detailed.txt"));
			
			// Sorted indexes
			String[] sortedGeneKeys = Utils.sortKeys(Global.geneIndex);
			String[] sortedBarcodeKeys = Utils.sortKeysByValues(Global.mappingBarcodeName);
			
			// Header
			bw_reads.write("Gene_id"); 
			bw_reads_detailed.write("Gene_id\tGene_name"); 
			for(String barcode:sortedBarcodeKeys) 
			{ 
				String mappedBarcode = Global.mappingBarcodeName.get(barcode);
				if(mappedBarcode != null) bw_reads.write("\t" + mappedBarcode); 
				if(mappedBarcode == null) mappedBarcode = "Unknown_Barcode";
				bw_reads_detailed.write("\t" + mappedBarcode); 
			}
			bw_reads.write("\n");  
			bw_reads_detailed.write("\n"); 
			
			// Actual values
			for(String gene:sortedGeneKeys)
			{
				String mappedGene = Global.mappingGeneIdGeneName.get(gene);
				if(mappedGene == null) mappedGene = "";
				bw_reads.write(gene);
				bw_reads_detailed.write(gene + "\t" + mappedGene); 
				for(String barcode:sortedBarcodeKeys) 
				{
					ResultStruct res = results.get(barcode);
					bw_reads.write("\t" + res.counts[Global.geneIndex.get(gene)]);
					bw_reads_detailed.write("\t" + res.counts[Global.geneIndex.get(gene)]);
				}
				bw_reads.write("\n");
				bw_reads_detailed.write("\n");
			}
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__no_feature\t__no_feature");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.noFeature);
			}
			bw_reads_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__ambiguous\t__ambiguous");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.ambiguous);
			}
			bw_reads_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__too_low_aQual\t__too_low_aQual");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.toolowAqual);
			}
			bw_reads_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__not_aligned\t__not_aligned");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.unmapped);
			}
			bw_reads_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__alignment_not_unique\t__alignment_not_unique");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.notUnique);
			}
			bw_reads_detailed.write("\n");
					
			bw_reads.close(); 
			bw_reads_detailed.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		System.out.println("\nCount matrix written in " + Parameters.outputFolder);
	}
}
