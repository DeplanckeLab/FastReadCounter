package model;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import gnu.trove.set.hash.THashSet;
import tools.DNAHashing;
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
	public THashSet<?>[] umis; // THashSet is used instead of HashSet because more efficient (especially in terms of RAM usage)
	
	//paired-end
	public int mateUnmapped;
	public int properlyPaired;
	public int countedPair = 0;
	public int countedUnique = 0;
	
	public ResultStruct(int nbGenes)
	{
		this.counts = new int[nbGenes]; // genes 
		if(Parameters.umi_dedup != UMIDedup.NONE) // I only generate UMI results with this option
		{
			this.umis = new THashSet<?>[nbGenes]; // umis
			for(int i = 0; i < this.umis.length; i++) this.umis[i] = new THashSet<Long>();
		}
	}
	
	@SuppressWarnings("unchecked")
	public void addUMI(int index, String umi)
	{
		((THashSet<Long>)this.umis[index]).add(DNAHashing.hashDNASequence(umi));
	}
	
	public static void createOutputDGE(HashMap<String, ResultStruct> results)
	{
		try
		{
			// Create the read count and transcript count matrices
			BufferedWriter bw_reads = new BufferedWriter(new FileWriter(Parameters.outputFolder + "counts.txt"));
			BufferedWriter bw_reads_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "counts.detailed.txt"));
			// Create the read count and transcript count matrices
			BufferedWriter bw_umis = null;
			BufferedWriter bw_umis_detailed = null;
			if(Parameters.use_bam_tags && Global.foundUTag != 0)
			{	
				bw_umis = new BufferedWriter(new FileWriter(Parameters.outputFolder + "umis.txt"));
				bw_umis_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "umis.detailed.txt"));
			}			
			
			// Sorted indexes
			String[] sortedGeneKeys = Utils.sortKeys(Global.geneIndex);
			String[] sortedBarcodeKeys = Utils.sortKeysByValues(Global.mappingBarcodeName);
			
			// Header
			bw_reads.write("Gene_id");
			if(bw_umis != null) bw_umis.write("Gene_id");
			bw_reads_detailed.write("Gene_id\tGene_name");
			if(bw_umis_detailed != null) bw_umis_detailed.write("Gene_id\tGene_name");
			for(String barcode:sortedBarcodeKeys) 
			{ 
				String mappedBarcode = Global.mappingBarcodeName.get(barcode);
				if(mappedBarcode == null) mappedBarcode = "Unknown_Barcode";
				bw_reads.write("\t" + mappedBarcode);
				if(bw_umis != null) bw_umis.write("\t" + mappedBarcode);
				bw_reads_detailed.write("\t" + mappedBarcode);
				if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + mappedBarcode);
			}
			bw_reads.write("\n"); 
			if(bw_umis != null) bw_umis.write("\n");
			bw_reads_detailed.write("\n");
			if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
			
			// Actual values
			for(String gene:sortedGeneKeys)
			{
				String mappedGene = Global.mappingGeneIdGeneName.get(gene);
				if(mappedGene == null) mappedGene = "";
				bw_reads.write(gene);
				if(bw_umis != null) bw_umis.write(gene);
				bw_reads_detailed.write(gene + "\t" + mappedGene); 
				if(bw_umis_detailed != null) bw_umis_detailed.write(gene + "\t" + mappedGene);
				for(String barcode:sortedBarcodeKeys) 
				{
					ResultStruct res = results.get(barcode);
					bw_reads.write("\t" + res.counts[Global.geneIndex.get(gene)]);
					if(bw_umis != null) bw_umis.write("\t" + res.umis[Global.geneIndex.get(gene)].size());
					bw_reads_detailed.write("\t" + res.counts[Global.geneIndex.get(gene)]);
					if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + res.umis[Global.geneIndex.get(gene)].size());
				}
				bw_reads.write("\n"); 
				if(bw_umis != null) bw_umis.write("\n");
				bw_reads_detailed.write("\n");
				if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
			}
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__no_feature\t__no_feature");
			if(bw_umis_detailed != null) bw_umis_detailed.write("__no_feature\t__no_feature");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.noFeature);
				if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + res.noFeature);
			}
			bw_reads_detailed.write("\n");
			if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__ambiguous\t__ambiguous");
			if(bw_umis_detailed != null) bw_umis_detailed.write("__ambiguous\t__ambiguous");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.ambiguous);
				if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + res.ambiguous);
			}
			bw_reads_detailed.write("\n");
			if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__too_low_aQual\t__too_low_aQual");
			if(bw_umis_detailed != null) bw_umis_detailed.write("__too_low_aQual\t__too_low_aQual");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.toolowAqual);
				if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + res.toolowAqual);
			}
			bw_reads_detailed.write("\n");
			if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__not_aligned\t__not_aligned");
			if(bw_umis_detailed != null) bw_umis_detailed.write("__not_aligned\t__not_aligned");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.unmapped);
				if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + res.unmapped);
			}
			bw_reads_detailed.write("\n");
			if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
			
			// Write complementing values (same as HTSeq-count)
			bw_reads_detailed.write("__alignment_not_unique\t__alignment_not_unique");
			if(bw_umis_detailed != null) bw_umis_detailed.write("__alignment_not_unique\t__alignment_not_unique");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.notUnique);
				if(bw_umis_detailed != null) bw_umis_detailed.write("\t" + res.notUnique);
			}
			bw_reads_detailed.write("\n");
			if(bw_umis_detailed != null) bw_umis_detailed.write("\n");
					
			bw_reads.close(); 
			if(bw_umis != null) bw_umis.close();
			bw_reads_detailed.close();
			if(bw_umis_detailed != null) bw_umis_detailed.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		System.out.println("\nRead/UMI count matrices written in " + Parameters.outputFolder);
	}
}
