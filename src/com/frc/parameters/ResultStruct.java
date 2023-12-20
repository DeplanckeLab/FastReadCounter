package com.frc.parameters;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import com.errors.ErrorMessage;
import com.tools.DNAHashing;
import com.tools.Utils;

import gnu.trove.set.hash.THashSet;

public class ResultStruct 
{
	public String barcode;

	public int foundGXTag = 0;
	public int foundUTag = 0;
	public long nbReads = 0;
	public long duplicates = 0;
	public int unmapped = 0;
	public int notUnique = 0;
	public int ambiguous = 0;
	public int mapped = 0;
	public int noFeature = 0;
	public int toolowAqual = 0;
	public int[] counts;
	public THashSet<?>[] umis; // THashSet is used instead of HashSet because more efficient (especially in terms of RAM usage)
	boolean isLong = false;
	
	//paired-end
	public int mateUnmapped = 0;
	public int properlyPaired = 0;
	public int countedPair = 0;
	public int countedUnique = 0;
	
	public ResultStruct() // Used for summary
	{
		this.counts = null;
		this.barcode = null;
		this.umis = null;
	}
	
	public ResultStruct(int nbGenes, String barcode)
	{
		this.counts = new int[nbGenes]; // genes
		this.barcode = barcode;
		this.umis = null;
	}
	
	public void addSummary(ResultStruct tmp)
	{
		this.nbReads += tmp.nbReads;
		this.unmapped += tmp.unmapped;
		this.notUnique += tmp.notUnique;
		this.ambiguous += tmp.ambiguous;
		this.mapped += tmp.mapped;
		this.noFeature += tmp.noFeature;
		this.toolowAqual += tmp.toolowAqual;
		this.duplicates += tmp.duplicates;
		this.foundGXTag += tmp.foundGXTag;
		this.foundUTag += tmp.foundUTag;
		
		//paired-end
		this.mateUnmapped += tmp.mateUnmapped;
		this.properlyPaired += tmp.properlyPaired;
		this.countedPair += tmp.countedPair;
		this.countedUnique += tmp.countedUnique;
	}
	
	public void add(ResultStruct tmp)
	{
		if(!this.barcode.equals(tmp.barcode)) new ErrorMessage("Cannot merge different barcodes");
		this.addSummary(tmp);
		if(tmp.isLong) this.isLong = true;
		
		// Count matrix
		if(this.counts.length != tmp.counts.length) new ErrorMessage("read count array length are different for sample " + this.barcode);
		for(int i = 0; i < this.counts.length; i++)
		{
			this.counts[i] += tmp.counts[i];
		}
		
		// UMI matrix
		if(this.umis != null && tmp.umis == null) new ErrorMessage("UMI array length are different for sample " + this.barcode);
		if(this.umis == null && tmp.umis != null) new ErrorMessage("UMI array length are different for sample " + this.barcode);
		if(this.umis != null && tmp.umis != null)
		{
			if(this.umis.length != tmp.umis.length) new ErrorMessage("UMI array length are different for sample " + this.barcode);
			for(int i = 0; i < this.umis.length; i++)
			{
				// TODO
				//if(this.isLong) (THashSet<Long>)(this.umis[i]).addAll((THashSet<Long>)tmp.umis[i]);
				//else (THashSet<Integer>)(this.umis[i]).addAll((THashSet<Integer>)tmp.umis[i]);
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	public void addUMI(int index, String umi)
	{
		if(this.umis == null) // Init here?
		{	
			this.umis = new THashSet<?>[this.counts.length]; // umis
			if(umi.length() <= 16) for(int i = 0; i < this.umis.length; i++) this.umis[i] = new THashSet<Integer>();
			else if(umi.length() <= 32) { isLong = true; for(int i = 0; i < this.umis.length; i++) this.umis[i] = new THashSet<Long>(); }
			else new ErrorMessage("Cannot procede UMIs greater than 32nt");
		}
		((THashSet<Integer>)this.umis[index]).add(DNAHashing.convertNuclStrToInt32(umi));
	}
	
	public static void createOutputDGE(HashMap<String, ResultStruct> results)
	{
		// Merge
		int foundUTag = 0;
		for(String barcode:results.keySet()) foundUTag += results.get(barcode).foundUTag;
		
		// Write output
		try
		{
			// Create the read count and transcript count matrices
			BufferedWriter bw_reads = new BufferedWriter(new FileWriter(Parameters.outputFolder + "counts.txt"));
			BufferedWriter bw_reads_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "counts.detailed.txt"));
			// Create the read count and transcript count matrices
			BufferedWriter bw_umis = null;
			BufferedWriter bw_umis_detailed = null;
			if(Parameters.use_bam_tags && foundUTag != 0)
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
