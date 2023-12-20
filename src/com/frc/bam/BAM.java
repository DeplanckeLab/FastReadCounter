package com.frc.bam;

import java.util.HashMap;
import java.util.HashSet;

import com.errors.ErrorMessage;
import com.frc.parameters.Global;
import com.frc.parameters.Parameters;
import com.frc.parameters.ResultStruct;
import com.frc.parameters.UMIDedup;
import com.tools.Utils;

import htsjdk.samtools.SAMRecord;

public class BAM 
{
	/**
	 * Using htsjdk to process the samRecord read from the BAM file created by the alignment tool
	 * 
	 */
	public static void processRecord(SAMRecord samRecord, HashMap<String, ResultStruct> results, HashMap<String, Read> pairedBuffer)
	{
		// Resolve barcode first
		String barcode = "Unknown";
		if(Parameters.doDemultiplexing)
		{
			barcode = (String)samRecord.getAttribute("CB");
			if(barcode == null || barcode.equals("-")) barcode = "Unknown"; // New in STAR, now unknown barcodes are labelled '-' instead of the CB tag being nonexistent.
		}

		// Check if barcode is in the list of barcodes
		ResultStruct res = results.get(barcode); // This should exist
		if(res == null) new ErrorMessage("Barcode " + barcode + " is found in the BAM file, but is not in your list of barcodes");
		res.nbReads++; // Count reads by barcode
		
		// Quantify read according to its tags
		if(Parameters.is_paired && !samRecord.getReadPairedFlag()) new ErrorMessage("You used the --paired option, but it seems that some reads are not paired?");
		if(samRecord.getDuplicateReadFlag()) res.duplicates++; // Just count this, but still process it
		if(samRecord.getReadUnmappedFlag()) res.unmapped++;
		else
		{
			boolean toProcess = true;
			if(samRecord.isSecondaryOrSupplementary())
			{
				res.notUnique++;
				if(!Parameters.keep_multiple_mapped_reads) toProcess = false;
			}
			if(toProcess)
			{
				// Paired-end
				if(Parameters.is_paired)
				{
					// If MateUnmapped flag, then the read cannot have the "ProperPair" flag
					if(samRecord.getMateUnmappedFlag()) res.mateUnmapped++;
					if(samRecord.getProperPairFlag()) res.properlyPaired++;
					if(samRecord.getMateUnmappedFlag())
					{
						res.countedUnique++;
						// Process this singleton read as if it was single-end (not stored)
						processSingleEndRead(samRecord, res);
					}
					else
					{
						String name = samRecord.getReadName();
						Read read1 = pairedBuffer.remove(samRecord.getReadName());
						if(read1 != null) // I run the counting only if two pairs are identified/available
						{
							res.countedPair++;
							if(samRecord.getMappingQuality() < Parameters.minAQual || read1.mapQ < Parameters.minAQual) res.toolowAqual++; // To match htseq-count. I do this first
							else
							{
								if(Parameters.use_bam_tags)
								{
									String gene_id_1 = read1.gene;
									String gene_id_2 = (String)samRecord.getAttribute("GX");
									String umi = null;
									if(Parameters.umi_dedup != UMIDedup.NONE)
									{
										umi = (String)samRecord.getAttribute("UB");
										if(umi == null) umi = (String)samRecord.getAttribute("UR"); // Non corrected
									}
									if(gene_id_1 != null && gene_id_1.equals("-")) gene_id_1 = null;
									if(gene_id_2 != null && gene_id_2.equals("-")) gene_id_2 = null;
									if(gene_id_1 != null && gene_id_2 != null && !gene_id_1.equals(gene_id_2)) { res.foundGXTag++; res.ambiguous++; }
									else if(gene_id_1 == null && gene_id_2 == null) res.noFeature++;
									else
									{
										res.foundGXTag++;
										String gene_id = gene_id_1;
										if(gene_id == null) gene_id = gene_id_2;
										res.mapped = res.mapped + 2;
										Integer indexGene = Global.geneIndex.get(gene_id); // From GTF
										if(indexGene == null) new ErrorMessage("ERROR: This gene " + gene_id + " is not in your GTF file. Please check again.");
										res.counts[indexGene]++;
										if(umi != null) { res.foundUTag++; res.addUMI(indexGene, umi); }
									}
								}
								else // Use positions
								{
									HashSet<String> overlappingGenes_1 = Utils.getOverlappingFeatures(read1.chr, read1.startV, read1.endV, read1.cigar, read1.negativeStrandFlag, read1.firstOfPair);
									HashSet<String> overlappingGenes_2 = Utils.getOverlappingFeatures(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag(), samRecord.getFirstOfPairFlag());
									
									if(overlappingGenes_1.size() == 0 && overlappingGenes_2.size() == 0) res.noFeature++;
									else
									{
										overlappingGenes_1.addAll(overlappingGenes_2);
										if(overlappingGenes_1.size() == 1)	
										{
											res.mapped++;
											String gene = overlappingGenes_1.iterator().next();
											res.counts[Global.geneIndex.get(gene)]++;
										}
										else res.ambiguous++;
									}
								}
							}
						}
						else 
						{
							String umi = null;
							if(Parameters.umi_dedup != UMIDedup.NONE)
							{
								umi = (String)samRecord.getAttribute("UB");
								if(umi == null) umi = (String)samRecord.getAttribute("UR"); // Non corrected
							}
							pairedBuffer.put(name, new Read(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag(), samRecord.getMappingQuality(), (String)samRecord.getAttribute("GX"), umi, samRecord.getFirstOfPairFlag()));
						}
					}
				}
				else
				{
					// Single-end
					processSingleEndRead(samRecord, res);
				}
			}
		}
	}
	
	private static void processSingleEndRead(SAMRecord samRecord, ResultStruct res)
	{
		if(samRecord.getMappingQuality() < Parameters.minAQual) res.toolowAqual++; // To match htseq-count. I do this first
		else
		{	
			if(Parameters.use_bam_tags)
			{
				String gene_id = (String)samRecord.getAttribute("GX");
				//String gene_name = (String)samRecord.getAttribute("GN"); // Not used
				String umi = null;
				if(Parameters.umi_dedup != UMIDedup.NONE)
				{
					umi = (String)samRecord.getAttribute("UB");
					if(umi == null) umi = (String)samRecord.getAttribute("UR"); // Non corrected
				}
				if(gene_id == null || gene_id.equals("-")) res.noFeature++;
				else 
				{	
					res.foundGXTag++;
					res.mapped++;
					Integer indexGene = Global.geneIndex.get(gene_id); // From GTF
					if(indexGene == null) new ErrorMessage("ERROR: This gene " + gene_id + " is not in your GTF file. Please check again.");
					res.counts[indexGene]++;
					if(umi != null) { res.foundUTag++; res.addUMI(indexGene, umi); }
				}
			}
			else // Use positions
			{
				HashSet<String> overlappingGenes = Utils.getOverlappingFeatures(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag());
				if(overlappingGenes.size() == 0) res.noFeature++;
				else if(overlappingGenes.size() == 1)	
				{
					res.mapped++;
					String gene = overlappingGenes.iterator().next();
					res.counts[Global.geneIndex.get(gene)]++;
				}
				else res.ambiguous++;
			}
		}
	}
}
