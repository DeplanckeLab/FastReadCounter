package model;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import tools.Utils;

public class BAM 
{
	/**
	 * Using Picard to read the reads from the BAM file created by the alignment tool
	 */
	public static void readBAM()
	{
		// TODO Handle paired reads? Count once each read (current) or count pairs and discard single reads?
		// TODO Handle UMI? Default in UB tag (corrected) UB:Z:GGTGCTTGTTACA or UR tag (non corrected) UR:Z:GAGTCGCGCTCAG
		
		// Init result struct
		HashMap<String, ResultStruct> results = new HashMap<String, ResultStruct>();
		for(String barcode:Parameters.barcodes) results.put(barcode, new ResultStruct(Global.geneIndex.size()));
		int foundGXTag = 0;
		
		// Start reading BAM file
		Long start = System.currentTimeMillis();
		System.out.println("\nReading the reads from the BAM file...");
		try
		{
			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = samReaderFactory.open(Parameters.inputBAMFile);
			SAMRecordIterator it = samReader.iterator();
			
			// Start reading the BAM file
			while(it.hasNext())
			{
				SAMRecord samRecord = it.next();
				Global.nbReads++;
				
				// Resolve barcode first
				String barcode = "Unknown";
				if(Parameters.doDemultiplexing)
				{
					barcode = (String)samRecord.getAttribute("CB");
					if(barcode == null || barcode.equals("-")) // New in STAR, now unknown barcodes are labelled '-' instead of the CB tag being nonexistent.
					{
						barcode = "Unknown";
						Global.notDemultiplexed++;
					}
				}
				
				// Check if barcode is in the list of barcodes
				ResultStruct res = results.get(barcode); // This should exist
				if(res == null) new ErrorMessage("Barcode " + barcode + " is found in the BAM file, but is not in your list of barcodes");
				res.nbReads++; // Count reads by barcode
				
				// Quantify read according to its tags
				if(samRecord.getReadUnmappedFlag()) { Global.unmapped++; res.unmapped++; }
				else if(samRecord.getMappingQuality() < Parameters.minAQual) { Global.toolowAqual++; res.toolowAqual++; }
				else
				{
					boolean toProcess = true;
					if(samRecord.getSupplementaryAlignmentFlag())
					{
						Global.notUnique++; res.notUnique++;
						if(!Parameters.keep_multiple_mapped_reads) toProcess = false;
					}
					if(toProcess)
					{
						if(Parameters.use_bam_tags)
						{
							String gene_id = (String)samRecord.getAttribute("GX");
							//String gene_name = (String)samRecord.getAttribute("GN"); // Not used
							if(gene_id == null || gene_id.equals("-")) { Global.noFeature++; res.noFeature++; }
							else 
							{	
								foundGXTag++;
								Global.mapped++; res.mapped++;
								Integer indexGene = Global.geneIndex.get(gene_id); // From GTF
								if(indexGene == null) new ErrorMessage("ERROR: This gene " + gene_id + " is not in your GTF file. Please check again.");
								res.counts[indexGene]++;
							}
						}
						else // Use positions
						{
							HashSet<String> overlappingGenes = Utils.getOverlappingFeatures(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag());
							if(overlappingGenes.size() == 0) { Global.noFeature++; res.noFeature++; }
							else if(overlappingGenes.size() == 1)	
							{
								Global.mapped++; res.mapped++;
								String gene = overlappingGenes.iterator().next();
								res.counts[Global.geneIndex.get(gene)]++;
							}
							else { Global.ambiguous++; res.ambiguous++; }
						}
					}
					

				}
				
				if(Global.nbReads % 10000000 == 0) System.out.println(Global.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			samReader.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		System.out.println(Global.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]\n");
		
		System.out.println(Global.mapped + " 'Mapped' reads (" + Utils.pcFormatter.format(((float)Global.mapped / Global.nbReads) * 100) + "%)");
		System.out.println(Global.ambiguous + " 'Ambiguous' reads (" + Utils.pcFormatter.format(((float)Global.ambiguous / Global.nbReads) * 100) + "%)");
		System.out.println(Global.noFeature + " 'No Features' reads (" + Utils.pcFormatter.format(((float)Global.noFeature / Global.nbReads) * 100) + "%)");
		System.out.println(Global.notUnique + " 'Not Unique' reads (" + Utils.pcFormatter.format(((float)Global.notUnique / Global.nbReads) * 100) + "%)");
		System.out.println(Global.unmapped + " 'Not Aligned' reads (" + Utils.pcFormatter.format(((float)Global.unmapped / Global.nbReads) * 100) + "%)");
		System.out.println(Global.toolowAqual + " 'Too Low aQual' reads (" + Utils.pcFormatter.format(((float)Global.toolowAqual / Global.nbReads) * 100) + "%)");
		if(Parameters.doDemultiplexing) System.out.println(Global.notDemultiplexed + " 'Aligned but not demultiplexed' reads (" + Utils.pcFormatter.format(((float)Global.notDemultiplexed / Global.nbReads) * 100) + "%)");
		
		if(Parameters.use_bam_tags && foundGXTag == 0) System.err.println("[Warning] You specified the --bamtag option, but no GX tag was found in the BAM file. Is this really a gene annotated BAM?");
		
		// Create the read count matrix
		System.out.println("Writing output file...");
		createOutputDGE(results);
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
			bw_reads_detailed.write("__alignment_not_unique\t__alignment_not_unique");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.notUnique);
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
			bw_reads_detailed.write("__no_feature\t__no_feature");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.noFeature);
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
			bw_reads_detailed.write("__too_low_aQual\t__too_low_aQual");
			for(String barcode:sortedBarcodeKeys) 
			{
				ResultStruct res = results.get(barcode);
				bw_reads_detailed.write("\t" + res.toolowAqual);
			}
			bw_reads_detailed.write("\n");
			
			bw_reads.close(); 
			bw_reads_detailed.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		System.out.println("Count matrix written in " + Parameters.outputFolder);
	}

}
