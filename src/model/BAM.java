package model;

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
	public static HashMap<String, ResultStruct> readBAM()
	{	
		// Init result struct
		HashMap<String, ResultStruct> results = new HashMap<String, ResultStruct>();
		for(String barcode:Parameters.barcodes) results.put(barcode, new ResultStruct(Global.geneIndex.size()));

		// Init temporary paired read buffer
		HashMap<String, Read> pairedBuffer = null;
		if(Parameters.is_paired) pairedBuffer = new HashMap<String, Read>();
		
		// Start reading BAM file
		Long start = System.currentTimeMillis();
		System.out.println("\nReading the reads from the BAM file provided: " + Parameters.inputBAMFile);
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
				if(Parameters.is_paired && !samRecord.getReadPairedFlag()) new ErrorMessage("You used the --paired option, but it seems that some reads are not paired?");
				if(samRecord.getDuplicateReadFlag()) Global.duplicates++; // Just count this
				if(samRecord.getReadUnmappedFlag()) { Global.unmapped++; res.unmapped++; }
				else
				{
					boolean toProcess = true;
					if(samRecord.isSecondaryOrSupplementary())
					{
						Global.notUnique++; res.notUnique++;
						if(!Parameters.keep_multiple_mapped_reads) toProcess = false;
					}
					if(toProcess)
					{
						// Paired-end
						if(Parameters.is_paired)
						{
							// If MateUnmapped flag, then the read cannot have the "ProperPair" flag
							if(samRecord.getMateUnmappedFlag()) { Global.mateUnmapped++; res.mateUnmapped++; }
							if(samRecord.getProperPairFlag()) { Global.properlyPaired++; res.properlyPaired++; }
							if(samRecord.getMateUnmappedFlag())
							{
								Global.countedUnique++;
								// Process this singleton read as if it was single-end (not stored)
								processSingleEndRead(samRecord, res);
							}
							else
							{
								String name = samRecord.getReadName();
								Read read1 = pairedBuffer.remove(samRecord.getReadName());
								if(read1 != null) // I run the counting only if two pairs are identified/available
								{
									Global.countedPair++;
									if(samRecord.getMappingQuality() < Parameters.minAQual || read1.mapQ < Parameters.minAQual) { Global.toolowAqual++; res.toolowAqual++; } // To match htseq-count. I do this first
									else
									{
										if(Parameters.use_bam_tags)
										{
											String gene_id_1 = read1.gene;
											String gene_id_2 = (String)samRecord.getAttribute("GX");
											if(gene_id_1 != null && gene_id_1.equals("-")) gene_id_1 = null;
											if(gene_id_2 != null && gene_id_2.equals("-")) gene_id_2 = null;
											if(gene_id_1 != null && gene_id_2 != null && !gene_id_1.equals(gene_id_2)) { Global.foundGXTag++; Global.ambiguous++; res.ambiguous++; }
											else if(gene_id_1 == null && gene_id_2 == null) { Global.noFeature++; res.noFeature++; }
											else
											{
												Global.foundGXTag++;
												String gene_id = gene_id_1;
												if(gene_id == null) gene_id = gene_id_2;
												Global.mapped = Global.mapped + 2; res.mapped = res.mapped + 2;
												Integer indexGene = Global.geneIndex.get(gene_id); // From GTF
												if(indexGene == null) new ErrorMessage("ERROR: This gene " + gene_id + " is not in your GTF file. Please check again.");
												res.counts[indexGene]++;
											}
										}
										else // Use positions
										{
											HashSet<String> overlappingGenes_1 = Utils.getOverlappingFeatures(read1.chr, read1.startV, read1.endV, read1.cigar, read1.negativeStrandFlag, read1.firstOfPair);
											HashSet<String> overlappingGenes_2 = Utils.getOverlappingFeatures(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag(), samRecord.getFirstOfPairFlag());
											
											if(overlappingGenes_1.size() == 0 && overlappingGenes_2.size() == 0) { Global.noFeature++; res.noFeature++; }
											else
											{
												overlappingGenes_1.addAll(overlappingGenes_2);
												if(overlappingGenes_1.size() == 1)	
												{
													Global.mapped++; res.mapped++;
													String gene = overlappingGenes_1.iterator().next();
													res.counts[Global.geneIndex.get(gene)]++;
												}
												else { Global.ambiguous++; res.ambiguous++;  }
											}
										}
									}
								}
								else 
								{
									pairedBuffer.put(name, new Read(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag(), samRecord.getMappingQuality(), (String)samRecord.getAttribute("GX"), samRecord.getFirstOfPairFlag()));
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
				
				if(Global.nbReads % 10000000 == 0) System.out.println(Global.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			samReader.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		System.out.println(Global.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		long total = Global.nbReads;
		
		System.out.println("\n[Global information]");
		System.out.println(Global.notUnique + " reads/pairs have secondary/supplement tag [multiple mappers] (" + Utils.pcFormatter.format(((float)Global.notUnique / total) * 100) + "%)");
		System.out.println(Global.duplicates + " reads/pairs are marked as duplicates (" + Utils.pcFormatter.format(((float)Global.duplicates / total) * 100) + "%)");
		
		if(Parameters.is_paired) 
		{
			System.out.println("\n[Paired-end specific]");
			System.out.println("In total, we retained " + (Global.countedPair + Global.countedUnique) + " pairs of reads and singletons (" + (Global.countedPair * 2 + Global.countedUnique) + " reads i.e. " + Utils.pcFormatter.format(((float)(Global.countedPair * 2 + Global.countedUnique) / total) * 100) + "% of total)");
			
			total = Global.countedPair + Global.countedUnique;
			// Global.mateUnmapped == Global.countedUnique
			//System.out.println(Global.mateUnmapped + " reads don't have aligned paired mate [singletons] (" + Utils.pcFormatter.format(((float)Global.mateUnmapped / total) * 100) + "%)");
			// Global.properlyPaired : I don't know if it is useful here?
			// Properly paired means that the read orientation goes F> ---gap --- < R and that the gap is roughly the expected size, 200-400 bp
			//System.out.println(Global.properlyPaired + " reads are properly paired (" + Utils.pcFormatter.format(((float)Global.properlyPaired / total) * 100) + "%)");
			System.out.println("\t" + Global.countedPair + " PAIRS of read (" + Global.countedPair * 2 + " reads = " + Utils.pcFormatter.format(((float)(Global.countedPair * 2) / (Global.nbReads - Global.notUnique)) * 100) + "% of retained reads)");
			System.out.println("\t" + Global.countedUnique + " SINGLETON reads (mate is not aligned) (" + Utils.pcFormatter.format(((float)Global.mateUnmapped / (Global.nbReads - Global.notUnique)) * 100) + "% of retained reads)");
		}
		
		System.out.println("\n[Summary]");
		System.out.println(Global.mapped + " reads/pairs pass all QCs and are included in the count table (" + Utils.pcFormatter.format(((float)Global.mapped / total) * 100) + "%)");
		System.out.println(Global.unmapped + " reads/pairs are not mapped [Unmapped] (" + Utils.pcFormatter.format(((float)Global.unmapped / total) * 100) + "%)");
		System.out.println(Global.ambiguous + " reads/pairs are aligned but map to multiple features [ambiguous] (" + Utils.pcFormatter.format(((float)Global.ambiguous / total) * 100) + "%)");
		System.out.println(Global.noFeature + " reads/pairs are aligned but do not map to any feature [no-feature] (" + Utils.pcFormatter.format(((float)Global.noFeature / total) * 100) + "%)");
		System.out.println(Global.toolowAqual + " reads/pairs have too low alignment quality [too low aQual] (" + Utils.pcFormatter.format(((float)Global.toolowAqual / total) * 100) + "%)");
		if(Parameters.doDemultiplexing) System.out.println(Global.notDemultiplexed + " 'Aligned but not demultiplexed' reads (" + Utils.pcFormatter.format(((float)Global.notDemultiplexed / total) * 100) + "%)");
		if(Parameters.use_bam_tags && Global.foundGXTag == 0) System.err.println("[Warning] You specified the --bamtag option, but no GX tag was found in the BAM file. Is this really a gene annotated BAM?");
	
		return results;
	}
	
	private static void processSingleEndRead(SAMRecord samRecord, ResultStruct res)
	{
		if(samRecord.getMappingQuality() < Parameters.minAQual) { Global.toolowAqual++; res.toolowAqual++; } // To match htseq-count. I do this first
		else
		{	
			if(Parameters.use_bam_tags)
			{
				String gene_id = (String)samRecord.getAttribute("GX");
				//String gene_name = (String)samRecord.getAttribute("GN"); // Not used
				if(gene_id == null || gene_id.equals("-")) { Global.noFeature++; res.noFeature++; }
				else 
				{	
					Global.foundGXTag++;
					Global.mapped++; res.mapped++;
					Integer indexGene = Global.geneIndex.get(gene_id); // From GTF
					if(indexGene == null) new ErrorMessage("ERROR: This gene " + gene_id + " is not in your GTF file. Please check again.");
					res.counts[indexGene]++;
				}
			}
			else // Use positions
			{
				HashSet<String> overlappingGenes = Utils.getOverlappingFeatures(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag(), samRecord.getFirstOfPairFlag());
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
}
