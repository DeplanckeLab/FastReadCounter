package com.frc.exec;

import java.io.IOException;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.HashMap;

import com.errors.ErrorMessage;
import com.frc.bam.BAM;
import com.frc.bam.Read;
import com.frc.parallel.CustomSamReader;
import com.frc.parameters.Global;
import com.frc.parameters.Parameters;
import com.frc.parameters.ResultStruct;
import com.tools.MemoryHandler;
import com.tools.Utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

public class JobDispatcher
{
	private static boolean heapMsg = false;
	private static HashMap<String, ResultStruct> results = null;
	private static HashMap<String, Read> pairedBuffer = null;
	
	private static void init()
	{
		// Init Main result struct
		results = new HashMap<String, ResultStruct>();
		for(String barcode:Parameters.barcodes) results.put(barcode, new ResultStruct(Global.geneIndex.size(), barcode));

		// Init temporary paired read buffer
		if(Parameters.is_paired) pairedBuffer = new HashMap<String, Read>();
	}
	
	public static HashMap<String, ResultStruct> readBAM()
	{
		init();
		Long start = System.currentTimeMillis();
		Long previous = start;
		
		// Reading file as one thread
		long nbReads = 0; // processed reads
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader samReader = samReaderFactory.open(Parameters.inputBAMFile);
		SAMRecordIterator it = samReader.iterator();
		
		// Start reading BAM file
		System.out.println("\nReading the reads from the BAM file provided: " + Parameters.inputBAMFile);
		while(it.hasNext())
		{
			SAMRecord samRecord = it.next();
			nbReads++;
			BAM.processRecord(samRecord, results, pairedBuffer);
			
			if(nbReads % 10000000 == 0) {System.out.println(nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "], Diff = " + Utils.toReadableTime(System.currentTimeMillis() - previous)); previous = System.currentTimeMillis();}
			if(nbReads % 100000 == 0)
			{
				if(MemoryHandler.isHeapFull() && !heapMsg)
				{
					System.err.println("[Warning] Java heap space is almost full! Current used Heap size = " + Utils.pcFormatter.format(MemoryHandler.getCurrentHeapSizePercent()) + "% of Total " + MemoryHandler.toReadableSize(MemoryHandler.maxHeap));
					System.err.println("[Warning] Program may crash unexpectedly or slow down drastically. Please consider rerunning the tool increasing the heap space. For e.g. to use 16Gb of heap space, use the following code:");
					System.err.println("[Warning]\tjava -Xmx16g -jar FastReadCounter.jar ....");
					heapMsg = true; // Avoid spamming the user...
				}
			}
		}
		try 
		{
			samReader.close();
		} 
		catch (IOException ioe) 
		{
			new ErrorMessage(ioe.getMessage());
		}
		System.out.println(nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		printResultSummary();
		return results;
	}
	
	public static HashMap<String, ResultStruct> readBAMParallel()
	{
		init();
		Long start = System.currentTimeMillis();
		
		long nbReads = 0; // processed reads
		ArrayList<Job> myjobs = new ArrayList<Job>(); // All threads
		String filePath = Parameters.inputBAMFile.getAbsolutePath();
		
		long size = Utils.getFileSize(filePath);
		long offset = size - size / Parameters.nbThreads; // Start by last block to know to which Read the thread needs to stop
		SAMRecord lastSAMRecord = null;
		
		// First, create the job pool
		for(int n_job = 1; n_job <= Parameters.nbThreads; n_job++)
		{
			//System.out.println("Creating Thread " + n_job);
			try
			{
				CustomSamReader mySAMReader = new CustomSamReader(filePath);
				if(n_job == Parameters.nbThreads) offset = 0; // Start at first read
				CloseableIterator<SAMRecord> it = mySAMReader.queryAtOffset(offset);
				
				if(it.hasNext()) 
				{
					SAMRecord samRecord = null; 
					if(n_job != Parameters.nbThreads) samRecord = it.next(); // Not for the first
					
					// Create part result struct
					HashMap<String, ResultStruct> result_tmp = new HashMap<String, ResultStruct>();
					for(String barcode:Parameters.barcodes) result_tmp.put(barcode, new ResultStruct(Global.geneIndex.size(), barcode));
					
					Job thread = new Job(it, mySAMReader, lastSAMRecord, result_tmp, pairedBuffer);
					thread.setDaemon(true);
					thread.start();
					
					lastSAMRecord = samRecord; // For next thread
					
					offset = offset - size / Parameters.nbThreads;
					myjobs.add(thread);
				}
			}
			catch(IOException ioe)
			{
				new ErrorMessage(ioe.getMessage());
			}
		}
		
		// Create Mama thread, waiting for the other to finish, and compiling the results
		while(true)
		{
			int nbFinished = 0;
			nbReads = 0;
			for(Job j:myjobs)
			{
				nbReads += j.nbReads;
				if(j.getState() == State.TERMINATED) nbFinished++;
			}
			if(nbFinished == myjobs.size()) break;
			System.out.println(nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			try
			{
				Thread.sleep(1000);
			}
			catch(InterruptedException ie)
			{
				new ErrorMessage(ie.getMessage());
			}
		}
		
		// Merge all results
		for(Job j:myjobs)
		{
			HashMap<String, ResultStruct> result_tmp = j.getResults();
			for(String barcode:result_tmp.keySet())
			{
				ResultStruct res = result_tmp.get(barcode);
				results.get(barcode).add(res);
			}
		}
		
		System.out.println(nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		printResultSummary();
		return results;
	}
	
	public static void printResultSummary()
	{
		ResultStruct summary = new ResultStruct();
		for(String barcode:results.keySet())
		{
			ResultStruct res = results.get(barcode);
			summary.addSummary(res);
		}	
		
		long total = summary.nbReads;
		
		System.out.println("\n[Global information]");
		System.out.println(summary.notUnique + " reads/pairs have secondary/supplement tag [multiple mappers] (" + Utils.pcFormatter.format(((float)summary.notUnique / total) * 100) + "%)");
		System.out.println(summary.duplicates + " reads/pairs are marked as duplicates (" + Utils.pcFormatter.format(((float)summary.duplicates / total) * 100) + "%)");
		
		if(Parameters.is_paired) 
		{
			System.out.println("\n[Paired-end specific]");
			System.out.println("In total, we retained " + (summary.countedPair + summary.countedUnique) + " pairs of reads and singletons (" + (summary.countedPair * 2 + summary.countedUnique) + " reads i.e. " + Utils.pcFormatter.format(((float)(summary.countedPair * 2 + summary.countedUnique) / total) * 100) + "% of total)");
			
			total = summary.countedPair + summary.countedUnique;
			// summary.mateUnmapped == summary.countedUnique
			//System.out.println(summary.mateUnmapped + " reads don't have aligned paired mate [singletons] (" + Utils.pcFormatter.format(((float)summary.mateUnmapped / total) * 100) + "%)");
			// summary.properlyPaired : I don't know if it is useful here?
			// Properly paired means that the read orientation goes F> ---gap --- < R and that the gap is roughly the expected size, 200-400 bp
			//System.out.println(summary.properlyPaired + " reads are properly paired (" + Utils.pcFormatter.format(((float)summary.properlyPaired / total) * 100) + "%)");
			System.out.println("\t" + summary.countedPair + " PAIRS of read (" + summary.countedPair * 2 + " reads = " + Utils.pcFormatter.format(((float)(summary.countedPair * 2) / (summary.nbReads - summary.notUnique)) * 100) + "% of retained reads)");
			System.out.println("\t" + summary.countedUnique + " SINGLETON reads (mate is not aligned) (" + Utils.pcFormatter.format(((float)summary.mateUnmapped / (summary.nbReads - summary.notUnique)) * 100) + "% of retained reads)");
		}
		
		System.out.println("\n[Summary]");
		System.out.println(summary.mapped + " reads/pairs pass all QCs and are included in the count table (" + Utils.pcFormatter.format(((float)summary.mapped / total) * 100) + "%)");
		System.out.println(summary.unmapped + " reads/pairs are not mapped [Unmapped] (" + Utils.pcFormatter.format(((float)summary.unmapped / total) * 100) + "%)");
		System.out.println(summary.ambiguous + " reads/pairs are aligned but map to multiple features [ambiguous] (" + Utils.pcFormatter.format(((float)summary.ambiguous / total) * 100) + "%)");
		System.out.println(summary.noFeature + " reads/pairs are aligned but do not map to any feature [no-feature] (" + Utils.pcFormatter.format(((float)summary.noFeature / total) * 100) + "%)");
		System.out.println(summary.toolowAqual + " reads/pairs have too low alignment quality [too low aQual] (" + Utils.pcFormatter.format(((float)summary.toolowAqual / total) * 100) + "%)");
		if(Parameters.doDemultiplexing) System.out.println(results.get("Unknown").nbReads + " 'Not demultiplexed' reads (" + Utils.pcFormatter.format(((float)results.get("Unknown").nbReads / total) * 100) + "%)");
		if(Parameters.use_bam_tags && summary.foundGXTag == 0) System.err.println("[Warning] You specified the --bamtag option, but no GX tag was found in the BAM file. Is this really a gene annotated BAM?");
		if(Parameters.use_bam_tags && summary.foundUTag == 0) System.err.println("[Warning] You specified the --bamtag option, but no UB or UR tag was found in the BAM file. The UMI count matrix will NOT be generated.");
	}
}
