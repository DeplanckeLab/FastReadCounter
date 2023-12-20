package com.frc.exec;

import java.util.HashMap;

import com.frc.bam.BAM;
import com.frc.bam.Read;
import com.frc.parallel.CustomSamReader;
import com.frc.parameters.ResultStruct;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

public class Job extends Thread
{
	private CloseableIterator<SAMRecord> it;
	private CustomSamReader samReader;
	private SAMRecord lastSAMRecord;
	private HashMap<String, ResultStruct> results;
	private HashMap<String, Read> pairedBuffer;
	
	public long nbReads = 0;
	
	public Job(CloseableIterator<SAMRecord> it, CustomSamReader samReader, SAMRecord lastSAMRecord, HashMap<String, ResultStruct> results, HashMap<String, Read> pairedBuffer)
	{
		this.it = it;
		this.samReader = samReader;
		this.lastSAMRecord = lastSAMRecord;
		this.results = results;
		this.pairedBuffer = pairedBuffer;
	}
	
	public boolean reachedSameRead(SAMRecord samRecord, SAMRecord samRecord2) 
	{
		// Note, I could add this in the equals() inherited function of SAMRecord
		if(samRecord == null) return false;
		if(samRecord2 == null) return false;
		if(!samRecord.getReadName().equals(samRecord2.getReadName())) return false;
		// If I stop here, it can fail in case of:
		// - Bad BAM files containing duplicated reads all over the place
		// - Ok BAM files with simply duplicated reads?
		if(samRecord.getAlignmentEnd() != samRecord2.getAlignmentEnd()) return false;
		if(samRecord.getAlignmentStart() != samRecord2.getAlignmentStart()) return false;
		// How could I know otherwise if it's exactly the same read? Should take the offset loci maybe?
		return true;
	}
		
	@Override
	public void run() 
	{
		// Start reading the BAM file
		while(it.hasNext())
		{
			SAMRecord samRecord = it.next();
			
			nbReads++;

			BAM.processRecord(samRecord, this.results, this.pairedBuffer);
			
			// TODO: OMG we should be careful with paired-end seq. Multiple mapped reads?
			if(reachedSameRead(samRecord, lastSAMRecord)) break; // I need to process it, so I break at the end
		}
		
		samReader.close();
	}
	
	public HashMap<String, ResultStruct> getResults()
	{
		return results;
	}
}
