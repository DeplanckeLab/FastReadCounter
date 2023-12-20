package com.frc.bam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class Read implements Comparable<Read>
{
	public String name;
	public String barcode;
	public String qualityBC;
	public String UMI;
	public String qualityUMI;
	public String gene;
	public SAMRecord samRecord;
	public String start;
	public String end;
	public boolean barcodeMatch = false;
	
	public Cigar cigar;
	public int startV;
	public int endV;
	public String chr;
	public boolean negativeStrandFlag;
	public int mapQ;
	public boolean firstOfPair = false;
	
	public Read(String chr, int start, int end, Cigar c, boolean strand, int mapQ, String gene, String UMI, boolean firstOfPair) 
	{
		this.chr = chr;
		this.startV = start;
		this.endV = end;
		this.cigar = c;
		this.negativeStrandFlag = strand;
		this.mapQ = mapQ;
		this.gene = gene;
		this.UMI = UMI;
		this.firstOfPair = firstOfPair;
	}
	
	@Override
	public int compareTo(Read r2) 
	{
		return this.name.compareTo(r2.name);
	}
}
