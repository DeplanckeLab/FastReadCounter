package model;

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
	
	@Override
	public int compareTo(Read r2) 
	{
		return this.name.compareTo(r2.name);
	}
}
