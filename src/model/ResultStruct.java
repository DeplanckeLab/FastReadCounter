package model;

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
	
	public ResultStruct(int nbGenes)
	{
		this.counts = new int[nbGenes]; // genes 
	}
}
