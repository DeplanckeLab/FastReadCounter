package com.frc.parameters;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

import com.errors.ErrorMessage;
import com.errors.WarningMessage;
import com.tools.Utils;

import htsjdk.tribble.annotation.Strand;

public class Parameters 
{
	public static final String currentVersion = "1.0";
	
	// Input parameters
	public static String outputFolder = null;
	public static File inputBAMFile = null;
	public static File inputGTFFile = null;
	public static File inputVCFFile = null;
	public static File inputBEDFile = null;
	public static Strand stranded = Strand.NONE;
	public static int minAQual = 10;
	public static HashSet<String> barcodes = null;
	public static String barcodeTag = "CB"; // Default barcode-corrected tag, CR is not corrected
	public static int barcodeLength = -1; // Default barcode-corrected tag, CR is not corrected
	public static boolean doDemultiplexing = false;
	public static boolean keep_multiple_mapped_reads = false;
	public static boolean use_bam_tags = false;
	public static boolean is_paired = false;
	public static UMIDedup umi_dedup = UMIDedup.NONE;
	public static int nbThreads = 1;
	
	public static void load(String[] args) throws Exception
	{
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) new ErrorMessage(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "--stranded":
					case "-s":
						i++;
						switch(args[i])
						{
						case "no":
							Parameters.stranded = Strand.NONE;
							break;
						case "yes":
							Parameters.stranded = Strand.FORWARD;
							break;
						case "reverse":
							Parameters.stranded = Strand.REVERSE;
							break;
						default:
							new ErrorMessage("The '-s' | '--stranded' option should be followed by one of the following parameters: [no, yes, reverse].");
						}
						break;
					case "--paired":
					case "-p":
						is_paired = true;
						break;
					case "--bam":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputBAMFile = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--bam' option should be followed by aligned BAM file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--gtf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputGTFFile = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--gtf' option should be followed by GTF file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--bamtag":
						use_bam_tags = true;
						break;
					case "--umi-dedup":
						i++;
						switch(args[i])
						{
						case "none":
							Parameters.umi_dedup =  UMIDedup.NONE;
							break;
						case "exact":
							Parameters.umi_dedup =  UMIDedup.EXACT;
							break;
						default:
							new ErrorMessage("The '--umi-dedup' option should be followed by one of the following parameters: [none, exact].");
						}
						break;
					case "--bed":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputBEDFile = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--bed' option should be followed by BED file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--vcf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputVCFFile = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '-vcf' option should be followed by VCF file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--barcodeFile":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							barcodes = Utils.readBarcodeFile(args[i]);
							doDemultiplexing = true;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--barcodeFile' option should be followed by a barcode file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--barcodeTag":
						i++;
						barcodeTag = args[i];
						break;
					case "-q":
						i++;
						try
						{
							minAQual = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorMessage("The '-q' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-t":
					case "--threads":
						i++;
						try
						{
							nbThreads = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorMessage("The '-t' | '--threads' option should be followed by an Integer. You entered " + args[i]);
						}
						if(nbThreads < 1) new ErrorMessage("The '-t' | '--threads' option should be followed by a positive Integer. You entered " + nbThreads);
						int nbLogicalThreads = Utils.getNbLogicalThreads();
						if(nbThreads >= nbLogicalThreads) new WarningMessage("There are " + nbLogicalThreads + " available logical threads on your machine. We would not recommend using more than " + (nbLogicalThreads - 1) + " (here you've set " + nbThreads + "). Optimal/maximum number should be ~nbLogicalThreads/2 = " + nbLogicalThreads / 2);
						break;
					case "--multiple-mapped":
						keep_multiple_mapped_reads = true;
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(inputBAMFile == null)
		{
			new ErrorMessage("Please use '--bam' option to specify the path of the aligned BAM file");
		}
		int c = 0;
		if(inputGTFFile != null) c++;
		if(inputBEDFile != null) c++;
		if(inputVCFFile != null) c++;;
		if(c == 0) new ErrorMessage("Please use '--gtf', '--bed', or '--vcf' option to specify how to count the features");
		if(c > 1) new ErrorMessage("You specified several features ('--gtf' '--bed' '--vcf'), while you should specify only ONE");
		if(inputVCFFile != null) 
		{
			new WarningMessage("VCF mode, [NOT Stranded] only.");
			Parameters.stranded = Strand.NONE;
		}
		System.out.println("[Parameters]");
		if(!use_bam_tags) 
		{
			if(Parameters.umi_dedup !=  UMIDedup.NONE) new ErrorMessage("You cannot specify --umi-dedup option without the --bamtag option, i.e. FRC can only count UMIs that are already annotated in the BAM file.");
			System.out.println("Stranded = " + Parameters.stranded);
		}
		else 
		{
			if(inputGTFFile == null) new ErrorMessage("The option --bamtag can only be used with option --gtf");
			System.out.println("Using BAM tags GN/GX for computing count matrix.");
			if(Parameters.umi_dedup !=  UMIDedup.NONE) System.out.println("UMI count matrix WILL BE generated using the " + Parameters.umi_dedup + " deduplication scheme");
			else System.out.println("UMI count matrix WILL NOT BE generated (use '--umi-dedup' option to do so)");
		}
		
		System.out.println("Min aQual = " + Parameters.minAQual);
		System.out.println("Multiple mapped reads will be " + ((keep_multiple_mapped_reads)?"KEPT":"DISCARDED"));
		if(Parameters.doDemultiplexing)
		{
			System.out.println("Samples WILL BE demultiplexed:");
			if(Parameters.barcodes.size() == 0) new ErrorMessage("\tNo barcodes were found in your barcode file.");
			System.out.println("\tBarcode file = " + Parameters.barcodes.size() + " barcodes of length " + Parameters.barcodeLength);
			System.out.println("\tBarcode tag = " + Parameters.barcodeTag);
		}
		else 
		{
			Global.mappingBarcodeName = new HashMap<String, String>();
			Parameters.barcodes = new HashSet<String>();
		}
		Parameters.barcodes.add("Unknown");
		Global.mappingBarcodeName.put("Unknown", "Unknown");
		
		System.out.println("Samples will be treated as " + (Parameters.is_paired?"PAIRED-END":"SINGLE-END"));
		
		if(outputFolder == null)
		{
			String path = inputBAMFile.getAbsolutePath();
			path = path.replaceAll("\\\\", "/");
			path = path.substring(0, path.lastIndexOf("/"));
			new WarningMessage("No output folder is specified, using the BAM folder as default. You can specify an output folder by using the '-o' option.");
			outputFolder = path;
		}
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
		System.out.println("Output folder = " + Parameters.outputFolder);
		if(nbThreads == 1) System.out.println("Run is NOT parallelized: 1 thread is used");
		else System.out.println("Run IS parallelized: " + nbThreads + " threads are used");
	}
	
	public static void printHelp()
	{
		System.out.println("FastReadCounter (FRC) " + Parameters.currentVersion + "\n\nOptions:");
		System.out.println("\t--bam %s \t\tPath of BAM file (do not need to be sorted or indexed)");
		System.out.println("\t--gtf %s \t\tPath of GTF file");
		System.out.println("\t--bed %s \t\tPath of BED file");
		System.out.println("\t--vcf %s \t\tPath of VCF file");
		System.out.println("\t--bamtag \t\t[Use with --gtf] If BAM file is already mapped to genes (STARsolo) and/or tagged with GN/GX tags, use these tags to produce count matrix instead of the positions in the GTF file");
		System.out.println("\t--barcodeFile %s \tPath of Barcode file for demultiplexing the BAM file using the barcode tag (default tag is 'CB')");
		System.out.println("\t--barcodeTag %s \tChanging the barcode tag to check for the --barcodeFile option (default tag is 'CB')");
		System.out.println("\t--multiple-mapped %s \tKeep multiple mapped read (default = discard them)");
		System.out.println("\t--umi-dedup %s \t\t[none, exact] How to dedup the UMIs (only if using --bamtag option, default = none)");
		System.out.println("\t-s %s \t\t\t[no, yes, reverse] Do you want to count only reads falling on same strand than feature? (default = no, unused with option --bamtag)");
		System.out.println("\t-q %i \t\t\tMinimum quality required for a read to be counted (default = 10)");
		System.out.println("\t-t | --threads %i \tNumber of threads to use (default = 1)");
		System.out.println("\t-o %s \t\t\tOutput folder (default = folder of BAM file)");
	}
}
