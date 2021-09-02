package model;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

import tools.Utils;

enum Strand{NO, YES, REVERSE};

public class Parameters 
{
	public static final String currentVersion = "1.0";
	
	// Input parameters
	public static String outputFolder = null;
	public static File inputBAMFile = null;
	public static File inputGTFFile = null;
	public static File inputVCFFile = null;
	public static File inputBEDFile = null;
	public static Strand stranded = Strand.NO;
	public static int minAQual = 10;
	public static HashSet<String> barcodes = null;
	public static String barcodeTag = "CB"; // Default barcode-corrected tag, CR is not corrected
	public static int barcodeLength = -1; // Default barcode-corrected tag, CR is not corrected
	public static boolean doDemultiplexing = false;
	public static boolean keep_multiple_mapped_reads = false;
	public static boolean use_bam_tags = false;
	public static boolean is_paired = false;
	
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
							Parameters.stranded = Strand.NO;
							break;
						case "yes":
							Parameters.stranded = Strand.YES;
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
			System.err.println("[Warning] VCF mode, [NOT Stranded] only.");
			Parameters.stranded = Strand.NO;
		}
		System.out.println("[Parameters]");
		if(!use_bam_tags) System.out.println("Stranded = " + Parameters.stranded);
		else 
		{
			if(inputGTFFile == null) new ErrorMessage("The option --bamtag can only be used with option --gtf");
			System.out.println("Using BAM tags GN/GX for computing count matrix.");
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
			System.err.println("[Warning] No output folder is specified, using the BAM folder as default. You can specify an output folder by using the '-o' option.");
			outputFolder = path;
		}
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
		System.out.println("Output folder = " + Parameters.outputFolder);
	}
	
	public static void printHelp()
	{
		System.out.println("FastReadCounter (FRC) " + Parameters.currentVersion + "\n\nOptions:");
		System.out.println("\t--bam %s \t[Required] Path of BAM file [do not need to be sorted or indexed].");
		System.out.println("\t--gtf %s \t[Required (or --bed or --vcf or --bamtag)] Path of GTF file.");
		System.out.println("\t--bed %s \t[Required (or --gtf or --vcf or --bamtag)] Path of BED file.");
		System.out.println("\t--vcf %s \t[Required (or --gtf or --bed or --bamtag)] Path of VCF file.");
		System.out.println("\t--bamtag \t[Optional (with --gtf)] BAM file is already mapped to genes (STARsolo) and/or tagged with GN/GX tags, use these tags to produce count matrix instead of the positions in the GTF file.");
		System.out.println("\t--barcodeFile %s \t[Optional] Path of Barcode file for demultiplexing the BAM file using the barcode tag (default tag is 'CB').");
		System.out.println("\t--barcodeTag %s \t[Optional] Changing the barcode tag to check for the --barcodeFile option (default tag is 'CB').");
		System.out.println("\t--multiple-mapped %s \tKeep multiple mapped read [default = discard them]");
		System.out.println("\t-s %s \t\t[no, yes, reverse] Do you want to count only reads falling on same strand than feature? [default = no].");
		System.out.println("\t-q %i \t\tMinimum quality required for a read to be counted [default = 10].");
		System.out.println("\t-o %s \t\tOutput folder [default = folder of BAM file]");
	}
}
