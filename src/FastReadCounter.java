import java.util.HashMap;

import model.BAM;
import model.BED;
import model.GTF;
import model.Parameters;
import model.ResultStruct;
import model.VCF;
import tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class FastReadCounter
{	
	// TODO If BAM is sorted, reduce RAM for UMI, by processing each gene separately
	// TODO For multiple reads. If the same read maps multiple times to the same gene, count it only once.
	// TODO Add the possibility to process multiple bams at once
	// TODO Add parallelization
	
	public static void main(String[] args) throws Exception
	{
		if(args.length < 1) Parameters.printHelp();
		else
		{
			Long start = System.currentTimeMillis();
			System.out.println("FastReadCounter (FRC) " + Parameters.currentVersion + "\n");
			
			// Reading args
			Parameters.load(args);
				
			// Reading regions
			Long time = System.currentTimeMillis();
			if(Parameters.inputGTFFile != null) GTF.readGTF();
			else if(Parameters.inputVCFFile != null) VCF.readVCF();
			else BED.readBED(); // At least one is not null
			System.out.println("Region reading DONE [" + Utils.toReadableTime(System.currentTimeMillis() - time) + "]");
						
			// Reading reads in BAM
			time = System.currentTimeMillis();
			HashMap<String, ResultStruct> results = BAM.readBAM();
			System.out.println("BAM reading DONE [" + Utils.toReadableTime(System.currentTimeMillis() - time) + "]");
			
			// Create the read count matrix
			time = System.currentTimeMillis();
			ResultStruct.createOutputDGE(results);
			System.out.println("Output file writing DONE [" + Utils.toReadableTime(System.currentTimeMillis() - time) + "]");
			
			System.out.println("\nFRC DONE [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		}
	}
}

