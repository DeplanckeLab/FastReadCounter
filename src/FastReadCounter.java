import java.util.HashMap;

import com.frc.bed.BED;
import com.frc.exec.JobDispatcher;
import com.frc.gtf.GTF;
import com.frc.parameters.Parameters;
import com.frc.parameters.ResultStruct;
import com.frc.vcf.VCF;
import com.tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class FastReadCounter
{	
	// TODO Add the possibility to process multiple bams at once
	// TODO Add parallelization (currently, bug, don't generate the exact same output, problem in sync?)
	// TODO Don't require a GTF file for --bamtag option
	// TODO Remove "Unkown" column when not detailed => Careful when only one sample (it's called Unkown)
	// TODO Add prefix for output, instead of folder
	// TODO Handle paired / multiple reads for paralellization (pairedBuffer should be created smartly)
	// TODO check merge UMI (pas fait)
	
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
			HashMap<String, ResultStruct> results;
			if(Parameters.nbThreads == 1) results = JobDispatcher.readBAM();
			else results = JobDispatcher.readBAMParallel();
			System.out.println("BAM reading DONE [" + Utils.toReadableTime(System.currentTimeMillis() - time) + "]");
			
			// Create the read count matrix
			time = System.currentTimeMillis();
			ResultStruct.createOutputDGE(results);
			System.out.println("Output file writing DONE [" + Utils.toReadableTime(System.currentTimeMillis() - time) + "]");
			
			System.out.println("\nFRC DONE [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		}
	}
}

