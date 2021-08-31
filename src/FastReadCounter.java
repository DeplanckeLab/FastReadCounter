import model.BAM;
import model.BED;
import model.GTF;
import model.Parameters;
import model.VCF;
import tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class FastReadCounter
{	
	public static void main(String[] args) throws Exception
	{
		if(args.length < 1) Parameters.printHelp();
		else
		{
			System.out.println("FastReadCounter (FRC) " + Parameters.currentVersion + "\n");
			Parameters.load(args);
			Long time = System.currentTimeMillis();
			if(Parameters.inputGTFFile != null) GTF.readGTF();
			else if(Parameters.inputVCFFile != null) VCF.readVCF();
			else BED.readBED(); // At least one is not null
			System.out.println(Utils.toReadableTime(System.currentTimeMillis() - time));
			time = System.currentTimeMillis();
			BAM.readBAM();
			System.out.println(Utils.toReadableTime(System.currentTimeMillis() - time));
		}
	}
}

