package model;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;

public class Forest 
{
	private static HashMap<String, IntervalTree> forest = null;
	
	public static void init() 
	{
		forest = new HashMap<>();
	}
	
	public static void addToTree(String chr, IntervalLabelled interval) 
	{
		IntervalTree tree = forest.get(chr);
		if(tree == null) tree = new IntervalTree();
		tree.insert(interval);
		forest.put(chr, tree);
	}
	
	public static HashSet<String> findOverlappingFeatures(String chr, int start, int end, boolean readNegativeStrandFlag)
	{
		HashSet<String> result = new HashSet<>();
		IntervalTree itree = forest.get(chr);
		if(itree != null)
		{
			List<Interval> found = itree.findOverlapping(new Interval(start, end));
			for(Interval i:found) 
			{
				IntervalLabelled g = (IntervalLabelled)i;
				if(Parameters.stranded == Strand.NO) result.add(g.featureName);
				else if(Parameters.stranded == Strand.REVERSE && g.readNegativeStrandFlag == readNegativeStrandFlag) result.add(g.featureName);
				else if(Parameters.stranded == Strand.YES && g.readNegativeStrandFlag != readNegativeStrandFlag) result.add(g.featureName);
			}
		}
		return result;
	}
}
