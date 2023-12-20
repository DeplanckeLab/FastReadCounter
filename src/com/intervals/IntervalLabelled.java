package com.intervals;

import htsjdk.tribble.index.interval.Interval;

public class IntervalLabelled extends Interval
{
	public String featureName;
	public boolean readNegativeStrandFlag;
	
	public IntervalLabelled(int start, int end, String featureName, boolean readNegativeStrandFlag) 
	{
		super(start, end);
		this.featureName = featureName;
		this.readNegativeStrandFlag = readNegativeStrandFlag;
	}
	
	@Override
	public String toString() {
		return super.toString() + " : " + featureName;
	}
}
