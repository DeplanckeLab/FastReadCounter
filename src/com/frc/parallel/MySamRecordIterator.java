package com.frc.parallel;
/*
The MIT License (MIT)


Copyright (c) 2021 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

import java.io.InputStream;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;

/**
 * This class and the rest of this package was taken from the great implementation of jvarkit from Pierre Lindenbaum: https://github.com/lindenb
 *
 */
public class MySamRecordIterator extends AbstractIterator<SAMRecord> implements CloseableIterator<SAMRecord> 
{
	/** bam record codec converting input stream to sam record */
    private final BAMRecordCodec bamRecordCodec;
    
    /** cleanup on close */
    private Runnable onClose = null;
    
    public MySamRecordIterator(final SAMFileHeader header,final InputStream in,final String url) 
    {
    	this.bamRecordCodec = new BAMRecordCodec(header);
    	this.bamRecordCodec.setInputStream(in, url);
    }
    
	@Override
	protected SAMRecord advance() 
	{
        final SAMRecord next = this.bamRecordCodec.decode();
        if(next==null) close();
        return next;
	}
	
	@Override
	public void close() 
	{
		if(onClose!=null) onClose.run();
	}
}
