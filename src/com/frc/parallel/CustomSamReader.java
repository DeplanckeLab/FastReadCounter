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

import java.io.Closeable;
import java.io.IOException;
import java.io.PushbackInputStream;
import java.nio.file.Paths;

import com.errors.ErrorMessage;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * This class and the rest of this package was taken from the great implementation of jvarkit from Pierre Lindenbaum: https://github.com/lindenb
 *
 */
/*public class CustomSamReader implements Closeable 
{
	private String path;
	private SAMFileHeader samFileHeader;
	private SeekableStream seekableStream;
	private SamReader samReader ;
	private BgzfBlockGuesser bgzfBlockGuesser; // class finding the next BGZF block
	public SortOrder sortOrder;
	
	public CustomSamReader(String path) throws IOException 
	{
		this.path = path;
		this.seekableStream = SeekableStreamFactory.getInstance().getStreamFor(Paths.get(path).toString());
		this.bgzfBlockGuesser= new BgzfBlockGuesser(this.seekableStream, this.path.toString());
		this.samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(this.seekableStream)); 
		this.samFileHeader = this.samReader.getFileHeader();
		this.sortOrder = this.samFileHeader.getSortOrder();
	}
	
	public SAMFileHeader getFileHeader() 
	{
		return samFileHeader;
	}

	public CloseableIterator<SAMRecord> queryAtOffset(final long start_offset) throws IOException 
	{
		if(start_offset >= this.seekableStream.length()) new ErrorMessage("Cannot create a reading thread with this offset (go beyond last file byte).");
		
		BgzfBlockGuesser.BgzfBlock bgzfBlock = this.bgzfBlockGuesser.guessNextBGZFPos(start_offset, seekableStream.length());
		if(bgzfBlock==null) new ErrorMessage("No bgz block found at "+start_offset);
		
		this.seekableStream.seek(bgzfBlock.pos);

		BlockCompressedInputStream bcis = new BlockCompressedInputStream(this.seekableStream);
		PushbackInputStream bpi = new PushbackInputStream(bcis, BamRecordGuesser.BUFFER_SIZE);
		BamRecordGuesser bamRecordGuesser = new BamRecordGuesser(this.samFileHeader);
		if(bamRecordGuesser.find(bpi)) return new MySamRecordIterator(this.samFileHeader, bpi, this.path.toString());
		return null;
	}
	
	@Override
	public void close() 
	{
		CloserUtil.close(this.seekableStream);
		this.bgzfBlockGuesser.close();
		try
		{
			this.samReader.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
	}
}*/

public class CustomSamReader implements Closeable 
{
	private String path;
	private SAMFileHeader samFileHeader;
	private SeekableStream seekableStream;
	private SamReader samReader ;
	private BgzfBlockGuesser bgzfBlockGuesser; // class finding the next BGZF block
	public SortOrder sortOrder;
	
	public CustomSamReader(final String path) throws IOException 
	{
		this.path = path;
		this.seekableStream = SeekableStreamFactory.getInstance().getStreamFor(Paths.get(path).toString());
		this.bgzfBlockGuesser= new BgzfBlockGuesser(this.seekableStream, this.path.toString());
		this.samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(this.seekableStream)); 
		this.samFileHeader = this.samReader.getFileHeader();
		this.sortOrder = this.samFileHeader.getSortOrder();
	}
	
	public SAMFileHeader getFileHeader() 
	{
		return samFileHeader;
	}

	/**
	 *  create a SAMRecord Iterator at the given byte offset
	 */
	public CloseableIterator<SAMRecord> queryAtOffset(final long start_offset) throws IOException 
	{
		if(start_offset >= this.seekableStream.length()) new ErrorMessage("Cannot create a reading thread with this offset (go beyond last file byte).");
		
		BgzfBlockGuesser.BgzfBlock bgzfBlock = this.bgzfBlockGuesser.guessNextBGZFPos(start_offset, seekableStream.length());
		if(bgzfBlock==null) new ErrorMessage("No bgz block found at "+start_offset);
		
		this.seekableStream.seek(bgzfBlock.pos);

		BlockCompressedInputStream bcis = new BlockCompressedInputStream(this.seekableStream);
		PushbackInputStream bpi = new PushbackInputStream(bcis, BamRecordGuesser.BUFFER_SIZE);
		BamRecordGuesser bamRecordGuesser = new BamRecordGuesser(this.samFileHeader);
		if(bamRecordGuesser.find(bpi)) return new MySamRecordIterator(this.samFileHeader, bpi, this.path.toString());
		return null;
	}
	
	@Override
	public void close() 
	{
		CloserUtil.close(this.seekableStream);
		this.bgzfBlockGuesser.close();
		try
		{
			this.samReader.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
	}
}