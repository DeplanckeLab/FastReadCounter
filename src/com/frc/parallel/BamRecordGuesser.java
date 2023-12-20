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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PushbackInputStream;
import java.util.Arrays;

import com.errors.ErrorMessage;
import com.tools.Utils;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.StringUtil;

/**
 * This class and the rest of this package was taken from the great implementation of jvarkit from Pierre Lindenbaum: https://github.com/lindenb
 *
 */
public class BamRecordGuesser 
{
	public static final int BUFFER_SIZE=100_000;
	private static final int FIXED_SAM_BLOCK_SIZE = 36;
	
	private final SAMFileHeader header;
	private final SAMSequenceDictionary dict;
	private final BAMRecordCodec bamRecordCodec;
	private final int max_flag;
	
	private final ByteArrayOutputStream consumed = new ByteArrayOutputStream(BUFFER_SIZE);
	private final BinaryCodec binaryCodec=new BinaryCodec();
	
	public BamRecordGuesser(final SAMFileHeader header) 
	{
		this.header = header;
		this.dict = Utils.extractRequired(this.header);
		this.bamRecordCodec  = new BAMRecordCodec(this.header);
		this.max_flag = Arrays.stream(SAMFlag.values()).mapToInt(F->F.intValue()).sum();
	}
	
	public boolean find(final PushbackInputStream pbis) throws IOException 
	{
		for(;;) 
		{
			if(canDecode(pbis)) return true;
			final int c=pbis.read(); //skip one byte
			if(c==-1) return false;
		}
	}

	@SuppressWarnings("unused")
	private boolean canDecode(final PushbackInputStream pbis) throws IOException
	{
		this.consumed.reset();
		final LimitInputStream limitIs = new LimitInputStream(pbis, BUFFER_SIZE);
		final TeeInputStream tee=new TeeInputStream(limitIs, this.consumed, false);
		this.binaryCodec.setInputStream(tee);
		int n=0;
		try 
		{
			final int  recordLength = this.binaryCodec.readInt();
			n+=4;
			if(recordLength<FIXED_SAM_BLOCK_SIZE) return false;
			if(recordLength>BUFFER_SIZE) return false;
			
	        final int referenceID =  this.binaryCodec.readInt();
	        n+=4;
	        if(referenceID!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) 
	        {
	        	if( referenceID<0 || referenceID>= this.dict.size()) return false; 
	        }
	        
	        final int coordinate = this.binaryCodec.readInt() + 1;
	        n+=4;
	        if(referenceID!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) 
	        {
		        if(coordinate<1 || coordinate> this.dict.getSequence(referenceID).getSequenceLength()) return false; 
		    }
	        
	        final short readNameLength = this.binaryCodec.readUByte();
	        n++;
	        if(readNameLength <= 0) return false;
	        
	        final short mappingQuality = this.binaryCodec.readUByte();
	        n++;
	        if(mappingQuality> SAMRecord.UNKNOWN_MAPPING_QUALITY) return false;
	
	        
	        final int bin = this.binaryCodec.readUShort();
	        n+=2;
	        final int cigarLen = this.binaryCodec.readUShort();
	        n+=2;
	        
	        final int flags = this.binaryCodec.readUShort();
	        n+=2;
	        if(flags > this.max_flag) return false;
	       
	        
	        final int readLen = this.binaryCodec.readInt();
	        if(readLen<0) return false;
	        n+=4;
	       
	        final int mateReferenceID = this.binaryCodec.readInt();
	        n+=4;
	        final int mateCoordinate = this.binaryCodec.readInt() + 1;
	        n+=4;
	        if(SAMFlag.READ_PAIRED.isSet(flags)) 
	        {
	        	if(!SAMFlag.MATE_UNMAPPED.isSet(flags))
	        	{
	        		if(mateReferenceID<0 || mateReferenceID>= this.dict.size()) return false; 
	        		if(mateCoordinate<1 || mateCoordinate> this.dict.getSequence(mateReferenceID).getSequenceLength()) return false; 
	        	}
	        }
	        
	        final int insertSize = binaryCodec.readInt();   
	        n+=4;
	        if(referenceID>=0) 
	        {
	        	int contigLen = this.dict.getSequence(referenceID).getSequenceLength();
	        	if(Math.abs(insertSize) > contigLen) return false;
	        }
	        
	        if(n!=FIXED_SAM_BLOCK_SIZE) throw new IllegalStateException(" "+n+" vs "+FIXED_SAM_BLOCK_SIZE);
	        
	        final int restLen = recordLength - FIXED_SAM_BLOCK_SIZE;
	        byte mRestOfBinaryData[]=new byte[restLen];
	        Utils.readFully(tee, mRestOfBinaryData);
	        final String s=StringUtil.bytesToString(mRestOfBinaryData, 0, readNameLength - 1);
			for(int i=0;i< s.length();i++) 
			{
				int ascii = s.charAt(i);
				if(ascii < 32 || ascii > 127) return false;
			}
		}
		catch(final Throwable err)
		{
			new ErrorMessage(err.getMessage());
		}
		finally
		{
			final byte array[]= this.consumed.toByteArray();
			pbis.unread(array);
			this.consumed.reset();
		}
		this.consumed.reset();
		try 
		{
			this.bamRecordCodec.setInputStream(tee);
			final SAMRecord rec = this.bamRecordCodec.decode();
			return rec!=null;
		}
		catch(final Throwable err) 
		{
			new ErrorMessage(err.getMessage());
		}
		finally
		{
			final byte array[]= consumed.toByteArray();
			pbis.unread(array);
			this.consumed.reset();
		}
		return false;
	}
}



