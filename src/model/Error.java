package model;

public class Error 
{
	public Error(String message) 
	{
		System.err.println(message);
		System.exit(-1);
	}
}
