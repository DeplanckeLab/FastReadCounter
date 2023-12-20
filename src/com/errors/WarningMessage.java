package com.errors;

public class WarningMessage 
{
	public String displayed_warning = "Warning";
	
	public WarningMessage(String warningMessage) 
	{
		if(warningMessage != null) this.displayed_warning = warningMessage;
		System.err.println("[Warning] " + this.displayed_warning);
	}
}
