Variables stored in separate files (CEOP formatted)

Filename
	
	Data_seperate_files_startdate(YYYYMMDD)_enddate(YYYYMMDD).zip

	e.g., Data_seperate_files_20050316_20050601.zip

	
Folder structure

	Networkname
		Stationname

		
Dataset Filename

	CSE_Network_Station_Variablename_depthfrom_depthto_startdate_enddate.ext

	CSE	- Continental Scale Experiment (CSE) acronym, if not applicable use Networkname
	Network	- Network abbreviation (e.g., OZNET)
	Station	- Station name (e.g., Widgiewa)
	Variablename - Name of the variable in the file (e.g., Soil-Moisture)
	depthfrom - Depth in the ground in which the variable was observed (upper boundary)
	depthto	- Depth in the ground in which the variable was observed (lower boundary)
	startdate -	Date of the first dataset in the file (format YYYYMMDD)
	enddate	- Date of the last dataset in the file (format YYYYMMDD)
	ext	- Extension .stm (Soil Temperature and Soil Moisture Data Set see CEOP standard)

	e.g., OZNET_OZNET_Widgiewa_Soil-Temperature_0.150000_0.150000_20010103_20090812.stm

	
File Content Sample

	2003/12/11 00:30 2003/12/11 00:40 OZNET      OZNET           Widgiewa         -35.09000   146.30600  121.00    0.15    0.15    28.30  U M

	UTC nominal date/time - yyyy/mm/dd HH:MM, where MM is 00 or 30, only
	UTC actual date/time - yyyy/mm/dd HH:MM
	CSE Identifier - Continental Scale Experiment (CSE) acronym, if not applicable use Networkname
	Network	- Network abbreviation (e.g., OZNET)
	Station	- Station name (e.g., Widgiewa)
	Latitude - Decimal degrees. South is negative. 
	Longitude -	Decimal degrees. West is negative.
	Elevation -	Meters above sea level
	Depth from - Depth in the ground in which the variable was observed (upper boundary)
	Depth to - Depth in the ground in which the variable was observed (lower boundary)
	Variable value
	ISMN Quality Flag
	Data Provider Quality Flag, if existing
	

For Definition of the CEOP Data Format see http://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/ceop_soils_format.html


Network Information

	COSMOS
		Abstract: A new project to measure soil moisture using a cosmic-ray technique. Currently, there are 67 stations deployed in 7 countries, which 59 are in USA, 1 in Germany, 1 in Switzerland, 1 in France, 1 in Brasil, 2 in Kenya, 1 in United Kingdom and 1 in Mexico.
		Continent: North_America
		Country: USA
		Stations: 109
		Status: running
		Data Range: from 2008-04-28
		Type: project
		Url: http://cosmos.hwr.arizona.edu/
		Reference: Zreda M., Desilets D., Ferré Ty P.A., Scott R.L., “Measuring soil moisture content non-invasively at intermediate spatial scale using cosmic-ray neutrons”, Geophysical Research Letters 35(21), 2008. ,   Zreda, M., W.J. Shuttleworth, X. Zeng, C. Zweck, D. Desilets, T. Franz, and R. Rosolem, 2012. COSMOS: the COsmic-ray Soil Moisture Observing System. Hydrology and Earth System Sciences 16, 4079-4099, doi: 10.5194/hess-16-4079-2012. 
		Variables: soil moisture, 
		Soil Moisture Depths: 0.00 - 0.04 m , 0.00 - 0.05 m , 0.00 - 0.06 m , 0.00 - 0.07 m , 0.00 - 0.08 m , 0.00 - 0.09 m , 0.00 - 0.10 m , 0.00 - 0.11 m , 0.00 - 0.12 m , 0.00 - 0.13 m , 0.00 - 0.14 m , 0.00 - 0.15 m , 0.00 - 0.16 m , 0.00 - 0.17 m , 0.00 - 0.18 m , 0.00 - 0.19 m , 0.00 - 0.20 m , 0.00 - 0.21 m , 0.00 - 0.22 m , 0.00 - 0.23 m , 0.00 - 0.24 m , 0.00 - 0.25 m , 0.00 - 0.26 m , 0.00 - 0.27 m , 0.00 - 0.28 m , 0.00 - 0.29 m , 0.00 - 0.30 m , 0.00 - 0.31 m , 0.00 - 0.32 m , 0.00 - 0.33 m , 0.00 - 0.34 m , 0.00 - 0.36 m , 0.00 - 0.37 m , 0.00 - 0.38 m , 0.00 - 0.39 m , 0.00 - 0.40 m , 0.00 - 0.41 m , 0.00 - 0.44 m , 0.00 - 0.55 m , 0.00 - 0.69 m , 0.00 - 0.90 m , 
		Soil Moisture Sensors: Cosmic-ray Probe, 

