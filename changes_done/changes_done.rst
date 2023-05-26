General:
---------
There was a Problem with reading data from the RISMA network. 
This data only contained soil moisture as variable and was of format  "header & files".
Further, this data was dowloaded from the ISMN webpage May 22, 23 at 09:00.


Problem:
----------
To my understanding, all sensor files of the "header & files" format should contain four columns: date/time, data, new flag, old flag.
For some RISMA sensor files the last column, old flag, contained strings such as:
    - "Frozen soil,OK,OK,OK,OK,OK" (see :example: `./RISMA_RISMA_MB11_sm_1.000000_1.000000_Hydraprobe-II-Sdi-12-A_19500101_20230521.stm`)
    - "DLT >= 1.5,OK,OK,OK,OK,OK" (see :example: `RISMA_RISMA_MB11_sm_0.000000_0.050000_Hydraprobe-II-Sdi-12-C_19500101_20230521.stm`)
    - "Out of WFV average range,NA,NA,OK,OK,OK" (see :example: `RISMA_RISMA_CEF_sm_0.000000_0.050000_Hydraprobe-II-Sdi-12-B_19500101_20230521.stm`) 

The method ``__read_csv()`` of the ``DataFile`` class could not parse these sensor files and threw following error: **pandas.errors.ParserError: Error tokenizing data. C error: Expected 5 fields in line 23, saw 9**.
This phenomenon only occured with some RISMA sensor files, other networks were fine.


Fix:
------
- Extended the ``__read_csv()`` method to take kezword arguments ``**kwargs`` as input
- Replaced lambda function ``readf()`` calling the ``pandas.read_csv()`` with a method of same name, that takes input arguments
- Defined additional default arguments: ``delim_whitespace=None, sep=None, low_memory=None``
- With a try and except block caught ParserError, and then call an adapted ``pd.read_csv()``, to fit the problems mentioned above
- Subsequently adapted the ``__read_format_header_values()`` method correctly to call the ``__read_csv()`` in its ``return`` statement
