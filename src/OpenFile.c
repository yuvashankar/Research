//This function opens a file, with all of the error checking

#include <stdio.h>
#include "edflib.h"


int OpenFile(const char* fileName, struct edf_hdr_struct *header)
{
	if(edfopen_file_readonly(fileName, header, EDFLIB_READ_ALL_ANNOTATIONS))
	{
		switch(header->filetype)
		{
			case EDFLIB_MALLOC_ERROR                : printf("\nmalloc error\n\n");
			                                        break;
			case EDFLIB_NO_SUCH_FILE_OR_DIRECTORY   : printf("\ncan not open file, no such file or directory\n\n");
			                                        break;
			case EDFLIB_FILE_CONTAINS_FORMAT_ERRORS : printf("\nthe file is not EDF(+) or BDF(+) compliant\n"
			                                               "(it contains format errors)\n\n");
			                                        break;
			case EDFLIB_MAXFILES_REACHED            : printf("\nto many files opened\n\n");
			                                        break;
			case EDFLIB_FILE_READ_ERROR             : printf("\na read error occurred\n\n");
			                                        break;
			case EDFLIB_FILE_ALREADY_OPENED         : printf("\nfile has already been opened\n\n");
			                                        break;
			default                                 : printf("\nunknown error\n\n");
			                                        break;
		}

	return(1);
	}
	return(0);
}