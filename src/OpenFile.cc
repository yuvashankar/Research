 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

//This function opens a file, with all of the error checking

#include <stdio.h>
#include "edflib.h"
#include "processEEG.h"
#include <wavelet.h>

/**
	\fn int OpenFile(const char* fileName, struct edf_hdr_struct *header)
	
	\brief Openes a .BDF file and allocates it to an edf_hdr_struct.

	\param fileName The name and location of the file to be opened
	\param header The pointer to the edf header structure

	\return 0 if file is opened successfully
	\return -1 if there is an error

*/
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

	return(-1);
	}
	return(0);
}