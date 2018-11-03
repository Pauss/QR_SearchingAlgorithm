/*
 ============================================================================
 Name        : Licenta_FBS_QR.c
 Author      : Paula Tanasa
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */
/*=========================================*/
/*includes*/
#include "Types.h"
#include "FileOperations.h"
/*=========================================*/
int main()
{

	FILE_DIM file_dim;
	FILE_ERRORS check_file;
	/*
	 * STEP1
	 * Get the Data for applying search algorithm.
	 * Extract a matrix m*n from file.
	 * Applying algorithm on a specific type system equation (complete info)*/   /*todo*/

	/* check validity of a file
	 * if valid extract matrix from file*/
	check_file = fileIsValid("Data_Invalid.txt", &file_dim);
	if( FILE_NO_ERROR == check_file)
	{

		printf("File is valid!\n");
	}
	else
	{
		/*
		 * Message to the output
		 * Algorithm can't be verified on a invalid set of data*/
		switch(check_file)
		{
		case FILE_NOT_FOUND:
		{
			printf("ERROR: File not found.\n");
			break;
		}
		case FILE_DIM_INVALID:
		{
			printf("ERROR: File dimension is not valid.\n");
			break;
		}
		case FILE_DATA_INVALID:
		{
			printf("ERROR: File data is not valid.\n");
			break;
		}
		default:
		{
			printf("ERROR wasn't detected.\n");
			break;
		}
	}

	}

	return 0;

}
