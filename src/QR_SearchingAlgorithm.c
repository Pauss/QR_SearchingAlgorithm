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
#include "MatrixComputations.h"
#include "QR_SearchingAlgorithm.h"
/*=========================================*/
/*private data*/
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
/*=========================================*/

/*=========================================*/
/*Description*/
/*=========================================*/
int main()
{

	T_FILE_DIM file_dim;
	T_FILE_ERRORS check_file;
	gsl_matrix* input_matrix;
	/*
	 * Get the Data for applying search algorithm.
	 * Extract a matrix m*n from file.
	 * Applying algorithm on a specific type system equation (complete info)*/   /*todo*/

	/* check validity of a file
	 * if valid extract matrix from file*/
	check_file = fileIsValid("Data_Invalid.txt", &file_dim);

	if( file_no_error == check_file)
	{

		/*File is valid*/
		input_matrix = gsl_matrix_alloc(file_dim.lines, file_dim.columns);
		convert_to_gsl(&file_dim, input_matrix);
		if (TRUE == compute_matrix_inverse(input_matrix)) {
			print_matrix(input_matrix);

			convert_to_gsl(&file_dim, input_matrix);
			print_matrix(input_matrix);
		}

		/*choose what strategy to use
		 * 1. Naive search with QR decomposition applied at each step.
		 * 2. Efficient search with QR decomposition applied only on first step. Save some time tho..
		 * 3. Even more nice strategy using GA (Genetic Algorithm)*/

		T_SEARCH_STRATEGIES strategy = naive_search;

		switch(strategy)
		{
		case naive_search:{


			/*TODO*/
			break;
		}
		case efficient_search:{
			/*TODO*/
			break;
		}
		case GA_search:{
			/*TODO*/
			break;
		}
		default: {
			/*TODO*/
			break;
		}

		}



	}
	else
	{
		/*
		 * Message to the output
		 * Algorithm can't be verified on a invalid set of data*/
		switch(check_file)
		{
		case file_not_found:
		{
			printf("ERROR: File not found.\n");
			break;
		}
		case file_dim_invalid:
		{
			printf("ERROR: File dimension is not valid.\n");
			break;
		}
		case file_data_invalid:
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
