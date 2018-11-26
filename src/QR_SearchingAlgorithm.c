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
#include "MatrixComputations.h"
#include "QR_SearchingAlgorithm.h"
#include "AlgorithmComputations.h"
#include "FileOperations.h"
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

	T_FILE_DIM* file_dim;
	T_FILE_ERRORS check_file = file_no_error;
	/*
	 * Get the Data for applying search algorithm.
	 * Extract a matrix m*n from file.
	 * Applying algorithm on a specific type system equation (complete info)*/   /*TODO*/

	/* check validity of a file
	 * if valid extract matrix from file*/
	check_file = fileIsValid("Data_Invalid.txt");

	if( file_no_error == check_file)
	{

		/*File is valid*/


		/* choose what strategy to use
		 * 1. Naive search with QR decomposition applied at each step.
		 * 2. Efficient search with QR decomposition applied only on first step. Save some time tho..
		 * 3. Even more nice strategy using GA (Genetic Algorithm)*/

		T_SEARCH_STRATEGIES strategy = naive_search;

		switch(strategy)
		{
		case naive_search:{


			/*TODO*/
			file_dim = get_file_dimensions();

			if (NULL != file_dim) {
				set_y_vector(file_dim);
				set_A_matrix(file_dim);
				compute_transitions_QR();
			}
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

		/* Message to the output
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
