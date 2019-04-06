/*
 ============================================================================
 Name        : QR_SearchingAlgorithm.c
 Author      : Paula Tanasa
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */
/*=========================================*/
/*includes*/
#include "QR_SearchingAlgorithm.h"
/*=========================================*/
/*private data*/
/*=========================================*/
/*global data*/
T_Intercept intercept = installed;
/*=========================================*/
/*private functions*/
/*=========================================*/
/*=========================================*/
/*Description*/
/*=========================================*/
int main()
{

	clock_t begin = clock();
	srand(time(NULL)*SEED_RAND);

	T_FILE_DIM* file_dim;
 	T_FILE_ERRORS check_file = file_no_error;
	/*
	 * Get the Data for applying search algorithm.
	 * Extract a matrix m*n from file.
	 * Applying algorithm on a specific type system equation (over_determined system)*/

	/* check validity of a file
	 * if valid extract matrix from file*/
	check_file = fileIsValid("house.txt"); //house //data_invalid

	if( file_no_error == check_file)
	{

		/*File is valid*/

		/* choose what strategy to use
		 * 1. Naive search with QR decomposition applied at each step.
		 * 2. Efficient search with QR decomposition applied only on first step. Save some time tho..
		 * 3. Even more nice strategy using GA (Genetic Algorithm)*/

		T_SEARCH_STRATEGIES strategy = GA_SA;//GA_search; naive_search //GA_SA // efficient_search

		//efficient method relevant
		T_EFFICIENT_METHOD columns_method = columns_transitions;//columns_transitions; columns_removal;

		//GA relevant
		T_SELECTION_METHOD selection_method = tournament;  //tournament; roulette_wheel
		T_OPERATOR_METHOD operator1 = flip; //flip //interchanging; interchanging_abs; reversing
		T_OPERATOR_METHOD operator2 = RRC; //_1point; uniform; RRC; _1point_simple; no_operator
		intercept = installed; //installed //not_installed

		file_dim = get_file_dimensions();

		if (NULL != file_dim) {
			set_y_vector(file_dim);
			set_A_matrix(file_dim);

			switch(strategy)
			{
			case naive_search:{

				printf("Performing Naive Search\n");

				naive_alg();

				break;
			}
			case efficient_search:{

				printf("Performing Efficient Search\n");

				efficient_alg(columns_method);

				break;
			}
			case GA_search:{
				printf("Performing GA Search\n");

				GA_naive_alg(selection_method, operator1, operator2);

				break;
			}
			case GA_SA:{
				printf("Performing GA Simulated Annealing Search\n");

				GA_simulated_annealing(operator1, operator2);

				break;
			}
			default: {
				/*TODO*/
				break;
			}

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

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("EXECUTION TIME %f", time_spent);

	return 0;

}
