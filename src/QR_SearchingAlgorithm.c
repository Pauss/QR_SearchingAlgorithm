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
static void read_inputs(CMD_INPUTS* cmd_data, char** argv);
/*Description*/
/*=========================================*/
int main(int argc, char **argv)
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

 	CMD_INPUTS cmd_data;

 	if (argc != NR_PARAMETERS) {
 		check_file = file_not_found;

	} else {
		read_inputs(&cmd_data, argv);
		check_file = fileIsValid(cmd_data.my_file);

	}

	if( file_no_error == check_file)
	{

		/*File is valid*/

		/* choose what strategy to use
		 * 1. Naive search with QR decomposition applied at each step.
		 * 2. Efficient search with QR decomposition applied only on first step. Save some time tho..
		 * 3. Even more nice strategy using GA (Genetic Algorithm)*/

		T_EFFICIENT_METHOD columns_method = columns_removal;//columns_transitions; columns_removal;
		intercept = installed; //installed //not_installed

		file_dim = get_file_dimensions();

		if (NULL != file_dim) {
			set_y_vector(file_dim);
			set_A_matrix(file_dim);
			set_model_elements();

			if(USE_GRAPHICS == 1)
			{
				if (fopen(OUTPUT_FILE, "r"))
				{
					fclose(fopen(OUTPUT_FILE, "w"));
				}
			}

			switch(cmd_data.strategy)
			{

			case naive_search:{

				printf("\nPerforming Naive Search\n");

				naive_alg();

				break;
			}
			case efficient_search:{

				printf("\nPerforming Efficient Search\n");

				efficient_alg(columns_method);

				break;
			}
			case GA_search:{
				printf("\nPerforming GA Search\n");

				GA_naive_alg(cmd_data.method, cmd_data.op1, cmd_data.op2);

				break;
			}
			case GA_search_BB:{
				printf("\nPerforming GA Building Blocks Search\n");

				GA_BB_alg(cmd_data.op1, cmd_data.op2);

				break;
			}
			case GA_SA:{
				printf("\nPerforming Simulated Annealing Search\n");

				GA_simulated_annealing(cmd_data.op1, cmd_data.op2);

				break;
			}
			case GA_HC:{
				printf("\nPerforming Hill Climbing Search\n");

				GA_hill_climbing(cmd_data.op1, cmd_data.op2);

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
			printf("\nERROR: File not found.\n");
			printf(cmd_data.my_file);
			break;
		}
		case file_dim_invalid:
		{
			printf("\nERROR: File dimension is not valid.\n");
			break;
		}
		case file_data_invalid:
		{
			printf("\nERROR: File data is not valid.\n");
			break;
		}
		default:
		{
			printf("\nERROR wasn't detected.\n");
			break;
		}
	}

	}

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\nEXECUTION TIME %f\n", time_spent);

	//fclose(stdout);
	//fclose(stderr);
	return 0;

}

static void read_inputs(CMD_INPUTS* cmd_data, char** argv) {
	uint8 i = 0;

	cmd_data->my_file = (char*) malloc(PATH_SIZE * sizeof(char));

	//lower case
	while (argv[1][i]) {
		putchar(tolower(argv[1][i]));
		i++;
	}

	//setting strategy from command line
	if (0 == strcmp(argv[1], "ga")) {
		cmd_data->strategy = GA_search;
	} else if (0 == strcmp(argv[1], "ga_bb")) {
		cmd_data->strategy = GA_search_BB;
	} else if (0 == strcmp(argv[1], "ns")) {
		cmd_data->strategy = naive_search;
	} else if (0 == strcmp(argv[1], "ga_sa")) {
		cmd_data->strategy = GA_SA;
	} else if (0 == strcmp(argv[1], "ga_hc")) {
		cmd_data->strategy = GA_HC;
	} else {
		printf("\nInvalid strategy!\n");
		cmd_data->strategy = invalid_strategy;
	}

	//setting operator1 from command line
	if (0 == strcmp(argv[2], "1")) {
		cmd_data->op1 = flip;
	} else if (0 == strcmp(argv[2], "2")) {
		cmd_data->op1 = interchanging;
	} else if (0 == strcmp(argv[2], "3")) {
		cmd_data->op1 = interchanging_abs;
	} else if (0 == strcmp(argv[2], "4")) {
		cmd_data->op1 = reversing;
	} else {
		printf("\nInvalid crossover operator!\n");
		cmd_data->strategy = invalid_strategy;
	}

	//setting operator2 from command line
	if (0 == strcmp(argv[3], "1")) {
		cmd_data->op2 = _1point;
	} else if (0 == strcmp(argv[3], "2")) {
		cmd_data->op2 = uniform;
	} else if (0 == strcmp(argv[3], "3")) {
		cmd_data->op2 = RRC;
	} else if (0 == strcmp(argv[3], "4")) {
		cmd_data->op2 = _1point_simple;
	} else {
		printf("\nInvalid mutation operator2!\n");
		cmd_data->strategy = invalid_strategy;
	}

	//setting selection method from command line, only for GA strategy
	if (0 == strcmp(argv[4], "1")) {
		cmd_data->method = roulette_wheel;
	} else if (0 == strcmp(argv[4], "2")) {
		cmd_data->method = tournament;
	} else {
		printf("\nInvalid selection method!\n");
		cmd_data->strategy = invalid_strategy;
	}

	strcpy(cmd_data->my_file, GENERATED_DATA_PATH);

	cmd_data->my_file = strcat(cmd_data->my_file, argv[5]);

}
