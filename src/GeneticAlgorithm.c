/*
 * GeneticAlgorithm.c
 *
 *  Created on: Mar 3, 2019
 *      Author: Paus
 */
/*=========================================*/
/*include*/
#include "GeneticAlgorithm.h"
/*=========================================*/
/*private data*/
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
/*=========================================*/


/*=========================================*/
/*This function compute a possible subset of ModelMatrix*/
/*=========================================*/
boolean get_random_model(gsl_vector* my_random_model, uint8 n, uint8 size_A)
{
	uint8 i = n;

	gsl_vector_set_zero(my_random_model);

	if (n >= my_random_model->size) {
		return FALSE;
	} else {
		while (i) {

			uint8 r = rand() % (size_A);

			//for getting n random columns
			while (gsl_vector_get(my_random_model, r) == 1) {
				r = rand() % (size_A);
			}

			gsl_vector_set(my_random_model, r, 1);
			i--;
		}
	}

	return TRUE;

}

void convert_submodels(gsl_vector* bit_model, gsl_vector* model)
{
	uint8 found = 0;

	model->size = bit_model->size;
	gsl_vector_set_zero(model);

	for (uint8 j = 0; j < bit_model->size; j++) {
		//if column was selected for sub-model set it with column number
		if (gsl_vector_get(bit_model, j) == 1) {

			gsl_vector_set(model, found, j);
			found++;
		}
	}

	model->size = found;

}

/*=========================================*/
/*This function applies an criterion over sub-models
 * AIC = n + n log 2π + n log(RSS/n) + 2(p + 1)*/
/*=========================================*/
boolean criterion(double RSS, uint8 n, uint8 k, double* result) {

	//result = (double) ((-2) * log(RSS)) + 2 * k;

	if(intercept)
	{
		k = k+1;
	}


	static double best_s = 10000;

	if(n)
	{
		*result = (double) (n + n * log(2 * M_PI) + n * log(RSS / n) + 2 * (k + 1));
	}

	if (*result <= best_s) {
		best_s = *result;

		return TRUE;

	}

	return FALSE;
}

/*=========================================*/
/*This function build and assign to an individual all corresponding values*/
/*=========================================*/
void build_individual (T_INDIVIDUAL* GA_individual, gsl_vector* rand_model, uint8 columns)
{

	GA_individual->columns = gsl_vector_alloc(rand_model->size);
	GA_individual->bit_columns = gsl_vector_alloc(rand_model->size);
	GA_individual->submodel = gsl_matrix_alloc(columns, columns);

	gsl_vector_memcpy(GA_individual->bit_columns, rand_model);

	convert_submodels(rand_model, GA_individual->columns);

	GA_individual->submodel = sub_model_matrix(GA_individual->columns);
	GA_individual->RSS = RSS_compute(GA_individual->submodel);
	GA_individual->fitness_value = 0;
	GA_individual->selection_probability = 0;
	GA_individual->selected = FALSE;
}

/*=========================================*/
/*This function implements generation of population*/
/*=========================================*/
void generate_population(T_INDIVIDUAL* GA_population, uint8 n, uint8 k)
{
	gsl_vector* random_model = gsl_vector_alloc(n);

	//generate population
	for (uint8 i = 0; i < NUMBER_OF_CHROMOSOMES; i++) {

		if (FALSE != get_random_model(random_model, k, n)) {

			build_individual(&GA_population[i], random_model, n);

		}

		//if size changes
		random_model->size = n;
	}

}

/*=========================================*/
/*This function implements Roultte Wheel selection*/
/*=========================================*/
void roulette_wheel(T_INDIVIDUAL* population, uint8 size) {
	double sum_fitness = 0;
	double sum_probabilities = 0;

	for (uint8 i = 0; i < size; i++) {
		sum_fitness += population[i].fitness_value;
	}

	for (uint8 i = 0; i < size; i++) {
		//it's used cumulative probability
		population[i].selection_probability = sum_probabilities + (double) (population[i].fitness_value / sum_fitness);
		sum_probabilities += (double) (population[i].fitness_value / sum_fitness);
	}

	for (uint8 i = 0; i < size - 1; i++) {

		double r = (double) (rand() / (double) (RAND_MAX));

		//r bigger than actual but less then next individual
		if (population[i].selection_probability < r
				&& r <= population[i + 1].selection_probability) {
			population[i].selected = TRUE;
		} else
			population[i].selected = FALSE;
	}

}

/*=========================================*/
/*This function implements a shuffle array elements method 'Fisher–Yates shuffle*/
/*=========================================*/
void shuffle_array(T_INDIVIDUAL* population, uint8 size) {
	T_INDIVIDUAL temp_population[size];

	for (uint8 i = size - 1; i > 0; i--) {
		uint8 j = rand() % i;

		copy_individual(temp_population, population, j, i);
		copy_individual(temp_population, population, i, j);

		copy_individual(population, temp_population, i, i);
		copy_individual(population, temp_population, j, j);

	}

	copy_population(population, temp_population, size);

}

uint8 get_index_of_BEST(T_INDIVIDUAL* population, uint8 size){

	uint8 best_index = 0;

	double min = population[best_index].fitness_value;

	for(uint8 index = 1; index< size; index++){
		if(population[index].fitness_value < min)
		{
			min = population[index].fitness_value;
			best_index = index;
		}
	}

	return best_index;
}


/*=========================================*/
/*This function implements Tournament selection*/
/*=========================================*/
void tournament_selection(T_INDIVIDUAL* population, uint8 size, uint8 k) {

	/*Steps:
	 * 1. Pick k random individuals from population
	 * 2. Select the best individual and add it to new population
	 * 3. repeat step 1 and 2 untile seletd deisred number of individuals
	*/

	gsl_vector* temp_mutation;
	uint8 model_size_n = population[0].columns->size;

	T_INDIVIDUAL temp_population[size];
	T_INDIVIDUAL temp_pool_population[2];
	T_INDIVIDUAL pool_population[size];

	for (uint8 i = 0; i < size; i += 2) {
		copy_population(temp_population, population, size);

		//get first element

		shuffle_array(temp_population, size);

		uint8 index1 = get_index_of_BEST(temp_population, k);

		//add it to temp pool
		copy_individual(temp_pool_population, temp_population, 0, index1);

		//get second element

		shuffle_array(temp_population, size);

		uint8 index2 = get_index_of_BEST(temp_population, k);

		//add it to pool
		copy_individual(temp_pool_population, temp_population, 1, index2);

		//apply muations
		new_population_computed(temp_pool_population, 2,
							model_size_n);

		//add new 2 individuals to generation
		copy_individual(pool_population, temp_pool_population, i, 0);
		copy_individual(pool_population, temp_pool_population, i+1, 1);

	}

	copy_population(population, pool_population);

}

/*=========================================*/
/*This function apply crossover on 2 individuals on a random position*/
/*=========================================*/
void crossover(gsl_vector* bit_individ_1, gsl_vector *bit_individ_2) {

	uint8 r = 0, i;
	while (r == 0 || r == bit_individ_1->size) {
		r = rand() % bit_individ_1->size;
	}

	for (i = r; i < bit_individ_1->size; i++) {
		double t = gsl_vector_get(bit_individ_1, i);
		gsl_vector_set(bit_individ_1, i, gsl_vector_get(bit_individ_2, i));
		gsl_vector_set(bit_individ_2, i, t);
	}
}

/*=========================================*/
/*This function counts number of bits set on 1*/
/*=========================================*/
int get_number_locus(gsl_vector *bit_individ) {
	int i, nr = 0;
	for (i = 0; i < bit_individ->size; i++) {
		if (gsl_vector_get(bit_individ, i) == 1)
			nr++;
	}

	return nr;
}

/*=========================================*/
/*This function apply mutation on individual on a n random positions, n also randomly choosed*/
/*=========================================*/
void flip_mutation(gsl_vector* bit_columns) {
	// get n random positions and change its value : if it's 0 become 1, and reverse;

	gsl_vector* aux = gsl_vector_alloc(bit_columns->size);

	uint8 n_random_positions = 0;

	while (n_random_positions == 0) {
		n_random_positions = rand() % (bit_columns->size);
	}

	while (n_random_positions) {
		gsl_vector_memcpy(aux, bit_columns);
		uint8 r = rand() % (aux->size);
		uint8 temp_el = gsl_vector_get(aux, r);
		temp_el ^= 1UL;
		gsl_vector_set(aux, r, temp_el);

		n_random_positions--;

	}

	gsl_vector_memcpy(bit_columns, aux);
	gsl_vector_free(aux);
}

/*=========================================*/
/*This function apply mutation on individuals on a 2 random positions, by switchig it's elements*/
/*=========================================*/
void interchanging_mutation (gsl_vector* bit_columns) {

	uint8 pos1 = 0;
	uint8 pos2 = 0;

	while(pos1 == pos2)
	{
		pos1 = rand() % (bit_columns->size);
		pos2 = rand() % (bit_columns->size);
	}

	gsl_vector_swap_elements(bit_columns, pos1, pos2);
}

/*=========================================*/
/*This function apply mutation on individuals that are after one random position, by reversing it's bit value*/
/*=========================================*/
void reversing_mutation (gsl_vector* bit_columns) {

	uint8 pos = 0;

	while (pos >= bit_columns->size - 1) {
		pos = rand() % (bit_columns->size);
	}

	for (uint8 i = pos + 1; i < bit_columns->size; i++) {
		uint8 temp_el = gsl_vector_get(bit_columns, i);
		temp_el ^= 1UL;
		gsl_vector_set(bit_columns, i, temp_el);
	}
}


/*=========================================*/
/*This function compute new population using selection method*/
/*=========================================*/
void new_population_computed(T_INDIVIDUAL* temp_population, uint8 new_size, uint8 n)
{
	gsl_vector* temp_mutation;

	//Apply mutation on selected individuals
	for (uint8 i = 0; i < new_size; i++) {
		temp_mutation = gsl_vector_alloc(temp_population[i].bit_columns->size);

		gsl_vector_memcpy(temp_mutation, temp_population[i].bit_columns);

		interchanging_mutation(temp_mutation);

		build_individual(&temp_population[i], temp_mutation, n);

	}
}

/*=========================================*/
/*This function copy values of a source population into destination population */
/*=========================================*/
void copy_population(T_INDIVIDUAL* dest, T_INDIVIDUAL* src, uint8 new_size) {
	for (uint8 i = 0; i < new_size; i++) {
		dest[i].RSS = src[i].RSS;
		dest[i].fitness_value = src[i].fitness_value;
		dest[i].selected = src[i].selected;
		dest[i].selection_probability = src[i].selected;

		dest[i].columns = gsl_vector_alloc(src[i].columns->size);
		dest[i].bit_columns = gsl_vector_alloc(src[i].bit_columns->size);
		dest[i].submodel = gsl_matrix_alloc(src[i].submodel->size1,
				src[i].submodel->size2);

		gsl_vector_memcpy(dest[i].columns, src[i].columns);
		gsl_vector_memcpy(dest[i].bit_columns, src[i].bit_columns);
		gsl_matrix_memcpy(dest[i].submodel, src[i].submodel);

	}
}

/*=========================================*/
/*This function copy values of a source population into destination population */
/*=========================================*/
void copy_individual(T_INDIVIDUAL* dest, T_INDIVIDUAL* src, uint8 index1,
		uint8 index2) {

	dest[index1].RSS = src[index2].RSS;
	dest[index1].fitness_value = src[index2].fitness_value;
	dest[index1].selected = src[index2].selected;
	dest[index1].selection_probability = src[index2].selected;

	dest[index1].columns = gsl_vector_alloc(src[index2].columns->size);
	dest[index1].bit_columns = gsl_vector_alloc(src[index2].bit_columns->size);
	dest[index1].submodel = gsl_matrix_alloc(src[index2].submodel->size1,
			src[index2].submodel->size2);

	gsl_vector_memcpy(dest[index1].columns, src[index2].columns);
	gsl_vector_memcpy(dest[index1].bit_columns, src[index2].bit_columns);
	gsl_matrix_memcpy(dest[index1].submodel, src[index2].submodel);

}

/*=========================================*/
/*This function compute new population using selection method*/
/*=========================================*/
void population_selected(T_INDIVIDUAL* population,
		T_INDIVIDUAL* temp_population, uint8 size, uint8* new_size) {

	uint8 temp_index = 0;

	for (uint8 i = 0; i < size; i++) {

		if (FALSE != population[i].selected) {

			copy_individual(temp_population, population, temp_index, i);

			temp_index++;
		}
	}

	*new_size = temp_index;

}

/*=========================================*/
/*This function prints to the console all data of an individual*/
/*=========================================*/
void print_individual(T_INDIVIDUAL individual)
{
	printf("\nColumns:");
	print_vector(individual.columns);
	printf("Fitness Value: %f\n", individual.fitness_value);
	printf("Selection Probability: %f\n", individual.selection_probability);
	printf("Is Selected: %d\n", individual.selected);
}

/*=========================================*/
/*This function calculate fitness function of each individual and return index of the best*/
/*=========================================*/
void fitness_func(T_INDIVIDUAL* population,uint8 model_size_n, uint8 model_size_k, double* result)
{

	for (uint8 i = 0; i < NUMBER_OF_CHROMOSOMES; i++) {

		criterion(population[i].RSS, model_size_n, model_size_k, result);

		//in each case update fitness value of each individual
		population[i].fitness_value = *result;

	}
}


/*The main components are the chromosome encoding
* the fitness function, selection, recombination and the evolution scheme.*/
//https://www.sciencedirect.com/science/article/pii/S0377042705000774
/*=========================================*/
/*This function implements a naive approach of genetic algorithm*/
/*=========================================*/
void naive_GA_alg(void)
{

	// Initialization of rand() function
	srand(time(NULL));

	gsl_matrix * main_model = get_A_matrix();

	Model_QR_components matrix_components;

	uint8 model_size_k = NUMBER_OF_GENES;
	uint8 model_size_n = main_model->size2;
	uint8 population_size = NUMBER_OF_CHROMOSOMES;
	uint8 new_size = 0;
	uint8 best_solution_index = 0;
	uint8 generation = 1;
	uint8 converge_value = 0;
	double temperature = TEMP;
	double result = 1;
	double MIN = (double) MIN_FITNESS;

	T_INDIVIDUAL GA_population[population_size];
	T_INDIVIDUAL new_GA_population[population_size];

	QR_decomposition(main_model, &matrix_components);

	//generate population
	generate_population(GA_population, model_size_n, model_size_k);

	//apply fitness to all individuals
	fitness_func(GA_population, model_size_n, model_size_k,
					&result);
	//get BEST
	best_solution_index = get_index_of_BEST(GA_population, population_size);

	if (GA_population[best_solution_index].fitness_value < MIN)

	{
		MIN = GA_population[best_solution_index].fitness_value;
	}

	//temperature = 0;

	while (temperature)

	{
		printf("\n==========#Generation %d#===========", generation);


		 /*select best individual for next generation (x %)
		 * 1. Roulette Wheel -> how % from SUM of all Results represents each Result)
		 * 2. Tournament Selection -> first selects two individuals with uniform probability -> chooses the one with the highest fitness.
		 * 3. Truncation Selection -> simply selects at random from the population having first eliminated K number of the least fit individuals
		 */

		tournament_selection(GA_population, population_size, TOURNAMENT_K);

		//apply fitness to all individuals
		fitness_func(GA_population, model_size_n, model_size_k,
						&result);
		//get BEST
		best_solution_index = get_index_of_BEST(GA_population, population_size);

		if (GA_population[best_solution_index].fitness_value < MIN)

		{
			MIN = GA_population[best_solution_index].fitness_value;
			temperature--;
		}
		else
		{
			converge_value ++;
		}

		/*roulette_wheel(GA_population, population_size);

		for (uint8 i = 0; i < population_size; i++) {
			printf("\nIndividual index: %d", i);
			print_individual(GA_population[i]);

		}

		population_selected(GA_population, new_GA_population, population_size,
				&new_size);

		copy_population(GA_population, new_GA_population, new_size);

		population_size = new_size;

		if (new_size) {

			printf("\n=============Selected==============\n");
			for (uint8 i = 0; i < population_size; i++) {
				print_individual(GA_population[i]);
			}

			//apply mutation and crossover on selected chromosomes => get next generation

			new_population_computed(GA_population, population_size,
					model_size_n);

			printf("\n==========Mutation Applied===========\n");
			for (uint8 i = 0; i < population_size; i++) {
				print_individual(GA_population[i]);
			}

		} else {
			printf("\nNo individual survived at generation %d!\n", generation);
			temperature = 0;
		}*/

		if (converge_value < CONVERGE) {

			generation++;
		} else {
			temperature = 0;
		}

	}

	printf("\nBEST\n");
	print_individual(GA_population[best_solution_index]);

}
