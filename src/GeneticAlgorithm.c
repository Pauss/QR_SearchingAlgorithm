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
/*This function implements generation of population*/
/*=========================================*/
void generate_population(T_INDIVIDUAL* GA_population, uint8 n, uint8 k)
{
	gsl_vector* random_model = gsl_vector_alloc(n);
	gsl_vector* bit_model = gsl_vector_alloc(n);

	//generate population
	for (uint8 i = 0; i < NUMBER_OF_CHROMOSOMES; i++) {

		if (FALSE != get_random_model(random_model, k, n)) {
			GA_population[i].columns = gsl_vector_alloc(random_model->size);
			GA_population[i].bit_columns = gsl_vector_alloc(random_model->size);
			GA_population[i].submodel = gsl_matrix_alloc(random_model->size,
					random_model->size);

			gsl_vector_memcpy(GA_population[i].bit_columns, random_model);

			convert_submodels(random_model, bit_model);

			GA_population[i].columns = gsl_vector_alloc(bit_model->size);

			gsl_vector_memcpy(GA_population[i].columns, bit_model);

			GA_population[i].submodel = sub_model_matrix(
					GA_population[i].columns);
			GA_population[i].RSS = RSS_compute(GA_population[i].submodel);
			GA_population[i].selection_probability = 0;
			GA_population[i].selected = FALSE;

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
		//double r = (double) rand() / (double) (RAND_MAX / 0.1);
		double r = (double) (rand() / (double)(RAND_MAX));

		printf("\nR = %lf\n",r);

		//r bigger than actual but less then next individual
		if (population[i].selection_probability < r
				&& r <= population[i + 1].selection_probability) {
			population[i].selected = TRUE;
		} else
			population[i].selected = FALSE;
	}

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

}

/*=========================================*/
/*This function apply crossover on 2 individuals on a random position*/
/*=========================================*/
void crossover(gsl_vector* bit_individ_1, gsl_vector *bit_individ_2) {

	uint8 r = 0,i;
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
	int i,nr=0;
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
void new_population_computed(T_INDIVIDUAL* temp_population, uint8 new_size)
{
	//Apply mutation on selected individuals
	for(uint8 i = 0; i< new_size; i++)
	{
		interchanging_mutation(temp_population[i].bit_columns);
	}
}

/*=========================================*/
/*This function compute new population using selection method*/
/*=========================================*/
void population_selected(T_INDIVIDUAL* population, T_INDIVIDUAL* temp_population, uint8 size, uint8* new_size)
{

	uint8 temp_index = 0;

	for (uint8 i = 0; i < size; i++) {

		if(FALSE != population[i].selected)
		{

			temp_population[temp_index].bit_columns =  gsl_vector_alloc(population[i].bit_columns->size);
			gsl_vector_memcpy(temp_population[temp_index].bit_columns, population[i].bit_columns);

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
	printf("\nChormosome Values\n");
	print_vector(individual.columns);
	printf("Fitness Value: %lf\n", individual.fitness_value);
	printf("Selection Probability: %lf\n", individual.selection_probability);
	printf("Is Selected: %d\n", individual.selected);
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
	double result = 1;
	uint8 best_solution_index = 0;

	T_INDIVIDUAL GA_population[population_size];
	T_INDIVIDUAL new_GA_population[population_size];

	QR_decomposition(main_model, &matrix_components);

	//generate population
	generate_population(GA_population, model_size_n, model_size_k);

	//apply fitness function on each individual and save best fitness
	for (uint8 i = 0; i < NUMBER_OF_CHROMOSOMES; i++) {
		if (criterion(GA_population[i].RSS, model_size_n, model_size_k,
				&result)) {
			best_solution_index = i;
		}

		//in each case update fitness value of each individual
		GA_population[i].fitness_value = result;

	}

	//print best individual
	print_individual(GA_population[best_solution_index]);

	/* select best individual for next generation (x %)
	 * 1. Roulette Wheel -> how % from SUM of all Results represents each Result)
	 * 2. Tournament Selection -> first selects two individuals with uniform probability -> chooses the one with the highest fitness.
	 * 3. Truncation Selection -> simply selects at random from the population having first eliminated K number of the least fit individuals
	 * 4.  */

	roulette_wheel(GA_population, population_size);

	for (uint8 i = 0; i < NUMBER_OF_CHROMOSOMES; i++) {
		printf("\nIndivid %d\n", i);
		print_individual(GA_population[i]);

	}

	/*Evolution
	 * 1. Replacement-with-elitism -> almost complete replacement except that the best one or two individuals from the source population are preserved in the successor population
	 *
	 */

	//apply mutation and crossover on selected chromosomes => get next generation
	population_selected(GA_population, new_GA_population, population_size,
			&new_size);


	for (uint8 i = 0; i < new_size; i++) {
		printf("\nIndivid %d\n", i);
		print_vector(new_GA_population[i].bit_columns);

	}

	if (new_size) {

		new_population_computed(new_GA_population, new_size);

		printf("\nAfter\n");
		for (uint8 i = 0; i < new_size; i++) {
			printf("\nIndivid %d\n", i);
			print_vector(new_GA_population[i].bit_columns);

		}
	}
	else
	{
		printf("\nNo individual survived!!\n");
	}

	//finish when best solution isn't updated for more than k iterations
	//OR
	// finish when temperature is 0


}