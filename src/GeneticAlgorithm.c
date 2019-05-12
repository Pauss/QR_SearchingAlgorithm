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
static void generate_population(T_INDIVIDUAL2* GA_population, uint16 size, uint16 n, uint16 k);
static void generate_individual(T_INDIVIDUAL2* individual, uint16 n, uint16 k);
static DATA_ERRORS convert_submodels(uint8* bit_model, uint16* model, uint16* size);
static void build_individual (T_INDIVIDUAL2* GA_individual, uint8* rand_model, uint16 columns);
static void new_population_computed(T_INDIVIDUAL2* temp_population, uint16 new_size, uint16 n, T_OPERATOR_METHOD op);

static void set_BEST(T_INDIVIDUAL2* best, T_INDIVIDUAL2* best_found, double* MIN);

static void copy_individual(T_INDIVIDUAL2* dest, T_INDIVIDUAL2* src);
static void copy_individual_into_population(T_INDIVIDUAL2* dest, T_INDIVIDUAL2* src, uint16 index1, uint16 index2);
static void copy_population(T_INDIVIDUAL2* dest, T_INDIVIDUAL2* src, uint16 new_size);

static void probability_selection(T_INDIVIDUAL2* population, uint16 size);
static void population_selected(T_INDIVIDUAL2* population,
		T_INDIVIDUAL2* temp_population, uint16 size, uint16* new_size);
static void selection_roulette_wheel(T_INDIVIDUAL2* population, uint16* size, uint16 model_size_n);

static DATA_ERRORS shuffle_array(T_INDIVIDUAL2* population, uint16 size);
static void selection_tournament(T_INDIVIDUAL2* population, uint16 size, uint16 k);

static boolean selection_building_blocks(T_INDIVIDUAL2* population, T_INDIVIDUAL2* schema, uint16* size, uint16 size_k);
static boolean individual_match_schema(T_INDIVIDUAL2* individual,T_INDIVIDUAL2* schema);

static boolean neighbor_acceptance(T_INDIVIDUAL2* current, T_INDIVIDUAL2* neighbor, double temperature);

static void crossover_1point_simple(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size);
static void crossover_1point(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size);
static void crossover_uniform(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size);
static void crossover_RRC(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size);

static void mutation_flip(uint8* bit_columns, uint16 size);
static void mutation_interchanging_abs(uint8* bit_columns, uint16 size);
static void mutation_interchanging (uint8* bit_columns, uint16 size);
static void mutation_reversing (uint8* bit_columns, uint16 size);

static uint16 get_index_of_BEST(T_INDIVIDUAL2* population, uint16 size);
static uint16 get_number_locus(uint8 *bit_individ, uint16 size);

static void print_individual(T_INDIVIDUAL2* individual);
static void print_population(T_INDIVIDUAL2* population, uint16 size);

static uint16 get_nr_genes(uint16 columns);

/*=========================================*/
/*This function compute a possible subset of ModelMatrix*/
/*=========================================*/
void get_random_model(uint8* my_random_model, uint16 k, uint16 n) {

	uint16 i = k;

	while (i) {

		uint16 r = rand() % (n);

		//for getting n random columns
		while (my_random_model[r] == 1) {
			r = rand() % (n);
		}

		my_random_model[r] = 1;
		i--;
	}

}

static DATA_ERRORS convert_submodels(uint8* bit_model, uint16* model, uint16* size) {
	uint16 found = 0;
	DATA_ERRORS ret_val = data_no_error;

	if (*size == 0) {

		ret_val = data_nr_of_columns_0;
	}

	else {

		for (uint8 j = 0; j < *size; j++) {
			//if column was selected for sub-model set it with column number
			if (bit_model[j] == 1) {
				model[found] = j;
				found++;
			}
		}

		if (found) {

			*size = found;
			//model = realloc(model, sizeof(uint16) * found);
		}

		else {
			ret_val = data_nr_of_columns_0;
		}

	}

	return ret_val;

}

/*=========================================*/
/*This function applies an criterion over sub-models
 * AIC = n + n log 2π + n log(RSS/n) + 2(p + 1)*/
/*=========================================*/
boolean criterion(double RSS, uint16 n, uint16 k, double* result) {

	//criterion based on AIC
	if (intercept) {
		k = k + 1;
	}

	static double best_s = MAX_FITNESS;

	if (n > 0 && RSS > 0 && k > 0) {

		*result = (double) (n + n * log(2 * M_PI) + n * log(RSS / n)
				+ 2 * (k + 1)); //(K+1)
	} else {
		*result = MAX_FITNESS;
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
static void build_individual (T_INDIVIDUAL2* GA_individual, uint8* rand_model, uint16 columns)
{

	DATA_ERRORS temp_error = data_no_error;
	gsl_matrix* temp_submodel;

	copy_vector8(GA_individual->bit_columns, rand_model, columns);

	temp_error = convert_submodels(rand_model, GA_individual->columns, &columns);

	if (temp_error == data_no_error) {

		GA_individual->size = columns;
		temp_submodel = submodel_matrix(GA_individual->columns, GA_individual->size);
		GA_individual->RSS = RSS_compute(temp_submodel);
		GA_individual->fitness_value = MAX_FITNESS;
		GA_individual->selection_probability = 0;
		GA_individual->selected = FALSE;


		gsl_matrix_free(temp_submodel);
	} else {

		GA_individual->fitness_value = MAX_FITNESS;
		GA_individual->selection_probability = 0;
		GA_individual->selected = FALSE;
		GA_individual->size = 0;

	}

}

/*=========================================*/
/*This function implements generation of population*/
/*=========================================*/
static void generate_population(T_INDIVIDUAL2* GA_population, uint16 size, uint16 n, uint16 k)
{
	uint8* random_model = (uint8*) calloc(n, sizeof(uint8));

	//generate population
	for (uint16 i = 0; i < size; i++) {

		individual_init(&GA_population[i]);
		generate_individual(&GA_population[i], n, k);
	}

	free(random_model);

}

/*=========================================*/
/*This function implements generation of an individual*/
/*=========================================*/
static void generate_individual(T_INDIVIDUAL2* individual, uint16 n, uint16 k)
{
	uint8* random_model = (uint8*) calloc(n, sizeof(uint8));

	if (n) {
		get_random_model(random_model, k, n);

		build_individual(individual, random_model, n);

		//if size changes
		//random_model = realloc(random_model, sizeof(uint8) * n);
	}
	free(random_model);
}

/*=========================================*/
/*This function compute probability of each individual to be selected in next generation*/
/*=========================================*/
static void probability_selection(T_INDIVIDUAL2* population, uint16 size) {
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
/*This function identifies if an individual match with a Schema*/
/*=========================================*/
static boolean individual_match_schema(T_INDIVIDUAL2* individual,T_INDIVIDUAL2* schema) {

	boolean ret_val = TRUE;

	if ((schema->size == 0) || (schema->columns == NULL)) {
		ret_val = FALSE;
	}

	if ((individual->size == 0) || (individual->columns == NULL)) {
		ret_val = FALSE;
	}

	if (schema->size > individual->size) {
		ret_val = FALSE;
	}

	else {

		for (uint16 i = 0; i < schema->size_bit; i++) {
			if (schema->bit_columns[i] == 1)
				if (individual->bit_columns[i] == 0) {
					ret_val = FALSE;
				}
		}
	}

	return ret_val;

}


/*=========================================*/
/*This function implements a shuffle array elements method 'Fisher–Yates shuffle*/
/*=========================================*/
static DATA_ERRORS shuffle_array(T_INDIVIDUAL2* population, uint16 size) {

	DATA_ERRORS ret_val = data_no_error;
	T_INDIVIDUAL2 temp_individual;

	individual_init(&temp_individual);

	if (size) {

		for (uint16 i = size - 1; i > 0; i--) {
			uint16 j = rand() % i;

			copy_individual(&temp_individual, &population[i]);
			copy_individual_into_population(population, population, i, j);
			copy_individual_into_population(population, &temp_individual, j , 0);
		}

	}

	else {
		ret_val = data_nr_of_individuals_0;
		printf("Shuffle array: error: %d", ret_val);
	}

	individual_dealloc(&temp_individual);

	return ret_val;
}
/*=========================================*/
/*This function gets index of best individual, best fitness output*/
/*=========================================*/
static uint16 get_index_of_BEST(T_INDIVIDUAL2* population, uint16 size){

	uint16 best_index = 0;

	double min = population[best_index].fitness_value;

	for(uint16 index = 1; index< size; index++){
		if(population[index].fitness_value < min)
		{
			min = population[index].fitness_value;
			best_index = index;
		}
	}

	return best_index;
}

static void set_BEST(T_INDIVIDUAL2* best, T_INDIVIDUAL2* best_found, double* MIN)
{
	if (best_found->fitness_value < *MIN)

	{
		*MIN = best_found->fitness_value;
		copy_individual(best, best_found);
	}
}

/*=========================================*/
/*This function implements Tournament selection*/
/*=========================================*/
static void selection_tournament(T_INDIVIDUAL2* population, uint16 size, uint16 k) {

	/*Steps:
	 * 1. Pick k random individuals from population (shuffled population)
	 * 2. Select the best individual and add it to new population
	 * 3. repeat step 1 and 2 until selected desired number of individuals
	 */

	T_INDIVIDUAL2 temp_population[size];

	T_INDIVIDUAL2 pool_population[size];

	for(uint16 i = 0 ; i< size; i++)
	{
		individual_init(&temp_population[i]);
		individual_init(&pool_population[i]);
	}

	for (uint16 i = 0; i < size; i++) {//size; i++) {
		copy_population(temp_population, population, size);

		if(shuffle_array(temp_population, size) == data_no_error)
		{
			//get best k individuals from tournament
			uint16 index1 = get_index_of_BEST(temp_population, k);

			copy_individual_into_population(pool_population, temp_population, i,
					index1);
		}

	}

}

/*=========================================*/
/*This function implements Roulette Wheel selection*/
/*=========================================*/
static void selection_roulette_wheel(T_INDIVIDUAL2* population, uint16* size, uint16 model_size_n)
{

	T_INDIVIDUAL2 temp_population[*size];
	T_INDIVIDUAL2 temp_pool[*size];

	uint16 temp_size = 0;
	uint16 temp_pool_size = *size;

	do {

		probability_selection(population, *size);

		population_selected(population, temp_pool, *size, &temp_pool_size);

		if (temp_pool_size)

		{
			for (uint8 i = 0; i < temp_pool_size; i++) {
				if (temp_size < *size) {

					individual_init(&temp_population[temp_size]);
					copy_individual_into_population(temp_population, temp_pool,
							temp_size, i);
					++temp_size;
				} else {
					break;
				}
			}

		}

	} while (temp_size < *size);

	copy_population(population, temp_population, temp_size);

}

/*=========================================*/
/*This function implements Building- Blocks Hypothesis*/
/*=========================================*/
static boolean selection_building_blocks(T_INDIVIDUAL2* population, T_INDIVIDUAL2* schemas, uint16* size, uint16 size_s)
{

	T_INDIVIDUAL2 temp_population[*size];
	uint16 temp_size = 0;
	boolean ret_val = FALSE;

	for (uint16 i = 0; i < *size; i++) {

		for(uint16 j = 0; j < size_s; j++)
		{

			if (FALSE != individual_match_schema(&population[i], &schemas[j])) {

				//TODO
				//use royal road functions

				individual_init(&temp_population[temp_size]);
				copy_individual_into_population(temp_population, population,
						temp_size, i);
				temp_size++;
				break;
			}

		}

	}

	if (temp_size > 0) {

			copy_population(population, temp_population, temp_size);
			*size = temp_size;
			ret_val = TRUE;

	}

	return ret_val;
}

/*=========================================*/
/*This function apply crossover on 2 individuals on a random single position */
/*=========================================*/
static void crossover_1point_simple(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size) {

	uint16 r = 0, i;

	while (r == 0 ) {
		r = rand() % size;
	}

	for (i = r; i < size; i++) {
		double t = bit_individ_1[i];
		bit_individ_1[i] = bit_individ_2[i];
		bit_individ_2[i] = t;
	}

}

/*=========================================*/
/*This function apply crossover on 2 individuals on a random single position but keeping parent's number of genes*/
/*=========================================*/
static void crossover_1point(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size) {

	uint16 r = 0, i;
	boolean condition = FALSE; //keeping same number of genes

	while (condition == FALSE) {

		r = rand() % (size);

		if (get_number_locus(bit_individ_1 , r)
				== get_number_locus(bit_individ_2, r)) {
			condition = TRUE;
		}

	}

	if (r) {
		for (i = r; i < size; i++) {
			double t = bit_individ_1[i];
			bit_individ_1[i] = bit_individ_2[i];
			bit_individ_2[i] = t;
		}
	} else {


		vector_swap(bit_individ_1, bit_individ_2, size);
		//gsl_vector_swap(bit_individ_1, bit_individ_2);

	}

}

/*=========================================*/
/*This function apply crossover on 2 individuals using uniform_crossover method*/
/*=========================================*/
static void crossover_uniform(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size) {

	for (uint8 i = 0; i < size; i++) {

		double temp1 = bit_individ_1[i];
		double temp2 = bit_individ_2[i];

		if (temp1 != temp2) {
			//random number between 0 and 1
			double r = (double) (rand() / (double) (RAND_MAX));

			if (r > REFERENCE_PROBABILITY) {
				bit_individ_1[i] = temp2;
				bit_individ_2[i] = temp1;

			}

		}

	}
}

/*=========================================*/
/*This function apply crossover on 2 individuals using random respectful crossover (RRC) method*/
/*Genes that are same in both parents, are keep in both children */
/*And the rest of them are randomly chosen between children*/
/*=========================================*/
static void crossover_RRC(uint8* bit_individ_1, uint8 *bit_individ_2, uint16 size) {


	for (uint8 i = 0; i < size; i++) {

		double temp1 = bit_individ_1[i];
		double temp2 = bit_individ_2[i];

		if (temp1 != temp2) {
			//random number between 0 and 1
			double r = (double) (rand() / (double) (RAND_MAX));

			if (r < REFERENCE_PROBABILITY) {
				bit_individ_1[i] = 1;
			} else {
				bit_individ_1[i] = 0;
			}

			r = (double) (rand() / (double) (RAND_MAX));

			if (r < REFERENCE_PROBABILITY) {
				bit_individ_2[i] = 1;
			} else {
				bit_individ_2[i] = 0;
			}

		}

	}
}

/*=========================================*/
/*This function counts number of bits set on 1*/
/*=========================================*/
static uint16 get_number_locus(uint8 *bit_individ, uint16 size) {
	uint16 i, nr = 0;
	for (i = 0; i < size; i++) {
		if (bit_individ[i] == 1)
			nr++;
	}

	return nr;
}

/*=========================================*/
/*This function apply mutation on individual on a n random positions, n also randomly chosen*/
/*=========================================*/
static void mutation_flip(uint8* bit_columns, uint16 size) {
	// get n random positions and change its value : if it's 0 become 1, and reverse;

	uint8* aux = (uint8*) calloc(size, sizeof(uint8));

	uint16 n_random_positions = 0;

	while (n_random_positions == 0) {
		n_random_positions = rand() % (size);
	}

	while (n_random_positions) {
		copy_vector8(aux, bit_columns, size);
		uint8 r = rand() % (size);
		uint8 temp_el = aux[r];
		temp_el ^= 1UL;
		aux[r] = temp_el;
		n_random_positions--;

	}

	copy_vector8(bit_columns, aux, size);
	free(aux);
}


/*=========================================*/
/*This function apply mutation on individuals on a 2 random positions, by switching it's elements*/
/*=========================================*/
static void mutation_interchanging_abs(uint8* bit_columns, uint16 size) {

	uint16 pos1 = 0;
	uint16 pos2 = 0;

	while((pos1 == pos2) || (bit_columns[pos1] == bit_columns[pos2]))
	{
		pos1 = rand() % (size);
		pos2 = rand() % (size);

	}

	vector_swap_elements(bit_columns, pos1, pos2);
}

/*=========================================*/
/*This function apply mutation on individuals on a 2 random positions, by switchig it's elements*/
/*=========================================*/
static void mutation_interchanging (uint8* bit_columns, uint16 size) {

	uint16 pos1 = 0;
	uint16 pos2 = 0;

	while((pos1 == pos2))
	{
		pos1 = rand() % (size);
		pos2 = rand() % (size);

	}

	vector_swap_elements(bit_columns, pos1, pos2);
}

/*=========================================*/
/*This function apply mutation on individuals that are after one random position, by reversing it's bit value*/
/*=========================================*/
static void mutation_reversing (uint8* bit_columns, uint16 size) {

	uint16 pos = 0;
	uint8 temp_el = 0;

	pos = rand() % (size - 1);

	for (uint8 i = pos + 1; i < size; i++) {
		temp_el = bit_columns[i];
		temp_el ^= 1UL;
		bit_columns[i] = temp_el;
	}
}


/*=========================================*/
/*This function compute new population using selection method*/
/*=========================================*/
static void new_population_computed(T_INDIVIDUAL2* temp_population, uint16 new_size, uint16 columns, T_OPERATOR_METHOD op)
{

	uint8* temp_mutation = (uint8*) calloc(columns, sizeof(uint8));
	uint8* temp_crossover1 = (uint8*) calloc(columns, sizeof(uint8));
	uint8* temp_crossover2 = (uint8*) calloc(columns, sizeof(uint8));

	if (op < _1point_simple) {


		//Apply mutation on selected individuals
		for (uint16 i = 0; i < new_size; i++) {

			copy_vector8(temp_mutation, temp_population[i].bit_columns, columns);

			switch (op) {
			case flip: {mutation_flip(temp_mutation, columns);break;}

			case interchanging: {mutation_interchanging(temp_mutation, columns);break;}

			case interchanging_abs: {mutation_interchanging_abs(temp_mutation, columns);break;}

			case reversing: {mutation_reversing(temp_mutation, columns);break;}

			default:break;}

			build_individual(&temp_population[i], temp_mutation, columns);


		}


	}

	else if (op < no_operator) {

		if(new_size % 2 == 1){
			new_size--;
		}

		if(new_size > 1)
		{
			//Apply crossover on selected individuals
			for (uint8 i = 0; i < new_size; i += 2) {

				copy_vector8(temp_crossover1, temp_population[i].bit_columns, columns);
				copy_vector8(temp_crossover2, temp_population[i+1].bit_columns, columns);

				switch (op) {

				case _1point_simple: {crossover_1point_simple(temp_crossover1, temp_crossover2, columns);break;}

				case _1point: {crossover_1point(temp_crossover1, temp_crossover2, columns);break;}

				case RRC: {crossover_RRC(temp_crossover1, temp_crossover2, columns);break;}

				case uniform: {crossover_uniform(temp_crossover1, temp_crossover2, columns);break;}

				default:break;}


				build_individual(&temp_population[i], temp_crossover1, columns);
				build_individual(&temp_population[i + 1], temp_crossover2, columns);


			}
		}

	}

	else {
		/*nothing to do*/
	}

	free(temp_mutation);
	free(temp_crossover1);
	free(temp_crossover2);

}

/*=========================================*/
/*This function copy values of a source population into destination population */
/*=========================================*/
static void copy_population(T_INDIVIDUAL2* dest, T_INDIVIDUAL2* src, uint16 new_size) {
	for (uint8 i = 0; i < new_size; i++) {

		dest[i].size = src[i].size;
		dest[i].size_bit = src[i].size_bit;
		dest[i].RSS = src[i].RSS;
		dest[i].fitness_value = src[i].fitness_value;
		dest[i].selected = src[i].selected;
		dest[i].selection_probability = src[i].selection_probability;

		copy_vector16(dest[i].columns, src[i].columns, src[i].size);
		copy_vector8(dest[i].bit_columns, src[i].bit_columns, src[i].size_bit);
	}
}

/*=========================================*/
/*This function copy values of a source individual into another individual */
/*=========================================*/
static void copy_individual(T_INDIVIDUAL2* dest, T_INDIVIDUAL2* src) {

	dest->size = src->size;
	dest->size_bit = src->size_bit;
	dest->RSS = src->RSS;
	dest->fitness_value = src->fitness_value;
	dest->selected = src->selected;
	dest->selection_probability = src->selection_probability;

	copy_vector16(dest->columns, src->columns, src->size);
	copy_vector8(dest->bit_columns, src->bit_columns, src->size_bit);

}

/*=========================================*/
/*This function copy an individual from a population to another population*/
/*=========================================*/
static void copy_individual_into_population(T_INDIVIDUAL2* dest, T_INDIVIDUAL2* src, uint16 index1, uint16 index2) {


	if ((NULL != src[index2].bit_columns)) {

		dest[index1].size = src[index2].size;
		dest[index1].size_bit = src[index2].size_bit;
		dest[index1].RSS = src[index2].RSS;
		dest[index1].fitness_value = src[index2].fitness_value;
		dest[index1].selected = src[index2].selected;
		dest[index1].selection_probability = src[index2].selection_probability;

		copy_vector16(dest[index1].columns, src[index2].columns, src[index2].size);
		copy_vector8(dest[index1].bit_columns, src[index2].bit_columns, src[index2].size_bit);

	} else {
		printf("\ncopy_individual_into_population : index out of bound\n");
	}

}

/*=========================================*/
/*This function compute new population using selection method*/
/*=========================================*/
static void population_selected(T_INDIVIDUAL2* population, T_INDIVIDUAL2* temp_population, uint16 size, uint16* new_size) {

	uint16 temp_index = 0;

	for (uint8 i = 0; i < size; i++) {

		if (FALSE != population[i].selected) {

			individual_init(&temp_population[temp_index]);

			copy_individual_into_population(temp_population, population, temp_index, i);

			temp_index++;
		}
	}

	*new_size = temp_index;

}

/*=========================================*/
/*This function prints to the console all data of an individual*/
/*=========================================*/
static void print_individual(T_INDIVIDUAL2* individual)
{
	printf("\nColumns:");
	if (individual->columns != NULL && individual->size > 0) {
		print_vector16(individual->columns, individual->size);
	}

	printf("Fitness Value: %f\n", individual->fitness_value);
}

/*=========================================*/
/*This function prints to the console all data of all individuals in a population*/
/*=========================================*/
static void print_population(T_INDIVIDUAL2* population, uint16 size) {

	if (size > 0) {
		for (uint16 i = 0; i < size; i++) {

			print_individual(&population[i]);
		}
	} else {
		printf("\nPopulation size is 0\n");
	}
}

/*=========================================*/
/*This function calculate fitness function of each individual and return index of the best*/
/*=========================================*/
void fitness_func(T_INDIVIDUAL2* individual, uint16 model_size_n, uint16 model_size_k, double* result)
{
	criterion(individual->RSS, model_size_n, individual->size, result);

	//in each case update fitness value of each individual
	individual->fitness_value = *result;
}

/*=========================================*/
/*This function implements a method for acceptance of a neighbor */
/*=========================================*/
static boolean neighbor_acceptance(T_INDIVIDUAL2* current, T_INDIVIDUAL2* neighbor, double temperature){

	boolean is_accepted = FALSE;

	double delta = current->fitness_value - neighbor->fitness_value;

	if (delta >= 0) {
		is_accepted = TRUE;
	} else {
		if (temperature) {

			double r = (double) (rand() / (double) (RAND_MAX));

			if ((double) exp(delta / temperature) > r) {

				is_accepted = TRUE;

			}

		}
	}

	return is_accepted;

}

/*=========================================*/
/*This function calculates number of genes by random or fixed number*/
/*=========================================*/
static uint16 get_nr_genes(uint16 columns)
{
	uint16 ret_val = 0;

	if (FIXED_NR_GENES == 0 && columns > 0) {

		while (ret_val == 0) {
			ret_val = rand() % columns;
		}
	}

	else if (FIXED_NR_GENES == 1) {

		ret_val = (uint16) PERCENTAJE(PERCENTAJE_OF_GENES, columns);
	} else {
		printf("\nNumber of columns is 0, set a bigger percentaje");
	}

	return ret_val;
}

/*=========================================*/
/*This function calculates mean fitness of a population*/
/*=========================================*/
static double get_mean_fitness(T_INDIVIDUAL2* population, uint16 size)
{
	double ret_val = 0;


	for (uint16 i = 0; i < size; i++) {

		ret_val += population[i].fitness_value;

	}

	if (size > 0 && ret_val > 0) {
		ret_val /= size;
	}

	return ret_val;
}

static void generate_schemas(T_INDIVIDUAL2* schemas, uint16 size)
{

	T_INDIVIDUAL2 temp_schema;
	uint16 columns_schema = 0;

	individual_init(&temp_schema);

	if(size > 0)
	{
		columns_schema = (uint16)(temp_schema.size_bit/size);
	}

	for (uint16 i = 0; i < size; i++) {

		generate_individual(&temp_schema, temp_schema.size_bit, columns_schema);

		individual_init(&schemas[i]);
					copy_individual_into_population(schemas, &temp_schema,
							i, 0);

	    individual_dealloc(&temp_schema);
		individual_init(&temp_schema);
	}

}

/*=========================================*/
/*This function generates a number of schemas*/
/*=========================================*/
static void build_schemas(T_INDIVIDUAL2* schemas, uint16* size_schemas, T_INDIVIDUAL2* population, uint16 size)
{
	double mean_fitness = 0;
	double mean_schema_fitness = 0;
	double result = 1;

	T_INDIVIDUAL2 temp_schema;

	uint16 columns_schema;
	uint16 individuals_s;
	uint16 temp_i = 0;
	uint16 max_search = CONVERGE;

	individual_init(&temp_schema);

	if(*size_schemas > 0)
	{
		columns_schema = (uint16)(temp_schema.size_bit/ *size_schemas);
	}

	if(size > 0)
	{
		mean_fitness = get_mean_fitness(population, size);
	}

	if(columns_schema > 0 )
	{
		while (temp_i < *size_schemas && max_search > 0) {

			mean_schema_fitness = 0;
			individuals_s = 0;

			generate_individual(&temp_schema, temp_schema.size_bit, columns_schema);
			fitness_func(&temp_schema, temp_schema.size_bit, temp_schema.size, &result);

			for (uint16 i = 0; i < size; i++) {
				if (FALSE != individual_match_schema(&population[i], &temp_schema)) {
					mean_schema_fitness += population[i].fitness_value;

					individuals_s++;
				}
			}

			if(mean_schema_fitness > 0 && individuals_s > 0)
			{
				mean_schema_fitness /= individuals_s;
			}


			if((mean_schema_fitness <= mean_fitness) && (mean_schema_fitness > 0))
			//if(temp_schema.fitness_value < mean_fitness && (double)temp_schema.fitness_value > (double)0 )
			{
				individual_init(&schemas[temp_i]);
				copy_individual_into_population(schemas, &temp_schema, temp_i, 0);
				temp_i++;

			}
			else
			{
				max_search--;

			}

		}

		*size_schemas = temp_i;

	}
	else
	{
		*size_schemas = 0;
	}

}

/*=========================================*/
/*This function combine 2 schemas*/
/*=========================================*/
static void combine_schemas(T_INDIVIDUAL2* schema1,T_INDIVIDUAL2* schema2)
{
	double mean_fitness = 0;
	double mean_schema_fitness = 0;
	double result = 1;
	uint16 temp_size = 0;
	uint8 model[schema1->size_bit];

	T_INDIVIDUAL2 temp_schema;

	individual_init(&temp_schema);

	for(uint16 i = 0; i < schema2->size_bit; i++)
	{
		if(schema1->bit_columns[i] == 1 || schema2->bit_columns[i])
		{
			//temp_schema.columns[temp_size++] = i;
			temp_size++;
			model[i] = 1;
		}
		else
		{
			model[i] = 0;
		}



	}

	temp_schema.size = temp_size;

	build_individual(&temp_schema, model, schema1->size_bit);

	fitness_func(&temp_schema, temp_schema.size_bit, temp_schema.size, &result);

	if(temp_schema.fitness_value <= schema1->fitness_value)
	{
		schema1->size = temp_size;
		copy_individual(schema1, &temp_schema);
	}

	if(temp_schema.fitness_value <= schema1->fitness_value)
	{
		schema2->size = temp_size;
		copy_individual(schema2, &temp_schema);
	}


}

/*=========================================*/
/*This function implements a naive approach of genetic algorithm*/
/*=========================================*/
void GA_naive_alg(T_SELECTION_METHOD method, T_OPERATOR_METHOD op1,  T_OPERATOR_METHOD op2)
{

	gsl_matrix* main_model = get_A_matrix();
	uint16 converge_value = INIT;
	uint16 generation = INIT;
	uint16 model_size_n = main_model->size2;
	uint16 model_size_k = get_nr_genes(main_model->size2);
	uint16 population_size = (uint16)PERCENTAJE(PERCENTAJE_OF_CHROMOSOMES, main_model->size1);
	uint16 best_solution_index = INIT;
	boolean b_converge = FALSE;
	boolean last_converge = FALSE;
	T_INDIVIDUAL2 Schema;
	individual_init(&Schema);
	generate_individual(&Schema, model_size_n, 2);

	double result = 1;
	double MIN = (double) MAX_FITNESS;

	if(model_size_k == 0 || population_size == 0)
	{
		printf("\nPercentaje value of genes or population size is invalid, set a bigger value!");
		converge_value = 0;
	}
	else {

		T_INDIVIDUAL2 GA_population[population_size];
		T_INDIVIDUAL2 best_solution;

		individual_init(&best_solution);

		//generate population
		generate_population(GA_population, population_size, model_size_n,
				model_size_k);

		//get fitness of all individuals
		for (uint16 i = 0; i < population_size; i++) {
			fitness_func(&GA_population[i], model_size_n, model_size_k,
					&result);
		}
		//get BEST
		best_solution_index = get_index_of_BEST(GA_population, population_size);

		if (GA_population[best_solution_index].fitness_value < MIN)

		{
			MIN = GA_population[best_solution_index].fitness_value;
			copy_individual(&best_solution,
					&GA_population[best_solution_index]);

		}

		//repeat until converge
		while (converge_value < CONVERGE)

		{

			//printf("\n==========#Generation %d#===========", ++generation);

			/*select best individual for next generation (x %)
			 * 1. Roulette Wheel -> how % from SUM of all Results represents each Result)
			 * 2. Tournament Selection -> first selects two individuals with uniform probability -> chooses the one with the highest fitness.
			 * 3. Truncation Selection -> simply selects at random from the population having first eliminated K number of the least fit individuals
			 */

			switch (method) {
			case tournament: selection_tournament(GA_population, population_size, PERCENTAJE(PERCENTAJE_OF_TOURNAMENT_K, model_size_k)); break;

			case roulette_wheel: selection_roulette_wheel(GA_population, &population_size, model_size_n); break;

			case building_blocks: if(FALSE != selection_building_blocks(GA_population, &Schema, &population_size, model_size_k))
								    {
										//selection_building_blocks(GA_population, &best_solution, &population_size, model_size_k);
								    } break;

			default: /*No method selected*/	break;

			}


			new_population_computed(GA_population, population_size,
					model_size_n, op2);
			new_population_computed(GA_population, population_size,
					model_size_n, op1);


			if (population_size) {

				//apply fitness to all individuals
				for (uint16 i = 0; i < population_size; i++) {
					fitness_func(&GA_population[i], model_size_n, model_size_k,
							&result);
				}

				//print_population(GA_population, population_size);

				//get BEST
				best_solution_index = get_index_of_BEST(GA_population,
						population_size);

				if (GA_population[best_solution_index].fitness_value < MIN)

				{
					MIN = GA_population[best_solution_index].fitness_value;
					copy_individual(&best_solution,
							&GA_population[best_solution_index]);

					b_converge = FALSE;
					//reset counter in case converge is not consecutively met
					converge_value = INIT;

				} else {
					b_converge = TRUE;

					if (last_converge) {
						converge_value++;
					}
				}

			} else {
				//if no individuals got selected in new population then search must be stopped
				converge_value = CONVERGE;
			}

			last_converge = b_converge;

		}

		printf("\nBEST");
		print_individual(&best_solution);

		individual_dealloc(&best_solution);
	}

}


/*=========================================*/
/*This function implements a naive approach of genetic algorithm with Building Blocks */
/*=========================================*/
void GA_BB_alg(T_OPERATOR_METHOD op1,  T_OPERATOR_METHOD op2)
{

	uint16 converge_value = INIT;
	uint16 generation = INIT;
	uint16 model_size_n = get_A_matrix()->size2;
	uint16 model_size_k = get_nr_genes(model_size_n);
	uint16 population_size = (uint16)PERCENTAJE(PERCENTAJE_OF_CHROMOSOMES, get_A_matrix()->size1);
	uint16 temp_population_size, old_population_size = population_size;
	uint16 nr_schemas = (uint16)PERCENTAJE(PERCENTAJE_OF_SCHEMAS, model_size_n);
	uint16 best_solution_index = INIT;
	uint16 best_solution_index2 = INIT;
	boolean b_converge = FALSE;
	boolean last_converge = FALSE;
	boolean correct_size = TRUE;
	uint16 temp_nr_schemas = nr_schemas;
	uint16 temp_nr_schemas_old = nr_schemas;


	double result = 1;
	double MIN = (double) MAX_FITNESS;

	if (model_size_k == 0 || population_size == 0) {
		printf(
				"\nPercentaje value of genes or population size is invalid, set a bigger value!");
		converge_value = 0;
		correct_size = FALSE;
	}

	if (nr_schemas == 0) {
		printf(
				"\nPercentaje value of schemas size is invalid, set a bigger value!");
		converge_value = 0;
		correct_size = FALSE;
	}

	if(FALSE != correct_size)
	{

		T_INDIVIDUAL2 GA_population[population_size];
		T_INDIVIDUAL2 GA_temp_population[temp_population_size];
		T_INDIVIDUAL2 GA_old_population[old_population_size];
		T_INDIVIDUAL2 best_solution;
		T_INDIVIDUAL2 Schemas[nr_schemas];
		T_INDIVIDUAL2 temp_Schemas[nr_schemas];
		uint16 min_size = INIT;
		uint16 max_size = INIT;

		individual_init(&best_solution);

		//generate population
		generate_population(GA_population, population_size, model_size_n,
				model_size_k);

		//get fitness of all individuals
		for (uint16 i = 0; i < population_size; i++) {
			fitness_func(&GA_population[i], model_size_n, model_size_k,
					&result);
		}

		//get BEST
		best_solution_index = get_index_of_BEST(GA_population, population_size);

		//set BEST
		set_BEST(&best_solution, &GA_population[best_solution_index], &MIN);

		build_schemas(Schemas, &temp_nr_schemas, GA_population, population_size);


		while (converge_value < CONVERGE)

		{
			generation++;

/*			if (FALSE
					!= selection_building_blocks(GA_population, Schemas,
							&population_size, temp_nr_schemas)) {
			}*/

			new_population_computed(GA_population, population_size,
					model_size_n, op2);
			new_population_computed(GA_population, population_size,
					model_size_n, op1);

			if (population_size) {

				//apply fitness to all individuals
				for (uint16 i = 0; i < population_size; i++) {
					fitness_func(&GA_population[i], model_size_n, model_size_k,
							&result);
				}

				//get BEST
				best_solution_index = get_index_of_BEST(GA_population,
						population_size);



				if (GA_population[best_solution_index].fitness_value < MIN)

				{
					set_BEST(&best_solution, &GA_population[best_solution_index], &MIN);

					b_converge = FALSE;
					//reset counter in case converge is not consecutively met
					converge_value = INIT;

				} else {
					b_converge = TRUE;

					if (last_converge) {
						converge_value++;
					}
				}

				temp_nr_schemas_old = temp_nr_schemas;
				temp_nr_schemas = nr_schemas;
				build_schemas(temp_Schemas, &temp_nr_schemas, GA_population,
						population_size);

				if(temp_nr_schemas_old < temp_nr_schemas )
				{
					min_size = temp_nr_schemas_old;
					max_size = temp_nr_schemas;
				}
				else
				{
					max_size = temp_nr_schemas_old;
					min_size = temp_nr_schemas;
				}


				for (uint16 i = 0; i < min_size - 1; i++) {

					combine_schemas(&temp_Schemas[i], &temp_Schemas[i + 1]);

				}

				for (uint16 i = 0; i < min_size; i++) {

					if (Schemas[i].fitness_value < temp_Schemas[i].fitness_value) {
						//Schema[i] doesn't change

					} else {
						copy_individual_into_population(Schemas, temp_Schemas,
								i, i);
					}

				}

				for(uint16 i= min_size; i< max_size; i++)
				{
					if(max_size == temp_nr_schemas_old)
					{
						/*nothing*/
					}
					else
					{
						copy_individual_into_population(Schemas, temp_Schemas,
								i, i);
					}
				}


				best_solution_index2 = get_index_of_BEST(Schemas,
						temp_nr_schemas);

				if (Schemas[best_solution_index2].fitness_value < MIN) {
					set_BEST(&best_solution, &Schemas[best_solution_index2], &MIN);

					b_converge = FALSE;
					//reset counter in case converge is not consecutively met
					converge_value = INIT;

				} else {
					b_converge = TRUE;

					if (last_converge) {
						converge_value++;
					}
				}


			} else {
				//if no individuals got selected in new population then search must be stopped
				converge_value = CONVERGE;
			}

			last_converge = b_converge;

		}

		printf("\nBEST");
		print_individual(&best_solution);

		individual_dealloc(&best_solution);

	}

}


/*=========================================*/
/*This function implements a constructor for individual*/
/*=========================================*/
void individual_init(T_INDIVIDUAL2* individual)
{
	uint16 columns = get_A_matrix()->size2;
	individual->size = 0;
	individual->size_bit = columns;
	individual->columns = (uint16*) calloc(columns, sizeof(uint16));
	individual->bit_columns = (uint8*) calloc(columns, sizeof(uint8));
	individual->RSS = 0;
	individual->fitness_value = MAX_FITNESS;
	individual->selected = FALSE;
	individual->selection_probability = 0;
}

/*=========================================*/
/*This function implements a de-constructor for individual*/
/*=========================================*/
void individual_dealloc(T_INDIVIDUAL2* individual)
{
	individual->size = 0;
	individual->size_bit = 0;
	individual->RSS = 0;
	individual->fitness_value = MAX_FITNESS;
	individual->selected = FALSE;
	individual->selection_probability = 0;
	free(individual->columns);
	free(individual->bit_columns);
}

/*=========================================*/
/*This function implements a naive approach of genetic algorithm*/
/*=========================================*/
void GA_simulated_annealing(T_OPERATOR_METHOD op1, T_OPERATOR_METHOD op2){

	uint16 model_size_n = get_A_matrix()->size2;
	uint16 model_size_k = get_nr_genes(model_size_n);
	uint16 interation = 0;

	T_INDIVIDUAL2 current_individual;
	T_INDIVIDUAL2 neighbor_individual;
	T_INDIVIDUAL2 best_individual;

	individual_init(&current_individual);
	individual_init(&neighbor_individual);
	individual_init(&best_individual);

	double result = 1;
	double MIN_fitness;
	double temperature = TEMP;

	if (model_size_k == 0) {
		printf("\nPercentaje value of genes is invalid, set a bigger value!");
		interation = NR_ITERATIONS;
	} else {

		generate_individual(&current_individual, model_size_n, model_size_k);

		fitness_func(&current_individual, model_size_n, model_size_k, &result);

		MIN_fitness = current_individual.fitness_value;

		copy_individual(&best_individual, &current_individual);

	}

	for (; interation < NR_ITERATIONS; interation++) {

		copy_individual(&neighbor_individual, &current_individual);

		new_population_computed(&neighbor_individual, 1, model_size_n, op1);

		fitness_func(&neighbor_individual, model_size_n, model_size_k, &result);

		temperature = (double) temperature * COOLING_RATE;

		if (neighbor_acceptance(&current_individual, &neighbor_individual,
				temperature)) {

			copy_individual(&current_individual, &neighbor_individual);

		}

		if (current_individual.fitness_value < MIN_fitness) {
			MIN_fitness = current_individual.fitness_value;

			copy_individual(&best_individual, &current_individual);

			printf("\n==========#Iteration  %d#===========", interation);
			print_individual(&neighbor_individual);
		}
	}


	printf("\nBEST");
	print_individual(&best_individual);

	individual_dealloc(&current_individual);
	individual_dealloc(&neighbor_individual);
	individual_dealloc(&best_individual);

}

