#ifndef QR_SEARCHINGALGORITHM_H_
#define QR_SEARCHINGALGORITHM_H_

/*=========================================*/
/*include*/
/*=========================================*/
/*define*/
/*=========================================*/
/*enumerations*/
typedef enum
{
	naive_search = 0,
	efficient_search,
	GA_search
}T_SEARCH_STRATEGIES;
/*=========================================*/
/*typedef*/
/*=========================================*/
/*external functions*/
void get_base_R(gsl_matrix* R, gsl_vector* x, gsl_matrix* R_result);
/*=========================================*/
#endif /* SRC_QR_SEARCHINGALGORITHM_H_ */
