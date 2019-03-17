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

typedef enum
{
	not_installed = 0,
	installed
}T_Intercept;

/*=========================================*/
/*typedef*/
/*=========================================*/
/*external functions*/
extern T_Intercept intercept;
/*=========================================*/
#endif /* SRC_QR_SEARCHINGALGORITHM_H_ */
