/*
 ============================================================================
 Name        : least_square_c.c
 Author      : Lucas Cesar
 Version     :
 Copyright   : 
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


/* Debug definitions start */
#define ACTIVATE_PRINT_ERRORS       1
#define ACTIVATE_PRINT_VALUES_L1    1
#define ACTIVATE_PRINT_VALUES_L2    1


#if (ACTIVATE_PRINT_VALUES_L1)
#define PRINT_L1(...)         printf(__VA_ARGS__)
#else
#define PRINT_L1(...)         ( (void) 0 )
#endif
#if (ACTIVATE_PRINT_VALUES_L2)
#define PRINT_DEBUG_L2(...)         printf(__VA_ARGS__)
#else
#define PRINT_DEBUG_L2(...)         ( (void) 0 )
#endif
#if (ACTIVATE_PRINT_ERRORS)
#define PRINT_ERRORS(...)         printf(__VA_ARGS__)
#else
#define PRINT_ERRORS(...)         ( (void) 0 )
#endif

/* Define Magical Numbers */
#define MAX_SIZE_MATRIX 3

/* Define Error codes */
#define ERR_OK 0
#define ERR_FAIL 1

/* Define types of method */
#define LINEAR 1
#define SQUARE 2
#define ROBUST 3
#define MULTIPLE 4

/* Define types of matrix */
#define TWO_DIMENSION 2
#define THREE_DIMENSION 3


/* Declare a struct that contains the values readed from files. */
struct auxiliar_values
{
	int num_of_values;
	double x[MAX_SIZE_MATRIX];
	double y;
	double weight;

};

struct base_matrices
{
	double x_t_x[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX];
	double inv_x_t_x[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX];
	double x_t_y[MAX_SIZE_MATRIX][1];
	double beta[MAX_SIZE_MATRIX][1];
};

void print_matrices(struct base_matrices* m)
{
	/* Print the results */
	PRINT_DEBUG_L2("\r\nMatrix XtX:\r\n");
	for(int i = 0; i < MAX_SIZE_MATRIX; i++)
	{
		for (int j = 0; j < MAX_SIZE_MATRIX; j++)
		{
			PRINT_DEBUG_L2("%.2lf",m->x_t_x[i][j]);
			if(j < MAX_SIZE_MATRIX - 1) PRINT_DEBUG_L2("\t");
			else PRINT_DEBUG_L2("\r\n");
		}
	}

	PRINT_DEBUG_L2("\r\nMatrix Xty:\r\n");
	for(int i = 0; i < MAX_SIZE_MATRIX; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			PRINT_DEBUG_L2("%.2lf",m->x_t_y[i][j]);
			if(j < MAX_SIZE_MATRIX - 1) PRINT_DEBUG_L2("\t");
			else PRINT_DEBUG_L2("\r\n");
		}
	}


	/* Print the results */
	PRINT_DEBUG_L2("\r\nInverse Matrix: \r\n");
	for(int i = 0; i < MAX_SIZE_MATRIX; i++)
	{
		for (int j = 0; j < MAX_SIZE_MATRIX; j++)
		{
			PRINT_DEBUG_L2 ("%.2lf",m->inv_x_t_x[i][j]);
			if(j < MAX_SIZE_MATRIX - 1) PRINT_DEBUG_L2("\t");
			else PRINT_DEBUG_L2("\r\n");
		}
	}

	/* Print the results */
	PRINT_DEBUG_L2("\r\nBeta: \r\n");
	for(int i = 0; i < MAX_SIZE_MATRIX; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			PRINT_DEBUG_L2 ("%.2lf",m->beta[i][j]);
			if(j < MAX_SIZE_MATRIX - 1) PRINT_DEBUG_L2("\t");
			else PRINT_DEBUG_L2("\r\n");
		}
	}
}

int inverse_matrix(double input[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX], double output[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX], int column, int row)
{
	double pivot = 0;
	double m = 0;
	double backup[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX];

	/* Defining output as Identity Matrix and fill backup with input */
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < column; j++)
		{
			if(i == j)
			{
				output[i][j] = 1;
			}
			else
			{
				output[i][j] = 0;
			}

			backup[i][j] = input[i][j];
		}
	}

	/* Iterate for each column in the matrix */
	for(int j = 0; j < column; j++)
	{
		pivot = input[j][j];

		for(int k = 0; k < column; k++)
		{
			input[j][k] = (input[j][k])/(pivot);
			output[j][k] = (output[j][k])/(pivot);
		}
		/* Iterate for each row of the matrix */
		for(int i = 0; i < row; i++)
		{
			if(i != j)
			{
				/* Find the coefficient */
				m = input[i][j];

				/* Iterate for each item on the respective row */
				for(int k = 0; k < column; k++)
				{
					/* Subtract the above row multiplied by the coefficient. */
					input[i][k] = (input[i][k]) - (m * input[j][k]);
					output[i][k] = (output[i][k]) - (m * output[j][k]);
				}
			}
		}
	}

	for(int i = 0; i < row; i++)
		{
			for(int j = 0; j < column; j++)
			{
				input[i][j] = backup[i][j];
			}
		}

	return ERR_OK;
}

int multiply_matrices(double first[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX], int row_fst, int column_fst, double second[MAX_SIZE_MATRIX][1], int row_sec, int column_sec, double result[MAX_SIZE_MATRIX][1])
{
	double sum = 0.0;

	/* Check function entries */
	if(first == 0 || second == 0 || result == 0)
	{
		PRINT_ERRORS("Invalid Entry!");
		return ERR_FAIL;
	}

	if(column_fst != row_sec)
	{
		PRINT_ERRORS("Invalid Entry! Column of First Matrix must be equal to row of Second Matrix!");
		return ERR_FAIL;
	}

	/* iterate for each column of the second matrix */
	for(int i = 0; i<column_sec; i++)
	{
		/* iterate for each row of the first matrix */
		for(int j = 0; j<row_fst; j++)
		{
			sum = 0.0;

			/* iterate for each column of the first matrix (equal to line on second matrix) */
			for(int h = 0; h<column_fst; h++)
			{
				sum += first[j][h]*second[h][i];
			}

			result[j][i] = sum;
		}
	}

	return ERR_OK;
}

int zero_matrices(struct base_matrices* m)
{
	/* Check function entries */
	if(m == 0)
	{
		PRINT_ERRORS("Invalid Entry!");
		return ERR_FAIL;
	}

	for(int i = 0; i < MAX_SIZE_MATRIX; i++)
	{
		/* Zero values from Beta */
		m->beta[i][0] = 0.0;
		m->x_t_y[i][0] = 0.0;

		for(int j = 0; j < MAX_SIZE_MATRIX; j++)
		{
			/* Zero values from Inverse of XtX, XtX and Xty */
			m->inv_x_t_x[i][j] = 0.0;
			m->x_t_x[i][j] = 0.0;
		}
	}
	return ERR_OK;
}

int input_t_output_matrix_step(int dimension, struct auxiliar_values* aux, double output[MAX_SIZE_MATRIX][1], bool is_weighted)
{
	/* Check all the entry values */
	if(dimension != TWO_DIMENSION && dimension != THREE_DIMENSION)
	{
		return ERR_FAIL;
	}

	/* Iterate for each Row and Column */
	for(int line = 0; line < dimension; line++)
	{
		if(is_weighted)
		{
			output[line][0] += aux->x[line] * aux->y * aux->weight;
		}
		else
		{
			output[line][0] += aux->x[line] * aux->y;
		}
	}

	return ERR_OK;
}

int gram_matrix_step(int dimension, struct auxiliar_values* aux, double gram_matrix[MAX_SIZE_MATRIX][MAX_SIZE_MATRIX], bool is_weighted)
{
	/* Check all the entry values */
	if(dimension != TWO_DIMENSION && dimension != THREE_DIMENSION)
	{
		return ERR_FAIL;
	}

	/* Iterate for each Row and Column */
	for(int line = 0; line < dimension; line++)
	{
		for(int column = 0; column < dimension; column++)
		{
			if(is_weighted)
			{
				gram_matrix[line][column] += aux->x[line] * aux->x[column] * aux->weight;
			}
			else
			{
				gram_matrix[line][column] += aux->x[line] * aux->x[column];
			}
		}
	}

	return ERR_OK;
}

int read_data_files(int type, char* file_name, struct auxiliar_values* a, struct base_matrices* m)
{
	FILE *fp;

	/* Check function entries */
	if(a == 0 || file_name == 0)
	{
		PRINT_ERRORS("Invalid Entry!");
		return ERR_FAIL;
	}

	/* Declare auxiliar variables to contain the values from text file */
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double w = 0.0;

	/* Initiate the variables with zero */
	a->num_of_values = 0;

	/* Open the file on read mode */
	fp = fopen(file_name, "r");

	/* Check how much values each line will contain according to dimensions value */
	if(type == LINEAR)
	{
		/* Read each line until find EOF value. For each line do the calculation of the values in the matrix. */
		do
		{
			x = 0.0;
			y = 0.0;

			fscanf(fp,"%lf;%lf", &x, &y);

			/* Fill auxiliary matrix values */
			a->x[0]=1;
			a->x[1]=x;
			a->y=y;

			PRINT_DEBUG_L2("x%d=%f y%d=%f\n\r",a->num_of_values, x, a->num_of_values, y);
			a->num_of_values += 1;

			input_t_output_matrix_step(TWO_DIMENSION, a, m->x_t_y, false);
			gram_matrix_step(TWO_DIMENSION, a, m->x_t_x, false);

		}while(!feof(fp));
	}
	else if(type == SQUARE)
	{
		/* Read each line until find EOF value. For each line do the calculation of the values in the matrix. */
		do
		{
			x = 0.0;
			y = 0.0;

			fscanf(fp,"%lf;%lf", &x, &y);

			/* Fill auxiliary matrix values */
			a->x[0]=1;
			a->x[1]=x;
			a->x[2]=x*x;
			a->y=y;

			PRINT_DEBUG_L2("x%d=%f y%d=%f\n\r",a->num_of_values, x, a->num_of_values, y);
			a->num_of_values += 1;

			input_t_output_matrix_step(THREE_DIMENSION, a, m->x_t_y, false);
			gram_matrix_step(THREE_DIMENSION, a, m->x_t_x, false);

		}while(!feof(fp));
	}
	else if(type == ROBUST)
	{
		/* Read each line until find EOF value. For each line do the calculation of the values in the matrix. */
		do
		{
			x = 0.0;
			y = 0.0;
			w = 0.0;

			fscanf(fp,"%lf;%lf;%lf", &w, &x, &y);

			/* Fill auxiliary matrix values */
			a->x[0]=1;
			a->x[1]=x;
			a->weight=w;
			a->y=y;

			PRINT_DEBUG_L2("x%d=%f y%d=%f\n\r",a->num_of_values, x, a->num_of_values, y);
			a->num_of_values += 1;

			input_t_output_matrix_step(TWO_DIMENSION, a, m->x_t_y, true);
			gram_matrix_step(TWO_DIMENSION, a, m->x_t_x, true);

		}while(!feof(fp));
	}
	else if(type == MULTIPLE)
	{
		/* Read each line until find EOF value. For each line do the calculation of the values in the matrix. */
		do
		{
			x = 0.0;
			y = 0.0;
			z = 0.0;

			fscanf(fp,"%lf;%lf;%lf", &x, &z, &y);

			/* Fill auxiliary matrix values */
			a->x[0]=1;
			a->x[1]=x;
			a->x[2]=z;
			a->y=y;

			PRINT_DEBUG_L2("x%d=%f y%d=%f\n\r",a->num_of_values, x, a->num_of_values, y);
			a->num_of_values += 1;

			input_t_output_matrix_step(THREE_DIMENSION, a, m->x_t_y, false);
			gram_matrix_step(THREE_DIMENSION, a, m->x_t_x, false);

		}while(!feof(fp));
	}
	else
	{
		PRINT_ERRORS("Invalid Operation!");
		fclose(fp);
		return ERR_FAIL;
	}

	/* Print the results */
	PRINT_DEBUG_L2("\n\rNumber of lines: %d\n\r", a->num_of_values);
	PRINT_DEBUG_L2("X value Sum    : %.2lf\n\r", a->x_sum);
	PRINT_DEBUG_L2("Y value Sum    : %.2lf\n\r", a->y_sum);
	PRINT_DEBUG_L2("X * Y value Sum: %.2lf\n\r", a->xy_sum);
	PRINT_DEBUG_L2("X * X value Sum: %.2lf\n\n\r", a->x_square_sum);

	fclose(fp);
	return ERR_OK;
}

int main(void)
{
	char key = 'a';
	char nome[50];
	int type = 0;
	int tries = 0;
	int test_value = 0;
	int test_value_2 = 0;


	struct auxiliar_values a;
	struct base_matrices matrices;

	PRINT_L1("-------------------\n\r");
	PRINT_L1("Least Square Method\n\r");
	PRINT_L1("-------------------\n\r");

	/* Initialize all the matrices values in zero */
	zero_matrices(&matrices);

	PRINT_L1("\n\rEnter the data source file:\n\r");
	scanf("%s",nome);
	PRINT_L1("\n\rFile: %s\n\n\r",nome);

	PRINT_L1("\n\rType the method:\n\r1 - Linear\n\r2 - Square\n\r3 - Robust (Weighted)\n\r4 - Two Variables\n\r");
	scanf("%d",&type);
	PRINT_L1("\n\rMethod: %d\n\r", type);

	while(type!=1 && type!=2 && type!=3 && type!=4 && tries<4)
	{
		PRINT_L1("\n\rWrong Entry! Type again.\n\r");
		fflush(stdin);
		scanf("%d",&type);

		tries++;
	}

	if(tries >= 4)
	{
		PRINT_ERRORS("\n\rWrong Entry!\n\r");
		sleep(5);
		return ERR_FAIL;
	}

	/* Read all data from the archive and prepare the data to populate the matrices */
	read_data_files(type, nome, &a, &matrices);

	if(type == LINEAR)
	{
		/* Reverse the matrix XtX */
		inverse_matrix(matrices.x_t_x, matrices.inv_x_t_x, 2, 2);

		/* Multiply Inverse XtX by Xty */
		multiply_matrices(matrices.inv_x_t_x, 2, 2, matrices.x_t_y, 2, 1, matrices.beta);

		PRINT_L1("\n\rEnter a test value:\n\r");
		scanf("%d",&test_value);
		PRINT_L1("\n\rResult: %f\n\r", matrices.beta[0][0]+test_value*matrices.beta[1][0]);
	}
	else if(type == SQUARE)
	{
		/* Reverse the matrix XtX */
		inverse_matrix(matrices.x_t_x, matrices.inv_x_t_x, 3, 3);

		/* Multiply Inverse XtX by Xty */
		multiply_matrices(matrices.inv_x_t_x, 3, 3, matrices.x_t_y, 3, 1, matrices.beta);

		PRINT_L1("\n\rEnter a test value:\n\r");
		scanf("%d",&test_value);
		PRINT_L1("\n\rResult: %lf\n\r", matrices.beta[0][0]+test_value*matrices.beta[1][0]+test_value*test_value*matrices.beta[2][0]);
	}
	else if(type == ROBUST)
	{
		/* Reverse the matrix XtX */
		inverse_matrix(matrices.x_t_x, matrices.inv_x_t_x, 2, 2);

		/* Multiply Inverse XtX by Xty */
		multiply_matrices(matrices.inv_x_t_x, 2, 2, matrices.x_t_y, 2, 1, matrices.beta);

		PRINT_L1("\n\rEnter a test value:\n\r");
		scanf("%d",&test_value);
		PRINT_L1("\n\rResult: %f\n\r", matrices.beta[0][0]+test_value*matrices.beta[1][0]);
	}
	else if(type == MULTIPLE)
	{
		/* Reverse the matrix XtX */
		inverse_matrix(matrices.x_t_x, matrices.inv_x_t_x, 3, 3);

		/* Multiply Inverse XtX by Xty */
		multiply_matrices(matrices.inv_x_t_x, 3, 3, matrices.x_t_y, 3, 1, matrices.beta);

		PRINT_L1("\n\rEnter a first test value:\n\r");
		scanf("%d",&test_value);
		PRINT_L1("\n\rEnter a second test value:\n\r");
		scanf("%d",&test_value);
		PRINT_L1("\n\rResult: %lf\n\r", matrices.beta[0][0]+test_value*matrices.beta[1][0]+test_value_2*matrices.beta[2][0]);
	}
	else
	{
		PRINT_ERRORS("Unknown Error!");
		return ERR_FAIL;
	}

	print_matrices(&matrices);

	PRINT_L1("\n\rPress \"q\" to quit!\n\r");

	do{
		scanf("%c",&key);
	}while(key!='q');

	return ERR_OK;
}
