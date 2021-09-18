#include <stdio.h>
#include <string.h>
#include "gf_cal.h"
#include "pgz_decoding.h"

unsigned char received_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char output_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

int pgz_decoding()
{
	unsigned char i = 0, j = 0, tmp = 0xFF, tmp_sum = 0xFF, lambda_root = 0;
	unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
	unsigned char lambda[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1];
	unsigned char omega[(CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN)];
	unsigned char omega_tmp[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1][(CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN - 1)];
	unsigned char err_location[CODEWORD_LEN];
	unsigned char lambda_dev[(CODEWORD_LEN - MESSAGE_LEN) / 2];
	unsigned char err_mag[CODEWORD_LEN];
	unsigned char codeword[CODEWORD_LEN];
	
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	memset(lambda, 0x0, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
	memset(omega, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN)));
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN); j++)
		{
			omega_tmp[i][j] = 0xFF;
		}
	}
	memset(err_location, 0x0, sizeof(unsigned char) * CODEWORD_LEN);
	memset(lambda_dev, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2));
	memset(err_mag, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(codeword, received_polynomial, sizeof(unsigned char) * CODEWORD_LEN);

	printf("Received:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", received_polynomial[i]);
	}
	printf("\n");

	/*compute the syndrome polynomial from received symbols*/
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			tmp = gf_multp(received_polynomial[j], (i + 1) * j);
			tmp_sum = gf_add(tmp, tmp_sum);
			//printf("%x %x\n", tmp, tmp_sum);
		}
		syndrome[i] = tmp_sum;
	}
	printf("Syndrome Polynomial:\n");
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		printf("%x ", syndrome[i]);
	}
	printf("\n");

	lambda[1] = gf_div(gf_add(gf_multp(syndrome[2], syndrome[1]), gf_multp(syndrome[0], syndrome[3])),
						gf_add(gf_multp(syndrome[1], syndrome[1]), gf_multp(syndrome[0], syndrome[2])));
	lambda[2] = gf_div(gf_add(gf_multp(syndrome[2], syndrome[2]), gf_multp(syndrome[1], syndrome[3])),
						gf_add(gf_multp(syndrome[1], syndrome[1]), gf_multp(syndrome[0], syndrome[2])));
	printf("Lambda Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		printf("%x ", lambda[i]);
	}
	printf("\n");

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN; j++)
		{
			omega_tmp[i][j] = gf_multp(lambda[i], syndrome[j]);
			//printf("%d %d: %x %x %x\n", i, j, lambda[i], syndrome[j], omega_tmp[i][j]);
		}
	}
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN; j++)
		{
			//printf("%d %d: %x %x\n", i, j, omega[i + j], omega_tmp[i][j]);
			omega[i + j] = gf_add(omega[i + j], omega_tmp[i][j]);
		}
	}
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		omega[i] = gf_multp(omega[i], gf_mod_single_term(i, CODEWORD_LEN - MESSAGE_LEN));
	}
	printf("Omega Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		printf("%x ", omega[i]);
	}
	printf("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		lambda_root = 0xFF;
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; j++)
		{
			//printf("%d %d: %x %x %x %x\n", i, j, lambda_root, lambda[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda[j], j * power_polynomial_table[i + 1][0]));
			lambda_root = gf_add(lambda_root, gf_multp(lambda[j], j * power_polynomial_table[i + 1][0]));
		}
		//printf("%d: %x %x\n", i, lambda_root, power_polynomial_table[i + 1][0]);
		if(0xFF == lambda_root)
		{
			err_location[GF_FIELD - 1 - i] = 1;/*inversion*/
		}
	}
	printf("Error Location:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", err_location[i]);
	}
	printf("\n");

	/*derivative of lambda*/
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2; i++)
	{
		if(0 != ((i + 1) % 2))
		{
			lambda_dev[i] = lambda[i + 1];
		}
	}
	printf("Lambda Derivative Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2; i++)
	{
		printf("%x ", lambda_dev[i]);
	}
	printf("\n");

	/*magnitude of error*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == err_location[i])
		{
			continue;
		}
		else
		{
			tmp = 0xFF;
			for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN - 1) / 2 + 1; j++)
			{
				//printf("%d %d: %x %x %x %x\n", i, j, tmp, omega[j], j * power_polynomial_table[i + 1][0], gf_multp(omega[j], j * power_polynomial_table[i + 1][0]));
				tmp = gf_add(tmp, gf_div(omega[j], j * power_polynomial_table[i + 1][0]));
			}
			tmp_sum = 0xFF;
			for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2; j++)
			{
				//printf("%d %d: %x %x %x %x\n", i, j, tmp_sum, lambda_dev[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda_dev[j], j * power_polynomial_table[i + 1][0]));
				tmp_sum = gf_add(tmp_sum, gf_div(lambda_dev[j], j * power_polynomial_table[i + 1][0]));
			}
			//printf("%d: %x %x\n", i, tmp, tmp_sum);
			err_mag[i] = gf_div(tmp, tmp_sum);
		}
	}
	printf("Error Magnitude:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", err_mag[i]);
	}
	printf("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		codeword[i] = gf_add(codeword[i], err_mag[i]);
	}
	printf("Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", codeword[i]);
	}
	printf("\n");
	
	return 0;
}
