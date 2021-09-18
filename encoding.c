#include <stdio.h>
#include <string.h>
#include "gf_cal.h"
#include "encoding.h"

unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1] =
{
	/*generated by MATLAB: genpoly = rsgenpoly(n, k, prim_poly)*/
	0x3,
	0x1,
	0x0,
	0x3,
	0x1
};

unsigned char message_polynomial[MESSAGE_LEN] = 
{
	/*power representation*/
	0x1,
	0x3,
	0x5
};

unsigned char encoded_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char error_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0x4,
	0x1,
	0xFF
};

unsigned char evaluation_encoding()
{
	unsigned char i = 0, j = 0;
	unsigned char tmp_power_result[MESSAGE_LEN], tmp_add_result = 0;

	unsigned char codeword[CODEWORD_LEN];

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_power_result[j] = gf_multp(message_polynomial[j], i * gf_location(j));
		}
		tmp_add_result = 0;
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_add_result = gf_pow2poly(tmp_power_result[j]) ^ tmp_add_result;
		}

		for(j = 0; j < GF_FIELD; j++)
		{
			codeword[i] = gf_poly2pow(tmp_add_result);
		}
	}

	printf("Evaluation Encoding Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", codeword[i]);
	}
	printf("\n");

	return 0;
}

unsigned char systematic_encoding()
{
	unsigned char i = 0;
	
	unsigned char clock_cycle = 0;
	unsigned char input_message = 0;
	unsigned char feed_back_reg = 0, reg_delay[CODEWORD_LEN - MESSAGE_LEN];
	unsigned char codeword[CODEWORD_LEN];

	/*for iterating calculation, tmp for C, not in HDL*/
	unsigned char feed_back_reg_prev = feed_back_reg, reg_delay_prev[CODEWORD_LEN - MESSAGE_LEN];

	memset(reg_delay, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	memset(codeword, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	if(0 == clock_cycle)
	{
		input_message = message_polynomial[MESSAGE_LEN - 1 - clock_cycle];
		feed_back_reg = gf_add(input_message, reg_delay[CODEWORD_LEN - MESSAGE_LEN - 1]);
		feed_back_reg_prev = feed_back_reg;
		memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	}
	for(clock_cycle = 1; clock_cycle < (CODEWORD_LEN - MESSAGE_LEN); clock_cycle++)
	{
		for(i = CODEWORD_LEN - MESSAGE_LEN - 1; i >0; i--)
		{
			reg_delay[i] = gf_add(gf_multp(feed_back_reg_prev, generator_polynomial[i]), reg_delay_prev[i - 1]);
		}
		reg_delay[0] = gf_multp(feed_back_reg_prev, generator_polynomial[0]);
		input_message = message_polynomial[MESSAGE_LEN - 1 - clock_cycle];		
		feed_back_reg = gf_add(input_message, reg_delay[CODEWORD_LEN - MESSAGE_LEN - 1]);

		feed_back_reg_prev = feed_back_reg;
		memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	}

	for(i = CODEWORD_LEN - MESSAGE_LEN; i < CODEWORD_LEN; i++)
	{
		codeword[i] = message_polynomial[i - (CODEWORD_LEN - MESSAGE_LEN)];
	}
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		codeword[i] = reg_delay[i];
	}
	memcpy(encoded_polynomial, codeword, sizeof(unsigned char) * CODEWORD_LEN);

	printf("Systematic Encoding Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", encoded_polynomial[i]);
	}
	printf("\n");

	return 0;
}
