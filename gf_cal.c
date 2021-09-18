#include <stdio.h>
#include <string.h>
#include "gf_cal.h"

unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x3},
	{0x4, 0x6},
	{0x5, 0x7},
	{0x6, 0x5}
};

unsigned char gf_pow2poly(unsigned char val_in_pow)
{
	unsigned char val_in_poly = 0;

	return power_polynomial_table[val_in_pow + 1][1];
}

unsigned char gf_poly2pow(unsigned char val_in_poly)
{
	unsigned char i = 0;
	unsigned char val_in_pow = 0;

	for(i = 0; i < GF_FIELD; i++)
	{
		if(power_polynomial_table[i][1] == val_in_poly)
		{
			val_in_pow = power_polynomial_table[i][0];
			break;
		}
	}

	return val_in_pow;
}

unsigned char gf_location(unsigned char val)
{
	unsigned char val_location = power_polynomial_table[val + 1][0];

	return val_location;
}

unsigned char gf_add(unsigned char a, unsigned char b)
{
	unsigned char i = 0;
	unsigned char sum_in_pow = 0;
	
	unsigned char sum_in_poly = gf_pow2poly(a) ^ gf_pow2poly(b);
	sum_in_pow = gf_poly2pow(sum_in_poly);

	return sum_in_pow;
}

unsigned char gf_multp(unsigned char a, unsigned char b)
{
	if((0xFF == a) || (0xFF == b))
	{
		return 0xFF;
	}

	unsigned char product_in_pow = (a + b) % (GF_FIELD - 1);

	return product_in_pow;
}

unsigned char gf_div(unsigned char a, unsigned char b)
{
	if(0xFF == a)
	{
		return 0xFF;
	}
	if(0xFF == b)
	{
		printf("div err.\n");
		return 0xFF;
	}

	//printf("div: %x %x\n", a, b);
	unsigned char quotient_in_pow = 0;
	if(a >= b)
	{
		quotient_in_pow = (a - b) % (GF_FIELD - 1);
	}
	else
	{
		quotient_in_pow = ((b / (GF_FIELD - 1) + 1) * (GF_FIELD - 1) + a - b) % (GF_FIELD - 1);
	}

	return quotient_in_pow;
}

unsigned char gf_mod_single_term(unsigned char a, unsigned char b)
{
	unsigned char remainder_in_pow = 0xFF;
	if(a >= b)
	{
		remainder_in_pow = 0xFF;
	}
	else
	{
		remainder_in_pow = 0;
	}

	return remainder_in_pow;
}
