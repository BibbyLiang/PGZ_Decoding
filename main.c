#include "gf_cal.h"
#include "pgz_decoding.h"
#include "encoding.h"

void main()
{
	unsigned char i = 0;

	systematic_encoding();

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
	}
	
	pgz_decoding();

	return;
}
