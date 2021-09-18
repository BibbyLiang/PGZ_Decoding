#ifndef GF_CAL_H
#define GF_CAL_H

#define GF_FIELD        8
#define MESSAGE_LEN     3
#define CODEWORD_LEN    GF_FIELD - 1    

extern unsigned char power_polynomial_table[GF_FIELD][2];

extern unsigned char gf_pow2poly(unsigned char val_in_pow);
extern unsigned char gf_poly2pow(unsigned char val_in_poly);
extern unsigned char gf_location(unsigned char val);
extern unsigned char gf_add(unsigned char a, unsigned char b);
extern unsigned char gf_multp(unsigned char a, unsigned char b);
extern unsigned char gf_div(unsigned char a, unsigned char b);
extern unsigned char gf_mod_single_term(unsigned char a, unsigned char b);

#endif
