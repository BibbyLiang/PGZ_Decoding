#ifndef PGZ_DECODING_H
#define PGZ_DECODING_H

#define GF_FIELD        8
#define MESSAGE_LEN     3
#define CODEWORD_LEN    GF_FIELD - 1    

extern unsigned char power_polynomial_table[GF_FIELD][2];
extern unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1];
extern unsigned char message_polynomial[MESSAGE_LEN];
extern unsigned char received_polynomial[CODEWORD_LEN];

extern int pgz_decoding();
#endif
