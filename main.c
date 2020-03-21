#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

#define SECURITY_BITS 512


/****************************************************************************
 * Crypto & Watermarking project
 *
 *
 * Quentin Bourgue
 * Paul-Henri Mignot
 *
 * Data Hiding in Homomorphically Encrypted Medical Image
 *
 ****************************************************************************/

/****************************************************************************
 *
 * Paillier 
 *
 ****************************************************************************/

void generatePrime(mpz_t output, unsigned int bits, gmp_randstate_t state) {
	do {
		mpz_urandomb(output, state, bits);
		mpz_nextprime(output, output);
	}
	while (mpz_sizeinbase(output, 2) != bits);
}

void computePhi(mpz_t p, mpz_t q, mpz_t output) {
	mpz_t pMinusOne, qMinusOne;
	mpz_inits(pMinusOne, qMinusOne, NULL);
	mpz_sub_ui(pMinusOne, p, 1);
	mpz_sub_ui(qMinusOne, q, 1);
	mpz_mul(output, pMinusOne, qMinusOne);
	mpz_clears(pMinusOne, qMinusOne, NULL);
}

void generateKeys(mpz_t p, mpz_t q, mpz_t e, mpz_t d, mpz_t N, gmp_randstate_t state) {
	mpz_t phi, gcd;
	mpz_inits(phi, gcd, NULL);
	computePhi(p, q, phi);
	mpz_mul(N, p, q);
	do {
		mpz_urandomm(e, state, phi);
		mpz_gcd(gcd, e, phi);
	} while (mpz_cmp_ui(gcd, 1) != 0);
	mpz_invert(d, e, phi);
	mpz_t test1, test2;
	mpz_inits(test1, test2, NULL);
	mpz_mul(test2, d, e);
	mpz_mod(test1, test2, phi);
	gmp_printf("e * d mod phi = %Zu", test1);
	gmp_printf("e = %Zu, d = %Zu, phi = %Zu\n", e, d, phi);
	mpz_clears(phi, gcd, NULL);
}


void paillierEncrypt(mpz_t m, mpz_t r, mpz_t N, mpz_t c){
	mpz_t plus1, Nsquare, rn, nm ;
	mpz_inits(plus1, Nsquare, rn, nm, NULL);
	mpz_add_ui(plus1, N, 1);
	mpz_pow_ui(Nsquare, N, 2);
	mpz_powm(nm, plus1, m,  Nsquare);
	mpz_powm(rn, r, N, Nsquare);
	mpz_mul(c, nm, rn);
	mpz_mod(c, c, Nsquare);
	mpz_clears(plus1, Nsquare, rn, nm, NULL);
}

void paillierDecrypt(mpz_t c, mpz_t N, mpz_t m, mpz_t phi){
	mpz_t rInvert, product, Nsquare, r, rminus, minusN, temp, minusOne;
	mpz_inits(rInvert, product, Nsquare, r, rminus, minusN, temp, minusOne, NULL);
	mpz_set_si(minusOne, -1);
	mpz_pow_ui(Nsquare, N, 2);
	mpz_powm(temp, N, minusOne, phi);
	mpz_powm(r, c, temp, N);
	mpz_mul_si(minusN, N,-1);
	mpz_powm(rminus, r, minusN, Nsquare);
	mpz_mul(product, c, rminus);
	mpz_mod(product, product, Nsquare);
	mpz_sub_ui(product, product, 1);
	mpz_div(m, product, N);
	mpz_clears(rInvert, product, Nsquare, r, rminus, minusN, temp, minusOne, NULL);
}

int testPaillier() {
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	mpz_t N, m, c, r, p, q, phi, o;
	mpz_t m2, m3, c2, c3, r2, r3, N2;
	mpz_inits(N, m, c, r, p, q, phi, o, NULL);
	mpz_inits(m2, m3, c2, c3, r2, r3, N2, NULL);
	generatePrime(p, SECURITY_BITS, state);
	generatePrime(q, SECURITY_BITS, state);
	mpz_mul(N, p, q);
	computePhi(p, q, phi);
	printf("m : ");
	gmp_scanf("%Zu", m);
	printf("\nr : ");
	gmp_scanf("%Zu", r);
	printf("m2 : ");
	gmp_scanf("%Zu", m2);
	printf("\nr2 : ");
	gmp_scanf("%Zu", r2);
	paillierEncrypt(m, r, N, c);
	gmp_printf("\nEncrypting %Zu gives %Zu\n", m, c);
	paillierEncrypt(m2, r2, N, c2);
	gmp_printf("\nEncrypting %Zu gives %Zu\n", m, c);
	//paillierDecrypt(c, N, m, phi);
	//gmp_printf("Decrypting %Zu gives %Zu\n", c, m);
	mpz_t mul_e;
        mpz_init(mul_e);
	mpz_mul(mul_e, c, c2);
	mpz_pow_ui(N2, N, 2);
	mpz_mod(mul_e, mul_e, N2);
	gmp_printf("\nMUL c et c2: %Zu", mul_e);
	paillierDecrypt(mul_e, N, o, phi);
	gmp_printf("Decrypting %Zu gives %Zu\n", mul_e, o);
	
	mpz_add(m, m, m2);
	mpz_mul(r, r, r2);
	paillierEncrypt(m, r, N, c);
	gmp_printf("\n m1+m2, r1+r2 Encrypting %Zu, %Zu gives %Zu\n", m,r, c);
        mpz_clears(N, m, c, r, p, q, phi, o, NULL);
	mpz_clears(m2, m3, c2, c3, r2, r3, N2, NULL);
	mpz_clear(mul_e);
        return 0;
}

/****************************************************************************
 *
 * Data Generation
 *
 * (Create an array of V elements, which are random integers of 8 bytes)
 ****************************************************************************/

void data_generation(mpz_t * data, unsigned long int V){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));

	mpz_t random_bytes;
	mpz_inits(random_bytes, NULL);

	unsigned long int i = 0;
	unsigned int bitcnt = 8;
	do {
		mpz_urandomb(random_bytes, state, bitcnt);
		mpz_init_set(data[i], random_bytes);
		gmp_printf("random_bytes n°%d: %Zu\n", i, random_bytes);
		i++;
	}
	while (i<V);

	mpz_clears(random_bytes, NULL);
}


/****************************************************************************
 *
 * Pre-watermarking 
 *
 * (Create an array of p elements, which are random bits)
 *
 * LSB
 ****************************************************************************/

void watermark_generation(mpz_t *watermark, unsigned int p){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));

	mpz_t random_bit;
	mpz_inits(random_bit, NULL);

	unsigned int j = 0;
	unsigned int bitcnt = 1;
	do {
		mpz_urandomb(random_bit, state, bitcnt);
		mpz_init_set(watermark[j], random_bit);
		gmp_printf("random_bit n°%d: %Zu\n",j, random_bit);
		j++;
	}
	while (j<p);

	mpz_clears(random_bit, NULL);
}

void pre_watermarking(mpz_t * watermark, mpz_t * data, int N, int p){
	// embed the watermark
	mpz_t pixel, two;
	mpz_inits(pixel, two, NULL);
	mpz_set_ui(two, 2);
	int even;
	for (int in = 0; in<N; in++){
		for (int jp = 0; jp<p; jp++){
			mpz_set(pixel,data[in*p+jp]);
			gmp_printf("data[%d,%d]=%Zu\n",in, jp, pixel);
			// TO DO: get the LSB
			//even = mpz_divisible_p(pixel,two);
			//gmp_printf("even=%Zu\n",even);

			// TO DO: embet the watermark on the LSB
		}
	}
	mpz_clears(pixel, two, NULL);
}


/****************************************************************************
 *
 * Encryption 
 *
 ****************************************************************************/



/****************************************************************************
 *
 * Message embedding 
 *
 ****************************************************************************/



/****************************************************************************
 *
 * Message extraction
 *
 ****************************************************************************/



/****************************************************************************
 *
 * Data extraction
 *
 ****************************************************************************/



/****************************************************************************
 *
 * Main 
 *
 ****************************************************************************/

int main(int argc, char* argv[]) {
	// Init parameters
	unsigned int p = 2;
	unsigned int N = 5;
	unsigned long int V=p*N;
	
	// Data Generation
	mpz_t * data;
	data = malloc(sizeof(mpz_t) * V);
	data_generation(data, V);

	// Pre-watermarking
	mpz_t * watermark;
	watermark = malloc(sizeof(mpz_t) * p);
	watermark_generation(watermark, p);
	pre_watermarking(watermark, data, N, p);
	
	// Encryption
	
	// Message embedding
	
	// Mesage extraction
	
	// Data extraction

	// test main

	// clear/free
	free(data);
	free(watermark);
	return 0;
}
