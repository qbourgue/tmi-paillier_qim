#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

#define SECURITY_BITS 512

void rsaEncrypt(mpz_t N, mpz_t e, mpz_t m, mpz_t c){
	mpz_powm(c, m, e, N);
}

void rsaDecrypt(mpz_t N, mpz_t d, mpz_t c, mpz_t m) {
	mpz_powm(m, c, d, N);
}

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



int testRSA()
{
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	mpz_t N, e, d, m, c, p, q;
	mpz_inits(N, e, d, m, c, p, q, NULL);
	generatePrime(p, SECURITY_BITS, state);
	generatePrime(q, SECURITY_BITS, state);
	generateKeys(p, q, e, d, N, state);
	gmp_scanf("%Zu", m);
	gmp_printf("N is %Zu, p is %Zu, q is %Zu\n", N, p, q);
	rsaEncrypt(N, e, m, c);
	gmp_printf("Encrypting %Zu gives %Zu\n", m, c);
	rsaDecrypt(N, d, c, m);
	gmp_printf("Decrypting %Zu gives %Zu\n", c, m);
	gmp_randclear(state);
	mpz_clears(N, e, d, m, c, p, q, NULL);
 	return 0;
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

/*
int qim_lsb(mpz_t x, mpz_t m){
	if m =0
		return (X%2)*2
	if m=1
		return (X%2)*2 + 1
        return 0;
}
*/


int main(int argc, char* argv[]) {
	//testRSA();
	testPaillier();
	return 0;
}




