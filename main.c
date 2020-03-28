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

// Baptiste and Vincent Decrypt
/*
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
*/

void paillierDecrypt(mpz_t cipher, mpz_t N, mpz_t output, mpz_t phi){
	mpz_t Nsquare, mu, dividende;
	mpz_inits(Nsquare, mu, dividende, NULL);
	mpz_pow_ui(Nsquare, N, 2);
	mpz_powm(dividende, cipher, phi, Nsquare);
	mpz_sub_ui(dividende, dividende, 1);
	mpz_div(output, dividende, N);
	// mu is the inv of phi modulo N
	mpz_invert(mu, phi, N);
	mpz_mul(output, output, mu);
	mpz_mod(output, output, N);
	mpz_clears(Nsquare, mu, dividende, NULL);
}
/*
void paillierDecrypt(mpz_t cipher, mpz_t N, mpz_t output, mpz_t phi){
	mpz_t Nsquare, dividende;
	mpz_inits(Nsquare, dividende, NULL);
	mpz_pow_ui(Nsquare, N, 2);
	mpz_powm(dividende, cipher, phi, Nsquare);
	mpz_sub_ui(dividende, dividende, 1);
	
	
	// (N exp phi mod N**2 - 1) / N
	
	mpz_div(output, dividende, N);
	mpz_mod(output, output, N);
	mpz_clears(Nsquare, dividende, NULL);
}
*/
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
	gmp_printf("\nEncrypting %Zu gives %Zu\n", m2, c2);
	//paillierDecrypt(c2, N, m2, phi);
	//gmp_printf("Decrypting c2 %Zu  __ gives __ %Zu\n", c2, m2);
	mpz_t mul_e;
        mpz_init(mul_e);
	mpz_mul(mul_e, c, c2);
	mpz_pow_ui(N2, N, 2);
	mpz_mod(mul_e, mul_e, N2);
	gmp_printf("\nMUL c et c2: %Zu \n", mul_e);
	paillierDecrypt(mul_e, N, o, phi);
	gmp_printf("Decrypting %Zu __ gives __ %Zu\n", mul_e, o);
	
	mpz_add(m, m, m2);
	mpz_mul(r, r, r2);
	paillierEncrypt(m, r, N, c);
	gmp_printf("\n Encrypting from m1+m2 %Zu, and r1+r2 %Zu gives %Zu\n", m,r, c);
	paillierDecrypt(c, N, m, phi);
	gmp_printf("Decrypting %Zu gives %Zu\n", c, m);
        mpz_clears(N, m, c, r, p, q, phi, o, NULL);
	mpz_clears(m2, m3, c2, c3, r2, r3, N2, NULL);
	mpz_clear(mul_e);
	gmp_randclear(state);
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
		//gmp_printf("random_bytes n°%d: %Zu\n", i, random_bytes);
		i++;
	}
	while (i<V);

	mpz_clears(random_bytes, NULL);
	gmp_randclear(state);
}


/****************************************************************************
 *
 * Pre-watermarking 
 *
 * (Create an array of p elements, which are random bits)
 *
 * (Embed the watermark on the data)
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
		//gmp_printf("random_bit n°%d: %Zu\n",j, random_bit);
		j++;
	}
	while (j<p);

	mpz_clears(random_bit, NULL);
	gmp_randclear(state);
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

			// get the parity of the pixel, 0 if odd(LSB=1), 1 if even(LSB=0)
			even = mpz_divisible_p(pixel,two);

			// embed the watermark on the LSB
			if ((mpz_cmp_ui(watermark[jp],1) != 0)&&(even==1)){ //Wi=0 and LSB=0 -> OK, nothing to do
			}
			else if ((mpz_cmp_ui(watermark[jp],1) != 0)&&(even==0)){ //Wi=0 and LSB=1 -> Insert Wi instead of the LSB --> substract one
				mpz_sub_ui(data[in*p+jp], data[in*p+jp], 1);
			}
			else if ((mpz_cmp_ui(watermark[jp],0) != 0)&&(even==1)){ //Wi=1 and LSB=0 -> Insert Wi instead of the LSB --> add one
				mpz_add_ui(data[in*p+jp], data[in*p+jp], 1);
			}
			else if ((mpz_cmp_ui(watermark[jp],0) != 0)&&(even==0)){ //Wi=1 and LSB=1 -> OK, nothing to do
			}
			else{
				printf("ERROR, THIS CASE SHOULD NOT HAPPENED !");
			}
		}
	}
	mpz_clears(pixel, two, NULL);
}


/****************************************************************************
 *
 * Encryption 
 *
 ****************************************************************************/

void paillier_init(mpz_t p1, mpz_t q1, mpz_t N1, mpz_t r, mpz_t phi){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	generatePrime(p1, SECURITY_BITS, state);
	generatePrime(q1, SECURITY_BITS, state);
	mpz_mul(N1, p1, q1);
	computePhi(p1, q1, phi);
	gmp_randclear(state);
}

void data_encryption(mpz_t * data, unsigned long int V, mpz_t * encrypted_data, mpz_t N1, mpz_t r, mpz_t phi){
	mpz_t encrypted_value, decrypted_value;
	mpz_inits(encrypted_value, decrypted_value, NULL);
	for (int i = 0; i<V; i++){
		paillierEncrypt(data[i], r, N1, encrypted_value);
		mpz_init_set(encrypted_data[i], encrypted_value);
		//paillierDecrypt(encrypted_value, N1, decrypted_value, phi); 
	}
	mpz_clears(encrypted_value, decrypted_value, NULL);
}


/****************************************************************************
 *
 * Message embedding 
 *
 ****************************************************************************/

void embed_message_to_pixel(mpz_t enc_emb_pixel, mpz_t enc_pixel, mpz_t embedding, mpz_t distortion, mpz_t N, gmp_randstate_t state){ 
	mpz_t r, enc_distortion, Nsquare;
	mpz_inits(r, enc_distortion, Nsquare, NULL);
	mpz_pow_ui(Nsquare, N, 2);
	int even;
	do {
		//Pick random r : 0<r<N
		do{
			mpz_urandomm(r, state, N);
		}
		while (mpz_cmp_ui(r,0)==0);
        	paillierEncrypt(distortion, r, N, enc_distortion);
		mpz_mul(enc_emb_pixel ,enc_pixel, enc_distortion);
		// The multiplication of Paillet encrypted value is modulo N squared.
	    	mpz_mod(enc_emb_pixel, enc_emb_pixel, Nsquare);
		even = mpz_divisible_ui_p(enc_emb_pixel,2);
	}
	while(mpz_cmp_ui(embedding,even) == 0);
        // when embeeding = 0, enc_emb_pixel must be even, so even must be 1;
        // when embeeding = 1, enc_emb_pixel must be odd, so even must be 0;
	mpz_clears(r, enc_distortion, Nsquare, NULL);
}

void message_embedding(mpz_t * enc_emb_data, mpz_t * enc_data, unsigned long int V, int p, int N, mpz_t N1, mpz_t * watermark){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	mpz_t enc_pixel, enc_emb_pixel, embedding_p, distortion_p;
	mpz_inits(enc_pixel, enc_emb_pixel, embedding_p, distortion_p, NULL);
	for (int in = 0; in<N; in++){
        for (int jp = 0; jp<p; jp++){
            mpz_set(enc_pixel,enc_data[in*p+jp]);
			mpz_set(embedding_p, watermark[in]);
			// distortion = embedding * delta / 2;
			mpz_set(distortion_p, watermark[in]);
			embed_message_to_pixel(enc_emb_pixel, enc_pixel, embedding_p, distortion_p, N1, state);
			mpz_init_set(enc_emb_data[in*p+jp], enc_emb_pixel);
		}
	}
	mpz_clears(enc_pixel, enc_emb_pixel, embedding_p, distortion_p, NULL);
	gmp_randclear(state);
}


/****************************************************************************
 *
 * Message extraction from encrypted domain
 *
 ****************************************************************************/

void message_extraction_enc(mpz_t * enc_emb_data, unsigned long int V, mpz_t * extracted_message, unsigned int p, unsigned int N){
	// TO DO: redundancy -> check for each subset that the message is the same
	mpz_t enc_pixel, two;
	mpz_inits(enc_pixel, two, NULL);
	mpz_set_ui(two, 2);
	unsigned int even;
	for (int in = 0; in<N; in++){ // define the message values
		mpz_set(enc_pixel,enc_emb_data[in*p]);
		// get the parity of the pixel, 0 if odd(LSB=1), 1 if even(LSB=0)
		even = mpz_divisible_p(enc_pixel,two);

		// if the encrypted pixel is even, the message is 0
		if (even==1){ mpz_init_set_ui(extracted_message[in], 0); }

		// if the encrypted pixel is odd, the message is 1
		else { mpz_init_set_ui(extracted_message[in], 1); }

		for (int jp = 1; jp<p; jp++){ // redundancy, check for each subset that the message value defined is correct
			mpz_set(enc_pixel,enc_emb_data[in*p+jp]);
			even = mpz_divisible_p(enc_pixel,two);
			if ((even==1) && (mpz_cmp_ui(extracted_message[in],0) == 0)){
				// redundancy: the message value i_n corresponds
			}
			// if the encrypted pixel is odd, the message is 1
			else if ((even==0) && (mpz_cmp_ui(extracted_message[in],1) == 0)) {
				// redundancy: the message value i_n corresponds
			}
			else {
				printf("PROBLEM !\n");
			}
		}
	}
	mpz_clears(enc_pixel, two, NULL);
}


/****************************************************************************
 *
 * Data decryption 
 *
 ****************************************************************************/

void data_decryption(mpz_t * encrypted_data, unsigned long int V, mpz_t * decrypted_data, mpz_t N1, mpz_t phi){
	mpz_t decrypted_value, modulo;
	mpz_inits(decrypted_value, modulo, NULL);
	mpz_set_ui(modulo, 256);
	for (int i = 0; i<V; i++){
		paillierDecrypt(encrypted_data[i], N1, decrypted_value, phi); 
		// If pixel decrypted is out of bound, reduce by two its value
		if (mpz_cmp_ui(decrypted_value, 256) == 0) {
			mpz_set_ui(decrypted_value, 254);
		}
		mpz_init_set(decrypted_data[i], decrypted_value);
	}
	mpz_clears(decrypted_value, modulo, NULL);

}


/****************************************************************************
 *
 * Message extraction from spatial domain
 *
 ****************************************************************************/

void message_extraction_spa(mpz_t * decrypted_data, unsigned long int V, mpz_t * extracted_message, mpz_t * watermark, unsigned int p, unsigned int N){
	// TO DO: redundancy -> check for each subset that the message is the same
	mpz_t spa_pixel, two;
	mpz_inits(spa_pixel, two, NULL);
	mpz_set_ui(two,2);
	unsigned int even;
	for (int in = 0; in<N; in++){
		mpz_set(spa_pixel,decrypted_data[in*p]);
		even = mpz_divisible_p(spa_pixel, two);

		// if the pixel in the spa. domain is the same than the one in the prewatermark -> the message is 0
		// get the parity of the pixel, 0 if odd(LSB=1), 1 if even(LSB=0)
		if ((mpz_cmp_ui(watermark[0],1) == 0)&&(even==1)){ //Wi=1 and LSB=0 -> the message is 1 
			mpz_init_set_ui(extracted_message[in], 1);
		}
		else if ((mpz_cmp_ui(watermark[0],1) == 0)&&(even==0)){ //Wi=1 and LSB=1 -> the message is 0 
			mpz_init_set_ui(extracted_message[in], 0);
		}
		else if ((mpz_cmp_ui(watermark[0],0) == 0)&&(even==1)){ //Wi=0 and LSB=0 -> the message is 0
			mpz_init_set_ui(extracted_message[in], 0);
		}
		else if ((mpz_cmp_ui(watermark[0],0) == 0)&&(even==0)){ //Wi=0 and LSB=1 -> the message is 1 
			mpz_init_set_ui(extracted_message[in], 1);
		}
		else{
			printf("ERROR, THIS CASE SHOULD NOT HAPPENED !");
		}
		for (int jp = 1; jp<p; jp++){
			mpz_set(spa_pixel,decrypted_data[in*p+jp]);
			even = mpz_divisible_p(spa_pixel, two);

			// if the pixel in the spa. domain is the same than the one in the prewatermark -> the message is 0
			// get the parity of the pixel, 0 if odd(LSB=1), 1 if even(LSB=0)
			if ((mpz_cmp_ui(watermark[jp],1) == 0)&&(even==1)){ //Wi=1 and LSB=0 -> the message is 1 
				mpz_init_set_ui(extracted_message[in], 1);
			}
			else if ((mpz_cmp_ui(watermark[jp],1) == 0)&&(even==0)){ //Wi=1 and LSB=1 -> the message is 0 
				mpz_init_set_ui(extracted_message[in], 0);
			}
			else if ((mpz_cmp_ui(watermark[jp],0) == 0)&&(even==1)){ //Wi=0 and LSB=0 -> the message is 0
				mpz_init_set_ui(extracted_message[in], 0);
			}
			else if ((mpz_cmp_ui(watermark[jp],0) == 0)&&(even==0)){ //Wi=0 and LSB=1 -> the message is 1 
				mpz_init_set_ui(extracted_message[in], 1);
			}
			else{
				printf("ERROR, THIS CASE SHOULD NOT HAPPENED !");
			}
		}
	}
	mpz_clears(spa_pixel, two, NULL);
}


/****************************************************************************
 *
 * DEBUG: display data 
 *
 ****************************************************************************/

void display_array(mpz_t *array, unsigned long int V){
	for (int i=0; i<V;i++){
		gmp_printf("array[%d]=%Zu\n",i,array[i]);
	}
	printf("\n");
}


/****************************************************************************
 *
 * Main 
 *
 ****************************************************************************/

int main(int argc, char* argv[]) {
	// Init parameters
	unsigned int p = 20;
	unsigned int N = 40;
	unsigned long int V=p*N;
	

	// Data Generation
	mpz_t * data;
	data = malloc(sizeof(mpz_t) * V);
	data_generation(data, V);


	// Pre-watermarking
	mpz_t * watermark;
	watermark = malloc(sizeof(mpz_t) * p);
	watermark_generation(watermark, p);

	printf("##########\nOriginal data:\n");
	display_array(data, V);
	printf("##########\nWatermark:\n");
	display_array(watermark, p);

	pre_watermarking(watermark, data, N, p);

	printf("##########\nPrewatermarked data:\n");
	display_array(data, V);
	

	// Encryption
	mpz_t * encrypted_data;
	encrypted_data = malloc(sizeof(mpz_t) * V);

	// Paillier parameters
	mpz_t N1, r, p1, q1, phi, o;
	mpz_inits(N1, r, p1, q1, phi, o, NULL);
	printf("r: ");
	gmp_scanf("%Zu",r);
	paillier_init(p1,q1,N1,r,phi);

		// Paillier encryption
	data_encryption(data, V, encrypted_data, N1, r, phi);
	//printf("##########\nEncrypted data:\n");
	//display_array(encrypted_data, V);


	// Message embedding
	mpz_t * enc_emb_data;
	enc_emb_data = malloc(sizeof(mpz_t) * V);
	mpz_t * message;
	message = malloc(sizeof(mpz_t) * N);
	watermark_generation(message, N);
	printf("##########\nMessage:\n");
	display_array(message, N);
	message_embedding(enc_emb_data, encrypted_data, V, p, N, N1, message);
	//printf("##########\nEncrypted emb. data:\n");
	//display_array(enc_emb_data, V);
	

	// Mesage extraction from encrypted domain
	mpz_t * extracted_message;
	extracted_message = malloc(sizeof(mpz_t) * N);
	message_extraction_enc(enc_emb_data, V, extracted_message, p, N);
	printf("##########\nExtracted message from encrypted domain:\n");
	display_array(extracted_message, N);
	

	// Data decryption 
	mpz_t * decrypted_data;
	decrypted_data = malloc(sizeof(mpz_t) * V);
	//data_decryption(encrypted_data, V, decrypted_data, N1, phi);
	data_decryption(enc_emb_data, V, decrypted_data, N1, phi);
	printf("##########\nDecrypted data:\n");
	display_array(decrypted_data, V);


	// Mesage extraction from spatial domain
	message_extraction_spa(decrypted_data, V, extracted_message, watermark, p, N);
	printf("##########\nExtracted message from spatial domain:\n");
	display_array(extracted_message, N);
	

	// clear/free
	mpz_clears(N1, r, p1, q1, phi, o, NULL);
	free(data);
	free(watermark);
	free(encrypted_data);
	free(message);
	free(extracted_message);
	return 0;
}
