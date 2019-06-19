typedef union {
	uint4 cl_vec;
	uint arr[4];
} md4_hash;

inline void index_to_plaintext(unsigned long index, unsigned char *plaintext) {
  __attribute__((xcl_pipeline_loop))
  for (int i = 7; i >= 0; i--) {
    // TODO: should verify that the compiler. If not, need to do division and
    // modulo by ourselves to optimize the circuit.
    //
    // Regardless, it may be more efficient to do the division once and then
    // calculate the remainder using multiplication and subtraction on the
    // quotient and dividend. Using the modulo operator followed by the
    // division operator, even with the same constant, likely causes the
    // synthesis tool to generate two divisions.
    unsigned long quotient = index / 95;

    // ASSUMPTION: charset is ASCII-95. In that case, instead of using a lookup
    // table we can translate the index to a printable character by adding 32,
    // or even better by ORing 0x20.
    //
    // remainder = n - (q * c), c = 95
    plaintext[i] = (index - quotient * 95) | 0x20;
    index = quotient;
  }
}


/*
 * MD4 OpenCL kernel based on Solar Designer's MD4 algorithm implementation at:
 * http://openwall.info/wiki/people/solar/software/public-domain-source-code/md4
 * This code is in public domain.
 *
 * This software is Copyright (c) 2010, Dhiru Kholia <dhiru.kholia at gmail.com>
 * and Copyright (c) 2012, magnum
 * and Copyright (c) 2015, Sayantan Datta <std2048@gmail.com>
 * and it is hereby released to the general public under the following terms:
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted.
 *
 * Useful References:
 * 1  nt_opencl_kernel.c (written by Alain Espinosa <alainesp at gmail.com>)
 * 2. http://tools.ietf.org/html/rfc1320
 * 3. http://en.wikipedia.org/wiki/MD4
 */

#undef MD4_LUT3 /* No good for this format, just here for reference */
#define HAVE_ANDNOT 1

/* The basic MD4 functions */
#if MD4_LUT3
#define F(x, y, z)	lut3(x, y, z, 0xca)
#elif USE_BITSELECT
#define F(x, y, z)	bitselect((z), (y), (x))
#elif HAVE_ANDNOT
#define F(x, y, z)	((x & y) ^ ((~x) & z))
#else
#define F(x, y, z)	(z ^ (x & (y ^ z)))
#endif

#if MD4_LUT3
#define G(x, y, z)	lut3(x, y, z, 0xe8)
#else
#define G(x, y, z)	(((x) & ((y) | (z))) | ((y) & (z)))
#endif

#if MD4_LUT3
#define H(x, y, z)	lut3(x, y, z, 0x96)
#define H2 H
#else
#define H(x, y, z)	(((x) ^ (y)) ^ (z))
#define H2(x, y, z)	((x) ^ ((y) ^ (z)))
#endif

/* The MD4 transformation for all three rounds. */
#define STEP(f, a, b, c, d, x, s)	  \
	(a) += f((b), (c), (d)) + (x); \
	(a) = rotate((a), (uint)(s)) //(a) = ((a << s) | (a >> (32 - s))) 

void md4_round1(__private md4_hash *hash, __private uint W, uchar s)
{
	hash->arr[0] = F(hash->arr[1], hash->arr[2], hash->arr[3]) + W;
	hash->cl_vec = rotate(hash->cl_vec, s);
}

void md4_round2(__private md4_hash *hash, __private uint W, uchar s)
{
	hash->arr[0] = G(hash->arr[1], hash->arr[2], hash->arr[3]) + W + 0x5a827999;
	hash->cl_vec = rotate(hash->cl_vec, s);
}

void md4_round3a(__private md4_hash *hash, __private uint W, uchar s)
{
	hash->arr[0] = H(hash->arr[1], hash->arr[2], hash->arr[3]) + W + 0x6ed9eba1;
	hash->cl_vec = rotate(hash->cl_vec, s);
}

void md4_round3b(__private md4_hash *hash, __private uint W, uchar s)
{
	hash->arr[0] = H2(hash->arr[1], hash->arr[2], hash->arr[3]) + W + 0x6ed9eba1;
	hash->cl_vec = rotate(hash->cl_vec, s);
}

// Since we're treating each MD4 round as its own task, instruct the
// synthesis tool to pipeline each call.
__attribute__((xcl_pipeline_workitems))
void md4_encrypt(__private md4_hash *hash, __private uint *W)
{
	hash->arr[0] = 0x67452301;
	hash->arr[1] = 0xefcdab89;
	hash->arr[2] = 0x98badcfe;
	hash->arr[3] = 0x10325476;

	/* Round 1 */
	md4_round1(hash, W[0], 3);
	md4_round1(hash, W[1], 7);
	md4_round1(hash, W[2], 11);
	md4_round1(hash, W[3], 19);
	md4_round1(hash, W[4], 3);
	md4_round1(hash, W[5], 7);
	md4_round1(hash, W[6], 11);
	md4_round1(hash, W[7], 19);
	md4_round1(hash, W[8], 3);
	md4_round1(hash, W[9], 7);
	md4_round1(hash, W[10], 11);
	md4_round1(hash, W[11], 19);
	md4_round1(hash, W[12], 3);
	md4_round1(hash, W[13], 7);
	md4_round1(hash, W[14], 11);
	md4_round1(hash, W[15], 19);

	/* Round 2 */
	md4_round2(hash, W[0], 3);
	md4_round2(hash, W[4], 5);
	md4_round2(hash, W[8], 9);
	md4_round2(hash, W[12], 13);
	md4_round2(hash, W[1], 3);
	md4_round2(hash, W[5], 5);
	md4_round2(hash, W[9], 9);
	md4_round2(hash, W[13], 13);
	md4_round2(hash, W[2], 3);
	md4_round2(hash, W[6], 5);
	md4_round2(hash, W[10], 9);
	md4_round2(hash, W[14], 13);
	md4_round2(hash, W[3], 3);
	md4_round2(hash, W[7], 5);
	md4_round2(hash, W[11], 9);
	md4_round2(hash, W[15], 13);

	/* Round 3 */
	md4_round3a(hash, W[0], 3);
	md4_round3b(hash, W[8], 9);
	md4_round3a(hash, W[4], 11);
	md4_round3b(hash, W[12], 15);
	md4_round3a(hash, W[2], 3);
	md4_round3b(hash, W[10], 9);
	md4_round3a(hash, W[6], 11);
	md4_round3b(hash, W[14], 15);
	md4_round3a(hash, W[1], 3);
	md4_round3b(hash, W[9], 9);
	md4_round3a(hash, W[5], 11);
	md4_round3b(hash, W[13], 15);
	md4_round3a(hash, W[3], 3);
	md4_round3b(hash, W[11], 9);
	md4_round3a(hash, W[7], 11);
	md4_round3b(hash, W[15], 15);

	hash->arr[0] += 0x67452301;
	hash->arr[1] += 0xefcdab89;
	hash->arr[2] += 0x98badcfe;
	hash->arr[3] += 0x10325476;
}

unsigned long ntlm_hash(unsigned char *plaintext, unsigned char *hash, unsigned int pos) {
  unsigned int key[16] = {0};
  unsigned int output[4];

  // Unrolling this loop enables the FPGA synthesizer to perform these
  // assignments in parallel
  /*for (int i = 0; i < 4; i++) {
    key[i] = plaintext[i * 2] | (plaintext[(i * 2) + 1] << 16);
  }*/

  key[0] = plaintext[0] | (plaintext[1] << 16);
  key[1] = plaintext[2] | (plaintext[3] << 16);
  key[2] = plaintext[4] | (plaintext[5] << 16);
  key[3] = plaintext[6] | (plaintext[7] << 16);
  key[4] = 0x80;
  key[14] = 0x80;

  md4_encrypt((md4_hash*)output, key);

  unsigned long ret = ((unsigned long)output[1]) << 32 | (unsigned long)output[0];
  return (ret + pos) % 6634204312890625UL;
}


/* TODO: specify array length in definition...somehow? */
__kernel void crackalack_fpga_ntlm8(__global unsigned long *g_start_indices, __global unsigned long *g_end_indices) {

  unsigned long index = g_start_indices[get_global_id(0)];
  unsigned char plaintext[8];
  unsigned char hash[8];


  for (unsigned int pos = 0; pos < 421999; pos++) {
    index_to_plaintext(index, plaintext);
    index = ntlm_hash(plaintext, hash, pos);
  }

  g_end_indices[get_global_id(0)] = index;
  return;
}
