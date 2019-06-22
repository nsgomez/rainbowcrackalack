// Modules for calculating remainder by multiplying by inverse of divisor
//
// Motivation: in index_to_plaintext, it'd be faster to unroll the loop and do
// modulus by 95^(i + 1). The lack of data dependencies allows the FPGA
// synthesizer to do all the assignments in parallel in one step. However, the
// modulus operator isn't optimized in XOCC, generating a suboptimal circuit
// that takes many more clock cycles to divide.
//
// We can optimize modulus by a constant using the multiplicative inverse,
// addition/subtraction and shifts.
// You can see https://www.hackersdelight.org/hdcodetxt/remuc.c.txt for
// an explanation, but the magic numbers in this section were derived by
// reversing x86 output on godbolt.org.
//
// Comparing the XOCC reports using this method as opposed to naive modulus
// using %, the hand-optimizations here use about 5x fewer flip-flops (FFs)
// and 1/3 the FPGA lookup tables (LUTs) while reducing the clock interval
// by about 25%. The reduced amount of sequential logic also lets us increase
// the clock frequency to about 450MHz as opposed to the default 250MHz, so
// the overall speed improvement is - theoretically - 2.25x.

#define DEFINE_DIV_A(divisor, magic, shifts) \
  __attribute__((always_inline)) \
  inline ulong div ## divisor(ulong x) { \
    ulong q = x * magic ## UL; \
    ulong hi = mul_hi(x, magic ## UL);\
    q = (x - hi) >> 1; \
    q = (q + hi) >> shifts; \
    return q; \
  }

#define DEFINE_DIV_B(divisor, magic, shifts) \
  __attribute__((always_inline)) \
  inline ulong div ## divisor(ulong x) { \
    return mul_hi(x, magic ## UL) >> shifts; \
  }

DEFINE_DIV_A(95, 0x58ED2308158ED231, 6);
DEFINE_DIV_B(9025, 0xE85F14E7CDD97ABD, 13);
DEFINE_DIV_A(857375, 0x391703EA2D9B97E7, 19);
DEFINE_DIV_B(81450625, 0xD2EC7934AAD9FFEF, 26);
DEFINE_DIV_A(7737809375, 0x1C3124A7F9102B07, 32);
DEFINE_DIV_B(735091890625, 0x17EE949A457A9FEB, 36);
DEFINE_DIV_A(69833729609375, 0x1F616A9FFC60F49, 45);
DEFINE_DIV_B(6634204312890625, 0x2B72345286B2DF8B, 50);

__attribute__((always_inline))
inline uchar mod95(ulong x) { return (x - (div95(x) * 95UL)) & 0xFF; }

//__attribute__((always_inline))
//__attribute__((xcl_pipeline_workitems))
//inline ulong mod6634204312890625(ulong x) { return x - (div6634204312890625(x) * 6634204312890625UL); }

__attribute__((always_inline))
inline void reduce_index(ulong index, ulong *reduced) {
  reduced[7] = index;
  reduced[6] = div95(index);
  reduced[5] = div9025(index);
  reduced[4] = div857375(index);
  reduced[3] = div81450625(index);
  reduced[2] = div7737809375(index);
  reduced[1] = div735091890625(index);
  reduced[0] = div69833729609375(index);
}

void reduced_indices_to_plaintext(const ulong *restrict reduced, uchar *restrict plaintext) {
  __attribute__((opencl_unroll_hint))
  for (uchar i = 0; i < 7; i++) {
    plaintext[i] = mod95(reduced[i]) + 32;
  }
}

__attribute__((xcl_pipeline_workitems))
void index_to_plaintext(ulong index, uchar *plaintext) {
  ulong reduced_indices[8];

  reduce_index(index, reduced_indices);
  reduced_indices_to_plaintext(reduced_indices, plaintext);
}

#if 0
void index_to_plaintext(ulong index, uchar *plaintext) {
  plaintext[7] = (uchar)32 + mod95(index);
  plaintext[6] = (uchar)32 + mod95(div95(index));
  plaintext[5] = (uchar)32 + mod95(div9025(index));
  plaintext[4] = (uchar)32 + mod95(div857375(index));
  plaintext[3] = (uchar)32 + mod95(div81450625(index));
  plaintext[2] = (uchar)32 + mod95(div7737809375(index));
  plaintext[1] = (uchar)32 + mod95(div735091890625(index));
  plaintext[0] = (uchar)32 + mod95(div69833729609375(index));

  // The code above is an expanded version of this loop. There are a couple of
  // tricks used to improve latency. First, we assume that the charset is
  // ASCII-95, so adding 32 to mod95(index) is the same as doing a lookup in an
  // ASCII-95 table but without the need to access memory. Second, iteratively
  // dividing `index` would introduce a data dependency in the loop, forcing
  // each assignment to `plaintext` to be done sequentially. We can each cell
  // of the plaintext output in parallel by unrolling the iterative divisions.
  // Since the assignments are now all a function of `index`, which doesn't
  // change, these can be done independently and the synthesizer parallelizes
  // the math.

  /*for (int i = 7; i >= 0; i--) {
    plaintext[i] = charset[index % 95];
    index = index / 95;
  }*/
}
#endif


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
#define F(x, y, z)	(((x) & (y)) ^ ((~(x)) & (z)))
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

// Using the OpenCL rotate() call usually causes the compiler and synthesizer
// to use a completely separate block to do the shift, which takes an extra
// clock cycle. An extra clock cycle in the MD4 module is a big bottleneck,
// since the MD4 loop is called 422,000 times and is completely sequential.

#define FAST_ROTATE(x, s) ((x) << (s)) | (((x) & (0xFFFFFFFFU << (32 - (s)))) >> (32 - (s)))

/* The MD4 transformation for all three rounds. */
#define STEP(f, a, b, c, d, x, s) \
  (a) = FAST_ROTATE((a) + f((b), (c), (d)) + (x), (s));


// NOTE: MD4 hashing in this kernel uses a sparse representation of the message
// block. Since NTLM-8 always uses 0x80 for key[4] and key[14] and all elements
// between key[5] and key[13] inclusive are zero, we only need four bytes for
// the key, which reduces the bandwidth needed for I/O in the MD4 component.

__attribute__((xcl_pipeline_workitems))
uint4 md4_encrypt(const uint4 W)
{
  uint4 hash;
  hash[0] = 0x67452301;
  hash[1] = 0xefcdab89;
  hash[2] = 0x98badcfe;
  hash[3] = 0x10325476;

  /* Round 1 */
  STEP(F, hash[0], hash[1], hash[2], hash[3], W[0], 3);
  STEP(F, hash[3], hash[0], hash[1], hash[2], W[1], 7);
  STEP(F, hash[2], hash[3], hash[0], hash[1], W[2], 11);
  STEP(F, hash[1], hash[2], hash[3], hash[0], W[3], 19);
  STEP(F, hash[0], hash[1], hash[2], hash[3], 0x80, 3);
  STEP(F, hash[3], hash[0], hash[1], hash[2], 0x00, 7);
  STEP(F, hash[2], hash[3], hash[0], hash[1], 0x00, 11);
  STEP(F, hash[1], hash[2], hash[3], hash[0], 0x00, 19);
  STEP(F, hash[0], hash[1], hash[2], hash[3], 0x00, 3);
  STEP(F, hash[3], hash[0], hash[1], hash[2], 0x00, 7);
  STEP(F, hash[2], hash[3], hash[0], hash[1], 0x00, 11);
  STEP(F, hash[1], hash[2], hash[3], hash[0], 0x00, 19);
  STEP(F, hash[0], hash[1], hash[2], hash[3], 0x00, 3);
  STEP(F, hash[3], hash[0], hash[1], hash[2], 0x00, 7);
  STEP(F, hash[2], hash[3], hash[0], hash[1], 0x80, 11);
  STEP(F, hash[1], hash[2], hash[3], hash[0], 0x00, 19);

  /* Round 2 */
  STEP(G, hash[0], hash[1], hash[2], hash[3], W[0] + 0x5a827999, 3);
  STEP(G, hash[3], hash[0], hash[1], hash[2], 0x80 + 0x5a827999, 5);
  STEP(G, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x5a827999, 9);
  STEP(G, hash[1], hash[2], hash[3], hash[0], 0x00 + 0x5a827999, 13);
  STEP(G, hash[0], hash[1], hash[2], hash[3], W[1] + 0x5a827999, 3);
  STEP(G, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x5a827999, 5);
  STEP(G, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x5a827999, 9);
  STEP(G, hash[1], hash[2], hash[3], hash[0], 0x00 + 0x5a827999, 13);
  STEP(G, hash[0], hash[1], hash[2], hash[3], W[2] + 0x5a827999, 3);
  STEP(G, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x5a827999, 5);
  STEP(G, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x5a827999, 9);
  STEP(G, hash[1], hash[2], hash[3], hash[0], 0x80 + 0x5a827999, 13);
  STEP(G, hash[0], hash[1], hash[2], hash[3], W[3] + 0x5a827999, 3);
  STEP(G, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x5a827999, 5);
  STEP(G, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x5a827999, 9);
  STEP(G, hash[1], hash[2], hash[3], hash[0], 0x00 + 0x5a827999, 13);

  /* Round 3 */
  STEP(H, hash[0], hash[1], hash[2], hash[3], W[0] + 0x6ed9eba1, 3);
  STEP(H2, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x6ed9eba1, 9);
  STEP(H, hash[2], hash[3], hash[0], hash[1], 0x80 + 0x6ed9eba1, 11);
  STEP(H2, hash[1], hash[2], hash[3], hash[0], 0x00 + 0x6ed9eba1, 15);
  STEP(H, hash[0], hash[1], hash[2], hash[3], W[2] + 0x6ed9eba1, 3);
  STEP(H2, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x6ed9eba1, 9);
  STEP(H, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x6ed9eba1, 11);
  STEP(H2, hash[1], hash[2], hash[3], hash[0], 0x80 + 0x6ed9eba1, 15);
  STEP(H, hash[0], hash[1], hash[2], hash[3], W[1] + 0x6ed9eba1, 3);
  STEP(H2, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x6ed9eba1, 9);
  STEP(H, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x6ed9eba1, 11);
  STEP(H2, hash[1], hash[2], hash[3], hash[0], 0x00 + 0x6ed9eba1, 15);
  STEP(H, hash[0], hash[1], hash[2], hash[3], W[3] + 0x6ed9eba1, 3);
  STEP(H2, hash[3], hash[0], hash[1], hash[2], 0x00 + 0x6ed9eba1, 9);
  STEP(H, hash[2], hash[3], hash[0], hash[1], 0x00 + 0x6ed9eba1, 11);
  STEP(H2, hash[1], hash[2], hash[3], hash[0], 0x00 + 0x6ed9eba1, 15);

  hash[0] = hash[0] + 0x67452301;
  hash[1] = hash[1] + 0xefcdab89;
  hash[2] = hash[2] + 0x98badcfe;
  hash[3] = hash[3] + 0x10325476;

  return hash;
}

__attribute__((xcl_pipeline_workitems))
inline ulong ntlm_hash(const uchar *plaintext) {
  uint4 key = 0;
  key[0] = plaintext[0] | (plaintext[1] << 16);
  key[1] = plaintext[2] | (plaintext[3] << 16);
  key[2] = plaintext[4] | (plaintext[5] << 16);
  key[3] = plaintext[6] | (plaintext[7] << 16);
  //key[4] = 0x80;
  //key[14] = 0x80;

  uint4 output = md4_encrypt(key);
  return ((ulong)output[1]) << 32 | (ulong)output[0];
}


/* TODO: specify array length in definition...somehow? */
__kernel
__attribute__((reqd_work_group_size(1, 1, 1)))
__attribute__((vec_type_hint(ulong)))
__attribute__((xcl_zero_global_work_offset))
void crackalack_fpga_ntlm8(__global const ulong *restrict g_start_indices, __global ulong *restrict g_end_indices)
{
  __private ulong ntlm, tmp;
  __private ulong index = g_start_indices[get_global_id(0)];
  __private uchar plaintext[8] __attribute__((xcl_array_partition(complete, 1)));

  __attribute__((xcl_pipeline_loop))
  for (uint pos = 0; pos < 421999; pos++) {
    index_to_plaintext(index, plaintext);
    ntlm = ntlm_hash(plaintext);
    tmp = div6634204312890625(ntlm);
    index = ntlm - (tmp * 6634204312890625UL);
  }

  g_end_indices[get_global_id(0)] = index;
  return;
}
