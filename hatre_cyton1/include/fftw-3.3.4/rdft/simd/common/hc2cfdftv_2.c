/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Tue Mar  4 13:51:49 EST 2014 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 2 -dit -name hc2cfdftv_2 -include hc2cfv.h */

/*
 * This function contains 5 FP additions, 6 FP multiplications,
 * (or, 3 additions, 4 multiplications, 2 fused multiply/add),
 * 9 stack variables, 1 constants, and 4 memory accesses
 */
#include "hc2cfv.h"

static void hc2cfdftv_2(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 2)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 2), MAKE_VOLATILE_STRIDE(8, rs)) {
	       V T1, T2, T4, T3, T5, T7, T6;
	       T1 = LD(&(Rp[0]), ms, &(Rp[0]));
	       T2 = LD(&(Rm[0]), -ms, &(Rm[0]));
	       T4 = LDW(&(W[0]));
	       T3 = VFMACONJ(T2, T1);
	       T5 = VZMULIJ(T4, VFNMSCONJ(T2, T1));
	       T7 = VCONJ(VMUL(LDK(KP500000000), VADD(T3, T5)));
	       T6 = VMUL(LDK(KP500000000), VSUB(T3, T5));
	       ST(&(Rm[0]), T7, -ms, &(Rm[0]));
	       ST(&(Rp[0]), T6, ms, &(Rp[0]));
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(1, 1),
     {TW_NEXT, VL, 0}
};

static const hc2c_desc desc = { 2, XSIMD_STRING("hc2cfdftv_2"), twinstr, &GENUS, {3, 4, 2, 0} };

void XSIMD(codelet_hc2cfdftv_2) (planner *p) {
     X(khc2c_register) (p, hc2cfdftv_2, &desc, HC2C_VIA_DFT);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 2 -dit -name hc2cfdftv_2 -include hc2cfv.h */

/*
 * This function contains 5 FP additions, 4 FP multiplications,
 * (or, 5 additions, 4 multiplications, 0 fused multiply/add),
 * 10 stack variables, 1 constants, and 4 memory accesses
 */
#include "hc2cfv.h"

static void hc2cfdftv_2(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 2)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 2), MAKE_VOLATILE_STRIDE(8, rs)) {
	       V T4, T6, T1, T3, T2, T5, T7, T8;
	       T1 = LD(&(Rp[0]), ms, &(Rp[0]));
	       T2 = LD(&(Rm[0]), -ms, &(Rm[0]));
	       T3 = VCONJ(T2);
	       T4 = VADD(T1, T3);
	       T5 = LDW(&(W[0]));
	       T6 = VZMULIJ(T5, VSUB(T3, T1));
	       T7 = VCONJ(VMUL(LDK(KP500000000), VSUB(T4, T6)));
	       ST(&(Rm[0]), T7, -ms, &(Rm[0]));
	       T8 = VMUL(LDK(KP500000000), VADD(T4, T6));
	       ST(&(Rp[0]), T8, ms, &(Rp[0]));
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(1, 1),
     {TW_NEXT, VL, 0}
};

static const hc2c_desc desc = { 2, XSIMD_STRING("hc2cfdftv_2"), twinstr, &GENUS, {5, 4, 0, 0} };

void XSIMD(codelet_hc2cfdftv_2) (planner *p) {
     X(khc2c_register) (p, hc2cfdftv_2, &desc, HC2C_VIA_DFT);
}
#endif				/* HAVE_FMA */
