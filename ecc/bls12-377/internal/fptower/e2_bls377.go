// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

package fptower

import (
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"
)

// Mul sets z to the E2-product of x,y, returns z
func (z *E2) Mul(x, y *E2) *E2 {
	var a, b, c fp.Element
	a.Add(&x.A0, &x.A1)
	b.Add(&y.A0, &y.A1)
	a.Mul(&a, &b)
	b.Mul(&x.A0, &y.A0)
	c.Mul(&x.A1, &y.A1)
	z.A1.Sub(&a, &b).Sub(&z.A1, &c)
	fp.MulBy5(&c)
	z.A0.Sub(&b, &c)
	return z
}

// Square sets z to the E2-product of x,x returns z
func (z *E2) Square(x *E2) *E2 {
	//algo 22 https://eprint.iacr.org/2010/354.pdf
	var c0, c2 fp.Element
	c0.Add(&x.A0, &x.A1)
	c2.Neg(&x.A1)
	fp.MulBy5(&c2)
	c2.Add(&c2, &x.A0)

	c0.Mul(&c0, &c2) // (x1+x2)*(x1+(u**2)x2)
	c2.Mul(&x.A0, &x.A1).Double(&c2)
	z.A1 = c2
	c2.Double(&c2)
	z.A0.Add(&c0, &c2)

	return z
}

// MulByNonResidue multiplies a E2 by (0,1)
func (z *E2) MulByNonResidue(x *E2) *E2 {
	a := x.A0
	b := x.A1 // fetching x.A1 in the function below is slower
	fp.MulBy5(&b)
	z.A0.Neg(&b)
	z.A1 = a
	return z
}

// MulByNonResidueInv multiplies a E2 by (0,1)^{-1}
func (z *E2) MulByNonResidueInv(x *E2) *E2 {
	//z.A1.MulByNonResidueInv(&x.A0)
	a := x.A1
	fiveinv := fp.Element{
		330620507644336508,
		9878087358076053079,
		11461392860540703536,
		6973035786057818995,
		8846909097162646007,
		104838758629667239,
	}
	z.A1.Mul(&x.A0, &fiveinv).Neg(&z.A1)
	z.A0 = a
	return z
}

// Inverse sets z to the E2-inverse of x, returns z
func (z *E2) Inverse(x *E2) *E2 {
	// Algorithm 8 from https://eprint.iacr.org/2010/354.pdf
	//var a, b, t0, t1, tmp fp.Element
	var t0, t1, tmp fp.Element
	a := &x.A0 // creating the buffers a, b is faster than querying &x.A0, &x.A1 in the functions call below
	b := &x.A1
	t0.Square(a)
	t1.Square(b)
	tmp.Set(&t1)
	fp.MulBy5(&tmp)
	t0.Add(&t0, &tmp)
	t1.Inverse(&t0)
	z.A0.Mul(a, &t1)
	z.A1.Mul(b, &t1).Neg(&z.A1)

	return z
}

// norm sets x to the norm of z
func (z *E2) norm(x *fp.Element) {
	var tmp fp.Element
	x.Square(&z.A1)
	tmp.Set(x)
	fp.MulBy5(&tmp)
	x.Square(&z.A0).Add(x, &tmp)
}

// MulBybTwistCurveCoeff multiplies by 1/(0,1)
func (z *E2) MulBybTwistCurveCoeff(x *E2) *E2 {

	var res E2
	res.A0.Set(&x.A1)
	res.A1.MulByNonResidueInv(&x.A0)
	z.Set(&res)

	return z
}

// expByCbrtExp is equivalent to z.Exp(x, 2fcd9602a41e52f994c7c00c11ebb13bcafbc5aac5e5ba91a943cc6a010800028f7c2405555555641d6aaaaaaaaaab)
//
// uses github.com/mmcloughlin/addchain v0.4.0 to generate a shorter addition chain
func (z *E2) expByCbrtExp(x *E2) *E2 {
	// addition chain:
	//
	//	_10      = 2*1
	//	_11      = 1 + _10
	//	_101     = _10 + _11
	//	_111     = _10 + _101
	//	_1001    = _10 + _111
	//	_1011    = _10 + _1001
	//	_1101    = _10 + _1011
	//	_1111    = _10 + _1101
	//	_10001   = _10 + _1111
	//	_10011   = _10 + _10001
	//	_10101   = _10 + _10011
	//	_10111   = _10 + _10101
	//	_11001   = _10 + _10111
	//	_11011   = _10 + _11001
	//	_11101   = _10 + _11011
	//	_11111   = _10 + _11101
	//	_111110  = 2*_11111
	//	_111111  = 1 + _111110
	//	_1010000 = _10001 + _111111
	//	i37      = ((_1010000 << 6 + _10101) << 5 + _10001) << 5
	//	i48      = ((_10011 + i37) << 3 + _101) << 5 + _11
	//	i67      = ((i48 << 8 + _10011) << 5 + _10001) << 4
	//	i82      = ((_1011 + i67) << 6 + _10111) << 6 + _1111
	//	i100     = ((i82 << 6 + _11011) << 7 + _10101) << 3
	//	i116     = ((_111 + i100) << 8 + _1101) << 5 + _101
	//	i137     = ((i116 << 9 + _11101) << 6 + _11101) << 4
	//	i153     = ((_1101 + i137) << 5 + _1001) << 8 + _101
	//	i173     = ((i153 << 6 + _1101) << 5 + _1101) << 7
	//	i187     = ((_11101 + i173) << 4 + _101) << 7 + _1001
	//	i208     = ((i187 << 7 + _10011) << 6 + _101) << 6
	//	i226     = ((_1011 + i208) << 5 + _111) << 10 + _111111
	//	i244     = ((i226 << 2 + _11) << 4 + _111) << 10
	//	i259     = ((_11001 + i244) << 7 + _10101) << 5 + _1111
	//	i279     = ((i259 << 6 + _101) << 7 + _11001) << 5
	//	i292     = ((_11101 + i279) << 3 + 1) << 7 + _111111
	//	i320     = ((i292 << 7 + _1111) << 4 + _101) << 15
	//	i332     = ((_10101 + i320) << 4 + _1111) << 5 + _11
	//	i356     = ((i332 << 4 + 1) << 12 + _1011) << 6
	//	i369     = ((_10101 + i356) << 4 + _1011) << 6 + _11
	//	i389     = ((i369 << 6 + _111) << 6 + _1001) << 6
	//	i401     = ((_11011 + i389) << 2 + _11) << 7 + 1
	//	i422     = ((i401 << 10 + _10001) << 4 + _1101) << 5
	//	i440     = ((1 + i422) << 10 + _10111) << 5 + _11011
	//	i458     = ((i440 << 5 + _10011) << 5 + _10111) << 6
	//	i471     = ((_1001 + i458) << 5 + _1111) << 5 + _1101
	//	i490     = ((i471 << 8 + _1001) << 2 + 1) << 7
	//	i507     = ((_111111 + i490) << 7 + _111111) << 7 + _10011
	//	i528     = ((i507 << 7 + _111) << 7 + _10011) << 5
	//	i543     = ((_10111 + i528) << 7 + _10011) << 5 + _11001
	//	i562     = ((i543 << 4 + _1111) << 10 + _10011) << 3
	//	i582     = ((1 + i562) << 10 + _111111) << 7 + _11001
	//	i598     = ((i582 << 6 + _11011) << 5 + _11011) << 3
	//	i612     = ((_111 + i598) << 4 + 1) << 7 + _111
	//	i637     = ((i612 << 8 + _10101) << 11 + _11111) << 4
	//	i657     = ((_111 + i637) << 10 + _111111) << 7 + _1111
	//	i676     = ((i657 << 5 + _1001) << 6 + _11001) << 6
	//	i690     = ((_11111 + i676) << 5 + _1111) << 6 + _1011
	//	i707     = ((i690 << 6 + _11001) << 6 + _10011) << 3
	//	i724     = ((_111 + i707) << 9 + _11) << 5 + _11
	//	i756     = ((i724 << 20 + _11001) << 5 + _10011) << 5
	//	i770     = ((_1011 + i756) << 6 + _11001) << 5 + _1101
	//	i793     = ((i770 << 9 + _1101) << 6 + _10101) << 6
	//	i808     = ((_10101 + i793) << 6 + _10101) << 6 + _10101
	//	i831     = ((i808 << 6 + _11001) << 10 + _11101) << 5
	//	i846     = ((_1101 + i831) << 6 + _10101) << 6 + _10101
	//	i866     = ((i846 << 6 + _10101) << 6 + _10101) << 6
	//	i881     = ((_10101 + i866) << 6 + _10101) << 6 + _10101
	//	return     i881 << 5 + _1011
	//
	// Operations: 749 squares 138 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var (
		t0  = new(E2)
		t1  = new(E2)
		t2  = new(E2)
		t3  = new(E2)
		t4  = new(E2)
		t5  = new(E2)
		t6  = new(E2)
		t7  = new(E2)
		t8  = new(E2)
		t9  = new(E2)
		t10 = new(E2)
		t11 = new(E2)
		t12 = new(E2)
		t13 = new(E2)
		t14 = new(E2)
		t15 = new(E2)
	)

	// Step 1: t8 = x^0x2
	t8.Square(x)

	// Step 2: t5 = x^0x3
	t5.Mul(x, t8)

	// Step 3: t14 = x^0x5
	t14.Mul(t8, t5)

	// Step 4: t6 = x^0x7
	t6.Mul(t8, t14)

	// Step 5: t9 = x^0x9
	t9.Mul(t8, t6)

	// Step 6: z = x^0xb
	z.Mul(t8, t9)

	// Step 7: t1 = x^0xd
	t1.Mul(t8, z)

	// Step 8: t7 = x^0xf
	t7.Mul(t8, t1)

	// Step 9: t13 = x^0x11
	t13.Mul(t8, t7)

	// Step 10: t4 = x^0x13
	t4.Mul(t8, t13)

	// Step 11: t0 = x^0x15
	t0.Mul(t8, t4)

	// Step 12: t12 = x^0x17
	t12.Mul(t8, t0)

	// Step 13: t3 = x^0x19
	t3.Mul(t8, t12)

	// Step 14: t11 = x^0x1b
	t11.Mul(t8, t3)

	// Step 15: t2 = x^0x1d
	t2.Mul(t8, t11)

	// Step 16: t8 = x^0x1f
	t8.Mul(t8, t2)

	// Step 17: t10 = x^0x3e
	t10.Square(t8)

	// Step 18: t10 = x^0x3f
	t10.Mul(x, t10)

	// Step 19: t15 = x^0x50
	t15.Mul(t13, t10)

	// Step 25: t15 = x^0x1400
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 26: t15 = x^0x1415
	t15.Mul(t0, t15)

	// Step 31: t15 = x^0x282a0
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 32: t15 = x^0x282b1
	t15.Mul(t13, t15)

	// Step 37: t15 = x^0x505620
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 38: t15 = x^0x505633
	t15.Mul(t4, t15)

	// Step 41: t15 = x^0x282b198
	for s := 0; s < 3; s++ {
		t15.Square(t15)
	}

	// Step 42: t15 = x^0x282b19d
	t15.Mul(t14, t15)

	// Step 47: t15 = x^0x505633a0
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 48: t15 = x^0x505633a3
	t15.Mul(t5, t15)

	// Step 56: t15 = x^0x505633a300
	for s := 0; s < 8; s++ {
		t15.Square(t15)
	}

	// Step 57: t15 = x^0x505633a313
	t15.Mul(t4, t15)

	// Step 62: t15 = x^0xa0ac6746260
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 63: t15 = x^0xa0ac6746271
	t15.Mul(t13, t15)

	// Step 67: t15 = x^0xa0ac67462710
	for s := 0; s < 4; s++ {
		t15.Square(t15)
	}

	// Step 68: t15 = x^0xa0ac6746271b
	t15.Mul(z, t15)

	// Step 74: t15 = x^0x282b19d189c6c0
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 75: t15 = x^0x282b19d189c6d7
	t15.Mul(t12, t15)

	// Step 81: t15 = x^0xa0ac6746271b5c0
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 82: t15 = x^0xa0ac6746271b5cf
	t15.Mul(t7, t15)

	// Step 88: t15 = x^0x282b19d189c6d73c0
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 89: t15 = x^0x282b19d189c6d73db
	t15.Mul(t11, t15)

	// Step 96: t15 = x^0x14158ce8c4e36b9ed80
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 97: t15 = x^0x14158ce8c4e36b9ed95
	t15.Mul(t0, t15)

	// Step 100: t15 = x^0xa0ac6746271b5cf6ca8
	for s := 0; s < 3; s++ {
		t15.Square(t15)
	}

	// Step 101: t15 = x^0xa0ac6746271b5cf6caf
	t15.Mul(t6, t15)

	// Step 109: t15 = x^0xa0ac6746271b5cf6caf00
	for s := 0; s < 8; s++ {
		t15.Square(t15)
	}

	// Step 110: t15 = x^0xa0ac6746271b5cf6caf0d
	t15.Mul(t1, t15)

	// Step 115: t15 = x^0x14158ce8c4e36b9ed95e1a0
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 116: t15 = x^0x14158ce8c4e36b9ed95e1a5
	t15.Mul(t14, t15)

	// Step 125: t15 = x^0x282b19d189c6d73db2bc34a00
	for s := 0; s < 9; s++ {
		t15.Square(t15)
	}

	// Step 126: t15 = x^0x282b19d189c6d73db2bc34a1d
	t15.Mul(t2, t15)

	// Step 132: t15 = x^0xa0ac6746271b5cf6caf0d28740
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 133: t15 = x^0xa0ac6746271b5cf6caf0d2875d
	t15.Mul(t2, t15)

	// Step 137: t15 = x^0xa0ac6746271b5cf6caf0d2875d0
	for s := 0; s < 4; s++ {
		t15.Square(t15)
	}

	// Step 138: t15 = x^0xa0ac6746271b5cf6caf0d2875dd
	t15.Mul(t1, t15)

	// Step 143: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba0
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 144: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba9
	t15.Mul(t9, t15)

	// Step 152: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba900
	for s := 0; s < 8; s++ {
		t15.Square(t15)
	}

	// Step 153: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba905
	t15.Mul(t14, t15)

	// Step 159: t15 = x^0x505633a3138dae7b65786943aeea4140
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 160: t15 = x^0x505633a3138dae7b65786943aeea414d
	t15.Mul(t1, t15)

	// Step 165: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829a0
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 166: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad
	t15.Mul(t1, t15)

	// Step 173: t15 = x^0x505633a3138dae7b65786943aeea414d680
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 174: t15 = x^0x505633a3138dae7b65786943aeea414d69d
	t15.Mul(t2, t15)

	// Step 178: t15 = x^0x505633a3138dae7b65786943aeea414d69d0
	for s := 0; s < 4; s++ {
		t15.Square(t15)
	}

	// Step 179: t15 = x^0x505633a3138dae7b65786943aeea414d69d5
	t15.Mul(t14, t15)

	// Step 186: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea80
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 187: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea89
	t15.Mul(t9, t15)

	// Step 194: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a754480
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 195: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a754493
	t15.Mul(t4, t15)

	// Step 201: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c0
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 202: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c5
	t15.Mul(t14, t15)

	// Step 208: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a754493140
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 209: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b
	t15.Mul(z, t15)

	// Step 214: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea89262960
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 215: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea89262967
	t15.Mul(t6, t15)

	// Step 225: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c00
	for s := 0; s < 10; s++ {
		t15.Square(t15)
	}

	// Step 226: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3f
	t15.Mul(t10, t15)

	// Step 228: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670fc
	for s := 0; s < 2; s++ {
		t15.Square(t15)
	}

	// Step 229: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff
	t15.Mul(t5, t15)

	// Step 233: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff0
	for s := 0; s < 4; s++ {
		t15.Square(t15)
	}

	// Step 234: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7
	t15.Mul(t6, t15)

	// Step 244: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc00
	for s := 0; s < 10; s++ {
		t15.Square(t15)
	}

	// Step 245: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc19
	t15.Mul(t3, t15)

	// Step 252: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c80
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 253: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c95
	t15.Mul(t0, t15)

	// Step 258: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192a0
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 259: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af
	t15.Mul(t7, t15)

	// Step 265: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc0
	for s := 0; s < 6; s++ {
		t15.Square(t15)
	}

	// Step 266: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc5
	t15.Mul(t14, t15)

	// Step 273: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e280
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 274: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299
	t15.Mul(t3, t15)

	// Step 279: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc5320
	for s := 0; s < 5; s++ {
		t15.Square(t15)
	}

	// Step 280: t15 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d
	t15.Mul(t2, t15)

	// Step 283: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e8
	for s := 0; s < 3; s++ {
		t15.Square(t15)
	}

	// Step 284: t15 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e9
	t15.Mul(x, t15)

	// Step 291: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf480
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 292: t15 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf
	t15.Mul(t10, t15)

	// Step 299: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f80
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 300: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f
	t15.Mul(t7, t15)

	// Step 304: t15 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f0
	for s := 0; s < 4; s++ {
		t15.Square(t15)
	}

	// Step 305: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5
	t14.Mul(t14, t15)

	// Step 320: t14 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8000
	for s := 0; s < 15; s++ {
		t14.Square(t14)
	}

	// Step 321: t14 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015
	t14.Mul(t0, t14)

	// Step 325: t14 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a80150
	for s := 0; s < 4; s++ {
		t14.Square(t14)
	}

	// Step 326: t14 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f
	t14.Mul(t7, t14)

	// Step 331: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be0
	for s := 0; s < 5; s++ {
		t14.Square(t14)
	}

	// Step 332: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3
	t14.Mul(t5, t14)

	// Step 336: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be30
	for s := 0; s < 4; s++ {
		t14.Square(t14)
	}

	// Step 337: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be31
	t14.Mul(x, t14)

	// Step 349: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be31000
	for s := 0; s < 12; s++ {
		t14.Square(t14)
	}

	// Step 350: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b
	t14.Mul(z, t14)

	// Step 356: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402c0
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 357: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5
	t14.Mul(t0, t14)

	// Step 361: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d50
	for s := 0; s < 4; s++ {
		t14.Square(t14)
	}

	// Step 362: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b
	t14.Mul(z, t14)

	// Step 368: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c0
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 369: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c3
	t14.Mul(t5, t14)

	// Step 375: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c0
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 376: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c7
	t14.Mul(t6, t14)

	// Step 382: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c0
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 383: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c9
	t14.Mul(t9, t14)

	// Step 389: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c7240
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 390: t14 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725b
	t14.Mul(t11, t14)

	// Step 392: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96c
	for s := 0; s < 2; s++ {
		t14.Square(t14)
	}

	// Step 393: t14 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f
	t14.Mul(t5, t14)

	// Step 400: t14 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b780
	for s := 0; s < 7; s++ {
		t14.Square(t14)
	}

	// Step 401: t14 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781
	t14.Mul(x, t14)

	// Step 411: t14 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0400
	for s := 0; s < 10; s++ {
		t14.Square(t14)
	}

	// Step 412: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411
	t13.Mul(t13, t14)

	// Step 416: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de04110
	for s := 0; s < 4; s++ {
		t13.Square(t13)
	}

	// Step 417: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d
	t13.Mul(t1, t13)

	// Step 422: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a0
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 423: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a1
	t13.Mul(x, t13)

	// Step 433: t13 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8400
	for s := 0; s < 10; s++ {
		t13.Square(t13)
	}

	// Step 434: t13 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417
	t13.Mul(t12, t13)

	// Step 439: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082e0
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 440: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb
	t13.Mul(t11, t13)

	// Step 445: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f60
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 446: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73
	t13.Mul(t4, t13)

	// Step 451: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee60
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 452: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee77
	t13.Mul(t12, t13)

	// Step 458: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc0
	for s := 0; s < 6; s++ {
		t13.Square(t13)
	}

	// Step 459: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc9
	t13.Mul(t9, t13)

	// Step 464: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b920
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 465: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f
	t13.Mul(t7, t13)

	// Step 470: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725e0
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 471: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed
	t13.Mul(t1, t13)

	// Step 479: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed00
	for s := 0; s < 8; s++ {
		t13.Square(t13)
	}

	// Step 480: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed09
	t13.Mul(t9, t13)

	// Step 482: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b424
	for s := 0; s < 2; s++ {
		t13.Square(t13)
	}

	// Step 483: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b425
	t13.Mul(x, t13)

	// Step 490: t13 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda1280
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 491: t13 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf
	t13.Mul(t10, t13)

	// Step 498: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095f80
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 499: t13 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf
	t13.Mul(t10, t13)

	// Step 506: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf80
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 507: t13 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf93
	t13.Mul(t4, t13)

	// Step 514: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc980
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 515: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987
	t13.Mul(t6, t13)

	// Step 522: t13 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c380
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 523: t13 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393
	t13.Mul(t4, t13)

	// Step 528: t13 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987260
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 529: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277
	t12.Mul(t12, t13)

	// Step 536: t12 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b80
	for s := 0; s < 7; s++ {
		t12.Square(t12)
	}

	// Step 537: t12 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93
	t12.Mul(t4, t12)

	// Step 542: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277260
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 543: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279
	t12.Mul(t3, t12)

	// Step 547: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc9872772790
	for s := 0; s < 4; s++ {
		t12.Square(t12)
	}

	// Step 548: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f
	t12.Mul(t7, t12)

	// Step 558: t12 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c00
	for s := 0; s < 10; s++ {
		t12.Square(t12)
	}

	// Step 559: t12 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c13
	t12.Mul(t4, t12)

	// Step 562: t12 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e098
	for s := 0; s < 3; s++ {
		t12.Square(t12)
	}

	// Step 563: t12 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e099
	t12.Mul(x, t12)

	// Step 573: t12 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf826400
	for s := 0; s < 10; s++ {
		t12.Square(t12)
	}

	// Step 574: t12 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f
	t12.Mul(t10, t12)

	// Step 581: t12 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f80
	for s := 0; s < 7; s++ {
		t12.Square(t12)
	}

	// Step 582: t12 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f99
	t12.Mul(t3, t12)

	// Step 588: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e640
	for s := 0; s < 6; s++ {
		t12.Square(t12)
	}

	// Step 589: t12 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65b
	t12.Mul(t11, t12)

	// Step 594: t12 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb60
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 595: t11 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7b
	t11.Mul(t11, t12)

	// Step 598: t11 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bd8
	for s := 0; s < 3; s++ {
		t11.Square(t11)
	}

	// Step 599: t11 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf
	t11.Mul(t6, t11)

	// Step 603: t11 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf0
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 604: t11 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf1
	t11.Mul(x, t11)

	// Step 611: t11 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def880
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 612: t11 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def887
	t11.Mul(t6, t11)

	// Step 620: t11 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def88700
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 621: t11 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def88715
	t11.Mul(t0, t11)

	// Step 632: t11 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a800
	for s := 0; s < 11; s++ {
		t11.Square(t11)
	}

	// Step 633: t11 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f
	t11.Mul(t8, t11)

	// Step 637: t11 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f0
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 638: t11 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f7
	t11.Mul(t6, t11)

	// Step 648: t11 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc00
	for s := 0; s < 10; s++ {
		t11.Square(t11)
	}

	// Step 649: t10 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f
	t10.Mul(t10, t11)

	// Step 656: t10 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f80
	for s := 0; s < 7; s++ {
		t10.Square(t10)
	}

	// Step 657: t10 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f
	t10.Mul(t7, t10)

	// Step 662: t10 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e0
	for s := 0; s < 5; s++ {
		t10.Square(t10)
	}

	// Step 663: t9 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e9
	t9.Mul(t9, t10)

	// Step 669: t9 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a40
	for s := 0; s < 6; s++ {
		t9.Square(t9)
	}

	// Step 670: t9 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a59
	t9.Mul(t3, t9)

	// Step 676: t9 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e9640
	for s := 0; s < 6; s++ {
		t9.Square(t9)
	}

	// Step 677: t8 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f
	t8.Mul(t8, t9)

	// Step 682: t8 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbe0
	for s := 0; s < 5; s++ {
		t8.Square(t8)
	}

	// Step 683: t7 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef
	t7.Mul(t7, t8)

	// Step 689: t7 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbc0
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 690: t7 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb
	t7.Mul(z, t7)

	// Step 696: t7 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2c0
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 697: t7 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d9
	t7.Mul(t3, t7)

	// Step 703: t7 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb640
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 704: t7 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653
	t7.Mul(t4, t7)

	// Step 707: t7 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b298
	for s := 0; s < 3; s++ {
		t7.Square(t7)
	}

	// Step 708: t6 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f
	t6.Mul(t6, t7)

	// Step 717: t6 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e00
	for s := 0; s < 9; s++ {
		t6.Square(t6)
	}

	// Step 718: t6 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e03
	t6.Mul(t5, t6)

	// Step 723: t6 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c060
	for s := 0; s < 5; s++ {
		t6.Square(t6)
	}

	// Step 724: t5 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063
	t5.Mul(t5, t6)

	// Step 744: t5 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c06300000
	for s := 0; s < 20; s++ {
		t5.Square(t5)
	}

	// Step 745: t5 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c06300019
	t5.Mul(t3, t5)

	// Step 750: t5 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c6000320
	for s := 0; s < 5; s++ {
		t5.Square(t5)
	}

	// Step 751: t4 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c6000333
	t4.Mul(t4, t5)

	// Step 756: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c0006660
	for s := 0; s < 5; s++ {
		t4.Square(t4)
	}

	// Step 757: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b
	t4.Mul(z, t4)

	// Step 763: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ac0
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 764: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad9
	t4.Mul(t3, t4)

	// Step 769: t4 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b20
	for s := 0; s < 5; s++ {
		t4.Square(t4)
	}

	// Step 770: t4 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d
	t4.Mul(t1, t4)

	// Step 779: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a00
	for s := 0; s < 9; s++ {
		t4.Square(t4)
	}

	// Step 780: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d
	t4.Mul(t1, t4)

	// Step 786: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968340
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 787: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968355
	t4.Mul(t0, t4)

	// Step 793: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d540
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 794: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d555
	t4.Mul(t0, t4)

	// Step 800: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968355540
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 801: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968355555
	t4.Mul(t0, t4)

	// Step 807: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d555540
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 808: t4 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d555555
	t4.Mul(t0, t4)

	// Step 814: t4 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968355555540
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 815: t3 = x^0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968355555559
	t3.Mul(t3, t4)

	// Step 825: t3 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d5555556400
	for s := 0; s < 10; s++ {
		t3.Square(t3)
	}

	// Step 826: t2 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d555555641d
	t2.Mul(t2, t3)

	// Step 831: t2 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83a0
	for s := 0; s < 5; s++ {
		t2.Square(t2)
	}

	// Step 832: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad
	t1.Mul(t1, t2)

	// Step 838: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb40
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 839: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55
	t1.Mul(t0, t1)

	// Step 845: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad540
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 846: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad555
	t1.Mul(t0, t1)

	// Step 852: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55540
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 853: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55555
	t1.Mul(t0, t1)

	// Step 859: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad555540
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 860: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad555555
	t1.Mul(t0, t1)

	// Step 866: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55555540
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 867: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55555555
	t1.Mul(t0, t1)

	// Step 873: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad555555540
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 874: t1 = x^0x505633a3138dae7b65786943aeea414d69d5124c52ce1fee0c9578a67a5f8f5002be3100b56c31c96f0208e8417dcee4bda12bf7e4c393b93cf82643f32def8871503ee1f8f4b2fbcb653e0318000ccd6cb41aaaaaaac83ad555555555
	t1.Mul(t0, t1)

	// Step 880: t1 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55555555540
	for s := 0; s < 6; s++ {
		t1.Square(t1)
	}

	// Step 881: t0 = x^0x14158ce8c4e36b9ed95e1a50ebba90535a75449314b387fb83255e299e97e3d400af8c402d5b0c725bc0823a105f73b92f684afdf930e4ee4f3e0990fccb7be21c540fb87e3d2cbef2d94f80c60003335b2d06aaaaaab20eb55555555555
	t0.Mul(t0, t1)

	// Step 886: t0 = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d555555641d6aaaaaaaaaaa0
	for s := 0; s < 5; s++ {
		t0.Square(t0)
	}

	// Step 887: z = x^0x282b19d189c6d73db2bc34a1d77520a6b4ea892629670ff7064abc533d2fc7a8015f18805ab618e4b781047420bee7725ed095fbf261c9dc9e7c1321f996f7c438a81f70fc7a597de5b29f018c000666b65a0d555555641d6aaaaaaaaaaab
	z.Mul(z, t0)

	return z
}
