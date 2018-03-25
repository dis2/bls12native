package bls12

type Fq12 [2]Fq6

func (t *Fq12) Set(e *Fq12) *Fq12 {
	*t = *e
	return t
}

func fq12MapFrb(c, a *Fq12, i int) {
	switch i {
	case 0:
		c.Set(a)
	case 1:
		fq2MapFrb(&c[0][0], &a[0][0], 1);
		fq2MapFrb(&c[1][0], &a[1][0], 1);
		fq2MapFrb(&c[0][1], &a[0][1], 1);
		fq2MapFrb(&c[1][1], &a[1][1], 1);
		fq2MapFrb(&c[0][2], &a[0][2], 1);
		fq2MapFrb(&c[1][2], &a[1][2], 1);
		fq2MulFrb(&c[1][0], &c[1][0], 1, 1);
		fq2MulFrb(&c[0][1], &c[0][1], 1, 2);
		fq2MulFrb(&c[1][1], &c[1][1], 1, 3);
		fq2MulFrb(&c[0][2], &c[0][2], 1, 4);
		fq2MulFrb(&c[1][2], &c[1][2], 1, 5);
	case 2:
		c[0][0].Set(&a[0][0])
		fq2MulFrb(&c[0][2], &a[0][2], 2, 1);
		fq2MulFrb(&c[0][1], &a[0][1], 2, 2);
		fq2Neg(&c[0][2], &c[0][2]);
		fq2MulFrb(&c[1][0], &a[1][0], 2, 1);
		fq2MulFrb(&c[1][2], &a[1][2], 2, 2);
		fq2MulFrb(&c[1][1], &a[1][1], 2, 3);
		fq2Neg(&c[1][2], &c[1][2]);
	case 3:
		fq2MapFrb(&c[0][0], &a[0][0], 1);
		fq2MapFrb(&c[1][0], &a[1][0], 1);
		fq2MapFrb(&c[0][1], &a[0][1], 1);
		fq2MapFrb(&c[1][1], &a[1][1], 1);
		fq2MapFrb(&c[0][2], &a[0][2], 1);
		fq2MapFrb(&c[1][2], &a[1][2], 1);
		fq2MulFrb(&c[0][1], &c[0][1], 3, 2);
		fq2MulFrb(&c[0][2], &c[0][2], 3, 4);
		fq2Neg(&c[0][2], &c[0][2]);
		fq2MulFrb(&c[1][0], &c[1][0], 3, 1);
		fq2MulFrb(&c[1][1], &c[1][1], 3, 3);
		fq2MulFrb(&c[1][2], &c[1][2], 3, 5);
		fq2Neg(&c[1][2], &c[1][2]);
	}
}


func fq12Conjugate(c,a *Fq12){
	c[0] = a[0]
	fq6Neg(&c[1], &a[1])
}

func fq12Exp(c, a *Fq12, s []uint64) {
	nw := len(s)
	base := *a
	var res Fq12
	res[0][0][0].SetOne()
	for i := 0; i < nw; i++ {
		w := s[i]
		for j := uint(0); j < 64; j++ {
			if w & 1 != 0 {
				fq12Mul(&res, &res, &base)
			}
			w >>= 1
			if w == 0 && i == nw-1 {
				break
			}
			fq12Sqr(&base, &base)
		}
	}
	*c = res
}

// assume negative sign of the exponent
func fq12ExpNeg(c, a *Fq12, s []uint64) {
	fq12Exp(c,a,s)
	fq12Conjugate(c,c)
}


//void fp12_mul_lazyr(fp12_t c, fp12_t a, fp12_t b) {
func fq12Mul(c,a,b *Fq12) {
	var u0, u1, u2, u3 LazyFq6
	var t0, t1 Fq6

	/* Karatsuba algorithm. */
	/* u0 = a_0 * b_0. */
	lfq6Mul(&u0, &a[0], &b[0]);
	/* u1 = a_1 * b_1. */
	lfq6Mul(&u1, &a[1], &b[1]);
	/* t1 = a_0 + a_1. */
	fq6Add(&t0, &a[0], &a[1]);
	/* t0 = b_0 + b_1. */
	fq6Add(&t1, &b[0], &b[1]);
	/* u2 = (&a_0 + a_1) * (&b_0 + b_1) */
	lfq6Mul(&u2, &t0, &t1);
	/* c_1 = u2 - a_0b_0 - a_1b_1. */
	for i := 0; i < 3; i++ {
		lfq2Add(&u3[i], &u0[i], &u1[i]);
		lfq2Sub(&u2[i], &u2[i], &u3[i]);
		fq2REDC(&c[1][i], &u2[i]);
	}
	/* c_0 = a_0b_0 + v * a_1b_1. */
	lfq2MulQNR(&u2[0], &u1[2]);
	u2[1][0] = u1[0][0]
	u2[1][1] = u1[0][1]
	u2[2][0] = u1[1][0]
	u2[2][1] = u1[1][1]
	for i := 0; i < 3; i++ {
		lfq2Add(&u2[i], &u0[i], &u2[i]);
		fq2REDC(&c[0][i], &u2[i]);
	}
}


//void fp12_sqr_cyc_lazyr(fp12_t c, fp12_t a) {
func fq12Sqr(c, a *Fq12) {
	var t0, t1 Fq2
	var u0, u1, u2, u3, u4 LazyFq2

	lfq2Sqr(&u2, &a[0][0]);
	lfq2Sqr(&u3, &a[1][1]);
	fq2Add(&t1, &a[0][0], &a[1][1]);

	lfq2MulQNR(&u0, &u3);
	lfq2Add(&u0, &u0, &u2);
	fq2REDC(&t0, &u0);

	lfq2Sqr(&u1, &t1);
	lfq2Add(&u2, &u2, &u3); // fp2_addd ??!
	lfq2Sub(&u1, &u1, &u2);
	fq2REDC(&t1, &u1);

	fq2Sub(&c[0][0], &t0, &a[0][0]);
	fq2Add(&c[0][0], &c[0][0], &c[0][0]);
	fq2Add(&c[0][0], &t0, &c[0][0]);

	fq2Add(&c[1][1], &t1, &a[1][1]);
	fq2Add(&c[1][1], &c[1][1], &c[1][1]);
	fq2Add(&c[1][1], &t1, &c[1][1]);

	lfq2Sqr(&u0, &a[0][1]);
	lfq2Sqr(&u1, &a[1][2]);
	fq2Add(&t0, &a[0][1], &a[1][2]);
	lfq2Sqr(&u2, &t0);

	lfq2Add(&u3, &u0, &u1);
	lfq2Sub(&u3, &u2, &u3);
	fq2REDC(&t0, &u3);

	fq2Add(&t1, &a[1][0], &a[0][2]);
	lfq2Sqr(&u3, &t1);
	lfq2Sqr(&u2, &a[1][0]);

	fq2MulQNR(&t1, &t0);
	fq2Add(&t0, &t1, &a[1][0]);
	fq2Add(&t0, &t0, &t0);
	fq2Add(&c[1][0], &t0, &t1);

	lfq2MulQNR(&u4, &u1);
	lfq2Add(&u4, &u0, &u4);
	fq2REDC(&t0, &u4);
	fq2Sub(&t1, &t0, &a[0][2]);

	lfq2Sqr(&u1, &a[0][2]);

	fq2Add(&t1, &t1, &t1);
	fq2Add(&c[0][2], &t1, &t0);

	lfq2MulQNR(&u4, &u1);
	lfq2Add(&u4, &u2, &u4);
	fq2REDC(&t0, &u4);
	fq2Sub(&t1, &t0, &a[0][1]);
	fq2Add(&t1, &t1, &t1);
	fq2Add(&c[0][1], &t1, &t0);

	lfq2Add(&u0, &u2, &u1);
	lfq2Sub(&u3, &u3, &u0);
	fq2REDC(&t0, &u3);
	fq2Add(&t1, &t0, &a[1][2]);
	fq2Dbl(&t1, &t1);
	fq2Add(&c[1][2], &t0, &t1);
}

func fq12Inv(c,a *Fq12) {
	var t0, t1 Fq6
	fq6Sqr(&t0, &a[0]);
	fq6Sqr(&t1, &a[1]);
	fq6MulAdjSqrt(&t1, &t1);
	fq6Sub(&t0, &t0, &t1);
	fq6Inv(&t0, &t0);

	fq6Mul(&c[0], &a[0], &t0);
	fq6Neg(&c[1], &a[1]);
	fq6Mul(&c[1], &c[1], &t0);
}

// static void pp_exp_b12(fp12_t c, fp12_t a) {
func fq12FinalExp(c,a *Fq12) {
	var t, t0, t1, t2, t3 Fq12
	/* First, compute c = a^(p^6 - 1). */
	/* t = a^{-1}. */
	fq12Inv(&t, a);
	/* c = a^(p^6). */
	fq12Conjugate(c, a);
	/* c = a^(p^6 - 1). */
	fq12Mul(c, c, &t);

	/* Second, compute c^(p^2 + 1). */
	/* t = c^(p^2). */
	fq12MapFrb(&t, c, 2)
//	fp12_frb(t, c, 1);
//	fp12_frb(t, t, 1);

	/* c = c^(p^2 + 1). */
	fq12Mul(c, c, &t);

	/* Now compute m^((p^4 - p^2 + 1) / r). */
	/* t0 = f^2. */
	fq12Sqr(&t0, c);

	/* t1 = f^x. */
	xp := []uint64{xbits}
	fq12ExpNeg(&t1, c, xp)

	/* t2 = f^(x^2). */
	fq12ExpNeg(&t2, &t1, xp);

	/* t1 = t2/(t1^2 * f). */
	fq12Conjugate(&t3, c);
	fq12Sqr(&t1, &t1);
	fq12Mul(&t1, &t1, &t3);
	fq12Conjugate(&t1, &t1);
	fq12Mul(&t1, &t1, &t2);

	/* t2 = t1^x. */
	fq12ExpNeg(&t2, &t1, xp)

	/* t3 = t2^x/t1. */
	fq12ExpNeg(&t3, &t2, xp);

	fq12Conjugate(&t1, &t1);
	fq12Mul(&t3, &t1, &t3);

	/* t1 = t1^(-p^3 ) * t2^(p^2). */
	fq12Conjugate(&t1, &t1);
	fq12MapFrb(&t1, &t1, 3);
	fq12MapFrb(&t2, &t2, 2);
	fq12Mul(&t1, &t1, &t2);

	/* t2 = f * f^2 * t3^x. */
	fq12ExpNeg(&t2, &t3, xp);

	fq12Mul(&t2, &t2, &t0);
	fq12Mul(&t2, &t2, c);

	/* Compute t1 * t2 * t3^p. */
	fq12Mul(&t1, &t1, &t2);
	fq12MapFrb(&t2, &t3, 1);
	fq12Mul(c, &t1, &t2);
}


func fq12MulSparse(c,a,b *Fq12) {
	var u0, u1, u2 LazyFq6
	var t0 Fq6
	/* t0 = a_0 * b_0. */
	lfq6MulSparse(&u0, &a[0], &b[0]);
	/* t1 = a_1 * b_1. */
	lfq2Mul(&u1[1], &a[1][2], &b[1][1]);
	lfq2MulQNR(&u1[0], &u1[1]);
	lfq2Mul(&u1[1], &a[1][0], &b[1][1]);
	lfq2Mul(&u1[2], &a[1][1], &b[1][1]);
	/* t2 = b_0 + b_1. */
	t0[0] = b[0][0]
	fq2Add(&t0[1], &b[0][1], &b[1][1]);
	/* c_1 = a_0 + a_1. */
	fq6Add(&c[1], &a[0], &a[1]);
	/* c_1 = (a_0 + a_1) * (&b_0 + b_1) */
	lfq6MulSparse(&u2, &c[1], &t0);
	for i := 0; i < 3; i++ {
		lfq2Sub(&u2[i], &u2[i], &u0[i]);
		lfq2Sub(&u2[i], &u2[i], &u1[i]);
	}
	fq2REDC(&c[1][0], &u2[0]);
	fq2REDC(&c[1][1], &u2[1]);
	fq2REDC(&c[1][2], &u2[2]);

	lfq2MulQNR(&u2[0], &u1[2]);
	lfq2Add(&u0[0], &u0[0], &u2[0]);
	lfq2Add(&u0[1], &u0[1], &u1[0]);
	lfq2Add(&u0[2], &u0[2], &u1[1]);
	/* c_0 = a_0b_0 + v * a_1b_1. */
	fq2REDC(&c[0][0], &u0[0]);
	fq2REDC(&c[0][1], &u0[1]);
	fq2REDC(&c[0][2], &u0[2]);
}

// line function evalautions, results go in l and r
func fq12MillerDbl(l *Fq12, r, q *G2, p *G1) {
	var t0, t1, t2, t3, t4, t5, t6 Fq2

	// TODO: make it LazyFq2
	//var u0, u1 LazyFq2

	/* A = x1^2. */
	fq2Sqr(&t0, &q.X);
	/* B = y1^2. */
	fq2Sqr(&t1, &q.Y);
	/* C = z1^2. */
	fq2Sqr(&t2, &q.Z);
	/* D = 3bC, for general b. */
	fq2Dbl(&t3, &t2);
	fq2Add(&t3, &t3, &t2);
	t4.GetB();
	fq2Mul(&t3, &t3, &t4);
	/* E = (&x1 + y1)^2 - A - B. */
	fq2Add(&t4, &q.X, &q.Y);
	fq2Sqr(&t4, &t4);
	fq2Sub(&t4, &t4, &t0);
	fq2Sub(&t4, &t4, &t1);

	/* F = (&y1 + z1)^2 - B - C. */
	fq2Add(&t5, &q.Y, &q.Z);
	fq2Sqr(&t5, &t5);
	fq2Sub(&t5, &t5, &t1);
	fq2Sub(&t5, &t5, &t2);

	/* G = 3D. */
	fq2Dbl(&t6, &t3);
	fq2Add(&t6, &t6, &t3);

	/* x3 = E * (&B - G). */
	fq2Sub(&r.X, &t1, &t6);
	fq2Mul(&r.X, &r.X, &t4);

	/* y3 = (&B + G)^2 -12D^2. */
	fq2Add(&t6, &t6, &t1);
	fq2Sqr(&t6, &t6);
	fq2Sqr(&t2, &t3);
	fq2Dbl(&r.Y, &t2);
	fq2Dbl(&t2, &r.Y);
	fq2Dbl(&r.Y, &t2);
	fq2Add(&r.Y, &r.Y, &t2);
	fq2Sub(&r.Y, &t6, &r.Y);

	/* z3 = 4B * F. */
	fq2Dbl(&r.Z, &t1);
	fq2Dbl(&r.Z, &r.Z);
	fq2Mul(&r.Z, &r.Z, &t5);

	const one = 1
	const zero = 0

	/* l00 = D - B. */
	fq2Sub(&l[one][one], &t3, &t1);

	/* l10 = (&3 * xp) * A. */
	fqMul(&l[one][zero][0], &p.X, &t0[0]);
	fqMul(&l[one][zero][1], &p.X, &t0[1]);

	/* l01 = F * (&-yp). */
	fqMul(&l[zero][zero][0], &t5[0], &p.Y);
	fqMul(&l[zero][zero][1], &t5[1], &p.Y);
	r.Affine = false
}

func fq12MillerAdd(l *Fq12, r, q *G2, p *G1) {
	const one = 1
	const zero = 0
	var t1, t2, t3, t4 Fq2
	var u1, u2 LazyFq2

	fq2Mul(&t1, &r.Z, &q.X);
	fq2Sub(&t1, &r.X, &t1);
	fq2Mul(&t2, &r.Z, &q.Y);
	fq2Sub(&t2, &r.Y, &t2);

	fq2Sqr(&t3, &t1);
	fq2Mul(&r.X, &t3, &r.X);
	fq2Mul(&t3, &t1, &t3);
	fq2Sqr(&t4, &t2);
	fq2Mul(&t4, &t4, &r.Z);
	fq2Add(&t4, &t3, &t4);

	fq2Sub(&t4, &t4, &r.X);
	fq2Sub(&t4, &t4, &r.X);
	fq2Sub(&r.X, &r.X, &t4);
	lfq2Mul(&u1, &t2, &r.X);
	lfq2Mul(&u2, &t3, &r.Y);
	lfq2Sub(&u2, &u1, &u2);
	fq2REDC(&r.Y, &u2);
	fq2Mul(&r.X, &t1, &t4);
	fq2Mul(&r.Z, &r.Z, &t3);

	fqMul(&l[one][zero][0], &t2[0], &p.X);
	fqMul(&l[one][zero][1], &t2[1], &p.X);
	fq2Neg(&l[one][zero], &l[one][zero]);

	lfq2Mul(&u1, &q.X, &t2);
	lfq2Mul(&u2, &q.Y, &t1);
	lfq2Sub(&u1, &u1, &u2);
	fq2REDC(&l[one][one], &u1);

	fqMul(&l[zero][zero][0], &t1[0], &p.Y);
	fqMul(&l[zero][zero][1], &t1[1], &p.Y);

	r.Affine = false
}

//pp_map_sim_oatep_k12
func fq12MillerLoop(r *Fq12, t []G2, q []G2, p []G1) {
	var l Fq12
	_p := make([]G1, len(p))
	_q := make([]G2, len(q))
	count := len(t)
	for i := 0; i < count; i++ {
		pi := &p[i]
		_q[i].Set(&q[i]).Neg(nil)
		pic := _p[i].Set(pi)
		// p.X*3
		fqDbl(&pic.X, &pi.X)
		fqAdd(&pic.X, &pic.X, &pi.X)
		fqNeg(&pic.Y, &pic.Y)
	}
	fq12MillerDbl(r, &t[0], &t[0], &_p[0])
	for i := 1; i < count; i++ {
		ti := &t[i]
		fq12MillerDbl(&l, ti, ti, &_p[i])
		fq12MulSparse(r, r, &l)
	}

	for i := 0; i < count; i++ {
		qi := &q[i]
		ti := &t[i]
		pi := &p[i]
		fq12MillerAdd(&l, ti, qi, pi);
		fq12MulSparse(r, r, &l)
	}
	for i := 61; i >= 0; i-- {
		fq12Sqr(r, r)
		for j := 0; j < len(t); j++ {
			tj := &t[j]
			qj := &q[j]
			pj := &p[j]
			_pj := &_p[j]
			fq12MillerDbl(&l, tj, tj, _pj);
			fq12MulSparse(r, r, &l)
			if (xbits & (1<<uint(i))) != 0 {
				fq12MillerAdd(&l, tj, qj, pj);
				fq12MulSparse(r, r, &l)
			}
		}
	}
}


