package bls12

type Fq6 [3]Fq2
type LazyFq6 [3]LazyFq2

func fq6Add(c, a, b *Fq6) {
	fq2Add(&c[0], &a[0], &b[0])
	fq2Add(&c[1], &a[1], &b[1])
	fq2Add(&c[2], &a[2], &b[2])
}

func fq6Sub(c, a, b *Fq6) {
	fq2Sub(&c[0], &a[0], &b[0])
	fq2Sub(&c[1], &a[1], &b[1])
	fq2Sub(&c[2], &a[2], &b[2])
}

func fq6Neg(c, a *Fq6) {
	fq2Neg(&c[0], &a[0])
	fq2Neg(&c[1], &a[1])
	fq2Neg(&c[2], &a[2])
}

func fq6Dbl(c, a *Fq6) {
	fq2Dbl(&c[0], &a[0])
	fq2Dbl(&c[1], &a[1])
	fq2Dbl(&c[2], &a[2])
}

func lfq6Sqr(c *LazyFq6,a *Fq6) {
	var u0, u1, u2, u3, u4, u5 LazyFq2
	var t0, t1, t2, t3 Fq2

	/* u0 = a_0^2 */
	lfq2Sqr(&u0, &a[0]);

	/* t1 = 2 * a_1 * a_2 */
	fq2Dbl(&t0, &a[1]);
	lfq2Mul(&u1, &t0, &a[2]);

	/* u2 = a_2^2. */
	lfq2Sqr(&u2, &a[2]);

	/* t4 = a_0 + a_2. */
	fq2Add(&t3, &a[0], &a[2]);

	/* u3 = (&a_0 + a_2 + a_1)^2. */
	fq2Add(&t2, &t3, &a[1]);
	lfq2Sqr(&u3, &t2);

	/* u4 = (&a_0 + a_2 - a_1)^2. */
	fq2Sub(&t1, &t3, &a[1]);
	lfq2Sqr(&u4, &t1);

	/* u4 = (&u4 + u3)/2. */
	lfq2Add(&u4, &u4, &u3);
	lfqHalve(&u4[0], &u4[0]);
	lfqHalve(&u4[1], &u4[1]);

	/* u3 = u3 - u4 - u1. */
	lfq2Add(&u5, &u1, &u4);
	lfq2Sub(&u3, &u3, &u5);

	/* c2 = u4 - u0 - u2. */
	lfq2Add(&u5, &u0, &u2);
	lfq2Sub(&c[2], &u4, &u5);

	/* c0 = u0 + u1 * E. */
	lfq2MulQNR(&u4, &u1);
	lfq2Add(&c[0], &u0, &u4);

	/* c1 = u3 + u2 * E. */
	lfq2MulQNR(&u4, &u2);
	lfq2Add(&c[1], &u3, &u4);
}


func lfq6Mul(c *LazyFq6, a,b *Fq6) {
	var u0, u1, u2, u3 LazyFq2
	var t0, t1 Fq2
	/* v0 = a_0b_0, v1 = a_1b_1, v2 = a_2b_2,
	 * t0 = a_1 + a_2, t1 = b_1 + b_2,
	 * u4 = u1 + u2, u5 = u0 + u1, u6 = u0 + u2 */
	lfq2Mul(&u0, &a[0], &b[0]);
	lfq2Mul(&u1, &a[1], &b[1]);
	lfq2Mul(&u2, &a[2], &b[2]);
	fq2Add(&t0, &a[1], &a[2]);
	fq2Add(&t1, &b[1], &b[2]);
	lfq2Add(&c[0], &u1, &u2);
	/* t2 (c_0) = v0 + E((a_1 + a_2)(b_1 + b_2) - v1 - v2) */
	lfq2Mul(&u3, &t0, &t1);
	lfq2Sub(&u3, &u3, &c[0]);
	lfq2MulQNR(&c[0], &u3);
	lfq2Add(&c[0], &c[0], &u0);
	/* c_1 = (a_0 + a_1)(b_0 + b_1) - v0 - v1 + Ev2 */
	fq2Add(&t0, &a[0], &a[1]);
	fq2Add(&t1, &b[0], &b[1]);
	lfq2Add(&c[1], &u0, &u1);
	lfq2Mul(&u3, &t0, &t1);
	lfq2Sub(&u3, &u3, &c[1]);
	lfq2MulQNR(&c[2], &u2);
	lfq2Add(&c[1], &u3, &c[2]);
	/* c_2 = (a_0 + a_2)(b_0 + b_2) - v0 + v1 - v2 */
	fq2Add(&t0, &a[0], &a[2]);
	fq2Add(&t1, &b[0], &b[2]);
	lfq2Add(&c[2], &u0, &u2);
	lfq2Mul(&u3, &t0, &t1);
	lfq2Sub(&u3, &u3, &c[2]);
	lfq2Add(&c[2], &u3, &u1);
}

func lfq6MulSparse(c *LazyFq6, a,b *Fq6) {
	var u0, u1, u2 LazyFq2
	var t0, t1 Fq2

	lfq2Mul(&u0, &a[0], &b[0]);
	lfq2Mul(&u1, &a[1], &b[1]);
	fq2Add(&t0, &a[0], &a[1]);
	fq2Add(&t1, &b[0], &b[1]);

	/* c_1 = (a_0 + a_1)(b_0 + b_1) - a_0b_0 - a_1b_1 */
	lfq2Mul(&u2, &t0, &t1);
	lfq2Sub(&u2, &u2, &u0);
	lfq2Sub(&c[1], &u2, &u1);

	/* c_0 = a_0b_0 + E a_2b_1 */
	lfq2Mul(&u2, &a[2], &b[1]);
	lfq2MulQNR(&c[0], &u2);
	lfq2Add(&c[0], &u0, &c[0]);

	/* c_2 = a_0b_2 + a_1b_1 */
	lfq2Mul(&u2, &a[2], &b[0]);
	lfq2Add(&c[2], &u1, &u2);
}


func fq6Mul(c,a,b *Fq6) {
	var t LazyFq6
	lfq6Mul(&t, a, b)
	fq2REDC(&c[0], &t[0])
	fq2REDC(&c[1], &t[1])
	fq2REDC(&c[2], &t[2])
}
func fq6Sqr(c,a *Fq6) {
	var t LazyFq6
	lfq6Sqr(&t, a)
	fq2REDC(&c[0], &t[0])
	fq2REDC(&c[1], &t[1])
	fq2REDC(&c[2], &t[2])
}


func fq6MulAdjSqrt(c,a *Fq6) {
	/* (a_0 + a_1 * v + a_2 * v^2) * v = a_2 + a_0 * v + a_1 * v^2 */
	t0 := a[0]
	fq2MulQNR(&c[0], &a[2]);
	c[2] = a[1]
	c[1] = t0
}

func fq6Inv(c,a *Fq6) {
	var v0,v1,v2,t0 Fq2
	/* v0 = a_0^2 - E * a_1 * a_2. */
	fq2Sqr(&t0, &a[0]);
	fq2Mul(&v0, &a[1], &a[2]);
	fq2MulQNR(&v2, &v0);
	fq2Sub(&v0, &t0, &v2);

	/* v1 = E * a_2^2 - a_0 * a_1. */
	fq2Sqr(&t0, &a[2]);
	fq2MulQNR(&v2, &t0);
	fq2Mul(&v1, &a[0], &a[1]);
	fq2Sub(&v1, &v2, &v1);

	/* v2 = a_1^2 - a_0 * a_2. */
	fq2Sqr(&t0, &a[1]);
	fq2Mul(&v2, &a[0], &a[2]);
	fq2Sub(&v2, &t0, &v2);

	fq2Mul(&t0, &a[1], &v2);
	fq2MulQNR(&c[1], &t0);

	fq2Mul(&c[0], &a[0], &v0);

	fq2Mul(&t0, &a[2], &v1);
	fq2MulQNR(&c[2], &t0);

	fq2Add(&t0, &c[0], &c[1]);
	fq2Add(&t0, &t0, &c[2]);
	fq2Inv(&t0, &t0);

	fq2Mul(&c[0], &v0, &t0);
	fq2Mul(&c[1], &v1, &t0);
	fq2Mul(&c[2], &v2, &t0);
}


