package bls12

type Fq2 [2]Fq
type LazyFq2 [2]LazyFq

var zero2 = Fq2{}

func (t *Fq2) Set(e *Fq2) *Fq2 {
	*t = *e
	return t
}

func (t *Fq2) GetB() *Fq2 {
	t[0].GetB()
	t[1].GetB()
	return t
}

func (t *Fq2) SetZero() *Fq2 {
	*t = Fq2{}
	return t
}

func (t *Fq2) Equal(a *Fq2) bool {
	return fqCmpN([]*Fq{&t[0], &t[1]}, []*Fq{&a[0], &a[1]}) == 1
}

func (t *Fq2) IsZero() bool {
	return t.Equal(&zero2)
}

func (t *Fq2) SetOne() *Fq2 {
	t[0].SetOne()
	t[1].SetZero()
	return t
}

func (t *Fq2) SetUInt64(n uint64) *Fq2 {
	t[0].SetUInt64(n)
	t[1].SetZero()
	return t
}


func (e *Fq2) Copy() (r *Fq2) {
	r = new(Fq2)
	*r = *e
	return r
}

func fq2AddGeneric(c,a,b *Fq2) {
	fqAdd(&c[0], &a[0], &b[0])
	fqAdd(&c[1], &a[1], &b[1])
}

func lfq2AddGeneric(c,a,b *LazyFq2) {
	lfqAdd(&c[0], &a[0], &b[0])
	lfqAdd(&c[1], &a[1], &b[1])
}

func lfq2SubGeneric(c,a,b *LazyFq2) {
	lfqSub(&c[0], &a[0], &b[0])
	lfqSub(&c[1], &a[1], &b[1])
}

func fq2SubGeneric(c,a,b *Fq2) {
	fqSub(&c[0], &a[0], &b[0])
	fqSub(&c[1], &a[1], &b[1])
}

func fq2NegGeneric(c,a *Fq2) {
	fqNeg(&c[0], &a[0])
	fqNeg(&c[1], &a[1])
}

func lfq2MulGeneric(c *LazyFq2,a,b *Fq2) {
	var t0,t1 Fq
	var lt2,lt0 LazyFq
	fqAdd(&t0, &a[0], &a[1]);
	fqAdd(&t1, &b[0], &b[1]);
	lfqMul(&c[0], &a[0], &b[0])
	lfqMul(&c[1], &a[1], &b[1]);
	lfqMul(&lt2, &t0, &t1);
	lfqAdd(&lt0, &c[0], &c[1]);
	lfqSub(&c[0], &c[0], &c[1]);
	lfqSub(&c[1], &lt2, &lt0)
}

func fq2REDCGeneric(c *Fq2, lc *LazyFq2) {
	fqREDC(&c[0], &lc[0])
	fqREDC(&c[1], &lc[1])
}

func fq2MulGeneric(c,a,b *Fq2) {
	var lc LazyFq2
	lfq2Mul(&lc, a, b)
	fq2REDC(c, &lc)
}

func fq2Inv(c, a *Fq2) {
	var t0, t1 Fq
	/* t0 = a_0^2, t1 = a_1^2. */
	fqSqr(&t0, &a[0]);
	fqSqr(&t1, &a[1]);
	/* t1 = 1/(a_0^2 + a_1^2). */
	fqAdd(&t0, &t0, &t1);
	fqInv(&t1, &t0);
	/* c_0 = a_0/(a_0^2 + a_1^2). */
	fqMul(&c[0], &a[0], &t1);
	/* c_1 = - a_1/(a_0^2 + a_1^2). */
	fqMul(&c[1], &a[1], &t1);
	fqNeg(&c[1], &c[1]);
}

func fq2Sqr(c,a *Fq2) {
	fq2Mul(c,a,a)
}

func lfq2Sqr(c *LazyFq2 ,a *Fq2) {
	lfq2Mul(c,a,a)
}

func fq2Dbl(c,a *Fq2) {
	fqDbl(&c[0], &a[0])
	fqDbl(&c[1], &a[1])
}

func lfq2MulQNR(a,c *LazyFq2) { // norh
	/* (a_0 + a_1 * u) * (1 + u) = (a_0 - a_1) + (a_0 + a_1) * u. */
	t := a[1]
	lfqAdd(&c[1], &a[0], &a[1]);
	lfqSub(&c[0], &a[0], &t);
}

func fq2MulQNR(a,c *Fq2) { // norh
	/* (a_0 + a_1 * u) * (1 + u) = (a_0 - a_1) + (a_0 + a_1) * u. */
	t := a[1]
	fqAdd(&c[1], &a[0], &a[1]);
	fqSub(&c[0], &a[0], &t);
}

func fq2MapFrb(c,a *Fq2, i int) {
	switch i%2 {
	case 0:
		*c = *a
	case 1:
		c[0] = a[0]
		fqNeg(&c[1], &a[1])
	}
}

func fq2MulFrb(c,a *Fq2, i,j int) {
	switch i {
	case 1:
		fq2Mul(c, a, &frbtab[j - 1]);
	case 2:
		fqMul(&c[0], &a[0], &frbtab2[j - 1]);
		fqMul(&c[1], &a[1], &frbtab2[j - 1]);
	default:
		fq2Mul(c, a, &frbtab3[j - 1]);
	}
}
