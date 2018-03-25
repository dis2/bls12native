package bls12

type G2 struct {
	X, Y, Z Fq2
	Affine  bool
	VarTime bool
}

func (r *G2) Set(p *G2) *G2 {
	r.X.Set(&p.X)
	r.Y.Set(&p.Y)
	r.Z.Set(&p.Z)
	r.Affine = p.Affine
	return r
}

func (r *G2) SetZero() *G2 {
	r.X.SetZero()
	r.Y.SetZero()
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (p *G2) Equal(q *G2) bool {
	var t Fq2
	t02 := [4]*Fq{&p.X[0], &p.X[1], &p.Y[0], &p.Y[1]}
	t13 := [4]*Fq{&q.X[0], &q.X[1], &q.Y[0], &q.Y[1]}

	// t0 = p.X * q.Z^2
	// t1 = q.X * p.Z^2
	// t2 = p.Y * q.Z^3
	// t3 = q.Y * p.Z^3
	// t0 == t1 && t2 = t3
	if !q.Affine {
		var tfq [2]Fq2
		t02[0], t02[1], t02[2], t02[3] = &tfq[0][0], &tfq[0][1], &tfq[1][0], &tfq[1][1]
		fq2Sqr(&t, &q.Z)
		fq2Mul(&tfq[0], &p.X, &t)
		fq2Mul(&t, &t, &q.Z)
		fq2Mul(&tfq[1], &p.Y, &t)
	}
	if !p.Affine {
		var tfq [2]Fq2
		t13[0], t13[1], t13[2], t13[3] = &tfq[0][0], &tfq[0][1], &tfq[1][0], &tfq[1][1]
		fq2Sqr(&t, &p.Z)
		fq2Mul(&tfq[0], &q.X, &t)
		fq2Mul(&t, &t, &p.Z)
		fq2Mul(&tfq[1], &q.Y, &t)
	}
	return fqCmpN(t02[:], t13[:]) == 1
}


func (r *G2) IsZero() bool {
	return r.X.IsZero() && r.Y.IsZero()
}

func (r *G2) SetXY(x, y *Fq2) *G2 {
	r.X.Set(x)
	r.Y.Set(y)
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G2) SetOne() *G2 {
	r.X[0] = g2x0
	r.X[1] = g2x1
	r.Y[0] = g2y0
	r.Y[1] = g2y1
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G2) G() *G2 {
	r.SetOne()
	return r
}

func (r *G2) SetAffine() *G2 {
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G2) Neg(p *G2) *G2 {
	if p == nil {
		p = r
	}
	r.Set(p)
	fq2Neg(&r.Y, &r.Y)
	return r
}

func (r *G2) Double(p *G2) (ret *G2) {
	ret = r
	if p == nil {
		p = r
	}
	var t0, t1, t2, t3 Fq2
	fq2Sqr(&t0, &p.X);
	fq2Add(&t2, &t0, &t0);
	fq2Add(&t0, &t2, &t0);
	fq2Sqr(&t3, &p.Y);
	fq2Mul(&t1, &t3, &p.X);
	fq2Add(&t1, &t1, &t1);
	fq2Add(&t1, &t1, &t1);
	fq2Sqr(&r.X, &t0);
	fq2Add(&t2, &t1, &t1);
	fq2Sub(&r.X, &r.X, &t2);
	fq2Mul(&r.Z, &p.Z, &p.Y);
	fq2Add(&r.Z, &r.Z, &r.Z);
	fq2Add(&t3, &t3, &t3);
	fq2Sqr(&t3, &t3);
	fq2Add(&t3, &t3, &t3);
	fq2Sub(&t1, &t1, &r.X);
	fq2Mul(&r.Y, &t0, &t1);
	fq2Sub(&r.Y, &r.Y, &t3);
	r.Affine = false
	return
}

func (r *G2) Add(p, q *G2) (ret *G2) {
	var t0, t1, t2, t3, t4, t5, t6 Fq2
	ret = r
	if p == nil {
		p = r
	}
	if p.IsZero() {
		*r = *q
		return
	}
	if q.IsZero() {
		*r = *p
		return
	}

	if *p == *q {
		if !p.Equal(q) {
			panic("point equal failed")
		}
	}

	// Do this early
	if p.Equal(q) {
		r.Double(p)
		return
	}

	if q.Affine {
		if !p.Affine {
			// t0 = z1^2.
			fq2Sqr(&t0, &p.Z)
			// t3 = U2 = x2 * z1^2.
			fq2Mul(&t3, &q.X, &t0)
			// t1 = S2 = y2 * z1^3.
			fq2Mul(&t1, &t0, &p.Z)
			fq2Mul(&t1, &t1, &q.Y)
			// t3 = H = U2 - x1.
			fq2Sub(&t3, &t3, &p.X)
			// t1 = R = 2 * (&S2 - y1).
			fq2Sub(&t1, &t1, &p.Y)
		} else {
			// H = x2 - x1.
			fq2Sub(&t3, &q.X, &p.X)
			// t1 = R = 2 * (&y2 - y1).
			fq2Sub(&t1, &q.Y, &p.Y)
		}
		// t2 = HH = H^2.
		fq2Sqr(&t2, &t3)
		// If E is zero.
		if t3.IsZero() {
			if t1.IsZero() {
				r.Double(p)
			} else {
				r.SetZero()
			}
			return
		}
		// t4 = I = 4*HH.
		fq2Mul(&t5, &t3, &t2)
		/* t4 = V = x1 * HH. */
		fq2Mul(&t4, &p.X, &t2)
		// x3 = R^2 - J - 2 * V.
		fq2Sqr(&r.X, &t1)
		fq2Sub(&r.X, &r.X, &t5)
		fq2Dbl(&t6, &t4)
		fq2Sub(&r.X, &r.X, &t6)
		// y3 = R * (&V - x3) - 2 * Y1 * J.
		fq2Sub(&t4, &t4, &r.X)
		fq2Mul(&t4, &t4, &t1)
		fq2Mul(&t1, &p.Y, &t5)
		fq2Sub(&r.Y, &t4, &t1)
		if (!p.Affine) {
			/* z3 = z1 * H. */
			fq2Mul(&r.Z, &p.Z, &t3);
		} else {
			/* z3 = H. */
			r.Z.Set(&t3)
		}

		p.Affine = false
		return
	}
	// t0 = z1^2.
	fq2Sqr(&t0, &p.Z)
	// t1 = z2^2.
	fq2Sqr(&t1, &q.Z)
	// t2 = U1 = x1 * z2^2.
	fq2Mul(&t2, &p.X, &t1)
	// t3 = U2 = x2 * z1^2.
	fq2Mul(&t3, &q.X, &t0)
	// t6 = z1^2 + z2^2.
	fq2Add(&t6, &t0, &t1)
	// t0 = S2 = y2 * z1^3.
	fq2Mul(&t0, &t0, &p.Z)
	fq2Mul(&t0, &t0, &q.Y)
	// t1 = S1 = y1 * z2^3.
	fq2Mul(&t1, &t1, &q.Z)
	fq2Mul(&t1, &t1, &p.Y)
	// t3 = H = U2 - U1.
	fq2Sub(&t3, &t3, &t2)
	// t0 = R = 2 * (&S2 - S1).
	fq2Sub(&t0, &t0, &t1)
	fq2Dbl(&t0, &t0)
	// If E is zero.
	if t3.IsZero() {
		if t0.IsZero() {
			r.Double(p)
		} else {
			r.SetZero()
		}
		return
	}
	// t4 = I = (&2*H)^2.
	fq2Dbl(&t4, &t3)
	fq2Sqr(&t4, &t4)
	// t5 = J = H * I.
	fq2Mul(&t5, &t3, &t4)
	// t4 = V = U1 * I.
	fq2Mul(&t4, &t2, &t4)
	// x3 = R^2 - J - 2 * V.
	fq2Sqr(&r.X, &t0)
	fq2Sub(&r.X, &r.X, &t5)
	fq2Dbl(&t2, &t4)
	fq2Sub(&r.X, &r.X, &t2)
	// y3 = R * (&V - x3) - 2 * S1 * J.
	fq2Sub(&t4, &t4, &r.X)
	fq2Mul(&t4, &t4, &t0)
	fq2Mul(&t1, &t1, &t5)
	fq2Dbl(&t1, &t1)
	fq2Sub(&r.Y, &t4, &t1)
	// z3 = (&(&z1 + z2)^2 - z1^2 - z2^2) * H.
	fq2Add(&r.Z, &p.Z, &q.Z)
	fq2Sqr(&r.Z, &r.Z)
	fq2Sub(&r.Z, &r.Z, &t6)
	fq2Mul(&r.Z, &r.Z, &t3)
	r.Affine = false
	return
}

func (r *G2) Mul(p *G2, s []uint64) {
	if p == nil {
		p = r
	}
	if p.VarTime {
		// TODO: wNAF
	}
	var t0, t1 G2
	t0.SetZero()
	t1.Set(p)

	// Platform code assumes the coords lie in sequence, ie only slice
	// length and first member address is consulted.
	t0sw := []*Fq{
		&t0.X[0], &t0.X[1], &t0.Y[0], &t0.Y[1], &t0.Z[0], &t0.Z[1],
	}
	t1sw := []*Fq{
		&t1.X[0], &t1.X[1], &t1.Y[0], &t1.Y[1], &t1.Z[0], &t1.Z[1],
	}
	for i := uint(len(s)*64 - 1); int(i) >= 0; i-- {
		w := (s[i/64] >> (i % 64)) & 1
		fqSwap(t0sw, t1sw, w^1)
		t0.Add(&t0, &t1)
		t1.Double(&t1)
		fqSwap(t0sw, t1sw, w^1)
	}
	r.ToAffine(&t0)
}

func (r *G2) ToAffine(p *G2) *G2 {
	if p == nil {
		p = r
	}
	if p.Affine {
		return r
	}
	var t0, t1 Fq2
	fq2Inv(&t1, &p.Z)
	fq2Sqr(&t0, &t1)
	fq2Mul(&r.X, &p.X, &t0)
	fq2Mul(&t0, &t0, &t1)
	fq2Mul(&r.Y, &p.Y, &t0)

	return r.SetAffine()
}

//void pp_dbl_k12_projc_lazyr(fp12_t l, ep2_t r, ep2_t q, ep_t p) {

