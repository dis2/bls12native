package bls12

type G1 struct {
	X, Y, Z Fq
	Affine  bool
	VarTime bool
}

func (r *G1) Set(p *G1) *G1 {
	r.X.Set(&p.X)
	r.Y.Set(&p.Y)
	r.Z.Set(&p.Z)
	r.Affine = p.Affine
	return r
}

func (r *G1) SetZero() *G1 {
	r.X.SetZero()
	r.Y.SetZero()
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G1) IsZero() bool {
	return r.X.IsZero() && r.Y.IsZero()
}

func (r *G1) SetXY(x, y *Fq) *G1 {
	r.X.Set(x)
	r.Y.Set(y)
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G1) SetOne() *G1 {
	r.X = g1x
	r.Y = g1y
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G1) SetAffine() *G1 {
	r.Z.SetOne()
	r.Affine = true
	return r
}

func (r *G1) ToAffine(p *G1) *G1 {
	if p == nil {
		p = r
	}
	if p.Affine {
		return r
	}

	var t0, t1 Fq

	fqInv(&t1, &p.Z)
	fqSqr(&t0, &t1)
	fqMul(&r.X, &p.X, &t0)
	fqMul(&t0, &t0, &t1)
	fqMul(&r.Y, &p.Y, &t0)

	return r.SetAffine()
}

func (p *G1) Equal(q *G1) bool {
	var t Fq
	t02 := [2]*Fq{&p.X, &p.Y}
	t13 := [2]*Fq{&q.X, &q.Y}

	// t0 = p.X * q.Z^2
	// t1 = q.X * p.Z^2
	// t2 = p.Y * q.Z^3
	// t3 = q.Y * p.Z^3
	// t0 == t1 && t2 = t3
	if !q.Affine {
		var tfq [2]Fq
		t02[0], t02[1] = &tfq[0], &tfq[1]
		fqSqr(&t, &q.Z)
		fqMul(t02[0], &p.X, &t)
		fqMul(&t, &t, &q.Z)
		fqMul(t02[1], &p.Y, &t)
	}
	if !p.Affine {
		var tfq [2]Fq
		t13[0], t13[1] = &tfq[0], &tfq[1]
		fqSqr(&t, &p.Z)
		fqMul(t13[0], &q.X, &t)
		fqMul(&t, &t, &p.Z)
		fqMul(t13[1], &q.Y, &t)
	}
	return fqCmpN(t02[:], t13[:]) == 1
}

func (r *G1) G() *G1 {
	r.SetOne()
	return r
}

func (r *G1) Neg(p *G1) *G1 {
	if p == nil {
		p = r
	}
	r.Set(p)
	fqNeg(&r.Y, &r.Y)
	return r
}

func (r *G1) Double(p *G1) (ret *G1) {
	ret = r
	if p == nil {
		p = r
	}
	// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
	var t0, t1, t2, t3 Fq
	fqSqr(&t0, &p.X)
	fqSqr(&t1, &p.Y)
	fqSqr(&t2, &t1)
	fqAdd(&t1, &t1, &p.X)
	fqSqr(&t1, &t1)
	fqSub(&t1, &t1, &t0)
	fqSub(&t1, &t1, &t2)
	fqDbl(&t1, &t1)
	fqDbl(&t3, &t0)
	fqAdd(&t0, &t3, &t0)
	fqSqr(&t3, &t0)
	fqMul(&r.Z, &p.Y, &p.Z)
	fqDbl(&r.Z, &r.Z)
	fqSub(&r.X, &t3, &t1)
	fqSub(&r.X, &r.X, &t1)
	fqSub(&r.Y, &t1, &r.X)
	fqMul(&r.Y, &r.Y, &t0)
	fqDbl(&t2, &t2)
	fqDbl(&t2, &t2)
	fqDbl(&t2, &t2)
	fqSub(&r.Y, &r.Y, &t2)
	r.Affine = false
	return
}

func (r *G1) Mul(p *G1, s []uint64) {
	if p == nil {
		p = r
	}
	if p.VarTime {
		// TODO: wNAF
	}
	var t0, t1 G1
	t0.SetZero()
	t1.Set(p)

	// Platform code assumes the coords lie in sequence, ie only slice
	// length and first member address is consulted.
	t0sw := []*Fq{
		&t0.X, &t0.Y, &t0.Z,
	}
	t1sw := []*Fq{
		&t1.X, &t1.Y, &t1.Z,
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

func (r *G1) Add(p, q *G1) (ret *G1) {
	var t0, t1, t2, t3, t4, t5, t6 Fq
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
	if q.Affine {
		if !p.Affine {
			// t0 = z1^2.
			fqSqr(&t0, &p.Z)
			// t3 = U2 = x2 * z1^2.
			fqMul(&t3, &q.X, &t0)
			// t1 = S2 = y2 * z1^3.
			fqMul(&t1, &t0, &p.Z)
			fqMul(&t1, &t1, &q.Y)
			// t3 = H = U2 - x1.
			fqSub(&t3, &t3, &p.X)
			// t1 = R = 2 * (&S2 - y1).
			fqSub(&t1, &t1, &p.Y)
			fqDbl(&t1, &t1)
		} else {
			// H = x2 - x1.
			fqSub(&t3, &q.X, &p.X)
			// t1 = R = 2 * (&y2 - y1).
			fqSub(&t1, &q.Y, &p.Y)
			fqDbl(&t1, &t1)
		}
		// t2 = HH = H^2.
		fqSqr(&t2, &t3)
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
		fqDbl(&t4, &t2)
		fqDbl(&t4, &t4)
		// t5 = J = H * I.
		fqMul(&t5, &t3, &t4)
		// t4 = V = x1 * I.
		fqMul(&t4, &p.X, &t4)
		// x3 = R^2 - J - 2 * V.
		fqSqr(&r.X, &t1)
		fqSub(&r.X, &r.X, &t5)
		fqDbl(&t6, &t4)
		fqSub(&r.X, &r.X, &t6)
		// y3 = R * (&V - x3) - 2 * Y1 * J.
		fqSub(&t4, &t4, &r.X)
		fqMul(&t4, &t4, &t1)
		fqMul(&t1, &p.Y, &t5)
		fqDbl(&t1, &t1)
		fqSub(&r.Y, &t4, &t1)
		if !p.Affine {
			// z3 = (&z1 + H)^2 - z1^2 - HH.
			fqAdd(&r.Z, &p.Z, &t3)
			fqSqr(&r.Z, &r.Z)
			fqSub(&r.Z, &r.Z, &t0)
			fqSub(&r.Z, &r.Z, &t2)
		} else {
			// z3 = 2 * H.
			fqDbl(&r.Z, &t3)
		}
		p.Affine = false
		return
	}
	// t0 = z1^2.
	fqSqr(&t0, &p.Z)
	// t1 = z2^2.
	fqSqr(&t1, &q.Z)
	// t2 = U1 = x1 * z2^2.
	fqMul(&t2, &p.X, &t1)
	// t3 = U2 = x2 * z1^2.
	fqMul(&t3, &q.X, &t0)
	// t6 = z1^2 + z2^2.
	fqAdd(&t6, &t0, &t1)
	// t0 = S2 = y2 * z1^3.
	fqMul(&t0, &t0, &p.Z)
	fqMul(&t0, &t0, &q.Y)
	// t1 = S1 = y1 * z2^3.
	fqMul(&t1, &t1, &q.Z)
	fqMul(&t1, &t1, &p.Y)
	// t3 = H = U2 - U1.
	fqSub(&t3, &t3, &t2)
	// t0 = R = 2 * (&S2 - S1).
	fqSub(&t0, &t0, &t1)
	fqDbl(&t0, &t0)
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
	fqDbl(&t4, &t3)
	fqSqr(&t4, &t4)
	// t5 = J = H * I.
	fqMul(&t5, &t3, &t4)
	// t4 = V = U1 * I.
	fqMul(&t4, &t2, &t4)
	// x3 = R^2 - J - 2 * V.
	fqSqr(&r.X, &t0)
	fqSub(&r.X, &r.X, &t5)
	fqDbl(&t2, &t4)
	fqSub(&r.X, &r.X, &t2)
	// y3 = R * (&V - x3) - 2 * S1 * J.
	fqSub(&t4, &t4, &r.X)
	fqMul(&t4, &t4, &t0)
	fqMul(&t1, &t1, &t5)
	fqDbl(&t1, &t1)
	fqSub(&r.Y, &t4, &t1)
	// z3 = (&(&z1 + z2)^2 - z1^2 - z2^2) * H.
	fqAdd(&r.Z, &p.Z, &q.Z)
	fqSqr(&r.Z, &r.Z)
	fqSub(&r.Z, &r.Z, &t6)
	fqMul(&r.Z, &r.Z, &t3)
	r.Affine = false
	return
}
