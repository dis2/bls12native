package bls12

type GT = Fq12

func (r *GT) Add(p, q *GT) (ret *GT) {
	fq12Mul(r,p,q)
	return r
}

func (r *GT) Pair(p *G1, q *G2) (ret *GT) {
	return r.Pairs([]*G1{p}, []*G2{q})
}

func (r *GT) Pairs(p []*G1, q []*G2) (ret *GT) {
	var t = make([]G2, len(q))
	_p := make([]G1, len(p))
	_q := make([]G2, len(q))
	for i, pi := range p {
		_p[i].ToAffine(pi)
		_q[i].ToAffine(q[i])
	}
	fq12MillerLoop(r, t, _q, _p)
	fq12Conjugate(r, r)
	fq12FinalExp(r, r)
	return r
}
