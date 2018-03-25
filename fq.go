package bls12

import "math/big"
//import "math/bits"

type Fq [6]uint64
type LazyFq [12]uint64

func parseInt(s string) (i *big.Int) {
	i = new(big.Int)
	i.SetString(s, 10)
	return
}

func (e *Fq) GetB() *Fq {
	*e = B
	return e
}

func parseFq(s string) (a Fq) {
	a.SetWords(parseInt(s).Bits())
	return
}

func parseFqMont(s string) (a Fq) {
	a.SetWords(parseInt(s).Bits())
	fqMul(&a, &a, &r2)
	return
}

// Set Fq to big.Int words.
func (e *Fq) SetWords(b []big.Word) {
	n := uint(len(b))
	if uint64(^uint(0))>>63 == 1 {
		for i := uint(0); i < n; i++ {
			e[i] = uint64(b[i])
		}
	} else {
		for i := uint(0); i < n; i++ {
			e[i/2] |= uint64(b[i]) << (32 * (i % 2))
		}
	}
}

// Get big.Int words representing the Fq.
func (e *Fq) GetWords() []big.Word {
	if uint64(^uint(0))>>63 == 1 {
		var b [6]big.Word
		for i := 0; i < 6; i++ {
			b[i] = big.Word(e[i])
		}
		return b[:]
	} else {
		var b [12]big.Word
		for i := uint(0); i < 12; i++ {
			b[i] = big.Word(e[i]) >> (32 * (i % 2))
		}
		return b[:]
	}
}

func (a *Fq) RawInt() (res *big.Int) {
	res = new(big.Int)
	res.SetBits(a.GetWords())
	return
}

func (e *Fq) FromInt(i *big.Int) *Fq {
	e.SetWords(i.Bits())
	fqMul(e, e, &r2)
	return e
}

func (e *Fq) SetUInt64(i uint64) *Fq {
	*e = Fq{i}
	fqMul(e, e, &r2)
	return e
}


func (e *Fq) ToInt() (res *big.Int) {
	res = new(big.Int)
	var t Fq
	fqREDC(&t, &LazyFq{e[0], e[1], e[2], e[3], e[4], e[5]})
	res.SetBits(t.GetWords())
	return
}

func (e *Fq) IsResidue() bool {
	var t Fq
	fqMul(&t, e, &qm1h)
	return t.IsOne()
}

// Guess all, or none: Branching on the full (degenerate) condition is fine
// provided the comparison leading to that is robust against time sampling
// of partial value.
func (e *Fq) Equal(a *Fq) bool {
	return fqCmp(e, a) == 1
}
func (e *Fq) IsOne() bool {
	return fqCmp(e, &one) == 1
}
func (e *Fq) IsZero() bool {
	return fqCmp(e, &zero) == 1
}

func (t *Fq) Set(e *Fq) *Fq {
	t[0], t[1], t[2], t[3], t[4], t[5] = e[0], e[1], e[2], e[3], e[4], e[5]
	return t
}

func (t *Fq) SetZero() *Fq {
	t[0], t[1], t[2], t[3], t[4], t[5] = 0, 0, 0, 0, 0, 0
	return t
}

func (t *Fq) SetOne() *Fq {
	t[0], t[1], t[2], t[3], t[4], t[5] = one[0], one[1], one[2], one[3], one[4], one[5]
	return t
}

func (e *Fq) Copy() (r *Fq) {
	t := new(Fq)
	t[0], t[1], t[2], t[3], t[4], t[5] = e[0], e[1], e[2], e[3], e[4], e[5]
	return t
}
func (e *LazyFq) Copy() (t *LazyFq) {
	t = new(LazyFq)
	*t = *e
	return
}

func fqReduceGeneric(a *Fq) {
	var b Fq
	var carry uint64
	for i, pi := range Q {
		ai := a[i]
		bi := ai - pi - carry
		b[i] = bi
		carry = (pi&^ai | (pi|^ai)&bi) >> 63
	}

	carry = -carry
	ncarry := ^carry
	for i := 0; i < 6; i++ {
		a[i] = (a[i] & carry) | (b[i] & ncarry)
	}
}

func lfqMulGeneric(cbuf *LazyFq, aa, bb *Fq) {
	var carry uint64

	for i := 0; i < 6; i++ {
		carry = 0
		b := aa[i]
		bh, bl := b>>32, b&mask32
		for j := 0; j < 6; j++ {
			c := bb[j]
			a := cbuf[j+i]

			ch, cl := c>>32, c&mask32
			ah, al := a>>32, a&mask32

			w := bl * cl
			x := bh * cl
			y := bl * ch

			r0 := (w & mask32) + al + (carry & mask32)

			z := bh * ch

			r1 := (r0 >> 32) + (w >> 32) + (x & mask32) + (y & mask32) + ah + (carry >> 32)
			r2 := (r1 >> 32) + (x >> 32) + (y >> 32) + (z & mask32)
			carry = (((r2 >> 32) + (z >> 32)) << 32) | (r2 & mask32)
			cbuf[i+j] = (r1 << 32) | (r0 & mask32)
		}
		cbuf[i+6] = carry
	}
}

func fqMulGeneric(c, a, b *Fq) {
	var cbuf LazyFq
	lfqMul(&cbuf, a, b)
	fqREDC(c, &cbuf)
}


func fqAddGeneric(c, a, b *Fq) {
	var carry uint64
	for i, ai := range a {
		bi := b[i]
		ci := ai + bi + carry
		c[i] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
	fqReduce(c)
}

func fqSubGeneric(c, a, b *Fq) {
	var t Fq
	fqNeg(&t, b)
	fqAdd(c, a, &t)
	fqReduce(c)
}

// BEWARE: abuf will be trashed, make a copy should you still need it.
func fqREDCGeneric(cbuf *Fq, abuf *LazyFq) {
	var carry, carry2 uint64
	for i := 0; i < 6; i++ {
		b := qInv64 * abuf[i]
		carry = 0
		bh, bl := b>>32, b&mask32
		for j := 0; j < 6; j++ {
			c := Q[j]
			a := abuf[i+j]

			ch, cl := c>>32, c&mask32
			ah, al := a>>32, a&mask32

			w := bl * cl
			x := bh * cl
			y := bl * ch

			r0 := (w & mask32) + al + (carry & mask32)

			z := bh * ch

			r1 := (r0 >> 32) + (w >> 32) + (x & mask32) + (y & mask32) + ah + (carry >> 32)
			r2 := (r1 >> 32) + (x >> 32) + (y >> 32) + (z & mask32)
			carry = (((r2 >> 32) + (z >> 32)) << 32) | (r2 & mask32)

			if j > 0 {
				abuf[i+j] = (r1 << 32) | (r0 & mask32)
			}
		}
		a := abuf[i+6]
		l := (a & mask32) + (carry2 & mask32) + (carry & mask32)
		h := (l >> 32) + (a >> 32) + (carry2 >> 32) + (carry >> 32)
		carry2 = h >> 32
		v := (h << 32) | (l & mask32)
		abuf[i+6] = v
	}
	cbuf[0] = abuf[6]
	cbuf[1] = abuf[7]
	cbuf[2] = abuf[8]
	cbuf[3] = abuf[9]
	cbuf[4] = abuf[10]
	cbuf[5] = abuf[11]
	fqReduce(cbuf)
}

func fqNegGeneric(c, a *Fq) {
	var carry uint64
	for i, pi := range Q {
		ai := a[i]
		ci := pi - ai - carry
		c[i] = ci
		carry = (ai&^pi | (ai|^pi)&ci) >> 63
	}
}

func lfqSubGeneric(c, a, b *LazyFq) {
	var carry uint64
	for i, pi := range a {
		bi := b[i]
		ti := pi - bi - carry
		c[i] = ti
		carry = (bi&^pi | (bi|^pi)&ti) >> 63
	}
	mask := -carry
	carry = 0
	for i, ai := range Q {
		ai &= mask
		bi := c[i+6]
		ci := ai + bi + carry
		c[i+6] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
}

func lfqAddGeneric(c, a, b *LazyFq) {
	var carry uint64
	for i, ai := range a {
		bi := b[i]
		ci := ai + bi + carry
		c[i] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
	var d Fq
	for i, pi := range Q {
		ai := c[i+6]
		bi := ai - pi - carry
		d[i] = bi
		carry = (pi&^ai | (pi|^ai)&bi) >> 63
	}

	carry = -carry
	ncarry := ^carry
	for i := 0; i < 6; i++ {
		c[i+6] = (c[i+6] & carry) | (d[i] & ncarry)
	}
}

func fqSqr(c, a *Fq) {
	fqMul(c, a, a)
}

func fqDbl(c, a *Fq) {
	fqAdd(c, a, a)
}

func lfqHalve(c, a *LazyFq) {
	mask := -(a[0] & 1)
	var carry uint64
	for i, ai := range Q {
		ai &= mask
		bi := c[i+6]
		ci := ai + bi + carry
		c[i+6] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
	for i := 6; i < 12; i++ {
		ai := a[i]
		ci := ai + carry
		c[i] = ci
		carry = (ai | ai&^ci) >> 63
	}
	for i := 11; i >= 0; i-- {
		nc := c[i] & 1
		c[i] >>= 1
		c[i] |= (carry<<63)
		carry = nc
	}
}

func fqExp(c, a *Fq, s []uint64) {
	nw := len(s)
	base := *a
	res := new(Fq).SetOne()
	for i := 0; i < nw; i++ {
		w := s[i]
		for j := uint(0); j < 64; j++ {
			if w & (1<<j) != 0 {
				fqMul(res, res, &base)
			}
			fqSqr(&base, &base)
		}
	}
	*c = *res
}
/*
func fqExp(c, a *Fq, s []uint64) {
	r := new(Fq).Set(a)
	sl := len(s)-1
	for i := sl * 64 + bits.Len64(s[sl]) - 2; i >= 0; i-- {
		fqSqr(r, r)
		if (s[i/64] & (1<<uint(i%64))) != 0 {
			fqMul(r, r, a)
		}
	}
	c.Set(r)
}*/

func fqInv(c, a *Fq) {
	fqExp(c, a, qm2[:])
}

func fqCmp(a, b *Fq) int {
	return fqCmpN([]*Fq{a}, []*Fq{b})
}

func fqCmpN(a,b []*Fq) int {
	var res uint64
	for j, ja := range a {
		jb := b[j]
		for i := 0; i < 6; i++ {
			res |= (ja[i] ^ jb[i])
		}
	}
	return int((((res & mask32) | (res >> 32)) - 1) >> 63)
}


func fqSwap(a, b []*Fq, swap uint64) {
	swap -= 1

	for i, ai := range a {
		bi := b[i]
		x0 := swap & (ai[0] ^ bi[0])
		ai[0] ^= x0
		bi[0] ^= x0
		x1 := swap & (ai[1] ^ bi[1])
		ai[1] ^= x1
		bi[1] ^= x1
		x2 := swap & (ai[2] ^ bi[2])
		ai[2] ^= x2
		bi[2] ^= x2
		x3 := swap & (ai[3] ^ bi[3])
		ai[3] ^= x3
		bi[3] ^= x3
		x4 := swap & (ai[4] ^ bi[4])
		ai[4] ^= x4
		bi[4] ^= x4
		x5 := swap & (ai[5] ^ bi[5])
		ai[5] ^= x5
		bi[5] ^= x5
	}
}

func lfqSub2(c, a, b *LazyFq) {
	var carry uint64
	for i, ai := range a {
		bi := b[i]
		ci := ai - bi - carry
		c[i] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
	mask := -carry
	carry = 0
	for i, ai := range Q {
		ai &= mask
		bi := c[i+6]
		ci := ai + bi + carry
		c[i+6] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
}


