package bls12

import "testing"
import "math/big"
import "crypto/rand"

func fqReduceTest(a *Fq) {
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

func lfqAddTest(c, a, b *LazyFq) {
	var carry uint64
	for i, ai := range a {
		bi := b[i]
		ci := ai + bi + carry
		c[i] = ci
		carry = (ai&bi | (ai|bi)&^ci) >> 63
	}
}

func lfqSubTest(c, a, b *LazyFq) {
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

func fqREDCTest(cbuf *Fq, abuf *LazyFq) {
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
	fqReduceTest(cbuf)
}


func rnd() (f Fq, i *big.Int) {
	lim := new(big.Int).SetInt64(1024)
	lim = Q.RawInt()
	i, _ = rand.Int(rand.Reader, lim)
	f.FromInt(i)
	return
}

/*
func TestSubSimple(t *testing.T) {
	c,a,b := LazyFq{},LazyFq{30}, LazyFq{20}
	lfqSub(&c, &a, &b)
	t.Logf("%v\n", c)
}*/

func TestFqInv(t *testing.T) {
	for i := 0; i < 100; i++ {
		a, ia := rnd()
		fqInv(&a, &a)
		ia.ModInverse(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed little fermat")
		}
	}
}

func TestFqAdd(t *testing.T) {
	for i := 0; i < 10000; i++ {
		a, ia := rnd()
		b, ib := rnd()
		fqAdd(&a, &a, &b)
		ia.Add(ia, ib)
		ia.Mod(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed addition")
		}
	}
}

func TestFqSub(t *testing.T) {
	for i := 0; i < 10000; i++ {
		a, ia := rnd()
		b, ib := rnd()
		fqSub(&a, &a, &b)
		ia.Sub(ia, ib)
		ia.Mod(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed substraction")
		}
	}
}

func TestLazyFqSub(t *testing.T) {
	a, ia := rnd()
	b, ib := rnd()
	c, ic := rnd()
	var la,lb LazyFq
	lfqMul(&la, &a, &b)
	lfqMul(&lb, &a, &c)
	lia := new(big.Int).Mul(ia, ib)
	lib := new(big.Int).Mul(ia, ic)

	fqREDC(&a, la.Copy())
	if a.ToInt().Cmp(lia.Mod(lia, Q.RawInt())) != 0 {
		t.Fatal("lazy prep failed")
	}

	fqREDC(&b, lb.Copy())
	if b.ToInt().Cmp(lib.Mod(lib, Q.RawInt())) != 0 {
		t.Fatal("lazy prep failed")
	}

	// lia == la
	// lib == lb
	for i := 0; i < 10000; i++ {
		lfqSubTest(&la, &la, &lb)
		lia.Sub(lia, lib)
		lia.Mod(lia, Q.RawInt())
		if i % 10 == 0 {
			fqREDC(&a, la.Copy())
			if a.ToInt().Cmp(lia) != 0 {
				t.Fatalf("failed lazy substraction %d %d", a.ToInt(), lia)
			}
		}
	}
}


func TestLazyFqAdd(t *testing.T) {
	a, ia := rnd()
	b, ib := rnd()
	c, ic := rnd()
	var la,lb LazyFq
	lfqMul(&la, &a, &b)
	lfqMul(&lb, &a, &c)
	lia := new(big.Int).Mul(ia, ib)
	lib := new(big.Int).Mul(ia, ic)

	fqREDC(&a, la.Copy())
	if a.ToInt().Cmp(lia.Mod(lia, Q.RawInt())) != 0 {
		t.Fatal("lazy prep failed")
	}

	fqREDC(&b, lb.Copy())
	if b.ToInt().Cmp(lib.Mod(lib, Q.RawInt())) != 0 {
		t.Fatal("lazy prep failed")
	}

	// lia == la
	// lib == lb
	// Add does not do long reduces
	// It is safe to do only 2^384 / Q (~8) additions before run out of
	for i := 0; i < 10000; i++ {
		lfqAdd(&la, &la, &lb)
		lia.Add(lia, lib)
		lia.Mod(lia, Q.RawInt())
		if i % 100 == 0 {
			fqREDC(&a, la.Copy())
			if a.ToInt().Cmp(lia) != 0 {
				t.Fatalf("failed lazy addition %d %d", a.ToInt(), lia)
			}
		}
	}
}


func TestFqNeg(t *testing.T) {
	for i := 0; i < 10000; i++ {
		a, ia := rnd()
		fqNeg(&a, &a)
		ia.Neg(ia)
		ia.Mod(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed neg")
		}
	}
}

func TestFqMul(t *testing.T) {

	for i := 0; i < 10000; i++ {
		a, ia := rnd()
		b, ib := rnd()
		fqMul(&a, &a, &b)
		ia.Mul(ia, ib)
		ia.Mod(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed mult")
		}
	}
}

func TestLazyFqMulREDC(t *testing.T) {

	for i := 0; i < 10000; i++ {
		a, ia := rnd()
		b, ib := rnd()
		var a2 LazyFq
		lfqMul(&a2, &a, &b)
		fqREDC(&a, &a2)
		ia.Mul(ia, ib)
		ia.Mod(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed mult")
		}
	}
}

func TestFqSquare(t *testing.T) {
	for i := 0; i < 10000; i++ {
		a, ia := rnd()
		fqSqr(&a, &a)
		ia.Mul(ia, ia)
		ia.Mod(ia, Q.RawInt())
		if a.ToInt().Cmp(ia) != 0 {
			t.Fatal("failed square")
		}
	}
}

func BenchmarkFqAdd(b *testing.B) {
	fa, _ := rnd()
	fb, _ := rnd()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fqAdd(&fa, &fa, &fb)
	}
}

func BenchmarkFqLazyAdd(b *testing.B) {
	fc := LazyFq{1}
	fa := LazyFq{0, 0, 0, 0, 0, 2}
	fb := LazyFq{0, 0, 0, 0, 0, 0, 0, 0, 05}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		lfqAdd(&fc, &fa, &fb)
	}
}

func BenchmarkFqLazySub(b *testing.B) {
	fc := LazyFq{1}
	fa := LazyFq{0, 0, 0, 0, 0, 2}
	fb := LazyFq{0, 0, 0, 0, 0, 0, 0, 0, 05}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		lfqSub(&fc, &fa, &fb)
	}
}


func BenchmarkFqSub(b *testing.B) {
	fa, _ := rnd()
	fb, _ := rnd()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fqSub(&fa, &fa, &fb)
	}
}

func BenchmarkFqNeg(b *testing.B) {
	a, _ := rnd()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fqNeg(&a, &a)
	}
}

func BenchmarkFqMul(b *testing.B) {
	fa, _ := rnd()
	fb, _ := rnd()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fqMul(&fa, &fa, &fb)
	}
}

func BenchmarkFqLazyMul(b *testing.B) {
	a, x := Fq{1, 2}, Fq{3, 4}
	var c LazyFq
	for i := 0; i < b.N; i++ {
		lfqMul(&c, &a, &x)
	}
}

func BenchmarkFqREDC(b *testing.B) {
	a, x := Fq{1, 2}, Fq{3, 4}
	var c LazyFq
	lfqMul(&c, &a, &x)
	for i := 0; i < b.N; i++ {
		fqREDC(&a, &c)
	}
}

func BenchmarkFqReduce(b *testing.B) {
	a := Fq{1, 2}
	for i := 0; i < b.N; i++ {
		fqReduce(&a)
	}
}

func BenchmarkFqInv(b *testing.B) {
	a := Fq{1, 2}
	for i := 0; i < b.N; i++ {
		fqInv(&a, &a)
	}
}
