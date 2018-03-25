package bls12

import "testing"
import "math/rand"


func TestG2First100(t *testing.T) {
	var a, b G2
	a.SetOne()
	b.SetOne()
	for i, n := range loadTest("g2first100") {
		if !n[0].Equal(&a.X[0]) || !n[1].Equal(&a.X[1]) || !n[2].Equal(&a.Y[0]) || !n[3].Equal(&a.Y[1]) {
			t.Fatalf("bad point at %d\n", i)
		}

		a.Add(&a, &b)
//		a.Double(&a)
		a.ToAffine(&a)
	}
}

func TestG2Double100(t *testing.T) {
	var a, b G2
	a.SetOne()
	b.SetOne()
	for i, n := range loadTest("g2dbl100") {
		if !n[0].Equal(&a.X[0]) || !n[1].Equal(&a.X[1]) || !n[2].Equal(&a.Y[0]) || !n[3].Equal(&a.Y[1]) {
			t.Fatalf("bad point at %d\n", i)
		}

		a.Double(&a)
		a.ToAffine(&a)
	}
}



func BenchmarkG2Double(b *testing.B) {
	a := new(G2).G()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Double(nil)
	}
}

func BenchmarkG2Add(b *testing.B) {
	p := new(G2).G()
	q := new(G2).G().Double(nil)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.Add(nil, q)
	}
}

func BenchmarkG2Mul(b *testing.B) {
	p := new(G2).G()
	s := []uint64{rand.Uint64(), rand.Uint64(), rand.Uint64(), rand.Uint64()}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.Mul(nil, s)
	}
}
