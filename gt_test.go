package bls12

import "testing"

func BenchmarkGTPair(b *testing.B) {
	p := new(G1).G()
	q := new(G2).G()
	gt := new(GT)
	for i := 0; i < b.N; i++ {
		gt.Pair(p,q)
	}
}
