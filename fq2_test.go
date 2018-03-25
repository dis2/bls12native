package bls12

import "testing"
//import "math/big"
//import "crypto/rand"

func TestFq2Mul100(t *testing.T) {
	var a, b Fq2
	for i, n := range loadTest("fq2mul100") {
		a.SetUInt64(uint64(i * 5))
		b.SetUInt64(uint64(i * 10))
		fq2Mul(&a, &a, &b)
		if a[0] != n[0] || a[1] != n[1] {
			t.Fatalf("bad mult at %d, %d(%d) != %d\n", i, a[0].ToInt(), a[1].ToInt(), n[0].ToInt())
		}
	}
}



