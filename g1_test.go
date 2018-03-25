package bls12

import "testing"
import "math/rand"
import "strings"
import "os"
import "bufio"


func loadTest(n string) (res [][]Fq) {
	f, _ := os.Open("testdata/"+n+".txt")
	r := bufio.NewReader(f)
	for {
		ln, err := r.ReadString('\n')
		if err != nil {
			break
		}
		var vv []Fq
		for _, v := range strings.Split(ln, " ") {
			vv = append(vv, parseFqMont(v))
		}
		res = append(res, vv)
	}
	return
}

func TestG1First100(t *testing.T) {
	var a, b G1
	a.SetOne()
	b.SetOne()

	for i, n := range loadTest("g1first100") {
		if !n[0].Equal(&a.X) || !n[1].Equal(&a.Y) {
			t.Fatalf("bad point at %d\n", i)
		}
		a.Add(&a,&b)
		a.ToAffine(&a)
	}
}

func BenchmarkG1Double(b *testing.B) {
	a := G1{Fq{1}, Fq{2}, Fq{}, false, false}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Double(nil)
	}
}

func BenchmarkG1Add(b *testing.B) {
	p := new(G1).G()
	q := new(G1).G().Double(nil)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.Add(nil, q)
	}
}

func BenchmarkG1Mul(b *testing.B) {
	p := new(G1).G()
	s := []uint64{rand.Uint64(), rand.Uint64(), rand.Uint64(), rand.Uint64()}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.Mul(nil, s)
	}
}
