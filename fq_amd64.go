// +build !noasm,!generic,amd64

package bls12

//go:noescape
func fqMul(c, a, b *Fq)

//go:noescape
func lfqMul(c *LazyFq, a, b *Fq)

//go:noescape
func fqAdd(c, a, b *Fq)

//go:noescape
func lfqAdd(c, a, b *LazyFq)

//go:noescape
func lfqSub(c, a, b *LazyFq)

//go:noescape
func fqREDC(c *Fq, a *LazyFq)

//go:noescape
func fqReduce(c *Fq)

//go:noescape
func fqSub(c, a, b *Fq)

//go:noescape
func fqNeg(c, a *Fq)
