// +build !amd64 generic noasm

package bls12

var lfqAdd = lfqAddGeneric
var lfqSub = lfqSubGeneric
var fqAdd = fqAddGeneric
var fqSub = fqSubGeneric
var fqNeg = fqNegGeneric
var fqREDC = fqREDCGeneric
var fqMul = fqMulGeneric
var lfqMul = lfqMulGeneric
var fqReduce = fqReduceGeneric
