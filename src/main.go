// main
package main

import (
	"fmt"
	"quantum"
)

var r1 quantum.QuantumReg

func main() {
	r1 = quantum.QuantumNewQureg(0, 1)
	quantum.QuantumHadamard(0, &r1)
	quantum.QuantumPrintQureg(r1)
	fmt.Println(quantum.QuantumBMeasure(0, &r1))
}
