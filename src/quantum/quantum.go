package quantum

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

type quantum_error struct {
	When time.Time
	What string
}

var QuantumStatus int = 0
var QuantumLambda float32 = 0.0
var QECType int = 0
var QECWidth int = 0

const numRegs = 4

func QuantumQECSetStatus(stype int, swidth int) {
	QECType = stype
	QECWidth = swidth
}

func QuantumQECGetStatus(ptype *int, pwidth *int) {
	if ptype != nil {
		*ptype = QECType
	}
	if pwidth != nil {
		*pwidth = QECWidth
	}
}

func QuantumQECEncode(type1 int, width1 int, reg *QuantumReg) {
	var lambda float32

	lambda = QuantumGetDecoherence()

	QuantumSetDecoherence(0)

	for i := 0; i < reg.width; i++ {
		if i == reg.width-1 {
			QuantumSetDecoherence(lambda)
		}
		if i < width1 {
			QuantumHadamard(reg.width+i, reg)
			QuantumHadamard(2*reg.width+i, reg)
			QuantumCNOT(reg.width+i, i, reg)
			QuantumCNOT(2*reg.width+i, i, reg)
		} else {
			QuantumCNOT(i, reg.width+i, reg)
			QuantumCNOT(i, 2*reg.width+i, reg)
		}
	}
	QuantumQECSetStatus(1, reg.width)
	reg.width *= 3
}

func QuantumQECDecode(type1 int, width1 int, reg *QuantumReg) {
	var swidth, a, b int
	var lambda float32

	lambda = QuantumGetDecoherence()

	QuantumSetDecoherence(0)

	swidth = reg.width / 3

	QuantumQECSetStatus(0, 0)

	for i := reg.width/3 - 1; i >= 0; i-- {
		if i == 0 {
			QuantumSetDecoherence(lambda)
		}
		if i < width1 {
			QuantumCNOT(2*swidth+i, i, reg)
			QuantumCNOT(swidth+i, i, reg)
			QuantumHadamard(2*swidth+i, reg)
			QuantumHadamard(swidth+i, reg)
		} else {
			QuantumCNOT(i, 2*swidth+i, reg)
			QuantumCNOT(i, swidth+i, reg)
		}
	}
	for i := 1; i <= swidth; i++ {
		a = QuantumBMeasure(swidth, reg)
		b = QuantumBMeasure(2*swidth-i, reg)

		if a == 1 && b == 1 && i-1 < QECWidth {
			QuantumSigmaZ(i-1, reg) // Z = HXH
		}
	}
}

func QuantumQECCounter(inc int, frequency int, reg *QuantumReg) int {
	var counter int = 0
	var freq int = 1 << 30

	if inc > 0 {
		counter += inc
	} else if inc < 0 {
		counter = 0
	}

	if frequency > 0 {
		freq = frequency
	}

	if counter >= freq {
		counter = 0
		QuantumQECDecode(QECType, QECWidth, reg)
		QuantumQECEncode(QECType, QECWidth, reg)
	}
	return counter
}

func QuantumCNOTft(control int, target int, reg *QuantumReg) {
	var lambda float32

	tmp := QECType
	QECType = 0

	lambda = QuantumGetDecoherence()
	QuantumSetDecoherence(0)

	QuantumCNOT(control, target, reg)
	QuantumCNOT(control+QECWidth, target+QECWidth, reg)
	QuantumSetDecoherence(lambda)
	QuantumCNOT(control+2*QECWidth, target+2*QECWidth, reg)

	QuantumQECCounter(1, 0, reg)

	QECType = tmp
}

func QuantumSigmaXft(target int, reg *QuantumReg) {
	var lambda float32
	var tmp int

	tmp = QECType
	QECType = 0

	lambda = QuantumGetDecoherence()
	QuantumSetDecoherence(0)

	QuantumSigmaX(target, reg)
	QuantumSigmaX(target+QECWidth, reg)
	QuantumSetDecoherence(lambda)
	QuantumSigmaX(target+2*QECWidth, reg)

	QuantumQECCounter(1, 0, reg)

	QECType = tmp
}

func QuantumToffolift(control1 int, control2 int, target int, reg *QuantumReg) {
	var c1, c2 int
	var mask uint32

	mask = (1 << uint(target)) +
		(1 << uint(target+QECWidth)) +
		(1 << uint(target+2*QECWidth))

	for i := 0; i < reg.size; i++ {
		c1 = 0
		c2 = 0
		if reg.node[i].state&(1<<uint(control1)) != 0 {
			c1 = 1
		}
		if reg.node[i].state&(1<<uint(control1+QECWidth)) != 0 {
			c1 ^= 1
		}
		if reg.node[i].state&(1<<uint(control1+2*QECWidth)) != 0 {
			c1 ^= 1
		}
		if reg.node[i].state&(1<<uint(control2)) != 0 {
			c2 = 1
		}
		if reg.node[i].state&(1<<uint(control2+QECWidth)) != 0 {
			c2 ^= 1
		}
		if reg.node[i].state&(1<<uint(control2+2*QECWidth)) != 0 {
			c2 ^= 1
		}
		if c1 == 1 && c2 == 1 {
			reg.node[i].state = reg.node[i].state ^ mask
		}
	}
	QuantumDecohere(reg)

	QuantumQECCounter(1, 0, reg)
}

func QuantumGetDecoherence() float32 {
	return QuantumLambda
}

func QuantumSetDecoherence(l float32) {
	if l != 0 {
		QuantumStatus = 1
		QuantumLambda = l
	} else {
		QuantumStatus = 0
	}
}

//This is a semi-quantum fulladder. It adds to b_in
//a c-number. Carry-in bit is c_in and carry_out is
//c_out. xlt-l and L are enablebits. See documentation
//for further information
func Muxfa(a int, bIn int, cIn int, cOut int, xltL int, L int, total int, reg *QuantumReg) {
	if a == 0 {
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumCNOT(bIn, cIn, reg)
	}

	if a == 3 {
		QuantumToffoli(L, cIn, cOut, reg)
		QuantumCNOT(L, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumCNOT(bIn, cIn, reg)
	}

	if a == 1 {
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumCNOT(bIn, cIn, reg)
	}

	if a == 2 {
		QuantumSigmaX(xltL, reg)
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumCNOT(bIn, cIn, reg)
		QuantumSigmaX(xltL, reg)
	}
}

func MuxfaInv(a int, bIn int, cIn int, cOut int, xltL int, L int, total int, reg *QuantumReg) {
	if a == 0 {
		QuantumCNOT(bIn, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
	}

	if a == 3 {
		QuantumCNOT(bIn, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumCNOT(L, cIn, reg)
		QuantumToffoli(L, cIn, cOut, reg)
	}

	if a == 1 {
		QuantumCNOT(bIn, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, bIn, reg)
	}

	if a == 2 {
		QuantumSigmaX(xltL, reg)
		QuantumCNOT(bIn, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumToffoli(bIn, cIn, cOut, reg)
		QuantumToffoli(L, xltL, bIn, reg)
		QuantumSigmaX(xltL, reg)
	}
}

//This is a semi-quantum halfadder. It adds to b_in
//a c-number. Carry-in bit is c_in and carry_out is
//not necessary. xlt-l and L are enablebits. See
//documentation for further information
func Muxha(a int, bIn int, cIn int, xltL int, L int, total int, reg *QuantumReg) {
	if a == 0 {
		QuantumCNOT(bIn, cIn, reg)
	}
	if a == 3 {
		QuantumCNOT(L, cIn, reg)
		QuantumCNOT(bIn, cIn, reg)
	}

	if a == 1 {
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumCNOT(bIn, cIn, reg)
	}

	if a == 2 {
		QuantumSigmaX(xltL, reg)
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumCNOT(bIn, cIn, reg)
		QuantumSigmaX(xltL, reg)
	}
}

func MuxhaInv(a int, bIn int, cIn int, xltL int, L int, total int, reg *QuantumReg) {
	if a == 0 {
		QuantumCNOT(bIn, cIn, reg)
	}

	if a == 3 {
		QuantumCNOT(bIn, cIn, reg)
		QuantumCNOT(L, cIn, reg)
	}

	if a == 1 {
		QuantumCNOT(bIn, cIn, reg)
		QuantumToffoli(L, xltL, cIn, reg)
	}

	if a == 2 {
		QuantumSigmaX(xltL, reg)
		QuantumCNOT(bIn, cIn, reg)
		QuantumToffoli(L, xltL, cIn, reg)
		QuantumSigmaX(xltL, reg)
	}
}

func Madd(a int, aInv int, width int, reg *QuantumReg) {
	var total, j int

	total = numRegs*width + 2
	for i := 0; i < width-1; i++ {
		if (1<<uint(i))&a != 0 {
			j = 1 << 1
		} else {
			j = 0
		}
		if (1<<uint(i))&aInv != 0 {
			j += 1
		}
		Muxfa(j, width+i, i, i+1, 2*width, 2*width+1, total, reg)
	}
	j = 0
	if (1<<uint(width-1))&a != 0 {
		j = 2
	}
	if (1<<uint(width-1))&aInv != 0 {
		j += 1
	}
	Muxha(j, 2*width-1, width-1, 2*width, 2*width+1, total, reg)
}

func MaddInv(a int, aInv int, width int, reg *QuantumReg) {
	var total, j int

	total = numRegs*width + 2
	j = 0

	if 1<<uint(width-1)&a != 0 {
		j = 2
	}
	if 1<<uint(width-1)&aInv != 0 {
		j += 1
	}

	MuxhaInv(j, width-1, 2*width-1, 2*width, 2*width+1, total, reg)

	for i := width - 2; i >= 0; i-- {
		if (1<<uint(i))&a != 0 {
			j = 1 << 1
		} else {
			j = 0
		}
		if (1<<uint(i))&aInv != 0 {
			j += 1
		}
		MuxfaInv(j, i, width+i, width+1+i, 2*width, 2*width+1, total, reg)
	}
}

func Addn(N int, a int, width int, reg *QuantumReg) {
	testSum(N-a, width, reg)
	Madd((1<<uint(width))+a-N, a, width, reg)
}

func AddnInv(N int, a int, width int, reg *QuantumReg) {
	QuantumCNOT(2*width+1, 2*width, reg)
	MaddInv(1<<uint(width)-a, N-a, width, reg)

	QuantumSwapTheLeads(width, reg)

	testSum(a, width, reg)
}

func AddModN(N int, a int, width int, reg *QuantumReg) {
	Addn(N, a, width, reg)
	AddnInv(N, a, width, reg)
}

func testSum(compare int, width int, reg *QuantumReg) {
	if compare&(1<<uint(width-1)) != 0 {
		QuantumCNOT(2*width-1, width-1, reg)
		QuantumSigmaX(2*width-1, reg)
		QuantumCNOT(2*width-1, 0, reg)
	} else {
		QuantumSigmaX(2*width-1, reg)
		QuantumCNOT(2*width-1, width-1, reg)
	}
	for i := width - 2; i > 0; i-- {
		if compare&1<<uint(i) != 0 {
			QuantumToffoli(i+1, width+i, i, reg)
			QuantumSigmaX(width+i, reg)
			QuantumToffoli(i+1, width+i, 0, reg)
		} else {
			QuantumSigmaX(width+i, reg)
			QuantumToffoli(i+1, width+i, i, reg)
		}
	}
	if compare&1 != 0 {
		QuantumSigmaX(width, reg)
		QuantumToffoli(width, 1, 0, reg)
	}
	QuantumToffoli(2*width+1, 0, 2*width, reg)

	if compare&1 != 0 {
		QuantumToffoli(width, 1, 0, reg)
		QuantumSigmaX(width, reg)
	}

	for i := 1; i <= width-2; i++ {
		if compare&1<<uint(i) != 0 {
			QuantumToffoli(i+1, width+i, 0, reg)
			QuantumSigmaX(width+i, reg)
			QuantumToffoli(i+1, width+i, i, reg)
		} else {
			QuantumToffoli(i+1, width+i, i, reg)
			QuantumSigmaX(width+i, reg)
		}
	}
	if compare&1<<uint(width-1) != 0 {
		QuantumCNOT(2*width-1, 0, reg)
		QuantumSigmaX(2*width-1, reg)
		QuantumCNOT(2*width-1, width-1, reg)
	} else {
		QuantumCNOT(2*width-1, width-1, reg)
		QuantumSigmaX(2*width-1, reg)
	}
}

/* Calculates the number of qubits required to store N */
func QuantumGetwidth(n int) int {
	var i int

	for i = 1; (1 << uint(i)) < n; i++ {
	}

	return i
}

/* Calculate the inverse modulus of N and C */
func QuantumInverseMod(n int, c int) int {
	var i int

	for i = 1; (i*c)%n != 1; i++ {
	}

	return i
}

/* Fractional approximation of a decimal value */
func QuantumFracApprox(a *int, b *int, width int) {

	var f, g float64
	var i int
	var num2, den2, num1, den1, num, den int = 0, 1, 1, 0, 0, 0

	f = float64(*a / *b)

	i = int(g + 5e-6)

	g = f

	g -= float64(i) - float64(5e-6)
	g = 1.0 / g

	tt := 1 << uint(width)
	if i*den1+den2 > int(tt) {
		return
	}

	num = i*num1 + num2
	den = i*den1 + den2

	num2 = num1
	den2 = den1
	num1 = num
	den1 = den

	for math.Abs(float64(num)/float64(den))-f > 1.0/(2*float64(tt)) {
		i = int(g + 0.000005)

		g -= float64(i) - float64(5e-6)
		g = 1.0 / g

		if i*den1+den2 > int(tt) {
			break
		}

		num = i*num1 + num2
		den = i*den1 + den2

		num2 = num1
		den2 = den1
		num1 = num
		den1 = den
	}

	*a = num
	*b = den
}

/* Calculate the greatest common divisor with Euclid's algorithm */
func QuantumGCD(u int, v int) int {
	var r int

	for v != 0 {
		r = u % v
		u = v
		v = r
	}

	return u
}

func QuantumMeasure(reg QuantumReg) int {
	var r float32
	r141272 := rand.New(rand.NewSource(time.Now().UnixNano()))

	r = r141272.Float32()

	for i := 0; i < reg.size; i++ {
		r -= quantum_prob(reg.node[i].amplitude)
		if 0 >= r {
			return int(reg.node[i].state)
		}
	}

	return -1
}

/* Measure a single bit of a quantum register. The bit measured is
   indicated by its position POS, starting with 0 as the least
   significant bit. The new state of the quantum register depends on
   the result of the measurement. */
func QuantumBMeasure(pos int, reg *QuantumReg) int {
	var result int = 0
	var pa float32 = 0.0
	var r float32
	var pos2 uint32
	var out QuantumReg
	r141272 := rand.New(rand.NewSource(time.Now().UnixNano()))

	pos2 = 1 << uint(pos)

	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&pos2 == 0 {
			pa += quantum_prob(reg.node[i].amplitude)
		}
	}

	r = r141272.Float32()

	if r > pa {
		result = 1
	}

	out = QuantumStateCollapse(pos, result, *reg)

	*reg = out

	return result
}

func QuantumDecohere(reg *QuantumReg) {
	var u, v, s, x, angle float64
	r141272 := rand.New(rand.NewSource(time.Now().UnixNano()))

	QuantumGateCounter(1)

	if QuantumStatus != 0 {
		nrands := make([]float64, reg.width)
		s = 1.0 // starting point. value not used.
		for i := 0; i < reg.width; i++ {
			for s >= 1.0 {
				u = 2*r141272.Float64() - 1
				v = 2*r141272.Float64() - 1
				s = u*u + v*v
			}
			x = u * math.Sqrt(-2*math.Log(s)/s)
			x *= math.Sqrt(2 * float64(QuantumLambda))
			nrands[i] = x / 2
		}
		for i := 0; i < reg.size; i++ {
			angle = 0.0
			for j := 0; j < reg.width; j++ {
				if reg.node[i].state&(1<<uint(j)) != 0 {
					angle += nrands[j]
				} else {
					angle -= nrands[j]
				}
			}
			reg.node[i].amplitude *= quantum_cexp(float32(angle))
		}
	}
}

func btoi(b bool) int {
	if b {
		return 1
	}
	return 0
}

func itob(i int) bool {
	return i != 0
}

func (e *quantum_error) Error() string {
	return fmt.Sprintf("at %v, %s", e.When, e.What)
}

type quantum_reg_node struct {
	amplitude complex64
	state     uint32
}

type QuantumReg struct {
	width, size, hashw int
	node               []quantum_reg_node
	hash               []int
}

type QuantumMatrix struct {
	rows, cols int
	t          []complex64
}

func MM(m QuantumMatrix, x int, y int) complex64 {
	return m.t[x+(y*m.cols)]
}

func quantum_cexp(a float32) complex64 {
	return complex(float32(math.Cos(float64(a))), float32(math.Sin(float64(a))))
}

func QuantumNewMatrix(cols int, rows int) QuantumMatrix {
	var m QuantumMatrix

	m.rows = rows
	m.cols = cols

	m.t = make([]complex64, cols*rows)

	return m
}

func QuantumMMult(A QuantumMatrix, B QuantumMatrix) (QuantumMatrix, error) {
	var C QuantumMatrix

	if A.cols != B.rows {
		return C, &quantum_error{
			time.Now(),
			"Error in matrix size",
		}
	}

	C = QuantumNewMatrix(B.cols, A.rows)

	for i := 0; i < B.cols; i++ {
		for j := 0; j < A.rows; j++ {
			for k := 0; k < B.rows; k++ {
				C.t[i+j*C.cols] += MM(A, k, j) * MM(B, i, k)
			}
		}
	}

	return C, nil
}

func quantum_hash64(key uint32, width int) uint {
	var k32 uint32

	k32 = (key & 0xFFFFFFFF) ^ (key >> 32)
	k32 *= 0x9e370001
	k32 = uint32(k32) >> (32 - uint32(width))
	return uint(k32)
}

func Quantum_get_state(a uint32, reg QuantumReg) int {
	var i uint

	if reg.hashw == 0 {
		return int(a)
	}

	i = quantum_hash64(a, reg.hashw)

	for reg.hash[i] != 0 {
		if reg.node[reg.hash[i]-1].state == a {
			return reg.hash[i] - 1
		}
		i++
		if i == (1 << uint(reg.hashw)) {
			i = 0
		}
	}

	return -1
}

func quantum_add_hash(a uint32, pos int, reg *QuantumReg) error {
	i := quantum_hash64(a, reg.hashw)

	mark := 0

	for reg.hash[i] != 0 {
		i += 1
		if i == (1 << uint(reg.hashw)) {
			if mark == 0 {
				i = 0
				mark = 1
			} else {
				return &quantum_error{
					time.Now(),
					"Hash Full",
				}
			}
		}
	}

	reg.hash[i] = pos + 1

	return nil
}

func quantum_reconstruct_hash(reg *QuantumReg) {
	if reg.hashw == 0 {
		return
	}

	for i := 0; i < (1 << uint(reg.hashw)); i++ {
		reg.hash[i] = 0
	}

	for i := 0; i < reg.size; i++ {
		quantum_add_hash(reg.node[i].state, i, reg)
	}
}

func QuantumNewQureg(initval uint32, width int) QuantumReg {
	var reg QuantumReg

	reg.width = width
	reg.size = 1
	reg.hashw = width + 2

	// Allocate memory for 1 base state

	reg.node = make([]quantum_reg_node, 1)
	reg.hash = make([]int, (1 << uint(reg.hashw)))

	reg.node[0].state = initval
	reg.node[0].amplitude = 1

	return reg
}

func quantum_prob(a complex64) float32 {
	return real(a)*real(a) + imag(a)*imag(a)
}

func quantum_destroy_hash(reg *QuantumReg) {
	reg.hash = make([]int, 0)
}

func quantum_delete_qureg(reg *QuantumReg) {
	if reg.hashw != 0 {
		quantum_destroy_hash(reg)
	}
	reg.node = make([]quantum_reg_node, 0)
}

func quantum_delete_qureg_hashpreserve(reg *QuantumReg) {
	reg.node = make([]quantum_reg_node, 0)
}

func QuantumPrintQureg(reg QuantumReg) {
	for i := 0; i < reg.size; i++ {
		fmt.Printf(" %v %vi|%v> (%v) (|", real(reg.node[i].amplitude),
			imag(reg.node[i].amplitude),
			reg.node[i].state, quantum_prob(reg.node[i].amplitude))
		for j := reg.width - 1; j >= 0; j-- {
			if j%4 == 3 {
				fmt.Printf(" ")
			}
			fmt.Printf("%v", btoi((((1 << uint(j)) & reg.node[i].state) > 0)))
		}
		fmt.Printf(">)\n")
	}
	fmt.Printf("\n")
}

func QuantumPrintExpn(reg QuantumReg) {
	for i := 0; i < reg.size; i++ {
		fmt.Printf("%v: %v\n", i, uint(reg.node[i].state)-uint(i)*(1<<uint(reg.width/2)))
	}
}

func QuantumAddscratch(bits int, reg *QuantumReg) {
	var l uint32
	reg.width += bits

	for i := 0; i < reg.size; i++ {
		l = reg.node[i].state << uint(bits)
		reg.node[i].state = l
	}
}

func QuantumStateCollapse(pos int, value int, reg QuantumReg) QuantumReg {
	var size int = 0
	var d float64 = 0.0
	var lpat, rpat, pos2 uint32
	var out QuantumReg

	pos2 = 1 << uint(pos)

	for i := 0; i < reg.size; i++ {
		if (itob(int(reg.node[i].state&pos2)) && itob(value)) ||
			(!(itob(int(reg.node[i].state & pos2))) && !itob(value)) {
			d += float64(quantum_prob(reg.node[i].amplitude))
			size++
		}
	}

	out.width = reg.width - 1
	out.size = size
	out.node = make([]quantum_reg_node, size)
	out.hashw = reg.hashw
	out.hash = reg.hash

	j := 0

	for i := 0; i < reg.size; i++ {
		if (itob(int(reg.node[i].state&pos2)) && itob(value)) ||
			(!(itob(int(reg.node[i].state & pos2))) && !itob(value)) {
			rpat = 0
			for k := 0; k < pos; k++ {
				rpat += 1 << uint(k)
			}

			rpat &= reg.node[i].state

			lpat = 0
			for k := 32*8 - 1; k > pos; k-- {
				lpat += 1 << uint(k)
			}

			lpat &= reg.node[i].state

			out.node[j].state = (lpat >> 1) | rpat
			out.node[j].amplitude = reg.node[i].amplitude * complex(1/float32(math.Sqrt(d)), 0)

			j++
		}
	}

	return out
}

func QuantumPrintHash(reg QuantumReg) {
	for i := 0; i < (1 << uint(reg.hashw)); i++ {
		if i != 0 {
			fmt.Printf("%v: %v %v\n", i, reg.hash[i]-1,
				reg.node[reg.hash[i]-1].state)
		}
	}
}

func QuantumGate1(target int, m QuantumMatrix, reg *QuantumReg) error {
	var iset int
	var j, k int = 0, 0
	var addsize, decsize int = 0, 0
	var t1 complex64
	var tnot complex64 = 0
	var l1 uint32
	var limit float32

	if (m.cols != 2) || (m.rows != 2) {
		return &quantum_error{
			time.Now(),
			"Error in the matrix size",
		}
	}

	if reg.hashw != 0 {
		for i := 0; i < reg.size; i++ {
			if Quantum_get_state(reg.node[i].state^(1<<uint(target)), *reg) == -1 {
				addsize++
			}
		}

		tt := make([]quantum_reg_node, reg.size+addsize)
		copy(tt, reg.node)
		reg.node = tt

		for i := 0; i < addsize; i++ {
			reg.node[i+reg.size].state = uint32(0)
			reg.node[i+reg.size].amplitude = 0
		}
	}

	done := make([]int, reg.size+addsize)

	k = reg.size
	l1 = 1 << uint(reg.width)
	limit = (1.0 / float32(l1)) * math.SmallestNonzeroFloat32

	for i := 0; i < reg.size; i++ {
		if done[i] == 0 {
			iset = int(reg.node[i].state & (1 << uint(target)))
			tnot = 0
			j = Quantum_get_state(reg.node[i].state^(1<<uint(target)), *reg)
			t1 = reg.node[i].amplitude
			if j >= 0 {
				tnot = reg.node[j].amplitude
			}

			if iset != 0 {
				reg.node[i].amplitude = m.t[2]*tnot + m.t[3]*t1
			} else {
				reg.node[i].amplitude = m.t[0]*t1 + m.t[1]*tnot
			}

			if j >= 0 {
				if iset != 0 {
					reg.node[j].amplitude = m.t[0]*tnot + m.t[1]*t1
				} else {
					reg.node[j].amplitude = m.t[2]*t1 + m.t[3]*tnot
				}
			} else {
				if m.t[1] == 0 && itob(iset) {
					break
				}
				if m.t[2] == 0 && !itob(iset) {
					break
				}

				reg.node[k].state = reg.node[i].state ^ (1 << uint(target))

				if iset != 0 {
					reg.node[k].amplitude = m.t[1] * t1
				} else {
					reg.node[k].amplitude = m.t[2] * t1
				}

				k++
			}
			if j >= 0 {
				done[j] = 1
			}
		}
	}
	reg.size += addsize

	if reg.hashw != 0 {
		j := 0
		for i := 0; i < reg.size; i++ {
			if quantum_prob(reg.node[i].amplitude) < limit {
				j++
				decsize++
			} else if j != 0 {
				reg.node[i-j].state = reg.node[i].state
				reg.node[i-j].amplitude = reg.node[i].amplitude
			}
		}

		if decsize != 0 {
			reg.size -= decsize
			tt := make([]quantum_reg_node, reg.size)
			copy(tt, reg.node)
			reg.node = tt
		}
	}
	QuantumDecohere(reg)
	return nil
}

func QuantumHadamard(target int, reg *QuantumReg) {
	var m QuantumMatrix

	m = QuantumNewMatrix(2, 2)

	m.t[0] = complex(float32(math.Sqrt(1.0/2.0)), 0)
	m.t[1] = complex(float32(math.Sqrt(1.0/2.0)), 0)
	m.t[2] = complex(float32(math.Sqrt(1.0/2.0)), 0)
	m.t[3] = complex(float32(-math.Sqrt(1.0/2.0)), 0)

	QuantumGate1(target, m, reg)
}

func QuantumWalsh(width int, reg *QuantumReg) {
	for i := 0; i < width; i++ {
		QuantumHadamard(i, reg)
	}
}

func QuantumRX(target int, gamma float32, reg *QuantumReg) {
	var m QuantumMatrix

	m = QuantumNewMatrix(2, 2)

	m.t[0] = complex(float32(math.Cos(float64(gamma)/2)), 0)
	m.t[1] = complex(0, float32(-math.Sin(float64(gamma)/2)))
	m.t[2] = complex(0, float32(-math.Sin(float64(gamma)/2)))
	m.t[3] = complex(float32(math.Cos(float64(gamma)/2)), 0)

	QuantumGate1(target, m, reg)
}

func QuantumRY(target int, gamma float32, reg *QuantumReg) {
	var m QuantumMatrix

	m = QuantumNewMatrix(2, 2)

	m.t[0] = complex(float32(math.Cos(float64(gamma)/2)), 0)
	m.t[1] = complex(float32(-math.Sin(float64(gamma)/2)), 0)
	m.t[2] = complex(float32(math.Sin(float64(gamma)/2)), 0)
	m.t[3] = complex(float32(math.Cos(float64(gamma)/2)), 0)

	QuantumGate1(target, m, reg)
}

func QuantumRZ(target int, gamma float32, reg *QuantumReg) {
	z := quantum_cexp(gamma / 2)

	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&(1<<uint(target)) != 0 {
			reg.node[i].amplitude *= z
		} else {
			reg.node[i].amplitude /= z
		}
	}

	QuantumDecohere(reg)
}

func QuantumPhaseScale(target int, gamma float32, reg *QuantumReg) {
	z := quantum_cexp(gamma / 2)

	for i := 0; i < reg.size; i++ {
		reg.node[i].amplitude *= z
	}

	QuantumDecohere(reg)
}

func QuantumPhaseKick(target int, gamma float32, reg *QuantumReg) {
	z := quantum_cexp(gamma / 2)

	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&(1<<uint(target)) != 0 {
			reg.node[i].amplitude *= z
		}
	}

	QuantumDecohere(reg)
}

func QuantumCondPhase(control int, target int, reg *QuantumReg) {
	tt := (1 << uint(control-target))
	z := quantum_cexp(math.Pi / float32(tt))

	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&(1<<uint(control)) != 0 {
			if reg.node[i].state&(1<<uint(target)) != 0 {
				reg.node[i].amplitude *= z
			}
		}
	}

	QuantumDecohere(reg)
}

func QuantumCondPhaseInv(control int, target int, reg *QuantumReg) {
	tt := (1 << uint(control-target))
	z := quantum_cexp(-math.Pi / float32(tt))

	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&(1<<uint(control)) != 0 {
			if reg.node[i].state&(1<<uint(target)) != 0 {
				reg.node[i].amplitude *= z
			}
		}
	}

	QuantumDecohere(reg)
}

func QuantumCondPhaseKick(control int, target int, gamma float32,
	reg *QuantumReg) {
	z := quantum_cexp(gamma)

	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&(1<<uint(control)) != 0 {
			if reg.node[i].state&(1<<uint(target)) != 0 {
				reg.node[i].amplitude *= z
			}
		}
	}

	QuantumDecohere(reg)
}

func QuantumCNOT(control int, target int, reg *QuantumReg) {
	var qec int

	QuantumQECGetStatus(&qec, nil)

	if qec != 0 {
		QuantumCNOTft(control, target, reg)
	} else {
		for i := 0; i < reg.size; i++ {
			if reg.node[i].state&(1<<uint(control)) != 0 {
				reg.node[i].state ^= 1 << uint(target)
			}
		}
		QuantumDecohere(reg)
	}
}

func QuantumSigmaZ(target int, reg *QuantumReg) {
	for i := 0; i < reg.size; i++ {
		if reg.node[i].state&(1<<uint(target)) != 0 {
			reg.node[i].amplitude *= complex(-1, 0)
		}
	}
	QuantumDecohere(reg)
}

func QuantumToffoli(control1 int, control2 int, target int, reg *QuantumReg) {
	var qec int

	QuantumQECGetStatus(&qec, nil)

	if qec != 0 {
		QuantumToffolift(control1, control2, target, reg)
	} else {
		for i := 0; i < reg.size; i++ {
			if reg.node[i].state&(1<<uint(control1)) != 0 {
				if reg.node[i].state&(1<<uint(control2)) != 0 {
					reg.node[i].state ^= (1 << uint(target))
				}
			}
		}
	}
	QuantumDecohere(reg)
}

func QuantumSigmaX(target int, reg *QuantumReg) {
	var qec int

	QuantumQECGetStatus(&qec, nil)

	if qec != 0 {
		QuantumSigmaXft(target, reg)
	} else {
		for i := 0; i < reg.size; i++ {
			reg.node[i].state ^= (1 << uint(target))
		}
	}

	QuantumDecohere(reg)
}

func QuantumSigmaY(target int, reg *QuantumReg) {
	for i := 0; i < reg.size; i++ {
		/* Flip the target bit of each basis state and multiply with
		+/- i */
		reg.node[i].state ^= (1 << uint(target))

		if reg.node[i].state&(1<<uint(target)) != 0 {
			reg.node[i].amplitude *= complex(0, 1)
		} else {
			reg.node[i].amplitude *= complex(0, -1)
		}
	}
	QuantumDecohere(reg)
}

func QuantumSwapTheLeads(width int, reg *QuantumReg) {
	var pat1, pat2 uint32
	var qec int
	var l uint32

	QuantumQECGetStatus(&qec, nil)

	if qec != 0 {
		for i := 0; i < width; i++ {
			QuantumCNOT(i, width+i, reg)
			QuantumCNOT(width+i, i, reg)
			QuantumCNOT(i, width+i, reg)
		}
	} else {
		for i := 0; i < reg.size; i++ {
			pat1 = reg.node[i].state % (1 << uint(width))
			pat2 = 0
			for j := 0; j < width; j++ {
				pat2 += reg.node[i].state & (1 << uint(width+j))
			}
			l = reg.node[i].state - (pat1 + pat2)
			l += (pat1 << uint(width))
			l += (pat2 >> uint(width))
			reg.node[i].state = l
		}
	}
}

func QuantumSwapTheLeadsOMULNControlled(control int, width int, reg *QuantumReg) {
	for i := 0; i < width; i++ {
		QuantumToffoli(control, width+i, 2*width+i+2, reg)
		QuantumToffoli(control, 2*width+i+2, width+i, reg)
		QuantumToffoli(control, width+i, 2*width+i+2, reg)
	}
}

func QuantumGateCounter(inc int) int {
	counter := 0

	if inc > 0 {
		counter += inc
	} else if inc < 0 {
		counter = 0
	}

	return counter
}
