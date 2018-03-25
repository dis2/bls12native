#define X0 R8
#define X1 R9
#define X2 R10
#define X3 R11
#define X4 R12
#define X5 R13

#define Y0 AX
#define Y1 BX
#define Y2 CX
#define Y3 DX
#define Y4 R14
#define Y5 R15



#define LO0 AX
#define LO1 CX
#define LO2 R8
#define LO3 R9
#define LO4 R10
#define LO5 R11

#define HI0 R12
#define HI1 R13
#define HI2 R14
#define HI3 R15
#define HI4 R12
#define HI5 R13





#define LO6 R10
#define LO7 R11
#define LO8 AX
#define LO9 CX
#define LO10 R8
#define LO11 R9

#define HI6 AX
#define HI7 CX
#define HI8 R8
#define HI9 R9
#define HI10 R10
#define HI11 R11



#define D0 0(DI)
#define D1 8(DI)
#define D2 16(DI)
#define D3 24(DI)
#define D4 32(DI)
#define D5 40(DI)
#define D6 48(DI)
#define D7 56(DI)
#define D8 64(DI)
#define D9 72(DI)
#define D10 80(DI)
#define D11 88(DI)

#define OUT0 0(BX)
#define OUT1 8(BX)
#define OUT2 16(BX)
#define OUT3 24(BX)
#define OUT4 32(BX)
#define OUT5 40(BX)
#define OUT6 48(BX)
#define OUT7 56(BX)
#define OUT8 64(BX)
#define OUT9 72(BX)
#define OUT10 80(BX)
#define OUT11 88(BX)


#define S0 0(SI)
#define S1 8(SI)
#define S2 16(SI)
#define S3 24(SI)
#define S4 32(SI)
#define S5 40(SI)
#define S6 48(SI)
#define S7 56(SI)
#define S8 64(SI)
#define S9 72(SI)
#define S10 80(SI)
#define S11 88(SI)

#define OP(op,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5) \
	op a1, b1 \
	op a2, b2 \
	op a3, b3 \
	op a4, b4 \ 
	op a5, b5

#define GET_C(dst) MOVQ c+0(FP), dst
#define GET_A(dst) MOVQ a+8(FP), dst
#define GET_B(dst) MOVQ b+16(FP), dst

#define REDUCE_X \
	MOVQ X0, Y0 \
	OP(MOVQ, 	X1, X2, X3, X4, X5, 	Y1, Y2, Y3, Y4, Y5) \
	SUBQ Q0, Y0 \
	OP(SBBQ, 	Q1, Q2, Q3, Q4, Q5, 	Y1, Y2, Y3, Y4, Y5) \
	CMOVQCC Y0, X0 \
	OP(CMOVQCC, 	Y1, Y2, Y3, Y4, Y5, 	X1, X2, X3, X4, X5)

#define MULINPUT(off) \
	MOVQ INPUT, AX \
	MULQ off*8(SI)

#define MULSTEP1(pos,oreg) \
	MOVQ DX, CARRY \
	MULINPUT(pos) \
	ADDQ CARRY, AX \
	ADCQ $0, DX \
	MOVQ AX, oreg

#define MULSTEP2(pos, ireg, oreg) \
	MOVQ DX, CARRY \
	MULINPUT(pos) \
	ADDQ CARRY, AX \
	ADCQ $0, DX \
	ADDQ ireg, AX \
	ADCQ $0, DX \
	MOVQ AX, oreg

#define REDCSTEP(pos,q,ireg,oreg) \
	MOVQ DX, CARRY \
	MOVQ INPUT, AX \
	MULQ ·Q+pos*8(SB) \
	ADDQ ireg, AX \
	ADCQ $0, DX \
	ADDQ CARRY, AX \
	ADCQ $0, DX \
	MOVQ AX, oreg

#define CARRY R15
#define INPUT R14
#define BASE CX

#define TMP R15
#define PCARRY BX
#define CPTR BX
