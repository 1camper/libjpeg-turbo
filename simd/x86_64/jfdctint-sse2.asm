;
; jfdctint.asm - accurate integer FDCT (64-bit SSE2)
;
; Copyright 2009 Pierre Ossman <ossman@cendio.se> for Cendio AB
; Copyright (C) 2009, 2016, D. R. Commander.
;
; Based on the x86 SIMD extension for IJG JPEG library
; Copyright (C) 1999-2006, MIYASAKA Masaru.
; For conditions of distribution and use, see copyright notice in jsimdext.inc
;
; This file should be assembled with NASM (Netwide Assembler),
; can *not* be assembled with Microsoft's MASM or any compatible
; assembler (including Borland's Turbo Assembler).
; NASM is available from http://nasm.sourceforge.net/ or
; http://sourceforge.net/project/showfiles.php?group_id=6208
;
; This file contains a slow-but-accurate integer implementation of the
; forward DCT (Discrete Cosine Transform). The following code is based
; directly on the IJG's original jfdctint.c; see the jfdctint.c for
; more details.
;
; [TAB8]

%include "jsimdext.inc"
%include "jdct.inc"

; --------------------------------------------------------------------------

%define CONST_BITS  13
%define PASS1_BITS  2

%define DESCALE_P1  (CONST_BITS - PASS1_BITS)
%define DESCALE_P2  (CONST_BITS + PASS1_BITS)

%if CONST_BITS == 13
F_0_298 equ  2446  ; FIX(0.298631336)
F_0_390 equ  3196  ; FIX(0.390180644)
F_0_541 equ  4433  ; FIX(0.541196100)
F_0_765 equ  6270  ; FIX(0.765366865)
F_0_899 equ  7373  ; FIX(0.899976223)
F_1_175 equ  9633  ; FIX(1.175875602)
F_1_501 equ 12299  ; FIX(1.501321110)
F_1_847 equ 15137  ; FIX(1.847759065)
F_1_961 equ 16069  ; FIX(1.961570560)
F_2_053 equ 16819  ; FIX(2.053119869)
F_2_562 equ 20995  ; FIX(2.562915447)
F_3_072 equ 25172  ; FIX(3.072711026)
%else
; NASM cannot do compile-time arithmetic on floating-point constants.
%define DESCALE(x, n)  (((x) + (1 << ((n) - 1))) >> (n))
F_0_298 equ DESCALE( 320652955, 30 - CONST_BITS)  ; FIX(0.298631336)
F_0_390 equ DESCALE( 418953276, 30 - CONST_BITS)  ; FIX(0.390180644)
F_0_541 equ DESCALE( 581104887, 30 - CONST_BITS)  ; FIX(0.541196100)
F_0_765 equ DESCALE( 821806413, 30 - CONST_BITS)  ; FIX(0.765366865)
F_0_899 equ DESCALE( 966342111, 30 - CONST_BITS)  ; FIX(0.899976223)
F_1_175 equ DESCALE(1262586813, 30 - CONST_BITS)  ; FIX(1.175875602)
F_1_501 equ DESCALE(1612031267, 30 - CONST_BITS)  ; FIX(1.501321110)
F_1_847 equ DESCALE(1984016188, 30 - CONST_BITS)  ; FIX(1.847759065)
F_1_961 equ DESCALE(2106220350, 30 - CONST_BITS)  ; FIX(1.961570560)
F_2_053 equ DESCALE(2204520673, 30 - CONST_BITS)  ; FIX(2.053119869)
F_2_562 equ DESCALE(2751909506, 30 - CONST_BITS)  ; FIX(2.562915447)
F_3_072 equ DESCALE(3299298341, 30 - CONST_BITS)  ; FIX(3.072711026)
%endif

; --------------------------------------------------------------------------
    SECTION     SEG_CONST

    alignz      32
    GLOBAL_DATA(jconst_fdct_islow_sse2)

EXTN(jconst_fdct_islow_sse2):

PW_F130_F054   times 4 dw  (F_0_541 + F_0_765),  F_0_541
PW_F054_MF130  times 4 dw  F_0_541, (F_0_541 - F_1_847)
PW_MF078_F117  times 4 dw  (F_1_175 - F_1_961),  F_1_175
PW_F117_F078   times 4 dw  F_1_175, (F_1_175 - F_0_390)
const_ref:
PW_MF060_MF089 times 4 dw  (F_0_298 - F_0_899), -F_0_899
PW_MF089_F060  times 4 dw -F_0_899, (F_1_501 - F_0_899)
PW_MF050_MF256 times 4 dw  (F_2_053 - F_2_562), -F_2_562
PW_MF256_F050  times 4 dw -F_2_562, (F_3_072 - F_2_562)
PD_DESCALE_P1  times 4 dd  1 << (DESCALE_P1 - 1)
PD_DESCALE_P2  times 4 dd  1 << (DESCALE_P2 - 1)
PW_DESCALE_P2X times 8 dw  1 << (PASS1_BITS - 1)

%define REL(c) rax + (c - const_ref)

    alignz      32

; --------------------------------------------------------------------------
    SECTION     SEG_TEXT
    BITS        64
;
; Perform the forward DCT on one block of samples.
;
; GLOBAL(void)
; jsimd_fdct_islow_sse2(DCTELEM *data)
;

%define m0  xmm0
%define m1  xmm1
%define m2  xmm2
%define m3  xmm3
%define m4  xmm4
%define m5  xmm5
%define m6  xmm6
%define m7  xmm7
%define m8  xmm8
%define m9  xmm9
%define m10 xmm10
%define m11 xmm11
%define m12 xmm12
%define m13 xmm13
%define m14 xmm14

%macro SWAP 2.nolist
%xdefine %%tmp m%1
%xdefine m%1 m%2
%xdefine m%2 %%tmp
%endmacro

%ifdef WIN64
%macro IFWIN64 1+
    %1
%endmacro
%else
%macro IFWIN64 1+.nolist
%endmacro
%endif

    align       32
    GLOBAL_FUNCTION(jsimd_fdct_islow_sse2)

EXTN(jsimd_fdct_islow_sse2):
%ifdef WIN64
    %define     data rcx
%else
    %define     data rdi
%endif
    ; ---- Pass 1: process rows.

    movaps      m0, XMMWORD [XMMBLOCK(0,0,data,SIZEOF_DCTELEM)] ; m0=(00 01 02 03 04 05 06 07)
    IFWIN64 mov         rdx, rsp
    movaps      m1, XMMWORD [XMMBLOCK(1,0,data,SIZEOF_DCTELEM)] ; m1=(10 11 12 13 14 15 16 17)
    movdqa      m4, m0                   ; transpose coefficients(phase 1)
    punpcklwd   m0, m1                   ; m0=(00 10 01 11 02 12 03 13)
    IFWIN64 movaps      [rdx + 8 + 0 * SIZEOF_XMMWORD], xmm6 ; shadow space
    movaps      m2, XMMWORD [XMMBLOCK(2,0,data,SIZEOF_DCTELEM)] ; m2=(20 21 22 23 24 25 26 27)
    IFWIN64 add         rsp, 8 - 8 * SIZEOF_XMMWORD
    punpckhwd   m4, m1                   ; m4=(04 14 05 15 06 16 07 17)
    movdqa      m5, m2                   ; transpose coefficients(phase 1)
    IFWIN64 movaps      [rdx + 8 - 7 * SIZEOF_XMMWORD], xmm7
    movaps      m3, XMMWORD [XMMBLOCK(3,0,data,SIZEOF_DCTELEM)] ; m3=(30 31 32 33 34 35 36 37)
    punpcklwd   m2, m3                   ; m2=(20 30 21 31 22 32 23 33)
    IFWIN64 movaps      [rdx + 8 - 6 * SIZEOF_XMMWORD], xmm8
    movaps      m6, XMMWORD [XMMBLOCK(4,0,data,SIZEOF_DCTELEM)] ; m6=( 4 12 20 28 36 44 52 60)
    IFWIN64 movaps      [rdx + 8 - 5 * SIZEOF_XMMWORD], xmm9
    punpckhwd   m5, m3                   ; m5=(24 34 25 35 26 36 27 37)
    SWAP  8, 2 ; movdqa      m8, m2      ; m8=(20 30 21 31 22 32 23 33)  ;  0  1  8  3  4  5  6  7    2  9 10 11 12 13 14
    SWAP  9, 5 ; movdqa      m9, m5      ; m9=(24 34 25 35 26 36 27 37)  ;  0  1  8  3  4  9  6  7    2  5 10 11 12 13 14
    movaps      m7, XMMWORD [XMMBLOCK(5,0,data,SIZEOF_DCTELEM)] ; m7=( 5 13 21 29 37 45 53 61)

    movdqa      m2, m6                   ; transpose coefficients(phase 1)
    punpcklwd   m6, m7                   ; m6=(40 50 41 51 42 52 43 53)
    movaps      m1, XMMWORD [XMMBLOCK(6,0,data,SIZEOF_DCTELEM)] ; m1=( 6 14 22 30 38 46 54 62)
    punpckhwd   m2, m7                   ; m2=(44 54 45 55 46 56 47 57)

    movaps      m3, XMMWORD [XMMBLOCK(7,0,data,SIZEOF_DCTELEM)] ; m3=( 7 15 23 31 39 47 55 63)
    movdqa      m5, m1                   ; transpose coefficients(phase 1)
    punpcklwd   m1, m3                   ; m1=(60 70 61 71 62 72 63 73)
    movdqa      m7, m6                   ; transpose coefficients(phase 2)
    punpckhwd   m5, m3                   ; m5=(64 74 65 75 66 76 67 77)

    movdqa      m3, m2                   ; transpose coefficients(phase 2)
    punpckldq   m6, m1                   ; m6=(40 50 60 70 41 51 61 71)

    IFWIN64 movaps      [rdx + 8 - 2 * SIZEOF_XMMWORD], xmm12

    punpckhdq   m7, m1                   ; m7=(42 52 62 72 43 53 63 73)

    lea         raxp, [rel const_ref]

    punpckldq   m2, m5                   ; m2=(44 54 64 74 45 55 65 75)

    IFWIN64 movaps      [rdx + 8 - 8 * SIZEOF_XMMWORD], xmm13

    punpckhdq   m3, m5                   ; m3=(46 56 66 76 47 57 67 77)

    SWAP  8, 1 ; movdqa      m1, m8      ; m1=(20 30 21 31 22 32 23 33)  ;  0  2  8  3  4  9  6  7    1  5 10 11 12 13 14
    SWAP  9, 5 ; movdqa      m5, m9      ; m5=(24 34 25 35 26 36 27 37)  ;  0  2  8  3  4  5  6  7    1  9 10 11 12 13 14
    SWAP  9, 7 ; movdqa      m9, m7      ; m9=(42 52 62 72 43 53 63 73)  ;  0  2  8  3  4  5  6  9    1  7 10 11 12 13 14
    SWAP  8, 2 ; movdqa      m8, m2      ; m8=(44 54 64 74 45 55 65 75)  ;  0  2  1  3  4  5  6  9    8  7 10 11 12 13 14

    movdqa      m7, m0                   ; transpose coefficients(phase 2)
    punpckldq   m0, m1                   ; m0=(00 10 20 30 01 11 21 31)
    movdqa      m2, m4                   ; transpose coefficients(phase 2)
    punpckldq   m4, m5                   ; m4=(04 14 24 34 05 15 25 35)

    IFWIN64 movaps      [rdx + 8 + 1 * SIZEOF_XMMWORD], xmm14

    punpckhdq   m7, m1                   ; m7=(02 12 22 32 03 13 23 33)
    movdqa      m1, m0                   ; transpose coefficients(phase 3)
    punpckhdq   m2, m5                   ; m2=(06 16 26 36 07 17 27 37)

    movdqa      m5, m2                   ; transpose coefficients(phase 3)
    punpcklqdq  m2, m3                   ; m2=(06 16 26 36 46 56 66 76)=data6

    movaps      m12, XMMWORD [REL(PW_F130_F054)]

    punpcklqdq  m0, m6                   ; m0=(00 10 20 30 40 50 60 70)=data0

    IFWIN64 movaps      [rdx + 8 - 4 * SIZEOF_XMMWORD], xmm10

    punpckhqdq  m1, m6                   ; m1=(01 11 21 31 41 51 61 71)=data1
    movdqa      m6, m1
    psubw       m1, m2                   ; m1=data1-data6=tmp6
    paddw       m6, m2                   ; m6=data1+data6=tmp1
    punpckhqdq  m5, m3                   ; m5=(07 17 27 37 47 57 67 77)=data7

    movdqa      m3, m0
    psubw       m0, m5                   ; m0=data0-data7=tmp7
    paddw       m3, m5                   ; m3=data0+data7=tmp0

    SWAP  9, 2 ; movdqa      m2, m9      ; m2=(42 52 62 72 43 53 63 73)  ;  0  2  7  3  4  5  6  9    8  1 10 11 12 13 14
    SWAP  8, 5 ; movdqa      m5, m8      ; m5=(44 54 64 74 45 55 65 75)  ;  0  2  7  3  4  8  6  9    5  1 10 11 12 13 14
    SWAP  8, 1 ; movdqa      m8, m1      ; m8=tmp6                       ;  0  5  7  3  4  8  6  9    2  1 10 11 12 13 14
    SWAP  9, 0 ; movdqa      m9, m0      ; m9=tmp7                       ;  1  5  7  3  4  8  6  9    2  0 10 11 12 13 14

    movdqa      m1, m7                   ; transpose coefficients(phase 3)
    punpcklqdq  m7, m2                   ; m7=(02 12 22 32 42 52 62 72)=data2

    movaps      m13, XMMWORD [REL(PW_F054_MF130)]

    movdqa      m0, m4                   ; transpose coefficients(phase 3)
    punpckhqdq  m1, m2                   ; m1=(03 13 23 33 43 53 63 73)=data3
    movdqa      m2, m1

    movaps      m14, XMMWORD [REL(PD_DESCALE_P1)]

    punpcklqdq  m4, m5                   ; m4=(04 14 24 34 44 54 64 74)=data4
    paddw       m1, m4                   ; m1=data3+data4=tmp3
    psubw       m2, m4                   ; m2=data3-data4=tmp4

    ; -- Even part

    movdqa      m4, m3
    punpckhqdq  m0, m5                   ; m0=(05 15 25 35 45 55 65 75)=data5
    movdqa      m5, m7
    paddw       m3, m1                   ; m3=tmp10
    psubw       m4, m1                   ; m4=tmp13
    paddw       m7, m0                   ; m7=data2+data5=tmp2
    psubw       m5, m0                   ; m5=data2-data5=tmp5
    movdqa      m0, m6
    movdqa      m1, m3
    paddw       m6, m7                   ; m6=tmp11
    psubw       m0, m7                   ; m0=tmp12
    movdqa      m7, m4                   ; m4=tmp13
    paddw       m3, m6                   ; m3=tmp10+tmp11
    psubw       m1, m6                   ; m1=tmp10-tmp11
    movdqa      m6, m4
    punpcklwd   m7, m0                   ; m0=tmp12
    psllw       m3, PASS1_BITS           ; m3=data0
    psllw       m1, PASS1_BITS           ; m1=data4

    ; (Original)
    ; z1 = (tmp12 + tmp13) * 0.541196100;
    ; data2 = z1 + tmp13 * 0.765366865;
    ; data6 = z1 + tmp12 * -1.847759065;
    ;
    ; (This implementation)
    ; data2 = tmp13 * (0.541196100 + 0.765366865) + tmp12 * 0.541196100;
    ; data6 = tmp13 * 0.541196100 + tmp12 * (0.541196100 - 1.847759065);

    IFWIN64 movaps      [rdx + 8 - 3 * SIZEOF_XMMWORD], xmm11

    punpckhwd   m6, m0
    movdqa      m4, m7
    pmaddwd     m7, m12                  ; m7=data2L
    movdqa      m0, m6
    pmaddwd     m6, m12                  ; m6=data2H

    SWAP 10, 3 ; movdqa      m10, m3     ; m10=data0                     ;  1  5  7 10  4  8  6  9    2  0  3 11 12 13 14
    SWAP 12, 1 ; movdqa      m12, m1     ; m12=data4                     ;  1 12  7 10  4  8  6 11    2  0  3  9  5 13 14
    SWAP  8, 3 ; movdqa      m3, m8      ; m3=tmp6                       ;  1 12  7  2 13  8  6 11   10  0  3  9  5  4 14
    SWAP  9, 1 ; movdqa      m1, m9      ; m1=tmp7                       ;  1  0  7  2 13  8  6 11   10 12  3  9  5  4 14

    pmaddwd     m4, m13                  ; m4=data6L
    paddd       m7, m14
    pmaddwd     m0, m13                  ; m0=data6H
    paddd       m6, m14
    paddd       m4, m14
    psrad       m7, DESCALE_P1
    paddd       m0, m14
    psrad       m6, DESCALE_P1
    psrad       m4, DESCALE_P1
    packssdw    m7, m6                   ; m7=data2
    psrad       m0, DESCALE_P1
    movdqa      m6, m2                   ; m2=tmp4

    movaps      m8, XMMWORD [REL(PW_MF078_F117)]

    packssdw    m4, m0                   ; m4=data6

    SWAP 11, 7 ; movdqa      m11, m7     ; m11=data2                     ;  1  5  7 10  4  8  6 11    2  0  3  9 12 13 14
    SWAP 13, 4 ; movdqa      m13, m4     ; m13=data6                     ;  1 12  7 10 13  8  6 11    2  0  3  9  5  4 14

    ; -- Odd part

    movdqa      m0, m5                   ; m5=tmp5
    paddw       m6, m3                   ; m6=z3

    movaps      m9, XMMWORD [REL(PW_F117_F078)]

    paddw       m0, m1                   ; m0=z4

    ; (Original)
    ; z5 = (z3 + z4) * 1.175875602;
    ; z3 = z3 * -1.961570560;  z4 = z4 * -0.390180644;
    ; z3 += z5;  z4 += z5;
    ;
    ; (This implementation)
    ; z3 = z3 * (1.175875602 - 1.961570560) + z4 * 1.175875602;
    ; z4 = z3 * 1.175875602 + z4 * (1.175875602 - 0.390180644);

    movdqa      m7, m6
    movdqa      m4, m6
    punpcklwd   m7, m0
    punpckhwd   m4, m0
    movdqa      m6, m7
    pmaddwd     m7, m8                   ; m7=z3L
    movdqa      m0, m4
    pmaddwd     m4, m8                   ; m4=z3H
    SWAP  8, 7 ; movdqa      m8, m7      ; m8=z3L                       ;  1  0  7  2 13  8  6 10   11 12  3  9  5  4 14
    pmaddwd     m6, m9                   ; m6=z4L
    movdqa      m7, m2
    pmaddwd     m0, m9                   ; m0=z4H
    SWAP  9, 4 ; movdqa      m9, m4      ; m9=z3H                       ;  1  0  7  2 12  8  6 10   11 13  3  9  5  4 14

    ; (Original)
    ; z1 = tmp4 + tmp7;  z2 = tmp5 + tmp6;
    ; tmp4 = tmp4 * 0.298631336;  tmp5 = tmp5 * 2.053119869;
    ; tmp6 = tmp6 * 3.072711026;  tmp7 = tmp7 * 1.501321110;
    ; z1 = z1 * -0.899976223;  z2 = z2 * -2.562915447;
    ; data7 = tmp4 + z1 + z3;  data5 = tmp5 + z2 + z4;
    ; data3 = tmp6 + z2 + z3;  data1 = tmp7 + z1 + z4;
    ;
    ; (This implementation)
    ; tmp4 = tmp4 * (0.298631336 - 0.899976223) + tmp7 * -0.899976223;
    ; tmp5 = tmp5 * (2.053119869 - 2.562915447) + tmp6 * -2.562915447;
    ; tmp6 = tmp5 * -2.562915447 + tmp6 * (3.072711026 - 2.562915447);
    ; tmp7 = tmp4 * -0.899976223 + tmp7 * (1.501321110 - 0.899976223);
    ; data7 = tmp4 + z3;  data5 = tmp5 + z4;
    ; data3 = tmp6 + z3;  data1 = tmp7 + z4;

    movdqa      m4, m2
    punpcklwd   m7, m1
    punpckhwd   m4, m1
    movdqa      m2, m7
    pmaddwd     m7, [REL(PW_MF060_MF089)] ; m7=tmp4L
    movdqa      m1, m4
    pmaddwd     m4, [REL(PW_MF060_MF089)] ; m4=tmp4H
    pmaddwd     m2, [REL(PW_MF089_F060)]  ; m2=tmp7L
    paddd       m7, m8                    ; m7=data7L
    pmaddwd     m1, [REL(PW_MF089_F060)]  ; m1=tmp7H

    paddd       m4, m9                   ; m4=data7H
    paddd       m2, m6                   ; m2=data1L
    paddd       m7, m14
    paddd       m1, m0                   ; m1=data1H

    paddd       m4, m14
    paddd       m2, m14
    psrad       m7, DESCALE_P1
    paddd       m1, m14
    psrad       m4, DESCALE_P1
    psrad       m2, DESCALE_P1
    psrad       m1, DESCALE_P1

    packssdw    m7, m4                   ; m7=data7
    movdqa      m4, m5
    packssdw    m2, m1                   ; m2=(01 11 21 31 41 51 61 71)

    movdqa      m1, m5
    punpcklwd   m4, m3
    movdqa      m5, m4
    punpckhwd   m1, m3
    movdqa      m3, m1
    pmaddwd     m4, [REL(PW_MF050_MF256)] ; m4=tmp5L
    pmaddwd     m5, [REL(PW_MF256_F050)]  ; m5=tmp6L
    pmaddwd     m1, [REL(PW_MF050_MF256)] ; m1=tmp5H
    pmaddwd     m3, [REL(PW_MF256_F050)]  ; m3=tmp6H

    paddd       m4, m6                   ; m4=data5L
    paddd       m5, m8                   ; m5=data3L
    paddd       m1, m0                   ; m1=data5H
    paddd       m3, m9                   ; m3=data3H

    SWAP 10, 6 ; movdqa      m6, m10     ; m6=(00 10 20 30 40 50 60 70)  ;  1  0  7  2 12  8  3 10   11 13  6  9  5  4 14
    SWAP 11, 0 ; movdqa      m0, m11     ; m0=(02 12 22 32 42 52 62 72)  ;  9  0  7  2 12  8  3 10   11 13  6  1  5  4 14

    paddd       m4, m14
    paddd       m5, m14
    paddd       m1, m14
    paddd       m3, m14
    psrad       m4, DESCALE_P1
    psrad       m5, DESCALE_P1
    psrad       m1, DESCALE_P1
    psrad       m3, DESCALE_P1

    packssdw    m4, m1                   ; m4=data5

    ; ---- Pass 2: process columns.

    movdqa      m1, m6                   ; transpose coefficients(phase 1)
    packssdw    m5, m3                   ; m5=(03 13 23 33 43 53 63 73)
    movdqa      m3, m0                   ; transpose coefficients(phase 1)
    punpcklwd   m6, m2                   ; m6=(00 01 10 11 20 21 30 31)
    punpckhwd   m1, m2                   ; m1=(40 41 50 51 60 61 70 71)
    SWAP 12, 2 ; movdqa      m2, m12     ; m2=col4                       ;  9  0  5  2 12  8  3 10   11 13  6  1  7  4 14
    punpcklwd   m0, m5                   ; m0=(02 03 12 13 22 23 32 33)
    SWAP 10, 0 ; movdqa      m10, m0     ; m10=(02 03 12 13 22 23 32 33) ;  6  0  5  2 12  8  3 10   11 13  9  1  7  4 14
    movdqa      m0, m2                   ; transpose coefficients(phase 1)
    punpckhwd   m3, m5                   ; m3=(42 43 52 53 62 63 72 73)
    SWAP 13, 5 ; movdqa      m5, m13     ; m5=col6                       ;  6  0  5  2 12  4  3 10   11 13  9  1  7  8 14
    SWAP 11, 3 ; movdqa      m11, m3     ; m11=(42 43 52 53 62 63 72 73) ;  6  0  5  1 12  4  3 10   11 13  9  2  7  8 14

    ; m2=(04 14 24 34 44 54 64 74), m5=(06 16 26 36 46 56 66 76)
    ; m4=(05 15 25 35 45 55 65 75), m7=(07 17 27 37 47 57 67 77)

    IFWIN64 movaps      xmm8, [rdx + 8 - 6 * SIZEOF_XMMWORD] ; m13

    movdqa      m3, m5                   ; transpose coefficients(phase 1)
    punpcklwd   m2, m4                   ; m2=(04 05 14 15 24 25 34 35)
    punpcklwd   m5, m7                   ; m5=(06 07 16 17 26 27 36 37)
    punpckhwd   m0, m4                   ; m0=(44 45 54 55 64 65 74 75)
    punpckhwd   m3, m7                   ; m3=(46 47 56 57 66 67 76 77)

    SWAP 12, 4                                                           ;  6  0  5  1  7  4  3 10   11 13  9  2 12  8 14
    movdqa      m4, m2                   ; transpose coefficients(phase 2)
    movdqa      m7, m0                   ; transpose coefficients(phase 2)
    punpckldq   m2, m5                   ; m2=(04 05 06 07 14 15 16 17)
    punpckldq   m0, m3                   ; m0=(44 45 46 47 54 55 56 57)
    punpckhdq   m4, m5                   ; m4=(24 25 26 27 34 35 36 37)
    SWAP 10, 5 ; movdqa      m5, m10     ; m5=(02 03 12 13 22 23 32 33)  ;  6  0  5  1  7  9  3 10   11 13  4  2 12  8 14
    SWAP 10, 4 ; movdqa      m10, m4     ; m10=(24 25 26 27 34 35 36 37) ;  6  0  5  2  4  9  3 10   11 13  7  1 12  8 14
    movdqa      m4, m6                   ; transpose coefficients(phase 2)
    punpckhdq   m7, m3                   ; m7=(64 65 66 67 74 75 76 77)
    SWAP 11, 3 ; movdqa      m3, m11     ; m3=(42 43 52 53 62 63 72 73)  ;  6  0  5  2  7  9  3 10   11 13  4  1 12  8 14
    SWAP 11, 0 ; movdqa      m11, m0     ; m11=(44 45 46 47 54 55 56 57) ;  1  0  5  2  4  9  3 10   11 13  7  6 12  8 14

    IFWIN64 movaps      xmm12, [rdx + 8 - 2 * SIZEOF_XMMWORD] ; m12

    movdqa      m0, m1                   ; transpose coefficients(phase 2)
    punpckldq   m6, m5                   ; m6=(00 01 02 03 10 11 12 13)
    punpckhdq   m4, m5                   ; m4=(20 21 22 23 30 31 32 33)
    punpckldq   m1, m3                   ; m1=(40 41 42 43 50 51 52 53)
    punpckhdq   m0, m3                   ; m0=(60 61 62 63 70 71 72 73)

    movdqa      m5, m6                   ; transpose coefficients(phase 3)
    punpcklqdq  m6, m2                   ; m6=(00 01 02 03 04 05 06 07)=data0
    movdqa      m3, m0                   ; transpose coefficients(phase 3)
    punpcklqdq  m0, m7                   ; m0=(60 61 62 63 64 65 66 67)=data6

    movaps      m14, [REL(PD_DESCALE_P2)]

    punpckhqdq  m5, m2                   ; m5=(10 11 12 13 14 15 16 17)=data1
    punpckhqdq  m3, m7                   ; m3=(70 71 72 73 74 75 76 77)=data7

    movdqa      m7, m6
    movdqa      m2, m5
    psubw       m5, m0                   ; m5=data1-data6=tmp6
    psubw       m6, m3                   ; m6=data0-data7=tmp7
    paddw       m7, m3                   ; m7=data0+data7=tmp0
    SWAP 11, 3 ; movdqa      m3, m11     ; m3=(44 45 46 47 54 55 56 57)  ;  7  0  5  6  4  9  3 10   11 13  1  2 12  8 14
    SWAP 11, 6 ; movdqa      m11, m6     ; m11=tmp7                      ;  7  0  5  6  4  1  2 10   11 13  9  3 12  8 14
    paddw       m2, m0                   ; m2=data1+data6=tmp1
    SWAP 10, 0 ; movdqa      m0, m10     ; m0=(24 25 26 27 34 35 36 37)  ;  7  0  5  2  4  9  3 10   11 13  1  6 12  8 14
    SWAP 10, 5 ; movdqa      m10, m5     ; m10=tmp6                      ;  7  0  5  6  4  1  3 10   11 13  9  2 12  8 14

    movdqa      m5, m4                   ; transpose coefficients(phase 3)
    punpcklqdq  m4, m0                   ; m4=(20 21 22 23 24 25 26 27)=data2
    movdqa      m6, m1                   ; transpose coefficients(phase 3)
    punpcklqdq  m1, m3                   ; m1=(40 41 42 43 44 45 46 47)=data4

    movaps      m8, [REL(PW_DESCALE_P2X)]

    punpckhqdq  m5, m0                   ; m5=(30 31 32 33 34 35 36 37)=data3

    movaps      m9, [REL(PW_F054_MF130)]

    punpckhqdq  m6, m3                   ; m6=(50 51 52 53 54 55 56 57)=data5

    movdqa      m3, m4
    movdqa      m0, m5
    paddw       m5, m1                   ; m5=data3+data4=tmp3
    paddw       m4, m6                   ; m4=data2+data5=tmp2
    psubw       m3, m6                   ; m3=data2-data5=tmp5
    psubw       m0, m1                   ; m0=data3-data4=tmp4

    ; -- Even part

    movdqa      m1, m7
    movdqa      m6, m2
    paddw       m7, m5                   ; m7=tmp10
    paddw       m2, m4                   ; m2=tmp11
    psubw       m1, m5                   ; m1=tmp13
    movdqa      m5, m7
    paddw       m7, m2                   ; m7=tmp10+tmp11
    psubw       m6, m4                   ; m6=tmp12
    movdqa      m4, m1                   ; m1=tmp13
    psubw       m5, m2                   ; m5=tmp10-tmp11

    paddw       m7, m8
    movdqa      m2, m1
    punpcklwd   m4, m6                   ; m6=tmp12
    paddw       m5, m8

    movaps      m8, [REL(PW_F130_F054)]

    psraw       m7, PASS1_BITS           ; m7=data0
    psraw       m5, PASS1_BITS           ; m5=data4

    ; (Original)
    ; z1 = (tmp12 + tmp13) * 0.541196100;
    ; data2 = z1 + tmp13 * 0.765366865;
    ; data6 = z1 + tmp12 * -1.847759065;
    ;
    ; (This implementation)
    ; data2 = tmp13 * (0.541196100 + 0.765366865) + tmp12 * 0.541196100;
    ; data6 = tmp13 * 0.541196100 + tmp12 * (0.541196100 - 1.847759065);
    punpckhwd   m2, m6
    movdqa      m1, m4
    pmaddwd     m4, m8                   ; m4=data2L
    movaps      XMMWORD [XMMBLOCK(0,0,data,SIZEOF_DCTELEM)], m7
    movdqa      m6, m2
    pmaddwd     m2, m8                   ; m2=data2H
    paddd       m4, m14
    pmaddwd     m1, m9                   ; m1=data6L
    pmaddwd     m6, m9                   ; m6=data6H

    paddd       m2, m14
    psrad       m4, DESCALE_P2
    paddd       m1, m14
    paddd       m6, m14
    psrad       m2, DESCALE_P2
    movaps      XMMWORD [XMMBLOCK(4,0,data,SIZEOF_DCTELEM)], m5
    psrad       m1, DESCALE_P2
    psrad       m6, DESCALE_P2

    packssdw    m4, m2                   ; m4=data2
    movdqa      m2, m0                   ; m0=tmp4
    ; -- Odd part

    SWAP 10, 7 ; movdqa      m7, m10     ; m7=tmp6                       ;  7  0  5  6  4  1  2  9   11 13 10  3 12  8 14
    SWAP 11, 5 ; movdqa      m5, m11     ; m5=tmp7                       ;  7  0  5  6  4  3  2  9   11 13 10  1 12  8 14

    movaps      m10, XMMWORD [REL(PW_MF078_F117)]

    packssdw    m1, m6                   ; m1=data6
    movdqa      m6, m3                   ; m3=tmp5

    movaps      XMMWORD [XMMBLOCK(2,0,data,SIZEOF_DCTELEM)], m4
    paddw       m2, m7                   ; m2=z3
    paddw       m6, m5                   ; m6=z4
    movaps      XMMWORD [XMMBLOCK(6,0,data,SIZEOF_DCTELEM)], m1

    ; (Original)
    ; z5 = (z3 + z4) * 1.175875602;
    ; z3 = z3 * -1.961570560;  z4 = z4 * -0.390180644;
    ; z3 += z5;  z4 += z5;
    ;
    ; (This implementation)
    ; z3 = z3 * (1.175875602 - 1.961570560) + z4 * 1.175875602;
    ; z4 = z3 * 1.175875602 + z4 * (1.175875602 - 0.390180644);

    movdqa      m4, m2
    movdqa      m1, m2

    movaps      m11, XMMWORD [REL(PW_F117_F078)]

    punpcklwd   m4, m6

    movaps      m9, XMMWORD [REL(PW_MF060_MF089)]

    punpckhwd   m1, m6
    movdqa      m2, m4
    pmaddwd     m4, m10                  ; m4=z3L
    SWAP  8, 4 ; movdqa      m8, m4      ; m8=z3L                        ;  7  0  5  6 11  3  2  9    4 13 10  1 12  8 14
    movdqa      m4, m0
    movdqa      m6, m1
    pmaddwd     m1, m10                  ; m1=z3H
    pmaddwd     m2, m11                  ; m2=z4L
    punpcklwd   m4, m5
    pmaddwd     m6, m11                  ; m6=z4H
    SWAP 11, 1 ; movdqa      m11, m1     ; m11=z3H                       ;  7  1  5  6 11  3  2  9    4 13 10  0 12  8 14

    ; (Original)
    ; z1 = tmp4 + tmp7;  z2 = tmp5 + tmp6;
    ; tmp4 = tmp4 * 0.298631336;  tmp5 = tmp5 * 2.053119869;
    ; tmp6 = tmp6 * 3.072711026;  tmp7 = tmp7 * 1.501321110;
    ; z1 = z1 * -0.899976223;  z2 = z2 * -2.562915447;
    ; data7 = tmp4 + z1 + z3;  data5 = tmp5 + z2 + z4;
    ; data3 = tmp6 + z2 + z3;  data1 = tmp7 + z1 + z4;
    ;
    ; (This implementation)
    ; tmp4 = tmp4 * (0.298631336 - 0.899976223) + tmp7 * -0.899976223;
    ; tmp5 = tmp5 * (2.053119869 - 2.562915447) + tmp6 * -2.562915447;
    ; tmp6 = tmp5 * -2.562915447 + tmp6 * (3.072711026 - 2.562915447);
    ; tmp7 = tmp4 * -0.899976223 + tmp7 * (1.501321110 - 0.899976223);
    ; data7 = tmp4 + z3;  data5 = tmp5 + z4;
    ; data3 = tmp6 + z3;  data1 = tmp7 + z4;

    movdqa      m1, m0
    movaps      m10, XMMWORD [REL(PW_MF089_F060)]

    movdqa      m0, m4
    pmaddwd     m4, m9                   ; m4=tmp4L
    punpckhwd   m1, m5
    pmaddwd     m0, m10                  ; m0=tmp7L
    movdqa      m5, m1
    paddd       m4, m8                   ; m4=data7L
    pmaddwd     m1, m9                   ; m1=tmp4H
    paddd       m0, m2                   ; m0=data1L
    pmaddwd     m5, m10                  ; m5=tmp7H

    paddd       m4, m14
    paddd       m1, m11                  ; m1=data7H
    paddd       m0, m14
    paddd       m5, m6                   ; m5=data1H

    psrad       m4, DESCALE_P2
    paddd       m1, m14
    psrad       m0, DESCALE_P2
    paddd       m5, m14

    movaps      m9, XMMWORD [REL(PW_MF050_MF256)]

    psrad       m1, DESCALE_P2
    psrad       m5, DESCALE_P2

    packssdw    m4, m1                   ; m4=data7
    movdqa      m1, m3
    packssdw    m0, m5                   ; m0=data1
    movdqa      m5, m3

    movaps      XMMWORD [XMMBLOCK(7,0,data,SIZEOF_DCTELEM)], m4

    punpcklwd   m1, m7

    IFWIN64 movaps      xmm11, [rdx + 8 - 3 * SIZEOF_XMMWORD] ; m4

    punpckhwd   m5, m7

    movaps      XMMWORD [XMMBLOCK(1,0,data,SIZEOF_DCTELEM)], m0

    movaps      m0, XMMWORD [REL(PW_MF256_F050)]

    movdqa      m3, m1
    pmaddwd     m1, m9                   ; m1=tmp5L
    movdqa      m7, m5
    pmaddwd     m5, m9                   ; m5=tmp5H
    pmaddwd     m3, m0                   ; m3=tmp6L
    IFWIN64 movaps      xmm13, [rdx + 8 - 8 * SIZEOF_XMMWORD] ; m9

    paddd       m1, m2                   ; m1=data5L
    pmaddwd     m7, m0                   ; m7=tmp6H

    IFWIN64 movaps      xmm10, [rdx + 8 - 4 * SIZEOF_XMMWORD] ; m10

    paddd       m5, m6                   ; m5=data5H
    paddd       m3, m8                   ; m3=data3L

    SWAP 11, 7                                                           ;  7  1  5  6 11  3  2  0    4 13 10  9 12  8 14

    paddd       m7, m11                  ; m7=data3H

    IFWIN64 movaps      xmm9, [rdx + 8 - 5 * SIZEOF_XMMWORD] ; m11

    paddd       m1, m14
    paddd       m5, m14

    IFWIN64 movaps      xmm7, [rdx + 8 - 7 * SIZEOF_XMMWORD] ; m0
    IFWIN64 mov         rsp, rdx

    paddd       m3, m14
    paddd       m7, m14
    psrad       m1, DESCALE_P2
    psrad       m5, DESCALE_P2

    IFWIN64 movaps      xmm14, [rdx + 8 + 1 * SIZEOF_XMMWORD] ; m14

    psrad       m3, DESCALE_P2
    psrad       m7, DESCALE_P2

    packssdw    m1, m5                   ; m1=data5
    packssdw    m3, m7                   ; m3=data3

    movaps      XMMWORD [XMMBLOCK(5,0,data,SIZEOF_DCTELEM)], m1
    movaps      XMMWORD [XMMBLOCK(3,0,data,SIZEOF_DCTELEM)], m3

    IFWIN64 movaps      xmm6, [rdx + 8 + 0 * SIZEOF_XMMWORD] ; m3
    ret

; For some reason, the OS X linker does not honor the request to align the
; segment unless we do this.
    align       32
