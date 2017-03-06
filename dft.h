/******************************************************************
 * @file         dft.h
 * @date         2015-11-19 22:03:08
 * @author       Yao
 * @e-mail       yao.jiang@tongji.edu.cn
 * @brief        Discrete Fourier Transfrom
 *
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <dos.h>
#include <windows.h>

#define size 1024
#define pi 3.141592654

typedef struct {
	short real;
	short imag;
}complex; 


int inv;
long npt;
complex x[size];
//WORD data[1024];//表示输入的一个数;
FILE *in, *out;

complex EE(complex b1, complex b2);  //复数相乘
void read_data();               //文件
void save_data();
void initdata(int n);
void dft();
void fft_dit();
void testEE();
void sin512();