/******************************************************************
 * @file         dft.c
 * @date         2015-11-19 22:32:46
 * @author       Yao
 * @e-mail       yao.jiang@tongji.edu.cn
 * @brief        Discrete Fourier Transfrom
 *
 *****************************************************************/

#include "dft.h"

/************************************************************************/
/* 两个复数相乘 ，复数相乘的话，进行截位  ,定点话计算                                */
/************************************************************************/

complex EE(complex b1, complex b2) {
	complex b3;
	int temp;
	int temp1;
	int temp2;
	int temp3;
	b1.real = b1.real & 0xffff; //将进行计算的数据进行，进行规格化
	b1.imag = b1.imag & 0xffff; //0xffff = 1111 1111 1111 1111

	//旋转因子的位数是8bit加上符号位9bit
	b2.real = b2.real & 0x01ff; //统一的有效数据是八位
	b2.imag = b2.imag & 0x01ff; //

	/*if (b1.real >= 32768) { //如果第八位是符号位为一位 0x0080 = 1000 0000 0000 0000
		b1.real = b1.real + 0x8000+0x8000;  //补码之后再加上符号位
	}
	if (b1.imag >= 32768) {
		b1.imag = b1.imag + 0x8000 + 0x8000;  
	}
	if (b2.real >= 256) {
		b2.real = (b2.real + 512)&0x01ff;  
	}
	if (b2.imag >= 256) {
		b2.imag = (b2.imag + 512)&0x01ff; 
	}*/

	//对于计算的数据进行截位处理
	temp = (((int)(b1.real*b2.real) & 0x01ffffff) >> 9) & 0xffff;  //0x1ffffff = 1 1111 1111 1111 1111 1111 1111
	temp1 = (((int)(b1.imag*b2.imag) & 0x01ffffff) >> 9) & 0xffff;  //这里的意思是取第14位到第七位
	temp2 = (((int)(b1.real*b2.imag) & 0x01ffffff) >> 9) & 0xffff;
	temp3 = (((int)(b2.real*b1.imag) & 0x01ffffff) >> 9) & 0xffff;
	b3.real = (temp - temp1)&0xffff;                     //每次计算完成之后保证计算的结果只有8bit数据
	b3.imag = (temp3 + temp2)&0xffff;
	return (b3);
}


/************************************************************************/
/* 初始化数据                                                            */
/************************************************************************/
void initdata(int n) {
	extern long npt;
	extern int inv;
	extern complex x[size];
	int i = 0;
	char creal[17] = { 0 };
	char cimag[17] = { 0 };
	short real;
	short imag;
	int k = 0;
	int m;
	errno_t err;
	err = fopen_s(&out, "E:\\Projects\\20151119-DFT-Project\\Debug\\inputdata.dat", "w");
	if (err)
	{
		printf("Cannot open file inputdata.dat\n");
		system("pause");
		exit(1);
	}

	switch (n) {
	case 0:
		npt = 0;
		inv = 0;
		for (i = 0;i <= size;++i) {
			x[i].real = 0.0;
			x[i].imag = 0.0;
		}
		break;
	case 1:
		npt = 128;
		for (i = 0;i < npt;++i) {
			x[i].real = sin(pi*i / 256);
			x[i].imag = 0.0;
		}
		break;
	case 2:
		npt = 128;
		for (i = 0;i < npt;++i) {
			x[i].real = cos(-2 * pi*i / npt);
			x[i].imag = sin(-2 * pi*i / npt);
		}
		break;
	case 3:
		npt = 1024;
		for (int j = 0;j < 8;++j) {
			for (i = 0;i < 128;++i) {
				x[k].real = i & 0x00ff;
				x[k].imag = 0 & 0x00ff;
				for (m = 0;m < 8;++m) //将旋转因子输出来
				{
					real = x[k].real;
					creal[7 - m] = ((real >> m) & 0x0001) + 48;
					imag = x[k].imag;
					cimag[7 - m] = ((imag >> m) & 0x0001) + 48;
				}
				fprintf(out, "%s%s\n", creal, cimag);
				k = k + 1;
			}
		}
		break;
	case 4:
		npt = 1024;
			for (i = 0;i < 1024;++i) {
				x[i].real = i & 0xffff;
				x[i].imag = 0 & 0xffff;
				for (m = 0;m < 16;++m) 
				{
					real = x[i].real;
					creal[15 - m] = ((real >> m) & 0x0001) + 48;
					imag = x[i].imag;
					cimag[15 - m] = ((imag >> m) & 0x0001) + 48;
				}
				fprintf(out, "%s%s\n", creal, cimag);
			}
	}
	fclose(out);
}



/************************************************************************/
/* fft_dit                                                              */
/************************************************************************/
void fft_dit() {
	extern complex x[size];
	extern long npt;
	extern int inv;
	short A = pow(2, 8) - 1;
	int f, m, LH, nm, i, j, k, L;
	double p, ps;
	int le, B, ip;
	complex w, t;
	LH = npt / 2;
	f = npt;
	for (m = 1;(f = f / 2) != 1;++m);//2的m次方=npt,m为级数
	nm = npt - 2;
	j = npt / 2;
	//变址计算
	for (i = 1;i <= nm;++i) {
		if (i < j) {
			t = x[j];
			x[j] = x[i];
			x[i] = t;
		}
		k = LH;
		while (j >= k) {
			j = j - k;
			k = k / 2;
		}
		j = j + k;
	}
	//fft_dit运算
	for (L = 0;L <= m - 1;++L) {           //每一级
		le = (int)pow(2.0, L + 1);         //具有相同旋转因子的蝶形运算的间隔---蝶距
		B = le / 2;                        //2的L次方，不同旋转因子的个数，蝶形运算两个间隔 同一蝶形运算两输入数据的距离 
		for (j = 0;j <= B - 1;++j) {
			p = pow(2.0, m - L - 1)*j;
			ps = 2 * pi *p / npt;          //旋转因子
			w.real = cos(ps)*A;
			w.imag = sin(ps)*A;
			if (w.real >= 0) {
				w.real = w.real & 0x01ff;  //直接变成补码正数的补码是其本身
			}
			else {
				w.real = (w.real + 512) & 0x01ff;//将负数转化为补码再加上符号位
			}
			if (w.imag >= 0) {
				w.imag = w.imag & 0x01ff;      //直接变成补码正数的补码是其本身
			}
			else {
				w.imag = (w.imag + 512) & 0x01ff;//将负数转化为补码再加上符号位
			}
			for (i = j;i <= npt - 1;i = i + le) {
				ip = i + B;
				t = EE(x[ip], w);

				x[ip].real = (x[i].real - t.real) & 0xffff;
				x[ip].imag = (x[i].imag - t.imag) & 0xffff;
				x[i].real = (x[i].real + t.real) & 0xffff;
				x[i].imag = (x[i].imag + t.imag) & 0xffff;
			}
		}
	}

}

void sin512() {
	int i, j;
	short real;
	short imag;
	short real1;
	short imag1;
	short A =255;
	char creal[10] = { 0 };
	char cimag[10] = { 0 };
	int m;
	int index;
	char cindex[7] = { 0 };
	errno_t err;
	err = fopen_s(&out, "E:\\1024data.dat", "w");
	if (err)
	{
		printf("Cannot open file 1024data.dat\n");
		system("pause");
		exit(1);
	}
	int k = 0;
	for (i = 0;i < 16;++i)
	{
		for (j = 0;j < 4;++j)
		{
			real = (short)(A*cos(-2 * 16 *i*j*pi / 1024));
			if (real >= 0) {
				real = real & 0x01ff;
			}
			else {
				real = (real + 512) & 0x01ff;//加上符号位就好了
			}
			imag = (short)(A*sin(-2 * 16 *i*j*pi / 1024));
			if (imag >= 0) {
				imag = imag & 0x01ff; 			}
			else {
				imag = (imag + 512) & 0x01ff;//将负数转化为补码再加上符号位
			}
			for (m = 0;m < 9;++m) //将旋转因子输出来
			{
				real1 = real;
				creal[8 - m] = ((real >> m) & 0x0001) + 48;
				imag1 = imag;
				cimag[8 - m] = ((imag >> m) & 0x0001) + 48;
			}
			for (m = 0;m < 6;++m) //将旋转因子输出来
			{
				index = k;
				cindex[5 - m] = ((index >> m) & 0x0001) + 48;
			}
			fprintf(out, "6'b%s:ROM_data<=18'b%s%s\n", cindex,creal, cimag);
			k = k + 1;
		}
	}
	fclose(out);
	printf("Sin 1024 计算完毕\n");
}


/************************************************************************/
/* 读取文件数据                                                         */
/************************************************************************/

void read_data() {
	extern long npt;
	extern complex x[size];
	int n;
	int num;

	//  [11/19/2015/22:43 Yao]
	errno_t err;
	for (n = 0;n < size;++n) {
		x[n].real = 0;
		x[n].imag =  0;
	}

	//  [11/19/2015/22:42 Yao] if ((in = fopen("list.dat", "rt")) == NULL)
	err = fopen_s(&in, "E:\\Projects\\20151119-DFT-Project\\Debug\\dft.dat", "rt");
	if(err)
	{
		printf("Cannot open the file list.dat\n");
		system("pause");
		exit(1);
	}
	fscanf_s(in, "%ld", &npt);
	for (n = 0;n < npt;++n) {
		fscanf_s(in, "%d %lf %lf",&num, &x[n].real, &x[n].imag);
	}
	fclose(in);
}

/************************************************************************/
/* 输出数据到文件中                                                      */
/************************************************************************/

void save_data() {
	long k;
	int k1;
	extern long npt;
	extern complex x[size];
	char creal[17] = { 0 };
	char cimag[17] = { 0 };
	int i;
	short real;
	short imag;
	errno_t err;
	err = fopen_s(&out, "E:\\Projects\\20151119-DFT-Project\\Debug\\dftout.dat", "w");
	if(err)
	{
		printf("Cannot open file dftout.dat\n");
		system("pause");
		exit(1);
	}
	for (k = 0;k < npt; ++k) {
		k1 = k ;
		for (i = 0;i < 16;++i) //将旋转因子输出来
		{
			real = x[k].real;
			creal[15 - i] = ((real >> i) & 0x0001) + 48;
			imag = x[k].imag;
			cimag[15 - i] = ((imag >> i) & 0x0001) + 48;
		}
		fprintf(out,"%s%s\n", creal, cimag);
		//fprintf(out, "%d%d\n", x[k].real, x[k].imag);
	}
	fclose(out);
}

/*
void save_data() {
	short k;
	char creal[11] = { 0 };
	int i;
	int n;
	char s[1024][19] = { 0 };
	errno_t err;
	err = fopen_s(&in, "E:\\Projects\\20151119-DFT-Project\\Debug\\w_im_bin_flow5.txt", "rt");
	if (err)
	{
		printf("Cannot open the file w_im_bin_flow5.txt\n");
		system("pause");
		exit(1);
	}
	for (n = 0;n < 1024;++n) {
		fscanf_s(in, "%s",s[n]);
	}
	fclose(in);
	err = fopen_s(&out, "E:\\Projects\\20151119-DFT-Project\\Debug\\dftout.dat", "w");
	if (err)
	{
		printf("Cannot open file dftout.dat\n");
		system("pause");
		exit(1);
	}
	for (k = 0;k < 1024; ++k) {
		for (i = 0;i < 10;++i) //将旋转因子输出来
		{
			k = k & 0x03ff;
			creal[9 - i] = ((k >> i) & 0x0001) + 48;
		}
		fprintf(out, "    10'b%s:ROM_data<=18'b%s\n", creal,s[i]);
	}
	fclose(out);
}*/
void testEE() {
	complex b1;
	complex b2;
	complex b3;
	int m = 0;
	char creal[9] = { 0 };
	char cimag[9] = { 0 };
	b1.real = 127;
	b1.imag = 0;
	b2.real = -127;
	b2.imag = 0;

	if (b1.real >= 0) {
		b1.real = b1.real & 0x00ff; //直接变成补码正数的补码是其本身
	}
	else {
		b1.real = (b1.real + 256) & 0x00ff;//将负数转化为补码再加上符号位
	}

	if (b2.real >= 0) {
		b2.real = b2.real & 0x00ff; //直接变成补码正数的补码是其本身
	}
	else {
		b2.real = (b2.real + 256) & 0x00ff;//将负数转化为补码再加上符号位
	}

	b3 = EE(b1, b2);
	for (m = 0;m < 8;++m) //将旋转因子输出来
	{
		creal[7 - m] = ((b3.real >> m) & 0x0001) + 48;
		cimag[7 - m] = ((b3.imag >> m) & 0x0001) + 48;
	}
	printf("%s%s\n", creal, cimag);
}
