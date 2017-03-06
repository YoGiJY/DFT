/******************************************************************
 * @file         main.c
 * @date         2015-11-19 22:32:08
 * @author       Yao
 * @e-mail       yao.jiang@tongji.edu.cn
 * @brief        Discrete Fourier Transfrom
 *
 *****************************************************************/

#include "dft.h"
int main() {
	extern long npt;
	extern int inv;
	//testEE();
    // ����c++���Ե����벢�ҺͲ����� [1/1/2016 YoGiJY]
	sin512();
	// ������������ [1/1/2016 YoGiJY]
	initdata(4);
	fft_dit();
	save_data();
	system("pause");
	return 0;
}