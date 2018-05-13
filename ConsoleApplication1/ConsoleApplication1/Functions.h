#pragma once
#include<iostream>
#include<math.h>
#define nHEIGHT 100
#define nWIDTH 100
#define Mask 2
#define Air 3
#define Surface 1
#define N 1000
#define k0  1.68*pow(10, -15)
#define Fr  6 * pow(10, 14)
#define Ee  0.1 * pow(10, -19)
#define kb  1.38065*pow(10, -23)
#define Ts  300
#define c  2 * pow(10, -23)
#define Yast  0.0032*pow(10,8)
#define Ji  4 * pow(10,16)
#define kD 0.3
#define shigema 0.2
#define nLONGTH 100
//取k0=1.68*10（-15），Fr=5*10^14,Ee=10^-20,Kb=1.38*10^23,Ts=300. 离子辅助刻蚀率YAST Ji离子流量
using namespace std;
int nMASKTHICK;
int nWINKUAN, nWINCHANG;
cCell aCA[nHEIGHT][nWIDTH][nLONGTH];  
double f(double seta,double y,double x) 
{
	return (cos(fabs(atan(y / x))-seta))*exp(-(seta-1.57)*(seta-1.57)/(2* shigema*shigema))/(shigema*sqrt(2*3.14));
}

void initialSystem()
{
	//一开始做长方形开口掩膜
	cin >> nMASKTHICK;
	cout << "请输入掩膜的宽度\n";
	cin >> nWINKUAN;
	cout << "请输入掩膜的长度\n";
	cin >> nWINCHANG;
	
	for (int z = 0; z < nMASKTHICK; z++)  // nMASKTHICK为掩膜厚度
	{
		for (int y = 0; y < nWIDTH; y++)
		{
			for (int x = 0; x < nLONGTH; x++)
			{
				aCA[z][y][x].dAmount = 0;
				aCA[z][y][x].nSpecie = Mask; 
			}
		}
	}  //将掩膜厚度以内的元胞设置为Mask
	for (int z = 0; z <= nMASKTHICK-1; z++)
	{
		for (int y = (nWIDTH-nWINKUAN)/2-1; y < (nWIDTH +nWINKUAN) / 2-1; y++) // 窗口参数
		{
			for (int x =0 ; x <= nWINCHANG - 1; x++)
			{
				aCA[z][y][x].nSpecie = Air;   aCA[z][y][x].dAmount = 0;  //将窗口位置的元胞设置为Air
			}
		}
	}
	for (int y = (nWIDTH - nWINKUAN) / 2 - 1; y < (nWIDTH + nWINKUAN) / 2 - 1; y++) 
		{
		for (int x = 0; x <= nWINCHANG - 1; x++)
			{
				aCA[nMASKTHICK][y][x].nSpecie = Surface;   aCA[nMASKTHICK][y][x].dAmount = Full;  //将窗口位置下面的元胞设置为Surface
			}
		}
	cout << "元胞系统初始化完成，刻蚀开始\n";
}
void XUNHUAN()
{
	for (int z = nMASKTHICK - 1; z < nHEIGHT;z++)
	{
		for (int y = 0; y < nWIDTH; y++)  // 求出表面元胞法向量（全部的细胞扫描一次的方法来找到表面元胞细胞）
		{
			for (int x = 0; x < nLONGTH; x++)
			{
				if (aCA[z][y][x].nSpecie == Surface)
				{
					double dNormal[4] = { 0 };
					int aX[2] = { 0,0 }, aY[2] = { 0,0 }, aZ[2] = { 0,0 };  //分别定义X、Y、Z轴的邻居
					if (aCA[z][y][x - 1].nSpecie==Air)aX[0] = 1;  //X轴负方向邻居状态
					if (aCA[z][y][x + 1].nSpecie == Air)aX[1] = 1;  //X轴正方向邻居状态
					if (aCA[z][y - 1][x].nSpecie == Air)aY[0] = 1;
					if (aCA[z][y + 1][x].nSpecie == Air)aY[1] = 1;
					if (aCA[z - 1][y][x].nSpecie == Air)aZ[0] = 1;
					if (aCA[z + 1][y][x].nSpecie == Air)aZ[1] = 1;
					dNormal[0] = aX[0] - aX[1];  //X轴由正负方向的邻居状态得出的向量
					dNormal[1] = aY[0] - aY[1];
					dNormal[2] = aZ[0] - aZ[1];
					//得到法向量的三个分量
					//接下来计算刻蚀速率
					dNormal[3] = (fabs(dNormal[0]) + fabs(dNormal[1]) + fabs(dNormal[2])) / sqrt(dNormal[0] * dNormal[0] + dNormal[1] * dNormal[1] + dNormal[2] * dNormal[2]);
					double v1 = 0.0000; double v2 = 0;
					double seta1 = 0, seta2 = 0;
					double m = 0;
					m = k0 * Fr*exp(-Ee / (kb*Ts));
					v1 = m * (dNormal[3]);   //化学刻蚀速率 
					double beta = 1.2*(double)nWINKUAN/(double)nWIDTH;
					v2 = v1*beta;
					double v = v1 + v2;
					aCA[z][y][x].dAmount = aCA[z][y][x].dAmount - v;
					if (aCA[z][y][x].dAmount <= 0)
					{
						aCA[z][y][x].nSpecie = Air;
						aCA[z][y][x].dAmount = 0;
						if (aCA[z][y + 1][x].nSpecie == Si)
						{
							aCA[z][y + 1][x].nSpecie = Surface;
						}
						if (aCA[z][y - 1][x].nSpecie == Si)
						{
							aCA[z][y - 1][x].nSpecie = Surface;
						}
						if (aCA[z][y][x - 1].nSpecie == Si)
						{
							aCA[z][y][x - 1].nSpecie = Surface;
						}
						if (aCA[z][y][x + 1].nSpecie == Si)
						{
							aCA[z][y][x + 1].nSpecie = Surface;
						}
						if (aCA[z+1][y][x].nSpecie == Si)
						{
							aCA[z+1][y][x].nSpecie = Surface;
						}
						if (aCA[z-1][y][x].nSpecie == Si)
						{
							aCA[z-1][y][x].nSpecie = Surface;
						}
					}
					if (aCA[z][y][x].dAmount > 0)
					{
						aCA[z][y][x].dAmount = aCA[z][y][x].dAmount + v1 * kD;
						/*cout << "坐标是是(" << x << "," << y << ")" << "的沉积后含量是" << aCA[y][x].dAmount << "\n";*/
					}
				}
			}
		}
	}
}