#include "stdafx.h"
#include<iostream>
#include <fstream>
#include<cmath>
#include"Classes.h"
#include"Functions.h"
#include <string>
using namespace std;
void main()
{
   /* int rgb = ((int)0<< 16 | (int)255 << 8 | (int)255);
	float frgb = *reinterpret_cast<float*>(&rgb);
	cout<<frgb;
	*/
	cout<<"请输入掩膜厚度\n";
	initialSystem();
	for (int i = 0; i < 100; i++)
	{
		XUNHUAN();
		cout << "第"<<i+1<<"次刻蚀/沉积循环完成\n";
	}
	ofstream outfile;
	string a = "result.txt";
	outfile.open(a);
	outfile << "# .PCD v0.7 - Point Cloud Data file format\nVERSION 0.7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH 1000000\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS 1000000\nDATA ascii\n";
	for (int z = 0; z < nHEIGHT; z++)
	{
		for (int y = 0; y < nWIDTH; y++)
		{
			for (int x = 0; x < nLONGTH; x++)
			{
				switch (aCA[z][y][x].nSpecie)
				{
				case Air:
					outfile << 0 << " " << 0<< " " << 0<< " ";
					outfile << 2.17624e-38;  //黄色
					break;
				case Si:
					outfile << (float)x << " " << (float)y << " " << (float)z << " ";
					outfile << 2.17624e-38;  //黄色
					break;
				case Surface:
					outfile << (float)x << " " << (float)y << " " << (float)z << " ";
					outfile << 9.18341e-41;  //青色
					break;
				case Mask:
					outfile << (float)x << " " << (float)y << " " << (float)z << " ";
					outfile <<2.18925e-38;   //红色
					break;
				}
				outfile << "\n";
			}
		}
		cout << "第" << z+1 << "层完成输出\n";
	}
	cout << "已完成输出\n";
	outfile.close();
}

