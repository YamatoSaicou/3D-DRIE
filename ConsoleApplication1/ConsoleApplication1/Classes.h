#pragma once
#define Si 4
#define Mask 2
#define Air 3
#define Full 0.5
class cCell {
public:
	double dAmount;  //Ԫ���������������������
	int nSpecie;
	cCell()
	{
		nSpecie = Si; dAmount = Full;
	}  //Ĭ������ΪSi��Ĭ������ΪFull
	cCell(int a, double b)
	{
		nSpecie = a; dAmount = b;
	}  //���캯��
};
