#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>

#define E 1.0e-12         /*定义全局变量相对误差限*/

void assign(double array[5][501])     /*将矩阵A转存为数组C[5][501]*/
{
	int i,j,k;
	//所有元素归零
	for(i=0;i<=4;i++)
	{
		for(j=0;j<=500;j++)
			array[i][j]=0;
	}
	//第0,4行赋值
	for(i=2;i<=500;i++)
	{   k=500-i;
		array[0][i]=-0.064;
	    array[4][k]=-0.064;
	}
	//第1,3行赋值
	for(i=1;i<=500;i++)
	{   k=500-i;
		array[1][i]=0.16;
	    array[3][k]=0.16;
	}
    //第2行赋值
    for(j=0;j<=500;j++)
    {   k=j+1;
		array[2][j]=(1.64-0.024*k)*sin((double)(0.2*k))-0.64*exp((double)(0.1/k));}
}

int max2(int a,int b)	      /*求两个整型数最大值*/
{
	if(a>b)                 
	return a;
    return b;		                           
}

int max3(int a,int b,int c)      /*求三整型数最大值*/
{   int t;
	if(a>b)
	t=a;                             
	else t=b;
	if(t<c) t=c;
	return(t);
}

int min2(int a,int b)          /*求两个整型数最小值的子程序*/
{
	if(a>b)                 
	return b;
    return a;		   
}


double mifa(double u[501],double array[5][501],double p)  /*带原点平移的幂法*/
{
	int i,j,k=0;                              /* u[501]为初始迭代向量*/
	double a,b,c=0;             /* array[5][501]为矩阵A的转存矩阵*/
	double y[501];                                 /*p为平移量*/
                               /*选用第一种迭代格式*/
    while(1)
	{
	k++;
	a=0;
	b=0;
    //求ηk-1
    for(i=0;i<=500;i++)                           
	{
		a=a+u[i]*u[i];
	}
	a=sqrt(a);
   //求yk-1
	for(i=0;i<=500;i++)
	{
		y[i]=u[i]/a;
	}
    //求uk
    for(i=0;i<=500;i++)
	{
		u[i]=0;
		for(j=max2(i-2,0);j<=min2(i+2,500);j++)
		{
			u[i]+=array[i-j+2][j]*y[j];
		}
		u[i]=u[i]-p*y[i];                   /*引入平移量*/
	}
    //求βk
	for(i=0;i<=500;i++)
	{
		b+=y[i]*u[i];
	}
	if(fabs((b-c)/b)<=E)                   /*达到精度水平，迭代终止*/
		break;
	    c=b;
	}
	printf("迭代次数为:%d ",k);
	return (b+p);                        /*直接返回A的特征值*/    
}

void LU(double array[5][501])             /*对矩阵A进行Doolittle分解*/
{                                    /*矩阵A转存在C[5][501]中*/
	int j,k,t;           /*分解结果L，U分别存在C[5][501]的上半部与下半部*/
	for(k=0;k<=500;k++)
	{
		for(j=k;j<=min2((k+2),500);j++)
		{
			for(t=max3(0,k-2,j-2);t<=(k-1);t++)
			{
				array[k-j+2][j]-=array[k-t+2][t]*array[t-j+2][j];
			}
		}
		if(k<500)
		for(j=k+1;j<=min2((k+2),500);j++)
		{
			for(t=max3(0,k-2,j-2);t<=(k-1);t++)
			{
				array[j-k+2][k]-=array[j-t+2][t]*array[t-k+2][k];
			}
			array[j-k+2][k]=array[j-k+2][k]/array[2][k];
		}
	}
}

double rmifa(double u[501],double array[5][501],double p)  
{                                              /*带原点平移的反幂法*/
	int i,j;
	double a,b,c=0;
	double y[501];
    //引入平移量
	for(i=0;i<=500;i++)
	{
		array[2][i]-=p;
	}
    //先将矩阵Doolittle分解
	LU(array);

	while(1)
	{
	a=0;
	b=0;
   //求ηk-1
	for(i=0;i<=500;i++)
	{
		a=a+u[i]*u[i];
	}
	a=sqrt(a);
   //求yk-1
	for(i=0;i<=500;i++)
	{
		y[i]=u[i]/a;
	}
    //回带过程，求解uk
    for(i=0;i<=500;i++)
	{
		u[i]=y[i];
	}
    //Ly=b
	for(i=1;i<=500;i++)
	{
		for(j=max2(0,(i-2));j<=(i-1);j++)
		{
			u[i]-=array[i-j+2][j]*u[j];
		}
		
	}
   //Ux=y
	u[500]=u[500]/array[2][500];

	for(i=499;i>=0;i--)
	{
		for(j=i+1;j<=min2((i+2),500);j++)
		{
			u[i]-=array[i-j+2][j]*u[j];
		}

		u[i]=u[i]/array[2][i];
	}
    //求βk
	for(i=0;i<=500;i++)
	{
		b+=y[i]*u[i];
	}
	if(fabs((b-c)/b)<=E)                      /*达到精度要求，迭代终止*/
		break;
	    c=b;
	}
	return (p+(1/b));             /*直接返回距离原点P最接近的A的特征值*/

}


void chushi(double a[])                   /*用随机数为初始迭代向量赋值*/
{
	int i;
    srand((int)time(0));
	for(i=0;i<=500;i++)
	{
		a[i]=(10.0*rand()/RAND_MAX);       /*生成0~10的随机数*/
	}

} 



void main()
{
	int i;
	double a1,a501,a,b,as;
	printf("     《数值分析》计算实习题目第一题\n");
	printf("       SY1406227     王斐\n");
   double u[501]; 
   double MatrixC[5][501];
   assign(MatrixC);

   //用带原点平移的幂法求解λ1，λ501
	chushi(u);
	a=mifa(u,MatrixC,0);
	chushi(u);
	b=mifa(u,MatrixC,a);
	if(a<0)
	{
		a1=a;
		a501=b;
	}else
	{
		a1=b;
		a501=a;
	}
	printf("λ1=%.12e\n",a1);
	printf("λ501=%.12e\n",a501);

    //用反幂法求λs
    chushi(u);
    as=rmifa(u,MatrixC,0);
	printf("λs=%.12e\n",as);
    //用带原点平移的反幂法求λik
	for(i=1;i<=39;i++)
	{
		a=a1+(i*(a501-a1))/40;
		assign(MatrixC);
		chushi(u);
		b=rmifa(u,MatrixC,a);
		printf("与μ%02d=%+.12e最接近的特征值λi%02d=%+.12e\n",i,a,i,b);
	}
    //求A的条件数
	a=fabs((a1/as));
	printf("A的（谱范数）条件数cond<A>2=%.12e\n",a);
    //求detA
	assign(MatrixC);
	LU(MatrixC);
	a=1;
	for(i=0;i<=500;i++)
	{
		a*=MatrixC[2][i];
	}
	printf("行列式detA=%.12e\n",a);

   //测试对不同迭代初始向量对λmax计算结果迭代次数的影响。
	printf("改变迭代初始向量，λmax迭代次数的测试如下：\n");
    assign(MatrixC);
	for(i=0;i<=50;i++)
	{
    chushi(u);
	a=mifa(u,MatrixC,0);
    printf("%03dλmax=%+e ",i,a);
	if(((i+1)%2)==0)
	printf("\n");
	}
}