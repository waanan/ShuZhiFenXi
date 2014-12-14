#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>

#define E 1.0e-12         /*����ȫ�ֱ�����������*/

void assign(double array[5][501])     /*������Aת��Ϊ����C[5][501]*/
{
	int i,j,k;
	//����Ԫ�ع���
	for(i=0;i<=4;i++)
	{
		for(j=0;j<=500;j++)
			array[i][j]=0;
	}
	//��0,4�и�ֵ
	for(i=2;i<=500;i++)
	{   k=500-i;
		array[0][i]=-0.064;
	    array[4][k]=-0.064;
	}
	//��1,3�и�ֵ
	for(i=1;i<=500;i++)
	{   k=500-i;
		array[1][i]=0.16;
	    array[3][k]=0.16;
	}
    //��2�и�ֵ
    for(j=0;j<=500;j++)
    {   k=j+1;
		array[2][j]=(1.64-0.024*k)*sin((double)(0.2*k))-0.64*exp((double)(0.1/k));}
}

int max2(int a,int b)	      /*���������������ֵ*/
{
	if(a>b)                 
	return a;
    return b;		                           
}

int max3(int a,int b,int c)      /*�������������ֵ*/
{   int t;
	if(a>b)
	t=a;                             
	else t=b;
	if(t<c) t=c;
	return(t);
}

int min2(int a,int b)          /*��������������Сֵ���ӳ���*/
{
	if(a>b)                 
	return b;
    return a;		   
}


double mifa(double u[501],double array[5][501],double p)  /*��ԭ��ƽ�Ƶ��ݷ�*/
{
	int i,j,k=0;                              /* u[501]Ϊ��ʼ��������*/
	double a,b,c=0;             /* array[5][501]Ϊ����A��ת�����*/
	double y[501];                                 /*pΪƽ����*/
                               /*ѡ�õ�һ�ֵ�����ʽ*/
    while(1)
	{
	k++;
	a=0;
	b=0;
    //���k-1
    for(i=0;i<=500;i++)                           
	{
		a=a+u[i]*u[i];
	}
	a=sqrt(a);
   //��yk-1
	for(i=0;i<=500;i++)
	{
		y[i]=u[i]/a;
	}
    //��uk
    for(i=0;i<=500;i++)
	{
		u[i]=0;
		for(j=max2(i-2,0);j<=min2(i+2,500);j++)
		{
			u[i]+=array[i-j+2][j]*y[j];
		}
		u[i]=u[i]-p*y[i];                   /*����ƽ����*/
	}
    //���k
	for(i=0;i<=500;i++)
	{
		b+=y[i]*u[i];
	}
	if(fabs((b-c)/b)<=E)                   /*�ﵽ����ˮƽ��������ֹ*/
		break;
	    c=b;
	}
	printf("��������Ϊ:%d ",k);
	return (b+p);                        /*ֱ�ӷ���A������ֵ*/    
}

void LU(double array[5][501])             /*�Ծ���A����Doolittle�ֽ�*/
{                                    /*����Aת����C[5][501]��*/
	int j,k,t;           /*�ֽ���L��U�ֱ����C[5][501]���ϰ벿���°벿*/
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
{                                              /*��ԭ��ƽ�Ƶķ��ݷ�*/
	int i,j;
	double a,b,c=0;
	double y[501];
    //����ƽ����
	for(i=0;i<=500;i++)
	{
		array[2][i]-=p;
	}
    //�Ƚ�����Doolittle�ֽ�
	LU(array);

	while(1)
	{
	a=0;
	b=0;
   //���k-1
	for(i=0;i<=500;i++)
	{
		a=a+u[i]*u[i];
	}
	a=sqrt(a);
   //��yk-1
	for(i=0;i<=500;i++)
	{
		y[i]=u[i]/a;
	}
    //�ش����̣����uk
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
    //���k
	for(i=0;i<=500;i++)
	{
		b+=y[i]*u[i];
	}
	if(fabs((b-c)/b)<=E)                      /*�ﵽ����Ҫ�󣬵�����ֹ*/
		break;
	    c=b;
	}
	return (p+(1/b));             /*ֱ�ӷ��ؾ���ԭ��P��ӽ���A������ֵ*/

}


void chushi(double a[])                   /*�������Ϊ��ʼ����������ֵ*/
{
	int i;
    srand((int)time(0));
	for(i=0;i<=500;i++)
	{
		a[i]=(10.0*rand()/RAND_MAX);       /*����0~10�������*/
	}

} 



void main()
{
	int i;
	double a1,a501,a,b,as;
	printf("     ����ֵ����������ʵϰ��Ŀ��һ��\n");
	printf("       SY1406227     ���\n");
   double u[501]; 
   double MatrixC[5][501];
   assign(MatrixC);

   //�ô�ԭ��ƽ�Ƶ��ݷ�����1����501
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
	printf("��1=%.12e\n",a1);
	printf("��501=%.12e\n",a501);

    //�÷��ݷ����s
    chushi(u);
    as=rmifa(u,MatrixC,0);
	printf("��s=%.12e\n",as);
    //�ô�ԭ��ƽ�Ƶķ��ݷ����ik
	for(i=1;i<=39;i++)
	{
		a=a1+(i*(a501-a1))/40;
		assign(MatrixC);
		chushi(u);
		b=rmifa(u,MatrixC,a);
		printf("���%02d=%+.12e��ӽ�������ֵ��i%02d=%+.12e\n",i,a,i,b);
	}
    //��A��������
	a=fabs((a1/as));
	printf("A�ģ��׷�����������cond<A>2=%.12e\n",a);
    //��detA
	assign(MatrixC);
	LU(MatrixC);
	a=1;
	for(i=0;i<=500;i++)
	{
		a*=MatrixC[2][i];
	}
	printf("����ʽdetA=%.12e\n",a);

   //���ԶԲ�ͬ������ʼ�����Ԧ�max����������������Ӱ�졣
	printf("�ı������ʼ��������max���������Ĳ������£�\n");
    assign(MatrixC);
	for(i=0;i<=50;i++)
	{
    chushi(u);
	a=mifa(u,MatrixC,0);
    printf("%03d��max=%+e ",i,a);
	if(((i+1)%2)==0)
	printf("\n");
	}
}