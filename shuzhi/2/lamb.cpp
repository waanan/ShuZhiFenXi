#include <stdio.h>
#include <math.h>

#define E 1.0e-12         /*定义全局变量相对误差限*/

void chushi(double array[10][10])  /*矩阵初始化*/
{
	int i,j,a,b;
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			a=i+1;
			b=j+1;
			if(a!=b)
				array[i][j]=sin(0.5*a+0.2*b);
			else
				array[i][j]=1.52*cos(a+1.2*b);
		}
	}
}

void nssj(double array[10][10])  /*矩阵拟上三角化*/
{
  int i,j,k;
   double d,c,h,t,temp1,temp2,u[10]={0},p[10],q[10];
	for(i=0;i<8;i++)
	{
		d=0;
		for(j=i+2;j<10;j++)
			d+=fabs(array[j][i]);
        if(d<=E)
        continue;

		d=0;
		for(j=i+1;j<10;j++)
		{
			d=d+array[j][i]*array[j][i];
		}
		d=sqrt(d);
		if(array[i+1][i]>0)
			c=-1*d;
		else
			c=d;
		h=c*c-c*array[i+1][i];
		for(j=i+1;j<10;j++)
		{
			u[j]=array[j][i];
		}
		u[i+1]=u[i+1]-c;
		for(j=0;j<10;j++)
		{
			temp1=0;
			temp2=0;
			for(k=i+1;k<10;k++)
			{
				temp1+=array[k][j]*u[k];
				temp2+=array[j][k]*u[k];
			}
			p[j]=temp1/h;
			q[j]=temp2/h;
		}
		temp1=0;
		for(j=i+1;j<10;j++)
		{
			temp1+=p[j]*u[j];
		}
		t=temp1/h;
		for(j=i+1;j<10;j++)
		{
			q[j]=q[j]-t*u[j];
		}
		for(j=0;j<10;j++)
		{
			for(k=i+1;k<10;k++)
			{
				array[j][k]=array[j][k]-q[j]*u[k];
			}
		}
		for(j=i+1;j<10;j++)
		{
			for(k=0;k<10;k++)
			{
				array[j][k]=array[j][k]-u[j]*p[k];
			} 
		}
		for(j=i+2;j<10;j++)
		{
			array[j][i]=0;
		}
	}


	//打印An-1
/*	printf("矩阵A 经过拟上三角化后所得矩阵An-1:\n");
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			if(j%3==0) printf("\n");
			printf("a[%d][%d]=%.12e    ",i,j,array[i][j]);
		}
		printf("\n");
	} */
}

//对矩阵进行QR分解
void QR_M(double M[10][10],double Q[10][10],int n)
{
	int i,j,column,row;
	double d,c,h,u[10]={0},w[10]={0},p[10]={10};
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
			Q[i][j]=1;
			else Q[i][j]=0;
		}
	}
	for(column=0;column<n-1;column++)
	{
		d=0;
		for(row=column+1;row<n;row++)
			d+=fabs(M[row][column]);
        if(d<=E)
        continue;
		d=0;
		for(row=column;row<n;row++)
		{
			d=d+M[row][column]*M[row][column];
		}
		d=sqrt(d);
		if(M[column][column]>0)
			c=-1*d;
		else
			c=d;
		h=c*c-c*M[column][column];
    	for(row=column;row<n;row++)
		{
			u[row]=M[row][column];
		}
		u[column]=u[column]-c;
		for(i=0;i<n;i++)
		{
			w[i]=0;
			for(j=column;j<n;j++)
			{
				w[i]=w[i]+Q[i][j]*u[j];
			}
		}
		for(i=0;i<n;i++)
		{
			for(j=column;j<n;j++)
			{
			Q[i][j]=Q[i][j]-w[i]*u[j]/h;
			}
		}
		for(i=0;i<n;i++)
		{
			p[i]=0;
			for(j=column;j<n;j++)
			{
				p[i]=p[i]+M[j][i]*u[j];
			}
			p[i]=p[i]/h;
		}
		for(i=column;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				M[i][j]=M[i][j]-u[i]*p[j];
			}
		}
		for(j=column+1;j<n;j++)
		{
			M[j][column]=0;
		}
	}

}

//解方程，得到特征值
int function(double M[10][10],double s[2][2],int m)
{
	int flag=0;
	double s_1,s_2,delt;
	s_1=M[m-1][m-1]+M[m][m];
	s_2=M[m-1][m-1]*M[m][m]-M[m-1][m]*M[m][m-1];
	delt=s_1*s_1-4*s_2;
	if(delt>0)
	{
		flag=1;
		s[0][0]=0.5*(s_1+sqrt(delt));
		s[1][0]=0.5*(s_1-sqrt(delt));
	}
	else
	{
		s[0][0]=0.5*s_1;
		s[0][1]=0.5*sqrt(-delt);
		s[1][0]=0.5*s_1;
		s[1][1]=-0.5*sqrt(-delt);
		flag=0;
	}
	return flag;
}

//双步位移QR迭代过程
void QR_d(double M[10][10],int m)
{
	int i,j,k;
	double s,t;
	double b[10][10]={0},Q[10][10]={0};
	s=M[m-2][m-2]+M[m-1][m-1];
	t=M[m-2][m-2]*M[m-1][m-1]-M[m-1][m-2]*M[m-2][m-1];
	for(i=0;i<=m-1;i++)
	{	
		for(j=0;j<=m-1;j++)
			for(k=0;k<=m-1;k++)
				b[i][j]+=M[i][k]*M[k][j];
	}
	for(i=0;i<=m-1;i++)
		for(j=0;j<=m-1;j++)
			b[i][j]-=s*M[i][j];
	for(i=0;i<=m-1;i++)
		b[i][i]+=t;
	QR_M(b,Q,m);
	for(i=0;i<=m-1;i++)
	{	
		for(j=0;j<=m-1;j++)
		{
			b[i][j]=0;
			for(k=0;k<=m-1;k++)
				b[i][j]+=Q[k][i]*M[k][j];
		}
	}
	for(i=0;i<=m-1;i++)
	{	
		for(j=0;j<=m-1;j++)
		{
			M[i][j]=0;
			for(k=0;k<=m-1;k++)
				M[i][j]+=b[i][k]*Q[k][j];
		}
	}
}

int sbwy_QR(double M[10][10],double Real[10],double NReal[10][2]) //带双步位移的QR迭代法
{
	int i,m=10,cond=1,count=0,countn=0,k,len=10000;
	double s[2][2],Q[10][10];
	while(cond!=7)
	{
		switch(cond)
		{
		case 1:
			if(fabs(M[m-1][m-2])<=1e-12)
				{
					Real[count]=M[m-1][m-1];
					count++;
					m--;
					cond=2;
				}
			else
					cond=3;
					break;
		case 2:
			if(m==1)
			{
				Real[count]=M[0][0];
				count++;
				cond=7;
			}
			else
				cond=1;
			break;
		case 3:
			if(m==2)
			{
				i=function(M,s,m-1); //调用求二个特征值函数
				if(i==1)
				{
					Real[count]=s[0][0];
					Real[count+1]=s[1][0];
					count+=2;
				}
				else
				{
					NReal[countn][0]=s[0][0];
					NReal[countn][1]=s[0][1];
					NReal[countn+1][0]=s[1][0];
					NReal[countn+1][1]=s[1][1];
					countn+=2;
				}
				cond=7;
				}
			else
				cond=4;
				break;
		case 4:
			if(fabs(M[m-2][m-3])<=E)
			{
				i=function(M,s,m-1); //调用求二特征值函数
				if(i==1)
				{
					Real[count]=s[0][0];
					Real[count+1]=s[1][0];
					count+=2;
				}
				else
				{
					NReal[countn][0]=s[0][0];
					NReal[countn][1]=s[0][1];
					NReal[countn+1][0]=s[1][0];
					NReal[countn+1][1]=s[1][1];
					countn+=2;
				}
				m-=2;
				cond=2;
			}
			else
			cond=5;
			break;
		case 5:
			if(k==len)
				cond=7;
			else
				cond=6;
			break;
		case 6:
			QR_d(M,m); //计算A(K+1)
			k++;
			cond=1;
			break;
		case 7:
			cond=7;
			break;
		}
	}
	return count;
}

//求特征值向量
void tzxl(double M[10][10],double A[10][10],double Real[10],int count)
{
	int i,j,k;
	double T_M[10][10],Q[10][10];
	for(k=0;k<count;k++)
	{
		for(i=0;i<10;i++)
		{
			for(j=0;j<10;j++)
			{
				T_M[i][j]=M[i][j];
			}
		}
		for(i=0;i<10;i++)
		{
			T_M[i][i]-=Real[k];
		}
		QR_M(T_M,Q,10);
		A[k][9]=1;
		for(i=8;i>=0;i--)
		{
			A[k][i]=0;
			for(j=i+1;j<10;j++)
				A[k][i]-=A[k][j]*T_M[i][j];
			A[k][i]/=T_M[i][i];
		}
	}
}
/*----------------------------------------------------------------
main 函数：用带双步位移的QR 法计算矩阵的特征值，先将矩阵变换为拟上三角矩阵；
计算实特征值对应的特征向量，使用QR分解法解齐次方程组。
------------------------------------------------------------------*/

void main()
{
	int i,j,k,m;
	double t,M[10][10],Q[10][10],A[10][10];
	double Real[10],NReal[10][2];
	chushi(M);
	nssj(M);


	/*  QR迭代法
	for(m=0;m<150;m++)
	{
	QR_M(M,Q,10);
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			t=0;
			for(k=i;k<10;k++)
			{
				t=t+M[i][k]*Q[k][j];
			}
			A[i][j]=t;
		}
	}
	QR_M(A,Q,10);
	for(i=0;i<10;i++)
		for(j=0;j<10;j++)
		{
			t=0;
			for(k=i;k<10;k++)
			{
				t=t+A[i][k]*Q[k][j];
			}
			M[i][j]=t;
		}
	}
	printf("矩阵An-1 QR迭代完成后所得矩阵R:\n");
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			if(j%3==0) printf("\n");
			printf("a[%d][%d]=%.12e    ",i,j,A[i][j]);
		}
		printf("\n");
	}

	printf("矩阵An-1 QR迭代完成后所得矩阵Q:\n");
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			if(j%3==0) printf("\n");
			printf("a[%d][%d]=%.12e    ",i,j,Q[i][j]);
		}
		printf("\n");
	}

	printf("矩阵An-1 QR迭代完成 RQ相乘:\n");
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			if(j%3==0) printf("\n");
			printf("a[%d][%d]=%.12e    ",i,j,M[i][j]);
		}
		printf("\n");
	}  */

	//双步位移法
	k=sbwy_QR(M,Real,NReal);
	printf("矩阵A 的实特征值分别为:\n");
	for(i=0;i<k;i++)
		printf("%.12e\n",Real[i]); //输出实特征值
	printf("\n");
	printf("矩阵A 的复数特征值分别为:\n");
	for(i=0;i<10-k;i++)
		printf("%.12e\t%.12ei\n",NReal[i][0],NReal[i][1]);
	
	//求特征向量
	chushi(M);
	tzxl(M,A,Real,k);

	printf("对应的特征向两为:\n");
	for(i=0;i<k;i++)
	{
		printf("T%d=\n",i);
		for(j=0;j<10;j++)
		{
			if(j%3==0) printf("\n");
			printf("%.12e    ",A[i][j]);
		}
		printf("\n");
	}  
}