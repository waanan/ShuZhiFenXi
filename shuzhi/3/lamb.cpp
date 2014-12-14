#include <stdio.h>
#include <math.h>

#define E 1.0e-12         /*定义全局变量相对误差限*/
#define E1 1.0e-7       /*定义全局变量相对误差限*/

//t u z 之间的数表
double tu_Z[6][6]={
	{-0.5, -0.34, 0.14, 0.94, 2.06, 3.5},
	{-0.42, -0.5, -0.26, 0.3, 1.18, 2.38},
	{-0.18, -0.5, -0.5, -0.18, 0.46, 1.42},
	{0.22, -0.34, -0.58, -0.5, -0.1, 0.62},
	{0.78, -0.02, -0.5, -0.66, -0.5, -0.02},
	{1.5, 0.46, -0.26, -0.66, -0.74, -0.5}
};

//列主元高斯消去法
void guass(double M[4][4],double A[4],int n)
{
	int i,j,k,Max;
	double Temp;
	for(i=0;i<n-1;i++)
	{
		Max= i;
		for(j=i;j<n;j++)
		{
			if(fabs(M[j][i])>fabs(M[Max][i]))
				Max=j;
		}
		if(Max!=i)
		{
			for(j=i;j<n;j++)
			{
				Temp=M[Max][j];
				M[Max][j]=M[i][j];
				M[i][j]=Temp;
			}
				Temp=A[Max];
				A[Max]=A[i];
				A[i]=Temp;
		}
		for(j=i+1;j<n;j++)
		{
			Temp=M[j][i]/M[i][i];
			M[j][i]=0;
			for(k=i+1;k<n;k++)
			{
				M[j][k]-=Temp*M[i][k];
			}
			A[j]-=Temp*A[i];
		}
	}
	A[n-1]/=M[n-1][n-1];
	for(i=n-2;i>=0;i--)
	{
		Temp=0;
		for(j=i+1;j<n;j++)
		{
			Temp+=A[j]*M[i][j];
		}
		A[i]-=Temp;
		A[i]/=M[i][i];
	}
}

//向量的无穷范数
double infi(double A[4])
{
	int i;
	double temp;
	temp=A[0];
	for(i=1;i<4;i++)
	{
		if(fabs(A[i])>fabs(temp))
			temp=A[i];
	}
	return fabs(temp);
}

//牛顿迭代法求非线性方程组的解 A[0~3]分别为t,u,v,w
void newton(double x,double y,double A[4])
{
	int i,j=0;
	double A1[4][4],detA[4];
	for(i=0;i<4;i++)
	{
		A[i]=0.5;
		detA[i]=0.5;
	}
	while(infi(detA)/infi(A)>E)
	{
		j++;
		if(j>400)
		{
			printf("牛顿迭代法计算失败");
			break;
		}
		A1[0][0]=-0.5*sin(A[0]);
		A1[0][1]=1.0;
		A1[0][2]=1.0;
		A1[0][3]=1.0;
		A1[1][0]=1.0;
		A1[1][1]=0.5*cos(A[1]);
		A1[1][2]=1.0;
		A1[1][3]=1.0;
		A1[2][0]=0.5;
		A1[2][1]=1.0;
		A1[2][2]=-sin(A[2]);
		A1[2][3]=1.0;
		A1[3][0]=1.0;
		A1[3][1]=0.5;
		A1[3][2]=1.0;
		A1[3][3]=cos(A[3]);
		detA[0]=-1.0*(0.5*cos(A[0])+A[1]+A[2]+A[3]-x-2.67);
		detA[1]=-1.0*(A[0]+0.5*sin(A[1])+A[2]+A[3]-y-1.07);
		detA[2]=-1.0*(0.5*A[0]+A[1]+cos(A[2])+A[3]-x-3.74);
		detA[3]=-1.0*(A[0]+0.5*A[1]+A[2]+sin(A[3])-y-0.79);
		guass(A1,detA,4);
	for(i=0;i<4;i++)
	{
		A[i]+=detA[i];
	}
	}
}

//插值法中间计算
double cha(double x,double y)
{
	int m,n,k,l,a,b;
	double t,u,tu[4],temp,su;
	newton(x,y,tu);
	t=tu[0];
	u=tu[1];
	m=(int)(t/0.2+0.5);
	n=(int)(u/0.4+0.5);
	if(m==0)  m=1;
	if(m==5)  m=4;
	if(n==0)  n=1;
    if(n==5)  n=4;
	temp=0.0;
	for(k=m-1;k<=m+1;k++)
	{
		for(l=n-1;l<=n+1;l++)
		{
		su=tu_Z[k][l];
		for(a=m-1;a<=m+1;a++)
		{
			if(a!=k)
			su*=(t-a*0.2)/(0.2*(k-a));
		}
		for(b=n-1;b<=n+1;b++)
		{
			if(b!=l)
				su*=(u-b*0.4)/(0.4*(l-b));
		}
		temp+=su;
		}
	}
	return temp;
}

//插值法求 x y z之间的对应关系
void chazhi(double xy_Z[11][21])
{
	int i,j;
	double x,y;
	for(i=0;i<11;i++)
		for(j=0;j<21;j++)
		{
			x=0.08*i;
			y=0.5+0.05*j;
			xy_Z[i][j]=cha(x,y);
			printf("x%d:%.2e    y%d:%.2e    f:%.12e\n",i,x,j,y,xy_Z[i][j]);
		}
}

//矩阵求逆
void reverse(double M[11][11],int n)
{
	int i,j,k,Max;
	double Temp;
	double T[21][21];
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			if(i!=j)
				T[i][j]=0;
			else
				T[i][j]=1.0;
		}


	for(i=0;i<n-1;i++)
	{
		Max= i;
		for(j=i;j<n;j++)
		{
			if(fabs(M[j][i])>fabs(M[Max][i]))
				Max=j;
		}
		if(Max!=i)
		{
			for(j=i;j<n;j++)
			{
				Temp=M[Max][j];
				M[Max][j]=M[i][j];
				M[i][j]=Temp;
			}
			for(j=0;j<n;j++)
			{
				Temp=T[Max][j];
				T[Max][j]=T[i][j];
				T[i][j]=Temp;
			}
		}
		for(j=i+1;j<n;j++)
		{
			Temp=M[j][i]/M[i][i];
			M[j][i]=0;
			for(k=i+1;k<n;k++)
			{
				M[j][k]-=Temp*M[i][k];
			}
			for(k=0;k<n;k++)
			{
				T[j][k]-=Temp*T[i][k];
			}
		}
	}

	for(i=n-1;i>=1;i--)
	{
		for(j=i-1;j>=0;j--)
		{
			Temp=M[j][i]/M[i][i];
			for(k=0;k<n;k++)
			{
				T[j][k]-=Temp*T[i][k];
			}
		}
	}

	for(i=0;i<n;i++)
	{
		Temp=M[i][i];
		for(j=0;j<n;j++)
		{
			M[i][j]=T[i][j]/Temp;
		}
	}
}

//矩阵拟合
double nihe(double xy_Z[11][21],double C[11][11],int n)
{
	int i,j,k,l;
	double det,m,B[11][11],G[21][21],T[11][11],Temp[11][11];
	for(i=0;i<11;i++)
	{
		for(j=0;j<n;j++)
		{
			B[i][j]=pow(0.08*i,j);
		}
	}
	for(i=0;i<21;i++)
	{
		for(j=0;j<n;j++)
		{
			G[i][j]=pow(0.5+0.05*i,j);
		}
	}
//Temp=UG 11*n
	for(i=0;i<11;i++)
	{
		for(j=0;j<n;j++)
		{
			m=0.0;
			for(k=0;k<21;k++)
				m+=xy_Z[i][k]*G[k][j];
			Temp[i][j]=m;
		}
	}
//T=GtG n*n
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			m=0.0;
			for(k=0;k<21;k++)
			{	m+=G[k][i]*G[k][j];}
			T[i][j]=m;
		}
	}
//T=GTG-1
	reverse(T,n);
//UG(GtG)-1
	for(i=0;i<11;i++)
	{	for(j=0;j<n;j++)
		{
			m=0.0;
			for(k=0;k<n;k++)
			{	m+=Temp[i][k]*T[k][j];}
			C[i][j]=m;
	}}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			m=0.0;
			for(k=0;k<11;k++)
				m+=B[k][i]*C[k][j];
			Temp[i][j]=m;
		}
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			m=0.0;
			for(k=0;k<11;k++)
				m+=B[k][i]*B[k][j];
			T[i][j]=m;
		}
	}
	reverse(T,n);
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			m=0.0;
			for(k=0;k<n;k++)
				m+=T[i][k]*Temp[k][j];
			C[i][j]=m;
		}
	}
//求误差
	det=0.0;
	for(i=0;i<11;i++)
	{	for(j=0;j<21;j++)
		{
			m=0.0;
			for(k=0;k<n;k++){
				for(l=0;l<n;l++)
				{
					m+=C[k][l]*pow(0.08*i,k)*pow(0.5+0.05*j,l);
			}}
			det+=pow(xy_Z[i][j]-m,2);
		}
	}
	return det;
}

void main()
{
	//freopen("answer.out","w",stdout);
	int i,j,k,m,n;
	double x,y,t1,t2,xy_Z[11][21];
	double det,C[11][11];
	chazhi(xy_Z);
	for(k=1;k<11;k++)
	{
		det=nihe(xy_Z,C,k+1);
		printf("k的值为  %d\n",k);
		printf("σ的值为  %.12e\n",det);
		if(det<E1)
		{
			printf("拟合过程结果如下:\n");
			printf("k的值为  %d\n",k);
			printf("σ的值为  %.12e\n",det);
			break;
		}
	}
	if(k>10)
		printf("拟合过程失败！\n");
	printf("拟合系数为:\n");
	for(i=0;i<k+1;i++)
		for(j=0;j<k+1;j++)
			printf("C %d %d为  %.12e\n",i,j,C[i][j]);

	for(i=1;i<=8;i++)
		for(j=1;j<=5;j++)
		{
			x=0.1*i;
			y=0.5+0.2*j;
			t1=cha(x,y);
			t2=0.0;
			for(m=0;m<k+1;m++){
				for(n=0;n<k+1;n++)
				{
					t2+=C[m][n]*pow(x,m)*pow(y,n);
				}}
			printf("x*:%.3e  y*:%.3e    f*:%.12e  p*:%.12e\n",x,y,t1,t2);

		}
			

  }