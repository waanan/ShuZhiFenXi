#include <stdio.h>
#include <malloc.h>
#include <math.h>

void guass_1(float *p,float *pa,int n);
void guass_2(float *p,float *pa,int n);

void main()
{
	int i,j,n;
    freopen("C:\\Users\\forln\\Desktop\\shuzhi\\number.in","r",stdin);
	scanf("%d",&n);
	while(n>0)
	{
		printf("当前矩阵为 %d 维。\n",n);
		float *p=(float*)malloc(sizeof(float)*n*(n+1)); 
		float *pa=(float*)malloc(sizeof(float)*n);
		float an;
        for(i=0;i<n*(n+1);i++)
			scanf("%f",&p[i]);

		guass_1(p,pa,n);       //guass消去法
		//guass_2(p,pa,n);         //guass列主元素消去法

        
		for(i=0;i<n;i++)
		{for(j=0;j<n+1;j++)
			{
				printf("%f  ",p[j+i*(n+1)]);
			}	
		printf("\n");
		}

		printf("计算所得结果为:");
		for(i=0;i<n;i++)
		printf(" %f ",pa[i]);
		printf("\n答案为：");
		for(i=0;i<n;i++)
		{
			scanf("%f",&an);
			printf(" %f ",an);
		}
		printf("\n\n");

		scanf("%d",&n);
		free(p);
		free(pa);
	}
	printf("\n");
}


void guass_b(float *p,float *pa,int n);
void guass_1(float *p,float *pa,int n)
{
	int i,j,k;
	int m=n+1;
	float Mik;
	for(i=0;i<n-1;i++)
	{

	 for(j=i+1;j<n;j++)
		{
			Mik=p[i+j*m]/p[i+i*m];
			printf("MIK %d %d %f\n",j+1,i+1,Mik);
			for(k=i+1;k<m;k++)
			{
				p[k+j*m]=p[k+j*m]-p[k+i*m]*Mik;
			}
		}
	 }
	guass_b(p,pa,n);
}
void guass_2(float *p,float *pa,int n)
{
	int i,j,k,max_line;
	int m=n+1;
	float Mik;
	for(i=0;i<n-1;i++)
	{
     max_line=i;
	 for(j=i+1;j<n;j++)
	 {
		 if(fabs(p[i+j*m])>fabs(p[i+max_line*m]))
			 max_line=j;
	 }
	 if(max_line!=i)
		 for(j=i;j<m;j++)
		 {
         Mik=p[j+i*m];
		 p[j+i*m]=p[j+max_line*m];
		 p[j+max_line*m]=Mik;
		 }
	 for(j=i+1;j<n;j++)
		{
			Mik=p[i+j*m]/p[i+i*m];
			printf("MIK %d %d %f\n",j+1,i+1,Mik);
			for(k=i+1;k<m;k++)
			{
				p[k+j*m]=p[k+j*m]-p[k+i*m]*Mik;
			}
		}
	 }
	guass_b(p,pa,n);
}
void guass_b(float *p,float *pa,int n)
{
	int m= n+1,i,j;
	float temp;
	pa[n-1]=p[n*m-1]/p[n*m-2];
	for(i=n-2;i>=0;i--)
	{   
		temp=p[(i+1)*m-1];
		for(j=m-2;j>i;j--)
		temp=temp-p[j+i*m]*pa[j];
		pa[i]=temp/p[i+i*m];
	}
}