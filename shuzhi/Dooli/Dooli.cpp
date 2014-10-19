#include <stdio.h>
#include <malloc.h>
#include <math.h>

void dooli_1(float *p,float *pa,int n);
void dooli_2(float *p,float *pa,int n);

void main()
{
	int i,j,n;
    freopen("C:\\Users\\forln\\Desktop\\shuzhi\\Dooli.in","r",stdin);
	scanf("%d",&n);
	while(n>0)
	{
		printf("当前矩阵为 %d 维。\n",n);
		float *p=(float*)malloc(sizeof(float)*n*n); 
		float *pa=(float*)malloc(sizeof(float)*n);
		float an;
        for(i=0;i<n*n;i++)
			scanf("%f",&p[i]);

		for(i=0;i<n;i++)
			scanf("%f",&pa[i]); 



		//dooli_1(p,pa,n);       //Doolittle
		dooli_2(p,pa,n);         //guass列主元素消去法

        
		for(i=0;i<n;i++)
		{for(j=0;j<n;j++)
			{
				printf("%f  ",p[j+i*n]);
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


void dooli_b(float *p,float *pa,int n);
void ch_b(float *pa,int *mx,int n);
void dooli_1(float *p,float *pa,int n)
{
	int i,j,k;
	float temp;
	for(i=1;i<n;i++)
		p[i*n]=p[i*n]/p[0];
	for(i=1;i<n-1;i++)
	{
		for(j=i;j<n;j++)
		{   
			temp=0;
			for(k=0;k<i;k++)
			{
				temp=temp+p[k+i*n]*p[j+k*n];
			}
			p[j+i*n]=p[j+i*n]-temp;
		}
		for(j=i+1;j<n;j++)
		{
			temp=0;
			for(k=0;k<i;k++)
			{
				temp=temp+p[k+j*n]*p[i+k*n];
			}
			p[i+j*n]=(p[i+j*n]-temp)/p[i+i*n];
		}
	}
	temp=0;
	for(k=0;k<n-1;k++)
	{
		temp=temp+p[k+(n-1)*n]*p[n-1+k*n];
	}
	p[n*n-1]=p[n*n-1]-temp;
	dooli_b(p,pa,n);
}
void dooli_2(float *p,float *pa,int n)
{
	int i,j,k,max_line;
	int *mx=(int*)malloc(sizeof(int)*(n-1));
	float temp;
	max_line=0;
	for(i=1;i<n;i++)
	{
		if(fabs(p[i*n])>fabs(p[max_line*n]))
			max_line=i;
	}
	mx[0]=max_line;
	if(max_line!=0)
		for(i=0;i<n;i++)
		{
			temp=p[i+max_line*n];
			p[i+max_line*n]=p[i];
			p[i]=temp;
		}

	for(i=1;i<n;i++)
		p[i*n]=p[i*n]/p[0];
	for(i=1;i<n-1;i++)
	{

		for(j=i;j<n;j++)
		{
			temp=0;
			for(k=0;k<i;k++)
			{
				temp=temp+p[k+j*n]*p[i+k*n];
			}
			p[i+j*n]=p[i+j*n]-temp;
		}
			
		max_line=i;
		for(j=i+1;j<n;j++)
		{
			if(fabs(p[i+j*n])>fabs(p[i+max_line*n]))
				max_line=j;
		}
		mx[i]=max_line;
		if(max_line!=i)
		for(j=0;j<n;j++)
		{
			temp=p[j+max_line*n];
			p[j+max_line*n]=p[j+i*n];
			p[j+i*n]=temp;
		}

		for(j=i+1;j<n;j++)
		{
			p[i+j*n]=p[i+j*n]/p[i+i*n];
			temp=0;
			for(k=0;k<i;k++)
			{
				temp=temp+p[k+i*n]*p[j+k*n];
			}
			p[j+i*n]=p[j+i*n]-temp;
		}
	}
	temp=0;
	for(k=0;k<n-1;k++)
	{
		temp=temp+p[k+(n-1)*n]*p[n-1+k*n];
	}
	p[n*n-1]=p[n*n-1]-temp;
	ch_b(pa,mx,n);
	free(mx);
	dooli_b(p,pa,n);
}
void dooli_b(float *p,float *pa,int n)
{
	int i,j;
	float temp;
	for(i=1;i<n;i++)
	{ 
		temp=0;
		for(j=0;j<i;j++)
			temp=temp+pa[j]*p[j+i*n];
		pa[i]=pa[i]-temp;
	}
	pa[n-1]=pa[n-1]/p[n*n-1];
	for(i=n-2;i>=0;i--)
	{
		temp=0;
		for(j=i+1;j<n;j++)
		temp=temp+pa[j]*p[j+i*n];
		pa[i]=(pa[i]-temp)/p[i+i*n];
	}
}
void ch_b(float *pa,int *mx,int n)
{
	int i;
	float temp;
	for(i=0;i<n-1;i++)
	{
		if(mx[i]!=i)
		{
			temp=pa[mx[i]];
			pa[mx[i]]=pa[i];
			pa[i]=temp;
		}
	}
}