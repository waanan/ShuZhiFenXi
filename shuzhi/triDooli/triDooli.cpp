#include <stdio.h>
#include <malloc.h>
#include <math.h>

void dooli(float *p,float *pa,int n,int r,int s);

void main()
{
	int i,j,n,r,s;
    freopen("C:\\Users\\forln\\Desktop\\shuzhi\\triDooli.in","r",stdin);
	scanf("%d",&n);
	scanf("%d",&r);
	scanf("%d",&s);
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



		dooli(p,pa,n,r,s);       //Doolittle

        
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
		scanf("%d",&r);
    	scanf("%d",&s);
		free(p);
		free(pa);
	}
	printf("\n");
}


void dooli_b(float *p,float *pa,int n,int r,int s);
void dooli(float *p,float *pa,int n,int r,int s)
{
	int i,j,k,min,max;
	float temp;
	for(i=1;i<=r;i++)
		p[i*n]=p[i*n]/p[0];
	for(i=1;i<n-1;i++)
	{
		min=i+s;
		if(n<min)
			min=n;
		for(j=i;j<=min;j++)
		{   
			max=0;
			if((i-r)>max)
				max=(i-r);
			if((j-s)>max)
				max=j-s;
			temp=0;
			for(k=max;k<i;k++)
			{
				temp=temp+p[k+i*n]*p[j+k*n];
			}
			p[j+i*n]=p[j+i*n]-temp;
		}
		min=i+r;
		if(n<min)
			min=n;
		for(j=i+1;j<=min;j++)
		{
			max=0;
			if((i-r)>max)
				max=i-r;
			if((j-s)>max)
				max=j-s;
			temp=0;
			for(k=max;k<i;k++)
			{
				temp=temp+p[k+j*n]*p[i+k*n];
			}
			p[i+j*n]=(p[i+j*n]-temp)/p[i+i*n];
		}
	}
	temp=0;
	max=0;
	if((n-1-r)>max)
			max=n-1-r;
	if((n-1-s)>max)
			max=n-1-s;
	for(k=max;k<n-1;k++)
	{
		temp=temp+p[k+(n-1)*n]*p[n-1+k*n];
	}
	p[n*n-1]=p[n*n-1]-temp;
	dooli_b(p,pa,n,r,s);
}

void dooli_b(float *p,float *pa,int n,int r,int s)
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
