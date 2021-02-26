#include<stdio.h>
#include<pthread.h>

struct params {int na,nb; double sum;};

void* bar(void* arg){
	struct params * p = (struct params *)arg;
	int n1 = (*p).na;
	int n2 = (*p).nb;
	double s = 0;
	for(int i=n1;i<n2;i++) s+=1.0/i;
	p->sum=s;
	printf("sum from %i to %i = %g\n",n1,n2,s);
}

int main(){
	int N=(int)1e6;
	int n1=1,n2=N/3,n3=N/3*2;
	double sum1=0,sum2=0,sum3=0;
	pthread_t t1,t2,t3;
	struct params p1 = {.na=n1,.nb=n2,.sum=sum1};
	struct params p2 = {.na=n2,.nb=n3,.sum=sum2};
	struct params p3 = {.na=n3,.nb=N,.sum=sum3};
	pthread_create(&t1,NULL,bar,(void*)&p1);
	pthread_create(&t2,NULL,bar,(void*)&p2);
	pthread_create(&t3,NULL,bar,(void*)&p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	double sum=p1.sum+p2.sum+p3.sum;
	printf("sum = %g\n",sum);
return 0;
}
