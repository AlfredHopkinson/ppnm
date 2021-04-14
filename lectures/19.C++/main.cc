#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<gc.h> //sudo apt install libgc-dev
#define SELF  (*this)
#define FOR(i) for(int i=0;i<size;i++)


struct vec {
	int size; double* data;
	
	vec(void):size(0),data(nullptr){}// default constructor

	vec(int n){ // parametrized constructor
		size=n;
		data=(double*)GC_MALLOC_ATOMIC(n*sizeof(double));
		}

	vec(vec& b){//copy constructor
		size=b.size;
		data=(double*)GC_MALLOC_ATOMIC(size*sizeof(double));
		FOR(i)SELF[i]=b[i];
		}

	vec(vec&& b){//move constructor
		size=b.size;
		data=(double*)GC_MALLOC_ATOMIC(size*sizeof(double));
		FOR(i)SELF[i]=b[i];
		}
					
	~vec(){/* GC_FREE(data);*/} // destructor

	double& operator[](int i){return data[i];} // v[i]=x;
	double& operator()(int i){return data[i];} // v(i)=x;
	void operator*=(double x){FOR(i)SELF[i]*=x;}
	vec operator+(vec b){ // a+b => (a).operator+(b)
		assert(size==b.size);
		vec v(size);
		FOR(i)v[i]=SELF[i]+b[i];
		return v;
		}
	void print(const char* s){
		printf("%s\n",s);
		FOR(i)printf("%9.3g ",SELF[i]);
		printf("\n");
	}
};

int main(){
	int n=3;
	vec v(n);
	vec u,w;
	v[1]=9;w[0]=7;
	u=v;
	u.print("u=");
return 0;
}
	
