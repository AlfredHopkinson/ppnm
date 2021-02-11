#include <stdio.h>
#include "hello.h"

void world(void);
char hello[] = "file scope";

void bar(void){
	char hello[] = "function scope";	
	printf("hello (function scope?):%s\n", hello);
	{
		char hello[] = "block scope";	
		printf("hello (block scope?):%s\n", hello);
	}
}

int main(){
	printf("hello, world\n");
	world();
	printf("hello (file scope?):%s\n", hello);
	bar();
return 0;
}

