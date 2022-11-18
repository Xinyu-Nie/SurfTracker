#include <stdio.h>
#include <unistd.h>
//#include <unistd.h>
 
int main()
 {
	int i = 0;
	char bar[102];
	const char *lable = "|/-\\";
	bar[0] = 0;
	int len = 555;
	int cnt = 0;
	while (i <= len)
	{
		printf("[%-100s][%d%%][%c]\r", bar, cnt, lable[i%4]);
		fflush(stdout);
		bar[cnt] = '#';
		if (i*1.0/len*100>cnt) cnt++;
		i++;
		bar[cnt] = 0;
		usleep(10000);
	}
	printf("\n");
	return 0;
}