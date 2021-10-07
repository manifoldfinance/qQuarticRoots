all: k.h
	$(CC) quartic.c -D KXVER=3 -Wall -fno-strict-aliasing -Wno-parentheses -g -O2 -shared -fPIC -o quartic.so -lm
clean:
	rm -f quartic.so
