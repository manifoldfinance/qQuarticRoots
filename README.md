# qQuarticRoots
C library for KDB+ quartic solver

## compile quartic.so
```bash
$ make
```
## load quartic.q to run quartic function
```q
q) \l quartic.q
q) quarticRoots[ 3.0; 6.0; -123.0; -126.0; 1080.0]
5 3 -4 -6f
```
