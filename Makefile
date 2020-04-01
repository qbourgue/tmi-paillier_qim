Programme : main.o
	gcc -g -Wall main.o -lgmp -o Programme

tp1m.o : main.c
	gcc -Wall -c main.c -lgmp -o main.o
