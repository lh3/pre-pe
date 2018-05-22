CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
PROG=pre-adna pre-lianti pre-dip-c pre-meta

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

pre-adna:kthread.o pre-adna.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

pre-lianti:kthread.o pre-lianti.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

pre-dip-c:kthread.o pre-dip-c.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

pre-meta:kthread.o pre-meta.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

clean:
		rm -fr gmon.out *.o ext/*.o a.out *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

pre-adna.o: kvec.h kstring.h khash.h kseq.h
pre-dip-c.o: kvec.h kseq.h
pre-lianti.o: kvec.h khash.h kseq.h
pre-meta.o: kvec.h kseq.h
