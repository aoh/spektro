OFLAGS=-O1
CFLAGS=-O2

spektro: spektro.c
	$(CC) $(CFLAGS) -o spektro spektro.c

spektro.c: spektro.scm 
	# fetch owl only if spektro.scm is updated and no owl yet
	test -x owl-lisp/bin/ol || make owl-lisp/bin/ol
	owl-lisp/bin/ol $(OFLAGS) -o spektro.c spektro.scm

install: spektro
	mkdir -p $(PREFIX)/bin
	install -m 755 spektro $(PREFIX)/bin/spektro

uninstall:
	-rm $(PREFIX)/bin/spektro

owl-lisp/bin/ol:
	# symlink owl-lisp here if you already have it to 
	# avoid a rebuild
	-git clone https://github.com/aoh/owl-lisp.git
	-cd owl-lisp && git pull
	cd owl-lisp && make

clean:
	-rm spektro.c spektro
