#!/usr/bin/ol --run

;;;
;;; Spektro -- a simple tool to assist in analysis of mass spectrometry runs
;;;

;; todo
; - add missing pieces to intervals
; - add config file for listing elements, masses and stuff

(import (owl args)) ;; command line argument parsing

(define (abort) (halt 1)) ;; nonzero program exit on errors


;;;
;;; Single goal solver 
;;; 

; find matching factors to satisfy the mass and charge with given elems
;  success: (report factors n) -> n-1
;  failure: n

(define (find-factors mass charge elems moves report n)
	(cond
		((eq? n 0) n)
		((eq? mass 0) 
			(if (eq? charge 0) (report moves n) n))
		((null? elems) n)
		((< mass 0) n)
		(else
			(bind (car elems)
				(lambda (name this-mass this-charge min)
					(let try-factors ((count min) (n n))
						(let ((next-mass (- mass (* count this-mass))))
							(if (< next-mass 0)
								n
								(try-factors (+ count 1)
									(find-factors next-mass 
										(- charge (* count this-charge)) 
										(cdr elems)
										(cons (cons count name) moves)
										report n))))))))))

; a solution is interesting if it is not a multiple of a smaller solution

(define (gcdl lst)
	(if (null? lst) 
		0
		(fold gcd (car lst) (cdr lst))))

(define (interesting? factors div)
	(= 1 (gcdl (cons div (map (o abs car) factors))))) 

(define (solve-goal mass charge elems n cook)
	(let loop ((divisor 1) (n n))
		(if (> n 0)
			(loop (+ divisor 1)
				(find-factors 
					(* mass divisor)
					(* charge divisor) 
					elems
					null
					(lambda (factors n)
						(cook mass charge factors divisor n))
					n))
			(begin
				(print "Enough goals computed for m/z " mass ". Use -n <limit> to get more if needed.")
				0))))


;;;
;;; Interval finder
;;;

(define (try-interval this step masses picked)
	(let 
		((picked (cons this picked))
		 (next (+ this step)))
		(cond
			((get masses next #false)
				(try-interval next step masses picked))
			((>= (length picked) 3)
				(reverse picked))
			(else #false))))

(define (walk-intervals lst step masses max)
	(cond
		((= step max)
			null)
		((try-interval (car lst) step masses null) =>
			(lambda (matched)
				(cons
					(cons step matched)
					(walk-intervals lst (+ step 1) masses max))))
		(else
			(walk-intervals lst (+ step 1) masses max))))

(define (find-intervals masses)
	(if (< (length masses) 3)
		(begin 
			(print "Too few masses for meaningful intervals: " masses)
			1)
		(let*
			((masses (sort < masses))
			 (max (car (reverse masses)))
			 (massff (fold (lambda (all this) (put all this this)) #empty (sort < masses))))
			(let loop ((lst masses))
				(if (> (length lst) 2)
					(append
						(walk-intervals lst 1 massff (- max (car lst)))
						(loop (cdr lst)))
					null)))))

(define (all-intervals masses)
	(let* 
		((all (find-intervals masses))
		 (lend (map (lambda (x) (cons (length x) x)) all))
		 (all (sort (lambda (a b) (< (car a) (car b))) lend))
		 (all (map cdr all)))
		(if (null? all)
			(print "No intervals found")
			(for-each 
				(lambda (lst)
					(print*
						(list "Interval with delta " (car lst) " and length " (length (cdr lst)) ": " (cdr lst))))
				all))))




;;;
;;; Formula formatting 
;;;

; signer = abs, lead or 

(define (format-formula factors signer)
	(cond
		((null? factors)
			null)
		((string? (car factors))
			(cons (car factors)
				(format-formula (cdr factors) 'lead)))
		((eq? (caar factors) 0)
			(format-formula (cdr factors) signer))
		((eq? signer 'abs)	; show the absolute value
			(let* 
				((this (car factors))
				 (n (abs (car this))))
				(cond
					((= n 1)
						(cons (cdr this)
							(format-formula (cdr factors) 'normal)))
					(else
						(ilist n
							(cdr this)
							(format-formula (cdr factors) 'normal))))))
		((eq? signer 'lead)	; show the sign in the number
			(let* ((this (car factors)) (n (car this)))
				(cond
					((= n 1)
						(cons (cdr this)
							(format-formula (cdr factors) 'normal)))
					((= n -1)
						(ilist "-" (cdr this)
							(format-formula (cdr factors) 'normal)))
					(else
						(ilist n
							(cdr this)
							(format-formula (cdr factors) 'normal))))))
		; emit a separate sign and continue as abs
		(else
			(let ((next (car factors)))
				(if (and (pair? next) (< (car next) 0))
					(cons " - "
						(format-formula factors 'abs))
					(cons " + "
						(format-formula factors 'abs)))))))

(define (format-solution mass charge factors divisor n)
	(print*
		(append
			(list "m/z " mass " = ")
			(append
				(format-formula factors 'lead)
				(if (> divisor 1)
					(list " / " divisor)
					null))))
	(- n 1))

(define (maybe-format-solution mass charge factors divisor n)
	(if (interesting? factors divisor)
		(format-solution mass charge factors divisor n)
		n))



;;;
;;; Series solver
;;;

(define (try-intervals lst-in)
	(let loop ((step 1) (this (car lst-in)) (lst (cdr lst-in)))
		(cond
			((null? lst)
				step)
			((eq? (+ this step) (car lst))
				(loop step (car lst) (cdr lst)))
			((> (+ this step) (car lst))
				(print "No interval applicable, stopping at " step)
				#false)
			(else
				(loop (+ step 1) (car lst-in) (cdr lst-in))))))

(define (get-interval lst)
	(if (< (length lst) 2)
		(begin
			(print "Series is too short")
			#false)
		(try-intervals lst)))

; factors = ((n . name) ...)
;; fixme, what is this then?

(define (adjust-elements elems factors steps)
	(fold 
		(lambda (elems factor)
			(map
				(lambda (elem)
					(bind elem
						(lambda (name mass charge min)
							(if (eq? name (cdr factor))
								(tuple name mass charge (div (- 0 (car factor)) steps))	
								; <- negative factors to allow removing originals in the iterations
								elem))))
				elems))
		elems factors))
		

(define (render-element elem)
   (if (and (tuple? elem) (= 4 (size elem))) ;; sanity check
      (lets ((name mass charge _ elem))
         (list->string
            (foldr render null
               (list "   name: " name ", mass: " mass ", charge: " charge))))
      elem))

; (cook mass charge factors divisor n) -> n'

(define (solve-series masses charge elems nequs)
	(let* 
		((masses (sort < masses))
		 (interval (get-interval masses))
		 (series-length (- (length masses) 1)))
		(if interval
			(begin
				(print "Solving mass series " masses " having interval " interval ".")
            (print "Total charge must be " charge ".")
            (print "Elements are:")
            (for-each (o print render-element) elems)
            (print "")
				(let 
					((n-left
						(solve-goal (car masses) charge elems nequs
							(lambda (mass charge factors divisor n)
								; may not be interesting but the extension can be
								(find-factors 
									(* interval divisor)
									0
									(adjust-elements elems factors series-length)
									null
									(lambda (chain-factors n)
										; check the full molecule for interestingness
										(if (interesting? (append chain-factors factors) divisor)
											(format-solution (car masses) charge 
												(append factors 
													(cons " + [" (append chain-factors (list "]"))))
												divisor
												n)
											n))
									n)))))
					(begin
						(if (= n-left 0)
							'ok
							(print "Unable to compute enough solutions for m/z " (car masses)))
						0)))
			(print "No repeating increment in " masses))))


;;;
;;; Command line processing 
;;;

(define (zap-elem fns args)
	(cond
		((null? args) null)
		((null? fns)
			(print "Too many arguments in element.")
			(abort))
		(else
			(cons ((car fns) (car args))
				(zap-elem (cdr fns) (cdr args))))))

(define (parse-element str)
	(let*
		((lst (c/,/ str))
		 (lst 
			(zap-elem
				(list
					(lambda (name-str) name-str)
					(lambda (mass-str) 
						(let ((mass (string->number mass-str 10)))
							(cond
								((not mass)	
									(print "Bad mass definition: " mass-str)
									(abort))
								((< mass 1)
									(print "Rather strange mass: " mass)
									(abort))
								(else mass))))
					(lambda (charge-str)
						(let ((charge (string->number charge-str 10)))
							(if charge charge
								(begin
									(print "Bad charge for element: " charge-str)
									(abort))))))
				lst))
		 (len (length lst)))
		(cond
			((eq? len 1)
				(print "No mass given for " (car lst))
				(abort))
			((eq? len 2)
				(print "No charge given for " (car lst))
				(abort))
			((eq? len 3)
				(tuple (car lst) (cadr lst) (caddr lst) 0))			; min element lisatty
			(else
				(print "Bad element definition: " str)
				(abort)))))

(define command-line-rules
   (cl-rules
      `((charge "-c" "--charge" cook ,(λ (arg) (string->number arg 10)) comment "total charge")
        (elem "-e" "--element" cook ,parse-element
            comment "define alement in form <name>,<mass>,[+-]<charge>"
            plural)
        (series "-s" "--series"
            comment "solve a mass-series by growing a submolecule")
        (count "-n" "--count" cook ,(λ (x) (string->number x 10)) default "10"
            comment "maximum number of unique solutions to find per given mass")
        (intervals "-i" "--intervals"
            comment "find interesting intervals in masses")
        ;; standard flags
        (about "-A" "--about")
        (version "-V" "--version")
        (help "-h" "--help"))))


(define usage-text "Usage: spektro [args] [mass] ...")

(define version "Spektro v0.1")

;;;
;;; Startup
;;;

(define (spektro args)
	(or
		(process-arguments (cdr args) command-line-rules usage-text
			(lambda (dict others)
				(let* 
					((elems (get dict 'elem null))
					 (elems (if (tuple? elems) (list elems) elems))
					 (charge (get dict 'charge 1)))
					(cond
						((getf dict 'help)
							(print usage-text)
							(print-rules command-line-rules)
							0)
						((getf dict 'about)
                     (print version " -- a simple tool to assist in analysis of mass-spectrometry results.")
                     (print "Copyright (c) Aki Helin, written at OUSPG")
                     (print "Some documentation is available at http://code.google.com/p/ouspg/wiki/Spektro.")
							0)
                  ((getf dict 'version)
                     (print version)
                     0)
						((null? others)
							(print "No masses?")
							0)
						((getf dict 'intervals)
							(let ((masses (map (λ (x) (string->number x 10)) others)))
								(if (all number? masses)
									(all-intervals masses)
									(begin
										(print "Cannot check intervals: bad masses: " others)
										1))))
						((null? elems)
							(print "No elements given. Hard to solve goals without them.")
							0)
						((getf dict 'series)
							(let ((masses (map (λ (x) (string->number x 10)) others)))
								(if (all number? masses)
									(solve-series masses charge elems (get dict 'count 10))
									(print "Bad series masses: " masses))))
						(else
							(for-each
								(lambda (goal)
									(let ((mass (string->number goal 10)))
										(cond
											((not mass)
												(print "Bad goal: " goal))
											((< mass 1)
												(print "Goal too small: " goal))
											(else 
												(fork-named (tuple 'mass mass)
													(lambda ()
														(solve-goal mass charge elems (get dict 'count 10) maybe-format-solution)))))))
								others))))))
		1))

; (spektro (list "-c" "0" "-s" "-e" "A,1,0" "5" "7" "9"))

; (dump (lambda (args) (set-signal-action 'halt) (spektro (cdr args))) "spektro.c")

spektro


