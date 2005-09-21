; This file provide some convenient ways to debug your program
; More specifically, in python mode of emacs, you can use
;     F9    to run a statement and forward to the next
;     F10  to run selected text
;     F11  to run the whole buffer.
;     F12  to switch between src and shell buffer.
;          start python-shell if necessary
;
; This file can be freely distriuted under GPL
;
; Author:  Bo Peng (ben.bob at gmail.com)
; Date:    Sep, 2004
; Last Modified: Oct, 14, 2004

(require `python-mode)

; ipython (with Prabhu updated python-mode.el )
; If you are using the usual python-mode, comment out the following
(setq py-shell-switch-buffers-on-execute nil)
; (require `ipython)

(defun goto-python-shell ()
 "Go to the python command window (start it if needed)"
 (interactive)
 (setq current-python-script-buffer (current-buffer))
 (py-shell)
 (end-of-buffer)
)

(defun goto-python-source ()
 "switch back to source window"
 (interactive)
 (switch-to-buffer-other-window current-python-script-buffer)
)

(defun goto-endof-quote ()
"goto end of quote"
  (interactive)
  (beginning-of-line)
  (while (not (looking-at "\"\"\""))
      (forward-line 1)
  )
  (forward-line 1)
)

(defun goto-endof-if-for ()
"goto end of if for"
  (interactive)
  (forward-line 1)
  (beginning-of-line)
    (while (looking-at " ")
      (forward-line 1)
      (beginning-of-line)
    )
   (forward-line 1)
)


;;; have to fix """ .... """ problem.
;;; 
;;; for my own use, hack this here.
;;; only for the case
;;;
;;;  some = """
;;;   ...
;;;  """
(defun py-goto-end-of-statement ()
  "Go to end of this statement."
 (let ((start (point)))
    (while (and
            (or (looking-at py-blank-or-comment-re)
                (py-in-literal))
            (not (eobp)))
       (forward-line 1))
  (end-of-line)
  (if (string= (concat (char-to-string (char-before) )
                       (char-to-string (char-before (- (point) 1))))
         "\"\"")
     (goto-endof-quote) 
  (progn
    (beginning-of-line)
    (if (or (or (or (looking-at "if" )   ; will look for the first ono-blank line
          (looking-at "for") )
          (looking-at "while") )
          (looking-at "def") )
       (goto-endof-if-for)
   (progn 
      (py-goto-beyond-final-line)
      (if (eobp)
        (progn (goto-char (point-max)) nil)
      t)
   ))))
  ) ; let
)


(defun py-execute-statement-and-step ()
 "select a statement, submit as a region and then step forward"
 (interactive)
 (beginning-of-line 1)
 (let ((beg (point)))
   ; (py-next-statement 1)
   (py-goto-end-of-statement)
   ; if last statement.
   (if (= (point) beg) (end-of-buffer))
   (py-execute-region beg (point))
 )
)

(defun py-execute-buffer-step-by-step ()
"execute the whole buffer but interactively"
(interactive)
 (save-excursion
     (goto-char (point-min) )
     (while (not (equal (point) (point-max)))
       (py-execute-statement-and-step)
     )
 )
)

(defun save-python-session ()
"execute the whole buffer but interactively"
(interactive)
 (save-excursion
    (py-shell)
   (write-region (point-min) (point-max) "userGuide.log")
 )
)


;  some key bindings
(define-key py-mode-map (quote [f9]) 'py-execute-statement-and-step)
(define-key py-mode-map (quote [f10]) `py-execute-region)
(define-key py-mode-map (quote [f11]) `py-execute-buffer-step-by-step)
(define-key py-mode-map (quote [f12]) `goto-python-shell)
(define-key py-shell-map (quote [f12]) `goto-python-source)
(define-key py-mode-map (quote [f8]) `save-python-session)

