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
(require `ipython)

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

(defun py-execute-statement-and-step ()
 "select a statement, submit as a region and then step forward"
 (interactive)
 (beginning-of-line 1)
 (let ((beg (point)))
   (py-next-statement 1)
   ; if last statement.
   (if (= (point) beg) (end-of-buffer ))
   (py-execute-region beg (point))
 )
)

;  some key bindings
(define-key py-mode-map (quote [f9]) 'py-execute-statement-and-step)
(define-key py-mode-map (quote [f10]) `py-execute-region)
(define-key py-mode-map (quote [f11]) `py-execute-buffer)
(define-key py-mode-map (quote [f12]) `goto-python-shell)
(define-key py-shell-map (quote [f12]) `goto-python-source)

