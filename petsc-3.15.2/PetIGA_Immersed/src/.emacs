;; .emacs

;;; uncomment this line to disable loading of "default.el" at startup
;; (setq inhibit-default-init t)

;; enable visual feedback on selections
;(setq transient-mark-mode t)

;; default to better frame titles
(setq frame-title-format
      (concat  "%b - emacs@" (system-name)))

;; default to unified diffs
(setq diff-switches "-u")

;; always end a file with a newline
;(setq require-final-newline 'query)

;;; uncomment for CJK utf-8 support for non-Asian users
;; (require 'un-define)


(global-set-key "\eOp" "0")
(global-set-key "\eOq" "1")
(global-set-key "\eOr" "2")
(global-set-key "\eOs" "3")
(global-set-key "\eOt" "4")
(global-set-key "\eOu" "5")
(global-set-key "\eOv" "6")
(global-set-key "\eOw" "7")
(global-set-key "\eOx" "8")
(global-set-key "\eOy" "9")
(global-set-key "\eOl" "+")
(global-set-key "\eOQ" "/")
(global-set-key "\eOR" "*")
(global-set-key "\eOS" "-")
(global-set-key "\eOn" ".")

;Comment macro
(fset 'comment
   "\C-ac\C-n")
(global-set-key (kbd "C-c f") 'comment)

;Uncomment macro
(fset 'uncomment
   "\C-a\C-d\C-n")
(global-set-key (kbd "C-c d") 'uncomment)
