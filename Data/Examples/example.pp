#PP 2.0
  This is an example of the new pp format 2.0
  We support comments;
  all lines starting with white space are considered comments
  The new format must start with a header #PP 2.0 

  Alignments look very much like clustal format. 
  Anchors can be specified by lines #A<n>; if anchors are given,
  there specification has to be complete (contiguous numbers <n> in 1..maxn, etc.).
  The structure annotation below is just a comment (it is ignored!)

seqA ACGTUG
seqB ACG-AC
     (((...
#A1  ...N..
#A2  ...1..

seqA ACGCGUAA
seqB ACGCGUAC
     ...)))..
#A1  .N......
#A2  .2......

#END

  Non-comment or non-blank lines between sections are disallowed; this could be confusing otherwise

#SECTION BASEPAIRS

#BPCUT 0.2
#STACKS

  In this section, all lines <i> <j> <p> [<p2>] specify base pair (and stack) probabilities;
  
1 6 0.12 0.99
1 7 0.50 0.98
1 12 0.8 0.9
2 11 0.7 0.69
3 10 0.8

#END

#SECTION INLOOP

#BPILCUT 0.001
#UILCUT 0.01

  In this section, all lines "<i> <j> : {<ix> <jx> <px>} ; { <kx> <px> }" specify in loop
  probabilities.
  
  We support multiline specifications, see below.

1 12 : 1 3 0.4 1 5 0.9 2 10 0.04 2 11 0.9 ; 2 0.7 3 0.8
2 11 : \
  3 10 0.4 \
  3 8 0.02 ; \
  3 0.07 4 0.06 10 0.01
3 10 : ; 4 0.1

  The specification of probabilities in the external loop looks like
  probabilities for loop of the pseudo base pair 0 <length>+1:
0 15 : 1 12 0.8 2 11 0.1 ; 13 0.99 14 0.99


#END
