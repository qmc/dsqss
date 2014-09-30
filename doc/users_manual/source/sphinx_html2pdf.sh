#!/bin/bash

page1=index
page2=intro
page3=algorithm
page4=install
page5=exe_method
page6=sample
page7=qa
page8=reference
page7=appendix
htmldir=../build/html
outfile=dsqss_manual

wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page1.html $page1.pdf
wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page2.html $page2.pdf
wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page3.html $page3.pdf
wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page4.html $page4.pdf
wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page5.html $page5.pdf
wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page6.html $page6.pdf
wkhtmltopdf --user-style-sheet pdf.css  $htmldir/$page7.html $page7.pdf

pdftk $page1.pdf $page2.pdf $page3.pdf $page4.pdf $page5.pdf $page6.pdf $page7.pdf cat output $outfile.pdf

rm *.html &>/dev/null

