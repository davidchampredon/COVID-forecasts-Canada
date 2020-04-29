
echo "LaTeX compilation for $1"

# Must compile twice for references
# Arguments for LaTeX: https://tex.stackexchange.com/questions/1492/passing-parameters-to-a-document
pdflatex "\def\prov{$1}\input{template-prov.tex}"
pdflatex "\def\prov{$1}\input{template-prov.tex}"

thedate=$(date +%Y_%m_%d)
# thedate='2020_04_18'  # <--- Manual overide

mv template-prov.pdf fcst_$1_$thedate.pdf

# Clean
rm template-prov*.aux
rm template-prov*.log
rm template-prov*.out
rm template-prov*.gz

