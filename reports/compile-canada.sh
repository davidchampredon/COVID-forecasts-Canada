
echo "LaTeX compilation for Canada"

# Must compile twice for references
pdflatex template-canada.tex
pdflatex template-canada.tex

thedate=$(date +%Y_%m_%d)
# thedate='2020_04_04'  # <--- Manual overide

mv template-canada.pdf fcst_Canada_$thedate.pdf

# Clean
rm template-canada*.aux
rm template-canada*.log
rm template-canada*.out
rm template-canada*.gz

