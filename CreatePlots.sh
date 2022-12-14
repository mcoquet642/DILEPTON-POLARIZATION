for TAU in `ls OUTPUT/DileptonRateTAU* | sed -e "s/OUTPUT\/DileptonRateTAU//" -e "s/.txt//" | sort -g`
do

echo "set terminal postscript eps enhanced color" >> next.plot
echo "set output 'DileptonRateTAU${TAU}.eps'" >> next.plot

echo "set xr [0:7]" >> next.plot
echo "set log y" >> next.plot
echo "set yr [1e-15:1e-5]" >> next.plot

echo "set xlabel 'Q'" >> next.plot
echo "set ylabel 'dN/d^4xdQ'" >> next.plot

echo "set st d lp" >> next.plot

INFOLINE=`head -n 1 OUTPUT/DileptonRateTAU${TAU}.txt`


echo "set label 'dN_{ch}/d{/Symbol h}=1900 | {/Symbol h}/s=0.16 | {/Symbol t}=${TAU} fm/c'  at screen 0.3,0.9" >> next.plot
echo "set label '${INFOLINE}' at screen 0.3,0.8" >> next.plot


STRING=""
STRING="$STRING 'OUTPUT/DileptonRateTAU${TAU}.txt' u 1:2 lc rgb 'red' dt ''  ti 'LO -- q+q -> l^{+}l^{-}',"
STRING="$STRING 'OUTPUT/DileptonRateTAU${TAU}.txt' u 1:4 lc rgb 'blue' dt '' ti 'NLO -- q+g -> l^{+}l^{-}',"
STRING="$STRING 'OUTPUT/DileptonRateTAU${TAU}.txt' u 1:3 lc rgb 'gray' dt '.' ti 'Equilibrium LO-- q+q -> l^{+}l^{-}',"

echo "p ${STRING}" >> next.plot

gnuplot next.plot
rm next.plot

ps2pdf -dEPSCrop DileptonRateTAU${TAU}.eps DileptonRateTAU${TAU}.pdf
rm DileptonRateTAU${TAU}.eps

done