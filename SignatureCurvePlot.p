# Gnuplot script file for plotting data in file "../output/test_geodict_import_gdt_import_test_13Dec05_17_11_31/Uoft.log"
      # This file is called   

      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "Temperature dependency of diffusion coefficient in electrolyte "
      set xlabel "SOC [-]"
      set ylabel "Cell voltage [V]"
#      set key 0.01,100
      set grid xtics lt 0 lw 1 
      set grid ytics lt 0 lw 1 
#      set label "Yield Point" at 0.003,260
#      set arrow from 0.0028,250 to 0.003,280
      set xr [0:0.95]
      set yr [0:1.5]
      set key left top
     plot "~/MA/BEST/Anode_A/A_Test_DiffCoeff/Output/A_Test_DiffCoeff_1C_20141014_10_14_41/Uoft.log" using ($1/3600):($3) title '298K' with lines, \
	  "~/MA/BEST/Anode_A/A_Test_DiffCoeff/Output/A_Test_DiffCoeff_1C_330K_20141014_16_08_38/Uoft.log" using ($1/3600):($3) title '330K' with lines, \
	  "~/MA/BEST/Anode_A/A_Test_DiffCoeff/Output/A_Test_DiffCoeff_1C_330K_DiffCoeffT_20141015_18_04_47/Uoft.log" using ($1/3600):($3) title '330K with dependence' with lines
     set size 1.0, 0.6
     set terminal postscript portrait enhanced  dashed lw 1 "Helvetica" 14 color
     set output "A_Test_DiffCoeff_1C.ps"
     replot
     set terminal x11
     set size 1,1
