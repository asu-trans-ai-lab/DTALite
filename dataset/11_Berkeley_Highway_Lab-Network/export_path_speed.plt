set title "Dynamic Speed Contour (Path 1 ) Unit: mph" 
set xlabel "Time Horizon"
set ylabel "Space (Node Sequence)"  offset -1
set xtics (" 0:00" 0 ," 0:10" 10 ," 0:20" 20 ," 0:30" 30 ," 0:40" 40 ," 0:50" 50 ," 1:00" 60 ," 1:10" 70 ," 1:20" 80 ," 1:30" 90 ," 1:40" 100 ," 1:50" 110 ," 2:00" 120 ) 
set ytics ("12" 0, "13" 2, "14" 4, "15" 7, "16" 10)
set xrange [0:121] 
set yrange [0:10] 
set palette defined (0 "white", 0.1 "red", 40 "yellow", 50 "green")
set pm3d map
splot 'C:\GitHub\dtalite_software_release\release\version_GMNS0.9\11_Berkeley_Highway_Lab-Network\export_path_speed.txt' matrix notitle
