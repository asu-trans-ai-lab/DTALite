set title "Space time trajectory diagram"
set xlabel "Time Horizon"
set ylabel "Space (distance)"  offset -1
set xtics (" 0:00" 0 ," 2:00" 120 ," 4:00" 240 ," 6:00" 360 ," 8:00" 480 ,"10:00" 600 ) 
set ytics (" " 0)
set xrange [0:661] 
set yrange [0:4.00] 
plot "agent1.txt" using 1:2 title 'agent 1'  with lines,\
"agent2.txt" using 1:2 title 'agent 2'  with lines,\
"agent3.txt" using 1:2 title 'agent 3'  with lines,\
"agent4.txt" using 1:2 title 'agent 4'  with lines
