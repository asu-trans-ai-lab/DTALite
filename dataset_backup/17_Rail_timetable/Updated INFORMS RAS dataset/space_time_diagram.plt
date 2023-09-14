set title "Space time trajectory diagram"
set xlabel "Time Horizon"
set ylabel "Space (distance)"  offset -1
set xtics (" 0:00" 0 ," 2:00" 120 ," 4:00" 240 ," 6:00" 360 ," 8:00" 480 ,"10:00" 600 ,"12:00" 720 ,"14:00" 840 ,"16:00" 960 ,"18:00" 1080 ,"20:00" 1200 ,"22:00" 1320 ,"24:00" 1440 ) 
set ytics (" " 0)
set xrange [0:1441] 
set yrange [0:80.00] 
plot "agent1.txt" using 1:2 title 'agent 1'  with lines,\
"agent2.txt" using 1:2 title 'agent 2'  with lines
