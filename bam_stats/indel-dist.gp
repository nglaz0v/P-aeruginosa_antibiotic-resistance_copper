
        set terminal png size 600,400 truecolor
        set output "bam_stats/indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "PAeruginosa.sorted" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	29860
2	3493
3	4053
4	522
5	481
6	978
7	166
8	163
9	459
10	159
11	70
12	133
13	1
14	49
15	52
16	58
17	16
18	12
19	0
20	0
21	0
22	0
24	0
25	0
end
1	33634
2	3749
3	3152
4	1079
5	335
6	1238
7	191
8	277
9	302
10	102
11	371
12	450
13	1
14	7
15	88
16	31
17	15
18	72
19	31
20	2
21	5
22	1
24	2
25	1
end
1	0.887792
2	0.931715
3	1.285850
4	0.483781
5	1.435821
6	0.789984
7	0.869110
8	0.588448
9	1.519868
10	1.558824
11	0.188679
12	0.295556
13	1.000000
14	7.000000
15	0.590909
16	1.870968
17	1.066667
18	0.166667
19	0.000000
20	0.000000
21	0.000000
22	0.000000
24	0.000000
25	0.000000
end
