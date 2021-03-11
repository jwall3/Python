gmx make_ndx -n index2.ndx -o index2.ndx <<EOF
!0
q
EOF

gmx editconf -f nvt.gro -o nvt_out.gro -n index2.ndx <<EOF
24
EOF

gmx trjconv -f md_xx.xtc -s md.tpr -o mdout.xtc -n index2.ndx -dt 100 <<EOF
24
EOF
