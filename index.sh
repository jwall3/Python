gmx make_ndx -n index2.ndx -o index2.ndx <<EOF
!0
q
EOF

gmx trjconv -f inxtc -s tpr -o outxtc -n index2.ndx -dt 1000 <<EOF
24
EOF
                         
