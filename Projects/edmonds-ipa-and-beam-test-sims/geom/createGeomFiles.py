# A script to generate the text files

radius_change_line="350"
new_radii=[ "350", "300", "250", "200" ]
length_change_line="500"
new_lengths=[ "0", "100", "200", "300", "400", "500"  ]

for radius in new_radii:
    for length in new_lengths:

        if radius=="350" and length=="500":
            continue;

        geom_file = open('protonAbsorber_cylindrical_r350mm_halfl500mm.txt', 'r');
        geom_filename="protonAbsorber_cylindrical_r"+radius+"mm_halfl"+length+"mm.txt"
        out_geom_file = open(geom_filename, 'w');

        for line in geom_file:
            if radius_change_line in line:
                newline = line.replace(radius_change_line, radius)
                out_geom_file.write(newline)
            elif length_change_line in line:
                newline = line.replace(length_change_line, length)
                out_geom_file.write(newline)
            else :
                out_geom_file.write(line)

        top_geom_file = open('geom_common_ipa_r350mm_halfl500mm.txt', 'r');
        top_geom_filename="geom_common_ipa_r"+radius+"mm_halfl"+length+"mm.txt"
        out_top_geom_file = open(top_geom_filename, 'w');

        for line in top_geom_file:
            if radius_change_line in line: # both radius and length should be on the same line
                newline = line.replace(radius_change_line, radius)
                newline = newline.replace(length_change_line, length)
                out_top_geom_file.write(newline)
            else :
                out_top_geom_file.write(line)

