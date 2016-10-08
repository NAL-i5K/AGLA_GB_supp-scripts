# fix_pseudogene.py
# =================
#requires https://github.com/hotdogee/gff3-py
from gff3 import Gff3
gff = Gff3('agla_v1_2-NALmod7.gff3')
type_map = {'exon': 'pseudogenic_exon', 'transcript': 'pseudogenic_transcript'}
pseudogenes = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'pseudogene']
for pseudogene in pseudogenes:
    # convert types
    for line in gff.descendants(pseudogene):
        if line['type'] in type_map:
            line['type'] = type_map[line['type']]
            print "changed line", line
gff.write('agla_v1_2-NALmod8.gff3')
