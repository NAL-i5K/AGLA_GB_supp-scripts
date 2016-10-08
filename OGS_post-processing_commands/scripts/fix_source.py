# fix_source.py
# =================

# script uses Han Lin's gff3 class to import a gff3, take a list of IDs that currently have source ManualCuration, and switch the source of these models (all parents and children) to I5K
# Requires: https://github.com/hotdogee/gff3-py
from gff3 import Gff3
gff = Gff3('agla_v1_2-NALmod3.gff3')
id_list = ['AGLA014663', 'AGLA003801', 'AGLA017751', 'AGLA003809','AGLA000919','AGLA000103']
source_map = {'ManualCuration': 'I5K'}

for feature_id in id_list:
    for feature in gff.features[feature_id]:
        for line in gff.descendants(feature):
            if line['source'] in source_map:
                line['source'] = source_map[line['source']]
            if feature['source'] in source_map:
                feature['source'] = source_map[feature['source']]

gff.write('agla_v1_2-NALmod4.gff3')


