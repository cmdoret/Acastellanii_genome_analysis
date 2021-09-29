"""
Uses biopython to plot a karyotype from input fasta and bed files.
cmdoret, 2020
"""
import sys
import re
from collections import defaultdict
import csv
from Bio import SeqIO
from Bio.Graphics import BasicChromosome
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib.units import cm
from matplotlib import cm as mcm
from matplotlib.colors import to_hex

logh = open(snakemake.log[0], "w")
# {chrom: length}
entries = {}
for rec in SeqIO.parse(snakemake.params["genome"], "fasta"):
    if len(rec.seq) > 100000:
        entries[rec.id] = len(rec.seq)
# all available colors
# cols = list(colors.getAllNamedColors().values())[8:]
cols = [to_hex(mcm.tab10(c)) for c in range(len(entries))]
# {chrom: {type: [features]}} where a feature is [start, end, type, color]
features = defaultdict(list)
type_id = 0
ftype_to_col = {}
with open(snakemake.input["features"]) as bedh:
    bed = csv.DictReader(
        bedh, delimiter="\t", fieldnames=["chrom", "start", "end", "type"]
    )
    for feat in bed:
        start, end = int(feat["start"]), int(feat["end"])
        if start > end:
            start, end = end, start
        try:
            if end > entries[feat["chrom"]]:
                continue
        except KeyError:
            continue
        if feat["type"] not in ftype_to_col.keys():
            ftype_to_col[feat["type"]] = cols[type_id]
            logh.write(f"{feat['type']}: {cols[type_id]}\n")
            type_id += 1
        features[feat["chrom"]].append(
            SeqFeature(
                location=FeatureLocation(int(feat["start"]), int(feat["end"])),
                type=feat["type"],
                qualifiers={"color": [ftype_to_col[feat["type"]]]},
            )
        )

max_len = max(entries.values())
telomere_length = 10000  # For illustration

chr_diagram = BasicChromosome.Organism(output_format="svg")
chr_diagram.page_size = (29.7 * cm, 21 * cm)  # A4 landscape

for name, length in entries.items():
    # features = [f for f in record.features if f.type == "tRNA"]
    chrom_feat = features[name]

    chrom_num = re.sub(r"^.*_([0-9]+)$", r"\1", name)
    cur_chromosome = BasicChromosome.Chromosome(chrom_num)
    # Set the scale to the MAXIMUM length plus the two telomeres in bp,
    # want the same scale used on all five chromosomes so they can be
    # compared to each other
    cur_chromosome.scale_num = max_len + 2 * telomere_length

    # Add an opening telomere
    start = BasicChromosome.TelomereSegment()
    start.scale = telomere_length
    cur_chromosome.add(start)

    # Add a body - again using bp as the scale length here.
    body = BasicChromosome.AnnotatedChromosomeSegment(length, chrom_feat)
    body.scale = length
    cur_chromosome.add(body)

    # Add a closing telomere
    end = BasicChromosome.TelomereSegment(inverted=True)
    end.scale = telomere_length
    cur_chromosome.add(end)

    # This chromosome is done
    chr_diagram.add(cur_chromosome)
chr_diagram.draw(snakemake.output[0], snakemake.params["title"])
logh.close()
