import os
import re
import twobitreader

class Flanking():

  def __init__(self, hg_dir):
    self.assemblies = {
        "hg19": {},
        "hg38": {}
    }
    self.assemblies["hg19"] = twobitreader.TwoBitFile(os.path.join(hg_dir, "hg19.2bit"))
    self.assemblies["hg38"] = twobitreader.TwoBitFile(os.path.join(hg_dir, "hg38.2bit"))


  def get(self, assembly_version, chrom_shortname, chrom_start, chrom_end):
    if(assembly_version == "GRCh38"):
      assembly_name = "hg38"
      assembly = self.assemblies[assembly_name]
    else:
      assembly_name = "hg19"
      assembly = self.assemblies[assembly_name]

    if(chrom_shortname == "MT"):
      chrom_shortname = "M"
    chrom_name = "chr" + str(chrom_shortname)

    chrom_record = assembly[chrom_name]
    flanking_bp = chrom_record[chrom_start-1:chrom_end].upper()

    return flanking_bp

  def list_chromosomes(self):
    chr_names = [str(i) for i in list(range(1, 23))] + ["X", "Y", "M"]
    chr_names = ["chr%s" % i for i in chr_names]

    return chr_names

  def get_chromosome_length(self, chr_name):
    return len(self.assemblies["hg19"][chr_name])

  def get_chromosome_lengths(self):
    chr_lens = {}
    for chr_name in self.list_chromosomes():
      chr_lens[chr_name] = self.get_chromosome_length(chr_name)
    return chr_lens

  def get_context(self, chromosome, pos, assembly, ref, alt, mutation_type):    
    if mutation_type != 'single base substitution':
      #print('"' + row['mutation_type'] + '"')
      return None

    nl = 1
    nr = 1
    ref_complements = { 'G': 'C', 'A': 'T' }
    if ref in ref_complements.keys():
      ref_c = ref_complements[ref]
      alt_c = Flanking.reverse_complement(alt)
      left_flank = self.get(assembly, chromosome, pos-nr, pos-1)
      right_flank = self.get(assembly, chromosome, pos+1, pos+nl)
      left_flank_c = Flanking.reverse_complement(right_flank)
      right_flank_c = Flanking.reverse_complement(left_flank)
      colname = ("%s[%c>%c]%s" % (left_flank_c, ref_c, alt_c, right_flank_c))
    else:
      left_flank = self.get(assembly, chromosome, pos-nl, pos-1)
      right_flank = self.get(assembly, chromosome, pos+1, pos+nr)
      colname = ("%s[%c>%c]%s" % (left_flank, ref, alt, right_flank))

    if 'N' in colname:
      return None

    return colname

  @staticmethod
  def reverse_complement(seq):
    trans = str.maketrans('ATGC', 'TACG')
    return seq.translate(trans)[::-1]


