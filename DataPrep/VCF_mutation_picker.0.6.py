import csv, sys
from collections import defaultdict

pheno_csv, vcf_file, mgo_file  = sys.argv[1:]

# CSV file sample:
#SampleID,Histology,Patient
#1A,N,1
#1C,M,1

pheno = csv.reader(open(pheno_csv, 'rU'))
pheno.next() # skip header.
dog2samples, sample2dog_info = defaultdict(lambda: defaultdict(list)), {}
for row in pheno:
        if len(row) != 3:
                print len(row), row
                sys.exit(0)
        sample, hist, patient = row
        assert sample not in sample2dog_info
        sample2dog_info[sample] = (patient, hist)
        dog2samples[patient][hist in ['N']].append(sample) # dog -> (true -> control samples, false -> case samples). TODO: Generalize/paramterize "normal" flag. Note that "'N'" appears below too.
        
### Variables ###
DP_cut=4 ### Must be greater than this value
GQ_cut=40 ### Must be greater than this value
IMPACTS=["HIGH", "MODERATE"]  ### SNPEFF mutation impacts to keep they have 4 possible values from least to most MODIFIER,LOW,MODERATE,HIGH

#KG
def score(control,case):
	diff=[]
	ctl=[]
	for i in control:
		if i[1]>=DP_cut and i[2]>=GQ_cut:
			if i in ctl:
				pass
			else:
				ctl.append(i[0])
	for cas in case:
		if len(ctl)==0:
			diff.append(".")
		elif cas[1]<DP_cut or cas[2]<GQ_cut:
			diff.append(".")
		elif cas[0] in ctl:
			diff.append("0")
		else:
			diff.append("1")
	return diff
					
col2sample = None
gene2summary = defaultdict(lambda: defaultdict(int))
for line in open(vcf_file):
        if line.startswith('#CHROM'):
                col2sample = dict(enumerate(line[:-1].split('\t')[9:]))
                assert len(set(col2sample.values())) == len(sample2dog_info)
                cols = sorted(col2sample.keys())
                for dog in sorted(dog2samples.keys()):
                        print '#%s:\t%s %s'%(dog,
                                              ','.join(['%s(%d)'%(col2sample[c], c) for c in cols if col2sample[c] in dog2samples[dog][True]]),
                                              ','.join(['%s(%d)'%(col2sample[c], c) for c in cols if col2sample[c] in dog2samples[dog][False]]))
                continue

        if line[0] == '#': continue

        ff = line[:-1].split('\t')
        ann = [p for p in ff[7].split(';') if p.startswith('ANN=')]
        if ann:
                assert len(ann) == 1
                annff = ann[0].split('|')
                if annff[2] in IMPACTS or (annff[2] == 'LOW' and annff[1].lower() == 'missense_variant'):
                        gene = annff[3]
                        dogs, dog2control, dog2case = set(), defaultdict(list), defaultdict(list)
                        for x, d in enumerate(ff[9:]):
                                val = ('./.', 0, 0)
                                if d != './.':
                                        dff = d.split(':')
                                        val = (dff[0], int(dff[2]), int(dff[3]))
                                
                                sample = col2sample[x]
                                dog, hist = sample2dog_info[sample]
                                dogs.add(dog)
                                if hist == 'N':
                                        dog2control[dog].append(val)
                                else:
                                        dog2case[dog].append((sample, val))
                        for dog in sorted(dogs):
                                scores = score(dog2control[dog], [v for sample, v in dog2case[dog]])
                                for sample, s in zip([sample for sample, v in dog2case[dog]], scores):
                                        gene2summary[gene][sample] += s == '1'
                                print '\t'.join([dog, gene, ff[0]+':'+ff[1]] + scores)

msamples = [s for s in sorted(sample2dog_info.keys()) if sample2dog_info[s][1] != 'N']

mut = open(mgo_file, 'w')
print >>mut, ','.join(['GENE']+msamples)
for g in sorted(gene2summary.keys()):
        if g.startswith('ENSCAF') or '_CANFA' in g: continue
        print >>mut, ','.join([g] + [str(gene2summary[g][s]) for s in msamples])
mut.close()
