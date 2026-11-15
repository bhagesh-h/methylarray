# Builds the statistics fixtures in this folder: two tested models over 2,000
# simulated CpGs whose t statistics correlate (rho 0.6), ten strong hits in the
# smoker model, and an empty table for the untested model. Random values.
#   python3 make_fixture.py
import gzip, math, random
random.seed(7)
header = "cpg_id\tcontrast\tt\tp_value\n"
def p_of(t):
    return math.erfc(abs(t) / math.sqrt(2))
rows = {"score": [], "smoker": []}
for i in range(1, 2001):
    a, b = random.gauss(0, 1), random.gauss(0, 1)
    t_score = a
    t_smoker = 0.6 * a + 0.8 * b + (8 if i <= 10 else 0)
    cpg = "cg%08d" % i
    rows["score"].append("%s\tscore\t%.5g\t%.5g\n" % (cpg, t_score, p_of(t_score)))
    rows["smoker"].append("%s\tyes-no\t%.5g\t%.5g\n" % (cpg, t_smoker, p_of(t_smoker)))
for model, lines in rows.items():
    with gzip.GzipFile("ewas_%s_stats.tsv.gz" % model, "wb", mtime=0) as f:
        f.write((header + "".join(lines)).encode())
with gzip.GzipFile("ewas_small_stats.tsv.gz", "wb", mtime=0) as f:
    f.write(header.encode())
