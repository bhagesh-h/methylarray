# Builds the simulated fixtures in this folder: 80 arrays (two failed sample QC)
# with clocks, cell lineages, two DNA methylation scores and phenotypes. The
# Horvath acceleration rises with BMI, metabolic syndrome is more common at high
# BMI, plate P3 shifts the Hannum acceleration, and the SNP outlier score rises
# with the Hannum acceleration, and ancestry PC1 with the Horvath acceleration.
# Random values.
#   python3 make_fixture.py
import random
random.seed(11)
n = 80
ids = ["S%03d" % i for i in range(1, n + 1)]
rows = []
for i, s in enumerate(ids):
    age = random.uniform(55, 80); men = random.random() < 0.5; bmi = random.gauss(30, 3)
    plate = "P%d" % (i % 4 + 1); chip = "2000%02d" % (i // 8); pos = "R0%dC01" % (i % 8 + 1)
    mets = "yes" if bmi + random.gauss(0, 2) > 30 else "no"
    gran = random.uniform(0.5, 0.7); mono = random.uniform(0.04, 0.1); b = random.uniform(0.03, 0.08); nk = random.uniform(0.02, 0.06)
    cd8 = random.uniform(0.05, 0.12); cd4 = max(0.01, 1 - gran - mono - b - nk - cd8)
    eaa_h = 0.6 * (bmi - 30) + random.gauss(0, 3)
    eaa_hn = random.gauss(0, 3) + (3 if plate == "P3" else 0)
    pace = 1 + 0.02 * (bmi - 30) + random.gauss(0, 0.08)
    rows.append(dict(sample=s, age=round(age, 1), men=men, bmi=round(bmi, 2), plate=plate, chip=chip, pos=pos, mets=mets,
                     gran=gran, mono=mono, b=b, nk=nk, cd8=cd8, cd4=cd4, eaa_h=eaa_h, eaa_hn=eaa_hn, pace=pace,
                     smoke=round(random.gauss(0, 1) + 0.3 * eaa_h / 3, 4), group=random.choice(["a", "b", "c"])))
def write(name, header, f):
    with open(name, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows: out.write("\t".join(str(v) for v in f(r)) + "\n")
write("samples.tsv", ["sample", "plate", "sentrix_id", "sentrix_position", "reported_sex", "reported_age"],
      lambda r: [r["sample"], r["plate"], r["chip"], r["pos"], "M" if r["men"] else "F", r["age"]])
write("samples_no_age.tsv", ["sample", "plate", "sentrix_id", "sentrix_position", "reported_sex", "reported_age"],
      lambda r: [r["sample"], r["plate"], r["chip"], r["pos"], "M" if r["men"] else "F", ""])
# SNP outlier scores from their own random stream (the other fixtures stay as they
# were): they rise with the Hannum acceleration, and three arrays lie just
# below the FAIL threshold of -4 (borderline contamination)
snp_random = random.Random(5)
with open("sample_qc.tsv", "w") as out:
    out.write("sample\tstatus\tsnp_outlier\n")
    for i, r in enumerate(rows):
        outlier = -4.2 if i in (10, 20, 30) else -5 + 0.04 * r["eaa_hn"] + snp_random.gauss(0, 0.08)
        out.write("%s\t%s\t%.3f\n" % (r["sample"], "FAIL" if i < 2 else "PASS", outlier))
write("clocks.tsv", ["sample", "age", "HorvathS2013", "HannumG2013", "DunedinPACE", "eaa_HorvathS2013", "eaa_HannumG2013"],
      lambda r: [r["sample"], r["age"], round(r["age"] + r["eaa_h"], 4), round(r["age"] + r["eaa_hn"], 4), round(r["pace"], 4), round(r["eaa_h"], 4), round(r["eaa_hn"], 4)])
write("cell_fractions.tsv", ["sample", "lineage_B", "lineage_CD4T", "lineage_CD8T", "lineage_NK", "lineage_Mono", "lineage_Gran"],
      lambda r: [r["sample"]] + [round(r[k], 5) for k in ("b", "cd4", "cd8", "nk", "mono", "gran")])
write("dnam_scores.tsv", ["sample", "dnam_smoking"], lambda r: [r["sample"], r["smoke"]])
write("phenotypes.tsv", ["sample", "bmi", "mets", "group"], lambda r: [r["sample"], r["bmi"], r["mets"], r["group"]])
# Ancestry PCs from their own random stream: PC1 rises with the Horvath acceleration
ancestry_random = random.Random(9)
with open("ancestry_pcs.tsv", "w") as out:
    out.write("sample\tancestry_pc1\tancestry_pc2\n")
    for r in rows:
        out.write("%s\t%.4f\t%.4f\n" % (r["sample"], 0.2 * r["eaa_h"] + ancestry_random.gauss(0, 0.8), ancestry_random.gauss(0, 1)))
# Genotype-robust clocks from their own random stream: the Horvath acceleration
# keeps most of its signal, the Hannum acceleration less
robust_random = random.Random(13)
with open("clocks_robust.tsv", "w") as out:
    out.write("sample\tage\tHorvathS2013\tHannumG2013\teaa_HorvathS2013\teaa_HannumG2013\n")
    for r in rows:
        h = 0.9 * r["eaa_h"] + robust_random.gauss(0, 1)
        hn = 0.5 * r["eaa_hn"] + robust_random.gauss(0, 2)
        out.write("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n" % (r["sample"], r["age"], r["age"] + h, r["age"] + hn, h, hn))
