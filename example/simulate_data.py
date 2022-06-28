#! /usr/bin/env python3
import msprime
import pathlib
from subprocess import run
import requests


bp_per_cm = int(1e6)
mu = 1e-8
chroms = list(range(1, 6))
chrlens = [(20 + i * 4) * bp_per_cm for i in range(1, 6)]

for i, (chrno, chrlen) in enumerate(zip(chroms, chrlens)):
    print(f"Simulating chr {chrno} via msprime")

    # simulate ancestry
    ts = msprime.sim_ancestry(
        samples=50,
        sequence_length=int(chrlen),
        recombination_rate=0.01 / bp_per_cm,
        population_size=10000,
        ploidy=2,
        random_seed=i + 1,
    )

    # simulate mutation
    mts = msprime.sim_mutations(ts, random_seed=i + 1, rate=mu, discrete_genome=True)

    # write vcf
    with open(f"{chrno}.vcf", "w") as f:
        mts.write_vcf(
            f, contig_id=f"{chrno}", individual_names=[f"ind{i:3d}" for i in range(50)]
        )

    # add a zero snp-density region
    cmd = f"""
        bgzip -f {chrno}.vcf; mv {chrno}.vcf.gz tmp_{chrno}.vcf.gz; bcftools index -f tmp_{chrno}.vcf.gz
        bcftools view -r {chrno}:1-5000000,{chrno}:8000000-{chrlen} tmp_{chrno}.vcf.gz -Oz -o {chrno}.vcf.gz
        rm tmp_{chrno}.vcf.gz*
    """
    run(cmd, shell=True)

    # write map file
    with open(f"{chrno}.map", "w") as f:
        bps = list(range(1, chrlen, chrlen // 100))
        cms = [bp / bp_per_cm for bp in bps]
        lines = []
        for bp, cm in zip(bps, cms):
            vid = "."
            line = f"{chrno} {vid} {cm} {bp}\n"  # sep=' '
            lines.append(line)
        f.writelines(lines)

    # download hapibd jar if needed
    hapibd_url = "https://faculty.washington.edu/browning/hap-ibd.jar"
    hapibd_jar = "hap-ibd.jar"
    if not pathlib.Path(hapibd_jar).exists():
        response = requests.get(hapibd_url)
        with open(hapibd_jar, "wb") as f:
            f.write(response.content)

    # run hapibd
    print("\tCalling IBD via hap-ibd")
    cmd = f"java -Xmx1g -jar {hapibd_jar} gt={chrno}.vcf.gz map={chrno}.map out={chrno}"
    _ = run(cmd, shell=True, capture_output=True)

    # clean up files
    pathlib.Path(f"{chrno}.log").unlink()
    pathlib.Path(f"{chrno}.hbd.gz").unlink()
