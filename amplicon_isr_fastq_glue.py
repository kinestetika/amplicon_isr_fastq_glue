import argparse
import time
import gzip
import os
import subprocess
import tempfile
from pathlib import Path
from collections import namedtuple


VERSION = "0.1"
START_TIME = time.monotonic()
LOG_FILE = Path('log.txt')

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

FastqSeq = namedtuple('FastqSeq', ['id', 'seq', 'quality'])
Sample = namedtuple('Sample', ['name', 'fwd_read_file', 'rev_read_file'])

def parse_arguments():
    parser = argparse.ArgumentParser(description='amplicon_isr_fastq_glue.py. (C) Marc Strous, 2022')
    parser.add_argument('-m', required=True, help='file with sample/read mappings (same format as metaamp)')
    parser.add_argument('-o', required=True, help='file with oligos (same format as metaamp)')
    parser.add_argument('-n', required=True, help='name of the analysis')
    parser.add_argument('-t', default=150, help='trim reads at this length (default 150)')
    parser.add_argument('-f', default=False, action="store_true", help='to overwrite previous files and analyses')
    return parser.parse_args()


def log(log_message, values=()):
    if len(values):
        final_msg = f'{format_runtime()} {log_message.format(*values)}'
    else:
        final_msg = f'{format_runtime()} {log_message}'
    print(final_msg)
    try:
        with open(LOG_FILE, 'a') as log_handle:
            log_handle.write(final_msg)
            log_handle.write('\n')
    except FileNotFoundError:
        pass


def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def run_external(exec, stdin=None, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, log_cmd=True):
    if log_cmd:
        log(exec)
    result = subprocess.run(exec.split(), stdout=stdout, stdin=stdin, stderr=stderr)
    if result.returncode != 0:
        raise Exception(f'Fatal: Error while trying to run "{exec}"')


def parse_mapping_file(file: Path) -> list[Sample]:
    log(f'Now parsing sample definitions from mapping file "{file}"...')
    sample_list = []
    sample_names = set()
    fastq_files = set()
    with open(file) as reader:
        for line in reader:
            line = line.strip()
            if line.startswith('#'):
                continue
            words = line.split('\t')
            if len(words) < 5:
                print(f'Warning: Skipping malformed line {line}')
            sample = Sample(words[0], Path(file.parent / words[3]), file.parent / Path(words[4]))
            if sample.name in sample_names:
                raise Exception(f'Fatal: Sample "{sample.name}" has non-unique name')
            if not sample.fwd_read_file.exists():
                raise Exception(f'Fatal: Sample "{sample.name}" fwd read file "{sample.fwd_read_file}" does not exist')
            if not sample.rev_read_file.exists():
                raise Exception(f'Fatal: Sample "{sample.name}" rev read file "{sample.rev_read_file}" does not exist')
            if sample.fwd_read_file in fastq_files:
                raise Exception(f'Fatal: Sample "{sample.name}" has non-unique fwd fastq file "{sample.fwd_read_file}"')
            if sample.rev_read_file in fastq_files:
                raise Exception(f'Fatal: Sample "{sample.name}" has non-unique rev fastq file "{sample.rev_read_file}"')
            sample_names.add(sample.name)
            fastq_files.add(sample.fwd_read_file)
            fastq_files.add(sample.rev_read_file)
            sample_list.append(sample)
    if not sample_list:
        raise Exception(f'Fatal: Zero samples parsed from mapping file {file}!')
    log(f'Successfully parsed {len(sample_list)} sample definitions.')
    return sample_list


def parse_oligo_file(file: Path) -> tuple:
    log(f'Now parsing oligos/primers from oligo file "{file}"...')
    fwd_primer, rev_primer = ('', '')
    with open(file) as reader:
        for line in reader:
            line = line.strip()
            if line.startswith('#'):
                continue
            words = line.split()
            if 'forward' == words[0]:
                fwd_primer = words[1]
            if 'reverse' == words[0]:
                rev_primer = ''.join(COMPLEMENT.get(base, base) for base in reversed(words[1]))
    if not fwd_primer or not rev_primer:
        raise Exception(f'Fatal: Failed to parse primers from {file}: fwd = "{fwd_primer}", rev = "{rev_primer}"!')
    log(f'Successfully parsed primerset: {fwd_primer} + (reverse complemented) {rev_primer}.')
    return fwd_primer, rev_primer


def count_reads_in_fastq_file(file: Path) -> int:
    count = 0
    if file.name.endswith('.gz'):
        with gzip.open(file, 'rt') as input:
            for line in input:
                count += 1
    else:
        with open(file) as input:
            for line in input:
                count += 1
    return count // 4


def count_reads_in_fasta_file(file: Path) -> int:
    count = 0
    if file.name.endswith('.gz'):
        with gzip.open(file, 'rt') as input:
            for line in input:
                if line.startswith('>'):
                    count += 1
    else:
        with open(file) as input:
            for line in input:
                if line.startswith('>'):
                    count += 1
    return count


class PairedReadFastqParser:
    def __init__(self, path1, path2):
        self.path1 = path1
        self.path2 = path2
        self.handle1 = None
        self.handle2 = None

    def __enter__(self):
        if str(self.path1).endswith('.gz'):
            self.handle1 = gzip.open(self.path1, 'rt')
        else:
            self.handle1 = open(self.path1)

        if str(self.path2).endswith('.gz'):
            self.handle2 = gzip.open(self.path2, 'rt')
        else:
            self.handle2 = open(self.path2)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle1.close()
        self.handle2.close()

    def __iter__(self):
        while True:
            l1 = self.handle1.readline()
            l2 = self.handle1.readline()
            l3 = self.handle1.readline()
            l4 = self.handle1.readline()
            l5 = self.handle2.readline()
            l6 = self.handle2.readline()
            l7 = self.handle2.readline()
            l8 = self.handle2.readline()
            if l1 and l2 and l3 and l4 and l5 and l6 and l7 and l8:
                yield (FastqSeq(l1.strip()[1:], l2.strip(), l4.strip()),
                       FastqSeq(l5.strip()[1:], l6.strip(), l8.strip()))
            else:
                break


def write_fastq_seq(handle, fastq_seq: FastqSeq):
    handle.write(f'@')
    handle.write(fastq_seq.id)
    handle.write('\n')
    handle.write(fastq_seq.seq)
    handle.write('\n')
    handle.write('+\n')
    handle.write(fastq_seq.quality)
    handle.write('\n')


def trim_fastq(fastq_seq: FastqSeq, pos) -> FastqSeq:
    return FastqSeq(fastq_seq.id, fastq_seq.seq[:pos], fastq_seq.quality[:pos])


def rev_compl_fastq(fastq_seq: FastqSeq) -> FastqSeq:
    return FastqSeq(fastq_seq.id,
                    ''.join(COMPLEMENT.get(base, base) for base in reversed(fastq_seq.seq)),
                    fastq_seq.quality[::-1])


def glue_fastq(fq1: FastqSeq, fq2: FastqSeq) -> FastqSeq:
    return FastqSeq(fq1.id,
                    fq1.seq + fq2.seq,
                    fq1.quality + fq2.quality)


def main():
    print(f'This is amplicon_isr_fastq_glue.py {VERSION}')
    args = parse_arguments()
    sample_list = parse_mapping_file(Path(args.m).absolute())
    (fwd_primer, rev_primer) = parse_oligo_file(Path(args.o).absolute())
    output_folder = Path(os.getcwd()) / args.n
    merged_fastq_dir = output_folder / 'merged_fastq'
    merged_fastq_dir.mkdir(exist_ok=True, parents=True)

    merged_noprimer_fastq_dir = output_folder / 'no_primer_fastq'
    merged_noprimer_fastq_dir.mkdir(exist_ok=True, parents=True)
    merged_noprimer_filtered_fasta_dir = output_folder / 'filtered_fasta'
    merged_noprimer_filtered_fasta_dir.mkdir(exist_ok=True, parents=True)
    merged_noprimer_filtered_dereplicated_fasta_dir = output_folder / 'dereplicated_fasta'
    merged_noprimer_filtered_dereplicated_fasta_dir.mkdir(exist_ok=True, parents=True)
    otu_fasta_dir = output_folder / 'otu_fasta'
    otu_fasta_dir.mkdir(exist_ok=True, parents=True)
    otu_abundance_dir = output_folder / 'otu_abundances'
    otu_abundance_dir.mkdir(exist_ok=True, parents=True)
    abundance_summary_file = output_folder / f'otu_table.txt'

    for sample in sample_list:
        log(f'=================================================')
        log(f'Started processing sample {sample.name}, now glue-ing fastq files, while trimming to {args.t} bp...')
        merged_fastq_file = merged_fastq_dir / f'{sample.name}.fq.gz'
        if not merged_fastq_file.exists() or args.f:
            with PairedReadFastqParser(sample.fwd_read_file, sample.rev_read_file) as reader, \
                    gzip.open(merged_fastq_file, 'wt') as writer:
                count = 0
                for fq1, fq2 in reader:
                    count += 1
                    fq1 = trim_fastq(fq1, args.t)
                    fq2 = rev_compl_fastq(trim_fastq(fq2, args.t))
                    fq3 = glue_fastq(fq1, fq2)
                    write_fastq_seq(writer, fq3)
                log(f'Wrote {count} reads to {merged_fastq_file}')
        else:
            reads_remaining = count_reads_in_fastq_file(merged_fastq_file)
            log(f'Keeping previous results (use -f to overwrite). {reads_remaining} reads in fastq file')

        log('Now removing primer sequences from reads using cutadapt...')
        merged_noprimer_fastq_file = merged_noprimer_fastq_dir / f'{sample.name}.no-primers.fastq'
        if not merged_noprimer_fastq_file.exists() or args.f:
            run_external(f'cutadapt -e 0 --quiet -j 1 -a ^{fwd_primer}...{rev_primer} --discard-untrimmed '
                         f'-o {merged_noprimer_fastq_file} {merged_fastq_file}')
        else:
            log('Keeping previous results (use -f to overwrite).')
        log(f'{count_reads_in_fastq_file(merged_noprimer_fastq_file)} reads in fastq file after cutting primers.')

        log('Now filtering reads for quality...')
        merged_noprimer_filtered_fasta_file = merged_noprimer_filtered_fasta_dir / f'{sample.name}.filtered.fasta'
        if not merged_noprimer_filtered_fasta_file.exists() or args.f:
            run_external(f'usearch -fastq_filter {merged_noprimer_fastq_file} -fastq_maxee 1 -fastaout {merged_noprimer_filtered_fasta_file}')
        else:
            log('Keeping previous results (use -f to overwrite).')
        log(f'{count_reads_in_fasta_file(merged_noprimer_filtered_fasta_file)} '
            f'reads in fasta file after filtering for quality.')

        log('Now dereplicating reads...')
        merged_noprimer_filtered_dereplicated_fasta_file = merged_noprimer_filtered_dereplicated_fasta_dir / f'{sample.name}.dereplicated.fasta'
        if not merged_noprimer_filtered_dereplicated_fasta_file.exists() or args.f:
            run_external(f'usearch -fastx_uniques {merged_noprimer_filtered_fasta_file} -fastaout {merged_noprimer_filtered_dereplicated_fasta_file} -sizeout')
        else:
            log('Keeping previous results (use -f to overwrite).')
        log(f'{count_reads_in_fasta_file(merged_noprimer_filtered_dereplicated_fasta_file)} '
            f'reads in fasta file after dereplication.')

    log(f'=================================================')
    log('Now creating otu\'s by clustering...')
    all_dereplicated_reads_file = merged_noprimer_filtered_dereplicated_fasta_dir / f'all.dereplicated.fasta'
    with open(all_dereplicated_reads_file, 'w') as writer:
        for sample in sample_list:
            merged_noprimer_filtered_dereplicated_fasta_file = merged_noprimer_filtered_dereplicated_fasta_dir / f'{sample.name}.dereplicated.fasta'
            with open(merged_noprimer_filtered_dereplicated_fasta_file) as reader:
                for line in reader:
                    writer.write(line)
    otu_fasta_file = otu_fasta_dir / f'all.otus.fasta'
    if not otu_fasta_file.exists() or args.f:
        run_external(f'usearch -cluster_otus {all_dereplicated_reads_file} -otus {otu_fasta_file} -minsize 2 -relabel OTU_')
    else:
        log('Keeping previous results (use -f to overwrite).')
    log(f'{count_reads_in_fasta_file(otu_fasta_file)} otu\'s created.')

    log('Now determining read counts (abundances) for each otu and writing otu table...')

    otu_table = {}
    for sample in sample_list:
        otu_abundance_file = otu_abundance_dir / f'{sample.name}.abundances.uc'
        merged_noprimer_filtered_fasta_file = merged_noprimer_filtered_fasta_dir / f'{sample.name}.filtered.fasta'
        if not otu_abundance_file.exists() or args.f:
            run_external(f'usearch -usearch_global {merged_noprimer_filtered_fasta_file} -threads 8 -db {otu_fasta_file} -strand plus -id 0.97 -uc {otu_abundance_file}')
        else:
            log('Keeping previous results (use -f to overwrite).')
        with open(otu_abundance_file) as reader:
            for line in reader:
                line = line.strip()
                if line.endswith('*'):
                    continue
                words = line.split('\t')

                try:
                    sample_abundances = otu_table[words[-1]]
                except:
                    sample_abundances = {}
                    otu_table[words[-1]] = sample_abundances

                prev_count = sample_abundances.get(sample.name, 0)
                sample_abundances[sample.name] = prev_count + 1

    otu_names = sorted(list(otu_table.keys()))
    with open(abundance_summary_file, 'w') as writer:
        writer.write('#OTU_ID')
        for sample in sample_list:
            writer.write(f'\t{sample.name}')
        writer.write('\n')
        for otu_name in otu_names:
            writer.write(otu_name)
            for sample in sample_list:
                try:
                    otu_counts = otu_table[otu_name][sample.name]
                except:
                    otu_counts = 0
                writer.write(f'\t{otu_counts}')

if __name__ == "__main__":
    main()
