from .misc import load_fasta 
import sys
import subprocess            

def mask_low_coverage_regions(assembly, reads, min_depth=1, threads=1, circularity=None):
    # align reads to assembly with minimap2
    assembly_base = assembly.split('/')[-1]
    work_dir = "/".join(assembly.split('/')[:-1])
    if work_dir == '':
        work_dir = '.'
    reads_base =reads.split('/')[-1]
    bamfile = assembly + "_" + reads_base + ".bam"
    sys.stderr.write("reads_base={}, bamfile(base)= {}\n".format(reads_base, bamfile))
    command = "minimap2 -t {} -a {} {} | samtools view -bS - | samtools sort -o {} -".format(threads, assembly, reads, bamfile)
    sys.stderr.write("run minimap: \n"+command+"\n")
    subprocess.run(command, shell=True)
     
    assembly_fasta = load_fasta(assembly)
    masked_assembly = {}
    original_assembly = {}
    for name, seq, _ in assembly_fasta:
        original_assembly[name] = seq
        masked_assembly[name] = ['N']*len(seq)

    # calculate depth at each position, mask positions below threshold
    command = ['samtools', 'depth', bamfile]
    sys.stderr.write("now find depth:\n"+" ".join(command)+"\n")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while process.poll() is None:
        contig, pos, depth = process.stdout.readline().rstrip().decode().split("\t")
        if pos:
            if int(depth) >= min_depth:
                pos = int(pos)-1 # 1-based positions ouput by samtools depth ??
                masked_assembly[contig][pos] = original_assembly[contig][pos]

    #print masked assembly to STDOUT
    for name, seq, _ in assembly_fasta:
        header = '>' + name
        if circularity and circularity[name]:
            header += ' circular=true'
        print(header)
        print(''.join(masked_assembly[name]))


if __name__ == '__main__':
    assembly, reads = sys.argv[1:3]
    min_depth = 5
    if len(sys.argv) > 3:
        min_depth = int(sys.argv[3])
    threads = 4
    if len(sys.argv) > 4:
        threads = int(sys.argv[4])
    sys.stderr.write("assembly={}, reads={}, min_depth={}, threads={}\n".format(assembly, reads, min_depth, threads))

    mask_low_coverage_regions(assembly, reads, min_depth=min_depth, threads=threads)
