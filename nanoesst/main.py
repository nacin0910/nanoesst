import argparse
import os
import subprocess
import csv
import sys
import glob
import shutil

# 动态获取当前脚本 (main.py) 所在的目录
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(PACKAGE_DIR)
DB_DIR = os.path.join(ROOT_DIR, "db")

# 定义 ESKAPEE 对应的标识和数据库绝对路径
PATHOGENS = {
    "Enterococcus faecium": {
        "prefix": "efaecium", 
        "fasta": os.path.join(DB_DIR, "Enterococcus_faecium.fasta"), 
        "db": os.path.join(DB_DIR, "Enterococcus_faecium.db")
    },
    "Staphylococcus aureus": {
        "prefix": "saureus", 
        "fasta": os.path.join(DB_DIR, "Staphylococcus_aureus.fasta"), 
        "db": os.path.join(DB_DIR, "Staphylococcus_aureus.db")
    },
    "Klebsiella pneumoniae": {
        "prefix": "kpneumoniae", 
        "fasta": os.path.join(DB_DIR, "Klebsiella_pneumoniae.fasta"), 
        "db": os.path.join(DB_DIR, "Klebsiella_pneumoniae.db")
    },
    "Acinetobacter baumannii": {
        "prefix": "abaumannii", 
        "fasta": os.path.join(DB_DIR, "Acinetobacter_baumannii.fasta"), 
        "db": os.path.join(DB_DIR, "Acinetobacter_baumannii.db")
    },
    "Pseudomonas aeruginosa": {
        "prefix": "paeruginosa", 
        "fasta": os.path.join(DB_DIR, "Pseudomonas_aeruginosa.fasta"), 
        "db": os.path.join(DB_DIR, "Pseudomonas_aeruginosa.db")
    },
    "Escherichia coli": {
        "prefix": "ecoli", 
        "fasta": os.path.join(DB_DIR, "Escherichia_coli.fasta"), 
        "db": os.path.join(DB_DIR, "Escherichia_coli.db")
    }
}

def run_cmd(cmd):
    """运行系统命令并打印"""
    print(f"[RUNNING] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def process_sample(fastq_path, sample_name, algorithm, syldb, krakendb, threads, outdir, samtools_f_flag, samtools_q_flag, map_all):
    print(f"\n{'='*40}\nProcessing Sample: {sample_name}\n{'='*40}")
    
    # 动态定义各步骤的输出路径
    out_fastplong = os.path.join(outdir, '1.fastplong')
    out_sylph = os.path.join(outdir, '2.1.sylph_out')
    out_kraken = os.path.join(outdir, '2.2.kraken_out')
    out_minimap = os.path.join(outdir, '3.minimap2')
    out_pymlst = os.path.join(outdir, '4.pymlst')
    
    # 创建目录
    for d in [out_fastplong, out_minimap, out_pymlst]:
        os.makedirs(d, exist_ok=True)
    
    # 第二步：fastplong
    clean_fastq = os.path.join(out_fastplong, f"{sample_name}_clean.fastq.gz")
    if not os.path.exists(clean_fastq):
        cmd_fastplong = f"fastplong -i {fastq_path} -o {clean_fastq}"
        run_cmd(cmd_fastplong)
    
    detected_pathogens = set()

    # 第三步：物种注释 (sylph 或 kraken2+bracken)
    if algorithm == 'sylph':
        os.makedirs(out_sylph, exist_ok=True)
        profile_tsv = os.path.join(out_sylph, f"{sample_name}_profile.tsv")
        cmd_sylph = f"sylph profile {syldb} {clean_fastq} -t {threads} > {profile_tsv}"
        run_cmd(cmd_sylph)
        
        if os.path.exists(profile_tsv):
            with open(profile_tsv, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    contig = row.get('Contig_name', '')
                    try:
                        abundance = float(row.get('Taxonomic_abundance', 0))
                    except ValueError:
                        abundance = 0.0
                    
                    for pathogen in PATHOGENS.keys():
                        if pathogen in contig and abundance > 1.0:
                            detected_pathogens.add(pathogen)

    elif algorithm == 'kraken':
        os.makedirs(out_kraken, exist_ok=True)
        kreport = os.path.join(out_kraken, f"{sample_name}.kreport")
        kraken_out = os.path.join(out_kraken, f"{sample_name}.kraken")
        bracken_out = os.path.join(out_kraken, f"{sample_name}.bracken")
        
        cmd_kraken = f"kraken2 --db {krakendb} --threads {threads} {clean_fastq} --report {kreport} --output {kraken_out}"
        run_cmd(cmd_kraken)
        
        cmd_bracken = f"bracken -d {krakendb} -i {kreport} -o {bracken_out}"
        run_cmd(cmd_bracken)
        
        if os.path.exists(bracken_out):
            with open(bracken_out, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    name = row.get('name', '')
                    try:
                        fraction = float(row.get('fraction_total_reads', 0))
                    except ValueError:
                        fraction = 0.0
                    
                    for pathogen in PATHOGENS.keys():
                        if pathogen == name and fraction > 0.01:
                            detected_pathogens.add(pathogen)

    # 处理 --map-all 参数
    if map_all:
        print("[INFO] --map-all flag is active. Mapping to all 6 ESKAPEE reference genomes regardless of abundance.")
        detected_pathogens = set(PATHOGENS.keys())
    else:
        if not detected_pathogens:
            print(f"[INFO] No ESKAPEE pathogens with >1% abundance found in {sample_name}.")
            return

    # 第四步与第五步：比对和 MLST
    for pathogen in detected_pathogens:
        info = PATHOGENS[pathogen]
        prefix = info['prefix']
        ref_fasta = info['fasta']
        mlst_db = info['db']
        
        print(f"\n[INFO] Mapping to: {pathogen}")
        
        # Minimap2 比对，并使用用户自定义的 filtering 参数
        bam_file = os.path.join(out_minimap, f"{sample_name}_{prefix}_aligned.bam")
        cmd_minimap = f"minimap2 -ax map-ont {ref_fasta} {clean_fastq} -t {threads} | samtools view -b -F {samtools_f_flag} -q {samtools_q_flag} - > {bam_file}"
        run_cmd(cmd_minimap)
        
        # 提取比对上的 reads，并使用 pigz 多线程压缩
        mapped_fastq = os.path.join(out_minimap, f"{sample_name}_{prefix}.fastq.gz")
        cmd_extract = f"samtools fastq {bam_file} | pigz -p {threads} > {mapped_fastq}"
        run_cmd(cmd_extract)
        
        # Pymlst 分型
        result_txt = os.path.join(out_pymlst, f"{sample_name}_{prefix}_result.txt")
        cmd_mlst = f"claMLST search2 {mlst_db} {mapped_fastq} --single > {result_txt}"
        run_cmd(cmd_mlst)

def summarize_mlst(processed_samples, outdir):
    """仅将本次运行涉及的 sample 结果整合到 CSV 中"""
    if not processed_samples:
        return
        
    print(f"\n{'='*40}\nSummarizing MLST Results for Current Run\n{'='*40}")
    output_csv = os.path.join(outdir, "ST_summary.csv")
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(['Sample', 'Pathogen', 'ST', 'Genes', 'Type'])
        
        for sample_name in processed_samples:
            for pathogen, info in PATHOGENS.items():
                file_path = os.path.join(outdir, f"4.pymlst/{sample_name}_{info['prefix']}_result.txt")
                
                if not os.path.exists(file_path):
                    continue
                
                try:
                    with open(file_path, 'r', encoding='utf-8') as in_f:
                        lines = [line for line in in_f.read().split('\n') if line.strip()]
                        if len(lines) < 2:
                            continue 
                        
                        header = lines[0].rstrip('\r\n').split('\t')
                        data = lines[1].rstrip('\r\n').split('\t')
                        
                        if len(header) >= 2 and len(data) >= 2:
                            st_value = data[1]
                            genes = "-".join(header[2:])
                            alleles = "-".join(data[2:])
                            writer.writerow([sample_name, pathogen, st_value, genes, alleles])
                except Exception as e:
                    print(f"[ERROR] Could not parse {file_path}: {e}")
                
    print(f"[INFO] Summary successfully written to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Nanopore ESKAPEE processing and ST typing pipeline.")
    parser.add_argument('mode', choices=['process', 'batch'], help="Mode to run: 'process' for single file, 'batch' for directory.")
    parser.add_argument('-i', '--input', required=True, help="Input fastq.gz file (process) or directory (batch).")
    parser.add_argument('-n', '--name', required=True, help="Sample name (process) or path to mapping file (batch).")
    parser.add_argument('-a', '--algorithm', choices=['sylph', 'kraken'], required=True, help="Taxonomic classification algorithm: 'sylph' or 'kraken'.")
    parser.add_argument('-syldb', help="Path to sylph database (.syldb). Required if -a is sylph.")
    parser.add_argument('-krakendb', help="Path to kraken database. Required if -a is kraken.")
    parser.add_argument('-t', '--threads', default=16, type=int, help="Number of threads (default: 16).")
    parser.add_argument('-o', '--outdir', default='.', help="Output directory path (default: current directory).")
    parser.add_argument('-f', '--force', action='store_true', help="Force overwrite. Clears the output directory if it exists.")
    parser.add_argument('--map-all', action='store_true', help="Map reads to all 6 ESKAPEE references regardless of abundance.")
    parser.add_argument('-F', '--filter-flag', default=4, type=int, help="Samtools view -F flag (default: 4).")
    parser.add_argument('-q', '--min-mapq', default=10, type=int, help="Samtools view -q min MAPQ score (default: 10).")
    
    args = parser.parse_args()
    
    # 检查算法参数关联
    if args.algorithm == 'sylph' and not args.syldb:
        parser.error("-syldb is required when using sylph algorithm (-a sylph).")
    if args.algorithm == 'kraken' and not args.krakendb:
        parser.error("-krakendb is required when using kraken algorithm (-a kraken).")
    
    # 检查内置 db 文件夹是否存在
    if not os.path.exists(DB_DIR):
        print(f"[ERROR] Internal 'db' folder not found at {DB_DIR}. Please check your installation.")
        sys.exit(1)

    # 目录覆盖安全逻辑
    if args.force:
        if args.outdir != '.' and args.outdir != './':
            if os.path.exists(args.outdir):
                print(f"[WARNING] --force is set. Clearing directory: {args.outdir}")
                shutil.rmtree(args.outdir)
        else:
            # 如果是当前目录，仅清除 pipeline 产生的特定文件夹，防止误删用户输入文件
            print("[WARNING] --force is set in current dir. Safely clearing pipeline subfolders only.")
            for subdir in ['1.fastplong', '2.1.sylph_out', '2.2.kraken_out', '3.minimap2', '4.pymlst']:
                if os.path.exists(subdir): shutil.rmtree(subdir)
    
    os.makedirs(args.outdir, exist_ok=True)

    processed_samples = []

    if args.mode == 'process':
        process_sample(args.input, args.name, args.algorithm, args.syldb, args.krakendb, args.threads, args.outdir, args.filter_flag, args.min_mapq, args.map_all)
        processed_samples.append(args.name)
        
    elif args.mode == 'batch':
        mapping = {}
        with open(args.name, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]
        
        # 【关键修复】按 barcode 长度降序排序，确保 P13 优先于 P1 匹配
        sorted_mapping = sorted(mapping.items(), key=lambda x: len(x[0]), reverse=True)
        
        if not os.path.isdir(args.input):
            print(f"[ERROR] Input must be a directory in batch mode: {args.input}")
            sys.exit(1)
            
        for f_name in os.listdir(args.input):
            if f_name.endswith('.fastq.gz'):
                for barcode, s_name in sorted_mapping:
                    if barcode in f_name:
                        fastq_path = os.path.join(args.input, f_name)
                        process_sample(fastq_path, s_name, args.algorithm, args.syldb, args.krakendb, args.threads, args.outdir, args.filter_flag, args.min_mapq, args.map_all)
                        processed_samples.append(s_name)
                        break
                        
    # 总结本次运行的样本
    summarize_mlst(processed_samples, args.outdir)

if __name__ == '__main__':
    main()
